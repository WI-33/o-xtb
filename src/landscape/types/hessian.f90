! This file is part of xtb.
!
! Copyright (C) 2026 xtb contributors
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

!> Hessian providers of the xtb landscape task suite (PLAN D3).
!>
!> The Hessian is the only quantity that every landscape task needs
!> (RFO/PRFO, IRC and frequency analysis), so the way it is produced is
!> abstracted behind `landscape_hessian_provider`.  The settings field
!> `hessian_source` selects the implementation through
!> `new_landscape_hessian_provider`; an unknown name is reported with
!> `env%error` and no provider is returned (never a silent fallback).
!>
!> Available sources:
!>
!>  * `numerical` -- `landscape_numerical_hessian_provider`, a self-built
!>    central difference Hessian.  The finite difference step and the SCC
!>    accuracy of the displaced single points are **explicit parameters**
!>    (`step` argument, `sccacc` component); the suite must not rely on the
!>    `$hess` control block, whose numeric keys are silently ignored on the
!>    non-`--copy` path (`docs/_landscape-recon/07-*.md`, section 4).
!>    Every displaced geometry is evaluated from the same wavefunction, so
!>    the whole differentiation sequence uses one restart policy (07-*.md,
!>    section 3.3) and the Hessian is deterministic.
!>  * `model` -- `landscape_model_hessian_provider`, an approximate model
!>    Hessian (Lindh-D2, `src/model_hessian/`) used as a cheap initial guess
!>    for RFO/BFGS.  It reuses the existing xtb model-Hessian implementation
!>    with the default force constants of `set%mhset`
!>    (`src/setparam.f90:267-278`) and does not evaluate the calculator at all
!>    (`n_eval = 0`).
!>  * `analytical` - not implemented; an analytic GFN2 Hessian only has to add
!>    another extension of `landscape_hessian_provider` and a `hessian_source`
!>    name.
!>
!> The numerical provider can evaluate its 6N displaced points in `jobs`
!> concurrent worker processes (`--l*-hess-jobs`, default 1 = the serial loop).
!> A worker is a `fork()` of the running program (`xtb_landscape_process`), so
!> it inherits the calculator, the geometry and the reference wavefunction; it
!> evaluates one Cartesian column, hands the two gradients back through a pipe
!> and leaves through `_exit()`.  The differentiation arithmetic always runs in
!> the parent, so the matrix does not depend on the worker count - and a worker
!> pool inside a *threaded* run is refused rather than nested (`fork` inherits
!> the OpenMP runtime).
!>
!> The `file` source reads the Hessian written by the existing xtb
!> implementation (`hessian`, TURBOMOLE layout) back in
!> (`landscape_read_hessian_file`), so an expensive numerical Hessian of a
!> preceding run can be reused by every landscape task.  The reader is strict
!> on purpose: a missing file, a wrong header, a truncated or oversized matrix
!> and an asymmetric or non-finite matrix are reported with `env%error`
!> instead of being silently replaced by zeros - `rdhess`
!> (`src/freq/io.f90:82`), which the xtb subcommands use, cannot do that: it
!> reads without `iostat`, so a damaged file terminates the program, and a
!> file without a `$hessian` header yields a silently zero matrix.
!>
!> Units: coordinates in bohr, energies in hartree, gradients in hartree/bohr
!> and Hessian matrix elements in hartree/bohr^2.  The returned Hessian is
!> dense and exactly symmetric.
module xtb_landscape_hessian
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
!$ use omp_lib, only : omp_set_num_threads, omp_get_max_threads
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   use xtb_type_restart, only : TRestart
   use xtb_type_calculator, only : TCalculator
   use xtb_type_data, only : scc_results
   use xtb_type_setvar, only : modhess_setvar
   use xtb_modelhessian_type, only : TModelHessian
   use xtb_modelhessian_lindh, only : newLindhD2ModelHessian
   use xtb_lin, only : lin
   use xtb_setparam, only : get_namespace
   use xtb_landscape_settings, only : landscape_log_failed, landscape_value_present, &
      & landscape_int_to_string, landscape_real_to_string
   use xtb_landscape_settings, only : landscape_lowercase
   use, intrinsic :: iso_c_binding, only : c_int, c_long
   use xtb_landscape_process, only : landscape_process_available, &
      & landscape_process_pipe, landscape_process_fork, landscape_process_close, &
      & landscape_process_read_int, landscape_process_write_int, &
      & landscape_process_read_real, landscape_process_write_real, &
      & landscape_process_wait, landscape_process_status_text, &
      & landscape_process_worker_exit
   implicit none(type, external)
   private

   public :: landscape_hessian_provider
   public :: landscape_numerical_hessian_provider
   public :: landscape_model_hessian_provider
   public :: landscape_hessian_file
   public :: landscape_hessian_source_known
   public :: new_landscape_hessian_provider
   public :: landscape_read_hessian_file, landscape_write_hessian_file

   !> Print level of the displaced single points of the numerical Hessian.
   !> Negative values silence the calculator, exactly like the numerical
   !> Hessian of xtb itself (`get_gradient_derivs` in
   !> `src/type/calculator.f90`).
   integer, parameter :: displacement_printlevel = -1

   !> Header record of a TURBOMOLE Hessian file (`wrhess`, `src/freq/io.f90`)
   character(len=*), parameter :: hessian_header = '$hessian'

   !> Record format of a TURBOMOLE Hessian matrix, identical to `wrhess`
   !> (`src/freq/io.f90:53`): five columns per record, prefixed by five blanks
   character(len=*), parameter :: hessian_format = '(5x,5f15.10)'

   !> Tolerance of the symmetry and finiteness check of a Hessian read from a
   !> file (hartree/bohr^2, relative to the largest matrix element)
   real(wp), parameter :: hessian_symmetry_tol = 1.0e-8_wp

   !> Maximum number of values parsed from one record of a Hessian file
   integer, parameter :: hessian_max_record = 256

   !> Abstract Hessian provider.
   !>
   !> `name` repeats the `hessian_source` string that selected the provider
   !> (for reports), `exact` is true for analytic implementations.
   type, abstract :: landscape_hessian_provider
      !> Provider name, one of the `hessian_source` values
      character(len=:), allocatable :: name
      !> True if the provider returns an exact (analytic) Hessian
      logical :: exact = .false.
      !> Number of energy/gradient evaluations performed by the last `compute`
      !> (0 for providers that do not displace the geometry); the engine adds
      !> it to its cost counter so the reported cost is not underestimated.
      integer :: n_eval = 0
      !> Print level of the progress lines of the provider: 0 silences them,
      !> anything above prints one line **before** an expensive differentiation
      !> starts.  It is set from the task settings by
      !> `new_landscape_hessian_provider`, so the same gate that echoes the
      !> steps of a run also explains the pause of a Hessian build (the first
      !> one and every periodic rebuild).  Providers that cost no evaluation
      !> (`model`, `file`) print nothing at all.
      integer :: printlevel = 0
   contains
      !> Compute a dense, symmetric Hessian (hartree/bohr^2)
      procedure(hessian_compute), deferred :: compute
   end type landscape_hessian_provider

   abstract interface
      subroutine hessian_compute(self, env, mol, chk, calc, hess, step, odlr, final_err)
         import :: landscape_hessian_provider, TEnvironment, TMolecule, &
            & TRestart, TCalculator, wp
         !> Provider instance
         class(landscape_hessian_provider), intent(inout) :: self
         !> Computational environment
         type(TEnvironment), intent(inout) :: env
         !> Molecular structure, restored on return
         type(TMolecule), intent(inout) :: mol
         !> Wavefunction data
         type(TRestart), intent(inout) :: chk
         !> Single point calculator
         class(TCalculator), intent(inout) :: calc
         !> Dense Hessian (3N, 3N) in hartree/bohr^2, symmetrised on return
         real(wp), intent(inout) :: hess(:, :)
         !> Finite difference step (bohr)
         real(wp), intent(in) :: step
         !> Use the ODLR approximation
         logical, intent(in) :: odlr
         !> Residual error estimate of the approximation (0 for exact data)
         real(wp), intent(out) :: final_err
      end subroutine hessian_compute
   end interface

   !> Self-built central difference Hessian.
   !>
   !> `hess(:, 3(i-1)+c) = (g(x + step*e_ic) - g(x - step*e_ic)) / (2 step)`.
   !> The calculator accuracy is set to `sccacc` for the displaced geometries
   !> and restored afterwards; the wavefunction of every displacement is reset
   !> from the incoming state, so the sequence never mixes restart policies.
   !>
   !> Cost: 6N energy/gradient evaluations (`src/hessian.F90`, PLAN D3).
   type, extends(landscape_hessian_provider) :: landscape_numerical_hessian_provider
      !> SCC accuracy used while building the Hessian (explicit parameter)
      real(wp) :: sccacc = 0.3_wp
      !> Number of processes that evaluate the displaced single points
      !> concurrently (`--l*-hess-jobs`).  The 3N Cartesian columns are
      !> independent: every point starts from the reference wavefunction, so
      !> the result cannot depend on the order in which the columns are
      !> evaluated and the arithmetic of each column is done by the parent
      !> either way.  `1` keeps the serial loop, which is the default and never
      !> touches the process module.
      integer :: jobs = 1
      !> Number of builds this provider instance has started; the progress
      !> line uses it to say `computing` for the first build and `rebuilding`
      !> for the periodic rebuild of the optimiser, the saddle search and the
      !> reaction path.
      integer :: n_builds = 0
   contains
      !> Central difference Hessian
      procedure :: compute => numerical_hessian_compute
   end type landscape_numerical_hessian_provider

   !> Model Hessian provider (Lindh-D2 initial guess, PLAN D3).
   !>
   !> Reuses the xtb model Hessian of `src/model_hessian/lindh.f90`
   !> (`newLindhD2ModelHessian`, the parameterisation selected by the default
   !> `set%mhset%model = p_modh_old`) with the default force constants
   !> (`src/setparam.f90:267-278`).  The construction is purely algebraic, so
   !> the provider never touches the calculator and reports `n_eval = 0`.
   type, extends(landscape_hessian_provider) :: landscape_model_hessian_provider
   contains
      !> Lindh-D2 model Hessian
      procedure :: compute => model_hessian_compute
   end type landscape_model_hessian_provider

   !> Hessian read from a file (`hessian_source = 'file'`).
   !>
   !> The file is the `hessian` file of xtb itself (TURBOMOLE layout, see
   !> `landscape_read_hessian_file`), which makes the result of a preceding
   !> numerical Hessian reusable without any algorithm knowing about it
   !> (PLAN section 5.2).  The provider does not evaluate the calculator
   !> (`n_eval = 0`) and marks the Hessian as exact data.
   type, extends(landscape_hessian_provider) :: landscape_hessian_file
      !> Name of the file holding the Hessian
      character(len=:), allocatable :: filename
   contains
      !> Hessian from the file
      procedure :: compute => file_hessian_compute
   end type landscape_hessian_file

contains

!> True if `source` names a known Hessian provider.  The name is matched case
!> insensitively.
pure function landscape_hessian_source_known(source) result(known)
   !> Value of the `hessian_source` setting
   character(len=*), intent(in) :: source
   !> True if the name can be resolved by `new_landscape_hessian_provider`
   logical :: known

   select case (landscape_lowercase(source))
   case ('numerical', 'model', 'file')
      known = .true.
   case default
      known = .false.
   end select

end function landscape_hessian_source_known

!> Construct the Hessian provider selected by `source`.
!>
!> On return `provider` is allocated if and only if the name is known;
!> an unknown name is reported with `env%error` and the caller must abort
!> (there is no default provider on this path, a silent fallback would make
!> the configured `hessian_source` a lie).  The name is matched case
!> insensitively, like `landscape_hessian_source_known`.
subroutine new_landscape_hessian_provider(env, source, provider, sccacc, filename, &
   & printlevel, jobs)
   character(len=*), parameter :: source_name = 'new_landscape_hessian_provider'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Value of the `hessian_source` setting
   character(len=*), intent(in) :: source
   !> Selected provider, unallocated for an unknown name
   class(landscape_hessian_provider), allocatable, intent(out) :: provider
   !> SCC accuracy of the numerical Hessian, defaults to xtb `accu_hess`
   real(wp), intent(in), optional :: sccacc
   !> File of the `file` source; required for that source, unused otherwise
   character(len=*), intent(in), optional :: filename
   !> Print level of the provider progress lines, 0 (silent) when absent
   integer, intent(in), optional :: printlevel
   !> Number of concurrent worker processes of the numerical Hessian, 1 when
   !> absent
   integer, intent(in), optional :: jobs

   real(wp) :: scc_accuracy
   integer :: n_jobs, n_threads
   character(len=200) :: buffer

   scc_accuracy = 0.3_wp
   if (present(sccacc)) scc_accuracy = sccacc
   n_jobs = 1
   if (present(jobs)) n_jobs = jobs

   if (n_jobs < 1) then
      call env%error("the number of concurrent Hessian workers must be at least 1", &
         & source_name)
      return
   end if

   ! H2 of report 220 section 4: `--l*-hessian-file` is accepted, echoed by the
   ! frequency task and silently ignored by every source but `file` (measured:
   ! even a file name that does not exist left the run at exit 0).  The name is
   ! reported here - the one place every task builds its provider through - so
   ! the notice cannot be forgotten by a task; an unknown source is refused
   ! below and gets no second message.
   if (present(filename)) then
      if (len_trim(filename) > 0 .and. landscape_hessian_source_known(source)) then
         if (landscape_lowercase(source) /= 'file') then
            call env%warning("--l*-hessian-file is ignored: the Hessian source '"// &
               & landscape_lowercase(source)//"' reads no file (only "// &
               & "'--l*-hessian file' does); the name has no effect on this run", &
               & source_name)
         end if
      end if
   end if

   select case (landscape_lowercase(source))
   case ('numerical')
      ! Refuse a worker count this build cannot deliver instead of ignoring it;
      ! `1` needs no process creation and always works.
      if (n_jobs > 1) then
         if (landscape_process_available() == 0_c_int) then
            call env%error("concurrent numerical Hessian workers (more than one "//&
               & "job) are not available in this build: process creation is "//&
               & "missing on this platform", source_name)
            return
         end if
         ! A worker inherits the OpenMP runtime of this process, and libgomp has
         ! no `fork` handler: a worker that inherited a team of more than one
         ! thread blocks in `gomp_team_barrier_wait_end` at its first parallel
         ! region (measured: `OMP_NUM_THREADS=8`, 4 workers, 30 atoms, hang;
         ! `OPENBLAS_NUM_THREADS` makes no difference).  A worker that pinned
         ! itself to one thread would not hang, but it would then also evaluate
         ! its points with one thread while the serial run of the same job
         ! count uses the ambient team - a difference at the level of
         ! `AGENTS.md` section 6.1 that costs the byte identity this option is
         ! required to keep.  The combination is therefore refused and not
         ! silently serialised or silently re-threaded: `jobs > 1` needs a
         ! **single-threaded run**, and `OMP_NUM_THREADS=1` is how that is said
         ! (a `-P`/`--parallel` value is checked here as well, because it sets
         ! the same team size).
         n_threads = 1
!$       n_threads = omp_get_max_threads()
         if (n_threads > 1) then
            write(buffer, '(a,i0,a,i0,a)') "--l*-hess-jobs ", n_jobs, &
               & " needs a single-threaded run (OMP_NUM_THREADS=1, and -P 1 "//&
               & "if -P/--parallel is given), but this run uses ", &
               & n_threads, " threads; set OMP_NUM_THREADS=1 or use "//&
               & "--l*-hess-jobs 1"
            call env%error(trim(buffer), source_name)
            return
         end if
      end if
      allocate(landscape_numerical_hessian_provider :: provider)
      select type (hp => provider)
      type is (landscape_numerical_hessian_provider)
         hp%sccacc = scc_accuracy
         hp%jobs = n_jobs
      end select
      provider%name = 'numerical'
      provider%exact = .false.
   case ('model')
      allocate(landscape_model_hessian_provider :: provider)
      provider%name = 'model'
      provider%exact = .false.
   case ('file')
      if (.not. present(filename)) then
         call env%error("the Hessian source 'file' requires a file name", source_name)
         return
      end if
      if (len_trim(filename) == 0) then
         call env%error("the Hessian source 'file' requires a non-empty file name", &
            & source_name)
         return
      end if
      allocate(landscape_hessian_file :: provider)
      select type (hp => provider)
      type is (landscape_hessian_file)
         hp%filename = trim(filename)
      end select
      provider%name = 'file'
      provider%exact = .true.
   case default
      call env%error("unknown hessian source: "//trim(source), source_name)
   end select

   if (present(printlevel) .and. allocated(provider)) then
      provider%printlevel = printlevel
   end if

end subroutine new_landscape_hessian_provider

!> Central difference Hessian of the attached calculator.
!>
!> The incoming geometry and wavefunction define the reference state of the
!> differentiation: every displaced single point starts from that state and the
!> calculator accuracy is swapped for `self%sccacc` for the whole sequence.
!> The geometry of `mol` and the accuracy of `calc` are restored before the
!> routine returns, also on the error path.
!>
!> The 3N Cartesian columns are independent of each other and of their order:
!> a column is one `chk%copy(wfn_ref)` plus two displaced single points, and
!> nothing of one column reaches the next.  That is what allows the columns to
!> be evaluated in `self%jobs` concurrent worker processes (`jobs > 1`,
!> `numerical_hessian_concurrent`) and it is why the resulting Hessian is the
!> same matrix for every worker count: the workers only move the *evaluation*
!> of the columns around, the differentiation arithmetic runs in this process,
!> one column at a time, exactly as in the serial loop.
subroutine numerical_hessian_compute(self, env, mol, chk, calc, hess, step, odlr, final_err)
   character(len=*), parameter :: source = 'landscape_numerical_hessian_compute'
   !> Provider instance
   class(landscape_numerical_hessian_provider), intent(inout) :: self
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure, restored on return
   type(TMolecule), intent(inout) :: mol
   !> Wavefunction data
   type(TRestart), intent(inout) :: chk
   !> Single point calculator
   class(TCalculator), intent(inout) :: calc
   !> Dense Hessian (3N, 3N), symmetrised on return
   real(wp), intent(inout) :: hess(:, :)
   !> Finite difference step (bohr)
   real(wp), intent(in) :: step
   !> Use the ODLR approximation
   logical, intent(in) :: odlr
   !> Residual error estimate (always zero for the full central difference)
   real(wp), intent(out) :: final_err

   real(wp), allocatable :: x0(:, :), xyz_plus(:, :), xyz_minus(:, :)
   real(wp), allocatable :: g_plus(:, :), g_minus(:, :)
   type(TRestart) :: wfn_ref
   type(scc_results) :: results
   real(wp) :: accuracy, egap, energy_plus, energy_minus
   real(wp) :: sigma(3, 3)
   integer :: icol, n
   logical :: failed
   character(len=160) :: buffer

   final_err = 0.0_wp
   self%n_eval = 0
   n = mol%n
   if (n <= 0) return

   if (odlr) then
      call env%error("the numerical Hessian provider always builds the full "//&
         & "dense Hessian and does not support the ODLR approximation", source)
      return
   end if
   ! A non-finite step passes the comparison below (`nan <= 0` is false) and
   ! reaches the displaced single points: the CLI parser accepts `nan`/`inf`
   ! for the Hessian step flags, so the value is named and refused here.
   if (.not. ieee_is_finite(step)) then
      write(buffer, '(a,es15.6,a)') "the finite difference step of the "//&
         & "numerical Hessian (--l*-hess-step) is ", step, ", not a finite number"
      call env%error(trim(buffer), source)
      return
   end if
   if (step <= 0.0_wp) then
      call env%error("the finite difference step of the numerical Hessian "//&
         & "must be positive", source)
      return
   end if

   allocate(x0(3, n), xyz_plus(3, n), xyz_minus(3, n))
   allocate(g_plus(3, n), g_minus(3, n))
   x0 = mol%xyz
   call wfn_ref%copy(chk)
   accuracy = calc%accuracy
   calc%accuracy = self%sccacc
   egap = 0.0_wp
   sigma = 0.0_wp

   ! Announce the pause before it starts: the loop below is 2 displaced single
   ! points per Cartesian component, i.e. 6N evaluations, all of them at
   ! `displacement_printlevel` and therefore silent.  The number is the exact
   ! cost of this implementation, not a guessed constant, and the periodic
   ! rebuild of a caller announces itself the same way (`n_builds > 0`).  A
   ! concurrent build names its worker count as well, so the setting is
   ! visible in the transcript and not only in the wall time; the serial line
   ! (the default) is unchanged.
   if (self%printlevel > 0) then
      if (self%n_builds == 0) then
         write(env%unit, '(1x,a,a,a,i0,a,es9.2,a)', advance='no') "computing the ", &
            & self%name, " Hessian: ", 6*n, &
            & " displaced single points (6N, step ", step, " bohr"
      else
         write(env%unit, '(1x,a,a,a,i0,a,es9.2,a)', advance='no') "rebuilding the ", &
            & self%name, " Hessian: ", 6*n, &
            & " displaced single points (6N, step ", step, " bohr"
      end if
      if (self%jobs > 1) then
         write(env%unit, '(a,i0,a)') ", ", self%jobs, " worker processes)"
      else
         write(env%unit, '(a)') ")"
      end if
   end if
   self%n_builds = self%n_builds + 1

   if (self%jobs > 1) then
      call numerical_hessian_concurrent(failed)
   else
      failed = .false.
      do icol = 1, 3*n
         call numerical_hessian_column(icol, failed)
         if (failed) exit
         hess(:, icol) = hess(:, icol) &
            & + reshape(g_plus - g_minus, [3*n])/(2.0_wp*step)
      end do
   end if

   call numerical_hessian_cleanup()
   if (failed) return
   hess = 0.5_wp*(hess + transpose(hess))

contains

   !> Evaluate one Cartesian column: the two displaced geometries of component
   !> `icol` (`3*(iat-1)+ic`, ascending atom by atom), each started from the
   !> reference wavefunction.  `g_plus` and `g_minus` hold the two gradients on
   !> return, `failed` is true when a single point logged an error.
   !>
   !> This is the only place a displaced point is evaluated, for the serial and
   !> the concurrent driver alike, so the two cannot drift apart; the
   !> differentiation arithmetic stays with the callers.
   subroutine numerical_hessian_column(icol, failed)
      !> Cartesian column, `1 .. 3N`
      integer, intent(in) :: icol
      !> True when a displaced single point reported an error
      logical, intent(out) :: failed

      integer :: iat, ic, nlog

      iat = (icol - 1)/3 + 1
      ic = icol - 3*(iat - 1)
      xyz_plus = x0
      xyz_plus(ic, iat) = x0(ic, iat) + step
      xyz_minus = x0
      xyz_minus(ic, iat) = x0(ic, iat) - step

      failed = .true.
      nlog = env%nLog
      mol%xyz = xyz_plus
      call chk%copy(wfn_ref)
      call calc%singlepoint(env, mol, chk, displacement_printlevel, .true., &
         & energy_plus, g_plus, sigma, egap, results)
      self%n_eval = self%n_eval + 1
      if (landscape_log_failed(env, nlog)) return

      mol%xyz = xyz_minus
      call chk%copy(wfn_ref)
      call calc%singlepoint(env, mol, chk, displacement_printlevel, .true., &
         & energy_minus, g_minus, sigma, egap, results)
      self%n_eval = self%n_eval + 1
      if (landscape_log_failed(env, nlog)) return

      failed = .false.
   end subroutine numerical_hessian_column

   !> Evaluate the 3N Cartesian columns in `self%jobs` concurrent worker
   !> processes.
   !>
   !> Every worker is a `fork()` of this process (`xtb_landscape_process`), so
   !> it inherits the calculator, the geometry and the reference wavefunction
   !> and pays no start-up cost; it evaluates its column with
   !> `numerical_hessian_column`, writes the two gradients into its own pipe
   !> and leaves through `_exit()`.  This process never writes to a worker's
   !> pipe and a worker never writes to a file: the only shared object is the
   !> pipe of that worker.
   !>
   !> Results are collected in launch order and applied column by column, so
   !> the collection order does not reach the matrix: each column is written
   !> exactly once.  A worker reports the number of points it attempted
   !> through one status word, which keeps the evaluation counter identical to
   !> the serial path (it counts a point that logged an error, like the serial
   !> loop).  One failed column fails the whole Hessian, exactly as in the
   !> serial loop; the workers still running are then closed and reaped, never
   !> detached.
   subroutine numerical_hessian_concurrent(failed)
      !> True when a worker reported a failure
      logical, intent(out) :: failed

      integer, allocatable :: q_fd(:), q_col(:)
      integer(c_int), allocatable :: q_pid(:)
      integer(c_int) :: fds(2), worker_pid, wait_status
      integer(c_int) :: report(1), status_word(1)
      integer(c_long) :: want, got, put, bytes
      integer :: head, count, next_col, icol, slot, iat, ic, ndone
      integer :: before, gone, ignored
      logical :: worker_failed

      failed = .false.
      allocate(q_fd(self%jobs), q_col(self%jobs), q_pid(self%jobs))
      head = 1
      count = 0
      next_col = 1
      ! One iteration starts **one** worker while the queue has room and
      ! collects one otherwise, so at most `self%jobs` workers exist.
      do while (.not. failed .and. (count > 0 .or. next_col <= 3*n))
         if (count < self%jobs .and. next_col <= 3*n) then
            slot = mod(head + count - 1, self%jobs) + 1
            if (landscape_process_pipe(fds) /= 0_c_int) then
               call env%error("the numerical Hessian could not create a worker "//&
                  & "pipe", source)
               failed = .true.
               cycle
            end if
            ! An inherited standard unit must not be re-flushed by a worker
            ! that dies abnormally (a Fortran runtime trace flushes): the
            ! parent's buffer is copied into the worker by `fork`, so it is
            ! emptied here first and a duplicated transcript is impossible.
            flush(env%unit)
            worker_pid = landscape_process_fork()
            if (worker_pid < 0_c_int) then
               ignored = int(landscape_process_close(fds(1)))
               ignored = int(landscape_process_close(fds(2)))
               call env%error("the numerical Hessian could not create a worker "//&
                  & "process; use --l*-hess-jobs 1", source)
               failed = .true.
               cycle
            end if
            if (worker_pid == 0_c_int) then
               ! This process is the worker: it evaluates its column, reports
               ! it through its own pipe and leaves.  The `_exit` below is
               ! final, so nothing of this frame ever returns to the caller.
               !
               ! A worker runs **one** thread.  The provider factory refuses a
               ! worker pool inside a threaded run (libgomp has no `fork`
               ! handler, and a one-thread worker could not reproduce the
               ! ambient team's arithmetic bit for bit); this pin states the
               ! invariant a second time, so a worker stays single-threaded
               ! even if the calculator changed the team size during the run.
!$             call omp_set_num_threads(1)
               !
               ! The read ends of this worker's siblings are inherited copies
               ! of the parent's; they are dropped here so that no leftover
               ! descriptor keeps another worker's pipe alive.
               do gone = 1, count
                  report = int(landscape_process_close(&
                     & q_fd(mod(head + gone - 2, self%jobs) + 1)))
               end do
               report = int(landscape_process_close(fds(1)))
               before = self%n_eval
               call numerical_hessian_column(next_col, worker_failed)
               if (worker_failed) then
                  ! The status word of a failed worker is the number of
                  ! displaced points it attempted, so the parent's counter
                  ! matches the serial path, which increments after the failing
                  ! call as well.
                  report = int(self%n_eval - before, c_int)
               else
                  report = 0_c_int
               end if
               bytes = int(storage_size(report(1))/8, c_long)
               put = landscape_process_write_int(fds(2), report, bytes)
               if (report(1) == 0_c_int) then
                  bytes = int(storage_size(0.0_wp)/8, c_long)*int(3*n, c_long)
                  put = landscape_process_write_real(fds(2), g_plus, bytes)
                  put = landscape_process_write_real(fds(2), g_minus, bytes)
               end if
               report = int(landscape_process_close(fds(2)))
               call landscape_process_worker_exit(0_c_int)
            end if
            ! This process is the parent: it keeps the read end only.  With the
            ! write end open here the read would never see the end of file of a
            ! worker that died before it reported.
            ignored = int(landscape_process_close(fds(2)))
            q_fd(slot) = fds(1)
            q_pid(slot) = worker_pid
            q_col(slot) = next_col
            count = count + 1
            next_col = next_col + 1
            cycle
         end if

         ! Collect the oldest worker: a blocking read of its status word, then
         ! of its two gradients.  One column is `48N + 4` bytes - two `3N`
         ! gradients of eight bytes each plus the four byte status word - so at
         ! N = 1366 atoms the payload exceeds the 64 KiB pipe buffer of Linux
         ! and the worker's write blocks until the parent has drained the pipe.
         ! That is the documented limit of this collector (N <= 1365 keeps one
         ! column in one write), not a deadlock: every pipe has exactly one
         ! writer, the parent closes the write end before it launches the next
         ! worker, and it always blocks on the pipe of the **oldest** worker and
         ! loops over short reads - so a blocked writer is never the pipe being
         ! waited for, and every worker makes progress in launch order.
         icol = q_col(head)
         want = int(storage_size(status_word(1))/8, c_long)
         got = landscape_process_read_int(q_fd(head), status_word, want)
         if (got == want) then
            if (status_word(1) == 0_c_int) then
               want = int(storage_size(0.0_wp)/8, c_long)*int(3*n, c_long)
               got = landscape_process_read_real(q_fd(head), g_plus, want)
               if (got == want) then
                  got = landscape_process_read_real(q_fd(head), g_minus, want)
               end if
               if (got /= want) status_word(1) = -1_c_int
            end if
         else
            status_word(1) = -1_c_int
         end if
         ignored = int(landscape_process_close(q_fd(head)))
         ignored = int(landscape_process_wait(q_pid(head), wait_status))
         head = mod(head, self%jobs) + 1
         count = count - 1

         if (status_word(1) == 0_c_int .and. wait_status == 0_c_int) then
            self%n_eval = self%n_eval + 2
            ! One column is written exactly once and the differentiation of a
            ! column is one expression, so the order in which the workers
            ! finish cannot reach the matrix.
            hess(:, icol) = hess(:, icol) &
               & + reshape(g_plus - g_minus, [3*n])/(2.0_wp*step)
         else
            failed = .true.
            iat = (icol - 1)/3 + 1
            ic = icol - 3*(iat - 1)
            if (status_word(1) > 0_c_int) then
               ndone = int(status_word(1))
               self%n_eval = self%n_eval + ndone
               write(buffer, '(a,i0,a,i0,a)') "a displaced single point of the "//&
                  & "numerical Hessian (atom ", iat, ", Cartesian component ", &
                  & ic, ") failed in its worker process; the Hessian is not "//&
                  & "available"
               call env%error(trim(buffer), source)
            else
               write(buffer, '(a)') "the worker process of the numerical "//&
                  & "Hessian failed (atom "//trim(landscape_int_to_string(iat))//&
                  & ", Cartesian component "//trim(landscape_int_to_string(ic))//&
                  & ", "//landscape_process_status_text(wait_status)//&
                  & ", status word "//trim(landscape_int_to_string(&
                     & int(status_word(1))))//")"
               call env%error(trim(buffer), source)
            end if
         end if
      end do

      ! Close and reap every worker that is still running; none is detached.
      do while (count > 0)
         ignored = int(landscape_process_close(q_fd(head)))
         ignored = int(landscape_process_wait(q_pid(head), wait_status))
         head = mod(head, self%jobs) + 1
         count = count - 1
      end do

   end subroutine numerical_hessian_concurrent

   !> Restore the geometry, the calculator accuracy and the wavefunction of
   !> the reference state.  The wavefunction is restored because the last
   !> displaced evaluation above (and the `chk%copy` inside the loop) would
   !> otherwise leak a displaced SCF state back to the caller, which matters
   !> for the `warm` restart policy.
   subroutine numerical_hessian_cleanup()
      mol%xyz = x0
      calc%accuracy = accuracy
      call chk%copy(wfn_ref)
   end subroutine numerical_hessian_cleanup

end subroutine numerical_hessian_compute

!> Lindh-D2 model Hessian (initial guess for RFO/BFGS).
!>
!> The model Hessian is built in redundant internal coordinates and contracted
!> to Cartesian by the existing xtb implementation
!> (`TModelHessian%compute`, `src/model_hessian/type.f90:41`); the provider
!> only supplies the geometry and the atomic numbers.  The construction is
!> algebraic and deterministic, so it costs no energy evaluation.
subroutine model_hessian_compute(self, env, mol, chk, calc, hess, step, odlr, final_err)
   !> Provider instance
   class(landscape_model_hessian_provider), intent(inout) :: self
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure
   type(TMolecule), intent(inout) :: mol
   !> Wavefunction data (unused, the model Hessian is geometry only)
   type(TRestart), intent(inout) :: chk
   !> Single point calculator (unused)
   class(TCalculator), intent(inout) :: calc
   !> Dense Hessian (3N, 3N), symmetrised on return
   real(wp), intent(inout) :: hess(:, :)
   !> Finite difference step (bohr), unused
   real(wp), intent(in) :: step
   !> Use the ODLR approximation, unused
   logical, intent(in) :: odlr
   !> Residual error estimate (0, the model is analytic in the geometry)
   real(wp), intent(out) :: final_err

   class(TModelHessian), allocatable :: model
   type(modhess_setvar) :: modh
   integer :: n

   final_err = 0.0_wp
   self%n_eval = 0
   n = mol%n
   if (n <= 0) return

   ! Defaults of `set%mhset` (`src/setparam.f90:267-278`): the 1995 Lindh
   ! parameterisation with the standard stretch/bend/torsion force constants.
   modh%kr = 0.4_wp
   modh%kf = 0.13_wp
   modh%kt = 0.0075_wp
   modh%ko = 0.0_wp
   modh%kd = 0.0_wp
   modh%kq = 0.0_wp
   modh%rcut = 70.0_wp
   modh%s6 = 20.0_wp

   allocate(model, source=newLindhD2ModelHessian(modh))
   hess = 0.0_wp
   call model%compute(env, mol%xyz, n, hess, mol%at)
   hess = 0.5_wp*(hess + transpose(hess))

end subroutine model_hessian_compute

!> Hessian of the attached file (`hessian_source = 'file'`).
!>
!> The dimension of the matrix has to match the molecule, because the
!> TURBOMOLE layout stores no dimension: a file that describes a different
!> molecule is only detectable through its number of values, which the reader
!> checks against `3n` and reports when it differs.  The attached calculator
!> and wavefunction are untouched, the provider never evaluates them.
subroutine file_hessian_compute(self, env, mol, chk, calc, hess, step, odlr, final_err)
   character(len=*), parameter :: source = 'landscape_file_hessian_compute'
   !> Provider instance
   class(landscape_hessian_file), intent(inout) :: self
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure (only the atom count is used)
   type(TMolecule), intent(inout) :: mol
   !> Wavefunction data (unused)
   type(TRestart), intent(inout) :: chk
   !> Single point calculator (unused)
   class(TCalculator), intent(inout) :: calc
   !> Dense Hessian (3N, 3N)
   real(wp), intent(inout) :: hess(:, :)
   !> Finite difference step (bohr), unused
   real(wp), intent(in) :: step
   !> Use the ODLR approximation, unused
   logical, intent(in) :: odlr
   !> Residual error estimate (0, the file is exact data)
   real(wp), intent(out) :: final_err

   logical :: ok

   final_err = 0.0_wp
   self%n_eval = 0
   if (odlr) then
      call env%error("the Hessian file provider returns the stored matrix and "//&
         & "does not support the ODLR approximation", source)
      return
   end if
   if (.not. allocated(self%filename)) then
      call env%error("the Hessian file provider has no file name", source)
      return
   end if
   if (size(hess, 1) /= 3*mol%n .or. size(hess, 2) /= 3*mol%n) then
      call env%error("the Hessian of the file provider must be a (3N, 3N) "//&
         & "matrix; the caller passed a different shape", source)
      return
   end if

   call landscape_read_hessian_file(env, self%filename, hess, ok)

end subroutine file_hessian_compute

!> Read a TURBOMOLE format Hessian from `filename`.
!>
!> The layout is the one of the `hessian` file written by xtb itself
!> (`wrhess`/`rdhess` in `src/freq/io.f90`): the record `$hessian` (a trailing
!> ` (projected)` is allowed) followed by the full square matrix, row by row,
!> five columns per record in the format `(5x,5f15.10)`; the optional closing
!> record `$end` is skipped.  The values are assigned to `hess(j, i)` for the
!> `j`-th entry of row `i`, so an asymmetric matrix (which the symmetric
!> writer never produces) is detected instead of being averaged silently.
!>
!> Failure modes, all reported with `env%error` and `ok = .false.`:
!>
!>  * the file does not exist or cannot be opened,
!>  * the first record is not a `$hessian` header,
!>  * a record holds a token that is not a number,
!>  * the number of values is not `(3n)^2`, i.e. the file describes a
!>    different molecule,
!>  * the matrix is not symmetric, or holds a non-finite value.
!>
!> `hess` is left untouched on failure.
subroutine landscape_read_hessian_file(env, filename, hess, ok)
   character(len=*), parameter :: source = 'landscape_read_hessian_file'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Name of the Hessian file
   character(len=*), intent(in) :: filename
   !> Dense Hessian (3n, 3n) in hartree/bohr^2
   real(wp), intent(out) :: hess(:, :)
   !> True if the matrix was read completely
   logical, intent(out) :: ok

   real(wp), allocatable :: work(:, :)
   real(wp) :: values(hessian_max_record), asym, scale
   character(len=1024) :: line, token
   character(len=128) :: iomsg
   integer :: unit, n3, nvals, irec, ios, itoken, ntoken
   integer :: istart, iend, irow, icol, i
   logical :: terminator, opened

   ok = .false.
   if (size(hess, 1) /= size(hess, 2) .or. size(hess, 1) <= 0) then
      call env%error("the Hessian target matrix must be a non-empty square "//&
         & "matrix", source)
      return
   end if
   n3 = size(hess, 1)

   ! `open_hessian_for_read` resolves the xtb file namespace (the same prefix
   ! `main_freq` and `wrhess` add to the files they write) and falls back to
   ! the plain name, so a Hessian written by a namespaced run is found both as
   ! `<namespace>.<name>` and as `<name>`.
   call open_hessian_for_read(unit, trim(filename), opened)
   if (.not. opened) then
      call env%error("cannot open Hessian file '"//trim(filename)//"'", source)
      return
   end if

   line = ''
   read(unit, '(a)', iostat=ios, iomsg=iomsg) line
   if (ios /= 0) then
      close(unit)
      call env%error("the Hessian file '"//trim(filename)//"' is empty or "//&
         & "unreadable ("//trim(iomsg)//")", source)
      return
   end if
   if (index(line, hessian_header) == 0) then
      close(unit)
      call env%error("the Hessian file '"//trim(filename)//"' is not a "//&
         & "TURBOMOLE Hessian: the first record does not hold '"//&
         & hessian_header//"'", source)
      return
   end if

   allocate(work(n3, n3), source=0.0_wp)
   nvals = 0
   irec = 0
   terminator = .false.
   do
      read(unit, '(a)', iostat=ios, iomsg=iomsg) line
      if (ios < 0) exit
      if (ios > 0) then
         close(unit)
         call env%error("cannot read record "//trim(landscape_int_to_string(irec + 2))//&
            & " of the Hessian file '"//trim(filename)//"' ("//trim(iomsg)//")", &
            & source)
         return
      end if
      irec = irec + 1
      if (index(line, '$end') > 0) then
         terminator = .true.
         exit
      end if
      ! Split the record into blank/comma separated tokens and parse each of
      ! them individually, so that a corrupt token is reported with its
      ! position instead of silently ending the list directed read.
      ntoken = 0
      i = 1
      do while (i <= len(line))
         if (is_separator(line(i:i))) then
            i = i + 1
            cycle
         end if
         istart = i
         do while (i <= len(line))
            if (is_separator(line(i:i))) exit
            i = i + 1
         end do
         iend = i - 1
         ntoken = ntoken + 1
         if (ntoken > hessian_max_record) then
            close(unit)
            call env%error("record "//trim(landscape_int_to_string(irec + 1))//" of the "//&
               & "Hessian file '"//trim(filename)//"' holds more than "//&
               & trim(landscape_int_to_string(hessian_max_record))//" values", source)
            return
         end if
         token = line(istart:iend)
         read(token, *, iostat=ios) values(ntoken)
         if (ios /= 0) then
            close(unit)
            call env%error("record "//trim(landscape_int_to_string(irec + 1))//" of the "//&
               & "Hessian file '"//trim(filename)//"' holds the token '"//&
               & trim(token)//"', which is not a number", source)
            return
         end if
      end do
      do itoken = 1, ntoken
         nvals = nvals + 1
         if (nvals > n3*n3) cycle
         irow = (nvals - 1)/n3 + 1
         icol = nvals - (irow - 1)*n3
         ! `rdhess` reads row `irow` into the column `icol` of the matrix, so
         ! the same assignment is used here
         work(icol, irow) = values(itoken)
      end do
   end do
   close(unit)

   if (terminator .and. nvals < n3*n3) then
      call env%error("the Hessian file '"//trim(filename)//"' ends after "//&
         & trim(landscape_int_to_string(nvals))//" values, but the matrix of 3N = "//&
         & trim(landscape_int_to_string(n3))//" needs "//trim(landscape_int_to_string(n3*n3)), &
         & source)
      return
   end if
   if (nvals /= n3*n3) then
      call env%error("the Hessian file '"//trim(filename)//"' holds "//&
         & trim(landscape_int_to_string(nvals))//" values, but the matrix of 3N = "//&
         & trim(landscape_int_to_string(n3))//" needs "//trim(landscape_int_to_string(n3*n3))//&
         & "; the file describes a different molecule", source)
      return
   end if

   do irow = 1, n3
      do icol = 1, n3
         if (.not. landscape_value_present(work(icol, irow))) then
            call env%error("the Hessian file '"//trim(filename)//"' holds a "//&
               & "non-finite matrix element at row "//trim(landscape_int_to_string(irow))//&
               & ", column "//trim(landscape_int_to_string(icol)), source)
            return
         end if
      end do
   end do

   asym = maxval(abs(work - transpose(work)))
   scale = max(1.0_wp, maxval(abs(work)))
   if (asym > hessian_symmetry_tol*scale) then
      call env%error("the Hessian read from '"//trim(filename)//"' is not "//&
         & "symmetric: max |h(i,j) - h(j,i)| = "//&
         & trim(landscape_real_to_string(asym, '(es16.8)'))//" Eh/bohr^2", source)
      return
   end if

   hess = work
   ok = .true.

end subroutine landscape_read_hessian_file

!> Write a dense Hessian in the TURBOMOLE format of xtb.
!>
!> The file holds the **complete** square matrix (row `i` carries the `n3`
!> elements `h(i, 1:n3)`), not the lower triangle alone: the packed array is
!> filled with `packed(lin(i,j)) = hess(i,j)` for `j <= i` and `lin(j,i)`
!> mirrors it for `j > i`, and the records are written five columns at a time -
!> exactly like `wrhess` (`src/freq/io.f90:53`) writes the `hessian` file of
!> `--hess`.  The writer is not `wrhess` itself because that routine cannot
!> report an unwritable target (`persistentEnv%io%writeFile` leaves the unit at
!> -1 and the following `write` terminates the program); the format is kept
!> identical, which the unit test `hessian-file-roundtrip` verifies by
!> comparing the output of both writers byte by byte.
subroutine landscape_write_hessian_file(env, filename, hess, ok)
   character(len=*), parameter :: source = 'landscape_write_hessian_file'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Name of the Hessian file
   character(len=*), intent(in) :: filename
   !> Dense Hessian (3n, 3n) in hartree/bohr^2
   real(wp), intent(in) :: hess(:, :)
   !> True if the matrix was written
   logical, intent(out) :: ok

   real(wp), allocatable :: packed(:)
   character(len=128) :: iomsg
   integer :: unit, n3, i, j, mincol, maxcol, ios
   logical :: opened

   ok = .false.
   if (size(hess, 1) /= size(hess, 2) .or. size(hess, 1) <= 0) then
      call env%error("the Hessian to write must be a non-empty square matrix", &
         & source)
      return
   end if
   n3 = size(hess, 1)

   call open_hessian_for_write(unit, trim(filename), opened)
   if (.not. opened) then
      call env%error("cannot write Hessian file '"//trim(filename)//"'", source)
      return
   end if

   allocate(packed(n3*(n3 + 1)/2), source=0.0_wp)
   do i = 1, n3
      do j = 1, i
         packed(lin(i, j)) = hess(i, j)
      end do
   end do

   write(unit, '(a)') hessian_header
   do i = 1, n3
      maxcol = 0
      do while (maxcol < n3)
         mincol = maxcol + 1
         maxcol = min(maxcol + 5, n3)
         write(unit, hessian_format, iostat=ios, iomsg=iomsg) &
            & (packed(lin(i, j)), j=mincol, maxcol)
         if (ios /= 0) then
            close(unit)
            call env%error("cannot write record "//trim(landscape_int_to_string(i))//&
               & " of the Hessian file '"//trim(filename)//"' ("//&
               & trim(iomsg)//")", source)
            return
         end if
      end do
   end do
   close(unit)
   ok = .true.

end subroutine landscape_write_hessian_file

!> Open a Hessian file for reading, namespaced name first.
!>
!> Mirrors `persistentEnv%io%readFile` (`src/type/iohandler.f90:202`): the
!> namespaced name is tried first, then the name as given, so a run with
!> `--namespace` finds both the file it wrote itself and the plain default
!> name.
subroutine open_hessian_for_read(unit, filename, opened)
   !> Unit of the opened file, -1 if it could not be opened
   integer, intent(out) :: unit
   !> Name of the Hessian file
   character(len=*), intent(in) :: filename
   !> True if the file was opened
   logical, intent(out) :: opened

   character(len=:), allocatable :: name
   integer :: ios

   unit = -1
   opened = .false.
   name = get_namespace(filename)
   open(newunit=unit, file=name, status='old', action='read', iostat=ios)
   if (ios == 0) then
      opened = .true.
      return
   end if
   unit = -1
   if (name == filename) return
   open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
   if (ios == 0) opened = .true.
   if (.not. opened) unit = -1

end subroutine open_hessian_for_read

!> Open a Hessian file for writing under the namespace of the run.
subroutine open_hessian_for_write(unit, filename, opened)
   !> Unit of the opened file, -1 if it could not be opened
   integer, intent(out) :: unit
   !> Name of the Hessian file
   character(len=*), intent(in) :: filename
   !> True if the file was opened
   logical, intent(out) :: opened

   character(len=:), allocatable :: name
   integer :: ios

   unit = -1
   name = get_namespace(filename)
   open(newunit=unit, file=name, action='write', iostat=ios)
   opened = ios == 0
   if (.not. opened) unit = -1

end subroutine open_hessian_for_write

!> True if the character separates the values of a Hessian record.
!>
!> Blank and comma are the separators of the format; tab and the carriage
!> return of a CRLF line ending are accepted as well, so a file produced on a
!> different platform is read instead of being rejected as corrupt.
pure function is_separator(char) result(separator)
   !> Character to test
   character(len=1), intent(in) :: char
   !> True for blank, comma, tab and carriage return
   logical :: separator

   separator = char == ' ' .or. char == ',' .or. char == achar(9) .or. &
      & char == achar(13)

end function is_separator

end module xtb_landscape_hessian
