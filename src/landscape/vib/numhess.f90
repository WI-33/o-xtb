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

!> Hessian production of the `vib` task (PLAN section 9.4).
!>
!> The task layer turns the `landscape_vib_settings` of a frequency analysis
!> into the three quantities every consumer of this task needs: the reference
!> energy and gradient of the geometry (a plain single point through the engine
!> cache) and the dense, symmetric Hessian of the configured provider
!> (`numerical`, `model` or `file`).  The provider itself is **not** reachable
!> from here by any other route: the engine is the only place that calls it
!> (PLAN section 4.2), so this module owns the engine instance and hands the
!> matrix on.
!>
!> Two properties are deliberate:
!>
!>  * the finite difference step (`engine%hess_step`) and the SCC accuracy of
!>    the displaced single points (`engine%hess_acc`) are taken from the
!>    settings object.  The xtb globals `set%step_hess`/`set%accu_hess` are
!>    **never** read: their `$hess` control block is silently ignored on the
!>    non-`--copy` path (`docs/_landscape-recon/07-*.md`).
!>  * the ODLR approximation (`src/o1numhess.f90:71`) is not implemented by any
!>    provider of this suite.  It is rejected here, before a single SCF, for
!>    every source - including `model` and `file`, whose `compute` would
!>    otherwise silently ignore the flag and return a full dense matrix.
!>
!> Units: coordinates in bohr, energies in hartree, gradients in hartree/bohr
!> and Hessian elements in hartree/bohr^2.
module xtb_landscape_vib_numhess
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   use xtb_type_restart, only : TRestart
   use xtb_type_calculator, only : TCalculator
   use xtb_landscape_settings, only : landscape_vib_settings
   use xtb_landscape_engine, only : landscape_engine, new_landscape_engine, &
      & landscape_hessian_provider, new_landscape_hessian_provider
   implicit none(type, external)
   private

   public :: landscape_vib_hessian

contains

!> Hessian, energy and gradient of `mol` for the `vib` task.
!>
!> On success `hess` holds the dense, symmetric Hessian of the configured
!> provider.  The energy and the gradient come from the engine evaluation of
!> the same geometry, which is also the reference state of a numerical
!> differentiation; they are returned because the frequency task needs the
!> energy for the thermochemistry and the gradient norm for its report.
subroutine landscape_vib_hessian(env, mol, chk, calc, settings, hess, energy, &
      & gradient, final_err, fail, n_eval, etemp)
   character(len=*), parameter :: source = 'landscape_vib_hessian'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure (geometry is not modified)
   type(TMolecule), intent(inout) :: mol
   !> Wavefunction data
   type(TRestart), intent(inout) :: chk
   !> Single point calculator
   class(TCalculator), intent(in) :: calc
   !> Frequency settings
   type(landscape_vib_settings), intent(in) :: settings
   !> Dense Hessian (3N, 3N) in hartree/bohr^2
   real(wp), intent(out) :: hess(:, :)
   !> Electronic energy of the geometry (hartree)
   real(wp), intent(out) :: energy
   !> Cartesian gradient of the geometry (hartree/bohr)
   real(wp), intent(out) :: gradient(:, :)
   !> Residual error estimate of the Hessian provider (0 for exact data)
   real(wp), intent(out) :: final_err
   !> Failure flag
   logical, intent(out) :: fail
   !> Number of energy/gradient evaluations (6N+1 for the numerical Hessian,
   !> 1 for the model and file sources)
   integer, intent(out), optional :: n_eval
   !> Effective electronic temperature read back from the calculator (kelvin);
   !> zero when the engine never touched the calculator (the provider failed
   !> before the first evaluation).  Recorded for the same reason as in the
   !> saddle and IRC tasks (M4 hard constraint 1): the value in effect, not
   !> the request.
   real(wp), intent(out), optional :: etemp

   type(landscape_engine) :: engine
   class(landscape_hessian_provider), allocatable :: hp
   logical :: hess_fail

   fail = .false.
   energy = 0.0_wp
   gradient = 0.0_wp
   hess = 0.0_wp
   if (present(etemp)) etemp = 0.0_wp
   final_err = 0.0_wp
   if (present(n_eval)) n_eval = 0

   ! The ODLR approximation of `src/o1numhess.f90` is not available in this
   ! suite.  Refusing it here covers every source; the numerical provider
   ! refuses it a second time as its own contract.
   if (settings%engine%odlr) then
      call env%error("the ODLR approximation of the numerical Hessian "//&
         & "(src/o1numhess.f90) is not implemented in the landscape suite", &
         & source)
      fail = .true.
      return
   end if

   call new_landscape_engine(engine, settings%engine)
   ! The factory reports an unknown source (and a missing file name) and
   ! returns without a provider; the caller must abort instead of falling back
   ! to a default backend.
   call new_landscape_hessian_provider(env, settings%hessian_source, hp, &
      & sccacc=settings%engine%hess_acc, filename=settings%hessian_file, &
      & printlevel=settings%printlevel, jobs=settings%engine%hess_jobs)
   if (.not. allocated(hp)) then
      fail = .true.
      call engine%destroy()
      return
   end if
   call engine%set_hessian_provider(hp)
   call engine%set_calculator(calc)

   call engine%energy_gradient(env, mol, chk, energy, gradient, fail)
   if (fail) then
      call env%error("landscape frequency analysis: the single point of the "//&
         & "reference geometry failed", source)
      call engine%destroy()
      return
   end if

   call engine%hessian(env, mol, chk, hess, hess_fail, step=settings%engine%hess_step, &
      & final_err=final_err)
   if (hess_fail) then
      call env%error("landscape frequency analysis: the "//&
         & trim(settings%hessian_source)//" Hessian provider failed", source)
      fail = .true.
      call engine%destroy()
      return
   end if

   ! Record the electronic temperature that is actually in effect.  The engine
   ! reads it back from the calculator (`engine_report`), so a `--etemp` value
   ! that never reached the engine (missing CLI sync) or a dead control field
   ! becomes visible on stdout instead of silently running at 300 K.  PLAN
   ! section 9.5 makes this record a gate for the M4 transition-state and IRC
   ! entries.  The same call also prints the cost counters of the run, and the
   ! unused notice follows it: it is derived from the builds the engine
   ! actually performed, so it agrees with the `Hessians` counter by
   ! construction and `--lvib-hess-step` can no longer be silently ignored for
   ! `--lvib-hessian model|file` (`37-stage2b-fix2.md` section T2).
   if (settings%printlevel > 0) then
      call engine%report(env%unit)
      call engine%report_unused_hess_step(env%unit)
   end if

   ! The engine adds the displaced points of the provider to its own counter,
   ! so the reported cost is the true 6N+1 of the numerical Hessian.
   if (present(n_eval)) n_eval = engine%n_eval
   ! The effective electronic temperature is read back from the calculator
   ! (M4 hard constraint 1); the frequency task stores it in its result so it
   ! reaches `landscape_<task>.json` as well as stdout.
   if (present(etemp)) etemp = engine%effective_etemp()
   call engine%destroy()

end subroutine landscape_vib_hessian

end module xtb_landscape_vib_numhess
