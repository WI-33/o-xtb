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

!> Transition state search task of the landscape suite (PLAN section 9.5).
!>
!> This is the task level entry point of the `--lts` run type.  It builds the
!> engine, wraps it into the injectable `landscape_energy_model`
!> (`xtb_landscape_tasks_opt::landscape_engine_model`) and runs one of the
!> searches:
!>
!>  * `xtb_landscape_saddle_prfo::landscape_saddle_optimize` for P-RFO and
!>    B-P-RFO, which follow an eigenvector of the (numerical) Hessian;
!>  * `xtb_landscape_saddle_dimer::landscape_dimer_optimize` for the dimer,
!>    which uses only first derivatives and therefore never asks the engine for
!>    a Hessian during the search;
!>  * `xtb_landscape_path_neb::landscape_band_optimize` for the band methods
!>    (`neb`, `cineb`), `xtb_landscape_path_autoneb::landscape_autoneb_optimize`
!>    for the incremental path building (`autoneb`) and
!>    `xtb_landscape_path_gsm::landscape_string_optimize` for the growing
!>    string (`string`, `gsm`).  The three chain-of-states methods share the
!>    input contract (reactant positional, product via `--lts-product`), the
!>    `path.xyz` format and the frequency confirmation of the candidate.
!>
!> After the search the stationary point is confirmed with a frequency analysis
!> (`xtb_landscape_vib_modes::landscape_vib_analysis`, the analysis of the M3
!> `--lvib` task) and exactly one imaginary mode is required.
!>
!> Output:
!>
!> | file                | content |
!> |---------------------|---------|
!> | `xtbts.xyz`         | converged geometry (energy and gradient norm in the comment) |
!> | `xtbts.log`         | method, frozen `etemp`, step table, frequency confirmation |
!> | `path.xyz`          | the band of `neb` / `cineb` / `string` / `gsm`, one frame per image |
!> | `path_autoneb.xyz`  | the path of AutoNEB (also its restart file) |
!> | `path_ts_list.xyz`  | every saddle candidate found, with its energy and imaginary mode count |
!>
!> M4 hard constraint 1 (PLAN section 9.5): the electronic temperature is
!> frozen for the whole run.  `settings%etemp` is the single source; it is
!> pushed into the engine control (and from there into the xTB calculator) and
!> the effective value is printed both from the settings and, read back from
!> the calculator, by `engine%report`.
!>
!> M4 hard constraint 2: the confirmation reports the number and value of the
!> imaginary modes, so a search that ended on a minimum (0) or on a higher
!> order saddle (>1) is reported and **fails** the task (H2 of
!> `29-m4-independent-review.md`): a non-saddle aborts with a non-zero exit
!> status after writing the geometry and the log for debugging.  A mode whose
!> |frequency| lies inside `0 +/- 3 cm^-1` is re-decided with a tighter finite
!> difference (`step = 2e-3` bohr, `sccacc = 0.1`) before its sign is trusted;
!> both spectra are printed (H1).  The dimer search applies the same rule to
!> its own curvature: a point whose gradient meets the criteria but whose
!> curvature along the dimer axis is not negative is not accepted
!> (`xtb_landscape_saddle_dimer`).
!>
!> Units: coordinates in bohr, energies in hartree, gradients in hartree/bohr,
!> frequencies in cm^-1, `etemp` in kelvin.
module xtb_landscape_tasks_ts
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : autokcal, autoaa
   use xtb_mctc_filetypes, only : fileType
   use xtb_io_writer, only : writeMolecule
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule, init
   use xtb_type_restart, only : TRestart
   use xtb_type_calculator, only : TCalculator
   use xtb_landscape_settings, only : landscape_ts_settings, landscape_vib_settings, &
      & landscape_opt_settings, landscape_engine_control, landscape_restart_name, &
      & landscape_ts_method_name, p_landscape_ts_prfo, p_landscape_ts_bprfo, &
      & p_landscape_ts_dimer, p_landscape_ts_cineb, p_landscape_ts_autoneb, &
      & landscape_ts_method_is_band, landscape_ts_method_is_string, landscape_log_failed, &
      & landscape_lowercase, landscape_int_to_string
   use xtb_landscape_results, only : landscape_result, landscape_vib_result
   use xtb_landscape_coords, only : geometry_pack, geometry_unpack, angstrom_to_bohr
   use xtb_landscape_engine, only : landscape_engine, new_landscape_engine, &
      & landscape_hessian_provider, new_landscape_hessian_provider
   use xtb_landscape_tasks_opt, only : landscape_engine_model
   use xtb_landscape_optimize, only : landscape_optimize_model
   use xtb_landscape_constraints, only : landscape_frozen_mask
   use xtb_landscape_saddle_prfo, only : landscape_saddle_optimize
   use xtb_landscape_saddle_dimer, only : landscape_dimer_optimize
   use xtb_landscape_vib_modes, only : landscape_vib_analysis
   use xtb_landscape_path_neb, only : landscape_band_optimize, landscape_band_type
   use xtb_landscape_path_autoneb, only : landscape_autoneb_optimize, &
      & landscape_autoneb_type
   use xtb_landscape_path_gsm, only : landscape_string_optimize, landscape_string_type
   use xtb_io_reader, only : readMolecule
   use xtb_mctc_filetypes, only : fileType, getFileType
   use xtb_landscape_json, only : landscape_json_document, &
      & landscape_json_document_init, landscape_json_add_result, &
      & landscape_json_add_core, &
      & landscape_json_add_counter, landscape_json_add_value, &
      & landscape_json_add_flag, landscape_json_write, landscape_json_status_ok, &
      & landscape_json_status_incomplete, landscape_json_status_error
   use xtb_setparam, only : get_namespace
   implicit none(type, external)
   private

   public :: landscape_ts
   public :: landscape_ts_curvature_band_hit
   public :: landscape_ts_curvature_band_cm1
   public :: landscape_ts_can_tighten
   public :: confirm_saddle   ! exposed for the curvature-band wiring test
   public :: landscape_ts_count_frames
   public :: landscape_ts_read_path

   !> Half width of the curvature band around zero (cm^-1) inside which the
   !> sign of a mode is limited by the resolution of the numerical Hessian of
   !> the confirmation.  Requirement 2 of the M4 hard constraints: a mode whose
   !> |frequency| falls inside `0 +/- 3 cm^-1` is recomputed with a tighter
   !> finite difference before its sign is trusted, because the sign is not
   !> meaningful at the default resolution there.
   real(wp), parameter :: landscape_ts_curvature_band_cm1 = 3.0_wp
   !> Tightened finite difference step (bohr) of the recomputation
   real(wp), parameter :: ts_tight_hess_step = 2.0e-3_wp
   !> Tightened SCC accuracy of the recomputation
   real(wp), parameter :: ts_tight_hess_acc = 0.1_wp

contains

!> Task level entry point: locate a first order saddle point.
!>
!> On return `mol%xyz` holds the converged geometry and `result` the search and
!> confirmation summary.  Every fatal condition (unknown method, a failed
!> search, an unconverged SCC) is reported through `env%error`, which the
!> caller turns into a non-zero exit status.
subroutine landscape_ts(env, mol, chk, calc, settings, result)
   character(len=*), parameter :: source = 'landscape_ts'
   !> Computational environment
   type(TEnvironment), intent(inout), target :: env
   !> Molecular structure, updated to the transition state
   type(TMolecule), intent(inout), target :: mol
   !> Wavefunction data
   type(TRestart), intent(inout), target :: chk
   !> Single point calculator
   class(TCalculator), intent(in) :: calc
   !> Transition state settings
   type(landscape_ts_settings), intent(in) :: settings
   !> Search result
   type(landscape_result), intent(inout) :: result

   type(landscape_engine), target :: engine
   type(landscape_engine_model) :: model
   class(landscape_hessian_provider), allocatable :: hp
   type(landscape_engine_control) :: ctrl
   real(wp), allocatable :: x0(:), xfinal(:), traj(:, :)
   integer :: n_scf_fail

   if (settings%printlevel > 0) then
      write(env%unit, '(1x,a,a,a)') "landscape: transition state search (", &
         & landscape_ts_method_name(settings%method), ")"
      if (.not. landscape_ts_method_is_band(settings%method)) then
         write(env%unit, '(1x,a,es15.6,a)') "initial trust radius:     ", &
            & settings%trust, " bohr"
         write(env%unit, '(1x,a,es15.6,a)') "maximum trust radius:     ", &
            & settings%trust_max, " bohr"
      end if
      write(env%unit, '(1x,a,i0)') "maximum number of steps:  ", settings%max_steps
      ! R5 of CLI redesign section 3.5: `e_tol` follows the convergence level
      ! (`set_conv_level`) and is part of the convergence test of the search
      ! itself.  The value comes from the settings object, not from the flag
      ! text.
      write(env%unit, '(1x,a,es15.6,a)') "energy tolerance:         ", &
         & settings%e_tol, " Eh"
      write(env%unit, '(1x,a,a)') "restart policy:           ", &
         & landscape_restart_name(settings%engine%restart_mode)
      ! The dimer search itself never builds a Hessian; the provider is still
      ! created because the frequency confirmation of the located point uses
      ! it, so the source is reported for every method.
      write(env%unit, '(1x,a,a)') "Hessian source:           ", &
         & trim(settings%hessian_source)
      ! CLI redesign sections 7.3 and 9.2: every flag that carries a unit must
      ! echo its effective value and unit from the settings object the task
      ! hands to the engine, not from the command line text.  Whether the
      ! step reaches a build is **not** decidable here (the path searches
      ! build their confirmation and their candidate Hessians at different
      ! points of the run, and a run that exhausts its budget never gets
      ! there), so the unused notice is not guessed from the method name: it
      ! is derived from the builds the engine actually performed
      ! (`engine%report_unused_hess_step`, printed with the cost counters at
      ! the end of the run) - see `37-stage2b-fix2.md` section T1.
      write(env%unit, '(1x,a,es15.6,a)') "Hessian step:             ", &
         & settings%engine%hess_step, " bohr"
      write(env%unit, '(1x,a,es15.6)') "Hessian SCC accuracy:     ", &
         & settings%engine%hess_acc
      if (settings%method /= p_landscape_ts_dimer) then
         ! MEDIUM-1 of `31-m4-review-round2.md`: the escape from a minimum is
         ! reported with the run so the displacement is visible in the log.
         write(env%unit, '(1x,a,f12.6,a)') "escape from a minimum:      ", &
            & settings%escape_step, " bohr along the lowest internal mode"
      end if
      ! D7 of the CLI redesign section 8.3: the escape walk moves the **input
      ! geometry** before the search, so the active default is echoed.  The
      ! notice is suppressed when `--lts-escape-step` named the length and is
      ! only printed for the two methods that use the escape walk
      ! (`xtb_landscape_saddle_prfo::escape_from_minimum` returns for every
      ! other method).
      if ((settings%method == p_landscape_ts_prfo .or. &
            & settings%method == p_landscape_ts_bprfo) .and. &
            & .not. settings%escape_step_was_set) then
         write(env%unit, '(a,f4.2,a)') "[NOTICE] --lts-escape-step default ", &
            & settings%escape_step, " bohr moves the input geometry off a "//&
            & "minimum"
         write(env%unit, '(a)') &
            & "         pass --lts-escape-step 0 to disable the escape walk"
      end if
      if (settings%method == p_landscape_ts_dimer) then
         write(env%unit, '(1x,a,es13.6,a)') "dimer separation:         ", &
            & settings%delta, " bohr"
         write(env%unit, '(1x,a,es13.6,a)') "dimer rotation probe:     ", &
            & settings%rot_dtheta, " rad"
         write(env%unit, '(1x,a,i0)') "dimer rotation iterations: ", &
            & settings%rot_max_iter
         write(env%unit, '(1x,a,a,a,i0)') "dimer initial axis:       ", &
            & trim(settings%n_init), " (seed ", settings%seed, ")"
      end if
      ! M4 hard constraint 1: the frozen electronic temperature is part of the
      ! result, not an implicit default.
      write(env%unit, '(1x,a,f12.3,a)') "electronic temperature:   ", &
         & settings%etemp, " K (frozen for the whole run)"
   end if

   ! `--lts-product` is a global option that a script may set for every run;
   ! the single ended methods (P-RFO, B-PRFO, dimer) ignore it, which must
   ! not happen silently.
   if (.not. landscape_ts_method_is_band(settings%method)) then
      if (allocated(settings%path%product_file)) then
         call env%warning("--lts-product is only used by the path methods "//&
            & "(`--lts-method neb|cineb|autoneb|string|gsm`) and is ignored "//&
            & "by '"//trim(landscape_ts_method_name(settings%method))//"'", &
            & 'landscape_ts')
      end if
   end if

   ! `$fix` is supported by the single ended searches (P-RFO, B-PRFO, dimer:
   ! their accepted step is projected onto the subspace `$fix` leaves free,
   ! `saddle/prfo.f90`, `saddle/dimer.f90`), but **not** by the band searches.
   ! A band image lies between the two endpoints and is rebuilt by
   ! interpolation, by the IDPP initial path and by the redistribution after
   ! every step (`path/idpp.f90`, `path/string.f90`, `path/gsm.f90`), so
   ! "which coordinate is fixed" has no defined meaning when the two endpoints
   ! disagree on it - `$fix atoms: 1` with a product whose atom 1 is elsewhere
   ! would move atom 1 of every interior image.  The run is refused instead of
   ! moved silently (docs/_simplify-recon/04-o2-and-silent-defects-fix.md
   ! section 5.2; PLAN section 6.9 open item).
   if (landscape_ts_method_is_band(settings%method)) then
      if (any(landscape_frozen_mask(mol%n))) then
         call env%error("the band searches (--lts-method neb|cineb|autoneb|"//&
            & "string|gsm) do not support $fix: their interior images are "//&
            & "interpolated between the two endpoints.  Use --lts-method "//&
            & "prfo, bprfo or dimer, or remove the $fix block", source)
         result%message = "the band searches do not support $fix"
         return
      end if
   end if

   ! The engine owns a copy of the calculator.  The frozen etemp is the one of
   ! the task settings, copied into the control object the engine pushes into
   ! the calculator, so a CLI `--etemp` reaches the SCF of every evaluation.
   ctrl = settings%engine
   ctrl%etemp = settings%etemp
   call new_landscape_engine(engine, ctrl)
   call new_landscape_hessian_provider(env, settings%hessian_source, hp, &
      & sccacc=settings%engine%hess_acc, filename=settings%hessian_file, &
      & printlevel=settings%printlevel)
   if (.not. allocated(hp)) then
      result%message = "transition state search failed"
      call engine%destroy()
      return
   end if
   call engine%set_hessian_provider(hp)
   call engine%set_calculator(calc)
   model%engine => engine
   model%env => env
   model%mol => mol
   model%chk => chk

   ! The chain-of-states searches (NEB / CI-NEB) are two more methods of the
   ! `--lts` run type: they reuse the engine, the Hessian provider and the
   ! frequency confirmation, but replace the single geometry search by a band
   ! of images between two fixed endpoints.
   if (landscape_ts_method_is_band(settings%method)) then
      call landscape_ts_band(env, mol, chk, engine, settings, result)
      if (settings%ljson) call write_json(env, mol, settings, result, engine%n_scf_fail)
      call engine%report(env%unit)
      call engine%report_unused_hess_step(env%unit)
      call engine%destroy()
      return
   end if

   x0 = geometry_pack(mol%xyz)
   select case (settings%method)
   case (p_landscape_ts_dimer)
      call landscape_dimer_optimize(env, model, settings, x0, result, xfinal, traj)
   case default
      call landscape_saddle_optimize(env, model, settings, x0, result, xfinal, traj)
   end select
   if (.not. allocated(xfinal)) then
      call env%error("the transition state search did not return a geometry", source)
      call engine%destroy()
      return
   end if

   call geometry_unpack(xfinal, mol%xyz)
   ! M1 (review 29): record the *effective* electronic temperature read back
   ! from the calculator, not the setting echoed from the request.  A value
   ! that never reached the calculator would otherwise be invisible in the
   ! result object (the M3 MEDIUM-B "dead field" defect).
   result%etemp = engine%effective_etemp()
   n_scf_fail = engine%n_scf_fail

   ! Frequency confirmation of the stationary point.  The Hessian is the one of
   ! the configured provider at the converged geometry; the analysis is the M3
   ! implementation, so the mode count and the 3N-6 / 3N-5 handling are
   ! identical to `--lvib`.
   if (settings%verify_ts) then
      call confirm_saddle(env, mol, chk, engine, settings, result)
   end if

   ! The evaluation counter is sampled *after* the confirmation: it is the
   ! cost of the complete run and equals the `landscape engine: energy
   ! evaluations` line printed below for the same run (M7 fix round, report
   ! 50 finding M-1 - the sampling point used to be before the confirmation,
   ! so the two numbers of one and the same stdout disagreed by the 6N finite
   ! difference evaluations of the confirmation Hessian, and the JSON
   ! counter was the search-only number).  `--lts-no-verify` therefore
   ! lowers the counter by exactly one 6N Hessian; the CLI check
   ! `landscape_e2e_check.py ts-counters` asserts that difference.
   result%n_cache_miss = engine%n_eval

   if (settings%printlevel > 0) then
      call result%write_log(env%unit)
      call result%conv%write_log(env%unit)
   end if

   if (settings%write_files) then
      call write_ts_files(env, mol, settings, result)
   end if

   call engine%report(env%unit)
   call engine%report_unused_hess_step(env%unit)
   call engine%destroy()

   if (settings%ljson) call write_json(env, mol, settings, result, n_scf_fail)

   if (n_scf_fail > 0) then
      call env%error("transition state search aborted: the SCC did not "//&
         & "converge at least once", source)
   end if

   if (.not. result%converged) then
      if (allocated(result%message)) then
         call env%error("transition state search did not converge: "//&
            & result%message, source)
      else
         call env%error("transition state search did not converge", source)
      end if
   end if

end subroutine landscape_ts

!> Write the structured result `landscape_ts.json`.
!>
!> The counter `method` is the `p_landscape_ts_*` identifier of the run
!> (1 `prfo`, 2 `bprfo`, 3 `dimer`, 4 `neb`, 5 `cineb`, 6 `autoneb`,
!> 7 `string`, 8 `gsm`); it is numeric because the document carries numbers
!> only, the mapping is the one of `landscape_ts_method_name`.  The band
!> methods fill `result%converged` from the band criterion and have no
!> convergence report (`result%conv`), so their `history` block is empty while
!> the counters and the confirmation are complete.
subroutine write_json(env, mol, settings, result, n_scf_fail)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecule of the run (the core block's geometry, the converged saddle)
   type(TMolecule), intent(in) :: mol
   !> Transition state settings
   type(landscape_ts_settings), intent(in) :: settings
   !> Search result
   type(landscape_result), intent(in) :: result
   !> Number of unconverged SCCs of the run
   integer, intent(in) :: n_scf_fail

   type(landscape_json_document) :: doc = landscape_json_document()
   logical :: ok
   integer :: status

   status = landscape_json_status_ok
   if (.not. result%converged) status = landscape_json_status_incomplete
   ! A converged search that did not end on a first order saddle (or a run
   ! with an unconverged SCC) has no usable result; the task reports both
   ! through `env` and exits non-zero.
   if (settings%verify_ts .and. .not. result%is_saddle) then
      status = landscape_json_status_error
   end if
   if (n_scf_fail > 0) status = landscape_json_status_error
   call landscape_json_document_init(doc, 'ts')
   call landscape_json_add_result(doc, result, status)
   call landscape_json_add_core(doc, mol, result%energy, settings%input_file)
   call landscape_json_add_counter(doc, 'method', settings%method)
   call landscape_json_add_flag(doc, 'verify_ts', settings%verify_ts)
   call landscape_json_add_value(doc, 'kcal_per_mol', result%energy*autokcal)
   call landscape_json_write(env, doc, ok)

end subroutine write_json

!> Confirm the stationary point with a frequency analysis.
!>
!> Exactly one imaginary mode is the definition of a first order saddle point.
!> A different count is not silently accepted: the result carries it (`n_imag`)
!> and `is_saddle`, and the task turns a non-saddle into a failure (H2 of
!> `29-m4-independent-review.md`): the three outcomes - first order saddle, a
!> minimum (0 imaginary modes) and a higher order saddle (>= 2) - are printed
!> explicitly and a non-saddle aborts the run with a non-zero exit status (the
!> geometry and the log are still written for debugging).
!>
!> M4 hard constraint 2 (H1): the sign of a mode whose |frequency| falls inside
!> `0 +/- 3 cm^-1` is not meaningful at the default finite difference
!> resolution, so such a Hessian is recomputed once with `step = 2e-3` bohr and
!> `sccacc = 0.1` and the two spectra are printed.  The conversion is explicit:
!> the Hessian eigenvalues are in hartree/bohr^2 and a curvature `lambda`
!> corresponds to `sqrt(|lambda|)*autorcm` cm^-1, so the band is
!> `(3/autorcm)^2 ~ 3.4e-7` hartree/bohr^2 wide in eigenvalue units.
subroutine confirm_saddle(env, mol, chk, engine, settings, result)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure at the converged geometry
   type(TMolecule), intent(inout) :: mol
   !> Wavefunction data
   type(TRestart), intent(inout) :: chk
   !> Engine
   type(landscape_engine), intent(inout) :: engine
   !> Transition state settings
   type(landscape_ts_settings), intent(in) :: settings
   !> Search result, updated with the mode count and the imaginary frequency
   type(landscape_result), intent(inout) :: result

   type(landscape_vib_settings) :: vib_settings
   type(landscape_vib_result) :: vib
   real(wp), allocatable :: hess(:, :)
   real(wp) :: imag_freq, step_used, acc_used
   logical :: fail, band_hit, can_tighten, tightened
   integer :: n3

   n3 = 3*mol%n
   vib_settings%project_tr = settings%project_tr
   vib_settings%linear_tol = settings%linear_tol
   vib_settings%imag_cutoff = settings%imag_cutoff

   ! --- first analysis with the configured finite difference parameters ---
   call build_confirmation_hessian(env, mol, chk, engine, hess, fail)
   if (fail) return
   call landscape_vib_analysis(env, mol, hess, vib_settings, vib, fail)
   if (fail) then
      call env%warning("the saddle confirmation failed; the imaginary mode "//&
         & "count is unknown", 'landscape_ts')
      return
   end if
   imag_freq = lowest_internal_frequency(vib, n3)

   ! M4 hard constraint 2: near-zero curvature is resolution limited.  Tighten
   ! the finite differences once and re-decide the sign.
   band_hit = landscape_ts_curvature_band_hit(vib%freq, vib%n_zero, &
      & landscape_ts_curvature_band_cm1)
   can_tighten = landscape_ts_can_tighten(settings%hessian_source)
   step_used = engine%ctrl%hess_step
   acc_used = engine%ctrl%hess_acc
   tightened = .false.
   if (band_hit) then
      if (settings%printlevel > 0) then
         write(env%unit, '(1x,a,f8.3,a)') "near-zero curvature: a mode lies "//&
            & "inside 0 +/- ", landscape_ts_curvature_band_cm1, " cm^-1"
      end if
      if (can_tighten .and. (abs(step_used - ts_tight_hess_step) > 1.0e-12_wp .or. &
            & abs(acc_used - ts_tight_hess_acc) > 1.0e-12_wp)) then
         if (settings%printlevel > 0) then
            write(env%unit, '(1x,a,es10.3,a,es10.3)') &
               & "  recomputing the Hessian with step ", ts_tight_hess_step, &
               & " bohr and sccacc ", ts_tight_hess_acc
         end if
         engine%ctrl%hess_step = ts_tight_hess_step
         engine%ctrl%hess_acc = ts_tight_hess_acc
         call build_confirmation_hessian(env, mol, chk, engine, hess, fail)
         if (.not. fail) then
            call landscape_vib_analysis(env, mol, hess, vib_settings, vib, fail)
            tightened = .not. fail
         end if
         fail = .false.
         if (tightened) then
            if (settings%printlevel > 0) then
               write(env%unit, '(1x,a,f16.2,a,es10.3)') &
                  & "  tightened re-decision: lowest internal mode ", &
                  & lowest_internal_frequency(vib, n3), " cm^-1 (was ", imag_freq, ")"
            end if
            imag_freq = lowest_internal_frequency(vib, n3)
         else
            call env%warning("the tightened Hessian of the confirmation could "//&
               & "not be built; the sign of the near-zero mode keeps the "//&
               & "default resolution", 'landscape_ts')
         end if
      else if (.not. can_tighten) then
         call env%warning("a near-zero curvature cannot be re-decided: the "//&
            & "Hessian source '"//trim(settings%hessian_source)//"' does not "//&
            & "support tighter finite differences, so the sign of the mode is "//&
            & "resolution limited", 'landscape_ts')
      end if
   end if

   result%n_imag = vib%n_imag
   result%imag_freq = imag_freq
   result%is_saddle = (vib%n_imag == 1)

   if (settings%printlevel > 0) then
      write(env%unit, '(1x,a,i0)') "saddle confirmation - imaginary modes: ", &
         & result%n_imag
      if (result%n_imag > 0) then
         write(env%unit, '(1x,a,f16.2,a)') "most negative frequency: ", &
            & imag_freq, " cm^-1"
      end if
      select case (result%n_imag)
      case (1)
         write(env%unit, '(1x,a)') "confirmation: first order saddle point"
      case (0)
         write(env%unit, '(1x,a)') "confirmation: minimum (0 imaginary "//&
            & "modes) - not a transition state"
      case default
         write(env%unit, '(1x,a,i0,a)') "confirmation: higher order saddle ("//&
            & "point with ", result%n_imag, " imaginary modes) - not a "//&
            & "transition state"
      end select
   end if

   if (.not. result%is_saddle) then
      ! H2: a converged search that did not land on a first order saddle is a
      ! failure, not a warning.  The geometry and the log are still written by
      ! the caller; this only carries the diagnosis to the exit status.
      result%converged = .false.
      if (result%n_imag == 0) then
         result%message = "the located stationary point is a minimum (0 "//&
            & "imaginary modes), not a first order saddle point"
      else
         result%message = "the located stationary point is a higher order "//&
            & "saddle ("//trim(landscape_int_to_string(result%n_imag))//& 
            & " imaginary modes), not a first order saddle point"
      end if
      call env%warning(trim(result%message), 'landscape_ts')
   end if

end subroutine confirm_saddle

!> Build a fresh confirmation Hessian through the engine provider.
!>
!> The engine's finite difference parameters (`ctrl%hess_step`, `ctrl%hess_acc`)
!> are used, so the caller can tighten them before the call.  A failed build is
!> reported with a warning and `fail = .true.`; the caller then keeps the
!> unknown mode count instead of aborting the whole run.
subroutine build_confirmation_hessian(env, mol, chk, engine, hess, fail)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure at the converged geometry
   type(TMolecule), intent(inout) :: mol
   !> Wavefunction data
   type(TRestart), intent(inout) :: chk
   !> Engine
   type(landscape_engine), intent(inout) :: engine
   !> Hessian (3N, 3N), allocated on return
   real(wp), allocatable, intent(out) :: hess(:, :)
   !> Failure flag
   logical, intent(out) :: fail

   integer :: n3

   fail = .false.
   n3 = 3*mol%n
   if (allocated(hess)) deallocate(hess)
   allocate(hess(n3, n3), source=0.0_wp)
   call engine%hessian(env, mol, chk, hess, fail)
   if (fail) then
      call env%warning("the Hessian of the saddle confirmation could not be "//&
         & "built; the imaginary mode count is unknown", 'landscape_ts')
   end if

end subroutine build_confirmation_hessian

!> Lowest internal frequency (cm^-1) of a confirmation result, zero if none.
pure function lowest_internal_frequency(vib, n3) result(freq)
   !> Confirmation result
   type(landscape_vib_result), intent(in) :: vib
   !> Number of Cartesian components (3N)
   integer, intent(in) :: n3
   !> Lowest internal mode (cm^-1)
   real(wp) :: freq

   if (vib%n_vib > 0 .and. allocated(vib%freq)) then
      freq = minval(vib%freq(vib%n_zero + 1:n3))
   else
      freq = 0.0_wp
   end if

end function lowest_internal_frequency

!> True if an internal mode lies inside `0 +/- band` cm^-1.
!>
!> The band is expressed in cm^-1 because that is the unit of the frequency
!> table; the underlying Hessian eigenvalues are related by
!> `|frequency| = sqrt(|lambda|)*autorcm`, so the band corresponds to
!> `(band/autorcm)^2` in hartree/bohr^2 (about 3.4e-7 for the default 3 cm^-1).
!> Only the internal modes are tested: the projected zero modes are exactly
!> zero by construction and must not trigger the recomputation.
pure function landscape_ts_curvature_band_hit(freq, n_zero, band) result(hit)
   !> Frequency table (cm^-1, imaginary modes negative)
   real(wp), intent(in) :: freq(:)
   !> Number of projected zero modes at the front of `freq`
   integer, intent(in) :: n_zero
   !> Half width of the band (cm^-1)
   real(wp), intent(in) :: band
   !> True if a mode lies inside the band
   logical :: hit

   integer :: i

   hit = .false.
   do i = n_zero + 1, size(freq)
      if (abs(freq(i)) <= band) then
         hit = .true.
         return
      end if
   end do

end function landscape_ts_curvature_band_hit

!> True if a Hessian source can be recomputed with tighter finite differences.
!>
!> The curvature band re-decision of `confirm_saddle` rebuilds the Hessian with
!> `step = 2e-3` bohr and `sccacc = 0.1`.  Those two parameters belong to the
!> numerical provider alone: the model Hessian (Lindh-D2) has no displacement
!> grid and the `file` provider returns one fixed matrix, so their curvature
!> sign stays resolution limited and the caller reports that instead of
!> pretending to re-decide.  The name is matched through
!> `landscape_lowercase`, like every other `hessian_source` consumer: the CLI
!> stores one spelling but a library caller may pass any (`--lts-hessian
!> NUMERICAL` used to be accepted by the validator and then silently skip the
!> tightening).
!>
!> The predicate is public because it is the gate of that behaviour and is
!> pinned by a unit test, exactly like `landscape_ts_curvature_band_hit`.
pure function landscape_ts_can_tighten(hessian_source) result(can_tighten)
   !> Value of the `hessian_source` setting
   character(len=*), intent(in) :: hessian_source
   !> True if the source supports a tighter recomputation
   logical :: can_tighten

   can_tighten = landscape_lowercase(hessian_source) == 'numerical'

end function landscape_ts_can_tighten

!> Write `xtbts.xyz` and `xtbts.log`.
!>
!> The names follow the `<prefix>ts.xyz` scheme of the minimisation task
!> (`xtbopt.xyz`), so a default prefix of `xtb` yields the `xtbts.xyz` /
!> `xtbts.log` pair of PLAN section 9.5.  The log holds the header, the step
!> table and the frequency confirmation, so the file alone documents the run.
subroutine write_ts_files(env, mol, settings, result)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure (converged geometry)
   type(TMolecule), intent(in) :: mol
   !> Transition state settings
   type(landscape_ts_settings), intent(in) :: settings
   !> Search result
   type(landscape_result), intent(in) :: result

   character(len=:), allocatable :: xyzname, logname, name
   integer :: xyzunit, unit, ierr

   ! M3 (review 29): the file names pass through `get_namespace`, exactly like
   ! the scan (`xtbscan.log`) and vib (`vib_normal_modes`) tasks, so a
   ! `--namespace ns --lts` run writes `ns.xtbts.xyz` / `ns.xtbts.log` instead
   ! of colliding with a plain run in the same directory.
   name = get_namespace(trim(settings%output_prefix)//'ts.xyz')
   xyzname = name
   open(newunit=xyzunit, file=xyzname, action='write', iostat=ierr)
   if (ierr /= 0) then
      call env%error("cannot open "//xyzname, 'landscape_ts')
      return
   end if
   call writeMolecule(mol, xyzunit, format=fileType%xyz, energy=result%energy, &
      & gnorm=result%gnorm)
   close(xyzunit)
   write(env%unit, '(1x,a,a)') "wrote ", xyzname

   logname = get_namespace(trim(settings%output_prefix)//'ts.log')
   open(newunit=unit, file=logname, action='write', iostat=ierr)
   if (ierr /= 0) then
      call env%error("cannot open "//logname, 'landscape_ts')
      unit = -1
      return
   end if
   write(unit, '(1x,a,a)') "landscape transition state search: ", &
      & landscape_ts_method_name(settings%method)
   ! the effective (read back) temperature, not the echoed setting (M1)
   write(unit, '(1x,a,f12.3,a)') "electronic temperature (frozen): ", &
      & result%etemp, " K"
   write(unit, '(1x,a,es15.6,a)') "initial trust radius: ", settings%trust, " bohr"
   write(unit, '(1x,a,a)') "Hessian source: ", trim(settings%hessian_source)
   if (allocated(result%message)) then
      write(unit, '(1x,a,a)') "termination: ", result%message
   end if
   call result%write_log(unit)
   call result%conv%write_log(unit)
   if (result%n_imag >= 0) then
      write(unit, '(1x,a,i0)') "confirmation imaginary modes: ", result%n_imag
      write(unit, '(1x,a,l1)') "is a first order saddle: ", result%is_saddle
      if (result%n_imag > 0) then
         write(unit, '(1x,a,f16.2,a)') "most negative frequency: ", &
            & result%imag_freq, " cm^-1"
      end if
   end if
   close(unit)
   write(env%unit, '(1x,a,a)') "wrote ", logname

end subroutine write_ts_files

!> Chain-of-states branch of the transition state search.
!>
!> All chain-of-states methods need two endpoints.  The reactant is the
!> geometry already read by xtb; the product is read from `--lts-product` (a
!> single frame xyz).  The input contract and every rejection path are
!> documented on `landscape_path_settings`; nothing is guessed and nothing is
!> silently ignored: a missing product, an unreadable file, a frame count
!> above one, a multi-frame input, a different number of atoms or a different
!> atom order are all `env%error`s.
!>
!> The method selects the search:
!>
!>  * `neb`, `cineb`: one band relaxation
!>    (`xtb_landscape_path_neb::landscape_band_optimize`);
!>  * `autoneb`: incremental path building
!>    (`xtb_landscape_path_autoneb::landscape_autoneb_optimize`) - the path is
!>    written to `path_autoneb.xyz`, which is also the restart file of a
!>    subsequent run (`--lts-autoneb-restart`);
!>  * `string`, `gsm`: growing string
!>    (`xtb_landscape_path_gsm::landscape_string_optimize`) - the path is
!>    written to `path.xyz` in the band format.
!>
!> After the search the **highest energy image** (the climbing image when one
!> was ignited) becomes the transition state candidate: it is written to
!> `mol%xyz` and confirmed with the M3 frequency analysis, exactly like the M4
!> methods (exactly one imaginary mode is required, a different count aborts
!> the run).  AutoNEB and the string additionally write `path_ts_list.xyz`,
!> the list of every saddle candidate the search produced (the stage
!> candidates of AutoNEB, the local maxima of the final string), each with
!> its energy and its imaginary mode count.
subroutine landscape_ts_band(env, mol, chk, engine, settings, result)
   character(len=*), parameter :: source = 'landscape_ts_band'
   !> Computational environment
   type(TEnvironment), intent(inout), target :: env
   !> Molecular structure (reactant on entry, candidate on return)
   type(TMolecule), intent(inout), target :: mol
   !> Wavefunction data
   type(TRestart), intent(inout), target :: chk
   !> Engine
   type(landscape_engine), intent(inout), target :: engine
   !> Transition state settings
   type(landscape_ts_settings), intent(in) :: settings
   !> Search result
   type(landscape_result), intent(inout) :: result

   type(landscape_engine_model) :: model
   type(landscape_band_type) :: band
   type(landscape_autoneb_type) :: auto
   type(landscape_string_type) :: string
   real(wp), allocatable :: reactant(:, :), product(:, :), images(:, :, :)
   real(wp), allocatable :: energies(:), forces(:, :, :), restart(:, :, :)
   real(wp), allocatable :: candidates(:, :, :), cand_energies(:)
   integer :: k
   logical :: fail, has_candidates

   model%engine => engine
   model%env => env
   model%mol => mol
   model%chk => chk

   call read_band_product(env, mol, settings, product, fail)
   if (fail) then
      result%message = "the product endpoint could not be read"
      return
   end if

   allocate(reactant(3, mol%n), source=mol%xyz)
   if (settings%path%relax_ends) then
      call relax_band_endpoints(env, model, settings, reactant, product)
      mol%xyz = reactant
   end if

   ! --- method dispatch --------------------------------------------------
   has_candidates = .false.
   k = 0
   select case (settings%method)
   case (p_landscape_ts_autoneb)
      ! restart reuse: the multi-frame path of an earlier run is adopted as
      ! the initial path, so the stages it already contains are not repeated
      ! (the semantics are documented on `xtb_landscape_path_autoneb`)
      if (allocated(settings%path%autoneb_restart_file)) then
         call landscape_ts_read_path(env, mol, settings%path%autoneb_restart_file, &
            & reactant, product, restart, fail)
         if (fail) then
            result%message = "the AutoNEB restart path could not be read"
            return
         end if
         call landscape_autoneb_optimize(env, model, settings%path, reactant, &
            & product, images, energies, forces, auto, initial_images=restart)
      else
         call landscape_autoneb_optimize(env, model, settings%path, reactant, &
            & product, images, energies, forces, auto)
      end if
      if (.not. allocated(images) .or. .not. allocated(auto%energies)) then
         result%message = auto%message
         return
      end if
      result%n_eval = auto%n_eval
      result%n_steps = auto%n_iter
      result%converged = auto%converged
      result%etemp = engine%effective_etemp()
      result%message = band_message(auto%converged, auto%message)
      ! A run that stopped at `max_paths` or `max_depth` did not build the
      ! requested path, even when the band of its last stage converged.  The
      ! cap message has to survive into the result (and therefore into
      ! `xtbts.log`), otherwise the run reports "converged" while it holds
      ! fewer images than requested (report 38 L-2).  It is reported as a
      ! warning, not as a failure: the band that was built is converged.
      if (auto%cap_reached) then
         result%message = auto%message
         call env%warning("the AutoNEB run did not reach the requested path: "//&
            & auto%message, source)
      end if
      if (settings%write_files) then
         call write_path_frames(env, mol, images, auto%energies, 'path_autoneb.xyz')
      end if
      if (settings%printlevel > 0) then
         call write_path_table(env, 'autoneb', auto%energies, auto%distances, &
            & auto%saddle_image, auto%climbing_image)
         call write_autoneb_table(env, auto)
      end if
      k = auto%saddle_image
      if (auto%climbing_image > 0) k = auto%climbing_image
      if (auto%n_stages > 0) then
         candidates = auto%stage_images(:, :, 1:auto%n_stages)
         cand_energies = auto%stage_energies(1:auto%n_stages)
         has_candidates = .true.
      end if
   case default
      if (landscape_ts_method_is_string(settings%method)) then
         call landscape_string_optimize(env, model, settings%path, reactant, &
            & product, images, energies, forces, string)
         if (.not. allocated(images)) then
            result%message = string%message
            return
         end if
         result%n_eval = string%n_eval
         result%n_steps = string%n_iter + string%n_climb_iter
         result%converged = string%converged
         result%etemp = engine%effective_etemp()
         result%message = band_message(string%converged, string%message)
         if (settings%write_files) then
            call write_path_frames(env, mol, images, string%energies, 'path.xyz')
         end if
         if (settings%printlevel > 0) then
            call write_path_table(env, 'string', string%energies, string%distances, &
               & string%highest_node, string%climbing_node)
            call write_string_table(env, string)
         end if
         k = string%highest_node
         call collect_string_candidates(images, string%energies, candidates, &
            & cand_energies)
         has_candidates = allocated(candidates)
      else
         call landscape_band_optimize(env, model, settings%path, reactant, product, &
            & images, energies, forces, band)
         if (.not. allocated(images) .or. .not. allocated(band%energies)) then
            result%message = band%message
            return
         end if
         result%n_eval = band%n_eval
         result%n_steps = band%n_iter
         result%converged = band%converged
         result%etemp = engine%effective_etemp()
         result%message = band_message(band%converged, band%message)
         if (settings%write_files) then
            call write_path_frames(env, mol, images, band%energies, 'path.xyz')
         end if
         if (settings%printlevel > 0) then
            call write_path_table(env, 'band', band%energies, band%distances, &
               & band%highest_image, band%climbing_image)
         end if
         ! transition state candidate: the climbing image when CI-NEB was
         ! run, otherwise the highest energy image of the band
         k = band%highest_image
         if (band%climbing_image > 0) k = band%climbing_image
      end if
   end select

   if (k < 1) then
      call env%error("the path has no internal node to use as a transition "//&
         & "state candidate", source)
      result%message = "no transition state candidate"
      return
   end if

   call landscape_ts_finish(env, model, mol, chk, engine, settings, images, k, &
      & candidates, cand_energies, has_candidates, result)

end subroutine landscape_ts_band

!> Message of a chain-of-states result: `converged` for a converged search.
pure function band_message(converged, message) result(text)
   !> Convergence flag
   logical, intent(in) :: converged
   !> Message of the search
   character(len=*), intent(in) :: message
   !> Message of the result object
   character(len=:), allocatable :: text

   if (converged) then
      text = "converged"
   else
      text = message
   end if

end function band_message

!> Common tail of the chain-of-states searches: the candidate geometry.
!>
!> The candidate `images(:, :, k)` is written to `mol%xyz` and evaluated; its
!> energy and gradient measures fill the result.  `path_ts_list.xyz` is written
!> first (AutoNEB and the string only), because the mode counting of the extra
!> candidates works on `mol%xyz` and restores it.  The frequency confirmation
!> and the `xtbts.*` files are the same for every chain-of-states method.
subroutine landscape_ts_finish(env, model, mol, chk, engine, settings, images, k, &
      & candidates, cand_energies, has_candidates, result)
   character(len=*), parameter :: source = 'landscape_ts_band'
   !> Computational environment
   type(TEnvironment), intent(inout), target :: env
   !> Potential
   class(landscape_engine_model), intent(inout) :: model
   !> Molecular structure (candidate on return)
   type(TMolecule), intent(inout), target :: mol
   !> Wavefunction data
   type(TRestart), intent(inout), target :: chk
   !> Engine
   type(landscape_engine), intent(inout), target :: engine
   !> Transition state settings
   type(landscape_ts_settings), intent(in) :: settings
   !> Path images
   real(wp), intent(in) :: images(:, :, :)
   !> Index of the candidate node
   integer, intent(in) :: k
   !> Extra saddle candidates (unallocated when there are none)
   real(wp), allocatable, intent(in) :: candidates(:, :, :)
   !> Energy of the extra candidates (hartree)
   real(wp), allocatable, intent(in) :: cand_energies(:)
   !> True if the extra candidates have to be written
   logical, intent(in) :: has_candidates
   !> Search result
   type(landscape_result), intent(inout) :: result

   real(wp), allocatable :: x0(:), gflat(:)
   real(wp) :: energy
   integer :: n_scf_fail
   logical :: fail

   mol%xyz = images(:, :, k)
   x0 = geometry_pack(mol%xyz)
   allocate(gflat(3*mol%n))
   call model%evaluate(x0, energy, gflat, fail)
   if (fail) then
      call env%error("the transition state candidate could not be evaluated", source)
      result%message = "the transition state candidate could not be evaluated"
      return
   end if
   result%energy = energy
   result%gnorm = sqrt(sum(gflat**2))
   result%gmax = maxval(abs(gflat))

   if (has_candidates .and. settings%write_files) then
      call write_ts_candidates(env, mol, chk, engine, settings, candidates, &
         & cand_energies, trim(landscape_ts_method_name(settings%method)))
   end if

   if (settings%verify_ts) then
      call confirm_saddle(env, mol, chk, engine, settings, result)
   end if

   ! Same sampling rule as the single-geometry search above (report 50
   ! finding M-1): the counter belongs after the confirmation Hessian so that
   ! the printed value, the JSON counter and the engine report of the run
   ! agree.
   result%n_cache_miss = engine%n_eval

   if (settings%printlevel > 0) then
      call result%write_log(env%unit)
      call result%conv%write_log(env%unit)
   end if

   if (settings%write_files) then
      call write_ts_files(env, mol, settings, result)
   end if

   n_scf_fail = engine%n_scf_fail
   if (n_scf_fail > 0) then
      call env%error("transition state search aborted: the SCC did not "//&
         & "converge at least once", source)
   end if
   if (.not. result%converged) then
      if (allocated(result%message)) then
         call env%error("transition state search did not converge: "//&
            & result%message, source)
      else
         call env%error("transition state search did not converge", source)
      end if
   end if

end subroutine landscape_ts_finish

!> Read and validate the product endpoint of a band search.
!>
!> Implements the input contract of `landscape_path_settings`:
!>
!>  * no `--lts-product`: a multi-frame input is named as such (the
!>    convenience path of using frame 2 is deliberately **not** implemented,
!>    a second parser is not worth the ambiguity), a single frame input is
!>    reported as a missing product;
!>  * `--lts-product` together with a multi-frame input: rejected, the two
!>    endpoints would be ambiguous;
!>  * the product must be a single frame xyz (any other extension, a missing
!>    file, an empty file and more than one frame are rejected);
!>  * the atom count and the atom order must match the reactant exactly.
!>
!> `env` collects a readable diagnostic for every branch; `fail` is set when
!> the endpoint could not be produced.
subroutine read_band_product(env, mol, settings, product, fail)
   character(len=*), parameter :: source = 'landscape_ts_band'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Reactant structure (atom list used for the validation)
   type(TMolecule), intent(in) :: mol
   !> Transition state settings
   type(landscape_ts_settings), intent(in) :: settings
   !> Product geometry (bohr), allocated on return
   real(wp), allocatable, intent(out) :: product(:, :)
   !> Failure flag
   logical, intent(out) :: fail

   type(TMolecule) :: prodmol
   character(len=:), allocatable :: fname
   integer :: unit, ierr, nframes, ninput, nlog, first, ia

   fail = .true.

   ! frame count of the positional input file, if it can be inspected
   ninput = -1
   if (allocated(settings%path%input_file)) then
      if (getFileType(settings%path%input_file) == fileType%xyz) then
         ninput = landscape_ts_count_frames(settings%path%input_file)
      end if
   end if

   if (.not. allocated(settings%path%product_file)) then
      if (ninput > 1) then
         call env%error("the input holds "//landscape_int_to_string(ninput)//&
            & " frames; pass the product with --lts-product", source)
      else
         call env%error("a band search needs a product geometry; pass it "//&
            & "with --lts-product <file>", source)
      end if
      return
   end if

   if (ninput > 1) then
      call env%error("the input holds "//landscape_int_to_string(ninput)//" frames "//&
         & "and --lts-product was given as well; the two endpoints are "//&
         & "ambiguous", source)
      return
   end if

   fname = settings%path%product_file
   if (getFileType(fname) /= fileType%xyz) then
      call env%error("the --lts-product geometry must be a single frame xyz "//&
         & "file, got '"//fname//"'", source)
      return
   end if
   nframes = landscape_ts_count_frames(fname)
   if (nframes < 0) then
      call env%error("cannot open the product geometry '"//fname//"'", source)
      return
   end if
   if (nframes /= 1) then
      call env%error("the product geometry '"//fname//"' holds "//&
         & landscape_int_to_string(nframes)//" frames; exactly one is required", source)
      return
   end if

   open(newunit=unit, file=fname, action='read', status='old', iostat=ierr)
   if (ierr /= 0) then
      call env%error("cannot open the product geometry '"//fname//"'", source)
      return
   end if
   nlog = env%nLog
   call readMolecule(env, prodmol, unit, fileType%xyz)
   close(unit)
   if (landscape_log_failed(env, nlog)) then
      return
   end if
   if (prodmol%n <= 0) then
      call env%error("the product geometry '"//fname//"' could not be read", source)
      return
   end if
   if (prodmol%n /= mol%n) then
      call env%error("the product geometry has "//landscape_int_to_string(prodmol%n)//&
         & " atoms but the reactant has "//landscape_int_to_string(mol%n), source)
      return
   end if
   first = 0
   do ia = 1, mol%n
      if (prodmol%at(ia) /= mol%at(ia) .and. first == 0) first = ia
   end do
   if (first /= 0) then
      call env%error("the atom order of the product geometry differs from "//&
         & "the reactant (first mismatch at atom "//landscape_int_to_string(first)//")", source)
      return
   end if

   allocate(product(3, mol%n), source=prodmol%xyz)
   if (settings%printlevel > 0) then
      write(env%unit, '(1x,a,a)') "band: product geometry ", fname
   end if
   fail = .false.

end subroutine read_band_product

!> Number of frames of an xyz file, `-1` if the file cannot be opened.
!>
!> The counter is deliberately not a parser: it reads the atom count of a
!> frame, skips the comment and the atom lines and looks for the next atom
!> count.  An empty file returns zero frames, which the caller reports as a
!> rejected product; a truncated file returns the number of complete frames
!> it contains.
function landscape_ts_count_frames(path) result(nframes)
   !> Path of the xyz file
   character(len=*), intent(in) :: path
   !> Number of complete frames, `-1` if the file cannot be read
   integer :: nframes

   character(len=512) :: line
   integer :: unit, iostat_, natom, i

   nframes = -1
   open(newunit=unit, file=path, action='read', status='old', iostat=iostat_)
   if (iostat_ /= 0) return
   nframes = 0
   do
      read(unit, *, iostat=iostat_) natom
      if (iostat_ /= 0) exit
      read(unit, '(a)', iostat=iostat_) line
      if (iostat_ /= 0) exit
      do i = 1, natom
         read(unit, '(a)', iostat=iostat_) line
         if (iostat_ /= 0) exit
      end do
      if (iostat_ /= 0) exit
      nframes = nframes + 1
   end do
   close(unit)

end function landscape_ts_count_frames

!> Relax the two endpoints of a band with the M1 minimiser.
!>
!> Only called with `landscape_path_settings%relax_ends`; the electronic
!> temperature, the SCC accuracy and the restart policy are those of the
!> surrounding run (the same engine control object is used), so the whole
!> band sequence is evaluated under one consistent policy.  The energy and
!> the gradient norm of both endpoints are printed **before and after** the
!> relaxation, so the effect of the option is visible in the log.
subroutine relax_band_endpoints(env, model, settings, reactant, product)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Potential
   class(landscape_engine_model), intent(inout) :: model
   !> Transition state settings
   type(landscape_ts_settings), intent(in) :: settings
   !> Reactant geometry (bohr), relaxed in place
   real(wp), intent(inout) :: reactant(:, :)
   !> Product geometry (bohr), relaxed in place
   real(wp), intent(inout) :: product(:, :)

   type(landscape_opt_settings) :: opt
   real(wp) :: e0, g0, e1, g1
   logical :: fail

   opt%max_steps = settings%max_steps
   opt%e_tol = settings%e_tol
   opt%g_max = settings%g_max
   opt%g_rms = settings%g_rms
   opt%d_max = settings%d_max
   opt%d_rms = settings%d_rms
   opt%printlevel = 0
   opt%write_files = .false.

   call relax_one_endpoint(env, model, opt, settings, reactant, "reactant", &
      & e0, g0, e1, g1, fail)
   if (fail) return
   call relax_one_endpoint(env, model, opt, settings, product, "product", &
      & e0, g0, e1, g1, fail)

end subroutine relax_band_endpoints

!> Relax one endpoint and print the energy and gradient norm before/after.
subroutine relax_one_endpoint(env, model, opt, settings, geometry, label, &
      & e_before, g_before, e_after, g_after, fail)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Potential
   class(landscape_engine_model), intent(inout) :: model
   !> Minimiser settings
   type(landscape_opt_settings), intent(in) :: opt
   !> Transition state settings (print level)
   type(landscape_ts_settings), intent(in) :: settings
   !> Endpoint geometry (bohr), relaxed in place
   real(wp), intent(inout) :: geometry(:, :)
   !> Label for the report
   character(len=*), intent(in) :: label
   !> Energy before the relaxation (hartree)
   real(wp), intent(out) :: e_before
   !> Gradient norm before the relaxation (hartree/bohr)
   real(wp), intent(out) :: g_before
   !> Energy after the relaxation (hartree)
   real(wp), intent(out) :: e_after
   !> Gradient norm after the relaxation (hartree/bohr)
   real(wp), intent(out) :: g_after
   !> Failure flag
   logical, intent(out) :: fail

   type(landscape_result) :: res
   real(wp), allocatable :: x0(:), xfinal(:), gflat(:)

   fail = .false.
   x0 = geometry_pack(geometry)
   allocate(gflat(3*size(geometry, 2)))
   call model%evaluate(x0, e_before, gflat, fail)
   if (fail) return
   g_before = sqrt(sum(gflat**2))

   call landscape_optimize_model(env, model, opt, x0, res, xfinal, &
      & frozen=landscape_frozen_mask(size(geometry, 2)))
   if (.not. allocated(xfinal)) then
      fail = .true.
      return
   end if
   call geometry_unpack(xfinal, geometry)
   call model%evaluate(xfinal, e_after, gflat, fail)
   if (fail) return
   g_after = sqrt(sum(gflat**2))

   if (settings%printlevel > 0) then
      write(env%unit, '(1x,a,a,a,f20.12,a,es12.5)') "band: endpoint relaxation (", &
         & label, ") before: E ", e_before, " Eh, |g| ", g_before
      write(env%unit, '(1x,a,a,a,f20.12,a,es12.5)') "band: endpoint relaxation (", &
         & label, ") after:  E ", e_after, " Eh, |g| ", g_after
   end if
   if (.not. res%converged) then
      call env%warning("the "//label//" endpoint optimisation did not "//&
         & "converge", 'landscape_ts_band')
   end if

end subroutine relax_one_endpoint

!> Write one frame per image of a path to an xyz file.
!>
!> Frame 1 is the reactant and the last frame the product, in the atom order
!> of the input; the comment line of every frame carries the image energy and
!> the frame number (`writeMolecule`).  The name passes through
!> `get_namespace`, exactly like the other landscape outputs, so
!> `--namespace ns` yields `ns.<name>`.
subroutine write_path_frames(env, mol, images, energies, name)
   character(len=*), parameter :: source = 'landscape_ts_band'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure (atom identities are copied from it)
   type(TMolecule), intent(in) :: mol
   !> Path images `(3, nat, n_images)`
   real(wp), intent(in) :: images(:, :, :)
   !> Energy of every image (hartree)
   real(wp), intent(in) :: energies(:)
   !> File name, relative to the namespace
   character(len=*), intent(in) :: name

   type(TMolecule) :: buffer
   character(len=:), allocatable :: fname
   integer :: unit, ierr, k, nframe

   if (size(energies) < size(images, 3)) then
      call env%error("the path has fewer energies than images", source)
      return
   end if
   fname = get_namespace(name)
   open(newunit=unit, file=fname, action='write', iostat=ierr)
   if (ierr /= 0) then
      call env%error("cannot open "//fname, source)
      return
   end if
   ! `init` instead of a derived type assignment (gfortran 12.2 stride bug,
   ! see the comment in `xtb_landscape_tasks_scan`)
   call init(buffer, mol%at, mol%xyz, chrg=mol%chrg, uhf=mol%uhf, &
      & lattice=mol%lattice, pbc=mol%pbc)
   nframe = 0
   do k = 1, size(images, 3)
      nframe = nframe + 1
      buffer%xyz = images(:, :, k)
      call writeMolecule(buffer, unit, format=fileType%xyz, &
         & energy=energies(k), number=nframe)
   end do
   close(unit)
   write(env%unit, '(1x,a,a,a,i0,a)') "wrote ", fname, " with ", nframe, &
      & " frames"

end subroutine write_path_frames

!> Write the saddle candidates of a search to `path_ts_list.xyz`.
!>
!> One frame per candidate; the comment line carries the energy and the
!> **imaginary mode count** of the candidate, measured with the same M3
!> frequency analysis as the confirmation of the final candidate (the Hessian
!> of the configured provider at the candidate geometry, tight finite
!> differences as configured).  A candidate whose mode count could not be
!> determined is written with `-1` and a warning; the file is a diagnostic, not
!> a decision, and it never aborts the run.
subroutine write_ts_candidates(env, mol, chk, engine, settings, candidates, &
      & cand_energies, label)
   character(len=*), parameter :: source = 'landscape_ts_band'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure (used as the carrier of the atom list)
   type(TMolecule), intent(inout) :: mol
   !> Wavefunction data
   type(TRestart), intent(inout) :: chk
   !> Engine
   type(landscape_engine), intent(inout) :: engine
   !> Transition state settings
   type(landscape_ts_settings), intent(in) :: settings
   !> Candidate geometries `(3, nat, n_candidates)`
   real(wp), allocatable, intent(in) :: candidates(:, :, :)
   !> Candidate energies (hartree)
   real(wp), allocatable, intent(in) :: cand_energies(:)
   !> Method name for the comment line
   character(len=*), intent(in) :: label

   character(len=:), allocatable :: fname
   integer :: unit, ierr, c, ia, n_imag
   integer :: n_cand

   if (.not. allocated(candidates)) return
   if (size(candidates, 3) < 1) return
   n_cand = size(candidates, 3)
   fname = get_namespace('path_ts_list.xyz')
   open(newunit=unit, file=fname, action='write', iostat=ierr)
   if (ierr /= 0) then
      call env%error("cannot open "//fname, source)
      return
   end if
   write(env%unit, '(1x,a)') "   candidate   energy/Eh    imaginary modes"
   do c = 1, n_cand
      call candidate_modes(env, mol, chk, engine, settings, candidates(:, :, c), &
         & n_imag)
      write(unit, '(i0)') mol%n
      if (c <= size(cand_energies)) then
         write(unit, '(a,f20.12,a,i0,a,a)') " energy: ", cand_energies(c), &
            & " Eh, imaginary modes: ", n_imag, " method: ", label
      else
         write(unit, '(a,i0,a,a)') " imaginary modes: ", n_imag, &
            & " method: ", label
      end if
      do ia = 1, mol%n
         write(unit, '(a4,1x,3f24.14)') mol%sym(ia), &
            & candidates(:, ia, c)*autoaa
      end do
      if (c <= size(cand_energies)) then
         write(env%unit, '(1x,i11,1x,f20.12,1x,i15)') c, cand_energies(c), n_imag
      end if
   end do
   close(unit)
   write(env%unit, '(1x,a,a,a,i0,a)') "wrote ", fname, " with ", n_cand, &
      & " candidate frames"

end subroutine write_ts_candidates

!> Imaginary mode count of one candidate geometry.
!>
!> The Hessian of the configured provider is built at the candidate geometry
!> through the same routine as the confirmation of the final candidate, so the
!> mode count and the projection are those of `--lvib`.  On return `mol%xyz`
!> holds the geometry it had on entry.  A failed Hessian or analysis leaves
!> `n_imag = -1` (unknown), which the caller writes out and reports.
subroutine candidate_modes(env, mol, chk, engine, settings, geometry, n_imag)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure (geometry carrier, restored on return)
   type(TMolecule), intent(inout) :: mol
   !> Wavefunction data
   type(TRestart), intent(inout) :: chk
   !> Engine
   type(landscape_engine), intent(inout) :: engine
   !> Transition state settings
   type(landscape_ts_settings), intent(in) :: settings
   !> Candidate geometry (bohr)
   real(wp), intent(in) :: geometry(:, :)
   !> Imaginary mode count, `-1` if it could not be determined
   integer, intent(out) :: n_imag

   type(landscape_vib_settings) :: vib_settings
   type(landscape_vib_result) :: vib
   real(wp), allocatable :: hess(:, :), saved(:, :)
   logical :: fail

   n_imag = -1
   if (size(geometry, 1) /= 3 .or. size(geometry, 2) /= mol%n) return
   if (.not. allocated(mol%sym)) return
   saved = mol%xyz
   mol%xyz = geometry
   vib_settings%project_tr = settings%project_tr
   vib_settings%linear_tol = settings%linear_tol
   vib_settings%imag_cutoff = settings%imag_cutoff
   call build_confirmation_hessian(env, mol, chk, engine, hess, fail)
   if (.not. fail) then
      call landscape_vib_analysis(env, mol, hess, vib_settings, vib, fail)
      if (.not. fail) n_imag = vib%n_imag
   end if
   mol%xyz = saved

end subroutine candidate_modes

!> Print the stage table of an AutoNEB run.
subroutine write_autoneb_table(env, auto)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> AutoNEB result
   type(landscape_autoneb_type), intent(in) :: auto

   integer :: s

   write(env%unit, '(1x,a)') "   stage    images       E(candidate)/Eh   barrier/kcal/mol"
   do s = 1, auto%n_stages
      write(env%unit, '(1x,i6,1x,i8,1x,f20.12,1x,f14.4)') s, &
         & auto%stage_n_images(s), auto%stage_energies(s), &
         & (auto%stage_energies(s) - auto%energies(1))*autokcal
   end do
   if (auto%cap_reached) then
      write(env%unit, '(1x,a,a)') "autoneb: budget limit reached: ", &
         & trim(auto%message)
   end if

end subroutine write_autoneb_table

!> Print the growth and convergence report of a string run.
subroutine write_string_table(env, string)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> String result
   type(landscape_string_type), intent(in) :: string

   write(env%unit, '(1x,a,i0,a,i0,a,es12.4,a)') "string: nodes ", string%n_nodes, &
      & ", growth steps ", string%n_grow, ", spacing ", string%spacing, " bohr"
   write(env%unit, '(1x,a,i0,a,i0,a,es12.4,a,es12.4)') "string: optimisation ", &
      & string%n_iter, ", climbing iterations ", string%n_climb_iter, &
      & ", f_rms_avg ", string%f_rms_avg, ", f_rms_max ", string%f_rms_max
   if (string%spacing_adjusted) then
      write(env%unit, '(1x,a)') "string: the node spacing condition enlarged "//&
         & "the requested growth step"
   end if

end subroutine write_string_table

!> Collect the saddle candidates of a converged string.
!>
!> Every interior node that is a strict local maximum of the energy profile is
!> a candidate; a monotone profile (no interior maximum) falls back to the
!> highest energy node, so the list is never empty for a path with an interior
!> node.
subroutine collect_string_candidates(images, energies, candidates, cand_energies)
   !> String nodes `(3, nat, n_nodes)`
   real(wp), intent(in) :: images(:, :, :)
   !> Node energies (hartree)
   real(wp), intent(in) :: energies(:)
   !> Candidate geometries, allocated on return
   real(wp), allocatable, intent(out) :: candidates(:, :, :)
   !> Candidate energies, allocated on return
   real(wp), allocatable, intent(out) :: cand_energies(:)

   integer, allocatable :: index(:)
   integer :: k, n, n_cand, best

   n = size(images, 3)
   allocate(index(n), source=0)
   n_cand = 0
   best = 0
   do k = 2, n - 1
      ! A block if, not `best == 0 .or. ...`: Fortran does not require the
      ! right operand to be short circuited, so the first iteration would read
      ! `energies(0)`, which a bounds checked build rejects.
      if (best == 0) then
         best = k
      else if (energies(k) > energies(best)) then
         best = k
      end if
      if (energies(k) > energies(k - 1) .and. energies(k) > energies(k + 1)) then
         n_cand = n_cand + 1
         index(n_cand) = k
      end if
   end do
   if (n_cand == 0) then
      if (best < 1) return
      n_cand = 1
      index(1) = best
   end if
   allocate(candidates(3, size(images, 2), n_cand), source=0.0_wp)
   allocate(cand_energies(n_cand), source=0.0_wp)
   do k = 1, n_cand
      candidates(:, :, k) = images(:, :, index(k))
      cand_energies(k) = energies(index(k))
   end do

end subroutine collect_string_candidates

!> Read the multi-frame path of an earlier AutoNEB run.
!>
!> The file is the `path_autoneb.xyz` of a previous run (or any xyz file with
!> the same frames); the reader is deliberately strict: every frame must have
!> exactly `mol%n` atoms, the first frame must be the reactant and the last
!> one the product (the two endpoints are fixed images, so a path with other
!> endpoints describes another reaction).  The geometries are converted from
!> the angstrom of the xyz file to the bohr of the suite; the atom order is
!> taken to be that of the reactant, which is the order a previous run of the
!> same system wrote.
subroutine landscape_ts_read_path(env, mol, path, reactant, product, frames, fail)
   character(len=*), parameter :: source = 'landscape_ts_read_path'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure (atom count and the reference endpoints)
   type(TMolecule), intent(in) :: mol
   !> Path of the xyz file
   character(len=*), intent(in) :: path
   !> Reactant geometry (bohr, `(3, nat)`)
   real(wp), intent(in) :: reactant(:, :)
   !> Product geometry (bohr, `(3, nat)`)
   real(wp), intent(in) :: product(:, :)
   !> Path frames (bohr), allocated on return
   real(wp), allocatable, intent(out) :: frames(:, :, :)
   !> Failure flag
   logical, intent(out) :: fail

   character(len=8) :: symbol
   integer :: unit, iostat_, natom, nframe, ifr, ia, i
   real(wp) :: xx, yy, zz
   real(wp), allocatable :: buffer(:, :, :)
   logical :: too_many

   fail = .true.
   open(newunit=unit, file=path, action='read', status='old', iostat=iostat_)
   if (iostat_ /= 0) then
      call env%error("cannot open the AutoNEB restart path '"//path//"'", source)
      return
   end if
   ! first pass: count the frames
   nframe = 0
   do
      read(unit, *, iostat=iostat_) natom
      if (iostat_ /= 0) exit
      read(unit, '(a)', iostat=iostat_) symbol
      if (iostat_ /= 0) exit
      do i = 1, natom
         read(unit, '(a)', iostat=iostat_) symbol
         if (iostat_ /= 0) exit
      end do
      if (iostat_ /= 0) exit
      nframe = nframe + 1
   end do
   close(unit)
   if (nframe < 3) then
      call env%error("the AutoNEB restart path '"//path//"' holds fewer than "//&
         & "three frames", source)
      return
   end if

   allocate(buffer(3, mol%n, nframe), source=0.0_wp)
   open(newunit=unit, file=path, action='read', status='old', iostat=iostat_)
   if (iostat_ /= 0) then
      call env%error("cannot open the AutoNEB restart path '"//path//"'", source)
      return
   end if
   do ifr = 1, size(buffer, 3)
      read(unit, *, iostat=iostat_) natom
      if (iostat_ /= 0) exit
      if (natom /= mol%n) then
         call env%error("the AutoNEB restart path '"//path//"' holds a frame "//&
            & "with a different number of atoms", source)
         close(unit)
         return
      end if
      read(unit, '(a)', iostat=iostat_) symbol
      if (iostat_ /= 0) exit
      too_many = .false.
      do ia = 1, mol%n
         read(unit, *, iostat=iostat_) symbol, xx, yy, zz
         if (iostat_ /= 0) then
            too_many = .true.
            exit
         end if
         buffer(1, ia, ifr) = angstrom_to_bohr(xx)
         buffer(2, ia, ifr) = angstrom_to_bohr(yy)
         buffer(3, ia, ifr) = angstrom_to_bohr(zz)
      end do
      if (too_many) exit
   end do
   close(unit)
   ! the endpoints are compared through their pairwise distances: the band
   ! kernel stores the product rotated into the frame of the reactant when the
   ! alignment is on, so a coordinate comparison would reject the driver's own
   ! restart file
   if (.not. geometry_shape_matches(buffer(:, :, 1), reactant, 1.0e-4_wp)) then
      call env%error("the AutoNEB restart path does not start at the reactant", source)
      return
   end if
   if (.not. geometry_shape_matches(buffer(:, :, size(buffer, 3)), product, &
         & 1.0e-4_wp)) then
      call env%error("the AutoNEB restart path does not end at the product", source)
      return
   end if
   frames = buffer
   fail = .false.

end subroutine landscape_ts_read_path

!> True if two geometries agree up to a rigid body motion.
!>
!> The check compares **all** pairwise atom distances, which is invariant under
!> a translation and a rotation and is therefore the right measure for the
!> endpoints of a path read back from a file (the band kernel may have stored
!> them in a rotated frame).  Atom order and count are compared implicitly.
pure function geometry_shape_matches(xyz, reference, tol) result(matches)
   !> Geometry (bohr)
   real(wp), intent(in) :: xyz(:, :)
   !> Reference geometry (bohr)
   real(wp), intent(in) :: reference(:, :)
   !> Tolerance (bohr)
   real(wp), intent(in) :: tol
   !> True if every pair distance agrees within `tol`
   logical :: matches

   real(wp) :: d1, d2
   integer :: a, b, nat

   matches = .false.
   if (size(xyz, 1) /= 3 .or. size(reference, 1) /= 3) return
   nat = size(xyz, 2)
   if (size(reference, 2) /= nat) return
   do a = 1, nat
      do b = a + 1, nat
         d1 = sqrt(sum((xyz(:, a) - xyz(:, b))**2))
         d2 = sqrt(sum((reference(:, a) - reference(:, b))**2))
         if (abs(d1 - d2) > tol) return
      end do
   end do
   matches = .true.

end function geometry_shape_matches

!> Print the image table of a path.
!>
!> `label` is the method family of the report line (`band`, `autoneb`,
!> `string`); the table itself has the same layout for every chain-of-states
!> search, so the reader of a log can compare the methods directly.
subroutine write_path_table(env, label, energies, distances, highest, climbing)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Method family of the summary lines
   character(len=*), intent(in) :: label
   !> Energy of every image (hartree)
   real(wp), intent(in) :: energies(:)
   !> Distance between consecutive images (bohr)
   real(wp), intent(in) :: distances(:)
   !> Highest energy internal image
   integer, intent(in) :: highest
   !> Climbing image (0 if none)
   integer, intent(in) :: climbing

   integer :: k, n_images

   n_images = size(energies)
   write(env%unit, '(1x,a)') "   image        energy/Eh       dE/kcal/mol   distance/bohr"
   do k = 1, n_images
      if (k == 1) then
         write(env%unit, '(1x,i6,1x,f20.12,1x,f14.4,1x,a)') k, &
            & energies(k), (energies(k) - energies(1))*autokcal, "reactant"
      else if (k == n_images) then
         write(env%unit, '(1x,i6,1x,f20.12,1x,f14.4,1x,a)') k, &
            & energies(k), (energies(k) - energies(1))*autokcal, "product"
      else if (k - 1 <= size(distances)) then
         write(env%unit, '(1x,i6,1x,f20.12,1x,f14.4,1x,f14.6)') k, &
            & energies(k), (energies(k) - energies(1))*autokcal, distances(k - 1)
      else
         write(env%unit, '(1x,i6,1x,f20.12,1x,f14.4)') k, &
            & energies(k), (energies(k) - energies(1))*autokcal
      end if
   end do
   if (climbing > 0) then
      write(env%unit, '(1x,a,i0)') trim(label)//": climbing image (TS candidate): ", &
         & climbing
   else
      write(env%unit, '(1x,a,i0)') trim(label)//": highest energy image (TS candidate): ", &
         & highest
   end if

end subroutine write_path_table

end module xtb_landscape_tasks_ts
