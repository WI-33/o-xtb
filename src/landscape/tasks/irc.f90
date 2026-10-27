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

!> Intrinsic reaction coordinate task of the landscape suite (PLAN section
!> 9.7).
!>
!> The task level entry point of the `irc` run type.  The **input geometry is
!> the transition state** (the usual xtb workflow
!> `xtb xtbts.xyz --gfn2 --lirc`): the geometry is confirmed with the M3
!> frequency analysis, **exactly one imaginary mode** is required, and the
!> mass weighted imaginary mode defines the initial direction of the two
!> integrations.  A geometry with zero or two or more imaginary modes aborts
!> with a non-zero exit status (the H2 lesson of the M4 round: not a warning).
!>
!> Direction contract: the transition vector is oriented with the M4 phase
!> convention (the largest absolute component becomes positive), which fixes
!> "forward" deterministically; `--lirc-swap` exchanges the two directions and
!> the choice is recorded in the result and in `irc.log`.  The component and
!> its atom are printed, so the direction is reproducible from the output.
!>
!> Output files (all through the xtb file namespace):
!>
!> | file                | content |
!> |---------------------|---------|
!> | `irc_forward.xyz`   | forward points, transition state first |
!> | `irc_backward.xyz`  | backward points, transition state first |
!> | `irc_path.xyz`      | both directions merged by signed arc length, highest point marked `TS` |
!> | `irc.log`           | method, frozen `etemp`, phase, per point table, endpoint assignment |
!>
!> The two endpoints are assigned to a minimum by an optional M1 relaxation
!> (`endpoint_relax`): `rmsd_to_min_fwd/bwd` is the **Kabsch aligned** RMS
!> deviation between the endpoint and the relaxed minimum (bohr; the endpoint is
!> rigid-body aligned to the minimum with the mass free Kabsch rotation of
!> `types/coords.f90` before the distance is taken) and `reached_min_fwd/bwd`
!> records whether the path terminated because it entered the basin (termination
!> by signal (b) or (c)) **and** the relaxation converged **and** the aligned
!> RMSD is below `rmsd_end` (`--lirc-rmsd-end`).  The RMSD criterion cannot be
!> applied without a relaxed minimum; with `--lirc-no-endpoint-relax` the
!> endpoint is therefore assigned by the termination signal alone and `irc.log`
!> states that the criterion was skipped.  The relaxation runs on
!> the same engine, i.e. under the same restart policy and the same frozen
!> `etemp` as the integration, and under the four convergence thresholds of
!> `--lirc-level` (its own step budget is `--lirc-max-steps`, like the inner
!> relaxations of the `ts` and `scan` tasks).
!>
!> Units: step length and arc lengths in `bohr*amu^1/2` (mass weighted
!> Cartesian), geometries in bohr, energies in hartree, gradients in
!> hartree/bohr, `etemp` in kelvin.
module xtb_landscape_tasks_irc
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_symbols, only : toSymbol
   use xtb_mctc_filetypes, only : fileType
   use xtb_io_writer, only : writeMolecule
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule, init
   use xtb_type_restart, only : TRestart
   use xtb_type_calculator, only : TCalculator
   use xtb_setparam, only : get_namespace
   use xtb_landscape_settings, only : landscape_irc_settings, landscape_vib_settings, &
      & landscape_opt_settings, landscape_engine_control, landscape_irc_method_name, &
      & landscape_restart_name, landscape_hessian_update_name, landscape_lowercase, &
      & landscape_real_to_string, &
      & p_landscape_irc_gs, p_landscape_irc_lqa, p_landscape_irc_hpc, &
      & p_landscape_irc_eulerpc, landscape_log_failed, landscape_int_to_string
   use xtb_landscape_results, only : landscape_result, landscape_irc_result, &
      & landscape_vib_result, p_landscape_irc_term_grad, &
      & p_landscape_irc_term_energy_rise, p_landscape_irc_term_limit_cycle, &
      & p_landscape_irc_term_max_steps, &
      & landscape_irc_term_name
   use xtb_landscape_coords, only : geometry_pack, geometry_unpack, vec_rms, vec_norm, &
      & align_geometry
   use xtb_landscape_engine, only : landscape_engine, new_landscape_engine, &
      & landscape_hessian_provider, new_landscape_hessian_provider
   use xtb_landscape_tasks_opt, only : landscape_engine_model
   use xtb_landscape_optimize, only : landscape_optimize_model
   use xtb_landscape_constraints, only : landscape_frozen_mask
   use xtb_landscape_vib_modes, only : landscape_vib_analysis
   use xtb_landscape_saddle_eigenfollow, only : landscape_transition_mode
   use xtb_landscape_irc_corrector, only : landscape_irc_masses, &
      & landscape_irc_hessian_to_weighted, landscape_irc_sqrt_mass, &
      & landscape_irc_select_direction, landscape_irc_path
   use xtb_landscape_irc_gs, only : landscape_irc_gs_walk
   use xtb_landscape_irc_lqa, only : landscape_irc_lqa_walk
   use xtb_landscape_irc_hpc, only : landscape_irc_hpc_walk
   use xtb_landscape_irc_eulerpc, only : landscape_irc_eulerpc_walk
   use xtb_landscape_json, only : landscape_json_document, &
      & landscape_json_document_init, landscape_json_add_counter, &
      & landscape_json_add_core, &
      & landscape_json_add_value, landscape_json_add_flag, &
      & landscape_json_add_series, landscape_json_write, &
      & landscape_json_status_ok, landscape_json_status_incomplete, &
      & landscape_json_status_error
   implicit none(type, external)
   private

   public :: landscape_irc
   !> Endpoint attribution rule of this task, public so that its criteria can be
   !> unit tested (`test/unit/test_landscape_irc.f90`)
   public :: endpoint_reached

contains

!> Task level entry point: follow the intrinsic reaction coordinate of a
!> transition state in both directions.
subroutine landscape_irc(env, mol, chk, calc, settings, result)
   character(len=*), parameter :: source = 'landscape_irc'
   !> Computational environment
   type(TEnvironment), intent(inout), target :: env
   !> Molecular structure (the transition state on entry)
   type(TMolecule), intent(inout), target :: mol
   !> Wavefunction data
   type(TRestart), intent(inout), target :: chk
   !> Single point calculator
   class(TCalculator), intent(in) :: calc
   !> IRC settings
   type(landscape_irc_settings), intent(in) :: settings
   !> Task result
   type(landscape_irc_result), intent(inout) :: result

   type(landscape_engine), target :: engine
   type(landscape_engine_model) :: model
   class(landscape_hessian_provider), allocatable :: hp
   type(landscape_engine_control) :: ctrl
   type(landscape_transition_mode) :: mode
   type(landscape_irc_path) :: path_fwd, path_bwd
   real(wp), allocatable :: hess(:, :), hess_work(:, :), mass(:), sqrt_mass(:)
   real(wp), allocatable :: hq(:, :), x_ts(:), v(:), x_end(:), x_min(:), g_ts(:)
   real(wp) :: energy, gnorm, rmsd_fwd, rmsd_bwd
   real(wp) :: e_min_fwd, e_min_bwd
   logical :: fail, ok, relaxed_fwd, relaxed_bwd, have_bwd
   integer :: n_scf_fail, nlog, isign_fwd, isign_bwd, atom, cart
   integer :: n_relax_fwd, n_relax_bwd
   character(len=:), allocatable :: unused_hess

   fail = .false.
   result%method = settings%method
   result%step_length = settings%step_length
   result%swapped = settings%swap
   result%n_imag = -1

   ! --- preconditions ----------------------------------------------------
   call check_irc_settings(env, settings, ok)
   if (.not. ok) then
      result%message = "invalid IRC settings"
      return
   end if

   if (settings%printlevel > 0) then
      write(env%unit, '(1x,a,a,a)') "landscape: intrinsic reaction coordinate (", &
         & landscape_irc_method_name(settings%method), ")"
      write(env%unit, '(1x,a,es13.6,a)') "step length:              ", &
         & settings%step_length, " bohr*amu^1/2"
      write(env%unit, '(1x,a,i0)') "maximum number of steps:  ", settings%max_steps
      ! R5 of CLI redesign section 3.5: `e_tol` follows the convergence level
      ! (`set_conv_level`) and reaches the endpoint relaxations through the
      ! inner minimisation settings (`opt%e_tol = settings%e_tol`).
      write(env%unit, '(1x,a,es15.6,a)') "energy tolerance:         ", &
         & settings%e_tol, " Eh"
      write(env%unit, '(1x,a,a)') "Hessian update:           ", &
         & landscape_hessian_update_name(settings%hessian_update)
      write(env%unit, '(1x,a,a)') "restart policy:           ", &
         & landscape_restart_name(settings%engine%restart_mode)
      ! The effective Hessian source, printed like the identical line of the
      ! scan, opt, ts and vib tasks.  It is an effective value, not a dangerous
      ! default of section 8.3, so it is printed regardless of whether
      ! `--lirc-hessian` named it (the default `numerical` is part of the
      ! answer: it decides whether the initial Hessian and the recalculations
      ! are displaced single points or the Lindh-D2 model).
      write(env%unit, '(1x,a,a)') "Hessian source:           ", &
         & trim(settings%hessian_source)
      ! CLI redesign sections 7.3 and 9.2: every flag that carries a unit must
      ! echo its effective value and unit, and a knob that no code path of the
      ! selected method reads must say so.  Both values come from the settings
      ! object the engine is built from, not from the command line text.  The
      ! predicate is the measured one of `35-stage2b-fix.md` section T1: every
      ! integrator builds the initial Hessian of the reaction path from the
      ! provider, and the provider ignores the finite difference step for
      ! every other Hessian source.
      write(env%unit, '(1x,a,es15.6,a)') "Hessian step:             ", &
         & settings%engine%hess_step, " bohr"
      write(env%unit, '(1x,a,es15.6)') "Hessian SCC accuracy:     ", &
         & settings%engine%hess_acc
      unused_hess = irc_unused_hessian_step(settings)
      if (len(unused_hess) > 0) then
         write(env%unit, '(1x,a,a,a)') "(not used: ", unused_hess, ")"
      end if
      write(env%unit, '(1x,a,f12.3,a)') "electronic temperature:   ", &
         & settings%etemp, " K (frozen for the whole run)"
      ! D2 of the CLI redesign section 8.3 asked for the active default
      ! integrator to be echoed, because `hpc` is not the `gs` of the older
      ! documentation (behaviour change of the M6 validation round).  That echo
      ! is the startup engine banner now
      ! (`xtb_prog_main::landscape_engine_banner`): it names the integrator of
      ! every reaction path run and says `(default)`/`(explicit)`.
   end if

   ! H2 of report 220 section 4: an option of this run type that only a subset
   ! of the integrators reads must not disappear without a word when the user
   ! selected an integrator outside that subset (see
   ! `irc_report_ignored_options`).
   call irc_report_ignored_options(env, settings)

   ! `$fix` is not implemented for the reaction path: the path advance happens
   ! in the mass weighted metric (`bohr*amu^1/2`), and the arc length column
   ! `s` is a contract of this task (M6: `s` is the mass weighted arc length,
   ! `docs/xtb-landscape-PLAN.md` section 5.3.7).  Restricting a step to the
   ! `$fix` subspace shortens the step actually taken, so `s` can no longer be
   ! accumulated from the *requested* step - and the LQA corrector iterates on
   ! a requested arc length (`Sigma_n` reaches `s_target` with a Newton step)
   ! which a projected step cannot reach, so its fixed point is not defined.
   ! Every method is refused with `$fix` rather than run in a metric that is
   ! neither the contract of this task nor the unconstrained one
   ! (docs/_simplify-recon/04-o2-and-silent-defects-fix.md section 5.3;
   ! PLAN section 6.9 open item).
   if (any(landscape_frozen_mask(mol%n))) then
      call env%error("the reaction path (--lirc) does not support $fix: its "//&
         & "path advance is defined in the mass weighted arc length, which a "//&
         & "fixed coordinate would have to redefine.  Use --lopt or --lts "//&
         & "(their searches project the step), or remove the $fix block", source)
      result%message = "the reaction path does not support $fix"
      return
   end if

   ! --- engine -----------------------------------------------------------
   ctrl = settings%engine
   ctrl%etemp = settings%etemp
   call new_landscape_engine(engine, ctrl)
   call new_landscape_hessian_provider(env, settings%hessian_source, hp, &
      & sccacc=settings%engine%hess_acc, filename=settings%hessian_file, &
      & printlevel=settings%printlevel, jobs=settings%engine%hess_jobs)
   if (.not. allocated(hp)) then
      result%message = "IRC failed"
      call engine%destroy()
      return
   end if
   call engine%set_hessian_provider(hp)
   call engine%set_calculator(calc)
   model%engine => engine
   model%env => env
   model%mol => mol
   model%chk => chk

   x_ts = geometry_pack(mol%xyz)
   ! --- transition state confirmation (exactly one imaginary mode) -------
   allocate(hess(3*mol%n, 3*mol%n), source=0.0_wp)
   call engine%hessian(env, mol, chk, hess, fail)
   if (fail) then
      result%message = "the transition state Hessian could not be built"
      call engine%destroy()
      return
   end if
   allocate(hess_work(3*mol%n, 3*mol%n), source=hess)
   call confirm_ts(env, mol, hess_work, settings, result)
   if (result%n_imag /= 1) then
      call engine%report(env%unit)
      call engine%destroy()
      return
   end if

   ! --- initial direction (mass weighted transition vector) --------------
   mass = landscape_irc_masses(mol)
   sqrt_mass = landscape_irc_sqrt_mass(mass)
   allocate(hq(3*mol%n, 3*mol%n), source=0.0_wp)
   call landscape_irc_hessian_to_weighted(sqrt_mass, hess, hq)
   call landscape_irc_select_direction(env, hq, settings%target_mode, mode, fail)
   if (fail) then
      result%message = "the initial direction could not be determined"
      call engine%destroy()
      return
   end if
   v = mode%vector
   result%phase_index = mode%phase_index
   result%phase_value = mode%phase_value
   if (settings%printlevel > 0) then
      atom = (mode%phase_index - 1)/3 + 1
      cart = mod(mode%phase_index - 1, 3) + 1
      write(env%unit, '(1x,a,i0,a,es13.6)') "transition mode index:     ", &
         & mode%index, " (eigenvalue ", mode%eigenvalue
      write(env%unit, '(1x,a)') "forward is defined by the largest absolute component "//&
         & "of the mass weighted transition vector"
      write(env%unit, '(1x,a,f12.6,a,i0,a,i0,a)') "phase component:          ", &
         & mode%phase_value, " (cartesian ", cart, " of atom ", atom, ")"
      if (settings%swap) then
         write(env%unit, '(1x,a)') "the two directions are exchanged (--lirc-swap)"
      end if
   end if
   isign_fwd = 1
   isign_bwd = -1
   if (settings%swap) then
      isign_fwd = -1
      isign_bwd = 1
   end if

   ! --- the two integrations --------------------------------------------
   call walk_direction(env, model, settings, mass, x_ts, hess, v, &
      & real(isign_fwd, wp), path_fwd, fail)
   if (fail) then
      result%message = "the forward IRC integration failed"
      call engine%destroy()
      return
   end if
   have_bwd = settings%both_directions
   if (have_bwd) then
      call walk_direction(env, model, settings, mass, x_ts, hess, v, &
         & real(isign_bwd, wp), path_bwd, fail)
      if (fail) then
         result%message = "the backward IRC integration failed"
         call engine%destroy()
         return
      end if
   end if

   ! --- endpoint assignment ---------------------------------------------
   rmsd_fwd = 0.0_wp
   rmsd_bwd = 0.0_wp
   e_min_fwd = 0.0_wp
   e_min_bwd = 0.0_wp
   n_relax_fwd = 0
   n_relax_bwd = 0
   relaxed_fwd = .false.
   relaxed_bwd = .false.
   if (settings%endpoint_relax .and. path_fwd%n > 0) then
      call relax_endpoint(env, model, path_fwd%xyz(:, path_fwd%n), settings, &
         & x_min, rmsd_fwd, e_min_fwd, n_relax_fwd, relaxed_fwd, fail)
      if (fail) then
         result%message = "the forward endpoint could not be relaxed"
         call engine%destroy()
         return
      end if
   end if
   if (settings%endpoint_relax .and. have_bwd .and. path_bwd%n > 0) then
      call relax_endpoint(env, model, path_bwd%xyz(:, path_bwd%n), settings, &
         & x_min, rmsd_bwd, e_min_bwd, n_relax_bwd, relaxed_bwd, fail)
      if (fail) then
         result%message = "the backward endpoint could not be relaxed"
         call engine%destroy()
         return
      end if
   end if

   ! --- fill the result --------------------------------------------------
   call fill_path(path_fwd, mol%n, result%xyz_fwd, result%e_fwd, result%s_fwd, &
      & result%gnorm_fwd, result%n_eval_fwd)
   result%n_fwd = path_fwd%n
   result%term_fwd = path_fwd%term
   if (have_bwd) then
      call fill_path(path_bwd, mol%n, result%xyz_bwd, result%e_bwd, result%s_bwd, &
         & result%gnorm_bwd, result%n_eval_bwd)
      result%n_bwd = path_bwd%n
      result%term_bwd = path_bwd%term
   end if

   allocate(g_ts(3*mol%n), source=0.0_wp)
   call model%evaluate(x_ts, energy, g_ts, fail)
   if (fail) then
      result%message = "the transition state could not be evaluated"
      call engine%destroy()
      return
   end if
   result%energy_ts = energy
   result%gnorm_ts = vec_norm(g_ts)
   result%rmsd_to_min_fwd = rmsd_fwd
   result%rmsd_to_min_bwd = rmsd_bwd
   result%e_min_fwd = e_min_fwd
   result%e_min_bwd = e_min_bwd
   result%n_relax_fwd = n_relax_fwd
   result%n_relax_bwd = n_relax_bwd
   result%n_pinned_fwd = path_fwd%n_pinned
   result%n_pinned_bwd = 0
   if (have_bwd) result%n_pinned_bwd = path_bwd%n_pinned
   result%reached_min_fwd = endpoint_reached(path_fwd%term, settings%endpoint_relax, &
      & relaxed_fwd, rmsd_fwd, settings%rmsd_end)
   result%reached_min_bwd = .false.
   if (have_bwd) result%reached_min_bwd = endpoint_reached(path_bwd%term, &
      & settings%endpoint_relax, relaxed_bwd, rmsd_bwd, settings%rmsd_end)
   result%etemp = engine%effective_etemp()

   if (settings%write_files) then
      nlog = env%nLog
      call write_irc_files(env, mol, settings, result, x_ts)
      if (landscape_log_failed(env, nlog)) then
         result%message = "the IRC output files could not be written"
         call engine%destroy()
         return
      end if
   end if

   if (settings%printlevel > 0) then
      call result%write_log(env%unit)
   end if

   n_scf_fail = engine%n_scf_fail
   call engine%report(env%unit)
   call engine%destroy()

   ! restore the transition state geometry: the model evaluation moved
   ! `mol%xyz` to the last evaluated point
   call geometry_unpack(x_ts, mol%xyz)

   if (n_scf_fail > 0) then
      call env%error("IRC aborted: the SCC did not converge at least once", source)
      result%message = "IRC aborted: unconverged SCC"
      return
   end if

   if (result%term_fwd == 1 .or. (have_bwd .and. result%term_bwd == 1)) then
      ! signal (a): at least one direction was cut off by the step budget
      call env%warning("an IRC direction was stopped by the step budget "//&
         & "before reaching a minimum", source)
   end if
   ! signal (d): the GS limit cycle guard truncated a stalled walk.  A warning
   ! is required (not optional): the endpoint is not a converged path point,
   ! it is one phase of the oscillation, so `reached_min_fwd` below is false
   ! unless an endpoint relaxation confirms the geometry independently.
   if (result%term_fwd == p_landscape_irc_term_limit_cycle) then
      call env%warning("the forward endpoint is one phase of the oscillation, "//&
         & "not a converged path point (the GS macro step stalled in a limit "//&
         & "cycle)", source)
   end if
   if (have_bwd .and. result%term_bwd == p_landscape_irc_term_limit_cycle) then
      call env%warning("the backward endpoint is one phase of the oscillation, "//&
         & "not a converged path point (the GS macro step stalled in a limit "//&
         & "cycle)", source)
   end if
   ! the message is the human readable verdict of the run; a truncated path
   ! must not claim completion (it reaches `--ljson` and `irc.log` through
   ! `landscape_irc_result_write`)
   if (result%term_fwd == p_landscape_irc_term_limit_cycle .or. &
         & (have_bwd .and. result%term_bwd == p_landscape_irc_term_limit_cycle)) then
      result%message = "IRC truncated by the limit cycle guard: the endpoint "//&
         & "is one phase of the oscillation, not a converged path point"
   else
      result%message = "IRC integration completed"
   end if

   if (settings%ljson) then
      call write_json(env, mol, settings, result, have_bwd, n_scf_fail)
   end if

end subroutine landscape_irc

!> Write the structured result `landscape_irc.json`.
!>
!> The `s` series is the arc length of the integrator in
!> `bohr*amu^1/2`; as in the printed table it is a chord summation for
!> `gs`/`hpc`/`eulerpc` and the integrated arc length for `lqa`, so the
!> series must not be compared across methods (M6 handover LOW-4).  The
!> energy series are the reference for the endpoint attribution: the last
!> point of each direction is the one the `reached_min_*` flag refers to.
!>
!> Status: an IRC of a first order saddle that produced both paths is `ok`
!> (or `incomplete` when a direction only stopped on the step budget).  A run
!> whose Gonzalez-Schlegel walk was truncated by the limit cycle guard
!> (`--lirc-osc-window`) is `incomplete` as well and `converged` is false: the
!> path was cut, the endpoint is one phase of the oscillation, so reporting
!> the run as a completed integration would be wrong.  A
!> run whose input is not a saddle, that could not build a direction or that
!> had an unconverged SCC is `error` - the task reports it through `env` and
!> the process exits non-zero.
subroutine write_json(env, mol, settings, result, have_bwd, n_scf_fail)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecule of the run (the core block's geometry, the transition state)
   type(TMolecule), intent(in) :: mol
   !> IRC settings
   type(landscape_irc_settings), intent(in) :: settings
   !> IRC result
   type(landscape_irc_result), intent(in) :: result
   !> True if the backward direction was integrated
   logical, intent(in) :: have_bwd
   !> Number of unconverged SCCs of the run
   integer, intent(in) :: n_scf_fail

   type(landscape_json_document) :: doc = landscape_json_document()
   logical :: ok, complete, truncated
   integer :: status, n_eval_fwd, n_eval_bwd

   n_eval_fwd = 0
   if (allocated(result%n_eval_fwd)) n_eval_fwd = sum(result%n_eval_fwd)
   n_eval_bwd = 0
   if (allocated(result%n_eval_bwd)) n_eval_bwd = sum(result%n_eval_bwd)

   ! the limit cycle guard cuts a stalled walk: the endpoint is not a
   ! converged path point, so the run is not `ok` and not `converged`
   truncated = result%term_fwd == p_landscape_irc_term_limit_cycle
   if (have_bwd) then
      truncated = truncated .or. result%term_bwd == p_landscape_irc_term_limit_cycle
   end if
   complete = result%is_saddle .and. (result%term_fwd /= 1) .and. .not. truncated
   if (have_bwd) complete = complete .and. (result%term_bwd /= 1)
   status = landscape_json_status_ok
   if (complete .and. result%n_fwd > 0) status = landscape_json_status_ok
   if (.not. complete) status = landscape_json_status_incomplete
   if (.not. result%is_saddle .or. result%n_fwd == 0) then
      status = landscape_json_status_error
   end if
   if (n_scf_fail > 0) status = landscape_json_status_error

   call landscape_json_document_init(doc, 'irc')
   doc%status = status
   ! a direction that stopped on the step budget is as unconverged as a
   ! truncated one: `status` already reports both as `incomplete`, and the
   ! documented meaning of this key is the verdict of the run, so it must not
   ! claim `converged` for a walk that ran out of budget (measured: a
   ! `--lirc-max-steps` run reported `status = incomplete` with
   ! `converged = true` before this line was extended)
   doc%converged = result%is_saddle .and. result%n_fwd > 0 .and. .not. truncated &
      & .and. result%term_fwd /= p_landscape_irc_term_max_steps
   if (allocated(result%message)) doc%message = result%message
   call landscape_json_add_core(doc, mol, result%energy_ts, settings%input_file)
   call landscape_json_add_counter(doc, 'method', settings%method)
   call landscape_json_add_counter(doc, 'n_fwd', result%n_fwd)
   call landscape_json_add_counter(doc, 'n_bwd', result%n_bwd)
   call landscape_json_add_counter(doc, 'n_imag', result%n_imag)
   call landscape_json_add_counter(doc, 'n_eval_fwd', n_eval_fwd)
   call landscape_json_add_counter(doc, 'n_eval_bwd', n_eval_bwd)
   call landscape_json_add_counter(doc, 'n_relax_fwd', result%n_relax_fwd)
   call landscape_json_add_counter(doc, 'n_relax_bwd', result%n_relax_bwd)
   call landscape_json_add_counter(doc, 'n_pinned_fwd', result%n_pinned_fwd)
   call landscape_json_add_counter(doc, 'n_pinned_bwd', result%n_pinned_bwd)
   call landscape_json_add_counter(doc, 'n_scf_fail', n_scf_fail)
   call landscape_json_add_counter(doc, 'term_fwd', result%term_fwd)
   call landscape_json_add_counter(doc, 'term_bwd', result%term_bwd)
   call landscape_json_add_value(doc, 'energy_ts', result%energy_ts)
   call landscape_json_add_value(doc, 'gnorm_ts', result%gnorm_ts)
   call landscape_json_add_value(doc, 'step_length', result%step_length)
   call landscape_json_add_value(doc, 'etemp', result%etemp)
   call landscape_json_add_value(doc, 'e_min_fwd', result%e_min_fwd)
   call landscape_json_add_value(doc, 'e_min_bwd', result%e_min_bwd)
   call landscape_json_add_value(doc, 'rmsd_to_min_fwd', result%rmsd_to_min_fwd)
   call landscape_json_add_value(doc, 'rmsd_to_min_bwd', result%rmsd_to_min_bwd)
   call landscape_json_add_flag(doc, 'is_saddle', result%is_saddle)
   call landscape_json_add_flag(doc, 'reached_min_fwd', result%reached_min_fwd)
   call landscape_json_add_flag(doc, 'reached_min_bwd', result%reached_min_bwd)
   call landscape_json_add_flag(doc, 'two_way', have_bwd)
   call landscape_json_add_flag(doc, 'swapped', result%swapped)
   if (allocated(result%e_fwd)) then
      call landscape_json_add_series(doc, 'energy_fwd', result%e_fwd)
      call landscape_json_add_series(doc, 's_fwd', result%s_fwd)
      call landscape_json_add_series(doc, 'gnorm_fwd', result%gnorm_fwd)
   end if
   if (allocated(result%e_bwd) .and. result%n_bwd > 0) then
      call landscape_json_add_series(doc, 'energy_bwd', result%e_bwd)
      call landscape_json_add_series(doc, 's_bwd', result%s_bwd)
      call landscape_json_add_series(doc, 'gnorm_bwd', result%gnorm_bwd)
   end if
   call landscape_json_write(env, doc, ok)

end subroutine write_json

!> Dispatch one direction to the integrator selected by the settings.
subroutine walk_direction(env, model, settings, mass, x0, hess, v, sign, path, fail)
   character(len=*), parameter :: source = 'landscape_irc_walk_direction'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Potential
   class(landscape_engine_model), intent(inout) :: model
   !> IRC settings
   type(landscape_irc_settings), intent(in) :: settings
   !> Atomic masses in amu (length 3N)
   real(wp), intent(in) :: mass(:)
   !> Transition state geometry, flat Cartesian (bohr)
   real(wp), intent(in) :: x0(:)
   !> Cartesian Hessian at the transition state (hartree/bohr^2)
   real(wp), intent(in) :: hess(:, :)
   !> Mass weighted transition vector (unit, phase fixed)
   real(wp), intent(in) :: v(:)
   !> Direction sign
   real(wp), intent(in) :: sign
   !> Resulting path
   type(landscape_irc_path), intent(inout) :: path
   !> Failure flag
   logical, intent(out) :: fail

   select case (settings%method)
   case (p_landscape_irc_gs)
      call landscape_irc_gs_walk(env, model, settings, mass, x0, hess, v, sign, &
         & path, fail)
   case (p_landscape_irc_lqa)
      call landscape_irc_lqa_walk(env, model, settings, mass, x0, hess, v, sign, &
         & path, fail)
   case (p_landscape_irc_hpc)
      call landscape_irc_hpc_walk(env, model, settings, mass, x0, hess, v, sign, &
         & path, fail)
   case (p_landscape_irc_eulerpc)
      call landscape_irc_eulerpc_walk(env, model, settings, mass, x0, hess, v, &
         & sign, path, fail)
   case default
      call env%error(trim(landscape_irc_method_name(settings%method))//&
         & " is not implemented yet", source)
      fail = .true.
   end select

end subroutine walk_direction

!> Confirm the transition state with the M3 frequency analysis.
!>
!> Exactly one imaginary mode is required; the three outcomes (first order
!> saddle, minimum, higher order saddle) are printed and a non-saddle is an
!> error (a non-zero exit status), not a warning.  `hess` is the Cartesian
!> Hessian of the transition state; the analysis projects it in place, so the
!> caller has to keep a pristine copy for the direction selection.
subroutine confirm_ts(env, mol, hess, settings, result)
   character(len=*), parameter :: source = 'landscape_irc_confirm_ts'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure (transition state)
   type(TMolecule), intent(in) :: mol
   !> Cartesian Hessian (hartree/bohr^2), projected in place
   real(wp), intent(inout) :: hess(:, :)
   !> IRC settings
   type(landscape_irc_settings), intent(in) :: settings
   !> IRC result, updated with the imaginary mode count
   type(landscape_irc_result), intent(inout) :: result

   type(landscape_vib_settings) :: vib_settings
   type(landscape_vib_result) :: vib
   logical :: fail

   vib_settings%project_tr = .true.
   vib_settings%imag_cutoff = -20.0_wp
   call landscape_vib_analysis(env, mol, hess, vib_settings, vib, fail)
   if (fail) then
      call env%error("the transition state confirmation failed", source)
      result%message = "the transition state confirmation failed"
      return
   end if
   result%n_imag = vib%n_imag
   result%is_saddle = (vib%n_imag == 1)

   if (settings%printlevel > 0) then
      write(env%unit, '(1x,a,i0)') "transition state imaginary modes: ", result%n_imag
      select case (vib%n_imag)
      case (1)
         write(env%unit, '(1x,a)') "confirmation: first order saddle point"
      case (0)
         write(env%unit, '(1x,a)') "confirmation: minimum (0 imaginary modes) - "//&
            & "not a transition state"
      case default
         write(env%unit, '(1x,a,i0,a)') "confirmation: higher order saddle (point "//&
            & "with ", vib%n_imag, " imaginary modes) - not a transition state"
      end select
   end if

   if (vib%n_imag /= 1) then
      if (vib%n_imag == 0) then
         result%message = "the input geometry is a minimum (0 imaginary modes), "//&
            & "not a transition state"
      else
         result%message = "the input geometry is a higher order saddle ("//&
            & trim(landscape_int_to_string(vib%n_imag))//" imaginary modes), not a "//&
            & "first order saddle point"
      end if
      call env%error(trim(result%message), source)
   end if

end subroutine confirm_ts

!> Report every IRC option the selected integrator does not read.
!>
!> The counterpart of `xtb_landscape_tasks_ts::ts_report_ignored_options` for
!> the reaction path (report 220, finding H2).  The four integrators of this
!> task share one settings type but not one set of knobs: the sphere
!> relaxation and the limit cycle guard belong to `gs`, the Euler substeps and
!> the surface fit to `lqa`/`hpc`, the modified Bulirsch-Stoer corrector and
!> the Hessian update to `hpc`/`eulerpc`, and the coarse phase of the
!> predictor-corrector to `eulerpc` alone.  A knob that no code path of the
!> selected integrator reads was accepted and then silently discarded, which
!> invalidates a study without a trace; each one is reported once and the run
!> continues (`env%warning` is the fatal error under `--strict` here, as
!> everywhere in the suite).
!>
!> The predicate is the effective value, compared against a default constructed
!> settings object, exactly like its transition state counterpart: the
!> settings record no per-option "was given" flag, and repeating a default
!> could not have changed anything.
!>
!> Not covered, on purpose: `--lirc-hessian-file` is reported by the provider
!> factory for every task, `--lirc-sccacc`/`--lirc-hess-step` by
!> `landscape_engine%report_unused_hess_step`, and `--lirc-step` outside the
!> GS window by `check_irc_settings`.
subroutine irc_report_ignored_options(env, settings)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> IRC settings of this run
   type(landscape_irc_settings), intent(in) :: settings

   type(landscape_irc_settings) :: dflt
   integer :: method

   method = settings%method

   ! --- the GS sphere relaxation and its guards --------------------------
   if (settings%max_micro_cycles /= dflt%max_micro_cycles) then
      if (.not. irc_method_reads(method, [p_landscape_irc_gs])) &
         & call report_ignored_irc_option(env, '--lirc-micro-cycles', method, &
         & 'the gs macro step only')
   end if
   if (settings%micro_step_thresh /= dflt%micro_step_thresh) then
      if (.not. irc_method_reads(method, [p_landscape_irc_gs])) &
         & call report_ignored_irc_option(env, '--lirc-micro-thresh', method, &
         & 'the gs macro step only')
   end if
   if (settings%osc_window /= dflt%osc_window) then
      if (.not. irc_method_reads(method, [p_landscape_irc_gs])) &
         & call report_ignored_irc_option(env, '--lirc-osc-window', method, &
         & 'the gs macro step only')
   end if
   if (settings%adaptive .neqv. dflt%adaptive) then
      if (.not. irc_method_reads(method, [p_landscape_irc_gs])) &
         & call report_ignored_irc_option(env, '--lirc-adaptive', method, &
         & 'the gs macro step only')
   end if

   ! --- the analytic propagation of lqa and hpc --------------------------
   if (settings%euler_n /= dflt%euler_n) then
      if (.not. irc_method_reads(method, [p_landscape_irc_lqa, p_landscape_irc_hpc])) &
         & call report_ignored_irc_option(env, '--lirc-euler-n', method, &
         & 'the lqa and hpc integrators only')
   end if

   ! --- the predictor-corrector ------------------------------------------
   if (settings%loose_cycles /= dflt%loose_cycles) then
      if (.not. irc_method_reads(method, [p_landscape_irc_eulerpc])) &
         & call report_ignored_irc_option(env, '--lirc-loose-cycles', method, &
         & 'the eulerpc integrator only')
   end if
   if (settings%dwi_n /= dflt%dwi_n) then
      if (.not. irc_method_reads(method, [p_landscape_irc_hpc, p_landscape_irc_eulerpc])) &
         & call report_ignored_irc_option(env, '--lirc-dwi-n', method, &
         & 'the hpc and eulerpc integrators only')
   end if
   if (settings%mbs_max_k /= dflt%mbs_max_k) then
      if (.not. irc_method_reads(method, [p_landscape_irc_hpc, p_landscape_irc_eulerpc])) &
         & call report_ignored_irc_option(env, '--lirc-mbs-max-k', method, &
         & 'the hpc and eulerpc integrators only')
   end if
   if (settings%mbs_points /= dflt%mbs_points) then
      if (.not. irc_method_reads(method, [p_landscape_irc_hpc, p_landscape_irc_eulerpc])) &
         & call report_ignored_irc_option(env, '--lirc-mbs-points', method, &
         & 'the hpc and eulerpc integrators only')
   end if
   if (settings%mbs_tol /= dflt%mbs_tol) then
      if (.not. irc_method_reads(method, [p_landscape_irc_hpc, p_landscape_irc_eulerpc])) &
         & call report_ignored_irc_option(env, '--lirc-mbs-tol', method, &
         & 'the hpc and eulerpc integrators only')
   end if
   if (settings%hessian_recalc /= dflt%hessian_recalc) then
      if (.not. irc_method_reads(method, [p_landscape_irc_hpc, p_landscape_irc_eulerpc])) &
         & call report_ignored_irc_option(env, '--lirc-recalc', method, &
         & 'the hpc and eulerpc integrators only')
   end if

end subroutine irc_report_ignored_options

!> True if `method` is one of the `readers`.
pure logical function irc_method_reads(method, readers)
   !> Selected IRC integrator
   integer, intent(in) :: method
   !> Integrators that read the option
   integer, intent(in) :: readers(:)

   irc_method_reads = any(readers == method)

end function irc_method_reads

!> Report one accepted IRC option the selected integrator ignores.
subroutine report_ignored_irc_option(env, flag, method, readers)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Spelling of the option
   character(len=*), intent(in) :: flag
   !> Selected integrator
   integer, intent(in) :: method
   !> Integrators that read the option, as prose
   character(len=*), intent(in) :: readers

   call env%warning(flag//" is not read by --lirc-method "// &
      & trim(landscape_irc_method_name(method))//" (it configures "//readers// &
      & "); the value has no effect on this run", 'landscape_irc')

end subroutine report_ignored_irc_option

!> Validate the IRC settings the task itself depends on.
subroutine check_irc_settings(env, settings, ok)
   character(len=*), parameter :: source = 'landscape_irc'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> IRC settings
   type(landscape_irc_settings), intent(in) :: settings
   !> True if every parameter is usable
   logical, intent(out) :: ok

   character(len=400) :: buffer

   ok = .true.
   if (settings%step_length <= 0.0_wp) then
      call env%error("the IRC step length must be positive", source)
      ok = .false.
   else if (settings%method == p_landscape_irc_gs) then
      ! A4 of report 220 section 4: the positivity check is the only one this
      ! value gets, although the GS macro step has a *measured* window.  The
      ! value is not refused - a user may have a reason - but the measured
      ! behaviour of the two ends of the window is stated (report 67: from
      ! 0.20 on the sphere relaxation stalls in a cancelling two cycle; at
      ! 0.40 the endpoint still sits 0.236 bohr from the nearest minimum, and
      ! the HPC fallback reaches a third well), so a study that used a larger
      ! step knows what it bought.  `hpc`/`eulerpc` are correct at 0.20 and
      ! are not reported.
      if (settings%step_length < 5.0e-2_wp .or. settings%step_length > 1.5e-1_wp) then
         buffer = "--lirc-step "// &
            & trim(adjustl(landscape_real_to_string(settings%step_length, '(f9.3)')))// &
            & " is outside the 0.05-0.15 bohr*amu^1/2 window measured for "//&
            & "--lirc-method gs: the walk stalls in a two cycle from 0.20 on "//&
            & "and a 0.40 step stops 0.236 bohr from the nearest minimum (the "//&
            & "HPC fallback reaches a third well); the requested value is used"
         call env%warning(trim(buffer), source)
      end if
   end if
   if (settings%euler_n < 2) then
      call env%error("the IRC predictor needs at least two Euler substeps", source)
      ok = .false.
   end if
   if (settings%max_micro_cycles < 1) then
      call env%error("the IRC sphere relaxation needs at least one micro cycle", source)
      ok = .false.
   end if
   if (settings%grad_tol < 0.0_wp) then
      call env%error("the IRC gradient threshold must not be negative", source)
      ok = .false.
   end if
   if (settings%e_rise_steps < 1) then
      call env%error("the IRC energy rise criterion needs at least one step", source)
      ok = .false.
   end if
   if (settings%osc_window < 0) then
      call env%error("the IRC limit cycle window must not be negative "//&
         & "(0 disables the guard)", source)
      ok = .false.
   end if
   if (settings%dwi_n < 2) then
      call env%error("the IRC DWI surface needs at least two sample points "//&
         & "(the current and the predicted point of a step)", source)
      ok = .false.
   end if
   if (settings%rmsd_end < 0.0_wp) then
      call env%error("the IRC endpoint RMSD threshold must not be negative", source)
      ok = .false.
   end if

end subroutine check_irc_settings

!> Relax one endpoint with the M1 minimiser and measure its RMSD.
!>
!> The reported `rmsd` is the **Kabsch aligned** RMS deviation between the
!> endpoint and the relaxed minimum: the endpoint geometry is translated to its
!> arithmetic centroid and rotated by the proper Kabsch rotation that best
!> reproduces the minimum (`align_geometry` / `kabsch_align` of
!> `types/coords.f90`, the mass free alignment of the M5 band), then the RMS
!> distance is taken.  Without the alignment a rigid body drift of the endpoint
!> would inflate the RMSD; the aligned value answers the assigned question, "is
!> this the same minimum?", and is reported as aligned in `irc.log`.
!>
!> The relaxation uses the **same engine model** as the integration, so the
!> restart policy and the frozen `etemp` of the IRC run are inherited (M6 hard
!> constraint 5: one restart policy per run), and it inherits the four
!> convergence thresholds of `--lirc-level` (`set_conv_level`) like the inner
!> relaxation of `tasks/ts.f90:1069-1075`: without the copy the level flag
!> would be parsed and then never read (MEDIUM-1 of
!> `45-m6-independent-review.md`).  `e_min` and `n_relax` report the relaxed
!> minimum energy and the macro steps the relaxation needed, which is what
!> makes the level flag observable.
subroutine relax_endpoint(env, model, x_end, settings, x_min, rmsd, e_min, &
      & n_relax, converged, fail)
   character(len=*), parameter :: source = 'landscape_irc_relax_endpoint'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Potential
   class(landscape_engine_model), intent(inout) :: model
   !> Endpoint geometry, flat Cartesian (bohr)
   real(wp), intent(in) :: x_end(:)
   !> IRC settings
   type(landscape_irc_settings), intent(in) :: settings
   !> Relaxed minimum (flat Cartesian, bohr)
   real(wp), allocatable, intent(out) :: x_min(:)
   !> Kabsch aligned RMS deviation between the endpoint and the minimum (bohr)
   real(wp), intent(out) :: rmsd
   !> Energy of the relaxed minimum (hartree)
   real(wp), intent(out) :: e_min
   !> Macro steps the relaxation used
   integer, intent(out) :: n_relax
   !> True if the relaxation converged
   logical, intent(out) :: converged
   !> Failure flag
   logical, intent(out) :: fail

   type(landscape_opt_settings) :: opt
   type(landscape_result) :: scratch
   real(wp), allocatable :: xyz_end(:, :), xyz_min(:, :)
   integer :: nat

   fail = .false.
   converged = .false.
   rmsd = 0.0_wp
   e_min = 0.0_wp
   n_relax = 0
   opt%write_files = .false.
   opt%printlevel = 0
   opt%engine = settings%engine
   opt%engine%etemp = settings%etemp
   ! the convergence thresholds of the relaxation are the task thresholds,
   ! which `--lirc-level` (`set_conv_level`) overwrites (the same copy as
   ! `tasks/ts.f90:1069-1075` and `scan/relax_scan.f90:577-585`)
   opt%max_steps = settings%max_steps
   opt%e_tol = settings%e_tol
   opt%g_max = settings%g_max
   opt%g_rms = settings%g_rms
   opt%d_max = settings%d_max
   opt%d_rms = settings%d_rms
   call landscape_optimize_model(env, model, opt, x_end, scratch, x_min, &
      & frozen=landscape_frozen_mask(size(x_end)/3))
   if (.not. allocated(x_min)) then
      call env%error("the IRC endpoint relaxation returned no geometry", source)
      fail = .true.
      return
   end if
   converged = scratch%converged
   e_min = scratch%energy
   n_relax = scratch%n_steps
   nat = size(x_end)/3
   allocate(xyz_end(3, nat), xyz_min(3, nat))
   call geometry_unpack(x_end, xyz_end)
   call geometry_unpack(x_min, xyz_min)
   call align_geometry(env, xyz_end, xyz_min)
   rmsd = vec_rms(geometry_pack(xyz_end) - geometry_pack(xyz_min))

end subroutine relax_endpoint

!> True if an endpoint counts as having reached a minimum.
!>
!> The path must have been stopped by signal (b) (a stationary point) or (c)
!> (the energy rose again after the minimum); a direction cut off by the step
!> budget is not assigned to a minimum.  Signal (d), the limit cycle guard,
!> is the **truncation** of a stalled walk: the endpoint is one phase of the
!> oscillation and is `step` scale away from the minimum, so the guard alone
!> does **not** assign it to a minimum.  It is accepted only when an endpoint
!> relaxation independently confirms the geometry, i.e. when the relaxation
!> was requested **and** converged **and** the Kabsch aligned RMSD to the
!> relaxed minimum is below `rmsd_end` - the same distance check as for (b)
!> and (c).  The RMSD criterion needs a relaxed minimum: with
!> `relax_requested = .false.` there is none, the criterion is therefore not
!> applied (the caller reports `rmsd = 0` and `irc.log` records that it was
!> skipped instead of silently claiming a satisfied distance check).
pure function endpoint_reached(term, relax_requested, relaxed, rmsd, rmsd_end) &
      & result(reached)
   !> Termination signal
   integer, intent(in) :: term
   !> True if the endpoint relaxation was requested
   logical, intent(in) :: relax_requested
   !> True if the endpoint relaxation converged
   logical, intent(in) :: relaxed
   !> Kabsch aligned RMSD between the endpoint and the relaxed minimum (bohr);
   !> only meaningful when `relax_requested` is true
   real(wp), intent(in) :: rmsd
   !> RMSD threshold below which the endpoint counts as reached the minimum
   real(wp), intent(in) :: rmsd_end
   !> True if the endpoint is assigned to a minimum
   logical :: reached

   reached = (term == p_landscape_irc_term_grad .or. &
      & term == p_landscape_irc_term_energy_rise)
   ! the limit cycle guard truncates a stalled walk: without an independent
   ! confirmation of the geometry the endpoint is **not** a minimum, only the
   ! relay above (the relaxation plus its RMSD) may accept it
   if (term == p_landscape_irc_term_limit_cycle) reached = relax_requested
   if (relax_requested) then
      reached = reached .and. relaxed .and. rmsd <= rmsd_end
   end if

end function endpoint_reached

!> Copy one path into the result arrays.
subroutine fill_path(path, nat, xyz, e, s, gnorm, n_eval)
   !> Path
   type(landscape_irc_path), intent(in) :: path
   !> Number of atoms
   integer, intent(in) :: nat
   !> Geometry array `(3, nat, n)`, allocated on return
   real(wp), allocatable, intent(out) :: xyz(:, :, :)
   !> Energies, allocated on return
   real(wp), allocatable, intent(out) :: e(:)
   !> Arc lengths, allocated on return
   real(wp), allocatable, intent(out) :: s(:)
   !> Gradient norms, allocated on return
   real(wp), allocatable, intent(out) :: gnorm(:)
   !> Evaluation counts, allocated on return
   integer, allocatable, intent(out) :: n_eval(:)

   integer :: k

   allocate(xyz(3, nat, path%n), source=0.0_wp)
   do k = 1, path%n
      call geometry_unpack(path%xyz(:, k), xyz(:, :, k))
   end do
   allocate(e(path%n), source=0.0_wp)
   allocate(s(path%n), source=0.0_wp)
   allocate(gnorm(path%n), source=0.0_wp)
   allocate(n_eval(path%n), source=0)
   if (path%n > 0) then
      e = path%e(1:path%n)
      s = path%s(1:path%n)
      gnorm = path%gnorm(1:path%n)
      n_eval = path%n_eval(1:path%n)
   end if

end subroutine fill_path

!> Write the four IRC output files.
subroutine write_irc_files(env, mol, settings, result, x_ts)
   character(len=*), parameter :: source = 'landscape_irc_write_files'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure (transition state)
   type(TMolecule), intent(in) :: mol
   !> IRC settings
   type(landscape_irc_settings), intent(in) :: settings
   !> IRC result
   type(landscape_irc_result), intent(in) :: result
   !> Transition state geometry, flat Cartesian (bohr)
   real(wp), intent(in) :: x_ts(:)

   real(wp), allocatable :: xyz_ts(:, :)

   allocate(xyz_ts(3, mol%n))
   call geometry_unpack(x_ts, xyz_ts)
   call write_direction_xyz(env, mol, trim(settings%output_prefix)//'_forward.xyz', &
      & xyz_ts, result%energy_ts, result%xyz_fwd, result%e_fwd, result%s_fwd)
   if (result%n_bwd > 0) then
      call write_direction_xyz(env, mol, trim(settings%output_prefix)//'_backward.xyz', &
         & xyz_ts, result%energy_ts, result%xyz_bwd, result%e_bwd, result%s_bwd)
   end if
   call write_path_xyz(env, mol, trim(settings%output_prefix)//'_path.xyz', &
      & xyz_ts, result)
   call write_irc_log(env, settings, result)

end subroutine write_irc_files

!> Write one direction to a multi-frame xyz file (transition state first).
subroutine write_direction_xyz(env, mol, name, xyz_ts, energy_ts, xyz, energies, s)
   character(len=*), parameter :: source = 'landscape_irc_write_direction_xyz'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure
   type(TMolecule), intent(in) :: mol
   !> File name relative to the namespace
   character(len=*), intent(in) :: name
   !> Transition state geometry `(3, nat)` (bohr)
   real(wp), intent(in) :: xyz_ts(:, :)
   !> Transition state energy (hartree)
   real(wp), intent(in) :: energy_ts
   !> Path geometries `(3, nat, n)` (bohr)
   real(wp), intent(in) :: xyz(:, :, :)
   !> Path energies (hartree)
   real(wp), intent(in) :: energies(:)
   !> Path arc lengths (bohr*amu^1/2)
   real(wp), intent(in) :: s(:)

   type(TMolecule) :: buffer
   character(len=:), allocatable :: fname
   integer :: unit, ierr, k, nframe

   if (size(s) /= size(xyz, 3)) then
      call env%error("the IRC path has fewer arc lengths than frames", source)
      return
   end if
   fname = get_namespace(name)
   open(newunit=unit, file=fname, action='write', iostat=ierr)
   if (ierr /= 0) then
      call env%error("cannot open "//fname, source)
      return
   end if
   ! `init` instead of a derived type assignment (gfortran 12.2 stride bug)
   call init(buffer, mol%at, mol%xyz, chrg=mol%chrg, uhf=mol%uhf, &
      & lattice=mol%lattice, pbc=mol%pbc)
   nframe = 0
   nframe = nframe + 1
   buffer%xyz = xyz_ts
   call writeMolecule(buffer, unit, format=fileType%xyz, energy=energy_ts, number=nframe)
   do k = 1, size(xyz, 3)
      nframe = nframe + 1
      buffer%xyz = xyz(:, :, k)
      call writeMolecule(buffer, unit, format=fileType%xyz, energy=energies(k), &
         & number=nframe)
   end do
   close(unit)
   write(env%unit, '(1x,a,a,a,i0,a)') "wrote ", fname, " with ", nframe, " frames"

end subroutine write_direction_xyz

!> Write the merged path, ordered by signed arc length and marking the TS.
subroutine write_path_xyz(env, mol, name, xyz_ts, result)
   character(len=*), parameter :: source = 'landscape_irc_write_path_xyz'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure
   type(TMolecule), intent(in) :: mol
   !> File name relative to the namespace
   character(len=*), intent(in) :: name
   !> Transition state geometry `(3, nat)` (bohr)
   real(wp), intent(in) :: xyz_ts(:, :)
   !> IRC result
   type(landscape_irc_result), intent(in) :: result

   character(len=:), allocatable :: fname
   real(wp) :: emax
   integer :: unit, ierr, k, iframe, its

   fname = get_namespace(name)
   open(newunit=unit, file=fname, action='write', iostat=ierr)
   if (ierr /= 0) then
      call env%error("cannot open "//fname, source)
      return
   end if

   ! find the highest energy frame (the transition state)
   emax = result%energy_ts
   its = 0
   do k = 1, result%n_bwd
      if (result%e_bwd(k) > emax) then
         emax = result%e_bwd(k)
         its = k
      end if
   end do
   do k = 1, result%n_fwd
      if (result%e_fwd(k) > emax) then
         emax = result%e_fwd(k)
         its = -k
      end if
   end do

   iframe = 0
   ! backward direction, reversed: signed arc length from the most negative up
   do k = result%n_bwd, 1, -1
      iframe = iframe + 1
      call write_path_frame(unit, mol, result%xyz_bwd(:, :, k), result%e_bwd(k), &
         & -result%s_bwd(k), (its == k), iframe)
   end do
   iframe = iframe + 1
   call write_path_frame(unit, mol, xyz_ts, result%energy_ts, 0.0_wp, its == 0, iframe)
   do k = 1, result%n_fwd
      iframe = iframe + 1
      call write_path_frame(unit, mol, result%xyz_fwd(:, :, k), result%e_fwd(k), &
         & result%s_fwd(k), (its == -k), iframe)
   end do
   close(unit)
   write(env%unit, '(1x,a,a,a,i0,a)') "wrote ", fname, " with ", iframe, " frames"

end subroutine write_path_xyz

!> Write one frame of `irc_path.xyz` with its signed arc length and TS marker.
subroutine write_path_frame(unit, mol, xyz, energy, s, is_ts, iframe)
   !> Output unit
   integer, intent(in) :: unit
   !> Molecular structure (atom identities)
   type(TMolecule), intent(in) :: mol
   !> Geometry `(3, nat)` (bohr)
   real(wp), intent(in) :: xyz(:, :)
   !> Energy (hartree)
   real(wp), intent(in) :: energy
   !> Signed arc length from the transition state (bohr*amu^1/2)
   real(wp), intent(in) :: s
   !> True for the transition state frame
   logical, intent(in) :: is_ts
   !> Frame index
   integer, intent(in) :: iframe

   character(len=8) :: marker
   integer :: i

   marker = ''
   if (is_ts) marker = ' TS'
   write(unit, '(i0)') mol%n
   write(unit, '(a,i0,a,f20.12,a,f16.6,a)') 'frame ', iframe, ' energy: ', &
      & energy, ' s: ', s, trim(marker)
   do i = 1, mol%n
      write(unit, '(1x,a2,1x,3(f18.10,1x))') trim(toSymbol(mol%at(i))), &
         & xyz(1, i), xyz(2, i), xyz(3, i)
   end do

end subroutine write_path_frame

!> Write the human readable IRC log.
subroutine write_irc_log(env, settings, result)
   character(len=*), parameter :: source = 'landscape_irc_write_log'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> IRC settings
   type(landscape_irc_settings), intent(in) :: settings
   !> IRC result
   type(landscape_irc_result), intent(in) :: result

   character(len=:), allocatable :: fname
   integer :: unit, ierr, k, atom, cart

   fname = get_namespace(trim(settings%output_prefix)//'.log')
   open(newunit=unit, file=fname, action='write', iostat=ierr)
   if (ierr /= 0) then
      call env%error("cannot open "//fname, source)
      return
   end if
   write(unit, '(1x,a,a)') "landscape intrinsic reaction coordinate: ", &
      & landscape_irc_method_name(settings%method)
   write(unit, '(1x,a,f12.3,a)') "electronic temperature (frozen): ", &
      & result%etemp, " K"
   write(unit, '(1x,a,es13.6,a)') "step length: ", result%step_length, &
      & " bohr*amu^1/2"
   write(unit, '(1x,a,i0)') "imaginary modes of the transition state: ", result%n_imag
   atom = (result%phase_index - 1)/3 + 1
   cart = mod(result%phase_index - 1, 3) + 1
   write(unit, '(1x,a,f14.8,a,i0,a,i0,a)') "phase component (forward): ", &
      & result%phase_value, " (cartesian ", cart, " of atom ", atom, ")"
   write(unit, '(1x,a,l1)') "directions exchanged: ", result%swapped
   write(unit, '(1x,a)') ""
   write(unit, '(1x,a)') "forward direction"
   write(unit, '(1x,a)') "   point          s/(bohr*amu^1/2)        energy/Eh      "//&
      & "|g|/(Eh/bohr)   E/g"
   do k = 1, result%n_fwd
      write(unit, '(1x,i8,1x,f20.10,1x,f20.12,1x,es15.6,1x,i6)') k, &
         & result%s_fwd(k), result%e_fwd(k), result%gnorm_fwd(k), &
         & result%n_eval_fwd(k)
   end do
   if (result%n_fwd > 0) then
      write(unit, '(1x,a,f10.2)') "mean energy/gradient evaluations per point: ", &
         & real(sum(result%n_eval_fwd(1:result%n_fwd)), wp)/real(result%n_fwd, wp)
   end if
   if (result%n_bwd > 0) then
      write(unit, '(1x,a)') ""
      write(unit, '(1x,a)') "backward direction"
      write(unit, '(1x,a)') "   point          s/(bohr*amu^1/2)        energy/Eh      "//&
         & "|g|/(Eh/bohr)   E/g"
      do k = 1, result%n_bwd
         write(unit, '(1x,i8,1x,f20.10,1x,f20.12,1x,es15.6,1x,i6)') k, &
            & result%s_bwd(k), result%e_bwd(k), result%gnorm_bwd(k), &
            & result%n_eval_bwd(k)
      end do
      write(unit, '(1x,a,f10.2)') "mean energy/gradient evaluations per point: ", &
         & real(sum(result%n_eval_bwd(1:result%n_bwd)), wp)/real(result%n_bwd, wp)
   end if
   write(unit, '(1x,a)') ""
   if (result%n_fwd > 0) then
      write(unit, '(1x,a,f20.12,a)') "forward endpoint energy:  ", &
         & result%e_fwd(result%n_fwd), " Eh"
      write(unit, '(1x,a,f16.6,a)') "forward arc length:       ", &
         & result%s_fwd(result%n_fwd), " bohr*amu^1/2"
   end if
   if (result%n_bwd > 0) then
      write(unit, '(1x,a,f20.12,a)') "backward endpoint energy: ", &
         & result%e_bwd(result%n_bwd), " Eh"
      write(unit, '(1x,a,f16.6,a)') "backward arc length:      ", &
         & result%s_bwd(result%n_bwd), " bohr*amu^1/2"
   end if
   write(unit, '(1x,a,a)') "termination forward:  ", &
      & landscape_irc_term_name(result%term_fwd)
   if (settings%both_directions) then
      write(unit, '(1x,a,a)') "termination backward: ", &
         & landscape_irc_term_name(result%term_bwd)
   else
      write(unit, '(1x,a)') "termination backward: not requested (--lirc-oneway)"
   end if
   write(unit, '(1x,a,a)') "arc length s is ", arc_length_definition(settings%method)
   write(unit, '(1x,a,f10.8)') "endpoint rmsd criterion: rmsd <= rmsd_end = ", &
      & settings%rmsd_end
   write(unit, '(1x,a)') "endpoint rmsd is Kabsch aligned (mass free proper "//&
      & "rotation onto the relaxed minimum)"
   if (settings%endpoint_relax) then
      write(unit, '(1x,a,l1,a,f12.6,a)') "forward reached minimum:  ", &
         & result%reached_min_fwd, " (rmsd ", result%rmsd_to_min_fwd, " bohr)"
      write(unit, '(1x,a,f20.12,a,i0,a)') "forward relaxed minimum:   ", &
         & result%e_min_fwd, " Eh (", result%n_relax_fwd, " relaxation steps)"
      if (settings%both_directions) then
         write(unit, '(1x,a,l1,a,f12.6,a)') "backward reached minimum: ", &
            & result%reached_min_bwd, " (rmsd ", result%rmsd_to_min_bwd, " bohr)"
         write(unit, '(1x,a,f20.12,a,i0,a)') "backward relaxed minimum:  ", &
            & result%e_min_bwd, " Eh (", result%n_relax_bwd, " relaxation steps)"
      end if
   else
      write(unit, '(1x,a,l1)') "forward reached minimum:  ", &
         & result%reached_min_fwd
      write(unit, '(1x,a)') "the rmsd_end criterion was skipped: no endpoint "//&
         & "relaxation was requested (--lirc-no-endpoint-relax)"
   end if
   write(unit, '(1x,a,i0,a,i0,a)') "points whose sphere relaxation was pinned "//&
      & "at the micro cycle budget: forward ", result%n_pinned_fwd, " of ", &
      & result%n_fwd, " points"
   if (settings%both_directions) then
      write(unit, '(1x,a,i0,a,i0,a)') "                            backward ", &
         & result%n_pinned_bwd, " of ", result%n_bwd, " points"
   end if
   close(unit)
   write(env%unit, '(1x,a,a)') "wrote ", fname

end subroutine write_irc_log

!> Human readable definition of the arc length `s` of an integrator.
!>
!> The four integrators measure `s` differently and the values are therefore
!> **not** comparable across methods (LOW-4 of the M6 review): `gs`, `hpc`
!> and `eulerpc` accumulate the **chord length** `|q_{k+1} - q_k|` (the
!> displacement actually taken), while `lqa` integrates the arc length of the
!> quadratic propagation along the path (`s_reached`, the value whose
!> termination the arc length integral defines).  On the same path the two
!> differ by ~5 % (measured: GS 4.0512 versus LQA 4.2500 bohr*amu^1/2 on the
!> GFN2 HCN/HNC transition state, where 0.05 + 42*0.10 = 4.25 is the exact
!> requested value), which is a definition difference, not a longer path.
pure function arc_length_definition(method) result(text)
   !> Integration method, one of `p_landscape_irc_*`
   integer, intent(in) :: method
   !> One sentence for `irc.log`
   character(len=:), allocatable :: text

   select case (method)
   case (p_landscape_irc_lqa)
      text = "the integrated path length of the quadratic propagation"
   case default
      text = "the accumulated chord length of the macro steps (not "//&
         & "comparable with the integrated length of --lirc-method lqa)"
   end select

end function arc_length_definition

!> Why the finite difference parameters of the numerical Hessian are not used.
!>
!> Empty when this run consumes `settings%engine%hess_step` and
!> `settings%engine%hess_acc`, otherwise a one line reason for the banner.
!> Every integrator builds the initial Hessian of the reaction path from the
!> provider, so a numerical source consumes the step; the provider ignores it
!> for every other source.  Measured on the bent HCN transition state: all four
!> integrators change `irc.log` with `--lirc-hess-step 0.2` (1 Hessian build),
!> while `--lirc-hessian model` aborts at the transition state confirmation
!> before the path is integrated and is byte identical with and without the
!> step - see `docs/_simplify-recon/35-stage2b-fix.md` section T1.
function irc_unused_hessian_step(settings) result(reason)
   !> Reaction path settings
   type(landscape_irc_settings), intent(in) :: settings
   !> Empty when the finite difference step reaches a provider, else the reason
   character(len=:), allocatable :: reason

   reason = ""
   if (landscape_lowercase(settings%hessian_source) /= 'numerical') then
      reason = "the Hessian source '"//trim(settings%hessian_source)//&
         & "' has no finite differences"
   end if

end function irc_unused_hessian_step

end module xtb_landscape_tasks_irc
