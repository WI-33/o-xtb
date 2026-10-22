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

!> Scan task of the landscape suite.
!>
!> This is the task level entry point of `--lscan`: it builds the engine,
!> wraps it into the `landscape_engine_model` of `xtb_landscape_tasks_opt`,
!> drives the 1D (`xtb_landscape_scan1d`) or the 2D
!> (`xtb_landscape_scan2d`) scan, writes the scan log and reports the
!> hysteresis and residual diagnostics.
!>
!> The scan log `xtbscan.log` (honouring the `--namespace` prefix) is written
!> frame by frame with `xtb_io_writer::writeMolecule`, exactly like the
!> existing relaxed scan (`src/scan_driver.f90:98,129`).  Every frame carries
!> the **de-biased** electronic energy (`Etot - Efix`), so the file is a drop
!> in replacement for the one produced by `--scan`, and its trajectory
!> semantics are preserved: frames are written in visit order, which for a
!> bidirectional scan includes the reverse branch.
!>
!> Units: coordinates in bohr, energies in hartree, gradients in hartree/bohr.
module xtb_landscape_tasks_scan
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : autokcal
   use xtb_mctc_filetypes, only : fileType
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule, init
   use xtb_type_restart, only : TRestart
   use xtb_type_calculator, only : TCalculator
   use xtb_io_writer, only : writeMolecule
   use xtb_setparam, only : get_namespace
   use xtb_landscape_settings, only : landscape_scan_settings, landscape_restart_name, &
      & landscape_scan_coord_name, landscape_missing_value, landscape_value_present, &
      & landscape_lowercase, landscape_opt_algorithm_from_name, &
      & landscape_opt_algorithm_name, p_landscape_opt_rfo, landscape_int_to_string
   use xtb_landscape_results, only : landscape_result
   use xtb_landscape_coords, only : geometry_pack, geometry_unpack
   use xtb_landscape_engine, only : landscape_engine, new_landscape_engine, &
      & landscape_hessian_provider, new_landscape_hessian_provider
   use xtb_landscape_tasks_opt, only : landscape_engine_model
   use xtb_landscape_scan_relax, only : landscape_scan_check
   use xtb_landscape_scan1d, only : landscape_scan_profile, landscape_scan1d, &
      & landscape_scan_profile_write
   use xtb_landscape_scan2d, only : landscape_scan_grid, landscape_scan2d, &
      & landscape_scan_grid_write
   use xtb_landscape_json, only : landscape_json_document, &
      & landscape_json_document_init, landscape_json_add_result, &
      & landscape_json_add_core, &
      & landscape_json_add_counter, landscape_json_add_value, &
      & landscape_json_add_flag, landscape_json_add_series, &
      & landscape_json_write, &
      & landscape_json_status_ok, landscape_json_status_incomplete, &
      & landscape_json_status_error
   implicit none(type, external)
   private

   public :: landscape_scan

contains

!> Run the scan described by `settings`.
!>
!> On return `mol%xyz` holds the geometry of the lowest energy scan point,
!> `result%energy` that de-biased energy and `result%converged` whether every
!> scan point converged.  Both kinds of incomplete point make it false: a
!> point whose single point failed carries no energy and is reported as an
!> error (exactly like an unconverged minimisation in
!> `xtb_landscape_tasks_opt`), and a point whose inner relaxation stopped
!> without meeting its thresholds still carries its energy but is reported as
!> a warning together with its per point residual.  The flag therefore tells
!> a caller that only inspects it the same thing as the `conv` column of the
!> scan table (`docs/_landscape-recon/17-m2-independent-review.md`
!> section 3.1).
!>
!> A point whose single point was discarded has **no** energy.  It is reported
!> as NaN in the scan table and left out of `xtbscan.log` (a missing frame is
!> honest, a frame carrying `energy: 0.0` or the all-zero geometry of a point
!> that never produced one is not, see
!> `docs/_landscape-recon/14-m2-validation-experiments.md` section 5.3).  When
!> no point at all carries an energy, `result%energy` is NaN as well and the
!> input geometry is left in place instead of being overwritten by a zero
!> frame.
subroutine landscape_scan(env, mol, chk, calc, settings, result)
   character(len=*), parameter :: source = 'landscape_scan'
   !> Computational environment
   type(TEnvironment), intent(inout), target :: env
   !> Molecular structure
   type(TMolecule), intent(inout), target :: mol
   !> Wavefunction data
   type(TRestart), intent(inout), target :: chk
   !> Single point calculator
   class(TCalculator), intent(in) :: calc
   !> Scan settings
   type(landscape_scan_settings), intent(in) :: settings
   !> Scan result
   type(landscape_result), intent(inout) :: result

   type(landscape_engine), target :: engine
   type(landscape_engine_model) :: model
   class(landscape_hessian_provider), allocatable :: hp
   type(landscape_scan_profile) :: profile
   type(landscape_scan_grid) :: grid
   real(wp), allocatable :: x0(:), xbest(:)
   real(wp) :: best_energy, best_residual
   integer :: n_fail, n_unconverged, n_valid, i
   logical :: fail
   integer :: inner_algorithm
   character(len=:), allocatable :: unused_hess

   if (settings%printlevel > 0) then
      write(env%unit, '(1x,a,a,a)') "landscape: ", trim(settings%mode), &
         & " scan"
      ! The inner minimiser of a relaxed scan is part of the settings and is
      ! handed to the constraint relaxation (`relax_scan.f90`), but no output
      ! file carries it (`xtbscan.log` keeps the `--scan` record layout) and
      ! the milestone above is mode only, so an ignored method was invisible.
      ! It is echoed here for the relaxed mode; a rigid scan builds the
      ! geometry and runs a single point, so it has no inner minimiser and
      ! stays quiet (stage 2 review L2).
      inner_algorithm = landscape_opt_algorithm_from_name(settings%method)
      if (landscape_lowercase(settings%mode) == 'relaxed') then
         write(env%unit, '(1x,a,a,a)') "  inner minimiser:         ", &
            & trim(landscape_opt_algorithm_name(inner_algorithm))
      end if
      do i = 1, settings%ndim
         write(env%unit, '(1x,a,i0,a,a,4(1x,i0))') "  coordinate ", i, ": ", &
            & trim(landscape_scan_coord_name(settings%coord(i)%ctype)), &
            & settings%coord(i)%atoms(1:4)
      end do
      write(env%unit, '(1x,a,i0)') "number of coordinates:    ", settings%ndim
      ! R5 of CLI redesign section 3.5: `e_tol` follows the convergence level
      ! of the scan (`set_conv_level`) and is the energy guard of the inner
      ! relaxation of every scan point (`relax_scan.f90`).  A rigid scan runs
      ! single points and never relaxes, so it says so instead of echoing a
      ! value no code path of that mode reads.
      write(env%unit, '(1x,a,es15.6,a)') "energy tolerance:         ", &
         & settings%e_tol, " Eh"
      if (landscape_lowercase(settings%mode) /= 'relaxed') then
         write(env%unit, '(1x,a)') &
            & "(not used: a rigid scan runs no inner relaxation)"
      end if
      write(env%unit, '(1x,a,a)') "restart policy:           ", &
         & landscape_restart_name(settings%engine%restart_mode)
      write(env%unit, '(1x,a,a)') "Hessian source:           ", &
         & trim(settings%hessian_source)
      ! CLI redesign sections 7.3 and 9.2: every flag that carries a unit must
      ! echo its effective value and unit, and a knob that no code path of the
      ! selected inner minimiser reads must say so.  Both values come from the
      ! settings object the inner minimisation is built from, not from the
      ! command line text.  The predicate is the measured one of
      ! `35-stage2b-fix.md` section T1: only an RFO inner minimiser builds a
      ! numerical Hessian (`optimize/driver.f90:183`), a rigid scan has no
      ! inner minimisation at all, and the provider ignores the finite
      ! difference step for every other Hessian source.
      write(env%unit, '(1x,a,es15.6,a)') "Hessian step:             ", &
         & settings%engine%hess_step, " bohr"
      write(env%unit, '(1x,a,es15.6)') "Hessian SCC accuracy:     ", &
         & settings%engine%hess_acc
      unused_hess = scan_unused_hessian_step(settings, inner_algorithm)
      if (len(unused_hess) > 0) then
         write(env%unit, '(1x,a,a,a)') "(not used: ", unused_hess, ")"
      end if
      ! D1 of the CLI redesign section 8.3: the suite force constant 0.5
      ! hartree/bohr^2 is ten times the xtb scan value 0.05 (`fcconstr`).  The
      ! notice prints the value **read from the settings object** and is
      ! suppressed when `--lscan-fc` named it; a rigid scan constructs the
      ! geometry and never touches the bias potential, so it stays quiet too.
      ! The mode is compared through `landscape_lowercase` because the
      ! predicate that accepts it (`landscape_scan_mode_known`) and
      ! `landscape_scan_check` both match case insensitively: a literal
      ! comparison would let `--lscan-mode RELAXED` run a relaxed scan that
      ! never echoed its force constant.
      if (landscape_lowercase(settings%mode) == 'relaxed' .and. &
            & .not. settings%constraint_fc_was_set) then
         write(env%unit, '(a,f4.2,a)') "[NOTICE] --lscan-fc default ", &
            & settings%constraint_fc, " hartree/bohr^2 differs from the "//&
            & "xtb scan force constant 0.05"
         write(env%unit, '(a)') &
            & "         pass --lscan-fc 0.05 for the legacy xtbscan.log "//&
            & "behaviour"
      end if
   end if

   ! All preconditions are checked before a single SCF is run, so an invalid
   ! request is reported instead of producing a partial scan.
   if (.not. landscape_scan_check(env, settings, mol%n)) then
      result%message = "invalid scan settings"
      return
   end if

   call new_landscape_engine(engine, settings%engine)
   call new_landscape_hessian_provider(env, settings%hessian_source, hp, &
      & sccacc=settings%engine%hess_acc, filename=settings%hessian_file, &
      & printlevel=settings%printlevel)
   if (.not. allocated(hp)) return
   call engine%set_hessian_provider(hp)
   call engine%set_calculator(calc)
   model%engine => engine
   model%env => env
   model%mol => mol
   model%chk => chk

   x0 = geometry_pack(mol%xyz)

   fail = .false.
   allocate(xbest(3*mol%n))
   ! The starting geometry is the fallback: a scan in which no point produced
   ! an energy must not report the all-zero coordinates of an unevaluated
   ! point as its result.
   xbest = x0
   if (settings%ndim == 2) then
      call landscape_scan2d(env, model, mol%at, settings, x0, grid, fail)
      n_fail = grid%n_fail
      n_unconverged = count(.not. grid%converged)
      n_valid = count(grid%energy_valid)
      best_energy = grid%min_energy
      best_residual = 0.0_wp
      if (n_valid > 0) then
         best_residual = grid%residual(grid%min_i, grid%min_j)
         xbest = grid%xyz(:, grid%min_visit)
      end if
      result%n_steps = grid%n_points
   else
      call landscape_scan1d(env, model, mol%at, settings, x0, profile, fail)
      n_fail = profile%n_fail
      n_unconverged = count(.not. profile%converged)
      n_valid = 0
      best_energy = landscape_missing_value()
      best_residual = 0.0_wp
      do i = 1, profile%n_points
         if (.not. profile%converged(i)) cycle
         if (.not. profile%energy_valid(i)) cycle
         n_valid = n_valid + 1
         ! Block if, not `n_valid == 1 .or. ...`: `best_energy` starts from the
         ! NaN sentinel, and Fortran does not require the right operand to be
         ! short circuited, so the comparison with that NaN would be evaluated
         ! on the first valid point and trap under -ffpe-trap=invalid.
         if (n_valid == 1) then
            best_energy = profile%energy(i)
            best_residual = profile%residual(i)
            xbest = profile%xyz(:, i)
         else if (profile%energy(i) < best_energy) then
            best_energy = profile%energy(i)
            best_residual = profile%residual(i)
            xbest = profile%xyz(:, i)
         end if
      end do
      result%n_steps = profile%n_points
   end if

   result%n_eval = engine%n_eval
   result%n_cache_miss = engine%n_eval
   ! M4 hard constraint 1: the effective electronic temperature is part of the
   ! result object (read back from the calculator), not only of the printed
   ! engine report.
   result%etemp = engine%effective_etemp()
   result%energy = best_energy
   result%gnorm = best_residual
   result%gmax = best_residual
   ! Every scan point, not merely every evaluated one: a point whose inner
   ! relaxation stopped short keeps its energy and only reaches the warning
   ! below, so it has to make the aggregate flag false as well.
   result%converged = (n_fail == 0) .and. (n_unconverged == 0)

   call geometry_unpack(xbest, mol%xyz)

   if (settings%printlevel > 0) then
      if (settings%ndim == 2) then
         call landscape_scan_grid_write(grid, env%unit, settings%coord(1)%ctype, &
            & settings%coord(2)%ctype)
      else
         call landscape_scan_profile_write(profile, env%unit, &
            & settings%coord(1)%ctype)
         write(env%unit, '(1x,a,es15.6,a)') "lowest energy: ", best_energy, " Eh"
      end if
      call engine%report(env%unit)
   end if

   if (settings%write_files) then
      if (settings%ndim == 2) then
         call write_scan_grid(env, mol, settings, grid)
      else
         call write_scan_profile(env, mol, settings, profile)
      end if
   end if

   call engine%destroy()

   if (n_fail > 0) then
      result%message = "some scan points failed"
      call env%error("landscape scan aborted: "//landscape_int_to_string(n_fail)//&
         & " of "//landscape_int_to_string(result%n_steps)//" scan points failed (no "//&
         & "energy; the failed points are the NaN rows of the scan table)", &
         & source)
   else if (n_unconverged == result%n_steps .and. result%n_steps > 0) then
      ! Not a single point relaxed its constraints.  The scan itself is not a
      ! failure of the machinery - a point may legitimately stop short of its
      ! convergence thresholds - but there is no minimum to report: both the
      ! one- and the two-dimensional minimum search count only converged points
      ! (`landscape_scan1d`, `landscape_scan2d`), so the energy stays at the
      ! missing value and the reference geometry is the input geometry (the
      ! `0.0` of an unevaluated point must never be reported as a result).
      result%message = "scan completed without a converged point"
      ! M7 close-out of PLAN line 1805: a result that is not available is a
      ! non-zero exit.  Before this round the branch only warned and the
      ! process exited 0 next to a NaN `TOTAL ENERGY` (the semantics recorded
      ! as INFO-3 in PLAN section 12.11); a caller of a scan chain cannot
      ! distinguish that from a valid result.  A scan with **some**
      ! unconverged points keeps exit 0 (see the next branch): it still has a
      ! converged minimum and the residuals make the caveat visible.
      call env%error("no scan point converged: the reported energy is "//&
         & "undefined (it stays NaN) and the reference geometry is the input "//&
         & "geometry; raise the step budget of the inner relaxation "//&
         & "(max_steps)", source)
   else if (n_unconverged > 0) then
      ! An inner relaxation that hit its step budget is not a hard failure - a
      ! scan point may legitimately converge loosely - but it has to be
      ! visible, because the constraint residual of such a point is non-zero
      ! (PLAN R-04).
      result%message = "scan completed with unconverged points"
      call env%warning(landscape_int_to_string(n_unconverged)//" scan point(s) did "//&
         & "not converge; check the reported residual", source)
   else
      result%message = "scan completed"
   end if

   ! The document is written independently of `write_files`: it is the
   ! regression artefact of the run, not a geometry file (see
   ! `landscape_settings%ljson`).
   if (settings%ljson) then
      if (settings%ndim == 2) then
         call write_json_grid(env, mol, settings, result, grid, n_fail, &
            & n_unconverged, n_valid)
      else
         call write_json_profile(env, mol, settings, result, profile, n_fail, &
            & n_unconverged, n_valid)
      end if
   end if

end subroutine landscape_scan

!> Write the structured result `landscape_scan.json` of a 1D scan.
!>
!> The `energy` series is the de-biased energy of every point in visit
!> order, i.e. the quantity of `xtbscan.log`; the counters are the point
!> counts the task prints (`n_points`, `n_valid`, `n_fail`,
!> `n_unconverged`) next to the aggregate `n_eval` of the result.
subroutine write_json_profile(env, mol, settings, result, profile, n_fail, &
      & n_unconverged, n_valid)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecule of the run (the core block's geometry)
   type(TMolecule), intent(in) :: mol
   !> Scan settings
   type(landscape_scan_settings), intent(in) :: settings
   !> Aggregate scan result
   type(landscape_result), intent(in) :: result
   !> 1D scan profile
   type(landscape_scan_profile), intent(in) :: profile
   !> Number of points whose single point failed
   integer, intent(in) :: n_fail
   !> Number of points whose inner relaxation stopped short
   integer, intent(in) :: n_unconverged
   !> Number of points that carry a usable energy
   integer, intent(in) :: n_valid

   type(landscape_json_document) :: doc = landscape_json_document()
   logical :: ok
   integer :: status

   call scan_json_status(result, n_fail, status)
   call landscape_json_document_init(doc, 'scan')
   call landscape_json_add_result(doc, result, status)
   call landscape_json_add_core(doc, mol, result%energy, settings%input_file)
   call landscape_json_add_counter(doc, 'n_points', profile%n_points)
   call landscape_json_add_counter(doc, 'n_forward', profile%n_forward)
   call landscape_json_add_counter(doc, 'n_valid', n_valid)
   call landscape_json_add_counter(doc, 'n_fail', n_fail)
   call landscape_json_add_counter(doc, 'n_unconverged', n_unconverged)
   call landscape_json_add_value(doc, 'kcal_per_mol', result%energy*autokcal)
   call landscape_json_add_flag(doc, 'bidirectional', profile%bidirectional)
   if (allocated(profile%energy)) then
      call landscape_json_add_series(doc, 'energy', profile%energy)
   end if
   if (allocated(profile%reached)) then
      call landscape_json_add_series(doc, 'reached', profile%reached)
   end if
   if (allocated(profile%residual)) then
      call landscape_json_add_series(doc, 'residual', profile%residual)
   end if
   call landscape_json_write(env, doc, ok)

end subroutine write_json_profile

!> Write the structured result `landscape_scan.json` of a 2D grid scan.
!>
!> The `energy` series holds the grid in **visit order** (the order the
!> serpentine sweep evaluated the points), so it lines up with the frames of
!> `xtbscan.log`; `n_points` is `n1*n2`.
subroutine write_json_grid(env, mol, settings, result, grid, n_fail, &
      & n_unconverged, n_valid)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecule of the run (the core block's geometry)
   type(TMolecule), intent(in) :: mol
   !> Scan settings
   type(landscape_scan_settings), intent(in) :: settings
   !> Aggregate scan result
   type(landscape_result), intent(in) :: result
   !> 2D scan grid
   type(landscape_scan_grid), intent(in) :: grid
   !> Number of points whose single point failed
   integer, intent(in) :: n_fail
   !> Number of points whose inner relaxation stopped short
   integer, intent(in) :: n_unconverged
   !> Number of points that carry a usable energy
   integer, intent(in) :: n_valid

   type(landscape_json_document) :: doc = landscape_json_document()
   real(wp), allocatable :: energy(:)
   logical :: ok
   integer :: status, i

   call scan_json_status(result, n_fail, status)
   call landscape_json_document_init(doc, 'scan')
   call landscape_json_add_result(doc, result, status)
   call landscape_json_add_core(doc, mol, result%energy, settings%input_file)
   call landscape_json_add_counter(doc, 'n_points', grid%n_points)
   call landscape_json_add_counter(doc, 'n1', grid%n1)
   call landscape_json_add_counter(doc, 'n2', grid%n2)
   call landscape_json_add_counter(doc, 'n_valid', n_valid)
   call landscape_json_add_counter(doc, 'n_fail', n_fail)
   call landscape_json_add_counter(doc, 'n_unconverged', n_unconverged)
   call landscape_json_add_value(doc, 'kcal_per_mol', result%energy*autokcal)
   if (allocated(grid%energy) .and. allocated(grid%visit_i)) then
      allocate(energy(grid%n_points))
      do i = 1, grid%n_points
         energy(i) = grid%energy(grid%visit_i(i), grid%visit_j(i))
      end do
      call landscape_json_add_series(doc, 'energy', energy)
   end if
   call landscape_json_write(env, doc, ok)

end subroutine write_json_grid

!> Verdict of a scan for the JSON document.
!>
!> A failed point makes the result unusable (status `error`; the task reports
!> it through `env` and exits non-zero), an all-unconverged scan has no
!> minimum either, and a partially unconverged scan is `incomplete` with a
!> usable minimum.
pure subroutine scan_json_status(result, n_fail, status)
   !> Aggregate scan result
   type(landscape_result), intent(in) :: result
   !> Number of points whose single point failed
   integer, intent(in) :: n_fail
   !> Task verdict
   integer, intent(out) :: status

   status = landscape_json_status_ok
   if (.not. result%converged) status = landscape_json_status_incomplete
   if (n_fail > 0) status = landscape_json_status_error
   if (result%n_steps > 0 .and. .not. result%converged) then
      if (.not. landscape_value_present(result%energy)) then
         status = landscape_json_status_error
      end if
   end if

end subroutine scan_json_status

!> Write the 1D scan log, one xyz frame per point with a de-biased energy.
!>
!> Only points that carry an energy are written.  A point whose single point
!> was discarded has neither an energy nor, necessarily, a meaningful
!> geometry (a point that failed before its first step keeps a zero frame), so
!> it is left out instead of producing a frame with a fabricated `0.0` energy.
subroutine write_scan_profile(env, mol, settings, profile)
   character(len=*), parameter :: source = 'landscape_write_scan_profile'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure (atom identities are copied from it)
   type(TMolecule), intent(in) :: mol
   !> Scan settings
   type(landscape_scan_settings), intent(in) :: settings
   !> Scan profile
   type(landscape_scan_profile), intent(in) :: profile

   type(TMolecule) :: buffer
   character(len=:), allocatable :: fname
   integer :: unit, ierr, iframe, nframe, nomitted

   fname = get_namespace('xtbscan.log')
   ! `status='replace'`: a scan that omits every frame (no point produced an
   ! energy) has to leave an empty log, not the frames of an earlier run - an
   ! open without `replace` keeps the old records until a new one is written.
   open(newunit=unit, file=fname, action='write', status='replace', iostat=ierr)
   if (ierr /= 0) then
      call env%error("cannot open "//fname, source)
      return
   end if
   ! The buffer is rebuilt with `init` instead of a derived type assignment,
   ! see the comment in `xtb_landscape_log::landscape_write_files`.
   call init(buffer, mol%at, mol%xyz, chrg=mol%chrg, uhf=mol%uhf, &
      & lattice=mol%lattice, pbc=mol%pbc)
   nframe = 0
   nomitted = 0
   do iframe = 1, profile%n_points
      if (.not. profile%energy_valid(iframe)) then
         nomitted = nomitted + 1
         cycle
      end if
      nframe = nframe + 1
      call geometry_unpack(profile%xyz(:, iframe), buffer%xyz)
      call writeMolecule(buffer, unit, format=fileType%xyz, &
         & energy=profile%energy(iframe), number=nframe)
   end do
   close(unit)
   write(env%unit, '(1x,a,a,a,i0,a,i0,a)') "wrote ", fname, " with ", nframe, &
      & " frame(s), ", nomitted, " failed point(s) omitted"

end subroutine write_scan_profile

!> Write the 2D scan log, one xyz frame per visited grid point with a
!> de-biased energy; failed grid points are omitted (see
!> `write_scan_profile`).
subroutine write_scan_grid(env, mol, settings, grid)
   character(len=*), parameter :: source = 'landscape_write_scan_grid'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure (atom identities are copied from it)
   type(TMolecule), intent(in) :: mol
   !> Scan settings
   type(landscape_scan_settings), intent(in) :: settings
   !> Scan grid
   type(landscape_scan_grid), intent(in) :: grid

   type(TMolecule) :: buffer
   character(len=:), allocatable :: fname
   integer :: unit, ierr, iframe, nframe, nomitted

   fname = get_namespace('xtbscan.log')
   ! See `write_scan_profile` for the `replace` of the previous log.
   open(newunit=unit, file=fname, action='write', status='replace', iostat=ierr)
   if (ierr /= 0) then
      call env%error("cannot open "//fname, source)
      return
   end if
   call init(buffer, mol%at, mol%xyz, chrg=mol%chrg, uhf=mol%uhf, &
      & lattice=mol%lattice, pbc=mol%pbc)
   ! The frames carry the energy of the point whose grid indices are read from
   ! the visit order, so the file is consistent with the serpentine traversal.
   nframe = 0
   nomitted = 0
   do iframe = 1, grid%n_points
      if (.not. grid%energy_valid(grid%visit_i(iframe), grid%visit_j(iframe))) then
         nomitted = nomitted + 1
         cycle
      end if
      nframe = nframe + 1
      call geometry_unpack(grid%xyz(:, iframe), buffer%xyz)
      call writeMolecule(buffer, unit, format=fileType%xyz, &
         & energy=grid%energy(grid%visit_i(iframe), grid%visit_j(iframe)), &
         & number=nframe)
   end do
   close(unit)
   write(env%unit, '(1x,a,a,a,i0,a,i0,a)') "wrote ", fname, " with ", nframe, &
      & " frame(s), ", nomitted, " failed point(s) omitted"

end subroutine write_scan_grid

!> Why the finite difference parameters of the numerical Hessian are not used.
!>
!> Empty when the scan consumes `settings%engine%hess_step` and
!> `settings%engine%hess_acc`, otherwise a one line reason for the banner.  A
!> rigid scan constructs the geometry and runs a single point, so it has no
!> inner minimisation and no Hessian; a relaxed scan builds a numerical
!> Hessian only when the inner minimiser is RFO, and the provider ignores the
!> finite difference step for every other Hessian source.  Measured on GFN2
!> water: `--lscan-method rfo --lscan-hessian numerical` changes
!> `xtbscan.log` with `--lscan-hess-step 0.2` (3 Hessian builds), while the
!> other five inner minimisers, `--lscan-mode rigid` and `--lscan-hessian
!> model` leave it byte identical (0 Hessian builds) - see
!> `docs/_simplify-recon/35-stage2b-fix.md` section T1.
function scan_unused_hessian_step(settings, algorithm) result(reason)
   !> Scan settings
   type(landscape_scan_settings), intent(in) :: settings
   !> Inner minimisation algorithm resolved from `settings%method`
   integer, intent(in) :: algorithm
   !> Empty when the finite difference step reaches a provider, else the reason
   character(len=:), allocatable :: reason

   reason = ""
   if (landscape_lowercase(settings%mode) /= 'relaxed') then
      reason = "a rigid scan builds no Hessian"
   else if (landscape_lowercase(settings%hessian_source) /= 'numerical') then
      reason = "the Hessian source '"//trim(settings%hessian_source)//&
         & "' has no finite differences"
   else if (algorithm /= p_landscape_opt_rfo) then
      reason = trim(landscape_opt_algorithm_name(algorithm))//&
         & " builds no numerical Hessian"
   end if

end function scan_unused_hessian_step

end module xtb_landscape_tasks_scan
