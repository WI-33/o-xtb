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

!> Frequency task of the landscape suite (PLAN section 9.4).
!>
!> The task driver is the task level entry point of the `vib` run type: it
!> builds the Hessian through `xtb_landscape_vib_numhess` (i.e. through the
!> configured `hessian_source`), analyses it (`xtb_landscape_vib_modes`) and
!> produces the output.
!>
!> Output files, all with the names and formats of the existing `--hess` run
!> so that existing tool chains keep working:
!>
!> | file                | written by |
!> |---------------------|------------|
!> | `hessian`           | this task (TURBOMOLE layout, non mass weighted, projected) |
!> | `vibspectrum`       | `main_freq` (reused) |
!> | `vib_normal_modes`  | this task (`writeNormalModesTurbomole`, the routine `main_freq` uses) |
!> | `g98.out`           | `main_freq` (reused) |
!> | `xtbhess.xyz`       | `main_freq` -> `distort`, only if imaginary modes exist |
!>
!> `write_files` of the settings controls the two files this suite writes
!> itself (`hessian`, `vib_normal_modes`); the reused xtb printout has no
!> output switch and always produces `vibspectrum`, `g98.out` and the
!> thermodynamic table.  All names go through the xtb file namespace, exactly
!> like the files of `--hess`.
!>
!> The result reports the reference energy, so the caller can write it back
!> into `res`: without that write-back the generic total energy line would
!> report the energy of an unrelated state (`--lopt` had exactly that bug).
!>
!> Units: Hessian elements in hartree/bohr^2, frequencies in cm^-1, energies
!> in hartree; the printed kcal/mol conversion is labelled in the output.
module xtb_landscape_tasks_vib
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : autokcal, autoamu
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   use xtb_type_restart, only : TRestart
   use xtb_type_calculator, only : TCalculator
   use xtb_setparam, only : get_namespace
   use xtb_io_writer_turbomole, only : writeNormalModesTurbomole
   use xtb_landscape_settings, only : landscape_vib_settings, landscape_restart_name, &
      & landscape_log_failed, landscape_lowercase, landscape_int_to_string
   use xtb_landscape_results, only : landscape_result, landscape_vib_result
   use xtb_landscape_coords, only : vec_norm, geometry_pack
   use xtb_landscape_hessian, only : landscape_write_hessian_file
   use xtb_landscape_constraints, only : landscape_frozen_mask
   use xtb_compliance, only : compliance_driver
   use xtb_landscape_vib_numhess, only : landscape_vib_hessian
   use xtb_landscape_vib_modes, only : landscape_vib_analysis
   use xtb_landscape_vib_thermo, only : landscape_vib_thermo_output, &
      & landscape_vib_set_thermo_settings, landscape_vib_check_thermo_settings
   use xtb_landscape_json, only : landscape_json_document, &
      & landscape_json_document_init, landscape_json_add_result, &
      & landscape_json_add_core, &
      & landscape_json_add_counter, landscape_json_add_value, &
      & landscape_json_add_flag, landscape_json_add_series, &
      & landscape_json_write, &
      & landscape_json_status_ok, landscape_json_status_incomplete
   implicit none(type, external)
   private

   public :: landscape_vib

contains

!> Task level entry point: harmonic frequencies and thermochemistry.
!>
!> On return `result` holds the reference energy and the summary of the
!> frequency analysis; `mol` is left at its input geometry (the imaginary mode
!> distortion of the reused printout restores it).
subroutine landscape_vib(env, mol, chk, calc, settings, result)
   character(len=*), parameter :: source = 'landscape_vib'
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
   !> Task result
   type(landscape_result), intent(inout) :: result

   type(landscape_vib_result) :: vib
   real(wp), allocatable :: hess(:, :), gradient(:, :)
   real(wp) :: energy, final_err, lowest, highest, etemp = 0.0_wp
   logical :: fail, ok
   integer :: nlog, n_eval

   fail = .false.
   allocate(hess(3*mol%n, 3*mol%n), gradient(3, mol%n))

   ! Merge the thermochemical parameters (priority: explicit command line >
   ! xcontrol `$thermo` > builtin default) and reject the values the reused
   ! RRHO model cannot handle, before any single point is evaluated.  The check
   ! reads the **merged** parameter state, so a bad value of the `$thermo`
   ! block (`scale=-1`, `temp=-5`, `imagthr=20`, ...) is refused exactly like
   ! the corresponding `--lvib-*` flag; validating only the settings object
   ! missed the block, because `set_thermo` checks nothing
   ! (`landscape_vib_check_thermo_settings`).
   call landscape_vib_set_thermo_settings(settings)
   call landscape_vib_check_thermo_settings(env, ok)
   if (.not. ok) then
      result%message = "frequency analysis failed"
      return
   end if

   if (settings%printlevel > 0) then
      write(env%unit, '(1x,a)') "landscape: harmonic frequencies"
      write(env%unit, '(1x,a,a)') "restart policy:           ", &
         & landscape_restart_name(settings%engine%restart_mode)
      write(env%unit, '(1x,a,a)') "Hessian source:           ", &
         & trim(settings%hessian_source)
      if (allocated(settings%hessian_file)) then
         write(env%unit, '(1x,a,a)') "Hessian file:             ", &
            & trim(settings%hessian_file)
      end if
      write(env%unit, '(1x,a,es15.6,a)') "Hessian step:             ", &
         & settings%engine%hess_step, " bohr"
      write(env%unit, '(1x,a,es15.6)') "Hessian SCC accuracy:     ", &
         & settings%engine%hess_acc
      write(env%unit, '(1x,a,l1)') "project translations/rotations: ", &
         & settings%project_tr
      ! R5 of CLI redesign section 3.5: `e_tol` belongs to the level-driven
      ! axes of the shared settings family.  The frequency task has no
      ! convergence loop and `--lvib` has no level of its own, so the value is
      ! echoed for completeness and marked as not consumed, exactly like the
      ! other knobs of this run type that no code path reads.
      write(env%unit, '(1x,a,es15.6,a)') "energy tolerance:         ", &
         & settings%e_tol, " Eh"
      write(env%unit, '(1x,a)') &
         & "(not used: the frequency task has no convergence loop)"
   end if

   ! A file based Hessian without a file name is a configuration error; it is
   ! reported here with the flag that is missing from the command line (the
   ! provider factory refuses it too, with a more generic message).  The source
   ! is compared case insensitively, so any spelling that the validator accepts
   ! reaches this guard (a literal comparison turned `--lvib-hessian FILE` into
   ! the generic factory message).
   if (landscape_lowercase(settings%hessian_source) == 'file') then
      if (.not. allocated(settings%hessian_file)) then
         result%message = "no Hessian file given"
         call env%error("the Hessian source 'file' needs --lvib-hessian-file", source)
         return
      end if
   end if

   call landscape_vib_hessian(env, mol, chk, calc, settings, hess, energy, &
      & gradient, final_err, fail, n_eval, etemp)
   if (fail) then
      result%message = "frequency analysis failed"
      return
   end if

   ! Legacy `compliance.dat` (CLI redesign section 3.9): the native `--hess`
   ! calls `compliance_driver` right after the numerical Hessian
   ! (`src/hessian.F90:331`), and the driver writes the file under its **bare**
   ! name (`src/compliance.f90:287`, no namespace).  The suite calls the same
   ! driver so that a `--lvib` run produces the file the `--hess` tool chain
   ! expects.
   !
   ! The **condition is deliberately not upstream's**: `src/hessian.F90:331`
   ! tests `freezeset%n == 0`, i.e. the `$fix freeze:` list, which is read only
   ! for the numerical Hessian and is a *different* list from the `$fix atoms:`
   ! list (`src/constrain_param.f90:335-357` writes `fixset` for `atoms` and
   ! `freezeset` for `freeze`; `man/xcontrol.7.adoc:212-216`).  Both directions
   ! are reachable: `$fix atoms: ...` leaves `freezeset%n == 0` and upstream
   ! would write the file, `$fix freeze: ...` sets it and upstream would not.
   ! The suite rejects the atoms of `fixset` here, because `fixset` /
   ! `landscape_frozen_mask` is the **single** `$fix` source of truth of the
   ! suite (`types/constraints.f90`).  The frequency analysis never reads that
   ! mask: the fixed atoms get a zero gradient at the calculator layer
   ! (`src/xtb/calculator.f90`), and the finite difference Hessian is built from
   ! gradient differences (`types/hessian.f90`), so the rows and columns of the
   ! frozen coordinates are zero - the compliance analysis of this run must not
   ! cover coordinates the Hessian does not contain.  Registered as a
   ! deviation in `docs/_simplify-recon/11-cli-redesign.md` section 3.9;
   ! `$fix freeze:` is not honoured by the suite at all, so a run with only a
   ! `freeze:` list takes the same branch as one without any `$fix` block.
   ! The suite refuses periodic systems before the task runs
   ! (`landscape_driver.f90`), so `mol%npbc == 0` documents the parity.
   if (mol%n > 1 .and. .not. any(landscape_frozen_mask(mol%n)) .and. &
      & mol%npbc == 0) then
      call compliance_driver(env, mol, hess)
   end if
   ! The electronic energy of the analysed geometry is part of the frequency
   ! result: the `--ohess` chain of `main.F90` fills the native frequency
   ! container from `result%vib` and `write_energy` prints the thermodynamic
   ! totals as `etot + htot`, so the electronic part has to travel with the
   ! thermal one (`landscape_vib_result%energy`).  It is filled together with
   ! the rest of `vib` below and not here: `landscape_vib_analysis` declares
   ! its result `intent(out)`, which resets **every** component to its default
   ! when it is called (measured: an assignment at this place produced an empty
   ! `fres` bridge and a `TOTAL ENERGY 0.000000000000` line).

   call landscape_vib_analysis(env, mol, hess, settings, vib, fail)
   if (fail) then
      result%message = "frequency analysis failed"
      return
   end if

   if (settings%write_files) then
      nlog = env%nLog
      call write_hessian_file(env, settings, hess)
      if (landscape_log_failed(env, nlog)) then
         result%message = "frequency analysis failed"
         return
      end if
   end if

   call landscape_vib_thermo_output(env, mol, chk, settings, vib, energy, ok)
   if (.not. ok) then
      result%message = "frequency analysis failed"
      return
   end if

   if (settings%write_files) then
      nlog = env%nLog
      call write_normal_modes(env, mol, vib)
      if (landscape_log_failed(env, nlog)) then
         result%message = "frequency analysis failed"
         return
      end if
   end if

   result%energy = energy
   result%gnorm = vec_norm(geometry_pack(gradient))
   result%gmax = result%gnorm
   result%n_eval = n_eval
   result%n_cache_miss = n_eval
   ! The imaginary mode count of the analysis belongs to the result as well:
   ! the printed table and `landscape_vib.json` must not report two different
   ! numbers (`-1` would be the "confirmation not performed" default of the
   ! saddle tasks).
   result%n_imag = vib%n_imag
   ! M4 hard constraint 1: the effective electronic temperature is part of the
   ! result object, not only of the printed engine report (the same read-back
   ! the saddle and IRC tasks perform).
   result%etemp = etemp
   result%converged = .true.
   if (settings%printlevel > 0) then
      write(env%unit, '(1x,a,i0,a,i0,a)') "modes: ", vib%n_vib, " internal, ", &
         & vib%n_zero, " projected zero modes"
      write(env%unit, '(1x,a,l1)') "linear molecule:          ", vib%linear
      write(env%unit, '(1x,a,i0)') "imaginary modes:          ", vib%n_imag
      if (vib%n_vib > 0) then
         lowest = vib%freq(vib%n_zero + 1)
         highest = vib%freq(3*mol%n)
         write(env%unit, '(1x,a,f12.2,a,f12.2,a)') "lowest mode: ", lowest, &
            & " cm^-1, highest mode: ", highest, " cm^-1"
      end if
      write(env%unit, '(1x,a,f20.12,a)') "energy:                   ", energy, " Eh"
      write(env%unit, '(1x,a,f20.12,a)') "energy:                   ", &
         & energy*autokcal, " kcal/mol"
      write(env%unit, '(1x,a,es15.6,a)') "gradient norm:            ", &
         & result%gnorm, " Eh/bohr"
      write(env%unit, '(1x,a,es15.6,a)') "zero point energy:        ", &
         & vib%zpe, " Eh"
      ! H(T) and G(T) are part of the summary as well: the reused
      ! thermodynamic table prints them, but the `--hess` property printout
      ! additionally reports `TOTAL ENTHALPY`, and the end-to-end chain
      ! (`test/unit/landscape_e2e_check.py`, mode `chain`) compares the two
      ! routes from stdout.  The values are the ones of the highest requested
      ! temperature (`freq_results%htot/%gtot`).
      write(env%unit, '(1x,a,f20.12,a)') "total enthalpy H(T):      ", &
         & vib%enthalpy, " Eh"
      write(env%unit, '(1x,a,f20.12,a)') "total free energy G(T):   ", &
         & vib%gibbs, " Eh"
      write(env%unit, '(1x,a,es15.6)') "Hessian residual:         ", final_err
   end if

   if (vib%n_imag > 0) then
      result%message = "frequency analysis completed with imaginary modes"
      call env%warning("the Hessian has "//trim(landscape_int_to_string(vib%n_imag))//&
         & " imaginary mode(s): the geometry is not a minimum", source)
   else
      result%message = "frequency analysis completed"
   end if

   ! The existing `--hess` always warns about a Hessian of a geometry that is
   ! not a stationary point (`src/hessian.F90:161-170`); the same threshold and
   ! the same wording keep the two paths consistent, so a user who runs `--lvib`
   ! instead of `--hess` gets the same reminder that the frequencies belong to
   ! a non-minimum structure.
   if (result%gnorm > 2.0e-3_wp) then
      write(env%unit, '(1x,a)') "!! INCOMPLETELY OPTIMIZED GEOMETRY !!"
      call env%warning("Hessian on incompletely optimized geometry!", source)
   end if

   ! The frequency analysis leaves the task through `result%vib` as well as
   ! through the shared result block: the `--ohess` chain of `main.F90` is a
   ! caller that needs the frequencies, the thermochemical totals and the
   ! electronic energy together (the result fields alone would not be enough).
   ! The energy comes from the single evaluation the Hessian provider made;
   ! recomputing it would shift the last digits (reusing the SCF restart
   ! wavefunction perturbs the next evaluation in the ~1e-7 hartree/bohr range).
   ! It is assigned here, after both calls with an `intent(out)` result.
   vib%energy = energy
   result%vib = vib

   if (settings%ljson) call write_json(env, mol, settings, result, vib)

end subroutine landscape_vib

!> Write the structured result `landscape_vib.json`.
!>
!> `zpe`, `enthalpy` and `gibbs` are the thermochemical quantities of the
!> reused model (`freq_results%zp/%htot/%gtot`, hartree); the end-to-end
!> validation compares them against the `--hess` route
!> (`test/unit/landscape_e2e_check.py`, mode `chain`).  `temperature` and
!> `freq_scal` are the **effective** values of the merged parameter state
!> (`landscape_vib_result%temperature`/`%freq_scal`), not the settings fields:
!> with a `$thermo temp=400` block the settings still carry the 298.15 K
!> default, so reporting them would describe the 400 K `gibbs` of the same
!> document as a 298.15 K result.  An imaginary mode does not
!> make the analysis a failure - it is a legitimate result of a geometry
!> that is not a minimum (the `--hess` path reports it the same way) - but it
!> is `incomplete` rather than `ok`, so a regression diff of the status field
!> catches a changed mode count.
subroutine write_json(env, mol, settings, result, vib)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecule of the run (the core block's geometry)
   type(TMolecule), intent(in) :: mol
   !> Frequency settings
   type(landscape_vib_settings), intent(in) :: settings
   !> Task result
   type(landscape_result), intent(in) :: result
   !> Frequency analysis
   type(landscape_vib_result), intent(in) :: vib

   type(landscape_json_document) :: doc = landscape_json_document()
   logical :: ok
   integer :: status

   status = landscape_json_status_ok
   if (vib%n_imag > 0) status = landscape_json_status_incomplete
   call landscape_json_document_init(doc, 'vib')
   call landscape_json_add_result(doc, result, status)
   call landscape_json_add_core(doc, mol, result%energy, settings%input_file)
   call landscape_json_add_counter(doc, 'n_atoms', size(vib%freq)/3)
   call landscape_json_add_counter(doc, 'n_modes', size(vib%freq))
   call landscape_json_add_counter(doc, 'n_zero', vib%n_zero)
   call landscape_json_add_counter(doc, 'n_vib', vib%n_vib)
   ! `n_imag` is already part of the shared result block, so it is not added a
   ! second time (a duplicate key would break a key based comparison).
   call landscape_json_add_value(doc, 'zpe', vib%zpe)
   call landscape_json_add_value(doc, 'enthalpy', vib%enthalpy)
   call landscape_json_add_value(doc, 'gibbs', vib%gibbs)
   ! The values of the **merged** parameter state (see the routine doc above).
   call landscape_json_add_value(doc, 'temperature', vib%temperature)
   call landscape_json_add_value(doc, 'freq_scal', vib%freq_scal)
   call landscape_json_add_flag(doc, 'linear', vib%linear)
   call landscape_json_add_flag(doc, 'project_tr', settings%project_tr)
   if (allocated(vib%freq)) then
      call landscape_json_add_series(doc, 'frequency', vib%freq)
   end if
   call landscape_json_write(env, doc, ok)

end subroutine write_json

!> Write the projected, non mass weighted Hessian as the `hessian` file.
!>
!> The name and the layout are the ones of the `hessian` file of `--hess`
!> (`wrhess`, `src/freq/io.f90:53`), so the file can be read back by
!> `--lvib --lvib-hessian-file hessian`.  The file namespace is applied once,
!> inside `landscape_write_hessian_file` (exactly like `wrhess` through
!> `open_file`, and like the read path of the file provider); applying it here
!> as well would produce `ns.ns.hessian` and break the write-read round trip.
subroutine write_hessian_file(env, settings, hess)
   character(len=*), parameter :: source = 'landscape_vib_write_hessian_file'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Frequency settings
   type(landscape_vib_settings), intent(in) :: settings
   !> Projected Hessian (hartree/bohr^2)
   real(wp), intent(in) :: hess(:, :)
   !> Success flag
   logical :: ok

   call landscape_write_hessian_file(env, 'hessian', hess, ok)
   if (.not. ok) then
      call env%error("the Hessian file could not be written", source)
      return
   end if
   if (settings%printlevel > 0) then
      write(env%unit, '(1x,a)') "wrote hessian, containing the non-mass-"//&
         & "weighted Hessian matrix in atomic units (Eh/bohr^2)"
   end if

end subroutine write_hessian_file

!> Write the mass weighted normal modes in the TURBOMOLE format.
!>
!> The routine is the one `main_freq` uses for `--hess`; it removes the mass
!> weighting from the modes itself, so the masses are passed in amu (`autoamu`
!> converts the atomic units of `TMolecule%atmass`) and the file is exactly the
!> file the existing run writes for `set%pr_nmtm`.
subroutine write_normal_modes(env, mol, vib)
   character(len=*), parameter :: source = 'landscape_vib_write_normal_modes'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure
   type(TMolecule), intent(in) :: mol
   !> Frequency result
   type(landscape_vib_result), intent(in) :: vib

   real(wp), allocatable :: atmass_amu(:)
   character(len=:), allocatable :: fname
   integer :: unit, ierr

   fname = get_namespace('vib_normal_modes')
   open(newunit=unit, file=fname, action='write', iostat=ierr)
   if (ierr /= 0) then
      call env%error("cannot write "//fname, source)
      return
   end if
   allocate(atmass_amu(mol%n), source=0.0_wp)
   atmass_amu = mol%atmass*autoamu
   call writeNormalModesTurbomole(unit, atmass_amu, vib%modes)
   close(unit)

end subroutine write_normal_modes

end module xtb_landscape_tasks_vib
