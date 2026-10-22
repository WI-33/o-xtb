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

!> Thermochemistry of the `vib` task.
!>
!> This is a **thin wrapper**, not an implementation: xtb already owns a
!> complete rigid-rotor-harmonic-oscillator (quasi-RRHO) model with the
!> Grimme free-rotor interpolation of the low frequency modes, the symmetry
!> number, the imaginary mode treatment and the multi-temperature table
!> (`print_thermo`, `src/main/property.F90:1268`, `thermodyn`,
!> `src/thermo.f90:79`).  Reimplementing any of it would be a second,
!> silently diverging spelling of the same physics (PLAN section 6, "explicitly
!> not done").
!>
!> What the wrapper does:
!>
!>  1. translate `landscape_vib_settings` into the parameter state the existing
!>     model reads (`set%thermotemp`, `set%nthermo`, `set%thermo_sthr`,
!>     `set%thermo_ithr`, `set%thermo_fscal`, `src/setparam.f90:284-291`) under
!>     the priority rule "explicit command line > xcontrol `$thermo` > builtin
!>     default" (`landscape_vib_set_thermo_settings`),
!>  2. fill the `freq_results` container of `xtb_type_data` from the
!>     `landscape_vib_result` (frequencies, mass weighted normal modes,
!>     reduced masses) and the reference energy,
!>  3. call `main_freq` (`src/main/property.F90:695`) - the routine behind
!>     `--hess` - which prints the frequency table, writes `vibspectrum` and
!>     `g98.out`, and reaches `print_thermo`/`thermodyn` for the thermodynamic
!>     table, and
!>  4. read the zero point energy, the imaginary mode count and the **effective**
!>     temperature out of the container.
!>
!> Steps 1 and 4 are also why `landscape_vib_result%temperature`/`%freq_scal`
!> exist: the structured result of the task has to carry the values that were
!> used, not the settings defaults the user may have overridden with a
!> `$thermo` block.
!>
!> `landscape_vib_check_thermo_settings` rejects the **merged** values the model
!> cannot handle; the task calls it before the first single point, so
!> `--lvib-temp 0` and `$thermo temp=-5` are refused by the same code.
!>
!> The `set` object is the only pre-existing mutable state the landscape suite
!> writes; it is written here and nowhere else, immediately before the call,
!> because the existing thermochemistry has no parameter argument.  The
!> settings object stays the source of truth - nothing in this suite reads a
!> thermochemical default from `set`.
!>
!> Not available on this path: IR and Raman intensities.  They are dipole
!> (and polarizability) *derivatives* that xtb obtains from its analytic
!> Hessian (`dipd` in `numhess`); the Hessian values of this suite come from
!> the abstract provider, which carries no such derivatives, so the intensity
!> columns of `vibspectrum` and `g98.out` are printed as zero.  The
!> frequencies and all thermochemical quantities are unaffected.
!>
!> Units: frequencies in cm^-1, temperatures in kelvin, `rotor_cutoff` and
!> `imag_cutoff` in cm^-1, energies in hartree.
module xtb_landscape_vib_thermo
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : autoamu
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   use xtb_type_restart, only : TRestart
   use xtb_type_data, only : freq_results
   use xtb_setparam, only : set
   use xtb_splitparam, only : atmass, init_split
   use xtb_propertyoutput, only : main_freq
   use xtb_landscape_settings, only : landscape_vib_settings, landscape_log_failed
   use xtb_landscape_results, only : landscape_vib_result
   implicit none(type, external)
   private

   public :: landscape_vib_thermo_output
   public :: landscape_vib_set_thermo_settings
   public :: landscape_vib_check_thermo_settings

contains

!> Frequency printout and thermochemistry of an analysed Hessian.
!>
!> The routine prints the frequency table, the reduced masses and the
!> thermodynamic functions, and writes `vibspectrum` and `g98.out` through the
!> reused `main_freq`.  It also reports the zero point energy in
!> `vib%zpe` (hartree) and cross-checks the imaginary mode count of the
!> analysis against the one the thermochemistry finds.
subroutine landscape_vib_thermo_output(env, mol, chk, settings, vib, energy, ok)
   character(len=*), parameter :: source = 'landscape_vib_thermo_output'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure
   type(TMolecule), intent(inout) :: mol
   !> Wavefunction data (only needed by the mode following of `main_freq`)
   type(TRestart), intent(inout) :: chk
   !> Frequency settings
   type(landscape_vib_settings), intent(in) :: settings
   !> Frequency result, receives the zero point energy
   type(landscape_vib_result), intent(inout) :: vib
   !> Electronic energy of the analysed geometry (hartree)
   real(wp), intent(in) :: energy
   !> True if the thermochemistry was produced
   logical, intent(out) :: ok

   type(freq_results) :: res
   integer :: nlog

   ok = .false.

   ! `print_thermo` reaches `axis2` (`src/axis_trafo.f90:59`), which reads the
   ! atomic masses from the global storage of `xtb_splitparam` **in amu**.
   ! The main program and the `thermo` submodule fill it the same way
   ! (`src/prog/main.F90:461`, `src/prog/thermo.f90:103`); the landscape suite
   ! runs without the main program, so it has to do it here - otherwise the
   ! reused routine reads unallocated memory.
   if (.not. allocated(atmass)) call init_split(mol%n)
   atmass = mol%atmass*autoamu

   call landscape_vib_set_thermo_settings(settings)

   ! The effective scaling factor is part of the result: `landscape_vib.json`
   ! reports it, and the scaling factor of the merged state is the one
   ! `print_thermo` applies (`set%thermo_fscal`).
   vib%freq_scal = set%thermo_fscal

   call res%allocate(mol%n)
   ! `n3true` is the number of frequencies the existing printout iterates
   ! over; xtb sets it to 3N unless atoms are frozen, which the landscape
   ! suite never does.
   res%n3true = 3*mol%n
   res%linear = vib%linear
   res%freq = vib%freq
   res%hess = vib%modes
   res%rmass = vib%rmass
   ! No dipole or polarizability derivatives are available from a Hessian
   ! provider (see the module documentation), so the intensity columns are
   ! zero rather than wrong.
   res%dipt = 0.0_wp
   res%polt = 0.0_wp
   res%etot = energy
   res%gnorm = 0.0_wp

   nlog = env%nLog
   call main_freq(env%unit, mol, chk%wfn, res)
   if (landscape_log_failed(env, nlog)) then
      call env%error("the frequency printout failed", source)
      return
   end if

   vib%zpe = res%zp
   ! The temperature **read back from the container**: `main_freq` stores the
   ! value it used for the table and for `gibbs` in `freq_results%temp`
   ! (`set%thermotemp(set%nthermo)`, `src/main/property.F90`), so the result
   ! field cannot diverge from the printout it describes.
   vib%temperature = res%temp
   ! H(T) and G(T) are read back so that the end-to-end chain can assert the
   ! two against the `--hess` route.  `print_thermo` stores the **thermal**
   ! contributions in `freq_results%htot/%gtot` (it prints the totals as
   ! `frqres%etot + frqres%htot` and `frqres%etot + frqres%gtot`,
   ! `src/main/property.F90:1315`), so the electronic energy is added here and
   ! the two fields of the result are the totals the `--hess` property
   ! printout reports as `TOTAL ENTHALPY` and `TOTAL FREE ENERGY`.
   vib%enthalpy = energy + res%htot
   vib%gibbs = energy + res%gtot
   if (res%nimag /= vib%n_imag) then
      call env%warning("the thermochemistry counted a different number of "//&
         & "imaginary modes than the mode analysis; a frequency scaling "//&
         & "factor shifts the imaginary cutoff", source)
   end if
   ok = .true.

end subroutine landscape_vib_thermo_output

!> Push the thermochemical settings into the parameter state of xtb.
!>
!> **Priority: explicit command line > xcontrol `$thermo` > builtin default.**
!> `set%thermotemp`/`set%nthermo`, `set%thermo_sthr`, `set%thermo_ithr` and
!> `set%thermo_fscal` already hold the `$thermo` value when the input file had
!> such a block (`src/set_module.f90:1960-2014`, called from
!> `xtb_constrain_param::read_userdata` before any task runs), and the builtin
!> default of `src/setparam.f90:284-294` otherwise.  A landscape CLI flag
!> therefore only overwrites the value when it was actually given; the
!> `*_was_set` flags of `landscape_vib_settings` carry that information,
!> because a settings field alone cannot distinguish "explicitly 298.15 K"
!> from "the default 298.15 K".
!>
!> The four keys are the complete set the `$thermo` block understands
!> (`src/set_module.f90:1981-2012`): `temp` (a list, `set%nthermo` entries),
!> `sthr`, `imagthr` and `scale`.  Nothing else of the block reaches the
!> thermochemistry, so nothing else has to be merged here.
subroutine landscape_vib_set_thermo_settings(settings)
   !> Frequency settings
   type(landscape_vib_settings), intent(in) :: settings

   integer :: ntemp

   if (settings%temperature_was_set) then
      ntemp = max(1, min(settings%n_temperature, size(set%thermotemp)))
      set%nthermo = ntemp
      set%thermotemp(1:ntemp) = settings%temperature(1:ntemp)
   end if
   if (settings%rotor_cutoff_was_set) set%thermo_sthr = settings%rotor_cutoff
   if (settings%imag_cutoff_was_set) set%thermo_ithr = settings%imag_cutoff
   if (settings%freq_scal_was_set) set%thermo_fscal = settings%freq_scal

end subroutine landscape_vib_set_thermo_settings

!> Reject thermochemical parameters the reused RRHO model cannot handle.
!>
!> The check runs on the **merged** parameter state, not on the settings object:
!> after `landscape_vib_set_thermo_settings` the values of `set%` are the
!> effective ones whatever their source (explicit command line, xcontrol
!> `$thermo` block or builtin default), so `--lvib-temp 0` and
!> `$thermo temp=-5` are caught by the same code.  Validating `settings%`
!> alone missed the block: `xtb_setparam::set_thermo` has no positivity check
!> of its own (`src/set_module.f90`, the `case('scale')` and `case('temp')`
!> branches), so `$thermo scale=-1` was multiplied into every frequency by
!> `print_thermo` and `$thermo imagthr=20` silently inverted the meaning of
!> the imaginary cutoff.  The reused model divides by the temperature
!> (`beta = 1/(kB T)`) and by the scaling factor, and it classifies a mode as
!> imaginary by `freq < imag_cutoff`, so those values produce Inf/NaN entries
!> or inverted flags instead of failing.
!>
!> The caller runs this **before** the first single point, so an unusable
!> parameter costs nothing (and is reported the same way from the command line
!> and from the API).
subroutine landscape_vib_check_thermo_settings(env, ok)
   character(len=*), parameter :: source = 'landscape_vib_thermo_settings'
   character(len=160) :: buffer
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> True if every effective parameter is usable
   logical, intent(out) :: ok

   integer :: ntemp, i

   ok = .true.
   ntemp = max(1, min(set%nthermo, size(set%thermotemp)))

   ! A non-finite value passes every comparison below (a comparison against
   ! NaN is false, a comparison against Inf is only caught on one side), and
   ! both routes accept it: `read(sec,*,iostat=)` of the command line parses
   ! `nan`/`inf`, and the `$thermo` block has no value check of its own.  In a
   ! release build that used to reach the RRHO model and print a table of NaN
   ! with exit status 0, so the whole-array check comes first here.
   do i = 1, ntemp
      if (.not. ieee_is_finite(set%thermotemp(i))) then
         write(buffer, '(a,es15.6,a)') "--lvib-temp ", set%thermotemp(i), &
            & " is not a finite temperature"
         call env%error(trim(buffer)//" (checked after the xcontrol `$thermo` "//&
            & "block and --lvib-temp are merged)", source)
         ok = .false.
      end if
   end do
   if (any(set%thermotemp(1:ntemp) <= 0.0_wp)) then
      call env%error("the thermochemical temperature must be positive "//&
         & "(checked after the xcontrol `$thermo` block and --lvib-temp are "//&
         & "merged)", source)
      ok = .false.
   end if
   if (.not. ieee_is_finite(set%thermo_fscal)) then
      write(buffer, '(a,es15.6,a)') "--lvib-freq-scal ", set%thermo_fscal, &
         & " is not a finite scaling factor"
      call env%error(trim(buffer)//" (checked after the xcontrol `$thermo` "//&
         & "block and --lvib-freq-scal are merged)", source)
      ok = .false.
   end if
   if (set%thermo_fscal <= 0.0_wp) then
      call env%error("the frequency scaling factor must be positive "//&
         & "(checked after the xcontrol `$thermo` block and --lvib-freq-scal "//&
         & "are merged)", source)
      ok = .false.
   end if
   if (.not. ieee_is_finite(set%thermo_ithr)) then
      write(buffer, '(a,es15.6,a)') "--lvib-imag-cutoff ", set%thermo_ithr, &
         & " is not a finite frequency cutoff"
      call env%error(trim(buffer)//" (checked after the xcontrol `$thermo` "//&
         & "block and --lvib-imag-cutoff are merged)", source)
      ok = .false.
   end if
   if (set%thermo_ithr > 0.0_wp) then
      call env%error("the imaginary frequency cutoff must not be positive "//&
         & "(checked after the xcontrol `$thermo` block and "//&
         & "--lvib-imag-cutoff are merged)", source)
      ok = .false.
   end if
   if (.not. ieee_is_finite(set%thermo_sthr)) then
      write(buffer, '(a,es15.6,a)') "--lvib-rotor-cutoff ", set%thermo_sthr, &
         & " is not a finite frequency cutoff"
      call env%error(trim(buffer)//" (checked after the xcontrol `$thermo` "//&
         & "block and --lvib-rotor-cutoff are merged)", source)
      ok = .false.
   end if
   if (set%thermo_sthr < 0.0_wp) then
      call env%error("the rotor cutoff of the quasi-RRHO treatment must not "//&
         & "be negative (checked after the xcontrol `$thermo` block and "//&
         & "--lvib-rotor-cutoff are merged)", source)
      ok = .false.
   end if

end subroutine landscape_vib_check_thermo_settings

end module xtb_landscape_vib_thermo
