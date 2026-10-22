! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
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

! Modified by the o-xtb landscape task suite, 2026.

module xtb_prog_main
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use xtb_mctc_accuracy, only: wp
   use xtb_mctc_constants, only: pi
   use xtb_mctc_io, only: stderr
   use xtb_mctc_timings
   use xtb_mctc_systools
   use xtb_mctc_convert
   use xtb_mctc_param
   use xtb_type_molecule
   use xtb_type_calculator
   use xtb_type_restart
   use xtb_tblite_restart, only: loadRestart, dumpRestart
   use xtb_type_param
   use xtb_type_data
   use xtb_type_environment, only: TEnvironment, init
   use xtb_prog_argparser
   use xtb_solv_state
   use xtb_setparam
   use xtb_sphereparam
   use xtb_scanparam
   use xtb_splitparam
   use xtb_fixparam
   use xtb_features, only: get_xtb_feature
   use xtb_constrain_param, only: read_userdata
   use xtb_shake, only: init_shake
   use xtb_gfnff_shake, only: gff_init_shake => init_shake
   use xtb_embedding, only: init_pcem
   use xtb_io_reader, only: readMolecule
   use xtb_io_writer, only: writeMolecule
   use xtb_mctc_filetypes, only: fileType, getFileType, generateFileMetaInfo, &
      & generateFileName
   use xtb_readin
   use xtb_printout
   use xtb_setmod
   use xtb_propertyoutput
   use xtb_io_writer_turbomole, only: writeResultsTurbomole
   use xtb_io_writer_orca, only: writeResultsOrca
   use xtb_io_writer_gaussian, only: writeResultsGaussianExternal
   use xtb_restart
   use xtb_readparam
   use xtb_scc_core, only: iniqshell
   use xtb_aespot, only: get_radcn
   use xtb_iniq, only: iniqcn
   use xtb_eeq
   use xtb_disp_ncoord, only: ncoord_gfn, dncoord_erf, dncoord_d3, ncoord_erf, &
      & ncoord_d3
   use xtb_basis
   use xtb_axis, only: axis3
   use xtb_hessian, only: numhess
   use xtb_dynamic, only: md
   use xtb_modef, only: modefollow
   use xtb_mdoptim, only: mdopt
   use xtb_screening, only: screen
   use xtb_xtb_calculator
   use xtb_gfnff_calculator
   use xtb_iff_calculator, only: TIFFCalculator
   use xtb_paramset
   use xtb_xtb_gfn0
   use xtb_xtb_gfn1
   use xtb_xtb_gfn2
   use xtb_main_setup
   use xtb_main_defaults, only: initDefaults
   use xtb_main_json, only: main_xtb_json, write_json_gfnff_lists
   use xtb_geoopt
   use xtb_metadynamic
   use xtb_biaspath
   use xtb_coffee
   use xtb_disp_dftd3param
   use xtb_disp_dftd4
   use xtb_gfnff_param, only: gff_print
   use xtb_gfnff_topology, only: TPrintTopo
   use xtb_gfnff_convert, only: struc_convert
   use xtb_scan
   use xtb_kopt
   use xtb_iff_iffprepare, only: prepare_IFF
   use xtb_iff_data, only: TIFFData
   use xtb_oniom, only: oniom_input, TOniomCalculator, calculateCharge
   use xtb_vertical, only: vfukui
   use xtb_tblite_calculator, only: TTBLiteCalculator, TTBLiteInput, &
         & newTBLiteWavefunction, get_ceh, num_grad_chrg
   use xtb_ptb_calculator, only: TPTBCalculator
   use xtb_solv_cpx, only: TCpcmx
   use xtb_dipro, only: get_jab, jab_input
   !> PTB related modules
   use xtb_main_json, only: main_ptb_json
   !> landscape task suite
   use xtb_landscape, only: run_landscape
   use xtb_landscape_settings, only: landscape_opt_settings, &
      & landscape_scan_settings, landscape_vib_settings, landscape_ts_settings, &
      & landscape_irc_settings, landscape_engine_control, &
      & landscape_opt_algorithm_from_name, landscape_conv_level_known, &
      & landscape_opt_algorithm_name, landscape_ts_method_name, &
      & landscape_irc_method_name, &
      & landscape_ts_method_from_name, &
      & landscape_ts_init_mode_known, landscape_ts_ninit_known, &
      & landscape_irc_method_from_name, &
      & landscape_scan_coord_from_name, landscape_scan_coord_known, &
      & landscape_scan_coord_natoms, landscape_scan_mode_known, &
      & p_landscape_scan_coord_distance, p_landscape_restart_off, &
      & p_landscape_restart_pin, p_landscape_restart_warm, p_landscape_ts_cineb, &
      & p_landscape_ts_autoneb, p_landscape_ts_string, p_landscape_ts_gsm, &
      & landscape_lowercase
   use xtb_landscape_hessian, only: landscape_hessian_source_known
   use xtb_landscape_results, only: landscape_result
   use xtb_landscape_json, only: landscape_json_write_run

   implicit none
   private

   public :: xtbMain

#ifdef _WIN32
   !> Switches the attached console to UTF-8 output, see
   !> `src/type/win32_console.c`: the labels below are UTF-8 (`·`, `²`, `α`,
   !> `⁻¹`), but a Windows console decodes them with its own code page, which
   !> misreads them on a GBK/CP936 system.  `_WIN32` is supplied by CMake's
   !> Windows block (gfortran's Fortran preprocessor defines no target macro),
   !> so on every other target no `xtb_console_utf8` symbol is referenced, and
   !> `win32_console.c` is not compiled either.
   interface
      subroutine c_console_utf8() bind(c, name='xtb_console_utf8')
      end subroutine c_console_utf8
   end interface
#endif

contains

   subroutine xtbMain(env, argParser)

      !> Source of errors in the main program unit
      character(len=*), parameter :: source = "prog_main"

      type(TEnvironment), intent(inout) :: env

      type(TArgParser), intent(inout) :: argParser

!! ========================================================================
!  use some wrapper types to bundle information together
      type(TMolecule) :: mol
      type(scc_results) :: res
      class(TCalculator), allocatable :: calc, cpxcalc
      type(freq_results) :: fres
      type(TRestart) :: chk
      type(chrg_parameter) :: chrgeq
      type(TIFFData), allocatable :: iff_data
      type(oniom_input) :: oniom
      type(jab_input) :: dipro
      type(TCpcmx) :: cpx
      type(TTBLiteInput) :: tblite
      type(landscape_opt_settings) :: lset
      type(landscape_scan_settings) :: lsset
      type(landscape_vib_settings) :: lvset
      type(landscape_ts_settings) :: ltset
      type(landscape_irc_settings) :: lircset
      type(landscape_result) :: lres
!  store important names and stuff like that in FORTRAN strings
      character(len=:), allocatable :: fname    ! geometry input file
      character(len=:), allocatable :: missfile ! positional file the parser did not classify
      character(len=:), allocatable :: xcontrol ! instruction file
      character(len=:), allocatable :: xrc      ! global instruction file
      character(len=:), allocatable :: fnv      ! parameter file
      character(len=:), allocatable :: tmpname  ! temporary string
      character(len=:), allocatable :: cdum     ! temporary string
      character(len=:), allocatable :: extension, basename, directory
      integer :: ftype

!! ========================================================================
!  default names for important files in xtb
      character(len=*), parameter :: p_fname_rc = '.xtbrc'
      character(len=*), parameter :: p_fname_param_gfn0 = 'param_gfn0-xtb.txt'
      character(len=*), parameter :: p_fname_param_gfn1 = 'param_gfn1-xtb.txt'
      character(len=*), parameter :: p_fname_param_gfn2 = 'param_gfn2-xtb.txt'
      character(len=*), parameter :: p_fname_param_gfnff = '.param_gfnff.xtb'
      character(len=*), parameter :: p_fname_param_ipea = 'param_ipea-xtb.txt'
      character(len=*), parameter :: p_fname_param_ptb = 'param_ptb.txt'

      integer :: gsolvstate
      integer :: i, j, k, l, idum
      integer :: ich, ictrl, iprop ! file handle
      real(wp) :: sigma(3, 3)
      real(wp), allocatable :: cn(:)
      real(wp), allocatable :: sat(:)
      real(wp), allocatable :: g(:, :)
      real(wp), allocatable :: fukui(:, :)
      real(wp) :: vec3(3)
      type(TxTBParameter) :: globpar
      real(wp), allocatable :: dcn(:, :, :)
      real(wp), allocatable :: dq(:, :, :)
      real(wp), allocatable :: dumdumdum(:, :, :)
      real(wp), allocatable :: q(:)
      real(wp), allocatable :: ql(:)
      real(wp), allocatable :: qr(:)
!! ------------------------------------------------------------------------
      integer, external :: ncore

!! ------------------------------------------------------------------------
      logical :: struc_conversion_done = .false.
      logical :: anyopt, anyhess
      ! logical for checking whether ODLR approximation can be used
      logical :: odlr_valid

!! ========================================================================
!  debugging variables for numerical gradient
      logical, parameter :: gen_param = .false.
      logical, parameter :: debug = .false.
      type(TRestart) :: wf0
      real(wp), allocatable :: coord(:, :), numg(:, :), gdum(:, :)
      real(wp) :: sdum(3, 3), nums(3, 3), eps(3, 3), latt(3, 3)

      real(wp), parameter :: step = 0.00001_wp, step2 = 0.5_wp / step ! for numerical gradient
      real(wp), parameter :: sstep = 1.0_wp * 10.0_wp**(-6), sstep2 = 0.5_wp / sstep ! for numerical sigma
      real(wp) :: er, el
      logical :: coffee ! if debugging gets really though, get a coffee
      integer :: iter_needed = 0

!! ------------------------------------------------------------------------
!  undocumented and unexplainable variables go here
      integer :: nFiles, iFile
      integer :: rohf, err
      real(wp) :: dum5, egap, etot, ipeashift
      real(wp) :: zero, t0, t1, w0, w1, etot2, g298
      real(wp) :: one, two
      real(wp) :: ea, ip
      real(wp) :: vomega
      real(wp) :: energy_gas
      parameter(zero=0.0_wp)
      parameter(one=1.0_wp)
      parameter(two=2.0_wp)
      logical :: ex, okbas
      !> Set by `geometry_optimization` when the native minimisation did not
      !> converge; read by the `--metaopt` guard below.  Explicitly initialised:
      !> the native optimisation block is not the only reader of the variable
      !> in every run type, and an undefined value here decided whether a
      !> metadynamics run was aborted (review HIGH-1 of stage 3b).
      logical :: epr, diff, murks = .false.
      logical :: exist
      logical :: lgrad, restart
      logical :: copycontrol
      logical :: newreader
      logical :: strict
      logical :: exitRun
      logical :: cold_fusion
      !> Landscape run type for the `$scan` diagnostic (empty when the run type
      !> is not a landscape task)
      character(len=16) :: lruntyp_flag
      !> Native ANCopt engine selected by `--method` (`p_engine_*`), 0 for a
      !> landscape method name and when the flag was not given; set by
      !> `parseArguments` (CLI redesign sections 4.7/4.8, R30/MEDIUM-1)
      integer :: landscape_native_engine = 0
      !> `--acc` was named on the command line (`intent(out)` of
      !> `parseArguments`, same bridge as `landscape_native_engine`).  A
      !> landscape run derives the SCC accuracy of its engine from the
      !> convergence level (`set_conv_level` of `xtb_landscape_settings`), so
      !> the global `set%acc` may overwrite it only when the user asked for
      !> that value explicitly (CLI redesign section 3.5, R5).  Without the
      !> guard the level-derived `acc` was silently replaced by the global
      !> default 1.0 and the axis was a dead knob.  The value is passed
      !> between the flag parser and the six synchronisation points of
      !> `xtbMain` instead of living in a variable of the program unit, which
      !> the project hard constraint "no module-level mutable variables"
      !> forbids.
      logical :: acc_was_set

!  OMP stuff
      integer :: TID, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
      integer :: nproc

      type(TPrintTopo) :: printTopo ! gfnff topology printout list

      ! A Windows console decodes this program's UTF-8 labels with its own
      ! output code page (GBK/CP936 on a Chinese-locale system, which turns
      ! `·` into 路); set it to UTF-8 before anything is printed.  The call is
      ! absent on every other target, and its result is ignored: a console
      ! nicety must never fail a run.  See
      ! `docs/_simplify-recon/194-windows-encoding.md`.
#ifdef _WIN32
      call c_console_utf8()
#endif

      xenv%home = env%xtbhome
      xenv%path = env%xtbpath

      ! ------------------------------------------------------------------------
      ! read the command line arguments
      call parseArguments(env, argParser, xcontrol, fnv, lgrad, &
         & restart, gsolvstate, strict, copycontrol, coffee, printTopo, oniom, dipro, tblite, lset, &
         & lsset, lvset, ltset, lircset, landscape_native_engine, acc_was_set)


      ! No solvation available for PTB!
      if (set%mode_extrun == p_ext_ptb) then
         if (allocated(set%solvInput%solvent)) then
            call env%error("Solvation is not implemented for PTB", source)
         endif
      end if

      ! Spin-polarization is only available in the tblite library
      if (set%mode_extrun /= p_ext_tblite .and. tblite%spin_polarized) then
         call env%error("Spin-polarization is only available with the tblite library! Try --tblite", source)
      end if

      ! If hessian (or ohess or bhess) is requested in combination with PTB, conduct GFN2-xTB + PTB hessian
      anyhess = (set%runtyp == p_run_hess) .or. (set%runtyp == p_run_ohess) .or. (set%runtyp == p_run_bhess)
      if (anyhess) then
         ! O1NumHess only supports Hessian calculation, nothing else
         odlr_valid = freezeset%n == 0 .and. set%mode_extrun /= p_ext_ptb
         if (set%o1numhess) then
            if (.not. odlr_valid) then
               call env%error("O1NumHess does not support frozen atoms/PTB", source)
            end if
         else if(set%mode_extrun == p_ext_ptb) then
            set%mode_extrun = p_ext_xtb
            set%ptbsetup%ptb_in_hessian = .true.
            call set_gfn(env, 'method', '2')
            call set_gfn(env, 'd4', 'true')
            tblite%method = "gfn2"
            set%ptbsetup%hessmethod = "GFN2-xTB"
         else
            if (set%elprop == p_elprop_alpha) then
               call env%error("Raman activities are not implemented with this method",source)
            end if
         endif
      endif 

      nFiles = argParser%countFiles()
      select case (nFiles)
      case (0)
         if (.not. coffee) then
          if (printTopo%warning) call env%error("Eventually the input file was given to wrtopo as an argument.", source)
            call env%error("No input file given, so there is nothing to do", source)
         else
            fname = 'coffee'
         end if
      case (1:)
         do iFile = 1, nFiles - 1
            call argParser%nextFile(fname)
            call positional_argument_error(env, fname, set%runtyp, .true., source)
         end do
         call argParser%nextFile(fname)
         ! A positional argument that is neither an option nor an existing file
         ! is reported as well.  Both cases are **errors**, not warnings (owner
         ! ruling, `docs/_simplify-recon/197`): the run used to go ahead with
         ! whichever file the parser recognised, so `oxtb a.xyz b.xyz` silently
         ! read `b.xyz` and `oxtb water.xyz --lvib normal` - a forgotten level
         ! keyword - exited 0 and computed the frequencies of the *unrelaxed*
         ! input geometry.  A command line the program cannot honour stops the
         ! run instead of producing a result nobody asked for; the message names
         ! the token and, when it is a convergence level, the flag that takes it.
         do
            call argParser%nextMissingFile(missfile)
            if (.not. allocated(missfile)) exit
            call positional_argument_error(env, missfile, set%runtyp, .false., source)
         end do
      end select

      ! The accumulated positional errors abort here, before anything is read:
      ! a run with an unusable command line must not spend a single evaluation.
      call env%checkpoint("Unusable positional argument on the command line")

      !> `argParser%nextFile` declares its output as `intent(out)`, so it
      !> deallocates the argument and only allocates it again when it finds a
      !> file argument.  An argument the parser does not classify as a file (a
      !> path that does not exist, for example) therefore leaves `fname`
      !> unallocated, and assigning from an unallocated deferred-length
      !> character is undefined behaviour - it segfaults in a checked build
      !> (`-fcheck=all`, `-O0`) while an optimised build happens to survive.
      !> Give it a defined value so the error path below can report the real
      !> problem.
      if (.not. allocated(fname)) fname = ''

      if (.not. allocated(xcontrol)) then
         if (copycontrol) then
            xcontrol = 'xtb.inp'
         else
            xcontrol = fname
         end if
      end if

      anyopt = ((set%runtyp == p_run_opt) .or. (set%runtyp == p_run_ohess) .or. &
         &   (set%runtyp == p_run_omd) .or. (set%runtyp == p_run_screen) .or. &
         &   (set%runtyp == p_run_metaopt))

      if (allocated(set%solvInput%cpxsolvent).and.anyopt) then
            call env%terminate("CPCM-X not implemented for geometry optimization. &
                 & Please use another solvation model for optimization instead.")
      endif     
 
      if ((set%mode_extrun == p_ext_ptb) .and. anyopt) call env%terminate("PTB not implemented for geometry optimization. &
         &Please use another method for optimization instead.")

      call env%checkpoint("Command line argument parsing failed")

      ! ------------------------------------------------------------------------
      ! read the detailed input file
      call rdcontrol(xcontrol, env, copy_file=copycontrol)

      call env%checkpoint("Reading '"//xcontrol//"' failed")

      ! ------------------------------------------------------------------------
      ! read dot-Files before reading the rc and after reading the xcontrol
      ! Total molecular charge
      call open_file(ich, '.CHRG', 'r')
      if (ich /= -1) then
         call getline(ich, cdum, iostat=err)
         if (err /= 0) then
            call env%error('.CHRG is empty!', source)
         else
            call set_chrg(env, cdum)
            call close_file(ich)
         end if
      end if

      call env%checkpoint("Reading charge from file failed")

      ! Number of unpaired electrons
      call open_file(ich, '.UHF', 'r')
      if (ich /= -1) then
         call getline(ich, cdum, iostat=err)
         if (err /= 0) then
            call env%error('.UHF is empty!', source)
         else
            call set_spin(env, cdum)
            call close_file(ich)
         end if
      end if

      call env%checkpoint("Reading multiplicity from file failed")

      ! efield read: tblite or gfnff and PTB only
      if (set%mode_extrun == p_ext_gfnff .or. set%mode_extrun == p_ext_ptb &
         & .or. set%mode_extrun == p_ext_tblite) then
         call open_file(ich, '.EFIELD', 'r')
         if (ich /= -1) then
            call getline(ich, cdum, iostat=err)
            if (err /= 0) then
               call env%error('.EFIELD is empty!', source)
            else
               call set_efield(env, cdum)
               ! Take electric field from file also for tblite if not already present
               if (set%mode_extrun == p_ext_tblite .and. .not.allocated(tblite%efield)) then
                  allocate(tblite%efield(3))
                  tblite%efield = set%efield
               end if
               call close_file(ich)
            end if
         end if
      end if

      ! If EFIELD is not zero when using xtb, print a warning
      if (((set%mode_extrun /= p_ext_ptb) .and. (set%mode_extrun /= p_ext_gfnff) &
         & .and. (set%mode_extrun /= p_ext_tblite)) &
         & .and. (sum(abs(set%efield)) /= 0.0_wp)) then
         call env%terminate("External electric field is not zero ('--efield' or file '.EFIELD'), &
            & but only supported via tblite or for GFN-FF and PTB")
      end if

      ! ------------------------------------------------------------------------
      ! read the xtbrc if you can find it (use rdpath directly instead of xfind)
      call rdpath(env%xtbpath, p_fname_rc, xrc, exist)
      if (exist) then
         ! The fallback rc file is not the file the caller selected, so its
         ! `$opt` block must not mark the parameter state as a user request:
         ! the repository ships an `.xtbrc` containing `$opt`, and every suite
         ! test runs with `XTBPATH` pointing at the repository root (CLI
         ! redesign section 3.8(c), `docs/_simplify-recon/42-stage3a-impl.md`).
         call rdcontrol(xrc, env, copy_file=.false., from_input=.false.)

         call env%checkpoint("Reading '"//xrc//"' failed")
      end if

      ! ------------------------------------------------------------------------
      ! FIXME: some settings that are still not automatic
      ! Make sure GFN0-xTB uses the correct exttyp
      if (set%gfn_method == 0) call set_exttyp('eht')
      rohf = 1 ! HS default
      egap = 0.0_wp
      ipeashift = 0.0_wp

      ! ========================================================================
      ! no user interaction up to now, time to show off!
      ! print the xtb banner with version number and compilation date
      ! making a fancy version of this is hard, x is difficult in ASCII art
      call xtb_header(env%unit)
      ! make sure you cannot blame us for destroying your computer
      call disclaimer(env%unit)
      ! how to cite this program
      call citation(env%unit)
      ! print current time
      call prdate('S')

      ! ------------------------------------------------------------------------
      ! get molecular structure
      if (coffee) then ! it's coffee time
         fname = 'caffeine'
         call get_coffee(mol)
         call generateFileMetaInfo(fname, directory, basename, extension)
      else
         call generateFileMetaInfo(fname, directory, basename, extension)
         ftype = getFileType(fname)
         call open_file(ich, fname, 'r')
         call readMolecule(env, mol, ich, ftype)
         call close_file(ich)
         if (mol%info%two_dimensional) then
            call env%warning("Two dimensional input structure detected", source)
         end if

         ! Special CT input file case
         if (mol%chrg /= 0.0_wp) then
            if (set%clichrg) then
               call env%warning("Charge in sdf/mol input was overwritten", source)
            else
               set%ichrg = nint(mol%chrg)
            end if
         end if

         call env%checkpoint("reading geometry input '"//fname//"' failed")
      end if

      if (mol%n == 1 .and. set%o1numhess) then
         call env%warning("O1NumHess not supported for single atoms. Using default semi-numerical Hessian.", source)
         set%o1numhess = .false.
      end if

      ! ------------------------------------------------------------------------
      ! initialize the global storage
      call init_fix(mol%n)
      call init_split(mol%n)
      call init_constr(mol%n, mol%at)
      call init_scan
      call init_walls
      call init_pcem
      if (set%runtyp == p_run_bhess) then
         call init_bhess(mol%n)
      else
         call init_metadyn(mol%n, metaset%maxsave)
      end if
      ! Initialize the atomic masses with the physical constants
      atmass = atomic_mass(mol%at) * autoamu ! from splitparam.f90
      call load_rmsdbias(rmsdset, mol%n, mol%at, mol%xyz)
      ! ------------------------------------------------------------------------
      ! CONSTRAINTS & SCANS
      ! now we are at a point that we can check for requested constraints
      call read_userdata(xcontrol, env, mol)

      ! ------------------------------------------------------------------------
      ! the xcontrol blocks the landscape suite does not execute
      !
      ! One drop point for the three guards of the CLI redesign (section 3.8).
      ! All three read `landscape_is_suite_run`, the single criterion that
      ! separates a suite run from a native one, so a native exit is never
      ! refused by mistake: `--opt --method ancopt|lancopt|fire|pbc_lbfgs` and
      ! `--ohess --method ancopt|...` keep the native xcontrol handling.  The
      ! block sits after `readMolecule` (`mol%npbc` is known), `rdcontrol`
      ! (`set%xcontrol_supplied` is known) and `read_userdata` (`nscan` is
      ! known), and before the native optimisation block, so no run is stopped
      ! after it has started working.
      !
      ! (a) `$scan` belongs to the ANCopt family: `relaxed_scan` executes it
      ! only on the `anyopt` path of this program (the `if (nscan > 0)` guard
      ! inside the `if (anyopt)` branch below).  No landscape run type executes
      ! it, and before this check a run such as `--lopt --input <file with
      ! $scan>` returned status 0 after a plain minimisation, wrote no
      ! `xtbscan.log` and printed no diagnostic: the caller asked for a scan
      ! and silently got a different task.  Since stage 3 the redeployed `--opt`
      ! is such a run type too, so the message points at the native exit
      ! (`--opt --method ancopt`) and not at `--opt` itself.
      !
      ! (b) the suite rejects periodic boundary conditions
      ! (`src/landscape/landscape_driver.f90`), but its message does not name
      ! the native exit; this guard is the one that does.
      !
      ! (c) the suite never reads the xcontrol `$opt` parameters, so a supplied
      ! `$opt` block is a dead knob for every run type that goes to the suite -
      ! refused instead of ignored (C7 of the redesign).
      if (landscape_is_suite_run(set%runtyp, landscape_native_engine)) then
         lruntyp_flag = ''
         select case (set%runtyp)
         case (p_run_opt)
            lruntyp_flag = '--opt'
         case (p_run_ohess)
            lruntyp_flag = '--ohess'
         case (p_run_lopt)
            lruntyp_flag = '--lopt'
         case (p_run_lscan)
            lruntyp_flag = '--lscan'
         case (p_run_lts)
            lruntyp_flag = '--lts'
         case (p_run_lirc)
            lruntyp_flag = '--lirc'
         case (p_run_lvib)
            lruntyp_flag = '--lvib'
         case default
            lruntyp_flag = ''
         end select
         if (nscan > 0) then
            call env%error(trim(lruntyp_flag)//" does not execute the xcontrol "//&
               & "$scan block: use --opt --method ancopt to run the native ANCopt "//&
               & "that performs it, or --lscan with --lscan-coord/--lscan-atoms "//&
               & "for a scan of the landscape suite", source)
         end if
         if (mol%npbc > 0) then
            call env%error(trim(lruntyp_flag)//" runs the landscape suite, which "//&
               & "does not support periodic boundary conditions: use --opt --method "//&
               & "ancopt for the native ANCopt (pbc_lbfgs / GFN-FF periodic) "//&
               & "optimisation", source)
         end if
         if (set%xcontrol_supplied) then
            call env%error("the xcontrol $opt block does nothing for "//&
               & trim(lruntyp_flag)//": the landscape suite does not read the "//&
               & "ANCopt parameters; remove the block, or use --opt --method "//&
               & "ancopt to run the native ANCopt (rf) that honours "//&
               & "engine=/optlevel=/hessian=/maxcycle=...", source)
         end if
      end if
      call env%checkpoint("the run type does not execute this xcontrol block")
      ! ------------------------------------------------------------------------
      ! get some memory
      allocate (cn(mol%n), sat(mol%n), g(3, mol%n), source=0.0_wp)
      set%periodic = mol%npbc > 0
      if (mol%npbc == 0) then
         if (set%do_cma_trafo) then
            allocate (coord(3, mol%n), source=0.0_wp)
            call axis3(1, mol%n, mol%at, mol%xyz, coord, vec3)
            mol%xyz = coord
            deallocate (coord)
         end if
      end if

      do i = 1, mol%n
         mol%z(i) = mol%at(i) - ncore(mol%at(i))
         ! lanthanides without f are treated as La
         if (mol%at(i) > 57 .and. mol%at(i) < 72) mol%z(i) = 3
      end do

      ! initialize time step for MD if requested autocomplete
      if (set%tstep_md < 0.0_wp) then
         set%tstep_md = (minval(atmass) / (atomic_mass(1) * autoamu))**(1.0_wp / 3.0_wp)
      end if

      mol%chrg = real(set%ichrg, wp)
   !! To assign charge
      mol%uhf = set%nalphabeta
      call initrand

      call setup_summary(env%unit, mol%n, fname, xcontrol, chk%wfn, xrc)

      ! The positional input file reaches the `provenance` block of the
      ! `--ljson` document through the settings objects, the only channel to
      ! the tasks.  It is assigned after the geometry was read, so a default
      ! input (`caffeine`) or a converted one is reported as the file the run
      ! actually used.  The suite settings are unused by a native run type and
      ! vice versa; both carry the name so the document is uniform.
      lset%input_file = fname
      lsset%input_file = fname
      lvset%input_file = fname
      ltset%input_file = fname
      lircset%input_file = fname

      ! Announce the engine that will drive the run.  It has to happen here:
      ! after every flag was parsed (so a default and an explicit engine are
      ! known) and before the first geometry is evaluated - the 2D => 3D
      ! converter below is already one, and the always-done single point is
      ! another.  The per-task echo this replaces ran after that single point.
      call landscape_engine_banner(env%unit, set%runtyp, &
         & landscape_native_engine, lset, lsset, ltset, lircset)

      ! ------------------------------------------------------------------------
      ! 2D => 3D STRUCTURE CONVERTER
      ! ------------------------------------------------------------------------
      if (mol%info%two_dimensional) then
         call struc_convert(env, restart, mol, chk, egap, set%etemp, set%maxscciter, &
                           &  set%optset%maxoptcycle, etot, g, sigma)
         struc_conversion_done = .true.
         mol%info%two_dimensional = .false.
      end if


      if (sum(abs(set%efield)) /= 0.0_wp) then
         write (env%unit, '(/,3x,a)') "--------------------------------------"
         write (env%unit, '(3x,a)') "--- external electric field / a.u. ---"
         write (env%unit, '(3x,3(a,f8.4))') "x = ", set%efield(1), " y = ", set%efield(2), " z = ", set%efield(3)
         write (env%unit, '(3x,a)') "--------------------------------------"
      end if

      ! initialize metadynamics
      call load_metadynamic(metaset, mol%n, mol%at, mol%xyz)

      ! restraining potential
      if (allocated(potset%xyz)) then
         if (lconstr_all_bonds) call constrain_all_bonds(mol%n, mol%at, potset%xyz)
         if (lconstr_all_angles) call constrain_all_angles(mol%n, mol%at, potset%xyz)
         if (lconstr_all_torsions) call constrain_all_torsions(mol%n, mol%at, potset%xyz)
         call setup_constrain_pot(mol%n, mol%at, potset%xyz)
      else
         if (lconstr_all_bonds) call constrain_all_bonds(mol%n, mol%at, mol%xyz)
         if (lconstr_all_angles) call constrain_all_angles(mol%n, mol%at, mol%xyz)
         if (lconstr_all_torsions) call constrain_all_torsions(mol%n, mol%at, mol%xyz)
         call setup_constrain_pot(mol%n, mol%at, mol%xyz)
      end if
      !  fragmentation for CMA constrain
      if (iatf1 == 0 .and. iatf2 == 0) then
         call ncoord_erf(mol%n, mol%at, mol%xyz, cn)
         call splitm(mol%n, mol%at, mol%xyz, cn)
      end if
      call splitprint(mol%n, mol%at, mol%xyz)

      if (set%verbose) then
         call fix_info(env%unit, mol%n, mol%at, mol%xyz)
         call pot_info(env%unit, mol%n, mol%at, mol%xyz)
      end if

      ! ------------------------------------------------------------------------
      ! write copy of detailed input
      if (copycontrol) then
         call open_set(ictrl, xcontrol)
         call write_set(ictrl)
         call close_set(ictrl)
      end if

      ! ------------------------------------------------------------------------
      ! if you have requested a define we stop here...
      if (set%define) then
         if (set%verbose) call main_geometry(env%unit, mol)
         call eval_define(set%veryverbose)
      end if
      call env%show('Please study the warnings concerning your input carefully')
      call raise('F', 'Please study the warnings concerning your input carefully')

      ! ========================================================================
      ! From here we switch to the method setup
      ! enable error on warnings
      if (strict) call mctc_strict
      env%strict = strict

      ! one last check on the input geometry
      call check_cold_fusion(env, mol, cold_fusion)
      if (cold_fusion) then
         call env%error("XTB REFUSES TO CONTINUE WITH THIS CALCULATION!")
         call env%terminate("Some atoms in the start geometry are *very* close")
      end if

      !> check if someone is still using GFN3...
      ! TODO: this is probably a bit outdated now
      if (set%gfn_method == 3) then
         call env%terminate('Wait for some months - for now, please use gfn_method=2!')
      end if

      ! ------------------------------------------------------------------------
      ! Print the method header and select the parameter file
      if (.not. allocated(fnv)) then
         select case (set%runtyp)
         case default
            call env%terminate('This is an internal error, please define your runtypes!')
         case (p_run_prescc,p_run_scc, p_run_grad, p_run_opt, p_run_hess, p_run_ohess, p_run_bhess, &
               p_run_md, p_run_omd, p_run_path, p_run_screen, &
               p_run_modef, p_run_mdopt, p_run_metaopt, p_run_lopt, p_run_lscan, p_run_lts, p_run_lvib, p_run_lirc)
            if (set%mode_extrun == p_ext_gfnff) then
               fnv = xfind(p_fname_param_gfnff)
            elseif (set%mode_extrun .eq. p_ext_ptb) then
                  fnv = xfind(p_fname_param_ptb)
            else
               if (set%gfn_method == 0) then
                  fnv = xfind(p_fname_param_gfn0)
               end if
               if (set%gfn_method == 1) then
                  fnv = xfind(p_fname_param_gfn1)
               end if
               if (set%gfn_method == 2) then
                  fnv = xfind(p_fname_param_gfn2)
               end if
            end if
         case (p_run_vip, p_run_vea, p_run_vipea, p_run_vfukui, p_run_vomega)
            if (set%mode_extrun .eq. p_ext_ptb) then
                  fnv = xfind(p_fname_param_ptb)
            else
               if (set%gfn_method == 0) then
                  fnv = xfind(p_fname_param_gfn0)
               end if
               if (set%gfn_method == 1) then
                  fnv = xfind(p_fname_param_gfn1)
               end if
               if (set%gfn_method == 2) then
                  fnv = xfind(p_fname_param_gfn2)
               end if
            end if
         end select
      end if

      !-------------------------------------------------------------------------
      ! Perform a precomputation of electronic properties for xTB-IFF
      if (set%mode_extrun == p_ext_iff) then
         allocate (iff_data)
         call prepare_IFF(env, mol, iff_data)
         call env%checkpoint("Could not generate electronic properties")
      end if

      ! ------------------------------------------------------------------------
      ! Obtain the parameter data
      call newCalculator(env, mol, calc, fnv, restart, set%acc, oniom, iff_data, tblite)
      call env%checkpoint("Could not setup single-point calculator")

      call initDefaults(env, calc, mol, gsolvstate)
      call env%checkpoint("Could not setup defaults")

      ! ------------------------------------------------------------------------
      ! initial guess, setup wavefunction
      select type (calc)
      type is (TxTBCalculator)
         call chk%wfn%allocate(mol%n, calc%basis%nshell, calc%basis%nao)

         ! Make sure number of electrons is initialized and multiplicity is consistent
         ! Default is restricted closed-shell singlet for even number of electrons
         ! and restricted open-shell doublet for odd number of electrons
         chk%wfn%nel = nint(sum(mol%z) - mol%chrg)
         if (mod(mol%uhf, 2) == mod(chk%wfn%nel, 2)) then
            ! Restricted closed-shell case
            chk%wfn%nopen = mol%uhf
         else
            ! Restricted open-shell case
            if (mol%uhf /= 0) then
               call env%terminate("Assigned number of unpaired electrons (flag '--uhf <int>' or <int> in file '.UHF') is not consistent with the total number of electrons")
            end if
            chk%wfn%nopen = mod(int(chk%wfn%nel), 2)
         end if 

         ! EN charges and CN
         if (set%gfn_method < 2) then
            call ncoord_d3(mol%n, mol%at, mol%xyz, cn)
         else
            call ncoord_gfn(mol%n, mol%at, mol%xyz, cn)
         end if
         if (mol%npbc > 0) then
            chk%wfn%q = real(set%ichrg, wp) / real(mol%n, wp)
         else
            if (set%guess_charges == p_guess_gasteiger) then
               call iniqcn(mol%n, mol%at, mol%z, mol%xyz, set%ichrg, 1.0_wp, chk%wfn%q, cn, set%gfn_method, .true.)
            else if (set%guess_charges == p_guess_goedecker) then
               call ncoord_erf(mol%n, mol%at, mol%xyz, cn)
               call goedecker_chrgeq(mol%n, mol%at, mol%xyz, real(set%ichrg, wp), cn, dcn, chk%wfn%q, dq, er, g, &
                                     .false., .false., .false.)
            else
               call ncoord_gfn(mol%n, mol%at, mol%xyz, cn)
               chk%wfn%q = real(set%ichrg, wp) / real(mol%n, wp)
            end if
         end if
         ! initialize shell charges from gasteiger charges
         call iniqshell(calc%xtbData, mol%n, mol%at, mol%z, calc%basis%nshell, chk%wfn%q, chk%wfn%qsh, set%gfn_method)
      type is (TTBLiteCalculator)
         call newTBLiteWavefunction(env, mol, calc, chk)
      end select

      ! get CEH charges !
      if (allocated(tblite%ceh)) then
         if (lset%ljson) call env%warning("--ljson writes nothing: the CEH "//&
            & "pre-calculation returns before the document stage", source)
         ! if numercal charges are requested !
         if (tblite%ceh%grad) then
            call num_grad_chrg(env,mol,tblite)
         ! only charges !
         else
            call get_ceh(env,mol,tblite)
         endif 

         ! stop calculation !
         call finalize_xtb(env, landscape_native_engine)
      end if

      ! ------------------------------------------------------------------------
      !> printout a header for the exttyp
      call calc%writeInfo(env%unit, mol)

      call delete_file('.sccnotconverged')

      call env%checkpoint("Setup for calculation failed")

      exist = .false.
      select type (calc)
      type is (TxTBCalculator)
         if (restart .and. calc%xtbData%level /= 0) then ! only in first run
            call readRestart(env, chk%wfn, 'xtbrestart', mol%n, mol%at, set%gfn_method, exist, .true.)
         end if
         ipeashift = calc%xtbData%ipeashift
      type is (TTBLiteCalculator)
         if (restart) then
            call loadRestart(env, chk, 'xtbrestart', exist)
            if (exist) write (env%unit, "(a)") "Wavefunction read from restart file"
         end if
         ipeashift = calc%ipeashift
      type is (TOniomCalculator)
         select type (xtb => calc%real_low)
         type is (TxTBCalculator)
            call chk%wfn%allocate(mol%n, xtb%basis%nshell, xtb%basis%nao)
            call newWavefunction(env, mol, xtb, chk)
         !! assigns only partial charges q and shell charges
            if (restart) then ! only in first run
               call readRestart(env, chk%wfn, 'xtbrestart', mol%n, mol%at, set%gfn_method, exist, .true.)
            end if
         end select
         if (.not. set%oniom_settings%fixed_chrgs) then
            set%oniom_settings%innerchrg = calculateCharge(calc, env, mol, chk)
         end if
         if (.not. set%oniom_settings%fixed_spin) then
            set%oniom_settings%innerspin = mol%uhf
         end if
      end select
      !-------------------------------------------------------------------------
      ! DIPRO calculation of coupling integrals for dimers
      if (dipro%diprocalc) then
         if (lset%ljson) call env%warning("--ljson writes nothing: the DIPRO "//&
            & "coupling calculation terminates before the document stage", source)
         call start_timing(11)
         call get_jab(env, tblite, mol, splitlist, dipro)
         call env%checkpoint("Something in your DIPRO calculation went wrong.")
         call stop_timing_run
         call stop_timing(11)
         write (*, '(A)') "----------------------------------------------------------"
         call prdate('E')
         write (*, '(A)') "----------------------------------------------------------"
         call prtiming(11, 'dipro')
         call terminate(0)
      end if

      ! ========================================================================
      ! the SP energy which is always done
      call start_timing(2)
      call calc%singlepoint(env, mol, chk, 2, exist, etot, g, sigma, egap, res)
      call stop_timing(2)
      select type (calc)
      type is (TGFFCalculator)
         gff_print = .false.
      type is (TxTBCalculator)
         if (restart) call writeRestart(env, chk%wfn, 'xtbrestart', set%gfn_method)
      type is (TTBLiteCalculator)
         if (restart) call dumpRestart(env, chk, 'xtbrestart')
      end select
      call env%checkpoint("Single point calculation terminated")

      ! write 2d => 3d converted structure
      if (struc_conversion_done) then
         call generateFileName(tmpname, 'gfnff_convert', extension, mol%ftype)
         write (env%unit, '(10x,a,1x,a,/)') &
            "converted geometry written to:", tmpname
         call open_file(ich, tmpname, 'w')
         call writeMolecule(mol, ich, energy=res%e_total, gnorm=res%gnorm)
         call close_file(ich)
      end if

      ! ========================================================================
      ! determine kopt for bhess including final biased geometry optimization
      if (set%runtyp == p_run_bhess) then
         call set_metadynamic(metaset, mol%n, mol%at, mol%xyz)
         call get_kopt(metaset, env, restart, mol, chk, calc, egap, set%etemp, set%maxscciter, &
            & set%optset%maxoptcycle, set%optset%optlev, etot, g, sigma, set%acc)
      end if

      ! ------------------------------------------------------------------------
      ! numerical gradient for debugging purposes
      if (debug) then
         !  generate a warning to keep release versions from calculating numerical gradients
         call env%warning('XTB IS CALCULATING NUMERICAL GRADIENTS, RESET DEBUG FOR RELEASE!')
         print'(/,"analytical gradient")'
         print *, g
         allocate (coord(3, mol%n), source=mol%xyz)
         allocate (numg(3, mol%n), gdum(3, mol%n), source=0.0_wp)
         wf0 = chk
         do i = 1, mol%n
            do j = 1, 3
               mol%xyz(j, i) = mol%xyz(j, i) + step
               chk = wf0
               call calc%singlepoint(env, mol, chk, 0, .true., er, gdum, sdum, egap, res)
               mol%xyz(j, i) = mol%xyz(j, i) - 2 * step
               chk = wf0
               call calc%singlepoint(env, mol, chk, 0, .true., el, gdum, sdum, egap, res)
               mol%xyz(j, i) = mol%xyz(j, i) + step
               numg(j, i) = step2 * (er - el)
            end do
         end do
         print'(/,"numerical gradient")'
         print *, numg
         print'(/,"difference gradient")'
         print *, g - numg
         deallocate (coord)
      end if

      ! numerical sigma (=volume*stressTensor) for debugging purposes
      if (debug .and. mol%npbc == 3) then
         !  generate a warning to keep release versions from calculating numerical gradients
         call env%warning('XTB IS CALCULATING NUMERICAL STRESS, RESET DEBUG FOR RELEASE!')
         print'(/,"analytical sigma")'
         print *, sigma
         if (.not. allocated(gdum)) allocate (gdum(3, mol%n), source=0.0_wp)
         allocate (coord(3, mol%n), source=mol%xyz)
         latt = mol%lattice
         nums = 0.0_wp
         !sdum = 0.0_wp
         wf0 = chk
         do j = 1, 3
            do i = 1, j
               ! Only eps_ij=step the rest equals zero: eps_(kl.ne.ij)=0
               eps = 0.0_wp
               eps(i, j) = sstep
               ! adjust position vectors and lattice to get er
               do k = 1, 3
                  mol%xyz(k, :) = kron(k, 1) * mol%xyz(1, :) + eps(k, 1) * mol%xyz(1, :) + &
                                &   kron(k, 2) * mol%xyz(2, :) + eps(k, 2) * mol%xyz(2, :) + &
                                &   kron(k, 3) * mol%xyz(3, :) + eps(k, 3) * mol%xyz(3, :)
                  mol%lattice(k, 1) = (kron(k, 1) + eps(k, 1)) * mol%lattice(1, 1) + &
                                   & (kron(k, 2) + eps(k, 2)) * mol%lattice(2, 1) + &
                                   & (kron(k, 3) + eps(k, 3)) * mol%lattice(3, 1)
                  mol%lattice(k, 2) = (kron(k, 1) + eps(k, 1)) * mol%lattice(1, 2) + &
                                   & (kron(k, 2) + eps(k, 2)) * mol%lattice(2, 2) + &
                                   & (kron(k, 3) + eps(k, 3)) * mol%lattice(3, 2)
                  mol%lattice(k, 3) = (kron(k, 1) + eps(k, 1)) * mol%lattice(1, 3) + &
                                   & (kron(k, 2) + eps(k, 2)) * mol%lattice(2, 3) + &
                                   & (kron(k, 3) + eps(k, 3)) * mol%lattice(3, 3)
               end do
               chk = wf0
               call calc%singlepoint(env, mol, chk, 0, .true., er, gdum, sdum, egap, res)

               ! reset coordinates and lattice
               mol%xyz = coord
               mol%lattice = latt
               ! adjust position vectors and lattice to get el
               do k = 1, 3
                  mol%xyz(k, :) = kron(k, 1) * mol%xyz(1, :) - eps(k, 1) * mol%xyz(1, :) + &
                                &  kron(k, 2) * mol%xyz(2, :) - eps(k, 2) * mol%xyz(2, :) + &
                                &  kron(k, 3) * mol%xyz(3, :) - eps(k, 3) * mol%xyz(3, :)
                  mol%lattice(k, 1) = (kron(k, 1) - eps(k, 1)) * mol%lattice(1, 1) + &
                                   & (kron(k, 2) - eps(k, 2)) * mol%lattice(2, 1) + &
                                   & (kron(k, 3) - eps(k, 3)) * mol%lattice(3, 1)
                  mol%lattice(k, 2) = (kron(k, 1) - eps(k, 1)) * mol%lattice(1, 2) + &
                                   & (kron(k, 2) - eps(k, 2)) * mol%lattice(2, 2) + &
                                   & (kron(k, 3) - eps(k, 3)) * mol%lattice(3, 2)
                  mol%lattice(k, 3) = (kron(k, 1) - eps(k, 1)) * mol%lattice(1, 3) + &
                                   & (kron(k, 2) - eps(k, 2)) * mol%lattice(2, 3) + &
                                   & (kron(k, 3) - eps(k, 3)) * mol%lattice(3, 3)
               end do
               chk = wf0
               call calc%singlepoint(env, mol, chk, 0, .true., el, gdum, sdum, egap, res)

               ! numerical sigma (=volume*stressTensor)
               nums(i, j) = sstep2 * (er - el)  ! divide by 2 times step size
               nums(j, i) = nums(i, j)  ! stress tensor is symmetric
               ! reset coordinates and lattice
               mol%xyz = coord
               mol%lattice = latt
            end do
         end do

         print'(/,"numerical sigma")'
         print *, nums
         print'(/,"difference sigma")'
         print *, sigma - nums
         deallocate (coord)
      end if

!---------------------------------------------!
! Geometry optimization(ANCopt,L_ANCopt,FIRE) !
!---------------------------------------------!
      ! Native exit of the `anyopt` family only: the redeployed `--opt`/`--ohess`
      ! keep their run type (`p_run_opt`/`p_run_ohess`), so `anyopt` is still
      ! true for them and this block would run a second, native optimisation
      ! next to the suite task (CLI redesign R24).  The two `anyopt`
      ! `env%terminate` guards above (CPCM-X, PTB) keep their behaviour on
      ! purpose: a native workaround is not available for those two either way.
      !
      ! The criterion is `landscape_is_suite_run`, not `landscape_native_engine
      ! /= 0`: `anyopt` also covers `--omd`/`--screen`/`--metaopt`, which are
      ! never redeployed and cannot take `--method` at all, so their
      ! `landscape_native_engine` is always 0 and an `/= 0` test would make this
      ! block unreachable for them -- `--omd` would run its dynamics on the
      ! input geometry (review HIGH-1 of stage 3b, regression of stage 3a).
      if (anyopt .and. .not. landscape_is_suite_run(set%runtyp, landscape_native_engine)) then

         if (mol%npbc > 0 .and. (set%mode_extrun == p_ext_gfnff &
                 & .or. set%mode_extrun == p_ext_mcgfnff)) then  ! if(npbc)
            deallocate (set%opt_engine)
            call set_opt(env, 'engine', 'pbc_lbfgs')  ! use lbfgs
         end if

         ! start optimization timer !
         call start_timing(3)

         ! calculation !
         call geometry_optimization &
            &     (env, mol, chk, calc, &
        &      egap,set%etemp,set%maxscciter,set%optset%maxoptcycle,etot,g,sigma,set%optset%optlev,.true.,.false.,murks, iter_needed)

         ! save results !
         res%e_total = etot
         res%gnorm = norm2(g)

         ! constrained optimization !
         if (nscan > 0) then
            call relaxed_scan(env, mol, chk, calc)
         end if

         ! in case of failure cretae xtblast geometry !
         if (murks) then
            call generateFileName(tmpname, 'xtblast', extension, mol%ftype)
            write (env%unit, '(/,a,1x,a,/)') &
               "last geometry written to:", tmpname
            call open_file(ich, tmpname, 'w')
            call writeMolecule(mol, ich, energy=res%e_total, gnorm=res%gnorm)
            call close_file(ich)
            call env%terminate("Geometry optimization failed")
         end if

         ! stop optimization timer !
         call stop_timing(3)

      end if

      ! ------------------------------------------------------------------------
      ! landscape suite (minimisation with the algorithms of src/landscape)
      !
      ! Since stage 3 the redeployed `--opt` and `--ohess` run the suite as
      ! well; `effective_landscape_task` is the single place that says so, the
      ! run types themselves stay `p_run_opt`/`p_run_ohess` (CLI redesign
      ! section 3.11).
      !
      ! `--ljson` is an output switch, not a run type.  Since the unified
      ! document of report 117 **every** run type writes a
      ! `landscape_<task>.json`: the suite tasks through their own drivers,
      ! the native run types through `landscape_json_write_run` at the common
      ! post-processing.  The old warning ("writes nothing without a landscape
      ! run type") is therefore obsolete and removed; it used to fire for the
      ! native default `--sp`, which is exactly the case that now writes
      ! `landscape_sp.json`.  The only run types that never reach the common
      ! post-processing are the tblite-gated DIPRO and CEH short-circuits,
      ! which report their own reason (see `src/prog/main.F90`).

      ! BC-5 of the CLI redesign: `--opt`/`--ohess` change their meaning, and
      ! that is announced on every run.  The notice goes to `env%unit` (stdout,
      ! like `[WARNING]`), is printed once per run and is not suppressed by
      ! `--silent`.  It deliberately does not name an algorithm: `--opt --method
      ! bfgs` is the suite with `bfgs`, so "lbfgs" would be wrong there.
      if (landscape_is_suite_run(set%runtyp, landscape_native_engine) .and. &
         & (set%runtyp == p_run_opt .or. set%runtyp == p_run_ohess)) then
         select case (set%runtyp)
         case (p_run_ohess)
            write(env%unit, '(a)') "[NOTICE] --ohess now runs the landscape "//&
               & "suite (optimisation followed by harmonic frequencies); pass "//&
               & "--method ancopt for the legacy ANCopt optimiser."
         case default
            write(env%unit, '(a)') "[NOTICE] --opt now runs the landscape "//&
               & "suite optimiser; pass --method ancopt for the legacy "//&
               & "ANCopt optimiser."
         end select
      end if

      if (effective_landscape_task(set%runtyp, landscape_native_engine) == p_run_lopt) then
         ! The CLI electronic-structure settings reach the engine only through
         ! its own control object: `xtb_single::singlepoint` forwards `acc` as
         ! an argument but drops it, and the engine writes both values into the
         ! calculator, so the CLI values have to be synced here
         ! (docs/_landscape-recon/12-*.md, R2-02;
         ! docs/_landscape-recon/23-m3-fix-report.md, MEDIUM-B).  The accuracy
         ! is synced **only** when `--acc` named it: a convergence level
         ! derives `engine%acc` from the upstream table (CLI redesign section
         ! 3.5, R5), and an unconditional copy of the global default 1.0
         ! turned that axis into a dead knob (`--lopt --lopt-level loose`
         ! reported 1.0 instead of 2.0).
         if (acc_was_set) lset%engine%acc = set%acc
         lset%engine%etemp = set%etemp
         ! `max_steps` of the run: the level multiplier (`normal` = 3) is the
         ! default, an explicit `--lopt-max-steps`/`--cycles` keeps its value.
         call lset%level_step_budget(mol%n)
         call start_timing(3)
         call run_landscape(env, mol, chk, calc, lset, lres)
         call stop_timing(3)
         ! Report the *converged* geometry: `res` feeds the TOTAL ENERGY line
         ! (`write_energy`) and xtbout.json, while `xtbopt.xyz` is written by
         ! the landscape task from `lres%energy`.  Without this write-back the
         ! reported energy is the one of the starting geometry.
         res%e_total = lres%energy
         res%gnorm = lres%gnorm

         ! The redeployed `--ohess` is the minimisation followed by the
         ! frequency analysis of the minimised geometry, i.e. the two calls the
         ! native `--ohess` makes (`main_geometry` + `main_freq`,
         ! `src/prog/main.F90`).  The native Hessian and frequency blocks are
         ! suppressed for it through `landscape_native_hessian_path`, so this
         ! chain is the only one that runs (CLI redesign section 3.11.1).
         if (landscape_ohess_chain(set%runtyp, landscape_native_engine)) then
            ! T2 of CLI redesign section 3.9: the frequency half runs only on a
            ! converged geometry.  A failed minimisation has already logged its
            ! single error (the `landscape_optimize` task); running the
            ! frequency analysis anyway would reuse the calculator the failed
            ! optimisation left behind and produce a misleading cascade
            ! (`Could not generate lattice points` / `Setup of Coulomb
            ! evaluator failed`).  The native `--ohess` continues after a
            ! failed geometry optimisation and still writes the frequency
            ! files; aborting cleanly here is a deliberate deviation (recorded
            ! in the design document).
            if (lres%converged) then
            if (acc_was_set) lvset%engine%acc = set%acc
            lvset%engine%etemp = set%etemp
            call run_landscape(env, mol, chk, calc, lvset, lres)
            call env%checkpoint("Landscape task terminated")
            ! Bridge the suite frequencies into the native frequency container.
            ! `anyhess` stays the criterion of the three readers (the property
            ! file, `xtbout.json` and the stdout `write_energy`), because the
            ! container is what they print - not the path that filled it.  The
            ! alternative, setting `anyhess` to false for the suite branch,
            ! would silently drop the `TOTAL ENTHALPY`/`TOTAL FREE ENERGY`
            ! lines and the thermochemistry keys of `xtbout.json`, and a half
            ! fix (container printed, not filled) would print
            ! `TOTAL ENERGY 0.000000000000` instead of failing.  Filling it
            ! structurally removes that state.
            call fres%deallocate()
            call fres%allocate(mol%n)
            fres%n3true = 3*mol%n
            fres%linear = lres%vib%linear
            fres%freq = lres%vib%freq
            fres%rmass = lres%vib%rmass
            fres%hess = lres%vib%modes
            ! no dipole or polarisability derivatives come from a Hessian
            ! provider, so the intensity columns are zero, exactly like
            ! `--lvib` (see `src/landscape/vib/thermo.f90`)
            fres%dipt = 0.0_wp
            fres%polt = 0.0_wp
            fres%etot = lres%vib%energy
            fres%gnorm = lres%gnorm
            ! `write_energy` prints `etot + htot`, and the suite fields hold
            ! the totals H(T)/G(T), so the thermal contribution is the
            ! difference (`src/main/property.F90`)
            fres%htot = lres%vib%enthalpy - lres%vib%energy
            fres%gtot = lres%vib%gibbs - lres%vib%energy
            fres%zp = lres%vib%zpe
            fres%temp = lres%vib%temperature
            fres%nimag = lres%vib%n_imag
            ! The frequency analysis evaluates the same geometry, so the
            ! reported energy must stay the converged one.
            res%e_total = lres%energy
            res%gnorm = lres%gnorm
            end if
         end if
         call env%checkpoint("Landscape task terminated")
      end if

      ! ------------------------------------------------------------------------
      ! landscape suite (scan with the algorithms of src/landscape)
      if (effective_landscape_task(set%runtyp, landscape_native_engine) == p_run_lscan) then
         if (acc_was_set) lsset%engine%acc = set%acc
         lsset%engine%etemp = set%etemp
         call lsset%level_step_budget(mol%n)
         call start_timing(3)
         call run_landscape(env, mol, chk, calc, lsset, lres)
         call stop_timing(3)
         ! Report the lowest scan point: `res` feeds the TOTAL ENERGY line and
         ! xtbout.json, while `xtbscan.log` is written by the scan task from
         ! the de-biased energies of all scan points.
         res%e_total = lres%energy
         res%gnorm = lres%gnorm
         call env%checkpoint("Landscape task terminated")
      end if

      ! ------------------------------------------------------------------------
      ! landscape suite (transition state search with src/landscape)
      if (effective_landscape_task(set%runtyp, landscape_native_engine) == p_run_lts) then
         ! M4 hard constraint 1: the electronic temperature is frozen for the
         ! whole run and recorded.  The settings carry it explicitly (rather
         ! than reading the inherited engine control) because the task copies
         ! it into the engine control it hands to the engine.
         ltset%etemp = set%etemp
         if (acc_was_set) ltset%engine%acc = set%acc
         call ltset%level_step_budget(mol%n)
         ! The band searches need the name of the positional input file: they
         ! use it only to detect and reject a multi-frame input (the reactant
         ! is the geometry already read, the product comes from
         ! `--lts-product`), see `landscape_path_settings`.
         ltset%path%input_file = fname
         ! `--lts-max-steps` is the iteration budget of the run; the band
         ! searches use the same budget, so the option does not have to be
         ! duplicated as a `--lts-path-max-steps`.
         ltset%path%max_steps = ltset%max_steps
         call start_timing(3)
         call run_landscape(env, mol, chk, calc, ltset, lres)
         call stop_timing(3)
         ! Report the converged stationary point: `res` feeds the TOTAL
         ! ENERGY line (`write_energy`) and xtbout.json, while `xtbts.xyz` is
         ! written by the landscape task from `lres%energy`.  Without this
         ! write-back the reported energy is the one of the starting geometry.
         res%e_total = lres%energy
         res%gnorm = lres%gnorm
         call env%checkpoint("Landscape task terminated")
      end if

      ! ------------------------------------------------------------------------
      ! landscape suite (intrinsic reaction coordinate of src/landscape)
      if (effective_landscape_task(set%runtyp, landscape_native_engine) == p_run_lirc) then
         ! M4/M6 hard constraint: the electronic temperature is frozen for the
         ! whole run and recorded.  The task copies it into the engine control
         ! it hands to the engine.
         lircset%etemp = set%etemp
         if (acc_was_set) lircset%engine%acc = set%acc
         call lircset%level_step_budget(mol%n)
         call start_timing(3)
         call run_landscape(env, mol, chk, calc, lircset, lres)
         call stop_timing(3)
         ! Report the transition state energy (the reference of the IRC):
         ! `res` feeds the TOTAL ENERGY line (`write_energy`) and xtbout.json,
         ! while the four `irc_*` files are written by the landscape task.
         res%e_total = lres%energy
         res%gnorm = lres%gnorm
         call env%checkpoint("Landscape task terminated")
      end if

      ! ------------------------------------------------------------------------
      ! landscape suite (harmonic frequencies and thermochemistry of src/landscape)
      if (effective_landscape_task(set%runtyp, landscape_native_engine) == p_run_lvib) then
         if (acc_was_set) lvset%engine%acc = set%acc
         lvset%engine%etemp = set%etemp
         call start_timing(3)
         call run_landscape(env, mol, chk, calc, lvset, lres)
         call stop_timing(3)
         ! Report the energy of the analysed geometry: `res` feeds the TOTAL
         ! ENERGY line (`write_energy`) and xtbout.json.  Without this
         ! write-back the reported energy would be the one of an unrelated
         ! state (the same defect the `--lopt` path had).
         res%e_total = lres%energy
         res%gnorm = lres%gnorm
         call env%checkpoint("Landscape task terminated")
      end if

      ! ------------------------------------------------------------------------
      ! automatic VIP and VEA single point (maybe after optimization)
      if (set%runtyp == p_run_vip .or. set%runtyp == p_run_vipea &
         & .or. set%runtyp == p_run_vomega) then
         call start_timing(2)
         call vip_header(env%unit)
         mol%chrg = mol%chrg + 1
         chk%wfn%nel = chk%wfn%nel - 1
         if (mod(chk%wfn%nel, 2) /= 0) chk%wfn%nopen = 1
         call calc%singlepoint(env, mol, chk, 1, exist, etot2, g, sigma, egap, res)
         ip = etot2 - etot - ipeashift
         write (env%unit, '(72("-"))')
         write (env%unit, '("empirical IP shift (eV):",f10.4)') &
            &                  autoev * ipeashift
         write (env%unit, '("delta SCC IP (eV):",f10.4)') autoev * ip
         write (env%unit, '(72("-"))')
         mol%chrg = mol%chrg - 1
         chk%wfn%nel = chk%wfn%nel + 1
         call stop_timing(2)
      end if

      if (set%runtyp == p_run_vea .or. set%runtyp == p_run_vipea &
         & .or. set%runtyp == p_run_vomega) then
         call start_timing(2)
         call vea_header(env%unit)
         mol%chrg = mol%chrg - 1
         chk%wfn%nel = chk%wfn%nel + 1
         if (mod(chk%wfn%nel, 2) /= 0) chk%wfn%nopen = 1
         call calc%singlepoint(env, mol, chk, 1, exist, etot2, g, sigma, egap, res)
         ea = etot - etot2 - ipeashift
         write (env%unit, '(72("-"))')
         write (env%unit, '("empirical EA shift (eV):",f10.4)') &
            &                  autoev * ipeashift
         write (env%unit, '("delta SCC EA (eV):",f10.4)') autoev * ea
         write (env%unit, '(72("-"))')

         mol%chrg = mol%chrg + 1
         chk%wfn%nel = chk%wfn%nel - 1
         call stop_timing(2)
      end if

      ! ------------------------------------------------------------------------
      ! vomega (electrophilicity) index
      if (set%runtyp == p_run_vomega) then
         write (env%unit, '(a)')
         write (env%unit, '(72("-"))')
         write (env%unit, '(a,1x,a)') &
            "Calculation of global electrophilicity index", &
            "(IP+EA)²/(8·(IP-EA))"
         vomega = (ip + ea)**2 / (8 * (ip - ea))
         write (env%unit, '("Global electrophilicity index (eV):",f10.4)') &
            autoev * vomega
         write (env%unit, '(72("-"))')
      end if

      ! ------------------------------------------------------------------------
      ! Fukui Index from Mulliken population analysis
      if (set%runtyp == p_run_vfukui) then
         allocate (fukui(3, mol%n))
         call vfukui(env, mol, chk, calc, fukui)
      end if

      ! ------------------------------------------------------------------------
      ! numerical hessian calculation
      !
      ! `--hess`/`--bhess` always take this path; `--ohess` takes it only in its
      ! native half (`--ohess --method ancopt|lancopt|fire|pbc_lbfgs`).  The
      ! redeployed `--ohess` gets its Hessian from the suite frequency chain,
      ! and running both would write `hessian`/`vibspectrum`/`g98.out` twice
      ! (CLI redesign section 3.11.1, T1).
      if (landscape_native_hessian_path(set%runtyp, landscape_native_engine)) then
         if (set%runtyp == p_run_bhess .and. set%mode_extrun /= p_ext_turbomole) then
            call generic_header(env%unit, "Biased Numerical Hessian", 49, 10)
         else if (set%runtyp == p_run_bhess .and. set%mode_extrun == p_ext_turbomole) then
            call generic_header(env%unit, "Biased Analytical TM Hessian", 49, 10)
         else if (set%mode_extrun == p_ext_turbomole) then
            call generic_header(env%unit, "Analytical TM Hessian", 49, 10)
         else
            call numhess_header(env%unit)
         end if
         if (mol%npbc > 0) then
            call env%error("Phonon calculations under PBC are not implemented", source)
         end if
         call start_timing(5)
         call numhess &
            &       (env, mol, chk, calc, &
            &        egap, set%etemp, set%maxscciter, etot, g, sigma, fres)
         call stop_timing(5)

         call env%checkpoint("Hessian calculation terminated")
      end if

      ! reset the gap, since it is currently not updated in ancopt and numhess
      if (allocated(chk%wfn%emo)) then
         res%hl_gap = chk%wfn%emo(chk%wfn%ihomo + 1) - chk%wfn%emo(chk%wfn%ihomo)
      end if

      ! CPCM-X post-SCF solvation
      if (allocated(calc%solvation)) then
         if (allocated(calc%solvation%cpxsolvent)) then
            select type (calc)
            type is (TxTBCalculator)
               call generic_header(env%unit, "CPCM-X post-SCF solvation evaluation", 49, 10)
               if (set%gfn_method /= 2) call env%warning("CPCM-X was parametrized for GFN2-xTB. &
                  &The results are probably inaccurate with other methods.")
               Call cpx%setup(env, calc%solvation%cpxsolvent)
               Call env%checkpoint("CPCM-X setup terminated")
               cpxcalc = calc
               deallocate (cpxcalc%solvation)
               call cpxcalc%singlepoint(env, mol, chk, 1, .false., energy_gas, g, sigma, egap, res)
          Call cpx%calc_solv(env, calc%solvation%cpxsolvent, energy_gas, 0.4_wp, 298.15_wp, 500, 0.0001_wp, res%e_total)
               Call cpx%print(set%verbose)
               Call env%checkpoint("CPCM-X post-SCF solvation evaluation terminated")
            type is (TGFFCalculator)
               call env%error("CPCM-X is not possible with a force field.", source)
            end select
         end if
      end if

      call env%checkpoint("Calculation terminated")

      ! ========================================================================
      ! PRINTOUT SECTION
      if (allocated(set%property_file)) then
         call open_file(iprop, set%property_file, 'w')
         if (iprop == -1) then
            iprop = env%unit
            deallocate (set%property_file)
         else
            write (env%unit, '(/,a)') "Property printout bound to '"//set%property_file//"'"
            if (allocated(cdum)) deallocate (cdum)
            call get_command(length=l)
            allocate (character(len=l) :: cdum)
            call get_command(cdum)
            write (iprop, '("command:  ''",a,"''")') cdum
            call rdvar('HOSTNAME', cdum, err)
            if (err == 0) &
               write (iprop, '("hostname: ''",a,"''")') cdum
            write (iprop, '("date:     ",a)') prtimestring('S')
         end if
      else
         iprop = env%unit
      end if

      call generic_header(iprop, 'Property Printout', 49, 10)
      if (lgrad) then
         call writeResultsTurbomole(mol, energy=etot, gradient=g, sigma=sigma, iter_number=iter_needed)
         if (allocated(basename)) then
            cdum = basename//'.engrad'
         else
            cdum = 'xtb-orca.engrad'
         end if
         call open_file(ich, cdum, 'w')
         call writeResultsOrca(ich, mol, etot, g)
         call close_file(ich)
      end if
      if (mol%ftype == fileType%gaussian) then
         if (allocated(basename)) then
            cdum = basename//'.EOu'
         else
            cdum = 'xtb-gaussian.EOu'
         end if
         call open_file(ich, cdum, 'w')
         call writeResultsGaussianExternal(ich, etot, res%dipole, g)
         call close_file(ich)
      end if

      if (set%periodic) then
         write (*, *) 'Periodic properties'
      else
         select type (calc)
         type is (TxTBCalculator)
            call main_property(iprop, env, mol, chk%wfn, calc%basis, calc%xtbData, res, &
               & calc%solvation, set%acc)
            call main_cube(set%verbose, mol, chk%wfn, calc%basis, res)
         type is (TTBLiteCalculator)
            call tblite_property(iprop, env, chk, calc, mol, res)
         type is (TGFFCalculator)
            call gfnff_property(iprop, mol%n, mol%xyz, calc%topo, chk%nlist)
         type is (TPTBCalculator)
            call ptb_property(iprop, env, chk, calc, mol, res)
         end select
      end if

      if (printTopo%any()) then
         select type (calc)
         type is (TGFFCalculator)
            call write_json_gfnff_lists(mol%n, res%e_total, res%gnorm, calc%topo, calc%neigh, chk%nlist, printTopo)
         end select
      end if
      if ((set%runtyp == p_run_opt) .or. (set%runtyp == p_run_ohess) .or. &
          (set%runtyp == p_run_omd) .or. (set%runtyp == p_run_screen) .or. &
          (set%runtyp == p_run_metaopt) .or. (set%runtyp == p_run_bhess)) then
         call main_geometry(iprop, mol)
      end if

      ! The frequency printout of the native Hessian path; the redeployed
      ! `--ohess` prints and stores its frequencies through the suite chain, so
      ! a second `main_freq` would repeat the whole thermodynamic table and
      ! overwrite the five files (CLI redesign section 3.11.1, T2).
      if (landscape_native_hessian_path(set%runtyp, landscape_native_engine)) then
         call generic_header(iprop, 'Frequency Printout', 49, 10)
         call main_freq(iprop, mol, chk%wfn, fres)
      end if

      if (set%pr_json) then
         select type (calc)
         type is (TxTBCalculator)
            call open_file(ich, 'xtbout.json', 'w')
            call main_xtb_json(ich, &
                               mol, chk%wfn, calc%basis, res, fres)
            call close_file(ich)
         type is (TPTBCalculator)
            call open_file(ich, 'xtbout.json', 'w')
            call main_ptb_json(ich, &
                               mol, chk%wfn, calc, res, fres)
            call close_file(ich)
         end select
      end if

      ! Unified machine readable document of a native run type (report 117).
      !
      ! The native run types (`--sp`, `--grad`, `--opt`, `--hess`, `--ohess`,
      ! `--md`, ...) have no `landscape_result`; `landscape_json_write_run`
      ! fills the same core block the suite tasks write, under the same
      ! `landscape_<task>.json` name, so one parser reads every task.  The
      ! document is written independently of `--json` (which still writes
      ! `xtbout.json` unchanged) and of the run type's own output files.  `res`
      ! holds the total energy and gradient norm of the run type specific
      ! block above, so the document describes the geometry `mol` ended on.
      if (lset%ljson .and. &
         & .not. landscape_is_suite_run(set%runtyp, landscape_native_engine)) then
         block
            logical :: ljson_ok
            call landscape_json_write_run(env, mol, res%e_total, res%gnorm, &
               & res%converged, landscape_runtyp_task_name(set%runtyp), fname, &
               & ljson_ok)
         end block
      end if

      if (allocated(set%property_file)) then
         if (iprop /= -1 .and. iprop /= env%unit) then
            call write_energy(iprop, res, fres, &
               & (set%runtyp == p_run_hess) .or. (set%runtyp == p_run_ohess) .or. (set%runtyp == p_run_bhess))
            call close_file(iprop)
         end if
      end if

      ! The second geometry writer of the native optimisers.  The run types
      ! that go to the landscape suite write `xtbopt.xyz` themselves
      ! (`xtb_landscape_log::landscape_write_files`), so for them this block
      ! would write the same file twice and print a second, contradictory
      ! `optimized geometry written to:` line (CLI rewrite read point 13).
      if (((set%runtyp == p_run_opt) .or. (set%runtyp == p_run_ohess) .or. &
          (set%runtyp == p_run_omd) .or. (set%runtyp == p_run_screen) .or. &
          (set%runtyp == p_run_metaopt) .or. (set%runtyp == p_run_bhess)) &
          .and. .not. landscape_is_suite_run(set%runtyp, landscape_native_engine)) then
         call generateFileName(tmpname, 'xtbopt', extension, mol%ftype)
         if (mol%ftype == 0) mol%ftype = getFileType(tmpname)
         write (env%unit, '(/,a,1x,a,/)') &
            "optimized geometry written to:", tmpname
         call open_file(ich, tmpname, 'w')
         call writeMolecule(mol, ich, energy=res%e_total, gnorm=res%gnorm)
         call close_file(ich)
      end if

      select type (calc)
      type is (TxTBCalculator)
         call write_energy(env%unit, res, fres, anyhess)
      type is (TOniomCalculator)
         call write_energy_oniom(env%unit, res, fres, anyhess)
      type is (TPTBCalculator)
      class default
         call write_energy_gff(env%unit, res, fres, anyhess)
      end select

      ! ------------------------------------------------------------------------
      !  xtb molecular dynamics
      if ((set%runtyp == p_run_md) .or. (set%runtyp == p_run_omd)) then
         if (metaset%maxsave > 0) then
            if (mol%npbc > 0) then
               call env%error("Metadynamic under PBC is not implemented", source)
            end if
            call metadyn_header(env%unit)
         else
            call md_header(env%unit)
         end if
         fixset%n = 0 ! no fixing for MD runs
         call start_timing(6)
         idum = 0
         select type (calc)
         class default
            if (set%shake_md) call init_shake(mol%n, mol%at, mol%xyz, chk%wfn%wbo)
         type is (TGFFCalculator)
            if (set%shake_md) call gff_init_shake(mol%n, mol%at, mol%xyz, calc%topo)
         end select
         call md &
            &     (env, mol, chk, calc, &
            &      egap, set%etemp, set%maxscciter, etot, g, sigma, 0, set%temp_md, idum)
         call env%checkpoint("Molecular dynamics terminated")
         call stop_timing(6)
      end if

      ! ------------------------------------------------------------------------
      !  metadynamics
      if (set%runtyp == p_run_metaopt) then
         if (mol%npbc > 0) then
            call env%warning("Metadynamic under PBC is not implemented", source)
         end if
         call metadyn_header(env%unit)
         ! check if ANCOPT already convered
         if (murks) then
            call env%error('Optimization did not converge, aborting', source)
         end if
         write (env%unit, '(1x,"output written to xtbmeta.log")')
         call open_file(ich, 'xtbmeta.log', 'w')
         call writeMolecule(mol, ich, fileType%xyz, energy=etot, gnorm=norm2(g))
         k = metaset%nstruc + 1
         call start_timing(6)
         do l = k, metaset%maxsave
            metaset%nstruc = l
            metaset%xyz(:, :, metaset%nstruc) = mol%xyz
            ! randomize structure to avoid zero RMSD
            do i = 1, mol%n
               do j = 1, 3
                  call random_number(er)
                  mol%xyz(j, i) = mol%xyz(j, i) + 1.0e-6_wp * er
               end do
            end do
            call geometry_optimization &
               &     (env, mol, chk, calc, &
               &      egap, set%etemp, set%maxscciter, set%optset%maxoptcycle, etot, g, sigma, &
               &      set%optset%optlev, set%verbose, .true., murks)
            if (.not. set%verbose) then
               write (env%unit, '("current energy:",1x,f20.8)') etot
            end if
            if (murks) then
               call close_file(ich)
               write (env%unit, '(/,3x,"***",1x,a,1x,"***",/)') &
                  "FAILED TO CONVERGE GEOMETRY OPTIMIZATION"
               call touch_file('NOT_CONVERGED')
            end if
            call writeMolecule(mol, ich, fileType%xyz, energy=etot, gnorm=norm2(g))
         end do
         call close_file(ich)
         call stop_timing(6)
      end if

      ! ------------------------------------------------------------------------
      !  path finder
      if (set%runtyp == p_run_path) then
         call rmsdpath_header(env%unit)
         if (mol%npbc > 0) then
            call env%warning("Metadynamics under PBC are not implemented", source)
         end if
         call start_timing(4)
         call bias_path(env, mol, chk, calc, egap, set%etemp, set%maxscciter, etot, g, sigma)
         call stop_timing(4)
      end if

      ! ------------------------------------------------------------------------
      !  screen over input structures
      if (set%runtyp == p_run_screen) then
         call start_timing(8)
         call screen(env, mol, chk, calc, egap, set%etemp, set%maxscciter, etot, g, sigma)
         call stop_timing(8)
      end if

      ! ------------------------------------------------------------------------
      !  mode following for conformer search
      if (set%runtyp == p_run_modef) then
         if (mol%npbc > 0) then
            call env%warning("Modefollowing under PBC is not implemented", source)
         end if
         call start_timing(9)
         call modefollow(env, mol, chk, calc, egap, set%etemp, set%maxscciter, etot, g, sigma)
         call stop_timing(9)
      end if

      ! ------------------------------------------------------------------------
      !  optimize along MD from xtb.trj for conformer searches
      if (set%runtyp == p_run_mdopt) then
         call start_timing(10)
         call mdopt(env, mol, chk, calc, egap, set%etemp, set%maxscciter, etot, g, sigma)
         call stop_timing(10)
      end if

      ! ------------------------------------------------------------------------
      !  to further speed up xtb calculations we dump our most important
      !  quantities in a restart file, so we can save some precious seconds
      select type (calc)
      type is (TxTBCalculator)
         if (restart) then
            call writeRestart(env, chk%wfn, 'xtbrestart', set%gfn_method)
         end if
      type is (TTBLiteCalculator)
         if (restart) call dumpRestart(env, chk, 'xtbrestart')
      end select

      ! ------------------------------------------------------------------------
      !  we may have generated some non-fatal errors, which have been saved,
      !  so we should tell the user, (s)he may want to know what went wrong
      call env%show("Runtime exception occurred")
      call raise('F', 'Some non-fatal runtime exceptions were caught,'// &
         &           ' please check:')

      ! ------------------------------------------------------------------------
      !  print all files xtb interacted with while running (for debugging mainly)
      if (set%verbose) then
         write (env%unit, '(a)')
         write (env%unit, '(72("-"))')
         call print_filelist(env%unit)
      end if

      ! ------------------------------------------------------------------------
      !  make some post processing afterward, show some timings and stuff
      call finalize_xtb(env, landscape_native_engine) 

   end subroutine xtbMain

!> Parse command line arguments and forward them to settings
   subroutine parseArguments(env, args, inputFile, paramFile, lgrad, &
         & restart, gsolvstate, strict, copycontrol, coffee, printTopo, oniom, dipro, tblite, lset, &
         & lsset, lvset, ltset, lircset, landscape_native_engine, acc_was_set)

      use xtb_mctc_global, only: persistentEnv

      !> Name of error producer
      character(len=*), parameter :: source = "prog_main_parseArguments"

      !> Calculation environment
      type(TEnvironment) :: env

      !> Command line argument parser
      type(TArgParser) :: args

      !> Detailed input file name
      character(len=:), allocatable, intent(out) :: inputFile

      !> Parameter file name
      character(len=:), allocatable, intent(out) :: paramFile

      !> Reference state for solvation free energies
      integer, intent(out) :: gsolvstate

      !> Restart calculation
      logical, intent(out) :: restart

      !> Handle warnings as errors
      logical, intent(out) :: strict

      !> Debugging with a lot of caffeine
      logical, intent(out) :: coffee

      !> topology printout list
      type(TPrintTopo), intent(out) :: printTopo

      !> Print the gradient to file
      logical, intent(out) :: lgrad

      !> Copy the detailed input file
      logical, intent(out) :: copycontrol

      !> Input for ONIOM model
      type(oniom_input), intent(out) :: oniom

      !> Input for DIPRO
      type(jab_input), intent(inout) :: dipro

      !> Stuff for second argument parser
!   integer  :: narg
!   character(len=p_str_length), dimension(p_arg_length) :: argv
!   type(TAtomList) :: atl
!   integer, allocatable :: list(:)

      !> Input for TBLite calculator
      type(TTBLiteInput), intent(out) :: tblite

      !> Landscape task settings (`--lopt`, `--lopt-engine`, `--lopt-level`)
      type(landscape_opt_settings), intent(inout) :: lset

      !> Landscape scan settings (`--lscan`, `--lscan-coord`, ...)
      type(landscape_scan_settings), intent(inout) :: lsset

      !> Landscape frequency settings (`--lvib`, `--lvib-hessian`, ...)
      type(landscape_vib_settings), intent(inout) :: lvset

      !> Landscape transition state settings (`--lts`, `--lts-method`, ...)
      type(landscape_ts_settings), intent(inout) :: ltset

      !> Landscape intrinsic reaction coordinate settings (`--lirc`, ...)
      type(landscape_irc_settings), intent(inout) :: lircset

      !> `--acc` was named on the command line.  It is the state the six
      !> synchronisation points of `xtbMain` read before they copy the global
      !> `set%acc` into a landscape engine control, so that an explicit
      !> `--acc` wins over the SCC accuracy a convergence level derives (CLI
      !> redesign section 3.5, R5).  Declared as a dummy - and not as a
      !> variable of the program unit - because the project hard constraint
      !> "no module-level mutable variables" forbids mutable program-unit
      !> state; reset at the top of the flag loop, so a second run in the same
      !> process (C API) starts clean.
      logical, intent(out) :: acc_was_set

!$    integer :: omp_get_num_procs
      integer :: nFlags
      integer :: idum, ndum
      logical :: ldum
      real(wp) :: ddum
      character(len=:), allocatable :: flag, sec
      logical :: exist
      integer :: iscan, idum2
      real(wp) :: scan_step_user(2)
      logical :: scan_step_given(2)

      !> Method name of the `--method` flag; only the string is stored while
      !> parsing, it is resolved after the loop (CLI redesign section 4.7)
      character(len=:), allocatable :: landscape_method_name

      !> Raw level keyword of `--opt`/`--ohess`, e.g. `tight`; empty when the
      !> flag was given without a level.  The keyword is kept because the
      !> integer `set%optset%optlev` cannot tell an unknown keyword from
      !> `normal` and cannot carry the suite-only `superloose` (CLI redesign
      !> section 3.5, finding F3 of the stage 3a fix round)
      character(len=:), allocatable :: landscape_opt_level_name

      !> Raw level keyword of the native `--bhess`/`--metaopt`, e.g. `tight`.
      !> The run type flags are mutually exclusive, so one pair of variables
      !> holds the keyword and the flag that supplied it; they are reported in
      !> the parser tail, where the effect of `--strict` is known, because
      !> `optlevel2int` maps an unknown keyword to `normal` silently
      character(len=:), allocatable :: landscape_native_level_name
      character(len=:), allocatable :: landscape_native_level_flag

      !> Effective suite convergence level the `--opt`/`--ohess` keyword maps
      !> to, printed on the banner of the redeployed run
      character(len=:), allocatable :: landscape_opt_conv_level

      !> `--cycles N` was given with a value the native parser accepted
      logical :: cycles_given

      !> Step budget value of `--cycles`
      integer :: cycles_value

      !> The task-specific step budget flags were given (CLI redesign 3.10b)
      logical :: lopt_max_steps_given, lscan_max_steps_given
      logical :: lts_max_steps_given, lirc_max_steps_given

      !> First run type requested on the command line and the flag that
      !> requested it; both are filled by `select_runtyp`, which reports a
      !> second, different run type as an error (CLI redesign section 17.1).
      character(len=:), allocatable :: selected_runtyp, selected_runtyp_flag

      !> Native ANCopt engine selected by `--method` (`p_engine_*`), 0 for a
      !> landscape method name and when the flag was not given (CLI redesign
      !> sections 4.7/4.8).  It is written before the flag loop because an
      !> `intent(out)` dummy is undefined on entry, and every run type - also
      !> the ones without `--method` - reads it afterwards.
      integer, intent(out) :: landscape_native_engine

      !> Unknown flags, blank separated; they are reported after the parsing
      !> loop, where the effect of `--strict` is known (section 17.5 5a).
      !> `message` is the report text of that block and of the `-P`/`--parallel`
      !> value report below
      character(len=:), allocatable :: unknown_flags, message
      integer :: unknown_count

      !> Raw value of a `-P`/`--parallel` argument that is not a positive
      !> integer (`0`, `-2`, `abc`); allocated and empty when the flag was
      !> given without a value, unallocated when the value is valid or the flag
      !> is absent.  The report happens in the parser tail, where the effect of
      !> `--strict` is known - the same lax/strict shape as the unknown option
      !> and the unknown native level there
      character(len=:), allocatable :: landscape_parallel_value

      !> Team size the `-P`/`--parallel` flag asked for, 0 when the flag is
      !> absent or its value was rejected; a value above the processor count is
      !> reported in the parser tail, where `nproc` is known
      integer :: landscape_parallel_team

      !> Raw value of an `-a`/`--acc` argument that is not a finite number
      !> (`nan`, `inf`); unallocated when the value is usable.  It is reported
      !> in the parser tail, where the effect of `--strict` is known - the same
      !> lax/strict shape as the `-P` value below - and the value is **not**
      !> applied, because an SCC accuracy of NaN turns every energy of the run
      !> into NaN instead of failing.
      character(len=:), allocatable :: landscape_acc_value

      !> `--lopt-recalc N` was given, and the value it carried.  The report of a
      !> negative value happens in the parser tail, where the effect of
      !> `--strict` is known - the same lax/strict shape as the `-P` value
      !> above - so the value has to survive the flag loop.
      logical :: lopt_recalc_given
      integer :: lopt_recalc_value

      !> Report text of the `--lopt-recalc` value diagnostic (the internal
      !> write needs a fixed length buffer, like `parallel_note`)
      character(len=128) :: recalc_note

      !> Report text of the `-P` over-subscription warning (the internal write
      !> needs a fixed length buffer; `message` is allocatable and cannot be
      !> written to before it is allocated)
      character(len=128) :: parallel_note

      !> Processor count of this machine, filled by the `-P`/`--parallel` case;
      !> 0 in a build without OpenMP, which therefore never reports an
      !> over-subscription
      integer :: nproc

      !> Exit status of the `--help` page: 0 for a page, 1 for an unknown
      !> topic (the topic lookup is in `src/xhelp.f90`)
      integer :: help_status

      interface
         !> The `--help` page and the `--help-values` table.  The definition is
         !> in `src/xhelp.f90`; the explicit interface is needed for the
         !> optional `values_only` argument of the last call site.
         subroutine help(iunit, topic, status, values_only)
            integer, intent(in) :: iunit
            character(len=*), intent(in) :: topic
            integer, intent(out) :: status
            logical, intent(in), optional :: values_only
         end subroutine help
      end interface

      set%gfn_method = 2
      dipro%diprocalc = .false.
      coffee = .false.
      strict = .false.
      unknown_count = 0
      restart = .true.
      copycontrol = .false.
      lgrad = .false.
      gsolvstate = solutionState%gsolv
      tblite%color = get_xtb_feature('color')
      scan_step_user = 0.0_wp
      scan_step_given = .false.

      ! `landscape_native_engine` is an `intent(out)` dummy and therefore
      ! undefined on entry; it must be defined on every return path, because
      ! the run types without `--method` read it as well (CLI redesign section
      ! 4.7, step 0).  The flag loop is the first statement sequence that can
      ! leave the subroutine through `env%error`, so the value is written here.
      landscape_native_engine = 0

      ! Local parse state of the `--opt <LEVEL>` and `--cycles N` bridges of
      ! the tail (CLI redesign sections 3.5 and 3.10b); a local of a subroutine
      ! without `save` is undefined on entry, so every flag is initialised
      ! here, before the loop can reach a `case`.  `acc_was_set` is the
      ! `intent(out)` dummy of the R5 rule; an `intent(out)` dummy is undefined
      ! on entry just like a local, and the flag loop has `return` paths (the
      ! early-type requests), so it is initialised here rather than only in the
      ! `--acc` case - otherwise those paths would hand an undefined value back
      ! to `xtbMain`.
      cycles_given = .false.
      cycles_value = 0
      lopt_max_steps_given = .false.
      lscan_max_steps_given = .false.
      lts_max_steps_given = .false.
      lirc_max_steps_given = .false.
      acc_was_set = .false.
      lopt_recalc_given = .false.
      lopt_recalc_value = 0
      landscape_parallel_team = 0
      nproc = 0

      nFlags = args%countFlags()
      call args%nextFlag(flag)
      do while (allocated(flag))
         if (len(flag) > 2 .and. flag(1:1) == '-' .and. flag(1:2) /= '--') then
            call env%warning("the use of '"//flag//"' is discouraged, "// &
               & "please use '-"//flag//"' next time", source)
            flag = '-'//flag
         end if
         select case (flag)
         case default
            ! Reported after the loop: `--strict` may stand after this flag,
            ! so whether an unknown option is a warning or an error cannot be
            ! decided here (CLI redesign section 17.5 5a).
            unknown_count = unknown_count + 1
            if (allocated(unknown_flags)) then
               unknown_flags = trim(unknown_flags)//' '//trim(flag)
            else
               unknown_flags = trim(flag)
            end if

         case ('-h', '--help')
            ! An optional topic may follow: an option spelling (with or
            ! without its dashes), a category key or any other word.  The
            ! parser returns a bare word through `nextArg` and the
            ! immediately following flag through `nextFlag`; a missing topic
            ! prints the grouped index and an unknown one is an error.
            call args%nextArg(sec)
            if (.not. allocated(sec)) call args%nextFlag(sec)
            if (allocated(sec)) then
               call help(env%unit, sec, help_status)
            else
               call help(env%unit, '', help_status)
            end if
            call terminate(help_status)

         case ('--help-values')
            ! Machine-readable dump of the legal value sets of the option
            ! table; the contract a completion generator codes against.
            call help(env%unit, '', help_status, values_only=.true.)
            call terminate(help_status)

         case ('--citation')
            call citation(env%unit)
            call terminate(0)

         case ('--license')
            call disclaimer(env%unit)
            call terminate(0)

         case ('--version')
            call xtb_header(env%unit)
            call terminate(0)

         case ('-v', '--verbose')
            set%verbose = .true.

         case ('-V', '--very-verbose')
            set%verbose = .true.
            set%veryverbose = .true.

         case ('--define')
            call set_define

         case ('-P', '--parallel')
!$          if (.false.) then
               call env%warning('Program compiled without threading support', source)
!$          end if
            ! Always remove next argument to keep argument parsing consistent
            call args%nextArg(sec)
            if (.not. allocated(sec)) then
               ! `-P` as the last argument used to be consumed and ignored: the
               ! missing value is reported in the tail instead of guessed
               landscape_parallel_value = ''
            else if (getValue(env, sec, idum)) then
               ! A parsed integer still reaches the OpenMP runtime, including a
               ! non-positive one, which the runtime maps to a single thread
               ! (measured, report 189 section 5).  Keeping the call keeps the
               ! accepted path byte identical and the value applied as given;
               ! only the report of a value that is not a positive integer is
               ! new, and it is a warning in a lax run and fatal with
               ! `--strict` (parser tail)
!$             nproc = omp_get_num_procs()
!$             call omp_set_num_threads(idum)
#ifdef WITH_MKL
!$             call mkl_set_num_threads(idum)
#endif
               if (idum < 1) then
                  landscape_parallel_value = sec
               else
                  landscape_parallel_team = idum
               end if
            else
               ! `getValue` reported the unparsable field itself; the flag and
               ! its accepted range are named in the tail
               landscape_parallel_value = sec
            end if

         case ('--blas-threads')
            !> Machine setting, not a run type selector: it deliberately does
            !> not go through `select_runtyp`, so `--sp --blas-threads 1` stays
            !> one single point instead of a reported run type conflict.
            call args%nextArg(sec)
            if (allocated(sec)) then
               if (getValue(env, sec, idum)) then
                  if (idum < 1) then
                     call env%error("--blas-threads needs at least one thread", source)
                  else
                     call set_blas_threads(env, idum)
                  end if
               else
                  ! The value was given but is not an integer (`getValue` reports
                  ! the unparsable field itself), so the message names the flag
                  ! instead of claiming that nothing was provided.
                  call env%error("--blas-threads needs an integer >= 1", source)
               end if
            else
               call env%error("--blas-threads needs an integer >= 1", source)
            end if

         case ('--restart')
            restart = .true.

         case ('--norestart')
            restart = .false.

         case ('--copy')
            copycontrol = .true.

         case ('--nocopy')
            copycontrol = .false.

         case ('--strict')
            strict = .true.

         case ('-I', '--input')
            call args%nextArg(inputFile)
            if (.not. allocated(inputFile)) then
               call env%error("Filename for detailed input is missing", source)
            end if

         case ('--namespace')
            call args%nextArg(persistentEnv%io%namespace)
            if (.not. allocated(persistentEnv%io%namespace)) then
               call env%error("Namespace argument is missing", source)
            end if

         case ('--vparam')
            call args%nextArg(paramFile)
            if (.not. allocated(paramFile)) then
               call env%error("Filename for --vparam is missing", source)
            else
               tblite%param = paramFile
            end if

         case ('--coffee')
            coffee = .true.

         case ('-a', '--acc')
            call args%nextArg(sec)
            if (allocated(sec)) then
               if (getValue(env, sec, ddum)) then
                  if (.not. ieee_is_finite(ddum)) then
                     ! `read(sec,*,iostat=)` parses `nan`/`inf`; a NaN
                     ! accuracy used to be stored and every SCC energy of the
                     ! run silently became NaN.  Reported in the parser tail
                     ! (below) and not applied, so the default of the level or
                     ! of `set%acc` stays in force.
                     landscape_acc_value = sec
                  else if (ddum < 1.e-4_wp) then
                     call env%warning("We cannot provide this level of accuracy, "//&
                        & "resetted accuracy to 0.0001", source)
                     set%acc = 1.e-4_wp
                  else if (ddum > 1.e+3_wp) then
                     call env%warning("We cannot provide this level of accuracy, "//&
                        & "resetted accuracy to 1000", source)
                     set%acc = 1.e+3_wp
                  else
                     set%acc = ddum
                  end if
                  ! A rejected value does not count as "named on the command
                  ! line", so the accuracy derived from the convergence level
                  ! stays in force; an accepted one wins over it (CLI redesign
                  ! section 3.5, R5).
                  if (.not. allocated(landscape_acc_value)) acc_was_set = .true.
               end if
               tblite%accuracy = set%acc
            else
               call env%error("Accuracy is not provided", source)
            end if

         case ('-c', '--chrg', '--charge')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_chrg(env, sec)
            else
               call env%error("Molecular charge is not provided", source)
            end if
         
         case('--ceh')
            if (get_xtb_feature('tblite')) then
               allocate(tblite%ceh)
               call select_runtyp(env, 'prescc', '--ceh', selected_runtyp, selected_runtyp_flag)
               ! check if numerical gradients should be calculated !
               call args%nextArg(sec)
               if (allocated(sec)) then
                  if (sec == 'grad') then
                     tblite%ceh%grad = .true.
                     ! check if step size is provided !
                     call args%nextArg(sec)
                     if (allocated(sec)) then
                        if (getValue(env, sec, ddum)) then
                           tblite%ceh%step = ddum
                        end if
                     end if

                  else
                     call env%warning("Unknown CEH option '"//sec//"' provided", source)
                  end if
               end if 
            else
               call env%error("CEH charges are only available through tblite library", source)
               return
            end if


         case ('-u', '--uhf')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_spin(env, sec)
            else
               call env%error("Number of unpaired electrons is not provided", source)
            end if

         case ("--efield")
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_efield(env, sec)
               ! Set electric field for tblite
               if (.not. allocated(tblite%efield)) then
                  allocate(tblite%efield(3))
               end if
               tblite%efield = set%efield
            else
               call env%error("Electric field is not provided", source)
            end if

         case ('--gfn')
            call args%nextArg(sec)
            if (allocated(sec)) then
               ! Only the four documented methods are accepted.  An unknown
               ! value used to fall through to the GFN2 default without an
               ! error, so the run silently computed a different method than
               ! the caller asked for.  The decoder is shared with the
               ! xcontrol `$gfn method` key, so the two routes cannot drift.
               select case (gfn_method_from_literal(sec))
               case (0)
                  call set_gfn(env, 'method', sec)
                  call set_exttyp('eht')
                  tblite%method = "gfn0"
               case (1)
                  call set_gfn(env, 'method', sec)
                  tblite%method = "gfn1"
               case (2)
                  call set_gfn(env, 'method', sec)
                  tblite%method = "gfn2"
               case (p_gfn_method_ff)
                  call set_gfn(env, 'method', sec)
                  call set_exttyp('ff')
                  tblite%method = "gfnff"
               case default
                  call env%error(gfn_method_error_message(sec), source)
               end select
            else
               call env%error("No method provided for --gfn option", source)
            end if

         case ('--gfn1')
            call set_gfn(env, 'method', '1')
            call env%warning("The use of '"//flag//"' is discouraged, "//&
               & "please use '--gfn 1' next time", source)
            tblite%method = "gfn1"

         case ('--gfn2')
            call set_gfn(env, 'method', '2')
            call set_gfn(env, 'd4', 'true')
            tblite%method = "gfn2"

         case ('--gfn0')
            call set_gfn(env, 'method', '0')
            call set_exttyp('eht')
            call env%warning("The use of '"//flag//"' is discouraged, "//&
               & "please use '--gfn 0' next time", source)

         case ('--gfnff')
            call set_exttyp('ff')

         case ('--gff')
            call set_exttyp('ff')

         case ('--mcgfnff')
            call set_exttyp('mcff')

         case ('--iff')
            call set_exttyp('iff')

         case ('--ptb')
            call set_exttyp('ptb')
            if (.not. get_xtb_feature('tblite')) then
               call ptb_feature_not_implemented(env)
               return
            endif

         case ('--tblite')
            if (get_xtb_feature('tblite')) then
               call set_exttyp('tblite')
            else
               call env%error("Compiled without support for tblite library", source)
               return
            end if

         case ('--color')
            call args%nextArg(sec)
            if (allocated(sec)) then
               select case (sec)
               case ('auto')
                  tblite%color = get_xtb_feature('color')
               case ('always')
                  tblite%color = .true.
               case ('never')
                  tblite%color = .false.
               case default
                  call env%warning("Unknown color option '"//sec//"' provided", source)
               end select
            else
               call env%error("No color scheme provided for --color option", source)
            end if

         case ('--spinpol')
            if (get_xtb_feature('tblite')) then
               tblite%spin_polarized = .true.
            else
           call env%error("Compiled without support for tblite library. This is required for spin-polarization", source)
               return
            end if

         case ('--dipro')
            if (get_xtb_feature('tblite')) then
               dipro%diprocalc = .true.
               call select_runtyp(env, 'scc', '--dipro', selected_runtyp, selected_runtyp_flag)
               call args%nextArg(sec)
               if (allocated(sec)) then
                  read (sec, '(f10.3)') dipro%othr
               else
                  dipro%othr = 0.1_wp
               end if
            else
               call env%error("Compiled without support for tblite library. This is required for DIPRO", source)
               return
            end if

         case ('--oniom')
            call set_exttyp('oniom')
            call args%nextArg(sec)

            if (.not. allocated(sec)) then ! handle no argument case !
               call env%error("No inner region is  provided for ONIOM", source)
               return
            end if
            call move_alloc(sec, oniom%first_arg)

            call args%nextArg(sec)
            if (.not. allocated(sec)) then
               call env%warning("No method is specified for ONIOM," &
                     &//achar(10)//" default gfn2:gfnff combination will be used", source)
               call move_alloc(oniom%first_arg, sec)
            end if

            inquire (file=sec, exist=exist)
            if (exist) then
               sec = read_whole_file(sec)
            end if
            call move_alloc(sec, oniom%second_arg)

         case ('--cut')
            call set_cut

         case ('--etemp')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_scc(env, 'temp', sec)
               !set etemp for tblite
               ldum = getValue(env, sec, tblite%etemp)
            else
               call env%error("Temperature in --etemp option is missing", source)
            end if

         case ('--esp')
            call select_runtyp(env, 'scc', '--esp', selected_runtyp, selected_runtyp_flag)
            call set_write(env, 'esp', 'true')

         case ('--stm')
            call select_runtyp(env, 'scc', '--stm', selected_runtyp, selected_runtyp_flag)
            call set_write(env, 'stm', 'true')

         case ('--cma')
            call set_cma

         case ('--tm')
            call set_exttyp('turbomole')

         case ('--enso')
            call set_enso_mode

         case ('--json')
            call set_write(env, 'json', 'true')
            Call setWRtopo("json", printTopo)

         case ('--ceasefiles')
            restart = .false.
            set%verbose = .false.
            set%ceasefiles = .true.
            call set_write(env, 'wiberg', 'false')
            call set_write(env, 'charges', 'false')
#ifdef _WIN32
            call set_opt(env, 'logfile', 'NUL')
#else
            call set_opt(env, 'logfile', '/dev/null')
#endif

         case ('--orca')
            call set_exttyp('orca')

         case ('--driver')
            call set_exttyp('driver')
            call args%nextArg(sec)
            if (allocated(sec)) then
               set%ext_driver%executable = sec
            end if

         case ('--mopac')
            call set_exttyp('mopac')

         case ('--pop')
            call set_write(env, 'mulliken', 'true')

         case ('--molden')
            call set_write(env, 'mos', 'true')

         case ('--dipole')
            call set_write(env, 'dipole', 'true')

         case ('--wbo')
            call set_write(env, 'wiberg', 'true')

         case ('--lmo')
            call set_write(env, 'mulliken', 'true')
            call set_write(env, 'lmo', 'true')

         case ('--ewin')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_siman(env, 'ewin', sec)
            else
               call env%error("Real argument for --ewin is missing", source)
            end if

         case ('--fod')
            call set_write(env, 'fod', 'true')
            call set_scc(env, 'temp', '5000.0')

         case ('--iterations', '--maxiterations')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_scc(env, 'maxiterations', sec)
               tblite%max_iter = set%maxscciter
            else
               call env%error("Integer argument for --iterations is missing", source)
            end if

         case ('--cycles')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_opt(env, 'maxcycle', sec)
               ! Keep the value for the suite bridge of the tail; `set_opt`
               ! rejects a negative or non-numeric value, so only an accepted
               ! value marks the flag as given (CLI redesign section 3.10b)
               read(sec, *, iostat=idum) idum2
               if (idum == 0 .and. idum2 >= 0) then
                  cycles_given = .true.
                  cycles_value = idum2
               end if
            else
               call env%error("Integer argument for --cycles is missing", source)
            end if

         case ('-g', '--gbsa')
            call args%nextArg(sec)
            if (allocated(sec)) then
               ! Read GBSA solvent name
               call set_gbsa(env, 'solvent', sec)
               call set_gbsa(env, 'alpb', 'false')
               call set_gbsa(env, 'kernel', 'still')
               ! Add solvation model also to tblite input
               if (.not. allocated(tblite%solvation)) then
                  allocate(tblite%solvation)
               end if
               if (allocated(tblite%solvation%solvation_model)) then
                  call env%error("Cannot specify multiple solvation models", source)
               end if
               tblite%solvation%solvation_model = "gbsa"
               tblite%solvation%solvent = sec
               ! Read possible reference state
               call args%nextArg(sec)
               if (allocated(sec)) then
                  if (sec == 'gsolv') then
                     gsolvstate = solutionState%gsolv
                     tblite%solvation%reference_state = "gsolv"
                  else if (sec == 'reference') then
                     gsolvstate = solutionState%reference
                     tblite%solvation%reference_state = "reference"
                  else if (sec == 'bar1M') then
                     gsolvstate = solutionState%mol1bar
                     tblite%solvation%reference_state = "bar1M"
                  else
                     call env%warning("Unknown reference state '"//sec//"'", source)
                  end if
               end if
            else
               call env%error("No solvent name provided for GBSA", source)
            end if

         case ('--alpb')
            call args%nextArg(sec)
            call set_gbsa(env, 'alpb', 'true')
            if (allocated(sec)) then
               ! Read ALPB solvent name
               call set_gbsa(env, 'solvent', sec)
               ! Add solvation model also to tblite input
               if (.not. allocated(tblite%solvation)) then
                  allocate(tblite%solvation)
               end if
               if (allocated(tblite%solvation%solvation_model)) then
                  call env%error("Cannot specify multiple solvation models", source)
               end if
               tblite%solvation%solvation_model = "alpb"
               tblite%solvation%solvent = sec
               ! Read possible reference state
               call args%nextArg(sec)
               if (allocated(sec)) then
                  if (sec == 'gsolv') then
                     gsolvstate = solutionState%gsolv
                     tblite%solvation%reference_state = "gsolv"
                  else if (sec == 'reference') then
                     gsolvstate = solutionState%reference
                     tblite%solvation%reference_state = "reference"
                  else if (sec == 'bar1M') then
                     gsolvstate = solutionState%mol1bar
                     tblite%solvation%reference_state = "bar1M"
                  else
                     call env%warning("Unknown reference state '"//sec//"'", source)
                  end if
               end if
            else
               call env%error("No solvent name provided for ALPB", source)
            end if

         case ('--gbe')
            call args%nextArg(sec)
            if (allocated(sec)) then
               if (.not. allocated(tblite%solvation)) then
                  allocate(tblite%solvation)
               end if
               if (allocated(tblite%solvation%solvation_model)) then
                  call env%error("Cannot specify multiple solvation models", source)
               end if
               tblite%solvation%solvation_model = "gbe"
               tblite%solvation%solvent = sec
            else
               call env%error("No solvent name or dielectric constant provided for GBE.", source)
            end if

         case ('--gb')
            call args%nextArg(sec)
            if (allocated(sec)) then
               if (.not. allocated(tblite%solvation)) then
                  allocate(tblite%solvation)
               end if
               if (allocated(tblite%solvation%solvation_model)) then
                  call env%error("Cannot specify multiple solvation models", source)
               end if
               tblite%solvation%solvation_model = "gb"
               tblite%solvation%solvent = sec
            else
               call env%error("No solvent name or dielectric constant provided for GB.", source)
            end if

         case ('--cosmo', '--tmcosmo')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_gbsa(env, 'solvent', sec)
               call set_gbsa(env, flag(3:), 'true')
               ! Add solvation model also to tblite input
               if (.not. allocated(tblite%solvation)) then
                  allocate(tblite%solvation)
               end if
               if (allocated(tblite%solvation%solvation_model)) then
                  call env%error("Cannot specify multiple solvation models", source)
               end if
               tblite%solvation%solvation_model = "cosmo"
               tblite%solvation%solvent = sec
               ! Read possible reference state
               call args%nextArg(sec)
               if (allocated(sec)) then
                  if (sec == 'gsolv') then
                     gsolvstate = solutionState%gsolv
                     tblite%solvation%reference_state = "gsolv"
                  else if (sec == 'reference') then
                     gsolvstate = solutionState%reference
                     tblite%solvation%reference_state = "reference"
                  else if (sec == 'bar1M') then
                     gsolvstate = solutionState%mol1bar
                     tblite%solvation%reference_state = "bar1M"
                  else
                     call env%warning("Unknown reference state '"//sec//"'", source)
                  end if
               end if
            else
               call env%error("No solvent name provided for COSMO", source)
            end if

         case ('--cpcmx')
            if (get_xtb_feature('cpcmx')) then
               call args%nextArg(sec)
               if (allocated(sec)) then
                  call set_gbsa(env, 'solvent', 'infinity')
                  call set_gbsa(env, 'cosmo', 'true')
                  call set_gbsa(env, 'cpcmx', sec)
               else
                  call env%error("No solvent name provided for CPCM-X", source)
               end if
            else
               call env%error("The CPCM-X library was not included in this version of xTB.", source)
            end if

         case ('--scc', '--sp')
            call select_runtyp(env, 'scc', '--scc', selected_runtyp, selected_runtyp_flag)

         case ('--vip')
            call set_gfn(env, 'method', '1')
            call select_runtyp(env, 'vip', '--vip', selected_runtyp, selected_runtyp_flag)

         case ('--vea')
            call set_gfn(env, 'method', '1')
            call select_runtyp(env, 'vea', '--vea', selected_runtyp, selected_runtyp_flag)

         case ('--vipea')
            call set_gfn(env, 'method', '1')
            call select_runtyp(env, 'vipea', '--vipea', selected_runtyp, selected_runtyp_flag)

         case ('--vomega')
            call set_gfn(env, 'method', '1')
            call select_runtyp(env, 'vomega', '--vomega', selected_runtyp, selected_runtyp_flag)

         case ('--vfukui')
            call select_runtyp(env, 'vfukui', '--vfukui', selected_runtyp, selected_runtyp_flag)

         case ('--alpha')
            call set_elprop('alpha')

         case ('--raman')
            ! `--raman` is a run type selection too (`set_elprop` runs the
            ! `hess` path), but the selection happens inside `set_elprop`, so
            ! the wrapper only records it (review finding M1 of
            ! `20-stage1-review.md`).
            call select_runtyp(env, 'hess', '--raman', selected_runtyp, &
               & selected_runtyp_flag, apply=.false.)
            call set_elprop('raman')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_raman(env,sec)
               call args%nextArg(sec)
               if (allocated(sec)) then
                  call set_raman(env,sec)
               endif
            endif


         case ('--grad')
            call select_runtyp(env, 'grad', '--grad', selected_runtyp, selected_runtyp_flag)
            lgrad = .true.

         case ('-o', '--opt')
            call select_runtyp(env, 'opt', '--opt', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_opt(env, 'optlevel', sec)
               ! Keep the raw keyword for the suite bridge of the tail; the
               ! integer above cannot tell `banana` from `normal` and cannot
               ! carry the suite-only `superloose` (section 3.5, F3)
               landscape_opt_level_name = sec
            end if

         case ('--hess')
            call select_runtyp(env, 'hess', '--hess', selected_runtyp, selected_runtyp_flag)

         case ('--md')
            call select_runtyp(env, 'md', '--md', selected_runtyp, selected_runtyp_flag)

         case ('--ohess')
            call select_runtyp(env, 'ohess', '--ohess', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_opt(env, 'optlevel', sec)
               ! Same bridge input as `--opt` (CLI redesign section 3.5)
               landscape_opt_level_name = sec
            end if

         case ('--bhess')
            call select_runtyp(env, 'bhess', '--bhess', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_opt(env, 'optlevel', sec)
               ! The native `optlevel2int` maps an unknown keyword to `normal`,
               ! so the keyword is reported in the parser tail (CLI section 3.5)
               landscape_native_level_name = sec
               landscape_native_level_flag = '--bhess'
            end if

         case ('--o1nh')
            set%o1numhess = .true.

         case ('--imagmin')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_hess(env, 'imagmin', '-'//sec)
            else
               call env%error("Imaginary significance floor is missing", source)
            end if

         case ('--imagmax')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_hess(env, 'imagmax', '-'//sec)
            else
               call env%error("Imaginary repair cutoff is missing", source)
            end if

         case ('--omd')
            call select_runtyp(env, 'omd', '--omd', selected_runtyp, selected_runtyp_flag)
            call set_opt(env, 'optlevel', '-1')

         case ('--siman')
            call select_runtyp(env, 'siman', '--siman', selected_runtyp, selected_runtyp_flag)
            call set_md(env, 'nvt', 'true')

         case ('--path')
            call select_runtyp(env, 'path', '--path', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_path(env, 'product', sec)
            end if

         case ('--screen')
            call select_runtyp(env, 'screen', '--screen', selected_runtyp, selected_runtyp_flag)

         case ('--gmd')
            call select_runtyp(env, 'gmd', '--gmd', selected_runtyp, selected_runtyp_flag)
            call env%error("This feature has been deprecated, I'm sorry.", source)

         case ('--modef')
            call select_runtyp(env, 'modef', '--modef', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_modef(env, 'mode', sec)
            end if

         case ('--mdopt')
            call select_runtyp(env, 'mdopt', '--mdopt', selected_runtyp, selected_runtyp_flag)

         case ('--ljson')
            ! Structured result of the landscape suite: every task writes
            ! `landscape_<task>.json` (`status` / `counters` / convergence
            ! `history`) independently of `--l*` and of the geometry files
            ! (PLAN section 5.5).  The option is valid for all five run
            ! types, so it is set on every settings object and read by
            ! whichever task the run type selects.
            lset%ljson = .true.
            lsset%ljson = .true.
            lvset%ljson = .true.
            ltset%ljson = .true.
            lircset%ljson = .true.

         case ('--lopt')
            ! Landscape suite: minimisation with the algorithms of
            ! `src/landscape/optimize` (see PLAN section 12).  The flag takes
            ! no positional argument; engine and convergence level are set by
            ! the dedicated flags below.
            call select_runtyp(env, 'lopt', '--lopt', selected_runtyp, selected_runtyp_flag)

         case ('--method')
            ! Unified method selector of the landscape suite (CLI redesign
            ! section 4).  Only the name is stored here: it is resolved after
            ! the parsing loop, where the run type is known regardless of the
            ! order of the flags (the run type is set inside the run type
            ! cases, so `--method bfgs --lopt` would not see `lopt` here).
            !
            ! `auto` is the default of the selected task; a name of the native
            ! ANCopt family (`ancopt`/`rf`, `lancopt`/`l-ancopt`, `fire`/
            ! `inertial`, `pbc_lbfgs`) selects the native engine and is legal
            ! with `--opt`/`--ohess` only (section 4.7/4.8).
            call args%nextArg(sec)
            if (allocated(sec)) then
               landscape_method_name = sec
            else
               call env%error("The --method flag requires a method name", source)
            end if

         case ('--lopt-engine')
            ! Minimisation algorithm of the landscape suite; alias of
            ! `--method` (CLI redesign section 4.3), the shared setter does the
            ! name lookup and the alias conflict check
            call select_runtyp(env, 'lopt', '--lopt-engine', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_landscape_opt_algorithm(env, lset, sec, '--lopt-engine', source)
            else
               call env%error("The --lopt-engine flag requires an algorithm name", source)
            end if

         case ('--lopt-level')
            ! Convergence level of the landscape suite (ORCA naming)
            call select_runtyp(env, 'lopt', '--lopt-level', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               if (landscape_conv_level_known(sec)) then
                  call lset%set_conv_level(sec)
               else
                  call env%error("Unknown landscape convergence level '"//sec//"'", source)
               end if
            else
               call env%error("The --lopt-level flag requires a level name", source)
            end if

         case ('--lopt-hessian')
            ! Hessian backend of the landscape suite.  `numerical` is the
            ! self-built central-difference Hessian, `model` the cheap
            ! Lindh-D2 model Hessian (`src/landscape/types/hessian.f90`).
            call select_runtyp(env, 'lopt', '--lopt-hessian', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               if (landscape_hessian_source_known(sec)) then
                  lset%hessian_source = landscape_lowercase(sec)
               else
                  call env%error("Unknown landscape Hessian source '"//sec//"'", source)
               end if
            else
               call env%error("The --lopt-hessian flag requires a source name", source)
            end if

         case ('--lopt-max-steps')
            ! Step budget of the minimisation (CLI redesign section 6.3).  The
            ! field existed with a fixed default of 256, it was simply not
            ! reachable from the command line; the mirror spelling of the scan
            ! is `--lscan-max-steps`.
            call select_runtyp(env, 'lopt', '--lopt-max-steps', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  lset%max_steps = idum2
                  ! Record the explicit budget: it wins over the one the
                  ! convergence level derives (`level_step_budget`).
                  lset%max_steps_source = '--lopt-max-steps'
                  lopt_max_steps_given = .true.
               else
                  call env%error("The --lopt-max-steps flag requires an integer", source)
               end if
            else
               call env%error("The --lopt-max-steps flag requires an integer", source)
            end if

         case ('--lopt-recalc')
            ! Hessian rebuild interval of the RFO step.  The field existed
            ! with a fixed default of 5 and no flag reached it; the saddle
            ! (`--lts-recalc`) and reaction path (`--lirc-recalc`) spellings
            ! are the same knob of their own settings object.  `0` disables
            ! the rebuild and an interval above the step budget simply never
            ! fires; a negative value is reported in the tail, where the
            ! effect of `--strict` is known.
            call select_runtyp(env, 'lopt', '--lopt-recalc', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  lset%recalc = idum2
                  lopt_recalc_value = idum2
                  lopt_recalc_given = .true.
               else
                  call env%error("The --lopt-recalc flag requires an integer", source)
               end if
            else
               call env%error("The --lopt-recalc flag requires an integer", source)
            end if

         case ('--lopt-hess-step')
            ! Finite difference step of the numerical Hessian (bohr, CLI
            ! redesign section 6.3); the minimisation had no knob for it, the
            ! field kept its default.  The scan mirror is `--lscan-hess-step`.
            call select_runtyp(env, 'lopt', '--lopt-hess-step', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  call set_landscape_hess_step(env, lset%engine, ddum, &
                     & '--lopt-hess-step', source)
               else
                  call env%error("The --lopt-hess-step flag requires a number", source)
               end if
            else
               call env%error("The --lopt-hess-step flag requires a number", source)
            end if

         case ('--lopt-sccacc')
            ! SCC accuracy of the displaced single points of the numerical
            ! Hessian (CLI redesign section 6.3); the scan mirror is
            ! `--lscan-sccacc`
            call select_runtyp(env, 'lopt', '--lopt-sccacc', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  lset%engine%hess_acc = ddum
                  lset%engine%hess_acc_source = '--lopt-sccacc'
               else
                  call env%error("The --lopt-sccacc flag requires a number", source)
               end if
            else
               call env%error("The --lopt-sccacc flag requires a number", source)
            end if

         case ('--lopt-restart')
            ! Wavefunction restart policy of the engine (CLI redesign section
            ! 6.3); the minimisation mirror of `--lscan-restart`
            call select_runtyp(env, 'lopt', '--lopt-restart', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               select case (landscape_lowercase(sec))
               case ('off')
                  lset%engine%restart_mode = p_landscape_restart_off
               case ('pin')
                  lset%engine%restart_mode = p_landscape_restart_pin
               case ('warm')
                  lset%engine%restart_mode = p_landscape_restart_warm
               case default
                  call env%error("Unknown landscape restart policy '"//sec//"'", source)
               end select
            else
               call env%error("The --lopt-restart flag requires off, pin or warm", source)
            end if

         case ('--lscan')
            ! Landscape suite: relaxed or rigid scan (PLAN section 12).  The
            ! flag takes no positional argument; the scanned coordinates and
            ! the grid are set by the dedicated `--lscan-*` flags below.
            call select_runtyp(env, 'lscan', '--lscan', selected_runtyp, selected_runtyp_flag)

         case ('--lscan-mode')
            ! Scan mode: relaxed or rigid
            call select_runtyp(env, 'lscan', '--lscan-mode', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               if (landscape_scan_mode_known(sec)) then
                  lsset%mode = landscape_lowercase(sec)
               else
                  call env%error("Unknown landscape scan mode '"//sec//"'", source)
               end if
            else
               call env%error("The --lscan-mode flag requires 'relaxed' or 'rigid'", source)
            end if

         case ('--lscan-method')
            ! Inner minimiser of a relaxed scan; alias of `--method` (CLI
            ! redesign section 4.3)
            call select_runtyp(env, 'lscan', '--lscan-method', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_landscape_scan_algorithm(env, lsset, sec, '--lscan-method', source)
            else
               call env%error("The --lscan-method flag requires an algorithm name", source)
            end if

         case ('--lscan-level')
            ! Convergence level of the inflated minimisation (ORCA naming)
            call select_runtyp(env, 'lscan', '--lscan-level', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               if (landscape_conv_level_known(sec)) then
                  call lsset%set_conv_level(sec)
               else
                  call env%error("Unknown landscape convergence level '"//sec//"'", source)
               end if
            else
               call env%error("The --lscan-level flag requires a level name", source)
            end if

         case ('--lscan-restart')
            ! Wavefunction restart policy of the scan point evaluations
            call select_runtyp(env, 'lscan', '--lscan-restart', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               select case (landscape_lowercase(sec))
               case ('off')
                  lsset%engine%restart_mode = p_landscape_restart_off
               case ('pin')
                  lsset%engine%restart_mode = p_landscape_restart_pin
               case ('warm')
                  lsset%engine%restart_mode = p_landscape_restart_warm
               case default
                  call env%error("Unknown landscape restart policy '"//sec//"'", source)
               end select
            else
               call env%error("The --lscan-restart flag requires off, pin or warm", source)
            end if

         case ('--lscan-ndim')
            ! Number of scanned coordinates (1 or 2)
            call select_runtyp(env, 'lscan', '--lscan-ndim', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  if (idum2 == 1 .or. idum2 == 2) then
                     lsset%ndim = idum2
                  else
                     call env%error("The --lscan-ndim flag accepts 1 or 2", source)
                  end if
               else
                  call env%error("The --lscan-ndim flag requires an integer", source)
               end if
            else
               call env%error("The --lscan-ndim flag requires 1 or 2", source)
            end if

         case ('--lscan-oneway')
            ! Scan only the forward branch (no hysteresis report)
            call select_runtyp(env, 'lscan', '--lscan-oneway', selected_runtyp, selected_runtyp_flag)
            lsset%bidirectional = .false.

         case ('--lscan-cold-start')
            ! Start every scan point from the input geometry
            call select_runtyp(env, 'lscan', '--lscan-cold-start', selected_runtyp, selected_runtyp_flag)
            lsset%warm_start = .false.

         case ('--lscan-fc')
            ! Force constant of the constraint bias (hartree/bohr^2)
            call select_runtyp(env, 'lscan', '--lscan-fc', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  lsset%constraint_fc = ddum
                  lsset%constraint_fc_was_set = .true.
               else
                  call env%error("The --lscan-fc flag requires a number", source)
               end if
            else
               call env%error("The --lscan-fc flag requires a number", source)
            end if

         case ('--lscan-expo')
            ! Exponent of the constraint bias potential
            call select_runtyp(env, 'lscan', '--lscan-expo', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  lsset%spring_expo = ddum
               else
                  call env%error("The --lscan-expo flag requires a number", source)
               end if
            else
               call env%error("The --lscan-expo flag requires a number", source)
            end if

         case ('--lscan-coord')
            ! Coordinate type of the first scanned coordinate
            call select_runtyp(env, 'lscan', '--lscan-coord', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               iscan = landscape_scan_coord_from_name(sec)
               if (iscan == 0) then
                  call env%error("Unknown landscape scan coordinate '"//sec//"'", source)
               else
                  lsset%coord(1)%ctype = iscan
               end if
            else
               call env%error("The --lscan-coord flag requires distance, angle or dihedral", source)
            end if

         case ('--lscan-atoms')
            ! 1-based atom indices of the first scanned coordinate
            call select_runtyp(env, 'lscan', '--lscan-atoms', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               call read_scan_atoms(env, sec, lsset%coord(1)%atoms, source)
            else
               call env%error("The --lscan-atoms flag requires a list of atom indices", source)
            end if

         case ('--lscan-step')
            ! Step of the first coordinate in angstrom (distance) or degree
            call select_runtyp(env, 'lscan', '--lscan-step', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  scan_step_user(1) = ddum
                  scan_step_given(1) = .true.
               else
                  call env%error("The --lscan-step flag requires a number", source)
               end if
            else
               call env%error("The --lscan-step flag requires a number", source)
            end if

         case ('--lscan-nstep')
            ! Number of steps of the first coordinate
            call select_runtyp(env, 'lscan', '--lscan-nstep', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  lsset%coord(1)%nstep = idum2
               else
                  call env%error("The --lscan-nstep flag requires an integer", source)
               end if
            else
               call env%error("The --lscan-nstep flag requires an integer", source)
            end if

         case ('--lscan-coord2')
            ! Coordinate type of the second scanned coordinate (sets ndim=2)
            call select_runtyp(env, 'lscan', '--lscan-coord2', selected_runtyp, selected_runtyp_flag)
            lsset%ndim = 2
            call args%nextArg(sec)
            if (allocated(sec)) then
               iscan = landscape_scan_coord_from_name(sec)
               if (iscan == 0) then
                  call env%error("Unknown landscape scan coordinate '"//sec//"'", source)
               else
                  lsset%coord(2)%ctype = iscan
               end if
            else
               call env%error("The --lscan-coord2 flag requires distance, angle or dihedral", source)
            end if

         case ('--lscan-atoms2')
            ! 1-based atom indices of the second scanned coordinate
            call select_runtyp(env, 'lscan', '--lscan-atoms2', selected_runtyp, selected_runtyp_flag)
            lsset%ndim = 2
            call args%nextArg(sec)
            if (allocated(sec)) then
               call read_scan_atoms(env, sec, lsset%coord(2)%atoms, source)
            else
               call env%error("The --lscan-atoms2 flag requires a list of atom indices", source)
            end if

         case ('--lscan-step2')
            ! Step of the second coordinate in angstrom (distance) or degree
            call select_runtyp(env, 'lscan', '--lscan-step2', selected_runtyp, selected_runtyp_flag)
            lsset%ndim = 2
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  scan_step_user(2) = ddum
                  scan_step_given(2) = .true.
               else
                  call env%error("The --lscan-step2 flag requires a number", source)
               end if
            else
               call env%error("The --lscan-step2 flag requires a number", source)
            end if

         case ('--lscan-nstep2')
            ! Number of steps of the second coordinate
            call select_runtyp(env, 'lscan', '--lscan-nstep2', selected_runtyp, selected_runtyp_flag)
            lsset%ndim = 2
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  lsset%coord(2)%nstep = idum2
               else
                  call env%error("The --lscan-nstep2 flag requires an integer", source)
               end if
            else
               call env%error("The --lscan-nstep2 flag requires an integer", source)
            end if

         case ('--lscan-max-steps')
            ! Step budget of the inner relaxation of every scan point
            call select_runtyp(env, 'lscan', '--lscan-max-steps', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  lsset%max_steps = idum2
                  lsset%max_steps_source = '--lscan-max-steps'
                  lscan_max_steps_given = .true.
               else
                  call env%error("The --lscan-max-steps flag requires an integer", source)
               end if
            else
               call env%error("The --lscan-max-steps flag requires an integer", source)
            end if

         case ('--lscan-hess-step')
            ! Finite difference step of the numerical Hessian of the scan
            ! points (bohr, CLI redesign sections 6.3 and 7.2); the canonical
            ! spelling, there is no older one for the scan
            call select_runtyp(env, 'lscan', '--lscan-hess-step', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  call set_landscape_hess_step(env, lsset%engine, ddum, &
                     & '--lscan-hess-step', source)
               else
                  call env%error("The --lscan-hess-step flag requires a number", source)
               end if
            else
               call env%error("The --lscan-hess-step flag requires a number", source)
            end if

         case ('--lscan-sccacc')
            ! SCC accuracy of the displaced single points of the numerical
            ! Hessian of the scan points (CLI redesign section 6.3)
            call select_runtyp(env, 'lscan', '--lscan-sccacc', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  lsset%engine%hess_acc = ddum
                  lsset%engine%hess_acc_source = '--lscan-sccacc'
               else
                  call env%error("The --lscan-sccacc flag requires a number", source)
               end if
            else
               call env%error("The --lscan-sccacc flag requires a number", source)
            end if

         case ('--lscan-hessian')
            ! Hessian backend of the scan points (`numerical`, `model` or
            ! `file`); needed for `--lscan-method rfo` and an inner BFGS
            ! start Hessian
            call select_runtyp(env, 'lscan', '--lscan-hessian', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               if (landscape_hessian_source_known(sec)) then
                  lsset%hessian_source = landscape_lowercase(sec)
               else
                  call env%error("Unknown landscape Hessian source '"//sec//"'", source)
               end if
            else
               call env%error("The --lscan-hessian flag requires a source name", source)
            end if

         case ('--lscan-hessian-file')
            ! File read by `--lscan-hessian file`
            call select_runtyp(env, 'lscan', '--lscan-hessian-file', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               lsset%hessian_file = sec
            else
               call env%error("The --lscan-hessian-file flag requires a file name", source)
            end if

         case ('--lts')
            ! Landscape suite: transition state search (src/landscape/saddle);
            ! the method, trust radius and Hessian backend are set by the
            ! dedicated `--lts-*` flags below.
            call select_runtyp(env, 'lts', '--lts', selected_runtyp, selected_runtyp_flag)

         case ('--lts-method')
            ! Transition state search method; alias of `--method` (CLI redesign
            ! section 4.3).  The shared setter carries the climbing image side
            ! effect, so `--lts --method cineb` and `--lts --lts-method cineb`
            ! are the same run.
            call select_runtyp(env, 'lts', '--lts-method', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_landscape_ts_algorithm(env, ltset, sec, '--lts-method', source)
            else
               call env%error("The --lts-method flag requires a method name", source)
            end if

         case ('--lts-trust')
            ! Initial trust radius of the P-RFO step (bohr)
            call select_runtyp(env, 'lts', '--lts-trust', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  ltset%trust = ddum
               else
                  call env%error("The --lts-trust flag requires a number", source)
               end if
            else
               call env%error("The --lts-trust flag requires a number", source)
            end if

         case ('--lts-trust-max')
            ! Upper bound of the trust radius (bohr)
            call select_runtyp(env, 'lts', '--lts-trust-max', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  ltset%trust_max = ddum
               else
                  call env%error("The --lts-trust-max flag requires a number", source)
               end if
            else
               call env%error("The --lts-trust-max flag requires a number", source)
            end if

         case ('--lts-hlow')
            ! Curvature floor of the downhill modes (hartree/bohr^2)
            call select_runtyp(env, 'lts', '--lts-hlow', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  ltset%hlow = ddum
               else
                  call env%error("The --lts-hlow flag requires a number", source)
               end if
            else
               call env%error("The --lts-hlow flag requires a number", source)
            end if

         case ('--lts-max-steps')
            ! Maximum number of macro iterations
            call select_runtyp(env, 'lts', '--lts-max-steps', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  ltset%max_steps = idum2
                  ltset%max_steps_source = '--lts-max-steps'
                  lts_max_steps_given = .true.
               else
                  call env%error("The --lts-max-steps flag requires an integer", source)
               end if
            else
               call env%error("The --lts-max-steps flag requires an integer", source)
            end if

         case ('--lts-level')
            ! Convergence thresholds of the search
            call select_runtyp(env, 'lts', '--lts-level', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               if (landscape_conv_level_known(sec)) then
                  call ltset%set_conv_level(sec)
               else
                  call env%error("Unknown landscape convergence level '"//sec//"'", source)
               end if
            else
               call env%error("The --lts-level flag requires a level name", source)
            end if

         case ('--lts-mode-overlap-min')
            ! Mode following guard threshold
            call select_runtyp(env, 'lts', '--lts-mode-overlap-min', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  ltset%mode_overlap_min = ddum
               else
                  call env%error("The --lts-mode-overlap-min flag requires a number", source)
               end if
            else
               call env%error("The --lts-mode-overlap-min flag requires a number", source)
            end if

         case ('--lts-alpha')
            ! Uphill/downhill trust radius fraction of the B-P-RFO
            call select_runtyp(env, 'lts', '--lts-alpha', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  ltset%alpha = ddum
               else
                  call env%error("The --lts-alpha flag requires a number", source)
               end if
            else
               call env%error("The --lts-alpha flag requires a number", source)
            end if

         case ('--lts-micro-cycles')
            ! Maximum number of restricted-step micro iterations
            call select_runtyp(env, 'lts', '--lts-micro-cycles', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  ltset%max_micro_cycles = idum2
               else
                  call env%error("The --lts-micro-cycles flag requires an integer", source)
               end if
            else
               call env%error("The --lts-micro-cycles flag requires an integer", source)
            end if

         case ('--lts-target-mode')
            ! Which negative mode is followed (1-based)
            call select_runtyp(env, 'lts', '--lts-target-mode', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  ltset%target_mode = idum2
               else
                  call env%error("The --lts-target-mode flag requires an integer", source)
               end if
            else
               call env%error("The --lts-target-mode flag requires an integer", source)
            end if

         case ('--lts-init-mode')
            ! Source of the initial transition vector
            call select_runtyp(env, 'lts', '--lts-init-mode', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               if (landscape_ts_init_mode_known(sec)) then
                  ltset%init_mode = landscape_lowercase(sec)
               else
                  call env%error("Unknown landscape initial mode '"//sec//"'", source)
               end if
            else
               call env%error("The --lts-init-mode flag requires 'h0' or 'force'", source)
            end if

         case ('--lts-delta')
            ! Dimer separation, the distance of each image from the midpoint
            ! (bohr)
            call select_runtyp(env, 'lts', '--lts-delta', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  ltset%delta = ddum
               else
                  call env%error("The --lts-delta flag requires a number", source)
               end if
            else
               call env%error("The --lts-delta flag requires a number", source)
            end if

         case ('--lts-dr')
            ! Finite difference step of the dimer translation probe (bohr)
            call select_runtyp(env, 'lts', '--lts-dr', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  ltset%dR = ddum
               else
                  call env%error("The --lts-dr flag requires a number", source)
               end if
            else
               call env%error("The --lts-dr flag requires a number", source)
            end if

         case ('--lts-rot-max-iter')
            ! Maximum number of dimer rotations per macro iteration
            call select_runtyp(env, 'lts', '--lts-rot-max-iter', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  ltset%rot_max_iter = idum2
               else
                  call env%error("The --lts-rot-max-iter flag requires an integer", source)
               end if
            else
               call env%error("The --lts-rot-max-iter flag requires an integer", source)
            end if

         case ('--lts-rot-alpha')
            ! Damping factor of the dimer rotation angle (Eq. 13)
            call select_runtyp(env, 'lts', '--lts-rot-alpha', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  ltset%rot_alpha = ddum
               else
                  call env%error("The --lts-rot-alpha flag requires a number", source)
               end if
            else
               call env%error("The --lts-rot-alpha flag requires a number", source)
            end if

         case ('--lts-rot-frms')
            ! RMS torque below which the dimer rotation is converged
            ! (hartree/bohr)
            call select_runtyp(env, 'lts', '--lts-rot-frms', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  ltset%rot_f_rms_th = ddum
               else
                  call env%error("The --lts-rot-frms flag requires a number", source)
               end if
            else
               call env%error("The --lts-rot-frms flag requires a number", source)
            end if

         case ('--lts-ninit')
            ! Source of the initial dimer axis
            call select_runtyp(env, 'lts', '--lts-ninit', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               if (landscape_ts_ninit_known(sec)) then
                  ltset%n_init = landscape_lowercase(sec)
               else
                  call env%error("Unknown landscape dimer initial axis '"//sec//"'", source)
               end if
            else
               call env%error("The --lts-ninit flag requires 'random' or 'force'", source)
            end if

         case ('--lts-seed')
            ! Seed of the random initial dimer axis (PLAN N10)
            call select_runtyp(env, 'lts', '--lts-seed', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  ltset%seed = idum2
               else
                  call env%error("The --lts-seed flag requires an integer", source)
               end if
            else
               call env%error("The --lts-seed flag requires an integer", source)
            end if

         case ('--lts-step-max')
            ! Hard bound of the dimer translation displacement (bohr)
            call select_runtyp(env, 'lts', '--lts-step-max', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  ltset%step_max = ddum
               else
                  call env%error("The --lts-step-max flag requires a number", source)
               end if
            else
               call env%error("The --lts-step-max flag requires a number", source)
            end if

         case ('--lts-translate-max-iter')
            ! Maximum number of quick-min displacement updates per dimer
            ! translation
            call select_runtyp(env, 'lts', '--lts-translate-max-iter', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  ltset%translate_max_iter = idum2
               else
                  call env%error("The --lts-translate-max-iter flag requires an integer", source)
               end if
            else
               call env%error("The --lts-translate-max-iter flag requires an integer", source)
            end if

         case ('--lts-recalc')
            ! Exact Hessian every N macro steps (0 disables)
            call select_runtyp(env, 'lts', '--lts-recalc', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  ltset%recalc = idum2
               else
                  call env%error("The --lts-recalc flag requires an integer", source)
               end if
            else
               call env%error("The --lts-recalc flag requires an integer", source)
            end if

         case ('--lts-no-verify')
            ! Skip the frequency confirmation of the stationary point
            call select_runtyp(env, 'lts', '--lts-no-verify', selected_runtyp, selected_runtyp_flag)
            ltset%verify_ts = .false.
         case ('--lts-no-rigid')
            ! Keep the rigid body components in the Cartesian step
            call select_runtyp(env, 'lts', '--lts-no-rigid', selected_runtyp, selected_runtyp_flag)
            ltset%remove_rigid = .false.

         case ('--lts-restart')
            ! Wavefunction restart policy of the engine
            call select_runtyp(env, 'lts', '--lts-restart', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               select case (landscape_lowercase(sec))
               case('off')
                  ltset%engine%restart_mode = p_landscape_restart_off
               case('pin')
                  ltset%engine%restart_mode = p_landscape_restart_pin
               case('warm')
                  ltset%engine%restart_mode = p_landscape_restart_warm
               case default
                  call env%error("The --lts-restart flag requires off, pin or warm", source)
               end select
            else
               call env%error("The --lts-restart flag requires off, pin or warm", source)
            end if

         case ('--lts-hessian')
            ! Hessian backend of the search
            call select_runtyp(env, 'lts', '--lts-hessian', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               if (landscape_hessian_source_known(sec)) then
                  ltset%hessian_source = landscape_lowercase(sec)
               else
                  call env%error("Unknown landscape Hessian source '"//sec//"'", source)
               end if
            else
               call env%error("The --lts-hessian flag requires a source name", source)
            end if

         case ('--lts-hessian-file')
            ! File read by `--lts-hessian file`
            call select_runtyp(env, 'lts', '--lts-hessian-file', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               ltset%hessian_file = sec
            else
               call env%error("The --lts-hessian-file flag requires a file name", source)
            end if

         case ('--lts-hess-step')
            ! Canonical spelling of the finite difference step of the
            ! confirmation Hessian (bohr, CLI redesign section 7.2);
            ! `--lts-step` is the older alias and both write the same field
            ! through the shared setter
            call select_runtyp(env, 'lts', '--lts-hess-step', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  call set_landscape_hess_step(env, ltset%engine, ddum, &
                     & '--lts-hess-step', source)
               else
                  call env%error("The --lts-hess-step flag requires a number", source)
               end if
            else
               call env%error("The --lts-hess-step flag requires a number", source)
            end if

         case ('--lts-step')
            ! Finite difference step of the confirmation Hessian (bohr); alias
            ! of `--lts-hess-step` (CLI redesign sections 6.2 and 7.2), the
            ! shared setter applies the alias conflict rule.  Keep the client
            ! side of M4 hard constraint 2 (H1) working.
            call select_runtyp(env, 'lts', '--lts-step', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  call set_landscape_hess_step(env, ltset%engine, ddum, &
                     & '--lts-step', source)
               else
                  call env%error("The --lts-step flag requires a number", source)
               end if
            else
               call env%error("The --lts-step flag requires a number", source)
            end if

         case ('--lts-sccacc')
            ! SCC accuracy of the displaced single points of the confirmation
            ! Hessian (H1)
            call select_runtyp(env, 'lts', '--lts-sccacc', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  ltset%engine%hess_acc = ddum
                  ltset%engine%hess_acc_source = '--lts-sccacc'
               else
                  call env%error("The --lts-sccacc flag requires a number", source)
               end if
            else
               call env%error("The --lts-sccacc flag requires a number", source)
            end if

         case ('--lts-step0')
            ! Upper bound of the first dimer translation (bohr)
            call select_runtyp(env, 'lts', '--lts-step0', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  ltset%step0 = ddum
               else
                  call env%error("The --lts-step0 flag requires a number", source)
               end if
            else
               call env%error("The --lts-step0 flag requires a number", source)
            end if

         case ('--lts-escape-step')
            ! Displacement along the lowest internal mode when the start
            ! geometry is a minimum (MEDIUM-1 of the review round 2); 0
            ! disables the escape.
            call select_runtyp(env, 'lts', '--lts-escape-step', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  ltset%escape_step = ddum
                  ltset%escape_step_was_set = .true.
               else
                  call env%error("The --lts-escape-step flag requires a number", source)
               end if
            else
               call env%error("The --lts-escape-step flag requires a number", source)
            end if

         case ('--lts-product')
            ! Product endpoint of a band search (`--lts-method neb|cineb`):
            ! a single frame xyz with the same atoms in the same order.
            call select_runtyp(env, 'lts', '--lts-product', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               ltset%path%product_file = sec
            else
               call env%error("The --lts-product flag requires a file name", source)
            end if

         case ('--lts-images')
            ! Number of band images, the two fixed endpoints included
            call select_runtyp(env, 'lts', '--lts-images', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  ltset%path%n_images = idum2
               else
                  call env%error("The --lts-images flag requires an integer", source)
               end if
            else
               call env%error("The --lts-images flag requires an integer", source)
            end if

         case ('--lts-kmin')
            ! Lower spring constant of the dynamic springs (hartree/bohr^2)
            call select_runtyp(env, 'lts', '--lts-kmin', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  ltset%path%k_min = ddum
               else
                  call env%error("The --lts-kmin flag requires a number", source)
               end if
            else
               call env%error("The --lts-kmin flag requires a number", source)
            end if

         case ('--lts-kmax')
            ! Upper spring constant (hartree/bohr^2)
            call select_runtyp(env, 'lts', '--lts-kmax', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  ltset%path%k_max = ddum
               else
                  call env%error("The --lts-kmax flag requires a number", source)
               end if
            else
               call env%error("The --lts-kmax flag requires a number", source)
            end if

         case ('--lts-kdecay')
            ! Decay factor of the dynamic springs (dimensionless)
            call select_runtyp(env, 'lts', '--lts-kdecay', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  ltset%path%k_decay = ddum
               else
                  call env%error("The --lts-kdecay flag requires a number", source)
               end if
            else
               call env%error("The --lts-kdecay flag requires a number", source)
            end if

         case ('--lts-no-dynamic-k')
            ! Use a single spring constant k_max for every image
            call select_runtyp(env, 'lts', '--lts-no-dynamic-k', selected_runtyp, selected_runtyp_flag)
            ltset%path%use_dynamic_k = .false.

         case ('--lts-climb')
            ! Switch on the climbing image refinement (CI-NEB)
            call select_runtyp(env, 'lts', '--lts-climb', selected_runtyp, selected_runtyp_flag)
            ltset%path%climb = .true.

         case ('--lts-climb-switch')
            ! Band force below which the climbing image is ignited
            ! (hartree/bohr); the default is the standard 0.5 eV/angstrom
            call select_runtyp(env, 'lts', '--lts-climb-switch', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  ltset%path%climb_switch_rms = ddum
               else
                  call env%error("The --lts-climb-switch flag requires a number", source)
               end if
            else
               call env%error("The --lts-climb-switch flag requires a number", source)
            end if

         case ('--lts-no-idpp')
            ! Straight linear interpolation instead of the IDPP initial path
            call select_runtyp(env, 'lts', '--lts-no-idpp', selected_runtyp, selected_runtyp_flag)
            ltset%path%ifidpp = .false.

         case ('--lts-no-align')
            ! Do not Kabsch align the images before every evaluation
            call select_runtyp(env, 'lts', '--lts-no-align', selected_runtyp, selected_runtyp_flag)
            ltset%path%align = .false.

         case ('--lts-relax-ends')
            ! Relax the two endpoints before the band is built
            call select_runtyp(env, 'lts', '--lts-relax-ends', selected_runtyp, selected_runtyp_flag)
            ltset%path%relax_ends = .true.

         case ('--lts-max-displ')
            ! Upper bound of one displacement component of a band step (bohr)
            call select_runtyp(env, 'lts', '--lts-max-displ', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  ltset%path%max_displ = ddum
               else
                  call env%error("The --lts-max-displ flag requires a number", source)
               end if
            else
               call env%error("The --lts-max-displ flag requires a number", source)
            end if

         case ('--lts-frms-avg')
            ! Convergence threshold of the average band force RMS
            call select_runtyp(env, 'lts', '--lts-frms-avg', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  ltset%path%f_rms_avg_tol = ddum
               else
                  call env%error("The --lts-frms-avg flag requires a number", source)
               end if
            else
               call env%error("The --lts-frms-avg flag requires a number", source)
            end if

         case ('--lts-frms-max')
            ! Convergence threshold of the largest per-image band force RMS
            call select_runtyp(env, 'lts', '--lts-frms-max', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  ltset%path%f_rms_max_tol = ddum
               else
                  call env%error("The --lts-frms-max flag requires a number", source)
               end if
            else
               call env%error("The --lts-frms-max flag requires a number", source)
            end if

         case ('--lts-autoneb-nstart')
            ! Number of images of the coarse AutoNEB start path
            call select_runtyp(env, 'lts', '--lts-autoneb-nstart', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  ltset%path%autoneb_nstart = idum2
               else
                  call env%error("The --lts-autoneb-nstart flag requires an integer", source)
               end if
            else
               call env%error("The --lts-autoneb-nstart flag requires an integer", source)
            end if

         case ('--lts-autoneb-max-paths')
            ! Hard upper bound of the number of images of an AutoNEB path
            call select_runtyp(env, 'lts', '--lts-autoneb-max-paths', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  ltset%path%max_paths = idum2
               else
                  call env%error("The --lts-autoneb-max-paths flag requires an integer", source)
               end if
            else
               call env%error("The --lts-autoneb-max-paths flag requires an integer", source)
            end if

         case ('--lts-autoneb-max-depth')
            ! Hard upper bound of the number of AutoNEB insertion stages
            call select_runtyp(env, 'lts', '--lts-autoneb-max-depth', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  ltset%path%max_depth = idum2
               else
                  call env%error("The --lts-autoneb-max-depth flag requires an integer", source)
               end if
            else
               call env%error("The --lts-autoneb-max-depth flag requires an integer", source)
            end if

         case ('--lts-autoneb-restart')
            ! Multi-frame path of an earlier AutoNEB run (`path_autoneb.xyz`),
            ! adopted as the initial path so completed stages are not repeated
            call select_runtyp(env, 'lts', '--lts-autoneb-restart', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               ltset%path%autoneb_restart_file = sec
            else
               call env%error("The --lts-autoneb-restart flag requires a file name", source)
            end if

         case ('--lts-string-nnodes')
            ! Target number of nodes of the growing string
            call select_runtyp(env, 'lts', '--lts-string-nnodes', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  ltset%path%string_nnodes = idum2
               else
                  call env%error("The --lts-string-nnodes flag requires an integer", source)
               end if
            else
               call env%error("The --lts-string-nnodes flag requires an integer", source)
            end if

         case ('--lts-string-grow-tol')
            ! Norm threshold of the growth indicator of Eq. 7 (hartree/bohr)
            call select_runtyp(env, 'lts', '--lts-string-grow-tol', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  ltset%path%string_grow_tol = ddum
               else
                  call env%error("The --lts-string-grow-tol flag requires a number", source)
               end if
            else
               call env%error("The --lts-string-grow-tol flag requires a number", source)
            end if

         case ('--lts-string-ds')
            ! Arc length increment of one growth step (Eq. 9, bohr)
            call select_runtyp(env, 'lts', '--lts-string-ds', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  ltset%path%string_ds = ddum
               else
                  call env%error("The --lts-string-ds flag requires a number", source)
               end if
            else
               call env%error("The --lts-string-ds flag requires a number", source)
            end if

         case ('--lts-string-nu2')
            ! Curvature floor of the node spacing condition (hartree/bohr^2)
            call select_runtyp(env, 'lts', '--lts-string-nu2', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  ltset%path%string_nu2 = ddum
               else
                  call env%error("The --lts-string-nu2 flag requires a number", source)
               end if
            else
               call env%error("The --lts-string-nu2 flag requires a number", source)
            end if

         case ('--lts-string-grow-steps')
            ! Maximum perpendicular relaxation iterations of one frontier node
            call select_runtyp(env, 'lts', '--lts-string-grow-steps', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  ltset%path%string_grow_steps = idum2
               else
                  call env%error("The --lts-string-grow-steps flag requires an integer", source)
               end if
            else
               call env%error("The --lts-string-grow-steps flag requires an integer", source)
            end if

         case ('--lts-no-string-spline')
            ! Central difference tangents instead of the cubic spline tangent
            call select_runtyp(env, 'lts', '--lts-no-string-spline', selected_runtyp, selected_runtyp_flag)
            ltset%path%string_spline = .false.

         case ('--lts-string-history')
            ! Print the growth history of the string
            call select_runtyp(env, 'lts', '--lts-string-history', selected_runtyp, selected_runtyp_flag)
            ltset%path%string_history = .true.

         case ('--lirc')
            ! Landscape suite: intrinsic reaction coordinate
            ! (src/landscape/irc); the method and the integrator parameters are
            ! set by the dedicated `--lirc-*` flags below.  The input geometry
            ! is the transition state.
            call select_runtyp(env, 'lirc', '--lirc', selected_runtyp, selected_runtyp_flag)

         case ('--lirc-method')
            ! IRC integration method (hpc default | lqa | eulerpc | gs); the default
            ! is `hpc` since the M6 validation round measured it both more accurate
            ! and 8.5x cheaper than `gs` (PLAN section 12.15, behaviour change).
            ! Alias of `--method` (CLI redesign section 4.3)
            call select_runtyp(env, 'lirc', '--lirc-method', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_landscape_irc_algorithm(env, lircset, sec, '--lirc-method', source)
            else
               call env%error("The --lirc-method flag requires a method name", source)
            end if

         case ('--lirc-step')
            ! Step length of the IRC integration (bohr*amu^1/2, mass weighted
            ! arc length)
            call select_runtyp(env, 'lirc', '--lirc-step', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  lircset%step_length = ddum
               else
                  call env%error("The --lirc-step flag requires a number", source)
               end if
            else
               call env%error("The --lirc-step flag requires a number", source)
            end if

         case ('--lirc-max-steps')
            ! Maximum number of macro steps of one direction
            call select_runtyp(env, 'lirc', '--lirc-max-steps', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  lircset%max_steps = idum2
                  lircset%max_steps_source = '--lirc-max-steps'
                  lirc_max_steps_given = .true.
               else
                  call env%error("The --lirc-max-steps flag requires an integer", source)
               end if
            else
               call env%error("The --lirc-max-steps flag requires an integer", source)
            end if

         case ('--lirc-level')
            ! Convergence thresholds of the inner endpoint relaxation (the four
            ! minimisation criteria; the path criteria are --lirc-grad-tol and
            ! --lirc-e-rise, and the relaxation step budget is --lirc-max-steps)
            call select_runtyp(env, 'lirc', '--lirc-level', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               if (landscape_conv_level_known(sec)) then
                  call lircset%set_conv_level(sec)
               else
                  call env%error("Unknown landscape convergence level '"//sec//"'", source)
               end if
            else
               call env%error("The --lirc-level flag requires a level name", source)
            end if

         case ('--lirc-swap')
            ! Exchange the forward and backward directions
            call select_runtyp(env, 'lirc', '--lirc-swap', selected_runtyp, selected_runtyp_flag)
            lircset%swap = .true.

         case ('--lirc-oneway')
            ! Integrate the forward direction only
            call select_runtyp(env, 'lirc', '--lirc-oneway', selected_runtyp, selected_runtyp_flag)
            lircset%both_directions = .false.

         case ('--lirc-target-mode')
            ! Which negative mode defines the initial direction (1-based)
            call select_runtyp(env, 'lirc', '--lirc-target-mode', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  lircset%target_mode = idum2
               else
                  call env%error("The --lirc-target-mode flag requires an integer", source)
               end if
            else
               call env%error("The --lirc-target-mode flag requires an integer", source)
            end if

         case ('--lirc-micro-cycles')
            ! Maximum number of GS micro cycles per macro step
            call select_runtyp(env, 'lirc', '--lirc-micro-cycles', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  lircset%max_micro_cycles = idum2
               else
                  call env%error("The --lirc-micro-cycles flag requires an integer", source)
               end if
            else
               call env%error("The --lirc-micro-cycles flag requires an integer", source)
            end if

         case ('--lirc-micro-thresh')
            ! Convergence threshold of the sphere relaxation (bohr*amu^1/2)
            call select_runtyp(env, 'lirc', '--lirc-micro-thresh', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  lircset%micro_step_thresh = ddum
               else
                  call env%error("The --lirc-micro-thresh flag requires a number", source)
               end if
            else
               call env%error("The --lirc-micro-thresh flag requires a number", source)
            end if

         case ('--lirc-euler-n')
            ! Number of Euler substeps of the analytic LQA propagation
            call select_runtyp(env, 'lirc', '--lirc-euler-n', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  lircset%euler_n = idum2
               else
                  call env%error("The --lirc-euler-n flag requires an integer", source)
               end if
            else
               call env%error("The --lirc-euler-n flag requires an integer", source)
            end if

         case ('--lirc-loose-cycles')
            ! EulerPC coarse phase: the first N macro steps use the plain
            ! Euler prediction (no corrector), 0 disables the coarse phase
            call select_runtyp(env, 'lirc', '--lirc-loose-cycles', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  lircset%loose_cycles = idum2
               else
                  call env%error("The --lirc-loose-cycles flag requires an integer", source)
               end if
            else
               call env%error("The --lirc-loose-cycles flag requires an integer", source)
            end if

         case ('--lirc-mbs-tol')
            ! Acceptance threshold of the mBS corrector (bohr*amu^1/2); 1e-5 is
            ! the MAPLE default, 1e-6 the stricter HPC paper value
            call select_runtyp(env, 'lirc', '--lirc-mbs-tol', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  lircset%mbs_tol = ddum
               else
                  call env%error("The --lirc-mbs-tol flag requires a number", source)
               end if
            else
               call env%error("The --lirc-mbs-tol flag requires a number", source)
            end if

         case ('--lirc-mbs-max-k')
            ! Maximum extrapolation level of the mBS corrector
            call select_runtyp(env, 'lirc', '--lirc-mbs-max-k', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  lircset%mbs_max_k = idum2
               else
                  call env%error("The --lirc-mbs-max-k flag requires an integer", source)
               end if
            else
               call env%error("The --lirc-mbs-max-k flag requires an integer", source)
            end if

         case ('--lirc-mbs-points')
            ! Base number of Euler substeps of the mBS corrector
            call select_runtyp(env, 'lirc', '--lirc-mbs-points', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  lircset%mbs_points = idum2
               else
                  call env%error("The --lirc-mbs-points flag requires an integer", source)
               end if
            else
               call env%error("The --lirc-mbs-points flag requires an integer", source)
            end if

         case ('--lirc-dwi-n')
            ! Capacity of the DWI surface fit; the interpolant is two point by
            ! construction, so every value >= 2 is behaviourally equivalent
            ! (`check_irc_settings` rejects values below 2)
            call select_runtyp(env, 'lirc', '--lirc-dwi-n', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  lircset%dwi_n = idum2
               else
                  call env%error("The --lirc-dwi-n flag requires an integer", source)
               end if
            else
               call env%error("The --lirc-dwi-n flag requires an integer", source)
            end if

         case ('--lirc-e-rise', '--lirc-erise-steps')
            ! Termination signal (c): number of consecutive rising energies
            call select_runtyp(env, 'lirc', '--lirc-e-rise', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  lircset%e_rise_steps = idum2
               else
                  call env%error("The --lirc-e-rise flag requires an integer", source)
               end if
            else
               call env%error("The --lirc-e-rise flag requires an integer", source)
            end if

         case ('--lirc-grad-tol')
            ! Termination signal (b): gradient norm threshold (hartree/bohr)
            call select_runtyp(env, 'lirc', '--lirc-grad-tol', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  lircset%grad_tol = ddum
               else
                  call env%error("The --lirc-grad-tol flag requires a number", source)
               end if
            else
               call env%error("The --lirc-grad-tol flag requires a number", source)
            end if

         case ('--lirc-osc-window')
            ! Limit cycle guard of the GS macro step: leave the path after this
            ! many consecutive macro steps that cancel each other (0 disables)
            call select_runtyp(env, 'lirc', '--lirc-osc-window', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  lircset%osc_window = idum2
               else
                  call env%error("The --lirc-osc-window flag requires an integer", source)
               end if
            else
               call env%error("The --lirc-osc-window flag requires an integer", source)
            end if

         case ('--lirc-adaptive')
            ! Opt-in adaptive step of the GS macro step: halve the effective
            ! step length when consecutive macro steps cancel (default off)
            call select_runtyp(env, 'lirc', '--lirc-adaptive', selected_runtyp, selected_runtyp_flag)
            lircset%adaptive = .true.

         case ('--lirc-recalc')
            ! Rebuild the exact Hessian every N macro steps (-1 disables)
            call select_runtyp(env, 'lirc', '--lirc-recalc', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) idum2
               if (idum == 0) then
                  lircset%hessian_recalc = idum2
               else
                  call env%error("The --lirc-recalc flag requires an integer", source)
               end if
            else
               call env%error("The --lirc-recalc flag requires an integer", source)
            end if

         case ('--lirc-no-endpoint-relax')
            ! Skip the M1 relaxation of the two endpoints
            call select_runtyp(env, 'lirc', '--lirc-no-endpoint-relax', selected_runtyp, selected_runtyp_flag)
            lircset%endpoint_relax = .false.

         case ('--lirc-endpoint-relax')
            ! Perform the M1 relaxation of the two endpoints (default; the
            ! explicit flag restores it after a --lirc-no-endpoint-relax)
            call select_runtyp(env, 'lirc', '--lirc-endpoint-relax', selected_runtyp, selected_runtyp_flag)
            lircset%endpoint_relax = .true.

         case ('--lirc-rmsd-end')
            ! Threshold of the endpoint attribution: the Kabsch aligned RMSD
            ! between the endpoint and its relaxed minimum has to be below this
            ! value (`endpoint_reached`, PLAN section 6.5)
            call select_runtyp(env, 'lirc', '--lirc-rmsd-end', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  lircset%rmsd_end = ddum
               else
                  call env%error("The --lirc-rmsd-end flag requires a number", source)
               end if
            else
               call env%error("The --lirc-rmsd-end flag requires a number", source)
            end if

         case ('--lirc-hessian')
            ! Hessian backend of the IRC (initial Hessian and recalculations)
            call select_runtyp(env, 'lirc', '--lirc-hessian', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               if (landscape_hessian_source_known(sec)) then
                  lircset%hessian_source = landscape_lowercase(sec)
               else
                  call env%error("Unknown landscape Hessian source '"//sec//"'", source)
               end if
            else
               call env%error("The --lirc-hessian flag requires a source name", source)
            end if

         case ('--lirc-hessian-file')
            ! File read by `--lirc-hessian file`
            call select_runtyp(env, 'lirc', '--lirc-hessian-file', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               lircset%hessian_file = sec
            else
               call env%error("The --lirc-hessian-file flag requires a file name", source)
            end if

         case ('--lirc-hess-step')
            ! Canonical spelling of the finite difference step of the
            ! numerical Hessian (bohr, CLI redesign section 7.2); the older
            ! `--lirc-step-hess` is kept as an alias and both write the same
            ! field through the shared setter.  The canonical name separates
            ! the FD displacement from `--lirc-step`, the arc length of the
            ! reaction path
            call select_runtyp(env, 'lirc', '--lirc-hess-step', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  call set_landscape_hess_step(env, lircset%engine, ddum, &
                     & '--lirc-hess-step', source)
               else
                  call env%error("The --lirc-hess-step flag requires a number", source)
               end if
            else
               call env%error("The --lirc-hess-step flag requires a number", source)
            end if

         case ('--lirc-step-hess')
            ! Finite difference step of the numerical Hessian (bohr); alias of
            ! `--lirc-hess-step` (CLI redesign sections 6.2 and 7.2)
            call select_runtyp(env, 'lirc', '--lirc-step-hess', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  call set_landscape_hess_step(env, lircset%engine, ddum, &
                     & '--lirc-step-hess', source)
               else
                  call env%error("The --lirc-step-hess flag requires a number", source)
               end if
            else
               call env%error("The --lirc-step-hess flag requires a number", source)
            end if

         case ('--lirc-sccacc')
            ! SCC accuracy of the displaced single points of the Hessian
            call select_runtyp(env, 'lirc', '--lirc-sccacc', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  lircset%engine%hess_acc = ddum
                  lircset%engine%hess_acc_source = '--lirc-sccacc'
               else
                  call env%error("The --lirc-sccacc flag requires a number", source)
               end if
            else
               call env%error("The --lirc-sccacc flag requires a number", source)
            end if

         case ('--lirc-restart')
            ! Wavefunction restart policy of the engine
            call select_runtyp(env, 'lirc', '--lirc-restart', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               select case (landscape_lowercase(sec))
               case('off')
                  lircset%engine%restart_mode = p_landscape_restart_off
               case('pin')
                  lircset%engine%restart_mode = p_landscape_restart_pin
               case('warm')
                  lircset%engine%restart_mode = p_landscape_restart_warm
               case default
                  call env%error("The --lirc-restart flag requires off, pin or warm", source)
               end select
            else
               call env%error("The --lirc-restart flag requires off, pin or warm", source)
            end if

         case ('--lvib')
            ! Landscape suite: harmonic frequencies and thermochemistry
            ! (src/landscape/vib); the Hessian backend and the thermochemical
            ! parameters are set by the dedicated `--lvib-*` flags below.
            call select_runtyp(env, 'lvib', '--lvib', selected_runtyp, selected_runtyp_flag)

         case ('--lvib-hessian')
            ! Hessian backend of the frequency task
            call select_runtyp(env, 'lvib', '--lvib-hessian', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               if (landscape_hessian_source_known(sec)) then
                  lvset%hessian_source = landscape_lowercase(sec)
               else
                  call env%error("Unknown landscape Hessian source '"//sec//"'", source)
               end if
            else
               call env%error("The --lvib-hessian flag requires a source name", source)
            end if

         case ('--lvib-hessian-file')
            ! File read by `--lvib-hessian file` (TURBOMOLE layout, the
            ! `hessian` file of a preceding run)
            call select_runtyp(env, 'lvib', '--lvib-hessian-file', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               lvset%hessian_file = sec
            else
               call env%error("The --lvib-hessian-file flag requires a file name", source)
            end if

         case ('--lvib-hess-step')
            ! Canonical spelling of the finite difference step of the
            ! numerical Hessian (bohr, CLI redesign section 7.2);
            ! `--lvib-step` is the older alias and both write the same field
            ! through the shared setter
            call select_runtyp(env, 'lvib', '--lvib-hess-step', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  call set_landscape_hess_step(env, lvset%engine, ddum, &
                     & '--lvib-hess-step', source)
               else
                  call env%error("The --lvib-hess-step flag requires a number", source)
               end if
            else
               call env%error("The --lvib-hess-step flag requires a number", source)
            end if

         case ('--lvib-step')
            ! Finite difference step of the numerical Hessian (bohr); alias of
            ! `--lvib-hess-step` (CLI redesign sections 6.2 and 7.2)
            call select_runtyp(env, 'lvib', '--lvib-step', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  call set_landscape_hess_step(env, lvset%engine, ddum, &
                     & '--lvib-step', source)
               else
                  call env%error("The --lvib-step flag requires a number", source)
               end if
            else
               call env%error("The --lvib-step flag requires a number", source)
            end if

         case ('--lvib-sccacc')
            ! SCC accuracy of the displaced single points of the Hessian
            call select_runtyp(env, 'lvib', '--lvib-sccacc', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  lvset%engine%hess_acc = ddum
                  lvset%engine%hess_acc_source = '--lvib-sccacc'
               else
                  call env%error("The --lvib-sccacc flag requires a number", source)
               end if
            else
               call env%error("The --lvib-sccacc flag requires a number", source)
            end if

         case ('--lvib-restart')
            ! Wavefunction restart policy of the engine (CLI redesign section
            ! 6.3); the frequency task mirror of `--lts-restart`
            call select_runtyp(env, 'lvib', '--lvib-restart', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               select case (landscape_lowercase(sec))
               case ('off')
                  lvset%engine%restart_mode = p_landscape_restart_off
               case ('pin')
                  lvset%engine%restart_mode = p_landscape_restart_pin
               case ('warm')
                  lvset%engine%restart_mode = p_landscape_restart_warm
               case default
                  call env%error("Unknown landscape restart policy '"//sec//"'", source)
               end select
            else
               call env%error("The --lvib-restart flag requires off, pin or warm", source)
            end if

         case ('--lvib-odlr')
            ! Request the ODLR approximation; no provider implements it, so
            ! the task reports the unsupported request instead of silently
            ! returning a full dense Hessian
            call select_runtyp(env, 'lvib', '--lvib-odlr', selected_runtyp, selected_runtyp_flag)
            lvset%engine%odlr = .true.

         case ('--lvib-no-project')
            ! Keep the translations and rotations in the Hessian
            call select_runtyp(env, 'lvib', '--lvib-no-project', selected_runtyp, selected_runtyp_flag)
            lvset%project_tr = .false.

         case ('--lvib-temp')
            ! Temperature of the thermodynamic table (kelvin)
            call select_runtyp(env, 'lvib', '--lvib-temp', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  lvset%temperature = 0.0_wp
                  lvset%temperature(1) = ddum
                  lvset%n_temperature = 1
                  lvset%temperature_was_set = .true.
               else
                  call env%error("The --lvib-temp flag requires a number", source)
               end if
            else
               call env%error("The --lvib-temp flag requires a number", source)
            end if

         case ('--lvib-rotor-cutoff')
            ! Rotor cutoff of the quasi-RRHO entropy treatment (cm^-1)
            call select_runtyp(env, 'lvib', '--lvib-rotor-cutoff', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  lvset%rotor_cutoff = ddum
                  lvset%rotor_cutoff_was_set = .true.
               else
                  call env%error("The --lvib-rotor-cutoff flag requires a number", source)
               end if
            else
               call env%error("The --lvib-rotor-cutoff flag requires a number", source)
            end if

         case ('--lvib-imag-cutoff')
            ! Imaginary frequency cutoff (cm^-1); needs the minus sign, so the
            ! parser has to be told that the value belongs to the flag
            call select_runtyp(env, 'lvib', '--lvib-imag-cutoff', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  lvset%imag_cutoff = ddum
                  lvset%imag_cutoff_was_set = .true.
               else
                  call env%error("The --lvib-imag-cutoff flag requires a number", source)
               end if
            else
               call env%error("The --lvib-imag-cutoff flag requires a number", source)
            end if

         case ('--lvib-freq-scal')
            ! Frequency scaling factor of the thermochemistry
            call select_runtyp(env, 'lvib', '--lvib-freq-scal', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               read(sec, *, iostat=idum) ddum
               if (idum == 0) then
                  lvset%freq_scal = ddum
                  lvset%freq_scal_was_set = .true.
               else
                  call env%error("The --lvib-freq-scal flag requires a number", source)
               end if
            else
               call env%error("The --lvib-freq-scal flag requires a number", source)
            end if

         case ('--lopt-hessian-file')
            ! File read by `--lopt-hessian file`
            call select_runtyp(env, 'lopt', '--lopt-hessian-file', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               lset%hessian_file = sec
            else
               call env%error("The --lopt-hessian-file flag requires a file name", source)
            end if

         case ('--metadyn')
            call select_runtyp(env, 'md', '--metadyn', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_metadyn(env, 'save', sec)
            end if
            call set_metadyn(env, 'static', 'false')

         case ('--metaopt')
            call select_runtyp(env, 'metaopt', '--metaopt', selected_runtyp, selected_runtyp_flag)
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_opt(env, 'optlevel', sec)
               ! Same tail report as `--bhess`
               landscape_native_level_name = sec
               landscape_native_level_flag = '--metaopt'
            end if

         case ('--nat')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_natom(env, sec)
            end if

         case ('--bias-input', '--gesc')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_metadyn(env, 'bias-input', sec)
            else
               call env%error("No input file for RMSD bias provided", source)
            end if

         case ('--wrtopo')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call setWRtopo(sec, printTopo)
               if (printTopo%warning) call env%error("A wrtopo argument has been misspelled.", source)
            else
               call env%error("The wrtopo keyword is missing an argument.", source)
            end if

         case ('--nocellopt')
            set%optcell = .false.

         end select
         call args%nextFlag(flag)
      end do

      ! Unknown flags.  They stay a warning, or become an error when `--strict`
      ! was given - the option used to be accepted silently, because the
      ! message is an `env%warning` and `env%strict` is only set after the
      ! parsing loop (`xtbMain`).  The report happens here and not inside the
      ! loop so that the decision does not depend on the position of `--strict`
      ! on the command line (CLI redesign section 17.5 5a).
      if (unknown_count > 0) then
         if (unknown_count == 1) then
            message = "Unknown option '"//unknown_flags//"' provided"
         else
            message = "Unknown options '"//unknown_flags//"' provided"
         end if
         if (strict) then
            call env%error(message, source)
         else
            call env%warning(message, source)
         end if
      end if

      ! The value of `-P`/`--parallel`.  A value that is not a positive integer
      ! used to run silently on one thread (`0`, `-2`, clamped by the OpenMP
      ! runtime) or on the default team (`abc`), and `-P` as the last argument
      ! was dropped without a word.  It is reported here for the same reason as
      ! the unknown option above: the decision needs the effect of `--strict`,
      ! which may stand after the flag.  A request above the processor count is
      ! the third silent case; it is only reported, and it stays a warning with
      ! `--strict` as well, because the value is applied as given rather than
      ! clamped - clamping it silently would be the defect class being fixed.
      if (allocated(landscape_parallel_value)) then
         if (len(landscape_parallel_value) == 0) then
            message = "-P/--parallel needs a positive integer, but no value "//&
               & "was given"
         else
            message = "-P/--parallel needs a positive integer, not '"//&
               & landscape_parallel_value//"'"
         end if
         if (strict) then
            call env%error(message, source)
         else
            call env%warning(message, source)
         end if
      else if (nproc > 0 .and. landscape_parallel_team > nproc) then
         write(parallel_note, '(a,i0,a,i0,a)') "-P/--parallel ", &
            & landscape_parallel_team, " exceeds the ", nproc, &
            & " processors of this machine"
         call env%warning(trim(parallel_note), source)
      end if

      ! `--lopt-recalc N` rebuilds the exact Hessian every N RFO steps, `0`
      ! disables the rebuild.  A negative value would not be a second spelling
      ! of "never": the step loop tests `recalc >= 1`, so it would degrade to 0
      ! without a word.  It is reported here - the decision needs the effect of
      ! `--strict`, which may stand after the flag - and the value is applied as
      ! given (the run rebuilds no Hessian), exactly like the `-P` value above
      ! rather than being silently clamped.
      if (lopt_recalc_given .and. lopt_recalc_value < 0) then
         write(recalc_note, '(a,i0,a)') "--lopt-recalc needs 0 (no rebuild) "//&
            & "or a positive number of steps, not ", lopt_recalc_value, &
            & "; the run rebuilds no Hessian"
         if (strict) then
            call env%error(trim(recalc_note), source)
         else
            call env%warning(trim(recalc_note), source)
         end if
      end if

      ! The value of `-a`/`--acc`.  `read(sec,*,iostat=)` accepts `nan` and
      ! `inf`, and the range check of the flag compares against both limits,
      ! so a non-finite value used to be stored as the SCC accuracy: every
      ! energy of the run became NaN without a diagnostic and the process
      ! still exited 0.  The value is reported here - the decision needs the
      ! effect of `--strict`, which may stand after the flag - and it is not
      ! applied, so the run continues with the default (or the level derived)
      ! accuracy.
      if (allocated(landscape_acc_value)) then
         message = "-a/--acc needs a finite accuracy, not '"//&
            & landscape_acc_value//"'; the value is ignored"
         if (strict) then
            call env%error(message, source)
         else
            call env%warning(message, source)
         end if
      end if

      ! The native level value of `--bhess`/`--metaopt` reaches
      ! `set%optset%optlev` through the shared `optlevel2int`
      ! (`src/set_module.f90`), which maps a keyword it does not know to
      ! `normal` - upstream behaviour kept for the shared parser.  A value the
      ! program does not know must be reported rather than silently replaced by
      ! the default, so the rejected keyword is reported here - the same
      ! lax/strict shape and the same position-independent decision as the
      ! unknown option above.
      if (allocated(landscape_native_level_name)) then
         if (.not. optlevel_keyword_known(landscape_native_level_name)) then
            message = "Unknown "//landscape_native_level_flag//" level '"//&
               & landscape_native_level_name//"': the native levels are "//&
               & "crude, sloppy, loose, lax, normal, tight, verytight and "//&
               & "extreme (or their numeric aliases -3 to 3)"
            if (strict) then
               call env%error(message, source)
            else
               call env%warning(message, source)
            end if
         end if
      end if

      ! The step sizes are given on the command line in angstrom (distance) or
      ! degree (angle/dihedral) and converted once, after every flag has been
      ! parsed: this removes the dependency on the order of `--lscan-coord`
      ! and `--lscan-step`.
      do iscan = 1, 2
         if (.not. scan_step_given(iscan)) cycle
         if (lsset%coord(iscan)%ctype == p_landscape_scan_coord_distance) then
            lsset%coord(iscan)%step = scan_step_user(iscan)*aatoau
         else
            lsset%coord(iscan)%step = scan_step_user(iscan)*pi/180.0_wp
         end if
      end do

      ! The unified method selector `--method` is resolved here, after every
      ! flag has been parsed: the run type is only known when the whole
      ! command line has been seen, because the run type flags select it inside
      ! their own case (CLI redesign section 4.7, step 2).  `--method` never
      ! calls `set_runtyp`, so it does not touch the "first run type wins" lock
      ! of `set_runtyp` - `--method bfgs --lopt` and `--lopt --method bfgs` are
      ! the same run.
      if (allocated(landscape_method_name)) then
         landscape_native_engine = landscape_native_engine_from_name(&
            & landscape_method_name)
         if (landscape_native_engine /= 0) then
            ! Native ANCopt engine family (section 4.8).  The engine is pinned
            ! explicitly, because `src/main/defaults.f90:42-46` would otherwise
            ! fill in `p_engine_lbfgs` for a system with more than 500 atoms
            ! and silently turn `--method ancopt` into L-ANCopt.
            select case (set%runtyp)
            case (p_run_opt, p_run_ohess)
               set%opt_engine = landscape_native_engine
            case default
               call env%error("--method "//landscape_method_name//&
                  & " selects a native optimizer engine, which is only "//&
                  & "available with --opt or --ohess", source)
            end select
         else if (landscape_lowercase(landscape_method_name) == 'auto') then
            ! `auto` is the default of the selected task; the settings objects
            ! already hold it, so there is nothing to do.  It is a no-op rather
            ! than a name lookup, so `--lopt-engine bfgs --method auto` keeps
            ! the explicit `bfgs` instead of reporting a conflict.
            continue
         else
            ! Landscape method name: dispatch on the run type (section 4.7,
            ! step 2).  A name of another task is an error that names the flag
            ! of that task (section 4.4), an unknown name is an error as well.
            select case (set%runtyp)
            case (p_run_lopt, p_run_opt, p_run_ohess)
               ! `p_run_opt`/`p_run_ohess` are here since stage 3 redeployed them
               ! to the suite: a landscape method name on `--opt`/`--ohess` now
               ! reaches the suite settings object the run really reads
               ! (`lset`), so the stage 2 refusal of this branch (R39 of the
               ! redesign) is gone.  `--method` never selects the run type, so
               ! `--opt` keeps its identity and `--opt --lopt` stays a conflict.
               if (landscape_method_task(landscape_method_name) /= 'opt') then
                  call landscape_method_mismatch(env, landscape_method_name, source)
               else
                  call set_landscape_opt_algorithm(env, lset, &
                     & landscape_method_name, '--method', source)
               end if
            case (p_run_lscan)
               if (landscape_method_task(landscape_method_name) /= 'opt') then
                  call landscape_method_mismatch(env, landscape_method_name, source)
               else
                  call set_landscape_scan_algorithm(env, lsset, &
                     & landscape_method_name, '--method', source)
               end if
            case (p_run_lts)
               if (landscape_method_task(landscape_method_name) /= 'ts') then
                  call landscape_method_mismatch(env, landscape_method_name, source)
               else
                  call set_landscape_ts_algorithm(env, ltset, &
                     & landscape_method_name, '--method', source)
               end if
            case (p_run_lirc)
               if (landscape_method_task(landscape_method_name) /= 'irc') then
                  call landscape_method_mismatch(env, landscape_method_name, source)
               else
                  call set_landscape_irc_algorithm(env, lircset, &
                     & landscape_method_name, '--method', source)
               end if
            case (p_run_lvib)
               call env%error("--method "//landscape_method_name//&
                  & " cannot be used with --lvib: the frequency task has no "//&
                  & "method selector (use --lvib-hessian for its Hessian "//&
                  & "backend)", source)
            case default
               call env%error("--method needs a landscape run type: add "//&
                  & "--lopt, --lscan, --lts, --lirc or --opt (or use "//&
                  & "--method ancopt with --opt)", source)
            end select
         end if
      end if

      ! ------------------------------------------------------------------------
      ! The upstream `--opt <LEVEL>` keyword and `--cycles N` are bridged onto
      ! the suite settings of the *redeployed* run types here
      ! (CLI redesign sections 3.5 and 3.10b; stage 3a fix round, findings B/C).
      !
      ! The suite reads neither `set%optset%optlev` nor
      ! `set%optset%maxoptcycle`, so without the bridge both flags are dead
      ! knobs for `--opt`/`--ohess` although they did reach the native ANCopt
      ! these run types used before stage 3.  The bridge sits in the run-type
      ! dispatch of the parse tail, where the run type and
      ! `landscape_native_engine` (resolved just above) are both known.  The
      ! priority rule of section 3.5 ("an explicit --lopt-level wins over
      ! --opt <LEVEL>") needs no code: `--lopt-level` selects the `lopt` run
      ! type and `--opt` the `opt` one, so `select_runtyp` already reports the
      ! two together as a conflict that names both flags.

      ! T2: `--opt <LEVEL>`/`--ohess <LEVEL>` set the level through the same
      ! field `--lopt-level` writes (`set_conv_level`), so
      ! `--opt tight` and `--lopt --lopt-level tight` are the same run.  The
      ! banner names the effective level, so the many-to-one map of the
      ! upstream names (and the `banana` refusal, F3) is never invisible.
      if ((set%runtyp == p_run_opt .or. set%runtyp == p_run_ohess) .and. &
         & landscape_native_engine == 0) then
         landscape_opt_conv_level = 'normal'
         if (allocated(landscape_opt_level_name)) then
            landscape_opt_conv_level = &
               & landscape_opt_level_to_conv_level(landscape_opt_level_name)
            if (len(landscape_opt_conv_level) == 0) then
               call env%error("Unknown --opt level '"// &
                  & landscape_opt_level_name//"': use crude, sloppy, loose, "//&
                  & "lax, normal, tight, verytight, extreme, superloose or "//&
                  & "vtight", source)
               landscape_opt_conv_level = 'normal'
            else
               call lset%set_conv_level(landscape_opt_conv_level)
               ! `--ohess` is one run of two stages: the minimisation and the
               ! frequency analysis of the minimised geometry, in the same
               ! process.  Both stages take the level, so the frequency half
               ! does not keep its own defaults and report the SCC accuracy of
               ! level `normal` while the minimisation ran at the requested
               ! one.  Measured before this write: `--ohess tight` printed
               ! `landscape engine: SCC accuracy 0.200000` for the
               ! minimisation and `1.000000` for the frequency analysis of
               ! the same run, i.e. the level reached only `lset`.  An
               ! explicit `--acc` keeps precedence over either half, because
               ! `xtbMain` overwrites `engine%acc` only when the parser
               ! reported `acc_was_set`.  The step budget of the level is
               ! derived for the minimisation only: the frequency task has no
               ! convergence loop and no budget (`--cycles` is refused for
               ! it, section 3.10b).
               if (set%runtyp == p_run_ohess) then
                  call lvset%set_conv_level(landscape_opt_conv_level)
               end if
            end if
         end if
         write(env%unit, '(a)') "[NOTICE] landscape convergence level: "// &
            & landscape_opt_conv_level
      else if (set%runtyp == p_run_opt .or. set%runtyp == p_run_ohess) then
         ! Native engine (`--opt`/`--ohess --method <engine>`): the level stays
         ! on `set%optset%optlev` and reaches `optlevel2int`, which maps an
         ! unknown keyword to `normal`.  Report it here, where the effect of
         ! `--strict` is known; the suite branch above has its own, wider set
         ! (it also knows `superloose`).
         if (allocated(landscape_opt_level_name)) then
            if (.not. optlevel_keyword_known(landscape_opt_level_name)) then
               message = "Unknown --opt level '"//landscape_opt_level_name//&
                  & "': the native levels are crude, sloppy, loose, lax, "//&
                  & "normal, tight, verytight and extreme (or their "//&
                  & "numeric aliases -3 to 3)"
               if (strict) then
                  call env%error(message, source)
               else
                  call env%warning(message, source)
               end if
            end if
         end if
      end if

      ! T3: `--cycles N` is the step budget of the run type that owns one.  The
      ! native run types, and `--opt`/`--ohess --method <native>`, keep the
      ! native semantics (`set%optset%maxoptcycle`, read in
      ! `src/optimizer.f90`), so they are the `case default`.  The bridge also
      ! records the source of the value, so the level-derived budget of
      ! `level_step_budget` (R5, section 3.5) cannot overwrite it.
      if (cycles_given) then
         select case (effective_landscape_task(set%runtyp, landscape_native_engine))
         case (p_run_lopt)
            call landscape_bridge_step_budget(env, cycles_value, &
               & lopt_max_steps_given, lset%max_steps, '--lopt-max-steps', source, &
               & lset%max_steps_source)
         case (p_run_lscan)
            call landscape_bridge_step_budget(env, cycles_value, &
               & lscan_max_steps_given, lsset%max_steps, '--lscan-max-steps', source, &
               & lsset%max_steps_source)
         case (p_run_lts)
            call landscape_bridge_step_budget(env, cycles_value, &
               & lts_max_steps_given, ltset%max_steps, '--lts-max-steps', source, &
               & ltset%max_steps_source)
         case (p_run_lirc)
            call landscape_bridge_step_budget(env, cycles_value, &
               & lirc_max_steps_given, lircset%max_steps, '--lirc-max-steps', source, &
               & lircset%max_steps_source)
         case (p_run_lvib)
            call env%error("--cycles has no meaning for --lvib: the "//&
               & "frequency task has no step budget", source)
         case default
            continue
         end select
      end if

      ! Some saddle methods imply the climbing image stage of the band
      ! (`cineb` and the three climbing node searches; `--lts-climb` does the
      ! same for `neb`).  The stage is applied **once**, to the value the
      ! flags settled on, so a repeated flag is deterministic:
      ! `--lts-method cineb --lts-method neb` and `--method cineb --method
      ! neb` are the same run (stage 2 review L3).  Applied inside the setter
      ! it was applied per flag and the first `cineb` left `climb` on.  The
      ! stage is only ever switched on, exactly like `--lts-climb`, so an
      ! explicit `--lts-climb` keeps working with any method.
      select case (ltset%method)
      case (p_landscape_ts_cineb, p_landscape_ts_autoneb, &
            & p_landscape_ts_string, p_landscape_ts_gsm)
         ltset%path%climb = .true.
      end select

   end subroutine parseArguments

!> Set the thread count of the linear algebra backend for this process.
!>
!> `--blas-threads` is the only caller, so without the flag xtb touches
!> nothing.  The pool cannot be reached through the environment from here: a
!> pthread OpenBLAS sizes it in an ELF constructor that runs before `main()`
!> (`gotoblas_init` -> `blas_get_cpu_number`, which reads
!> `OPENBLAS_NUM_THREADS`/`OMP_NUM_THREADS` and caches the result), so a
!> `setenv` call from the program is too late; the runtime C entry point is the
!> reliable handle.  Measured on OpenBLAS 0.3.21: `setenv` after startup left
!> `openblas_get_num_threads()` at 8, `openblas_set_num_threads` moved it to 1.
!>
!> That entry point is reached by two routes, in this order:
!>
!> 1. A **link-time** binding through the weak references in
!>    `src/prog/blas_threads.c`.  It is the only route that reaches a statically
!>    linked BLAS: a fully static executable has no `.dynsym` at all, so the
!>    entry points live in `.symtab` only and are invisible to `dlsym` on the
!>    global scope, even though the code is in the image.  The reference is
!>    weak so that the route stays optional - if nothing extracts the archive
!>    members defining the entry points the symbol stays NULL, the bridge in
!>    `src/prog/blas_threads.c` reports 0, no call goes through a NULL pointer
!>    and route 2 takes over.
!> 2. `dlopen(NULL, RTLD_NOW)` + `dlsym`, which only inspects the objects the
!>    loader has already mapped and therefore serves a dynamically linked BLAS
!>    (and oneMKL whatever its linkage).  Loading a BLAS library by name is
!>    deliberately avoided: on a machine whose BLAS is MKL that would pull in a
!>    second, unused OpenBLAS.
!>
!> `openblas_set_num_threads` is tried first, then `mkl_set_num_threads`; a
!> backend that offers neither (the reference BLAS, or a statically linked BLAS
!> whose entry points nothing extracted) is reported with a warning and the
!> run continues.  When the matching read-back entry point is present
!> (`openblas_get_num_threads`, or `mkl_get_max_threads`, the oneMKL function
!> `threadpoolctl` reads MKL's thread count with) the effective value is read
!> back and printed, which is the falsifier of the flag: a setter that never
!> ran cannot print the requested value.  A backend offering the setter but no
!> getter prints `(..., not readable)` rather than inventing a number.
!>
!> **Not verified against a real oneMKL**: no MKL is installed on the
!> development machine, so both MKL names are resolved by symbol name only and
!> the branch is checked against the oneMKL documentation, not against the
!> library (it was run once with stub symbols on a reference-BLAS build,
!> `docs/_simplify-recon/82-blas-threads-fix.md` section 2.1).  The OpenBLAS
!> branch and the no-symbol warning are measured
!> (`docs/_simplify-recon/80-blas-threads-validation.md`).
subroutine set_blas_threads(env, nthr)
   use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_char, c_funptr, &
      & c_null_char, c_null_ptr, c_associated, c_f_procpointer

   !> Name of error producer
   character(len=*), parameter :: source = 'prog_main_set_blas_threads'

   !> Calculation environment, only for the output unit
   type(TEnvironment), intent(inout) :: env

   !> Requested number of BLAS threads (>= 1)
   integer, intent(in) :: nthr

   !> Runtime entry points of the two supported backends.
   abstract interface
      !> `void openblas_set_num_threads(int)` / `void mkl_set_num_threads(int)`
      subroutine c_set_num_threads(n) bind(C)
         import :: c_int
         integer(c_int), value :: n
      end subroutine c_set_num_threads

      !> `int openblas_get_num_threads(void)` / `int mkl_get_max_threads(void)`
      function c_get_num_threads() bind(C) result(n)
         import :: c_int
         integer(c_int) :: n
      end function c_get_num_threads
   end interface

   interface
      !> Link-time bridge to the backend entry points, `xtb_blas_set_threads` in
      !> `src/prog/blas_threads.c`.  Returns 0 when this backend has no
      !> link-time binding, which is the cue to use the runtime lookup below.
      function c_blas_set_threads(backend, n) bind(C, name='xtb_blas_set_threads') result(status)
         import :: c_int
         integer(c_int), value :: backend
         integer(c_int), value :: n
         integer(c_int) :: status
      end function c_blas_set_threads

      !> `xtb_blas_get_threads` of the same file: returns 0 when the getter has
      !> no link-time binding, otherwise the reported value is written to `n`.
      function c_blas_get_threads(backend, n) bind(C, name='xtb_blas_get_threads') result(status)
         import :: c_int
         integer(c_int), value :: backend
         integer(c_int) :: n
         integer(c_int) :: status
      end function c_blas_get_threads

#ifndef _WIN32
      !> `void *dlopen(const char *, int)`; `c_null_ptr` is the NULL the
      !> interface reserves for "the current process" (an empty string would
      !> try to open a file of that name).
      function c_dlopen(path, mode) bind(C, name='dlopen') result(handle)
         import :: c_ptr, c_int
         type(c_ptr), value :: path
         integer(c_int), value :: mode
         type(c_ptr) :: handle
      end function c_dlopen

      !> `void *dlsym(void *, const char *)`, returned as a function pointer
      function c_dlsym(handle, symbol) bind(C, name='dlsym') result(addr)
         import :: c_ptr, c_char, c_funptr
         type(c_ptr), value :: handle
         character(kind=c_char), dimension(*), intent(in) :: symbol
         type(c_funptr) :: addr
      end function c_dlsym
#endif
   end interface

   !> The two backends, numbered as `src/prog/blas_threads.c` numbers them
   integer(c_int), parameter :: backend_openblas = 1
   integer(c_int), parameter :: backend_mkl = 2

   !> `RTLD_NOW` of glibc's `dlfcn.h`
   integer(c_int), parameter :: rtld_now = 2

   type(c_ptr) :: handle
   type(c_funptr) :: addr
   procedure(c_set_num_threads), pointer :: set_threads
   procedure(c_get_num_threads), pointer :: get_threads
   !> Backend that accepted the request, 0 while none has
   integer(c_int) :: backend
   integer(c_int) :: effective
   logical :: readable

   backend = 0_c_int
   effective = -1_c_int
   readable = .false.

   !> Route 1: the link-time binding.  This is the only route that works in a
   !> static image, where `dlsym` cannot see the entry points at all.
   if (c_blas_set_threads(backend_openblas, int(nthr, c_int)) /= 0_c_int) then
      backend = backend_openblas
   else if (c_blas_set_threads(backend_mkl, int(nthr, c_int)) /= 0_c_int) then
      backend = backend_mkl
   end if

   !> Route 2: the runtime symbol lookup, for a dynamically linked BLAS and for
   !> any backend route 1 could not bind.  `dlfcn.h` is POSIX; MinGW-w64 and
   !> MSVC have neither it nor `dlsym`, so the whole route is compiled in only
   !> where the dynamic loader exists (`_WIN32` is supplied by CMake, see the
   !> Windows block of `CMakeLists.txt`).  There route 1 is the only route and
   !> the flag reports the warning below rather than failing to link.
#ifndef _WIN32
   handle = c_null_ptr
   if (backend == 0_c_int) then
      handle = c_dlopen(c_null_ptr, rtld_now)
      if (.not. c_associated(handle)) then
         call env%warning("the current process cannot be opened for symbol lookup; "// &
            & "--blas-threads has no effect", source)
         return
      end if
      ! The handle is deliberately never passed to `dlclose`: it is the current
      ! process image, which must stay mapped, so the lookup handle lives for the
      ! rest of the process on purpose (not a leak to be "fixed").

      addr = c_dlsym(handle, 'openblas_set_num_threads'//c_null_char)
      if (c_associated(addr)) then
         call c_f_procpointer(addr, set_threads)
         call set_threads(int(nthr, c_int))
         backend = backend_openblas
      else
         addr = c_dlsym(handle, 'mkl_set_num_threads'//c_null_char)
         if (c_associated(addr)) then
            call c_f_procpointer(addr, set_threads)
            call set_threads(int(nthr, c_int))
            backend = backend_mkl
         end if
      end if
   end if
#else
   handle = c_null_ptr
#endif

   if (backend == 0_c_int) then
      call env%warning("no OpenBLAS or MKL symbol found (not exported by a "// &
         & "statically linked BLAS); --blas-threads has no effect", source)
      return
   end if

   !> Read the effective value back through the same bridge first, because that
   !> is the getter a static image also carries; the symbol lookup is the
   !> fallback for a getter only the loader exposes.  Without a getter the flag
   !> cannot be falsified, so the run says so instead of inventing a number.
   if (c_blas_get_threads(backend, effective) /= 0_c_int) then
      readable = .true.
#ifndef _WIN32
   else if (c_associated(handle)) then
      if (backend == backend_openblas) then
         addr = c_dlsym(handle, 'openblas_get_num_threads'//c_null_char)
      else
         addr = c_dlsym(handle, 'mkl_get_max_threads'//c_null_char)
      end if
      if (c_associated(addr)) then
         call c_f_procpointer(addr, get_threads)
         effective = get_threads()
         readable = .true.
      end if
#endif
   end if

   if (backend == backend_openblas) then
      if (readable) then
         write(env%unit, '(a,i0,a,i0,a)') 'blas threads: ', nthr, &
            & ' (OpenBLAS, read back ', int(effective), ')'
      else
         write(env%unit, '(a,i0,a)') 'blas threads: ', nthr, ' (OpenBLAS)'
      end if
   else
      if (readable) then
         write(env%unit, '(a,i0,a,i0,a)') 'blas threads: ', nthr, &
            & ' (MKL, read back ', int(effective), ')'
      else
         write(env%unit, '(a,i0,a)') 'blas threads: ', nthr, &
            & ' (MKL, not readable)'
      end if
   end if

end subroutine set_blas_threads

!> Select a run type from one command line flag.
!>
!> `set_runtyp` ("src/set_module.f90") accepts a run type only once and reports
!> a second call into the never printed `persistentEnv`, so conflicting run
!> types used to be dropped without any message ("--lopt --lscan" ran the
!> first one, "--lts --lopt-engine bfgs" ignored the sub-option).  This wrapper
!> keeps `selected_typ` (the run type of the first request, empty before the
!> first one) and `selected_flag` (the flag that asked for it) as parser state
!> and turns a request for a *different* run type into an error that names both
!> flags.  A repeated request of the same run type is a no-op because the task
!> sub-options (`--lopt-engine`, `--lscan-method`, ...) select their own run
!> type.
!>
!> `apply = .false.` only records the choice.  `--raman` selects `hess` through
!> `set_elprop` (`src/set_module.f90`), which calls `set_runtyp` itself; a
!> second call from here would hit the run type lock and report it.  The
!> `--raman` path therefore calls `set_runtyp` exactly once, as it always did,
!> while the conflict detection sees the choice (review finding M1 of
!> `docs/_simplify-recon/20-stage1-review.md`).
!>
!> The state is passed in and out explicitly instead of living in a `save`
!> variable: a hidden lock is exactly the defect this wrapper works around, and
!> it would survive the next call of the parser in the same process.
subroutine select_runtyp(env, typ, flag, selected_typ, selected_flag, apply)

   !> Name of error producer
   character(len=*), parameter :: source = 'prog_main_select_runtyp'

   !> Calculation environment for the conflict report
   type(TEnvironment), intent(inout) :: env

   !> Internal run type name as understood by `set_runtyp`
   character(len=*), intent(in) :: typ

   !> Command line flag that requests this run type
   character(len=*), intent(in) :: flag

   !> Run type of the first request, empty until the first request
   character(len=:), allocatable, intent(inout) :: selected_typ

   !> Flag of the first request
   character(len=:), allocatable, intent(inout) :: selected_flag

   !> Apply the run type through `set_runtyp`; `.false.` only records it
   logical, intent(in), optional :: apply

   logical :: do_apply

   do_apply = .true.
   if (present(apply)) do_apply = apply

   if (.not. allocated(selected_typ)) then
      if (do_apply) call set_runtyp(typ)
      selected_typ = typ
      selected_flag = flag
   else if (selected_typ /= typ) then
      call env%error("run type conflict: "//selected_flag//" already selected the '"// &
         & selected_typ//"' run type, "//flag//" would select '"//typ// &
         & "'; only one run type can run (use --ohess, or run two jobs)", source)
   end if

end subroutine select_runtyp

!> Native ANCopt engine selected by a `--method` name.
!>
!> The value is the `p_engine_*` identifier of the native optimizer
!> (`src/setparam.f90:77-80`) and 0 for a landscape method name.  The alias
!> spellings are the keys of the xcontrol `$opt engine=` value
!> (`src/set_module.f90:1886-1890`) with **two deliberate differences**: the
!> bare `lbfgs` is **not** listed here (it is the landscape minimiser of
!> section 4.1 / R21) and `lancopt` **is** listed although it is not an
!> xcontrol key - `$opt engine=lancopt` is rejected there (it reports
!> `engine 'lancopt' is not implemented` and selects nothing), only
!> `lbfgs`/`l-ancopt` do.  The CLI spelling exists because `lbfgs` is taken by
!> the suite (R9); the native L-ANCopt is reachable through `lancopt` and
!> through `l-ancopt`, which both routes accept.  Keeping the two apart is
!> what makes `--lopt --method lbfgs` the suite minimisation and
!> `--opt --method lancopt` the native engine (CLI redesign section 4.7).
!>
!> The earlier claim that this table reproduces `src/set_module.f90:1888-1891`
!> verbatim (stage 2 implementation report section 1) was **inaccurate**: it
!> named `set_module.f90:1889` as the source of `lancopt`, but that line is the
!> `lbfgs|l-ancopt` entry.  Corrected in
!> `docs/_simplify-recon/28-stage2-fix-report.md` section T6.
!>
!> The name is matched case insensitively.
pure function landscape_native_engine_from_name(name) result(engine)
   !> Method name from the command line
   character(len=*), intent(in) :: name
   !> Native engine identifier, 0 for a landscape method name
   integer :: engine

   select case (landscape_lowercase(name))
   case ('ancopt', 'rf')
      engine = p_engine_rf
   case ('lancopt', 'l-ancopt')
      engine = p_engine_lbfgs
   case ('fire', 'inertial')
      engine = p_engine_inertial
   case ('pbc_lbfgs')
      engine = p_engine_pbc_lbfgs
   case default
      engine = 0
   end select

end function landscape_native_engine_from_name

!> Lower case name of a native optimizer engine identifier.
!>
!> The spelling is one of the `--method` names `landscape_native_engine_from_name`
!> accepts (`ancopt`/`rf`, `lancopt`/`l-ancopt`, `fire`/`inertial`,
!> `pbc_lbfgs`), so the engine the banner announces is an engine the user can
!> ask for.  `ancopt` is the canonical spelling of `p_engine_rf`, `lancopt` of
!> `p_engine_lbfgs`: the bare `lbfgs` belongs to the suite minimiser.
!>
!> An identifier outside `p_engine_*` yields `unknown`; the only caller passes
!> the value `--method` resolved.
pure function landscape_native_engine_name(engine) result(name)
   !> Native engine identifier
   integer, intent(in) :: engine
   !> Lower case engine name
   character(len=:), allocatable :: name

   select case (engine)
   case (p_engine_rf)
      name = 'ancopt'
   case (p_engine_lbfgs)
      name = 'lancopt'
   case (p_engine_inertial)
      name = 'fire'
   case (p_engine_pbc_lbfgs)
      name = 'pbc_lbfgs'
   case default
      name = 'unknown'
   end select

end function landscape_native_engine_name

!> Suite convergence level an upstream `--opt <LEVEL>` keyword maps to.
!>
!> The upstream ANCopt level and the suite level are the same axis with eight
!> respectively five values, so the bridge is many-to-one (CLI redesign section
!> 3.5).  The keyword is matched case insensitively; an **empty** result marks a
!> keyword the suite does not know and the caller reports it instead of falling
!> back to `normal` (upstream `optlevel2int` drops an unknown keyword to
!> `normal`, finding F3 of the stage 3a fix round).  The function lives here and
!> not in `set_module.f90`, so the shared `optlevel2int` keeps its behaviour for
!> the native run types and its xcontrol users.
!>
!> The keywords are the eight upstream names, their numeric aliases (`-3`..`3`)
!> and the suite-only `superloose`; `vtight` is the second spelling of
!> `verytight`.
pure function landscape_opt_level_to_conv_level(name) result(level)
   !> Upstream `--opt` level keyword
   character(len=*), intent(in) :: name
   !> Suite convergence level name, empty for an unknown keyword
   character(len=:), allocatable :: level

   select case (landscape_lowercase(name))
   case ('crude', 'sloppy', 'superloose', '-3', '-2')
      level = 'superloose'
   case ('loose', 'lax', '-1')
      level = 'loose'
   case ('normal', '0')
      level = 'normal'
   case ('tight', '1')
      level = 'tight'
   case ('verytight', 'vtight', 'extreme', '2', '3')
      level = 'vtight'
   case default
      level = ''
   end select

end function landscape_opt_level_to_conv_level

!> Apply the `--cycles N` step budget of one landscape task.
!>
!> `--cycles` is the lowest-priority spelling (CLI redesign section 3.10b): it
!> sets the budget when the task-specific flag (`--lopt-max-steps`, ...) is
!> absent, and a **different** value through both flags is refused instead of
!> one silently overriding the other.  The same value through both is
!> idempotent.
subroutine landscape_bridge_step_budget(env, cycles, flag_given, budget, flag, source, budget_source)
   !> Computational environment for the conflict report
   type(TEnvironment), intent(inout) :: env
   !> Value of `--cycles`
   integer, intent(in) :: cycles
   !> The task-specific step budget flag was given
   logical, intent(in) :: flag_given
   !> Step budget of the task settings, overwritten when `--cycles` wins
   integer, intent(inout) :: budget
   !> Name of the task-specific flag
   character(len=*), intent(in) :: flag
   !> Name of error producer
   character(len=*), intent(in) :: source
   !> Flag that set `budget`, set when `--cycles` wins; it keeps the
   !> level-derived budget of `landscape_settings%level_step_budget` away
   !> from a value the user asked for (CLI redesign section 3.5, R5)
   character(len=:), allocatable, intent(inout), optional :: budget_source
   character(len=160) :: message

   if (.not. flag_given) then
      budget = cycles
      if (present(budget_source)) budget_source = '--cycles'
   else if (budget /= cycles) then
      write(message, '(a,i0,a,a,a,i0,a)') "--cycles ", cycles, &
         & " conflicts with ", flag, " ", budget, &
         & ": give only one step budget"
      call env%error(trim(message), source)
   end if

end subroutine landscape_bridge_step_budget

!> Landscape task executed by a run type, i.e. the answer to "which task does
!> this process run on the landscape side?".
!>
!> The run types stay what they are: a redeployment changes **who** executes a
!> run, not what it is (CLI redesign section 3.11, decision 1).  `--opt` and
!> `--ohess` without a `--method` are the suite's minimisation; with a native
!> engine they keep the native path and are returned unchanged here.  Every
!> other run type is returned unchanged as well, so the function is the single
!> place that knows about the redeployment and the callers can only ask about
!> one task at a time (`effective_landscape_task(...) == p_run_lopt`, ...).
!>
!> `landscape_is_suite_run` is the derived predicate that tells a suite run
!> from a native one; together they are the **only** criterion of the guards,
!> the suppression of the native blocks, the dispatch and `--cycles`/`--ljson`.
!> A second criterion such as `native_engine == p_engine_rf` (the H-2 mistake
!> of the redesign review) would classify native L-ANCopt and FIRE as suite
!> runs and report the three guards against them.
!>
!> `runtyp` is `set%runtyp`, `native_engine` is `landscape_native_engine`
!> (0 = the suite).
pure function effective_landscape_task(runtyp, native_engine) result(task)
   !> Run type of the process
   integer, intent(in) :: runtyp
   !> Native engine selected by `--method`, 0 for the suite
   integer, intent(in) :: native_engine
   !> Landscape task of the run
   integer :: task

   if ((runtyp == p_run_opt .or. runtyp == p_run_ohess) .and. native_engine == 0) then
      ! `--ohess` runs the minimisation first; its frequency analysis is the
      ! second stage of the same task (`landscape_ohess_chain`), not a task of
      ! its own.
      task = p_run_lopt
   else
      task = runtyp
   end if

end function effective_landscape_task

!> Short task name of a run type for the `task` key of the JSON document.
!>
!> The native run types need a stable, lower case name for the
!> `landscape_<task>.json` file and for the `task` key, the way the suite
!> tasks have one (`sp`, `opt`, `scan`, `ts`, `irc`, `vib`).  The mapping is
!> the one the `--help`/man run type table uses; an unknown run type (which
!> the dispatch rejects before this point anyway) yields `run`, so the writer
!> never produces an empty file name.
pure function landscape_runtyp_task_name(runtyp) result(name)
   !> Run type of the process
   integer, intent(in) :: runtyp
   !> Lower case task name
   character(len=:), allocatable :: name

   select case (runtyp)
   case (p_run_scc)
      name = 'sp'
   case (p_run_prescc)
      name = 'prescc'
   case (p_run_grad)
      name = 'grad'
   case (p_run_opt)
      name = 'opt'
   case (p_run_hess)
      name = 'hess'
   case (p_run_ohess)
      name = 'ohess'
   case (p_run_bhess)
      name = 'bhess'
   case (p_run_md)
      name = 'md'
   case (p_run_omd)
      name = 'omd'
   case (p_run_path)
      name = 'path'
   case (p_run_screen)
      name = 'screen'
   case (p_run_modef)
      name = 'modef'
   case (p_run_mdopt)
      name = 'mdopt'
   case (p_run_metaopt)
      name = 'metaopt'
   case (p_run_vip)
      name = 'vip'
   case (p_run_vea)
      name = 'vea'
   case (p_run_vipea)
      name = 'vipea'
   case (p_run_vfukui)
      name = 'vfukui'
   case (p_run_vomega)
      name = 'vomega'
   case default
      name = 'run'
   end select

end function landscape_runtyp_task_name

!> TRUE if the landscape suite executes this run type.
!>
!> This is the single criterion of the three xcontrol guards of
!> `src/prog/main.F90` (the `$scan` block, periodic systems and the `$opt`
!> block), of the `--ljson` message and of the suppressions of the native
!> geometry writer and the timing label.
pure function landscape_is_suite_run(runtyp, native_engine) result(yes)
   !> Run type of the process
   integer, intent(in) :: runtyp
   !> Native engine selected by `--method`, 0 for the suite
   integer, intent(in) :: native_engine
   !> TRUE if a landscape task runs
   logical :: yes

   select case (effective_landscape_task(runtyp, native_engine))
   case (p_run_lopt, p_run_lscan, p_run_lts, p_run_lirc, p_run_lvib)
      yes = .true.
   case default
      yes = .false.
   end select

end function landscape_is_suite_run

!> TRUE for the redeployed `--ohess`, the minimisation **followed by** the
!> frequency analysis.
!>
!> The native `--ohess` makes the same two calls (`main_geometry` +
!> `main_freq`); the suite chain is the counterpart of them for the redeployed
!> run type, and the condition of the `fres` bridge that keeps the console
!> output of `--ohess` unchanged.
pure function landscape_ohess_chain(runtyp, native_engine) result(yes)
   !> Run type of the process
   integer, intent(in) :: runtyp
   !> Native engine selected by `--method`, 0 for the suite
   integer, intent(in) :: native_engine
   !> TRUE if the suite runs the minimisation and the frequency analysis
   logical :: yes

   yes = runtyp == p_run_ohess .and. native_engine == 0

end function landscape_ohess_chain

!> Announce the engine of the run in one framed box, before any geometry is
!> evaluated.
!>
!> The engine used to be echoed by each task's print block, which runs after
!> the always-done single point: the run was already evaluating the surface
!> under an engine the user had not seen named, and two of the four echoes
!> (`--lopt-engine`, `--lirc-method`) only existed as `[NOTICE]` lines for the
!> **builtin default**, so an explicitly named engine was silent.  This is the
!> single announcement instead: one box, printed for every run type that has an
!> engine, whether the engine is the builtin default or was named on the
!> command line - the trailing `(default)`/`(explicit)` says which, so the two
!> cases are never confused.
!>
!> The label names the choice the run really makes: the minimisation engine of
!> `--lopt`/`--opt`/`--ohess` and of the relaxed scan's inner loop, the saddle
!> search of `--lts`, the integrator of `--lirc`, and the native `--method`
!> engine of `--opt`/`--ohess`.  A run type with no engine prints nothing: a
!> single point, a Hessian, MD (`--sp`, `--hess`, `--md`, ...), `--lvib` (it
!> has a Hessian backend, not an engine) and a rigid scan (it constructs the
!> geometry instead of minimising it).
!>
!> The box goes to `iunit`, which is the standard output of the run: the two
!> existing wordings are forbidden in a native run (no `landscape: ` in
!> `xtb/landscape-method-native`, no `[NOTICE]` in
!> `xtb/landscape-opt-native-exit`), and `generic_header` contains neither.
subroutine landscape_engine_banner(iunit, runtyp, native_engine, &
      & optset, scanset, tsset, ircset)
   !> I/O unit
   integer, intent(in) :: iunit
   !> Run type of the process
   integer, intent(in) :: runtyp
   !> Native ANCopt engine selected by `--method`, 0 when the suite runs
   integer, intent(in) :: native_engine
   !> Settings of the minimisation task (the redeployed `--opt`/`--ohess` too)
   type(landscape_opt_settings), intent(in) :: optset
   !> Settings of the scan task
   type(landscape_scan_settings), intent(in) :: scanset
   !> Settings of the saddle search task
   type(landscape_ts_settings), intent(in) :: tsset
   !> Settings of the reaction path task
   type(landscape_irc_settings), intent(in) :: ircset

   !> Box label and engine name of the line, and whether the engine was named
   !> on the command line or is the builtin default
   character(len=:), allocatable :: label, name, origin

   if (native_engine /= 0) then
      label = "Engine"
      name = landscape_native_engine_name(native_engine)
      origin = "explicit"
   else
      select case (runtyp)
      case (p_run_lopt, p_run_opt, p_run_ohess)
         label = "Engine"
         name = landscape_opt_algorithm_name(optset%algorithm)
         if (optset%algorithm_was_set) then
            origin = "explicit"
         else
            origin = "default"
         end if
      case (p_run_lscan)
         if (landscape_lowercase(scanset%mode) /= "relaxed") return
         label = "Engine"
         name = trim(scanset%method)
         if (allocated(scanset%method_source)) then
            origin = "explicit"
         else
            origin = "default"
         end if
      case (p_run_lts)
         label = "Saddle search"
         name = landscape_ts_method_name(tsset%method)
         if (allocated(tsset%method_source)) then
            origin = "explicit"
         else
            origin = "default"
         end if
      case (p_run_lirc)
         label = "IRC integrator"
         name = landscape_irc_method_name(ircset%method)
         if (allocated(ircset%method_source)) then
            origin = "explicit"
         else
            origin = "default"
         end if
      case default
         return
      end select
   end if

   write(iunit, '(a)')
   call generic_header(iunit, label//": "//name//" ("//origin//")", 49, 10)

end subroutine landscape_engine_banner

!> Report a positional command line argument the run cannot use.
!>
!> Two shapes: an additional **existing** input file (a run reads one input
!> file, the last positional one) and a token that is neither an option nor an
!> existing file (`normal`, `nope.xyz`).  Both were warnings until the owner
!> ruling recorded in `docs/_simplify-recon/197`, and the warning was not
!> harmless: `oxtb water.xyz --lvib normal` - a level keyword passed as a
!> positional argument - exited 0 and computed the frequencies of the
!> *unrelaxed* input geometry, because the token was dropped and the run went
!> ahead.  The caller raises the errors in the loop that walks the arguments
!> and aborts with `env%checkpoint` right after it, before anything is read.
!>
!> The message names the token, and a token that is a convergence level of the
!> run type also names the flag that takes it: `normal` is a legal value of
!> `--lopt-level`, and a user who passed it wants a level, not a file.
subroutine positional_argument_error(env, token, runtyp, duplicate, source)
   !> Computational environment for the error report
   type(TEnvironment), intent(inout) :: env
   !> The positional argument
   character(len=*), intent(in) :: token
   !> Run type of the process, for the level-keyword hint
   integer, intent(in) :: runtyp
   !> TRUE for an additional existing input file, FALSE for an unusable token
   logical, intent(in) :: duplicate
   !> Name of error producer
   character(len=*), intent(in) :: source

   !> Error text, and the level flag of the run type (empty when it has none)
   character(len=:), allocatable :: message, flag

   if (duplicate) then
      message = "Input file '"//token//"' is not used: a run reads one "// &
         & "input file, the last positional file of the command line"
   else
      message = "Input file '"//token//"' is neither an option nor an "// &
         & "existing file, so the run cannot use it"
      flag = landscape_level_flag(runtyp)
      if (len(flag) > 0 .and. landscape_conv_level_known(token)) then
         message = message//" (did you mean "//flag//" "// &
            & landscape_lowercase(token)//"?)"
      end if
   end if
   call env%error(message, source)

end subroutine positional_argument_error

!> Flag that takes the convergence level of a suite run type.
!>
!> Empty when the run type has no level keyword, which is why the hint of
!> `positional_argument_error` stays silent for them: `--lvib` has no level
!> flag at all, and `--opt`/`--ohess` take their level as **their own**
!> argument (`--opt tight`), so naming a separate flag for them would point at
!> a different run type.  Every other suite run type writes its level through
!> the flag returned here (`src/prog/main.F90`, the `--l*-level` branches).
pure function landscape_level_flag(runtyp) result(flag)
   !> Run type of the process
   integer, intent(in) :: runtyp
   !> Spelling of the level flag, empty when the run type has none
   character(len=:), allocatable :: flag

   select case (runtyp)
   case (p_run_lopt)
      flag = '--lopt-level'
   case (p_run_lscan)
      flag = '--lscan-level'
   case (p_run_lts)
      flag = '--lts-level'
   case (p_run_lirc)
      flag = '--lirc-level'
   case default
      flag = ''
   end select

end function landscape_level_flag

!> TRUE if the **native** numerical Hessian and frequency path of this run
!> executes.
!>
!> `--hess` and `--bhess` are native, and `--ohess` is native in its redeployed
!> half only (`--ohess --method ancopt|lancopt|fire|pbc_lbfgs`); the redeployed
!> `--ohess` gets its Hessian from the suite frequency chain instead.  Without
!> this distinction the same run would produce the Hessian twice and write
!> `hessian`/`vibspectrum`/`g98.out` twice (CLI redesign section 3.11.1, T1/T2).
!> Every branch assigns, including `case default`.
pure function landscape_native_hessian_path(runtyp, native_engine) result(yes)
   !> Run type of the process
   integer, intent(in) :: runtyp
   !> Native engine selected by `--method`, 0 for the suite
   integer, intent(in) :: native_engine
   !> TRUE if the native Hessian/frequency block runs
   logical :: yes

   select case (runtyp)
   case (p_run_hess, p_run_bhess)
      yes = .true.
   case (p_run_ohess)
      yes = native_engine /= 0
   case default
      yes = .false.
   end select

end function landscape_native_hessian_path

!> Task a landscape method name belongs to.
!>
!> Returns `'opt'` for the six minimisers (used by `--lopt` and `--lscan`),
!> `'ts'` for the eight saddle searches and `'irc'` for the four reaction path
!> integrators; an empty string marks an unknown name.  The answer is derived
!> from the three name lookups of `xtb_landscape_settings`, so the set of
!> valid names has exactly one spelling.  Every name belongs to exactly one
!> task (section 4.2: `gs` is the IRC integrator, not a saddle search; `gsm`
!> is the string search, not the IRC method).
function landscape_method_task(name) result(task)
   !> Method name from the command line
   character(len=*), intent(in) :: name
   !> Task of the name, empty for an unknown name
   character(len=:), allocatable :: task

   if (landscape_opt_algorithm_from_name(name) /= 0) then
      task = 'opt'
   else if (landscape_ts_method_from_name(name) /= 0) then
      task = 'ts'
   else if (landscape_irc_method_from_name(name) /= 0) then
      task = 'irc'
   else
      task = ''
   end if

end function landscape_method_task

!> Report a `--method` name that does not belong to the selected run type.
!>
!> The message names the run type flag of the task the name belongs to, so a
!> mistyped `--method dimer` on `--lopt` points at `--lts` (CLI redesign
!> section 4.4).  An unknown name and a `--method` without a run type are
!> reported here as well, because both are the same kind of user error.
subroutine landscape_method_mismatch(env, name, source)
   !> Computational environment for the error report
   type(TEnvironment), intent(inout) :: env
   !> Method name from the command line
   character(len=*), intent(in) :: name
   !> Name of error producer
   character(len=*), intent(in) :: source

   select case (landscape_method_task(name))
   case ('opt')
      call env%error(name//" is a minimisation method; use it with --lopt "//&
         & "or --lscan (the landscape suite was selected with --method)", source)
   case ('ts')
      call env%error(name//" is a transition-state method; use it with --lts "//&
         & "(or --lts-method "//name//")", source)
   case ('irc')
      call env%error(name//" is an intrinsic reaction coordinate method; use "//&
         & "it with --lirc", source)
   case default
      call env%error("unknown landscape method '"//trim(name)//"': ask for "//&
         & "--method auto for the task default, or name a minimiser "//&
         & "(sd, cg, lbfgs, bfgs, diis, rfo), a saddle search (prfo, bprfo, "//&
         & "dimer, neb, cineb, autoneb, string, gsm) or an IRC integrator "//&
         & "(gs, hpc, lqa, eulerpc)", source)
   end select

end subroutine landscape_method_mismatch

!> Select the landscape minimisation algorithm of the `opt` task.
!>
!> Shared setter of `--method` and `--lopt-engine` (CLI redesign section 4.3):
!> it resolves the name, records which flag selected it and refuses a second,
!> **different** selection through **another** flag.  The same value is
!> idempotent and a repeated flag keeps the last value, so the two aliases
!> compose in any order.  The run type is **not** touched here - the setter is
!> pure in that respect, which is what lets the parse-tail normalisation call
!> it for `--lopt` without disturbing the "first run type wins" lock.
subroutine set_landscape_opt_algorithm(env, settings, name, owner, source)
   !> Computational environment for the error report
   type(TEnvironment), intent(inout) :: env
   !> Landscape minimisation settings
   type(landscape_opt_settings), intent(inout) :: settings
   !> Algorithm name from the command line
   character(len=*), intent(in) :: name
   !> Flag that selected the name (used by the conflict message)
   character(len=*), intent(in) :: owner
   !> Name of error producer
   character(len=*), intent(in) :: source

   integer :: algorithm

   algorithm = landscape_opt_algorithm_from_name(name)
   if (algorithm == 0) then
      call env%error("Unknown landscape minimisation algorithm '"//name//"'", source)
      return
   end if
   if (allocated(settings%algorithm_source) .and. &
         & settings%algorithm_source /= owner .and. settings%algorithm /= algorithm) then
      call env%error("conflicting landscape methods: "//owner//" "//&
         & landscape_lowercase(name)//" and "//settings%algorithm_source//" "//&
         & landscape_opt_algorithm_name(settings%algorithm), source)
      return
   end if
   settings%algorithm = algorithm
   settings%algorithm_source = owner
   settings%algorithm_was_set = .true.

end subroutine set_landscape_opt_algorithm

!> Select the inner minimiser of the `scan` task.
!>
!> Shared setter of `--method` and `--lscan-method` (CLI redesign section 4.3,
!> same conflict rule as `set_landscape_opt_algorithm`).  The scan stores the
!> name, not the identifier, so it is stored lower case - the task resolves it
!> again with `landscape_opt_algorithm_from_name`.
subroutine set_landscape_scan_algorithm(env, settings, name, owner, source)
   !> Computational environment for the error report
   type(TEnvironment), intent(inout) :: env
   !> Landscape scan settings
   type(landscape_scan_settings), intent(inout) :: settings
   !> Algorithm name from the command line
   character(len=*), intent(in) :: name
   !> Flag that selected the name (used by the conflict message)
   character(len=*), intent(in) :: owner
   !> Name of error producer
   character(len=*), intent(in) :: source

   integer :: algorithm, previous

   algorithm = landscape_opt_algorithm_from_name(name)
   if (algorithm == 0) then
      call env%error("Unknown landscape scan minimiser '"//name//"'", source)
      return
   end if
   if (allocated(settings%method_source)) then
      previous = landscape_opt_algorithm_from_name(settings%method)
      if (settings%method_source /= owner .and. previous /= algorithm) then
         call env%error("conflicting landscape methods: "//owner//" "//&
            & landscape_lowercase(name)//" and "//settings%method_source//" "//&
            & trim(settings%method), source)
         return
      end if
   end if
   settings%method = landscape_lowercase(name)
   settings%method_source = owner

end subroutine set_landscape_scan_algorithm

!> Select the saddle search method of the `ts` task.
!>
!> Shared setter of `--method` and `--lts-method` (CLI redesign section 4.3,
!> same conflict rule).  The setter stores the method **only**; the climbing
!> image stage that `cineb`, `autoneb`, `string` and `gsm` imply is applied
!> once to the final value in `parseArguments`, after the flag loop (see the
!> comment there): a side effect inside the setter would be applied per flag,
!> so `--lts-method cineb --lts-method neb` kept `climb = .true.` while the
!> equivalent `--method cineb --method neb` did not (stage 2 review L3).
subroutine set_landscape_ts_algorithm(env, settings, name, owner, source)
   !> Computational environment for the error report
   type(TEnvironment), intent(inout) :: env
   !> Landscape transition state settings
   type(landscape_ts_settings), intent(inout) :: settings
   !> Method name from the command line
   character(len=*), intent(in) :: name
   !> Flag that selected the name (used by the conflict message)
   character(len=*), intent(in) :: owner
   !> Name of error producer
   character(len=*), intent(in) :: source

   integer :: method

   method = landscape_ts_method_from_name(name)
   if (method == 0) then
      call env%error("Unknown landscape transition state method '"//name//"'", source)
      return
   end if
   if (allocated(settings%method_source) .and. &
         & settings%method_source /= owner .and. settings%method /= method) then
      call env%error("conflicting landscape methods: "//owner//" "//&
         & landscape_lowercase(name)//" and "//settings%method_source//" "//&
         & landscape_ts_method_name(settings%method), source)
      return
   end if
   settings%method = method
   settings%method_source = owner

end subroutine set_landscape_ts_algorithm

!> Select the IRC integration method of the `lirc` task.
!>
!> Shared setter of `--method` and `--lirc-method` (CLI redesign section 4.3,
!> same conflict rule).
subroutine set_landscape_irc_algorithm(env, settings, name, owner, source)
   !> Computational environment for the error report
   type(TEnvironment), intent(inout) :: env
   !> Landscape reaction path settings
   type(landscape_irc_settings), intent(inout) :: settings
   !> Method name from the command line
   character(len=*), intent(in) :: name
   !> Flag that selected the name (used by the conflict message)
   character(len=*), intent(in) :: owner
   !> Name of error producer
   character(len=*), intent(in) :: source

   integer :: method

   method = landscape_irc_method_from_name(name)
   if (method == 0) then
      call env%error("Unknown landscape IRC method '"//name//"'", source)
      return
   end if
   if (allocated(settings%method_source) .and. &
         & settings%method_source /= owner .and. settings%method /= method) then
      call env%error("conflicting landscape methods: "//owner//" "//&
         & landscape_lowercase(name)//" and "//settings%method_source//" "//&
         & landscape_irc_method_name(settings%method), source)
      return
   end if
   settings%method = method
   settings%method_source = owner

end subroutine set_landscape_irc_algorithm

!> Set the finite difference step of the numerical Hessian of a task (bohr).
!>
!> Shared setter of the two spellings of every task that owns a Hessian step
!> (CLI redesign sections 6.3 and 7.2): `--lts-step` and `--lts-hess-step`,
!> `--lirc-step-hess` and `--lirc-hess-step`, `--lvib-step` and
!> `--lvib-hess-step`; `--lopt-hess-step` and `--lscan-hess-step` are the only
!> spelling of their task.  The alias rule is the one of the method selectors
!> (section 4.3): the same value through either spelling is equivalent, a
!> different value through **another** flag is an error that names the second
!> flag, and a repeated flag keeps the last value.  The value is compared
!> numerically, so `--lts-step 0.01 --lts-hess-step 1.0e-2` is the same
!> request and not a conflict.
subroutine set_landscape_hess_step(env, ctrl, step, owner, source)
   !> Computational environment for the error report
   type(TEnvironment), intent(inout) :: env
   !> Engine control that carries the finite difference step
   type(landscape_engine_control), intent(inout) :: ctrl
   !> Finite difference step from the command line (bohr)
   real(wp), intent(in) :: step
   !> Flag that set the step (used by the conflict message)
   character(len=*), intent(in) :: owner
   !> Name of error producer
   character(len=*), intent(in) :: source

   character(len=32) :: requested, previous

   if (allocated(ctrl%hess_step_source)) then
      if (ctrl%hess_step_source /= owner .and. ctrl%hess_step /= step) then
         write(requested, '(es14.6)') step
         write(previous, '(es14.6)') ctrl%hess_step
         call env%error("conflicting Hessian finite difference steps: "//&
            & owner//" "//trim(adjustl(requested))//" bohr and "//&
            & ctrl%hess_step_source//" "//trim(adjustl(previous))//" bohr", &
            & source)
         return
      end if
   end if
   ctrl%hess_step = step
   ctrl%hess_step_source = owner

end subroutine set_landscape_hess_step

!> Read up to four 1-based atom indices from a comma or blank separated list.
!>
!> Used by the `--lscan-atoms` and `--lscan-atoms2` flags.  The remaining
!> entries of `atoms` are zero; the scan task rejects a coordinate whose atom
!> list is incomplete (`xtb_landscape_scan_relax::landscape_scan_check`).
   subroutine read_scan_atoms(env, spec, atoms, source)
      !> Computational environment for the error report
      type(TEnvironment), intent(inout) :: env
      !> Atom index list
      character(len=*), intent(in) :: spec
      !> Parsed indices, unused entries are zero
      integer, intent(out) :: atoms(4)
      !> Source of the error report
      character(len=*), intent(in) :: source

      character(len=32) :: token
      integer :: nfound, pos, stop, ierr, value, last

      atoms = 0
      nfound = 0
      last = len_trim(spec)
      pos = 1
      do while (pos <= last .and. nfound < size(atoms))
         do while (pos <= last)
            if (spec(pos:pos) == ',' .or. spec(pos:pos) == ' ') then
               pos = pos + 1
            else
               exit
            end if
         end do
         if (pos > last) exit
         stop = pos
         do while (stop <= last)
            if (spec(stop:stop) == ',' .or. spec(stop:stop) == ' ') exit
            stop = stop + 1
         end do
         token = spec(pos:stop - 1)
         read(token, *, iostat=ierr) value
         if (ierr /= 0) then
            call env%error("cannot read the atom index '"//trim(token)//"'", source)
            return
         end if
         nfound = nfound + 1
         atoms(nfound) = value
         pos = stop
      end do

   end subroutine read_scan_atoms

!> kronecker delta
   function kron(i, j) result(res_kronij)
      integer, intent(in) :: i, j
      real(wp) :: res_kronij

      res_kronij = 0.0_wp
      if (i == j) res_kronij = 1.0_wp

   end function kron
   function read_whole_file(fname) result(list)
      character(len=*), intent(in) :: fname
      character(len=:), allocatable :: list
      integer :: io, stat
      character(len=:), allocatable :: line
      open (newunit=io, file=fname, iostat=stat)
      call getline(io, list, stat)
      do while (stat == 0)
         call getline(io, line, stat)
         if (stat == 0) list = list//","//line
      end do
      close (io, iostat=stat)
   end function read_whole_file

! set booleans for requested topology list printout
   subroutine setWRtopo(sec, printTopo)
      ! command line argument
      character(len=*), intent(in) :: sec
      ! type holds booleans of to be printed topology lists
      type(TPrintTopo), intent(inout) :: printTopo
      ! seperator for lists is ","
      character, parameter :: sep = ","
      ! current and old position of seperator
      integer :: curr_pos, old_pos
      integer :: lenSec, i

      curr_pos = 0
      old_pos = 0
      lenSec = len(sec)
      do i = 1, lenSec
         curr_pos = scan(sec(curr_pos + 1:lenSec), sep) + old_pos
         if (curr_pos /= old_pos) then
            call selectList(sec(old_pos + 1:curr_pos - 1), printTopo)
         else
            call selectList(sec(old_pos + 1:lenSec), printTopo)
            exit
         end if
         old_pos = curr_pos
      end do

   end subroutine setWRtopo

   subroutine selectList(secSplit, printTopo)
      ! part of command line argument
      character(len=*), intent(in) :: secSplit
      ! holds booleans of to be printed topology lists
      type(TPrintTopo), intent(inout) :: printTopo

      select case (secSplit)
      case ("etot")
         printTopo%etot = .true.
      case ("gnorm")
         printTopo%gnorm = .true.
      case ("nb")
         printTopo%nb = .true.
      case ("bpair")
         printTopo%bpair = .true.
      case ("alist")
         printTopo%alist = .true.
      case ("blist")
         printTopo%blist = .true.
      case ("tlist")
         printTopo%tlist = .true.
      case ("vtors")
         printTopo%vtors = .true.
      case ("vbond")
         printTopo%vbond = .true.
      case ("vangl")
         printTopo%vangl = .true.
      case ("hbbond")
         printTopo%hbbond = .true.
      case ("eeq")
         printTopo%eeq = .true.
      case ("json")
         printTopo%etot = .true.
         printTopo%gnorm = .true.
         printTopo%nb = .true.
         printTopo%bpair = .true.
         printTopo%alist = .true.
         printTopo%blist = .true.
         printTopo%tlist = .true.
         printTopo%vtors = .true.
         printTopo%vbond = .true.
         printTopo%vangl = .true.
         printTopo%hbbond = .true.
         printTopo%eeq = .true.
      case default
         printTopo%warning = .true.
      end select
   end subroutine selectList

   subroutine ptb_feature_not_implemented(env)
      !> Computational environment
      type(TEnvironment), intent(inout) :: env

      call env%error("PTB not available without 'tblite'. Compiled without support for 'tblite' library.")
      call env%error("Rebuild with the CMake option -DWITH_TBLITE=ON to use it.")
   end subroutine ptb_feature_not_implemented

   !>  make some post processing afterward, show some timings and stuff
   !>
   !> `native_engine` is the `--method` selector of the run (`xtbMain`); the
   !> timing label of slot 3 depends on whether the landscape suite or the
   !> native ANCopt family ran, and `landscape_native_engine` is a local of
   !> `xtbMain`, so it has to be passed in (CLI redesign read point 14).
   subroutine finalize_xtb(env, native_engine)
    
      !> Calculation environment
      type(TEnvironment), intent(in) :: env

      !> Native engine selected by `--method`, 0 for the landscape suite
      integer, intent(in) :: native_engine

      write (env%unit, '(a)')
      write (env%unit, '(72("-"))')
      call stop_timing_run
      call stop_timing(1)
      call prdate('E')
      write (env%unit, '(72("-"))')
      call prtiming(1, 'total')
      if (.not. set%runtyp == p_run_prescc) &
         call prtiming(2, 'SCF')
      ! The label names the timer that actually ran: the redeployed
      ! `--opt`/`--ohess` go through the landscape optimiser, so their timer slot
      ! 3 must not be announced as `ANC optimizer` (CLI redesign read point 14).
      if (((set%runtyp == p_run_opt) .or. (set%runtyp == p_run_ohess) .or. &
         &   (set%runtyp == p_run_omd) .or. (set%runtyp == p_run_metaopt)) .and. &
         &   .not. landscape_is_suite_run(set%runtyp, native_engine)) then
         call prtiming(3, 'ANC optimizer')
      end if
      if (effective_landscape_task(set%runtyp, native_engine) == p_run_lopt) then
         call prtiming(3, 'landscape optimizer')
      end if
      if (set%runtyp == p_run_lscan) then
         call prtiming(3, 'landscape scan')
      end if
      if (set%runtyp == p_run_lvib) then
         call prtiming(3, 'landscape frequencies')
      end if
      if (set%runtyp == p_run_lirc) then
         call prtiming(3, 'landscape IRC')
      end if
      if (set%runtyp == p_run_path) then
         call prtiming(4, 'path finder')
      end if
      ! The same criterion as the suppression of the native Hessian block
      ! (`landscape_native_hessian_path`): the redeployed `--ohess` takes its
      ! Hessian from the suite frequency chain, so slot 5 never runs for it and
      ! the label must not claim an `analytical hessian` that did not run (the
      ! raw run-type test left the label behind with a bogus `0.000 sec` - fix
      ! round T1 of `docs/_simplify-recon/45-stage3a-fix.md`).
      if (landscape_native_hessian_path(set%runtyp, native_engine)) then
         if (set%mode_extrun /= p_ext_turbomole) then
            call prtiming(5, 'analytical hessian')
         else
            call prtiming(5, 'numerical hessian')
         end if
      end if
      if ((set%runtyp == p_run_md) .or. (set%runtyp == p_run_omd) .or. &
          (set%runtyp == p_run_metaopt)) then
         call prtiming(6, 'MD')
      end if
      if (set%runtyp == p_run_screen) then
         call prtiming(8, 'screen')
      end if
      if (set%runtyp == p_run_modef) then
         call prtiming(9, 'mode following')
      end if
      if (set%runtyp == p_run_mdopt) then
         call prtiming(10, 'MD opt.')
      end if

      write (env%unit, '(a)')
      call terminate(0)

   end subroutine finalize_xtb

end module xtb_prog_main
