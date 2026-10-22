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

subroutine citation(iunit)
integer,intent(in) :: iunit
write(iunit,'(3x,a)') &
   "Cite this work as:", &
   "* C. Bannwarth, E. Caldeweyher, S. Ehlert, A. Hansen, P. Pracht,",&
   "  J. Seibert, S. Spicher, S. Grimme, WIREs Comput. Mol. Sci., 2020, 11,",&
   "  e01493. DOI: 10.1002/wcms.1493",&
   "",&
   "for GFN2-xTB:",&
   "* C. Bannwarth, S. Ehlert and S. Grimme., J. Chem. Theory Comput., 2019,",&
   "  15, 1652-1671. DOI: 10.1021/acs.jctc.8b01176",&
   "for GFN1-xTB:",&
   "* S. Grimme, C. Bannwarth, P. Shushkov, J. Chem. Theory Comput., 2017,",&
   "  13, 1989-2009. DOI: 10.1021/acs.jctc.7b00118", &
   "for GFN0-xTB:",&
   "* P. Pracht, E. Caldeweyher, S. Ehlert, S. Grimme, ChemRxiv, 2019, preprint.",&
   "  DOI: 10.26434/chemrxiv.8326202.v1",&
   "for GFN-FF:",&
   "* S. Spicher and S. Grimme, Angew. Chem. Int. Ed., 2020, 59, 15665-15673.",&
   "  DOI: 10.1002/anie.202004239", &
   "",&
   "for ALPB and GBSA implicit solvation:",&
   "* S. Ehlert, M. Stahn, S. Spicher, S. Grimme, J. Chem. Theory Comput.,", &
   "  2021, 17, 4250-4261. DOI: 10.1021/acs.jctc.1c00471", &
   "",&
   "for ddCOSMO and CPCM-X implicit solvation:",&
   "* M. Stahn, S. Ehlert, S. Grimme, J. Phys. Chem. A,", &
   "  2023, XX, XXXX-XXXX. DOI: 10.1021/acs.jpca.3c04382", &
   "",&
   "for DFT-D4:",&
   "* E. Caldeweyher, C. Bannwarth and S. Grimme, J. Chem. Phys., 2017,",&
   "  147, 034112. DOI: 10.1063/1.4993215", &
   "* E. Caldeweyher, S. Ehlert, A. Hansen, H. Neugebauer, S. Spicher,", &
   "  C. Bannwarth and S. Grimme, J. Chem. Phys., 2019, 150, 154122.", &
   "  DOI: 10.1063/1.5090222", &
   "* E. Caldeweyher, J.-M. Mewes, S. Ehlert and S. Grimme, Phys. Chem. Chem. Phys.",&
   "  2020, 22, 8499-8512. DOI: 10.1039/D0CP00502A",&
   "",&
   "for sTDA-xTB:",&
   "* S. Grimme and C. Bannwarth, J. Chem. Phys., 2016, 145, 054103.",&
   "  DOI: 10.1063/1.4959605",&
   "",&
   "in the mass-spec context:",&
   "* V. Asgeirsson, C. Bauer and S. Grimme, Chem. Sci., 2017, 8, 4879.",&
   "  DOI: 10.1039/c7sc00601b",&
   "* J. Koopman and S. Grimme, ACS Omega 2019, 4, 12, 15120-15133.",&
   "  DOI: 10.1021/acsomega.9b02011",&
   "",&
   "for metadynamics refer to:",&
   "* S. Grimme, J. Chem. Theory Comput., 2019, 155, 2847-2862", &
   "  DOI: 10.1021/acs.jctc.9b00143", &
   "",&
   "for SPH calculations refer to:",&
   "* S. Spicher and S. Grimme, J. Chem. Theory Comput., 2021, 17, 1701-1714", &
   "  DOI: 10.1021/acs.jctc.0c01306", &
   "",&
   "for ONIOM refer to:",&
   "* C. Plett, A. Katbashev, S. Ehlert, S. Grimme, M. Bursch,",&
   "  Phys. Chem. Chem. Phys., 2023, 25, 17860-17868. DOI: 10.1039/D3CP02178E",&
   "",&
   "for DIPRO refer to:",&
   "* J. Kohn, N. Gildemeister, S. Grimme, D. Fazzi, A. Hansen,",&
   "  J. Chem. Phys. 159, 144106 (2023). DOI: 10.1063/5.0167484",&
   "",&
   "for PTB refer to:",&
   "* S. Grimme, M. Mueller, A. Hansen, J. Chem. Phys., 2023, 158, 124111.",&
   "  DOI: 10.1063/5.0137838",&
   "",&
   "for O1NumHess refer to:",&
   "* B. Wang, S. Luo, Z. Wang, W. Liu,",&
   "  J. Chem. Theory Comput. 21, 10893-10909.",&
   "  DOI: 10.1021/acs.jctc.5c01354",&
   "",&
   "with help from (in alphabetical order)",&
   "P. Atkinson, C. Bannwarth, F. Bohle, G. Brandenburg, E. Caldeweyher,", &
   "M. Checinski, S. Dohm, S. Ehlert, S. Ehrlich, M. Friede, T. Froitzheim,", &
   "I. Gerasimov, C. Hölzer, A. Katbashev, J. Kohn, J. Koopman, C. Lavigne,", &
   "S. Lehtola, F. März, M. Müller, F. Musil, H. Neugebauer, J. Pisarek,", &
   "C. Plett, P. Pracht, F. Pultar, J. Seibert, L. M. Seidler, P. Shushkov, S. Spicher,", &
   "M. Stahn, M. Steiner, T. Strunk, J. Stückrath, T. Rose, and J. Unsleber", &
   ""
end subroutine citation

subroutine help(iunit, topic, status, values_only)
   use xtb_features, only : get_xtb_feature
   !> Write the `--help` page, the full text of one option, one category of
   !> the grouped index, or (with `values_only`) the machine readable value
   !> table of the `--help-values` flag.
   implicit none
   !> Unit to write to
   integer, intent(in) :: iunit
   !> Requested topic: an option spelling with or without its dashes, a
   !> category key, or an empty string for the grouped index
   character(len=*), intent(in) :: topic
   !> 0 when a page was written, 1 for an unknown topic
   integer, intent(out) :: status
   !> Write the `--help-values` table instead of a help page
   logical, intent(in), optional :: values_only

   !> One documented command line option of the page.
   type :: help_option
      !> The spellings exactly as the page prints them; this string starts
      !> with `-` so the help coverage gate reads it as an option line
      character(len=:), allocatable :: spellings
      !> Category key, see cat_key
      character(len=:), allocatable :: category
      !> One line gloss of the grouped index
      character(len=:), allocatable :: summary
      !> The full text of the option, one 4 space indented line each
      character(len=160), allocatable :: detail(:)
      !> `/spelling/` keys of every spelling, lowercased, for the lookup
      character(len=:), allocatable :: aliases
      !> The same keys with their original case; matched before `aliases` so
      !> the one letter pair `-v`/`-V` stays distinguishable
      character(len=:), allocatable :: aliases_exact
      !> Legal values of the option, one per element; unallocated when the
      !> option is a flag or takes a free value.  Printed as `Values: ...`
      !> by the topic page and dumped by `--help-values`
      character(len=16), allocatable :: values(:)
   end type help_option

   !> Number of categories and of options of the page
   integer, parameter :: n_cat = 17, n_opt = 223
   !> Legal values of every option with a constrained value set.  One named
   !> array per set, referenced by the `values=` argument of `add`, so a set
   !> lives in exactly one place: `--help <topic>` prints it as a uniform
   !> `Values: ...` line and `--help-values` dumps it for a generator.  The
   !> order is the one the prose lists the values in.
   character(len=16), parameter :: v_optlevel(10) = [character(len=16) :: &
      "crude", "sloppy", "loose", "lax", "normal", "tight", "verytight", &
      "vtight", "extreme", "superloose"]
   !> The upstream level keywords of the native `optlevel2int` bridge, the set
   !> `--bhess`/`--metaopt` reach: the eight upstream names plus the `vtight`
   !> spelling, and **without** the suite-only `superloose` (the native bridge
   !> silently maps an unknown keyword to `normal`).
   character(len=16), parameter :: v_optlevel_native(9) = [character(len=16) :: &
      "crude", "sloppy", "loose", "lax", "normal", "tight", "verytight", &
      "vtight", "extreme"]
   character(len=16), parameter :: v_level(5) = [character(len=16) :: &
      "superloose", "loose", "normal", "tight", "vtight"]
   character(len=16), parameter :: v_minimiser(6) = [character(len=16) :: &
      "sd", "cg", "lbfgs", "bfgs", "diis", "rfo"]
   character(len=16), parameter :: v_saddle(8) = [character(len=16) :: &
      "prfo", "bprfo", "dimer", "neb", "cineb", "autoneb", "string", "gsm"]
   character(len=16), parameter :: v_irc(4) = [character(len=16) :: &
      "gs", "hpc", "lqa", "eulerpc"]
   character(len=16), parameter :: v_method(26) = [character(len=16) :: &
      "auto", "sd", "cg", "lbfgs", "bfgs", "diis", "rfo", "prfo", "bprfo", &
      "dimer", "neb", "cineb", "autoneb", "string", "gsm", "gs", "hpc", &
      "lqa", "eulerpc", "ancopt", "rf", "lancopt", "l-ancopt", "fire", &
      "inertial", "pbc_lbfgs"]
   character(len=16), parameter :: v_hessian(3) = [character(len=16) :: &
      "numerical", "model", "file"]
   character(len=16), parameter :: v_restart(3) = [character(len=16) :: &
      "off", "pin", "warm"]
   character(len=16), parameter :: v_coord(3) = [character(len=16) :: &
      "distance", "angle", "dihedral"]
   character(len=16), parameter :: v_scanmode(2) = [character(len=16) :: &
      "relaxed", "rigid"]
   character(len=16), parameter :: v_initmode(2) = [character(len=16) :: &
      "h0", "force"]
   character(len=16), parameter :: v_ninit(2) = [character(len=16) :: &
      "random", "force"]
   character(len=16), parameter :: v_ndim(2) = [character(len=16) :: "1", "2"]
   character(len=16), parameter :: v_gfn(4) = [character(len=16) :: &
      "0", "1", "2", "ff"]
   character(len=16), parameter :: v_color(3) = [character(len=16) :: &
      "auto", "always", "never"]
   character(len=16), parameter :: v_wrtopo(13) = [character(len=16) :: &
      "etot", "gnorm", "nb", "bpair", "alist", "blist", "tlist", "vtors", &
      "vbond", "vangl", "hbbond", "eeq", "json"]
   character(len=16), parameter :: v_alpb(25) = [character(len=16) :: &
      "acetone", "acetonitrile", "aniline", "benzaldehyde", "benzene", &
      "ch2cl2", "chcl3", "cs2", "dioxane", "dmf", "dmso", "ethanol", &
      "ether", "ethylacetate", "furane", "hexadecane", "hexane", "methanol", &
      "nitromethane", "octanol", "woctanol", "phenol", "toluene", "thf", &
      "water"]
   character(len=16), parameter :: v_gbsa(14) = [character(len=16) :: &
      "acetone", "acetonitrile", "benzene", "ch2cl2", "chcl3", "cs2", "dmf", &
      "dmso", "ether", "water", "methanol", "hexane", "thf", "toluene"]
   !> Category keys, in print order
   character(len=16), parameter :: cat_key(n_cat) = [character(len=16) :: &
      "runtypes", &
      "suite", &
      "levels", &
      "minimisation", &
      "scan", &
      "saddle", &
      "band", &
      "irc", &
      "frequencies", &
      "hamiltonian", &
      "electronic", &
      "external", &
      "solvation", &
      "analysis", &
      "output", &
      "workflow", &
      "performance"]
   !> Category headings, aligned with cat_key
   character(len=64), parameter :: cat_title(n_cat) = [character(len=64) :: &
      "Native run types", &
      "Landscape suite tasks", &
      "Method selector and level keywords", &
      "Minimisation options (--lopt)", &
      "Scan options (--lscan)", &
      "Saddle search options (--lts)", &
      "Band and chain-of-states options (--lts)", &
      "IRC options (--lirc)", &
      "Frequency and thermochemistry options (--lvib)", &
      "Method and Hamiltonian", &
      "Molecule and electronic-structure settings", &
      "External engines and embedding", &
      "Solvation and environment", &
      "Analysis and property requests", &
      "Output and machine-readable", &
      "Input, workflow and diagnostics", &
      "Performance and machine settings"]
   type(help_option) :: opts(n_opt)
   integer :: n, i, icat
   character(len=:), allocatable :: norm, norm_lower

   !> Basic ANSI SGR sequences (8/16 colour only, reset after each run);
   !> they are decoration and are inserted into, never added to, the text
   character(len=*), parameter :: csi = achar(27) // '['
   character(len=*), parameter :: sgr_reset = csi // '0m'
   character(len=*), parameter :: sgr_heading = csi // '1;36m'
   character(len=*), parameter :: sgr_spelling = csi // '1m'
   character(len=*), parameter :: sgr_usage = csi // '1;32m'
   !> Whether the page is coloured; decided once here and passed to every
   !> writer (no module-level state, PLAN section 5)
   logical :: use_color

   use_color = colour_wanted()
   n = 0
   call add("--scc, --sp", "runtypes", &
      "single point calculation", &
      [character(len=160) :: "    performs a single point calculation"])
   call add("--grad", "runtypes", &
      "gradient calculation", &
      [character(len=160) :: "    performs a gradient calculation"])
   call add("-o, --opt [LEVEL]", "runtypes", &
      "minimise with the landscape suite", &
      [character(len=160) :: &
      "    optimize the geometry with the landscape task suite (stage 3); the legacy", &
      "    ANCopt is reachable as --opt --method ancopt, and the level keyword follows", &
      "    the upstream names crude, sloppy, loose, lax, normal (default), tight,", &
      "    verytight, extreme, the spelling vtight and the suite-only superloose; the", &
      "    run echoes the level it uses"], values=v_optlevel)
   call add("--cycles [INT]", "runtypes", &
      "maximum number of optimisation cycles", &
      [character(len=160) :: "    maximum number of optimization cycles."])
   call add("--hess", "runtypes", &
      "numerical Hessian of the input geometry", &
      [character(len=160) :: "    perform a numerical hessian calculation on input geometry"])
   call add("--ohess [LEVEL]", "runtypes", &
      "suite minimisation then harmonic analysis", &
      [character(len=160) :: &
      "    minimize the geometry with the landscape task suite and run its harmonic", &
      "    frequency analysis (the legacy ancopt(3) plus numerical Hessian path is", &
      "    --ohess --method ancopt)"], values=v_optlevel)
   call add("--bhess [LEVEL]", "runtypes", &
      "biased numerical Hessian on an ANCopt geometry", &
      [character(len=160) :: &
      "    perform a biased numerical hessian calculation on an ancopt(3) optimized geometry", &
      "", &
      "    Molecular Hessian runs with more than one atom, no frozen atoms, and no periodicity", &
      "    automatically print compliance constants, relaxed force constants, and bond local-mode", &
      "    frequencies, and write compliance.dat."], values=v_optlevel_native)
   call add("--o1nh", "runtypes", &
      "O1NumHess numerical Hessian (ODLR)", &
      [character(len=160) :: "    perform the numerical hessian calculation using the ODLR approximation (O1NumHess)"])
   call add("--imagmin REAL", "runtypes", &
      "O1NumHess imaginary repair: lower bound (cm^-1)", &
      [character(len=160) :: &
      "    magnitudes below this value (cm⁻¹) are ignored in the O1NumHess", &
      "    imaginary frequency repair (default 5)"])
   call add("--imagmax REAL", "runtypes", &
      "O1NumHess imaginary repair: upper bound (cm^-1)", &
      [character(len=160) :: &
      "    deepest imaginary frequency (cm⁻¹) considered in the O1NumHess", &
      "    imaginary frequency repair (default 200)"])
   call add("--md", "runtypes", &
      "molecular dynamics from the input geometry", &
      [character(len=160) :: "    molecular dynamics simulation on start geometry"])
   call add("--metadyn [INT]", "runtypes", &
      "metadynamics with int snapshots", &
      [character(len=160) :: &
      "    meta dynamics simulation on start geometry, saving int snapshots of the", &
      "    trajectory to bias the simulation"])
   call add("--bias-input, --gesc FILE", "runtypes", &
      "RMSD bias file for --metadyn", &
      [character(len=160) :: "    file with the RMSD bias of the metadynamics run (--metadyn)"])
   call add("--omd", "runtypes", &
      "MD from an ANCopt-minimised geometry", &
      [character(len=160) :: &
      "    molecular dynamics simulation on ancopt(3) optimized geometry,", &
      "    a loose optimization level will be chosen"])
   call add("--metaopt [LEVEL]", "runtypes", &
      "ANCopt, then metadynamics for other minima", &
      [character(len=160) :: &
      "    call ancopt(3) to perform a geometry optimization, then try to find other", &
      "    minimas by meta dynamics"], values=v_optlevel_native)
   call add("--path [FILE]", "runtypes", &
      "metadynamics path to a product structure", &
      [character(len=160) :: &
      "    use meta dynamics to calculate a path from the input geometry to the given", &
      "    product structure"])
   call add("--modef INT", "runtypes", &
      "vibrational mode following (mode INT)", &
      [character(len=160) :: &
      "    modefollowing algorithm. INT specifies the mode that should be used for the", &
      "    modefollowing."])
   call add("--screen", "runtypes", &
      "optimisation loop over an ensemble", &
      [character(len=160) :: &
      "    optimisation loop over an ensemble: QMDFF generates it when a file", &
      "    'solvent' exists, otherwise it is read from the file xtbscreen.xyz"])
   call add("--mdopt", "runtypes", &
      "re-optimise an MD trajectory (xtb.trj)", &
      [character(len=160) :: &
      "    re-optimise the structures of the molecular dynamics trajectory xtb.trj", &
      "    at the current optimisation level; the energies are written to", &
      "    xtb_ensemble.xyz"])
   ! `--siman` and `--gmd` are deliberately absent from this table: both are
   ! parsed by `src/prog/main.F90` but are not registered run types, so a run
   ! aborts in `set_runtyp` with an honest message; advertising them in
   ! `--help` only promises a capability the build does not have.  The two
   ! spellings are on the help/man gates' documented allowlist for that reason.
   call add("--ewin REAL", "runtypes", &
      "energy window of the conformational search (Eh)", &
      [character(len=160) :: "    energy window of the conformational search (default 20.0)"])
   call add("--ljson", "suite", &
      "write landscape_<TASK>.json for every task", &
      [character(len=160) :: &
      "    write the structured result landscape_<TASK>.json for every task: native", &
      "    (--sp, --grad, --opt, --hess, --ohess, --md, ...) and landscape (--lopt,", &
      "    --lscan, --lts, --lirc, --lvib).  The document carries a uniform core", &
      "    block (task, status, converged, energy, n_atoms, geometry, provenance)", &
      "    plus the task specific counters and convergence history, so one parser", &
      "    reads every task; it is meant for regression comparison and is written", &
      "    independently of the other output files (format 2)"])
   call add("--lopt", "suite", &
      "minimise with the landscape suite", &
      [character(len=160) :: &
      "    minimize the geometry with the landscape task suite (src/landscape); the", &
      "    algorithm is selected with --lopt-engine and the thresholds with --lopt-level", &
      "    an xcontrol $scan block is not executed by this run type and aborts the run", &
      "    (use --opt --method ancopt for that block, or --lscan with --lscan-coord)"])
   call add("--lscan", "suite", &
      "scan one or two internal coordinates", &
      [character(len=160) :: &
      "    scan an internal coordinate (distance, angle or dihedral) with the landscape", &
      "    task suite (src/landscape); a relaxed scan relaxes each scan point under a", &
      "    constraint potential, a rigid scan only moves the geometry. The result is", &
      "    written to xtbscan.log with the de-biased energy of every frame, together", &
      "    with a report of the constraint residual and, for a bidirectional scan, of", &
      "    the hysteresis between the forward and the reverse branch"])
   call add("--lts", "suite", &
      "transition-state or band search", &
      [character(len=160) :: &
      "    transition state search with the landscape task suite (src/landscape);", &
      "    writes xtbts.xyz and xtbts.log and confirms the imaginary mode count"])
   call add("--lirc", "suite", &
      "intrinsic reaction coordinate", &
      [character(len=160) :: &
      "    intrinsic reaction coordinate with the landscape task suite (src/landscape);", &
      "    the input geometry is the transition state and must have exactly one imaginary", &
      "    mode; writes irc_forward.xyz, irc_backward.xyz, irc_path.xyz and irc.log"])
   call add("--lvib", "suite", &
      "frequencies, modes and thermochemistry", &
      [character(len=160) :: &
      "    harmonic frequencies, normal modes and thermochemistry with the landscape task", &
      "    suite (src/landscape); writes vibspectrum, vib_normal_modes, g98.out and hessian", &
      "    with the names and formats of --hess", &
      "    an xcontrol $thermo block (temp, sthr, imagthr, scale) is honoured unless", &
      "    the matching --lvib-temp, --lvib-rotor-cutoff, --lvib-imag-cutoff or", &
      "    --lvib-freq-scal option is given: explicit command line > $thermo > default"])
   call add("-3, -2, -1, 0, 1, 2, 3", "levels", &
      "numeric aliases of the level keywords", &
      [character(len=160) :: &
      "    numeric aliases of the level keywords.  For --opt and --ohess they map", &
      "    onto the landscape levels (-3/-2 superloose, -1 loose, 0 normal, 1 tight,", &
      "    2/3 vtight); the native run types --bhess, --omd and --metaopt use", &
      "    their own table (-3 crude, -2 sloppy, 3 extreme)"])
   call add("--method NAME", "levels", &
      "algorithm of the current task (canonical selector)", &
      [character(len=160) :: &
      "    select the method of the current landscape run type; the value domain is", &
      "    filtered by the run type and an unknown or mismatched name aborts the run.", &
      "    --lopt and --lscan accept the minimisers sd, cg, lbfgs, bfgs, diis, rfo;", &
      "    --lts accepts the saddle searches prfo, bprfo, dimer, neb, cineb, autoneb,", &
      "    string, gsm; --lirc accepts the integrators gs, hpc, lqa, eulerpc; auto keeps", &
      "    the default of the task.  With --opt and --ohess the native engines ancopt", &
      "    (alias rf), lancopt (alias l-ancopt), fire (alias inertial) and pbc_lbfgs are", &
      "    accepted; a landscape name selects the suite that the redeployed --opt runs.", &
      "    --method is the canonical name of --lopt-engine, --lscan-method,", &
      "    --lts-method and --lirc-method: the same value is idempotent, a different", &
      "    value through two of the flags aborts the run.  The name is not", &
      "    case-sensitive"], values=v_method)
   call add("--lopt-engine NAME", "minimisation", &
      "minimiser: sd, cg, lbfgs, bfgs, diis or rfo", &
      [character(len=160) :: &
      "    minimisation algorithm of the landscape suite, one of sd, cg, lbfgs,", &
      "    bfgs, diis or rfo (default lbfgs).  Alias of --method.  The name is", &
      "    not case-sensitive"], values=v_minimiser)
   call add("--lopt-level LEVEL", "minimisation", &
      "convergence level of the minimisation", &
      [character(len=160) :: &
      "    convergence level of the landscape suite (thresholds, energy tolerance, SCC", &
      "    accuracy and the step budget), one of superloose, loose,", &
      "    normal (default), tight or vtight.  The level name is not case-sensitive"], &
      values=v_level)
   call add("--lopt-hessian SOURCE", "minimisation", &
      "Hessian backend: numerical, model or file", &
      [character(len=160) :: &
      "    Hessian backend of the landscape suite (minimisation and Hessian", &
      "    based steps), one of numerical (default, central differences of the", &
      "    electronic gradient), model (Lindh-D2 model Hessian) or file (read the", &
      "    Hessian written by an earlier run, named by --lopt-hessian-file)"], &
      values=v_hessian)
   call add("--lopt-max-steps N", "minimisation", &
      "step budget of the minimisation", &
      [character(len=160) :: &
      "    step budget of the minimisation; the default is derived from the", &
      "    convergence level (3 steps per atom at level normal, clamped to", &
      "    200..10000), --cycles and the level change it"])
   call add("--lopt-recalc N", "minimisation", &
      "rebuild the exact Hessian every N RFO steps", &
      [character(len=160) :: &
      "    rebuild the exact Hessian through the provider every N RFO steps of", &
      "    the minimisation (default 5); 0 disables the rebuild and keeps the", &
      "    Hessian update alone, an interval above the step budget never fires, and", &
      "    a negative value is reported and behaves like 0.  Only --lopt-engine rfo", &
      "    refreshes its Hessian at all; the saddle search and the reaction path", &
      "    have their own --lts-recalc and --lirc-recalc"])
   call add("--lopt-hess-step VALUE", "minimisation", &
      "FD step of the numerical Hessian (bohr)", &
      [character(len=160) :: &
      "    finite difference step of the numerical Hessian in bohr (default 0.005);", &
      "    only --lopt-engine rfo builds a numerical Hessian; the other five", &
      "    engines (including the default lbfgs) never do, and the value is", &
      "    ignored for --lopt-hessian model or file.  The run echoes the", &
      "    effective value and marks it as unused in that case.  The canonical", &
      "    spelling, mirrored by --lscan-hess-step and by --lts-hess-step /", &
      "    --lirc-hess-step / --lvib-hess-step for the other landscape tasks"])
   call add("--lopt-sccacc VALUE", "minimisation", &
      "SCC accuracy of the FD Hessian points", &
      [character(len=160) :: &
      "    SCC accuracy of the displaced single points of the numerical Hessian", &
      "    (default 0.3), consumed by the same engine and Hessian source as", &
      "    --lopt-hess-step"])
   call add("--lopt-restart MODE", "minimisation", &
      "wavefunction restart policy: off, pin or warm", &
      [character(len=160) :: "    wavefunction restart policy of the engine, off (default), pin or warm"], &
      values=v_restart)
   call add("--lopt-hessian-file NAME", "minimisation", &
      "Hessian file read by --lopt-hessian file", &
      [character(len=160) :: "    file read by --lopt-hessian file (TURBOMOLE Hessian layout)"])
   call add("--lscan-coord TYPE", "scan", &
      "type of the first coordinate: distance, angle, dihedral", &
      [character(len=160) :: &
      "    type of the first scanned coordinate, one of distance, angle or dihedral", &
      "    (the type name is not case-sensitive)"], values=v_coord)
   call add("--lscan-atoms LIST", "scan", &
      "atom indices of the first coordinate", &
      [character(len=160) :: &
      "    1-based atom indices of the first coordinate, comma or blank separated", &
      "    (two for a distance, three for an angle, four for a dihedral)"])
   call add("--lscan-step VALUE", "scan", &
      "step size of the first coordinate", &
      [character(len=160) :: &
      "    step size of the first coordinate in angstrom (distance) or degree", &
      "    (angle/dihedral)"])
   call add("--lscan-nstep N", "scan", &
      "number of steps of the first coordinate", &
      [character(len=160) :: "    number of steps of the first coordinate"])
   call add("--lscan-coord2, --lscan-atoms2, --lscan-step2, --lscan-nstep2", "scan", &
      "2nd coordinate", &
      [character(len=160) :: &
      "    same as above for the second scanned coordinate; any of them switches the", &
      "    scan to a two-dimensional Cartesian product grid"])
   call add("--lscan-mode MODE", "scan", &
      "scan mode: relaxed or rigid", &
      [character(len=160) :: "    scan mode, relaxed (default) or rigid (not case-sensitive)"], &
      values=v_scanmode)
   call add("--lscan-method NAME", "scan", &
      "inner minimiser of a relaxed scan", &
      [character(len=160) :: &
      "    inner minimiser of a relaxed scan, one of sd, cg, lbfgs (default), bfgs,", &
      "    diis or rfo.  Alias of --method.  The name is not case-sensitive"], &
      values=v_minimiser)
   call add("--lscan-level LEVEL", "scan", &
      "convergence level of the inner minimisation", &
      [character(len=160) :: &
      "    convergence level of the inner minimisation (thresholds, energy tolerance,", &
      "    SCC accuracy and the step budget), one of superloose, loose,", &
      "    normal (default), tight or vtight"], values=v_level)
   call add("--lscan-ndim N", "scan", &
      "number of scanned coordinates, 1 or 2", &
      [character(len=160) :: "    number of scanned coordinates, 1 (default) or 2"], &
      values=v_ndim)
   call add("--lscan-restart MODE", "scan", &
      "wavefunction policy of the scan points", &
      [character(len=160) :: &
      "    wavefunction policy of the scan points, off (default), pin or warm", &
      "    (warm_start reuses only the geometry, never the wavefunction)"], &
      values=v_restart)
   call add("--lscan-oneway", "scan", &
      "forward branch only (no hysteresis report)", &
      [character(len=160) :: "    scan the forward branch only (no reverse branch, no hysteresis report)"])
   call add("--lscan-cold-start", "scan", &
      "every scan point starts from the input geometry", &
      [character(len=160) :: "    start every scan point from the input geometry instead of the previous point"])
   call add("--lscan-fc VALUE", "scan", &
      "force constant of the constraint potential", &
      [character(len=160) :: "    force constant of the constraint potential in hartree/bohr^2 (default 0.5)"])
   call add("--lscan-expo VALUE", "scan", &
      "exponent of the constraint potential", &
      [character(len=160) :: "    exponent of the constraint potential (default 2.0)"])
   call add("--lscan-max-steps N", "scan", &
      "step budget of every inner relaxation", &
      [character(len=160) :: &
      "    step budget of the inner relaxation of every scan point; the default", &
      "    is derived from the convergence level (3 steps per atom at level", &
      "    normal, clamped to 200..10000), --cycles and the level change it"])
   call add("--lscan-hess-step VALUE", "scan", &
      "FD step of the scan-point Hessian (bohr)", &
      [character(len=160) :: &
      "    finite difference step of the numerical Hessian of the scan points in", &
      "    bohr (default 0.005); only an inner minimiser of --lscan-method rfo", &
      "    builds a numerical Hessian (the default lbfgs never does), a rigid", &
      "    scan has no inner minimisation at all, and the value is ignored for", &
      "    --lscan-hessian model or file.  The run echoes the effective value", &
      "    and marks it as unused in those cases"])
   call add("--lscan-sccacc VALUE", "scan", &
      "SCC accuracy of the scan-point Hessian", &
      [character(len=160) :: &
      "    SCC accuracy of the displaced single points of the numerical Hessian", &
      "    of the scan points (default 0.3), consumed by the same inner", &
      "    minimiser and Hessian source as --lscan-hess-step"])
   call add("--lscan-hessian SOURCE", "scan", &
      "Hessian backend of the scan points", &
      [character(len=160) :: "    Hessian backend of the scan points, one of numerical (default), model or file"], &
      values=v_hessian)
   call add("--lscan-hessian-file NAME", "scan", &
      "Hessian file read by --lscan-hessian file", &
      [character(len=160) :: "    file read by --lscan-hessian file (TURBOMOLE Hessian layout)"])
   call add("--lts-method NAME", "saddle", &
      "search method: prfo, bprfo, dimer, neb, ...", &
      [character(len=160) :: &
      "    search method, one of prfo (default), bprfo (restricted step), dimer", &
      "    (first derivatives only, no Hessian), neb, cineb (climbing image NEB),", &
      "    autoneb (incremental path building) or string / gsm (growing string).", &
      "    Alias of --method.  The method name is not case-sensitive"], &
      values=v_saddle)
   call add("--lts-trust VALUE", "saddle", &
      "initial trust radius (bohr)", &
      [character(len=160) :: "    initial trust radius of the step in bohr (default 0.37795)"])
   call add("--lts-trust-max VALUE", "saddle", &
      "upper bound of the trust radius (bohr)", &
      [character(len=160) :: "    upper bound of the trust radius in bohr (default 1.8897)"])
   call add("--lts-hlow VALUE", "saddle", &
      "curvature floor of the downhill modes", &
      [character(len=160) :: "    curvature floor of the downhill modes in hartree/bohr^2 (default 0.01)"])
   call add("--lts-max-steps N", "saddle", &
      "maximum number of macro iterations", &
      [character(len=160) :: &
      "    maximum number of macro iterations; the default is derived from the", &
      "    convergence level (3 steps per atom at level normal, clamped to", &
      "    200..10000), --cycles and the level change it"])
   call add("--lts-level LEVEL", "saddle", &
      "convergence level of the search", &
      [character(len=160) :: &
      "    convergence level (thresholds, energy tolerance, SCC accuracy and the step", &
      "    budget), one of superloose, loose, normal, tight, vtight"], values=v_level)
   call add("--lts-mode-overlap-min VALUE", "saddle", &
      "mode-following guard: reselect below overlap", &
      [character(len=160) :: "    mode following guard: reselect the transition vector below this overlap (default 0.3)"])
   call add("--lts-alpha VALUE", "saddle", &
      "fraction of the trust radius on the uphill subspace", &
      [character(len=160) :: "    fraction of the trust radius assigned to the uphill subspace by bprfo (default 0.5)"])
   call add("--lts-micro-cycles N", "saddle", &
      "restricted-step micro iterations per macro step", &
      [character(len=160) :: "    maximum number of restricted step micro iterations (default 50)"])
   call add("--lts-target-mode N", "saddle", &
      "which negative mode is followed", &
      [character(len=160) :: "    which negative mode is followed, 1 = most negative (default 1)"])
   call add("--lts-init-mode MODE", "saddle", &
      "initial transition vector: h0 or force", &
      [character(len=160) :: &
      "    initial transition vector, h0 (most negative Hessian mode, default) or force", &
      "    (largest Cartesian gradient component direction)"], values=v_initmode)
   call add("--lts-escape-step VALUE", "saddle", &
      "escape walk from a minimum (0 disables)", &
      [character(len=160) :: &
      "    displacement along the lowest internal mode in bohr (default 0.2), applied", &
      "    when the start geometry is a minimum so P-RFO has a direction to ascend", &
      "    along; 0 disables the escape"])
   call add("--lts-delta VALUE", "saddle", &
      "dimer separation (bohr)", &
      [character(len=160) :: &
      "    dimer separation, the distance of each image from the midpoint in bohr", &
      "    (default 0.018897, i.e. 0.01 angstrom)"])
   call add("--lts-dr VALUE", "saddle", &
      "FD step of the dimer translation probe", &
      [character(len=160) :: "    finite difference step of the dimer translation probe in bohr (default 0.0018897)"])
   call add("--lts-rot-max-iter N", "saddle", &
      "maximum dimer rotations per macro iteration", &
      [character(len=160) :: "    maximum number of dimer rotations per macro iteration (default 5)"])
   call add("--lts-rot-alpha VALUE", "saddle", &
      "damping of the dimer rotation angle", &
      [character(len=160) :: "    damping of the dimer rotation angle of Eq. 13, 1 = undamped (default 0.5)"])
   call add("--lts-rot-frms VALUE", "saddle", &
      "RMS torque of dimer rotation convergence", &
      [character(len=160) :: "    RMS torque below which the dimer rotation converged, in hartree/bohr (default 2.65e-4)"])
   call add("--lts-ninit MODE", "saddle", &
      "initial dimer axis: random or force", &
      [character(len=160) :: &
      "    initial dimer axis, random (default, seeded by --lts-seed) or force", &
      "    (largest Cartesian gradient component direction)"], values=v_ninit)
   call add("--lts-seed N", "saddle", &
      "seed of the random initial dimer axis", &
      [character(len=160) :: "    seed of the random initial dimer axis (default 42)"])
   call add("--lts-step-max VALUE", "saddle", &
      "maximum dimer translation step (bohr)", &
      [character(len=160) :: "    maximum dimer translation step in bohr (default 0.28346)"])
   call add("--lts-step0 VALUE", "saddle", &
      "upper bound of the first dimer translation", &
      [character(len=160) :: &
      "    upper bound of the first dimer translation in bohr (default 0.37795);", &
      "    the dR curvature probe can shorten the actual first displacement"])
   call add("--lts-translate-max-iter N", "saddle", &
      "quick-min updates per dimer translation", &
      [character(len=160) :: "    maximum number of quick-min displacement updates per dimer translation (default 1)"])
   call add("--lts-recalc N", "saddle", &
      "recompute the exact Hessian every N steps", &
      [character(len=160) :: "    recompute the exact Hessian every N macro steps, 0 disables (default 1)"])
   call add("--lts-no-verify", "saddle", &
      "skip the frequency confirmation", &
      [character(len=160) :: "    skip the frequency confirmation of the located stationary point"])
   call add("--lts-no-rigid", "saddle", &
      "keep the rigid-body components in the step", &
      [character(len=160) :: "    keep the rigid body (translation/rotation) components in the Cartesian step"])
   call add("--lts-restart MODE", "saddle", &
      "wavefunction restart policy: off, pin or warm", &
      [character(len=160) :: "    wavefunction restart policy of the engine, off (default), pin or warm"], &
      values=v_restart)
   call add("--lts-hessian SOURCE", "saddle", &
      "Hessian backend: numerical, model or file", &
      [character(len=160) :: "    Hessian backend, one of numerical (default), model (Lindh-D2) or file"], &
      values=v_hessian)
   call add("--lts-hessian-file NAME", "saddle", &
      "Hessian file read by --lts-hessian file", &
      [character(len=160) :: "    file read by --lts-hessian file (TURBOMOLE Hessian layout)"])
   call add("--lts-hess-step VALUE", "saddle", &
      "FD step of the confirmation Hessian (bohr)", &
      [character(len=160) :: &
      "    finite difference step of the numerical confirmation Hessian in bohr", &
      "    (default 0.005); the canonical name of --lts-step"])
   call add("--lts-step VALUE", "saddle", &
      "alias of --lts-hess-step", &
      [character(len=160) :: &
      "    alias of --lts-hess-step: finite difference step of the numerical", &
      "    confirmation Hessian in bohr (default 0.005);", &
      "    a mode inside 0 +/- 3 cm^-1 is re-decided with 0.002 bohr and sccacc 0.1"])
   call add("--lts-sccacc VALUE", "saddle", &
      "SCC accuracy of the confirmation Hessian", &
      [character(len=160) :: "    SCC accuracy of the displaced single points of the confirmation Hessian (default 0.3)"])
   call add("--lts-product NAME", "band", &
      "product endpoint of a chain-of-states search", &
      [character(len=160) :: &
      "    product endpoint of a chain-of-states search (neb, cineb, autoneb, string,", &
      "    gsm): a single frame xyz with the same atoms in the same order"])
   call add("--lts-images N", "band", &
      "number of images of a band", &
      [character(len=160) :: &
      "    number of images of a band, the two fixed endpoints included (default 9);", &
      "    for autoneb it is the requested final number of images"])
   call add("--lts-kmin VALUE, --lts-kmax VALUE", "band", &
      "dynamic spring constant bounds", &
      [character(len=160) :: &
      "    lower/upper spring constant of the dynamic springs in hartree/bohr^2", &
      "    (defaults 1e-2 and 1e-1, calibrated on the analytic LEPS surface)"])
   call add("--lts-kdecay VALUE", "band", &
      "decay factor of the dynamic spring formula", &
      [character(len=160) :: "    decay factor of the dynamic spring formula (default 0.5)"])
   call add("--lts-no-dynamic-k", "band", &
      "single spring constant k_max for every image", &
      [character(len=160) :: "    use the single spring constant k_max for every image"])
   call add("--lts-climb", "band", &
      "climbing image stage (CI-NEB)", &
      [character(len=160) :: &
      "    switch on the climbing image stage (CI-NEB); implied by --lts-method", &
      "    cineb, autoneb, string and gsm"])
   call add("--lts-climb-switch VALUE", "band", &
      "band force that ignites the climbing image", &
      [character(len=160) :: &
      "    band force (hartree/bohr) below which the climbing image is ignited", &
      "    (default 9.7234e-3, the standard 0.5 eV/angstrom)"])
   call add("--lts-no-idpp", "band", &
      "linear interpolation instead of IDPP", &
      [character(len=160) :: "    straight linear interpolation instead of the IDPP initial path"])
   call add("--lts-no-align", "band", &
      "do not Kabsch align the images", &
      [character(len=160) :: "    do not Kabsch align the images before every band evaluation"])
   call add("--lts-relax-ends", "band", &
      "relax the two endpoints before the search", &
      [character(len=160) :: "    relax the two endpoints with the minimiser before the search"])
   call add("--lts-max-displ VALUE", "band", &
      "upper bound of one displacement component", &
      [character(len=160) :: &
      "    upper bound of one displacement component of a band step in bohr", &
      "    (default 0.37795)"])
   call add("--lts-frms-avg VALUE, --lts-frms-max VALUE", "band", &
      "band convergence thresholds", &
      [character(len=160) :: &
      "    band convergence thresholds of the average and the largest per-image", &
      "    force RMS in hartree/bohr (defaults 3e-4 and 1e-3)"])
   call add("--lts-autoneb-nstart N", "band", &
      "AutoNEB: images of the coarse start path", &
      [character(len=160) :: "    number of images of the coarse start path of an AutoNEB run (default 5)"])
   call add("--lts-autoneb-max-paths N", "band", &
      "AutoNEB: upper bound of the number of images", &
      [character(len=160) :: &
      "    hard upper bound of the number of images of an AutoNEB path (default 10);", &
      "    a request above it is an error"])
   call add("--lts-autoneb-max-depth N", "band", &
      "AutoNEB: upper bound of the insertion stages", &
      [character(len=160) :: "    hard upper bound of the number of image insertion stages (default 5)"])
   call add("--lts-autoneb-restart NAME", "band", &
      "AutoNEB: restart from an earlier path", &
      [character(len=160) :: &
      "    multi-frame xyz path of an earlier AutoNEB run (its path_autoneb.xyz),", &
      "    adopted as the initial path so the completed stages are not repeated"])
   call add("--lts-string-nnodes N", "band", &
      "string: target number of nodes", &
      [character(len=160) :: "    target number of nodes of the growing string, endpoints included (default 9)"])
   call add("--lts-string-grow-tol VALUE", "band", &
      "string: norm threshold of the growth indicator", &
      [character(len=160) :: "    norm threshold of the growth indicator of Eq. 7 in hartree/bohr (default 5e-3)"])
   call add("--lts-string-ds VALUE", "band", &
      "string: arc length increment of one growth step", &
      [character(len=160) :: "    arc length increment of one growth step of the string in bohr (default 0.5)"])
   call add("--lts-string-nu2 VALUE", "band", &
      "string: curvature floor nu2 of the spacing", &
      [character(len=160) :: &
      "    curvature floor nu2 of the node spacing condition R > |gradV.tau|/(2 nu2)", &
      "    in hartree/bohr^2 (default 0.05)"])
   call add("--lts-string-grow-steps N", "band", &
      "string: perpendicular relaxation iterations", &
      [character(len=160) :: "    maximum perpendicular relaxation iterations of one frontier node (default 200)"])
   call add("--lts-no-string-spline", "band", &
      "string: central difference tangents", &
      [character(len=160) :: "    central difference tangents of the string instead of the cubic spline tangent"])
   call add("--lts-string-history", "band", &
      "string: print the growth history", &
      [character(len=160) :: "    print the growth history (position and perpendicular force of every node)"])
   call add("--lirc-method NAME", "irc", &
      "integration method: hpc, lqa, eulerpc or gs", &
      [character(len=160) :: &
      "    integration method, one of hpc (Hessian based predictor-corrector, default:", &
      "    most accurate and 1.98 E/g per point), lqa (Page-McIver analytic quadratic", &
      "    propagation, cheapest at 1.00 E/g), eulerpc (Euler predictor-corrector, use a", &
      "    step of 0.20 or less) or gs (Gonzalez-Schlegel sphere relaxation, compatible", &
      "    reference; 1.9-2.1 E/g per point and the largest path deviation, validated", &
      "    for steps up to 0.15).", &
      "    Alias of --method.  The method name is not case-sensitive"], &
      values=v_irc)
   call add("--lirc-step VALUE", "irc", &
      "mass weighted arc length step (bohr*amu^1/2)", &
      [character(len=160) :: &
      "    step length in bohr*amu^1/2 (mass weighted arc length, default 0.10; use 0.05", &
      "    for a tight path). The clean range is 0.05-0.15 for --lirc-method gs, whose", &
      "    macro step stalls in a two cycle from 0.20 on (--lirc-osc-window); hpc and", &
      "    eulerpc are correct at 0.20. The guard only detects the cancelling (limit", &
      "    cycle) subclass of the stall: a large step can still end far from the", &
      "    minimum (measured 0.236 bohr at 0.40 forward)"])
   call add("--lirc-max-steps N", "irc", &
      "maximum number of macro steps per direction", &
      [character(len=160) :: &
      "    maximum number of macro steps of one direction; the default is", &
      "    derived from the convergence level (3 steps per atom at level", &
      "    normal, clamped to 200..10000), --cycles and the level change it"])
   call add("--lirc-level LEVEL", "irc", &
      "convergence level of the endpoint relaxation", &
      [character(len=160) :: &
      "    convergence level of the endpoint relaxation (thresholds, energy tolerance,", &
      "    SCC accuracy and the step budget), one of superloose, loose,", &
      "    normal, tight, vtight (the IRC path itself is left by the gradient and", &
      "    energy-rise termination signals, see --lirc-grad-tol/--lirc-e-rise)"], &
      values=v_level)
   call add("--lirc-swap", "irc", &
      "exchange the forward and backward directions", &
      [character(len=160) :: "    exchange the forward and backward directions defined by the phase convention"])
   call add("--lirc-oneway", "irc", &
      "integrate the forward direction only", &
      [character(len=160) :: "    integrate the forward direction only"])
   call add("--lirc-target-mode N", "irc", &
      "which negative mode defines the direction", &
      [character(len=160) :: "    which negative mode defines the initial direction, 1 = most negative (default 1)"])
   call add("--lirc-micro-cycles N", "irc", &
      "Gonzalez-Schlegel micro cycles per macro step", &
      [character(len=160) :: "    maximum number of Gonzalez-Schlegel micro cycles per macro step (default 20)"])
   call add("--lirc-micro-thresh VALUE", "irc", &
      "sphere-relaxation convergence threshold", &
      [character(len=160) :: "    convergence threshold of the sphere relaxation in bohr*amu^1/2 (default 1e-3)"])
   call add("--lirc-euler-n N", "irc", &
      "Euler/Simpson substeps of the LQA propagation", &
      [character(len=160) :: "    Euler/Simpson substeps of the analytic LQA propagation (default 5000)"])
   call add("--lirc-loose-cycles N", "irc", &
      "EulerPC coarse phase: N Euler macro steps", &
      [character(len=160) :: &
      "    EulerPC coarse phase: use the plain Euler prediction for the first N macro", &
      "    steps, 0 disables the coarse phase (default 3)"])
   call add("--lirc-mbs-tol VALUE", "irc", &
      "acceptance threshold of the mBS corrector", &
      [character(len=160) :: &
      "    acceptance threshold of the modified Bulirsch-Stoer corrector in", &
      "    bohr*amu^1/2 (default 1e-5; the HPC paper uses 1e-6)"])
   call add("--lirc-mbs-max-k N", "irc", &
      "maximum extrapolation level of the mBS corrector", &
      [character(len=160) :: "    maximum extrapolation level of the mBS corrector (default 15)"])
   call add("--lirc-mbs-points N", "irc", &
      "base Euler substeps of the mBS corrector", &
      [character(len=160) :: "    base number of Euler substeps of the mBS corrector (default 20)"])
   call add("--lirc-dwi-n N", "irc", &
      "capacity of the DWI surface fit", &
      [character(len=160) :: &
      "    capacity of the DWI surface fit (default 4). The interpolant is two point by", &
      "    construction (current and predicted point), so any N >= 2 is equivalent;", &
      "    N < 2 is rejected"])
   call add("--lirc-e-rise N", "irc", &
      "stop after N consecutive rising energies", &
      [character(len=160) :: &
      "    termination signal (c): stop after N consecutive rising energies (default 3). A", &
      "    large N does not guarantee running to the step budget: the GS limit cycle", &
      "    guard (--lirc-osc-window) is an independent signal (d)"])
   call add("--lirc-erise-steps N", "irc", &
      "alias of --lirc-e-rise", &
      [character(len=160) :: "    alias of --lirc-e-rise N"])
   call add("--lirc-grad-tol VALUE", "irc", &
      "stop when the RMS gradient falls below this", &
      [character(len=160) :: &
      "    termination signal (b): stop when the RMS gradient per Cartesian component", &
      "    falls below this value in hartree/bohr (default 1e-3)"])
   call add("--lirc-osc-window N", "irc", &
      "GS limit cycle guard: N cancelling steps", &
      [character(len=160) :: &
      "    limit cycle guard of the GS macro step: leave the path after N consecutive", &
      "    macro steps that cancel each other (default 8, 0 disables the guard). The", &
      "    guard is bitwise inert on every validated path and truncates the stalled", &
      "    87 macro step walk of Mueller-Brown at step 0.20 to 13 steps. The endpoint", &
      "    it returns is not a converged path point - it is one phase of the", &
      "    oscillation, about one step length away from the minimum - so the run is", &
      "    reported as incomplete and reached_min stays false unless an endpoint", &
      "    relaxation confirms the geometry. Only --lirc-method gs is affected"])
   call add("--lirc-adaptive", "irc", &
      "adaptive GS step length (opt-in)", &
      [character(len=160) :: &
      "    opt-in adaptive step length of the GS macro step (default off): halve the", &
      "    effective step length (floor: a sixteenth of --lirc-step) at the first", &
      "    cancelling macro step pair, so a large step recovers accuracy instead of", &
      "    being truncated at the first cancellation (the guard still applies once", &
      "    the floor is reached). It changes the validated small step paths and is", &
      "    therefore a re-baseline; every reduction is printed. Only --lirc-method", &
      "    gs is affected; the other methods ignore it"])
   call add("--lirc-recalc N", "irc", &
      "rebuild the exact Hessian every N steps", &
      [character(len=160) :: &
      "    rebuild the exact Hessian every N macro steps, -1 disables (default -1). Pure", &
      "    Bofill updates are sufficient for short paths; long or strongly curved paths", &
      "    drift and need an explicit --lirc-recalc 5 .. 10"])
   call add("--lirc-no-endpoint-relax", "irc", &
      "do not relax the two endpoints", &
      [character(len=160) :: "    do not relax the two endpoints to assign them to a minimum"])
   call add("--lirc-endpoint-relax", "irc", &
      "relax the two endpoints (default)", &
      [character(len=160) :: &
      "    relax the two endpoints to assign them to a minimum (default; the explicit", &
      "    flag restores it after a --lirc-no-endpoint-relax)"])
   call add("--lirc-rmsd-end VALUE", "irc", &
      "endpoint-to-minimum RMSD threshold (bohr)", &
      [character(len=160) :: &
      "    Kabsch aligned RMSD between an endpoint and its relaxed minimum below which", &
      "    the endpoint counts as having reached a minimum, in bohr (default 0.18897).", &
      "    The criterion needs a relaxed endpoint and is skipped (and reported as", &
      "    skipped) with --lirc-no-endpoint-relax"])
   call add("--lirc-hessian SOURCE", "irc", &
      "Hessian backend: numerical, model or file", &
      [character(len=160) :: "    Hessian backend, one of numerical (default), model (Lindh-D2) or file"], &
      values=v_hessian)
   call add("--lirc-hessian-file NAME", "irc", &
      "Hessian file read by --lirc-hessian file", &
      [character(len=160) :: "    file read by --lirc-hessian file (TURBOMOLE Hessian layout)"])
   call add("--lirc-hess-step VALUE", "irc", &
      "FD step of the numerical Hessian (bohr)", &
      [character(len=160) :: &
      "    finite difference step of the numerical Hessian in bohr (default 0.005);", &
      "    the canonical name of --lirc-step-hess (not to be confused with", &
      "    --lirc-step, the arc length of the reaction path)"])
   call add("--lirc-step-hess VALUE", "irc", &
      "alias of --lirc-hess-step", &
      [character(len=160) :: &
      "    alias of --lirc-hess-step: finite difference step of the numerical", &
      "    Hessian in bohr (default 0.005)"])
   call add("--lirc-sccacc VALUE", "irc", &
      "SCC accuracy of the Hessian points", &
      [character(len=160) :: "    SCC accuracy of the displaced single points of the Hessian (default 0.3)"])
   call add("--lirc-restart MODE", "irc", &
      "wavefunction restart policy: off, pin or warm", &
      [character(len=160) :: "    wavefunction restart policy of the engine, off (default), pin or warm"], &
      values=v_restart)
   call add("--lvib-hessian SOURCE", "frequencies", &
      "Hessian backend: numerical, model or file", &
      [character(len=160) :: &
      "    Hessian backend of the frequency task, one of numerical (default, central", &
      "    differences of the electronic gradient), model (Lindh-D2) or file (read the", &
      "    Hessian of an earlier run)"], values=v_hessian)
   call add("--lvib-hessian-file NAME", "frequencies", &
      "Hessian file read by --lvib-hessian file", &
      [character(len=160) :: "    file read by --lvib-hessian file (TURBOMOLE Hessian layout)"])
   call add("--lvib-hess-step VALUE", "frequencies", &
      "FD step of the numerical Hessian (bohr)", &
      [character(len=160) :: &
      "    finite difference step of the numerical Hessian in bohr (default 0.005);", &
      "    the canonical name of --lvib-step"])
   call add("--lvib-step VALUE", "frequencies", &
      "alias of --lvib-hess-step", &
      [character(len=160) :: &
      "    alias of --lvib-hess-step: finite difference step of the numerical", &
      "    Hessian in bohr (default 0.005)"])
   call add("--lvib-sccacc VALUE", "frequencies", &
      "SCC accuracy of the FD Hessian points", &
      [character(len=160) :: "    SCC accuracy of the displaced single points of the numerical Hessian (default 0.3)"])
   call add("--lvib-restart MODE", "frequencies", &
      "wavefunction restart policy: off, pin or warm", &
      [character(len=160) :: &
      "    wavefunction restart policy of the engine, off (default), pin or warm.", &
      "    The finite difference sequence resets the wavefunction to the state it", &
      "    entered the Hessian with before every displacement, so warm and pin", &
      "    hand the provider the same start and produce a byte identical", &
      "    `hessian` file; only off is distinguishable (it records the cold start", &
      "    wavefunction).  The two spellings are kept, the equivalence is", &
      "    asserted by the xtb/landscape-lvib-restart-warm check"], &
      values=v_restart)
   call add("--lvib-odlr", "frequencies", &
      "request the ODLR approximation (not implemented)", &
      [character(len=160) :: "    request the ODLR approximation; not implemented, the run reports the request"])
   call add("--lvib-no-project", "frequencies", &
      "do not project translations and rotations", &
      [character(len=160) :: "    do not project the translations and rotations out of the Hessian"])
   call add("--lvib-temp VALUE", "frequencies", &
      "temperature of the thermodynamic table (K)", &
      [character(len=160) :: "    temperature of the thermodynamic table in K (default 298.15)"])
   call add("--lvib-rotor-cutoff VALUE", "frequencies", &
      "quasi-RRHO rotor cutoff (cm^-1)", &
      [character(len=160) :: "    rotor cutoff of the quasi-RRHO entropy treatment in cm^-1 (default 50.0)"])
   call add("--lvib-imag-cutoff VALUE", "frequencies", &
      "imaginary frequency cutoff (cm^-1)", &
      [character(len=160) :: "    imaginary frequency cutoff in cm^-1 (default -20.0, needs the minus sign)"])
   call add("--lvib-freq-scal VALUE", "frequencies", &
      "frequency scaling factor of the thermochemistry", &
      [character(len=160) :: "    frequency scaling factor of the thermochemistry (default 1.0)"])
   call add("--ceh <GRAD> [REAL]", "hamiltonian", &
      "CEH charges and numerical gradients", &
      [character(len=160) :: &
      "    calculate CEH (Charge-Extended Hückel model) charges and write them to the ceh.charges file", &
      "    optionally, calculate numerical gradients and write them to the ceh.charges.numgrad file", &
      "    with an adjustable step size for the numerical gradients."])
   call add("--gfn INT", "hamiltonian", &
      "GFN level: 0, 1, 2 or ff", &
      [character(len=160) :: &
      "    specify parametrisation of GFN-xTB: 0 (GFN0-xTB), 1 (GFN1-xTB), 2", &
      "    (GFN2-xTB, the default) or ff (GFN-FF, the same as --gfnff); any", &
      "    other value is refused with an error"], values=v_gfn)
   call add("--gfn0", "hamiltonian", &
      "GFN0-xTB Hamiltonian (a discouraged spelling)", &
      [character(len=160) :: &
      "    use the GFN0-xTB (EHT) Hamiltonian, the same as --gfn 0; a discouraged", &
      "    spelling, the run prints a warning and asks for --gfn 0"])
   call add("--gfn1", "hamiltonian", &
      "GFN1-xTB Hamiltonian (a discouraged spelling)", &
      [character(len=160) :: &
      "    use the GFN1-xTB Hamiltonian, the same as --gfn 1; a discouraged", &
      "    spelling, the run prints a warning and asks for --gfn 1"])
   call add("--gfn2", "hamiltonian", &
      "GFN2-xTB Hamiltonian (the default)", &
      [character(len=160) :: "    use the GFN2-xTB Hamiltonian (the default), the same as --gfn 2"])
   call add("--gfnff, --gff", "hamiltonian", &
      "GFN-FF force field", &
      [character(len=160) :: "    specify parametrisation of GFN-FF"])
   call add("--mcgfnff", "hamiltonian", &
      "MC-GFN-FF force field", &
      [character(len=160) :: "    use the MC-GFN-FF force field"])
   call add("--iff", "hamiltonian", &
      "IFF intermolecular force field", &
      [character(len=160) :: &
      "    use the IFF (intermolecular force field) calculator; the fragment", &
      "    sizes are taken from --nat"])
   call add("--nat LIST", "hamiltonian", &
      "atoms per IFF fragment", &
      [character(len=160) :: "    number of atoms per fragment of the IFF calculation, comma separated"])
   call add("--tblite,", "hamiltonian", &
      "use the tblite library", &
      [character(len=160) :: "    use tblite library as implementation for xTB, please note that solvation is not yet implemented within tblite"])
   call add("--dipro [REAL]", "hamiltonian", &
      "DIPRO coupling integrals (tblite only)", &
      [character(len=160) :: &
      "    calculate DIPRO coupling integrals (tblite builds only); the optional", &
      "    REAL is the coupling threshold (default 0.1)"])
   call add("--ptb,", "hamiltonian", &
      "density tight-binding method PTB", &
      [character(len=160) :: &
      "    performs single-point calculation with the density tight-binding method PTB.", &
      "    Provides electronic structure and properties, such as, e.g., atomic charges, bond orders, and dipole moments,", &
      "    but does not provide any energy-related properties, such as, e.g., total energy, nuclear gradients, or vibrational frequencies."])
   call add("--spinpol,", "hamiltonian", &
      "spin polarisation (tblite only)", &
      [character(len=160) :: "    enables spin-polarization for xTB methods (tblite required)"])
   call add("--vparam FILE", "hamiltonian", &
      "vTB parameter file", &
      [character(len=160) :: "    Parameter file for vTB calculation"])
   call add("-c, --chrg, --charge INT", "electronic", &
      "molecular charge", &
      [character(len=160) :: "    specify molecular charge as INT, overrides .CHRG file and xcontrol option"])
   call add("-u, --uhf INT", "electronic", &
      "number of unpaired electrons", &
      [character(len=160) :: &
      "    specify number of unpaired electrons as INT, overrides .UHF file and xcontrol option", &
      "    (default: 0 (singlet) for even electron number, 1 (doublet) for odd electron number)"])
   call add("--efield REAL,REAL,REAL", "electronic", &
      "static electric field in Cartesian coordinates", &
      [character(len=160) :: &
      "    static electric field in Cartesian coordinates, overrides '.EFIELD' file,", &
      "    works only via tblite for the xTB Hamiltonians, or with GFN-FF and PTB."])
   call add("--etemp REAL", "electronic", &
      "electronic temperature (default 300 K)", &
      [character(len=160) :: "    electronic temperature (default = 300K)"])
   call add("-a, --acc REAL", "electronic", &
      "SCC accuracy, lower is better", &
      [character(len=160) :: "    accuracy for SCC calculation, lower is better (default = 1.0)"])
   call add("--iterations, --maxiterations INT", "electronic", &
      "maximum SCC iterations per single point", &
      [character(len=160) :: &
      "    maximum number of SCC iterations per single point calculation", &
      "    (default = 250)"])
   call add("--tm", "external", &
      "external Turbomole energy and gradient", &
      [character(len=160) :: "    use the external Turbomole program for the energy and the gradient"])
   call add("--orca", "external", &
      "external ORCA energy and gradient", &
      [character(len=160) :: "    use the external ORCA program for the energy and the gradient"])
   call add("--mopac", "external", &
      "external MOPAC energy and gradient", &
      [character(len=160) :: "    use the external MOPAC program for the energy and the gradient"])
   call add("--driver [EXECUTABLE]", "external", &
      "external driver program", &
      [character(len=160) :: &
      "    use an external driver program for the energy and the gradient; the", &
      "    optional argument names the executable to run (there is no built-in default)"])
   call add("--oniom METHOD LIST", "external", &
      "subtractive embedding (ONIOM)", &
      [character(len=160) :: &
      "    use subtractive embedding via ONIOM method. 'METHOD' is given as 'inner:outer'", &
      "    where 'inner' can be 'orca', 'turbomole', 'gfn2', 'gfn1', or 'gfnff' and", &
      "    'outer' can be 'gfn2', 'gfn1', or 'gfnff'.", &
      "    The inner region is given as a comma separated indices directly in the commandline", &
      "    or in a file with each index on a separate line."])
   call add("--cut", "external", &
      "ONIOM: cut out the inner region only", &
      [character(len=160) :: &
      "    ONIOM: cut the inner region out of the input and stop without a", &
      "    calculation"])
   call add("--alpb SOLVENT [STATE]", "solvation", &
      "ALPB implicit solvation", &
      [character(len=160) :: &
      "    analytical linearized Poisson-Boltzmann (ALPB) model,", &
      "    available solvents are acetone, acetonitrile, aniline, benzaldehyde,", &
      "    benzene, ch2cl2, chcl3, cs2, dioxane, dmf, dmso, ethanol, ether,", &
      "    ethylacetate, furane, hexadecane, hexane, methanol, nitromethane, octanol,", &
      "    woctanol, phenol, toluene, thf, water.", &
      "    The solvent input is not case-sensitive. The Gsolv", &
      "    reference state can be chosen as reference, bar1M, or gsolv (default)."], &
      values=v_alpb)
   call add("-g, --gbsa SOLVENT [STATE]", "solvation", &
      "GBSA implicit solvation", &
      [character(len=160) :: &
      "    generalized born (GB) model with solvent accessable surface (SASA) model,", &
      "    available solvents are acetone, acetonitrile, benzene, CH2Cl2,", &
      "    CHCl3, CS2, DMF (only GFN2-xTB), DMSO, ether, H2O, methanol,", &
      "    n-hexane (only GFN2-xTB), THF and toluene.", &
      "    The solvent input is not case-sensitive.", &
      "    The Gsolv reference state can be chosen as reference, bar1M, or gsolv (default)."], &
      values=v_gbsa)
   call add("--gbe SOLVENT/EPSILON", "solvation", &
      "GBe solvation for finite epsilon (tblite)", &
      [character(len=160) :: &
      "    generalized Born for finite epsilon (GBe) solvation model (tblite required),", &
      "    includes only the electrostatic solvation contribution (no SASA model),", &
      "    available solvents are all solvents that are available for alpb.", &
      "    Additionally, the dielectric constant can be set manually."])
   call add("--gb SOLVENT/EPSILON", "solvation", &
      "GB solvation (tblite)", &
      [character(len=160) :: &
      "    generalized Born (GB) solvation model (tblite required),", &
      "    includes only the electrostatic solvation contribution (no SASA model),", &
      "    available solvents are all solvents that are available for alpb.", &
      "    Additionally, the dielectric constant can be set manually."])
   call add("--cosmo SOLVENT/EPSILON", "solvation", &
      "ddCOSMO solvation", &
      [character(len=160) :: &
      "    domain decomposition conductor-like screening model (ddCOSMO),", &
      "    available solvents are all solvents that are available for alpb.", &
      "    Additionally, the dielectric constant can be set manually or an ideal conductor", &
      "    can be chosen by setting epsilon to infinity."])
   call add("--tmcosmo SOLVENT/EPSILON", "solvation", &
      "ddCOSMO with the TM convention", &
      [character(len=160) :: "    same as --cosmo, but uses TM convention for writing the .cosmo files."])
   call add("--cpcmx SOLVENT", "solvation", &
      "CPCM-X solvation", &
      [character(len=160) :: &
      "    extended conduction-like polarizable continuum solvation model (CPCM-X),", &
      "    available solvents are all solvents included in the Minnesota Solvation Database."])
   call add("--cma", "solvation", &
      "centre-of-mass / principal-axis transform", &
      [character(len=160) :: &
      "    shifts molecule to center of mass and transforms cartesian coordinates into the", &
      "    coordinate system of the principle axis (not affected by ‘isotopes’-file)."])
   call add("--esp", "analysis", &
      "electrostatic potential on a VdW grid", &
      [character(len=160) :: "    calculate electrostatic potential on VdW-grid"])
   call add("--stm", "analysis", &
      "STM image", &
      [character(len=160) :: "    calculate STM image"])
   call add("--pop", "analysis", &
      "Mulliken population analysis", &
      [character(len=160) :: "    requests printout of Mulliken population analysis"])
   call add("--molden", "analysis", &
      "write a molden file", &
      [character(len=160) :: "    requests printout of molden file"])
   call add("--dipole", "analysis", &
      "dipole printout", &
      [character(len=160) :: "    requests dipole printout"])
   call add("--raman", "analysis", &
      "Raman spectrum (GFN2-xTB + PTB)", &
      [character(len=160) :: &
      "    requests Raman spectrum calculation via combination of GFN2-xTB and PTB", &
      "    using the temperature 'REAL' (default 298.15 K) and the wavelength of the", &
      "    incident laser which must be given in nm 'REAL' (default 514 nm)."])
   call add("--alpha", "analysis", &
      "static molecular dipole polarizabilities", &
      [character(len=160) :: "    requests the extension of electrical properties to static molecular dipole polarizabilities. "])
   call add("--wbo", "analysis", &
      "Wiberg bond order printout", &
      [character(len=160) :: "    requests Wiberg bond order printout"])
   call add("--lmo", "analysis", &
      "orbital localisation", &
      [character(len=160) :: "    requests localization of orbitals"])
   call add("--fod", "analysis", &
      "FOD calculation", &
      [character(len=160) :: "    requests FOD calculation"])
   call add("--enso", "analysis", &
      "machine-readable thermochemistry xtb_enso.json", &
      [character(len=160) :: &
      "    write the machine readable thermochemistry xtb_enso.json in a", &
      "    frequency run"])
   call add("--vip", "analysis", &
      "ionisation potential (GFN1)", &
      [character(len=160) :: &
      "    performs calculation of ionisation potential. This needs the param_ipea-xtb.txt", &
      "    parameters and a GFN1 Hamiltonian."])
   call add("--vea", "analysis", &
      "electron affinity (GFN1)", &
      [character(len=160) :: &
      "    performs calculation of electron affinity. This needs the param_ipea-xtb.txt", &
      "    parameters and a GFN1 Hamiltonian."])
   call add("--vipea", "analysis", &
      "IP and EA together (GFN1)", &
      [character(len=160) :: &
      "    performs calculation of electron affinity and ionisation potential.", &
      "    This needs the param_ipea-xtb.txt parameters and a GFN1 Hamiltonian."])
   call add("--vfukui", "analysis", &
      "Fukui indices", &
      [character(len=160) :: "    performs calculation of Fukui indices."])
   call add("--vomega", "analysis", &
      "electrophilicity index (GFN1)", &
      [character(len=160) :: &
      "    performs calculation of electrophilicity index. This needs the param_ipea-xtb.txt", &
      "    parameters and a GFN1 Hamiltonian."])
   call add("--color SCHEME", "output", &
      "colour scheme of the tblite output", &
      [character(len=160) :: &
      "    colour scheme of the tblite output, one of auto (follow the build", &
      "    feature, default), always or never (tblite builds only)"], values=v_color)
   call add("--namespace STRING", "output", &
      "namespace of every file of the run", &
      [character(len=160) :: &
      "    give this oxtb(1) run a namespace. All files, even temporary ones, will be named", &
      "    according to STRING (might not work everywhere)."])
   call add("--json", "output", &
      "write the xtbout.json file", &
      [character(len=160) :: "    write xtbout.json file"])
   call add("--help-values", "output", &
      "machine-readable dump of the option value sets", &
      [character(len=160) :: &
      "    print the legal value set of every documented option as a stable,", &
      "    line-oriented table, one record per option in --help order, with a", &
      "    TAB between the three fields", &
      "      category <spellings> <values>", &
      "    `spellings` is the spelling line of --help (a token that starts", &
      "    with a dash is an option spelling, any other token is its argument", &
      "    placeholder) and `values` is the space separated legal value set,", &
      "    empty when the option has none; lines starting with `#` are", &
      "    comments and fix this format (format 1).  A completion generator", &
      "    codes against this table"])
   call add("--wrtopo LIST", "output", &
      "comma separated topology printouts", &
      [character(len=160) :: &
      "    comma separated topology printouts, from etot, gnorm, nb, bpair, alist,", &
      "    blist, tlist, vtors, vbond, vangl, hbbond, eeq and json (json selects", &
      "    all of them); a misspelled keyword aborts the run"], values=v_wrtopo)
   call add("--citation", "output", &
      "print citation and terminate", &
      [character(len=160) :: "    print citation and terminate"])
   call add("--license", "output", &
      "print license and terminate", &
      [character(len=160) :: "    print license and terminate"])
   call add("--version", "output", &
      "print the version banner and terminate", &
      [character(len=160) :: "    print the version banner and terminate"])
   call add("-v, --verbose", "output", &
      "be more verbose", &
      [character(len=160) :: "    be more verbose (not supported in every unit)"])
   call add("-V, --very-verbose", "output", &
      "be very verbose", &
      [character(len=160) :: "    be very verbose (this also enables the verbose output)"])
   call add("--ceasefiles", "output", &
      "reduce output and files written", &
      [character(len=160) :: "    reduce the amount of output and files written"])
   call add("--strict", "output", &
      "turn all warnings into hard errors", &
      [character(len=160) :: "    turns all warnings into hard errors"])
   call add("-h, --help", "output", &
      "show this help page", &
      [character(len=160) :: "    show help page"])
   call add("--coffee", "workflow", &
      "run the built-in caffeine molecule", &
      [character(len=160) :: &
      "    ignore the input geometry and run the built-in caffeine molecule", &
      "    (a debugging helper)"])
   call add("-I, --input FILE", "workflow", &
      "xcontrol(7) instruction file", &
      [character(len=160) :: "    use FILE as input source for xcontrol(7) instructions"])
   call add("--[no]copy", "workflow", &
      "copy the xcontrol file at startup", &
      [character(len=160) :: "    copies the xcontrol file at startup (default = true)"])
   call add("--[no]restart", "workflow", &
      "restart from xtbrestart", &
      [character(len=160) :: "    restarts calculation from xtbrestart (default = true)"])
   call add("--define", "workflow", &
      "check the input and terminate", &
      [character(len=160) :: "    performs automatic check of input and terminate"])
   call add("--nocellopt", "workflow", &
      "keep a periodic cell fixed", &
      [character(len=160) :: &
      "    do not optimise the lattice of a periodic system, the cell is kept", &
      "    fixed (default: the cell is optimised)"])
   call add("-P, --parallel INT", "performance", &
      "number of OpenMP threads of the run", &
      [character(len=160) :: &
      "    number of OpenMP threads of the run, applied before the first", &
      "    parallel region so that it wins over OMP_NUM_THREADS; INT must be a", &
      "    positive integer - any other value is reported (and fatal with", &
      "    --strict) - and a value above the processor count is reported as", &
      "    well; the value itself is applied as given, never clamped"])
   call add("--blas-threads N", "performance", &
      "thread count of the linear algebra backend", &
      [character(len=160) :: &
      "    set the thread count of the linear algebra backend (OpenBLAS or MKL)", &
      "    at run time and print the value read back; without the flag the", &
      "    backend is left untouched (default)"])

   if (present(values_only)) then
      if (values_only) then
         call write_values_dump(iunit)
         status = 0
         return
      end if
   end if

   if (len_trim(topic) == 0) then
      call write_index(iunit, use_color)
      status = 0
      return
   end if

   norm = strip_dashes(topic)
   norm_lower = to_lower(norm)
   do icat = 1, n_cat
      if (norm_lower == cat_key(icat)) then
         call write_category(iunit, icat, use_color)
         status = 0
         return
      end if
   end do
   do i = 1, n
      if (index(opts(i)%aliases_exact, "/" // norm // "/") > 0) then
         call write_option(iunit, i, use_color)
         status = 0
         return
      end if
   end do
   do i = 1, n
      if (index(opts(i)%aliases, "/" // norm_lower // "/") > 0) then
         call write_option(iunit, i, use_color)
         status = 0
         return
      end if
   end do

   call write_unknown(iunit, topic, use_color)
   status = 1

contains

   !> Whether the page should be coloured.  `CLICOLOR_FORCE` set to a value
   !> other than `0` forces colour even into a pipe and wins over everything;
   !> otherwise `NO_COLOR` set and non-empty disables it; otherwise colour
   !> follows the output terminal (`color_support`, the same no-color.org
   !> precedence the tail of the page and `man 1 oxtb` document).
   function colour_wanted() result(use_color)
      logical :: use_color
      character(len=:), allocatable :: value
      logical :: is_set, force

      call read_env("CLICOLOR_FORCE", value, is_set)
      force = is_set .and. trim(value) /= "0"
      if (force) then
         use_color = .true.
      else
         call read_env("NO_COLOR", value, is_set)
         if (is_set) then
            use_color = .false.
         else
            use_color = get_xtb_feature("color")
         end if
      end if
   end function colour_wanted

   !> Read an environment variable.  `is_set` is false when it is unset or
   !> empty, which is what the colour conventions mean by "set".
   subroutine read_env(name, value, is_set)
      character(len=*), intent(in) :: name
      character(len=:), allocatable, intent(out) :: value
      logical, intent(out) :: is_set
      integer :: length, stat

      call get_environment_variable(name, length=length, status=stat)
      is_set = stat == 0 .and. length > 0
      if (is_set) then
         allocate(character(len=length) :: value)
         call get_environment_variable(name, value)
      else
         value = ""
      end if
   end subroutine read_env

   !> Wrap `text` in one SGR run when colour is on and return it unchanged
   !> otherwise, so the escape sequences are purely additive and stripping
   !> them reproduces the plain page byte for byte.
   function styled(text, sgr, use_color) result(res)
      character(len=*), intent(in) :: text, sgr
      logical, intent(in) :: use_color
      character(len=:), allocatable :: res

      if (use_color) then
         res = sgr // text // sgr_reset
      else
         res = text
      end if
   end function styled

   !> Append one option to opts; called once per documented option.
   subroutine add(spellings, category, summary, detail, values)
      character(len=*), intent(in) :: spellings, category, summary
      character(len=*), intent(in) :: detail(:)
      !> Legal values of the option; absent for a flag or a free value
      character(len=16), intent(in), optional :: values(:)
      n = n + 1
      opts(n)%spellings = spellings
      opts(n)%category = category
      opts(n)%summary = summary
      opts(n)%detail = detail
      if (present(values)) opts(n)%values = values
      call append_aliases(opts(n)%aliases, opts(n)%aliases_exact, spellings)
   end subroutine add

   !> Strip the leading dashes and the surrounding blanks of a topic.
   function strip_dashes(text) result(res)
      character(len=*), intent(in) :: text
      character(len=:), allocatable :: res
      res = trim(adjustl(text))
      do while (len(res) > 0 .and. res(1:1) == '-')
         res = res(2:)
      end do
   end function strip_dashes

   !> Lowercase a string; Fortran has no intrinsic for this.
   function to_lower(text) result(res)
      character(len=*), intent(in) :: text
      character(len=:), allocatable :: res
      integer :: i, code
      res = text
      do i = 1, len(res)
         code = iachar(res(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) &
            res(i:i) = achar(code + 32)
      end do
   end function to_lower

   !> Collect the `/key/` lookup keys of every spelling of `spellings`.
   subroutine append_aliases(aliases, aliases_exact, spellings)
      character(len=:), allocatable, intent(inout) :: aliases, aliases_exact
      character(len=*), intent(in) :: spellings
      character(len=:), allocatable :: token
      integer :: i
      aliases = ""
      aliases_exact = ""
      token = ""
      do i = 1, len(spellings)
         if (spellings(i:i) == ',' .or. spellings(i:i) == ' ') then
            call append_alias(aliases, aliases_exact, token)
            token = ""
         else
            token = token // spellings(i:i)
         end if
      end do
      call append_alias(aliases, aliases_exact, token)
   end subroutine append_aliases

   !> Append the key of one token; `--[no]x` yields `x` and `nox`.
   subroutine append_alias(aliases, aliases_exact, token)
      character(len=:), allocatable, intent(inout) :: aliases, aliases_exact
      character(len=*), intent(in) :: token
      character(len=:), allocatable :: key, key_exact
      integer :: pos
      if (len(token) == 0) return
      if (token(1:1) /= '-') return
      key_exact = strip_dashes(token)
      key = to_lower(key_exact)
      pos = index(key, '[no]')
      if (pos > 0) then
         aliases = aliases // "/" // key(1:pos-1) // key(pos+4:) // "/"
         aliases = aliases // "/" // key(1:pos-1) // "no" // key(pos+4:) // "/"
         aliases = aliases // "/" // key // "/"
         aliases_exact = aliases_exact // "/" // key_exact(1:pos-1) // &
            key_exact(pos+4:) // "/"
         aliases_exact = aliases_exact // "/" // key_exact(1:pos-1) // &
            "no" // key_exact(pos+4:) // "/"
         aliases_exact = aliases_exact // "/" // key_exact // "/"
      else
         aliases = aliases // "/" // key // "/"
         aliases_exact = aliases_exact // "/" // key_exact // "/"
      end if
   end subroutine append_alias

   !> The grouped index: preamble, one line per option under its category,
   !> then the machine settings, output conventions and man pointer.
   subroutine write_index(iunit, use_color)
      integer, intent(in) :: iunit
      logical, intent(in) :: use_color
      integer :: icat, i
      write(iunit,'(a)') &
         "Usage: oxtb [OPTIONS] <GEOMETRY> [OPTIONS]", &
         "", &
         "<GEOMETRY> may be provided as valid TM coordinate file (*coord in Bohr),", &
         "in xmol format (*.xyz in Ångström), sdf or mol file format, PDB format", &
         "genFormat input or Vasp's POSCAR format.", &
         "", &
         "For a full list visit: https://grimme-lab.github.io/mctc-lib/page/index.html", &
         "", &
         "Full text of one option: oxtb --help <OPTION>  (e.g. oxtb --help --method)", &
         "One category of the index: oxtb --help <CATEGORY>; the categories are:", &
         "runtypes, suite, levels, minimisation, scan, saddle, band, irc,", &
         "frequencies, hamiltonian, electronic, external, solvation, analysis,", &
         "output, workflow, performance", &
         ""
      do icat = 1, n_cat
         write(iunit,'(a)') ''
         write(iunit,'(a)') styled(trim(cat_title(icat)), sgr_heading, use_color)
         do i = 1, n
            if (opts(i)%category == cat_key(icat)) &
               call write_entry(iunit, i, use_color)
         end do
      end do
      call write_tail(iunit)
   end subroutine write_index

   !> One category of the index, the heading and its option lines.
   subroutine write_category(iunit, icat, use_color)
      integer, intent(in) :: iunit, icat
      logical, intent(in) :: use_color
      integer :: i
      write(iunit,'(a)') styled(trim(cat_title(icat)), sgr_heading, use_color)
      do i = 1, n
         if (opts(i)%category == cat_key(icat)) &
            call write_entry(iunit, i, use_color)
      end do
   end subroutine write_category

   !> One index line: the spellings and the one line gloss.
   subroutine write_entry(iunit, i, use_color)
      integer, intent(in) :: iunit, i
      logical, intent(in) :: use_color
      write(iunit,'(a)') "  " // &
         styled(trim(opts(i)%spellings), sgr_spelling, use_color) // "  " // &
         trim(opts(i)%summary)
   end subroutine write_entry

   !> The full text of one option: spellings, every detail line, the uniform
   !> `Values:` line when the option has a value set, blank.
   subroutine write_option(iunit, i, use_color)
      integer, intent(in) :: iunit, i
      logical, intent(in) :: use_color
      integer :: j
      write(iunit,'(a)') styled(trim(opts(i)%spellings), sgr_spelling, use_color)
      do j = 1, size(opts(i)%detail)
         write(iunit,'(a)') trim(opts(i)%detail(j))
      end do
      if (allocated(opts(i)%values)) then
         if (size(opts(i)%values) > 0) then
            write(iunit,'(a)') "Values: " // join_values(opts(i)%values, ", ")
         end if
      end if
      write(iunit,'(a)') ''
   end subroutine write_option

   !> Join a value set with `separator`; an empty set gives an empty string.
   pure function join_values(values, separator) result(res)
      !> Value set
      character(len=16), allocatable, intent(in) :: values(:)
      !> Separator between two values
      character(len=*), intent(in) :: separator
      !> Joined list
      character(len=:), allocatable :: res
      integer :: k
      res = ""
      if (.not. allocated(values)) return
      do k = 1, size(values)
         if (k > 1) res = res // separator
         res = res // trim(values(k))
      end do
   end function join_values

   !> The `--help-values` table: one tab separated record per documented
   !> option, in `--help` order, plus `#` comment lines that fix the format so
   !> a completion generator never has to guess.
   subroutine write_values_dump(iunit)
      integer, intent(in) :: iunit
      integer :: i, n_valued
      n_valued = 0
      do i = 1, n
         if (allocated(opts(i)%values)) then
            if (size(opts(i)%values) > 0) n_valued = n_valued + 1
         end if
      end do
      write(iunit,'(a)') &
         "# oxtb --help-values format 1", &
         "# one record per documented option, in --help order; one TAB between", &
         "# the three fields:", &
         "#   field 1 category, field 2 spellings, field 3 values", &
         "# field 2 is the spelling line of `oxtb --help`: comma separated", &
         "# tokens, a token that starts with a dash is an option spelling and", &
         "# any other token is its argument placeholder; field 3 is the space", &
         "# separated legal value set, empty when the option has none.", &
         "# lines starting with # are comments"
      write(iunit,'(a,i0,a,i0)') "# options: ", n, " with values: ", n_valued
      do i = 1, n
         write(iunit,'(a,a,a,a,a)') trim(opts(i)%category), achar(9), &
            trim(opts(i)%spellings), achar(9), join_values(opts(i)%values, " ")
      end do
   end subroutine write_values_dump

   !> The closing prose: machine settings, output conventions, man pointer.
   subroutine write_tail(iunit)
      integer, intent(in) :: iunit
      write(iunit,'(a)') &
         "", &
         "Useful Maschine Settings: ", &
         "", &
         "export OMP_NUM_THREADS=<NCORE>,1", &
         "export OMP_STACKSIZE=4G", &
         "ulimit -s unlimited", &
         "", &
         "A threaded BLAS is the largest performance trap on multi core machines:", &
         "the OpenMP threads and the BLAS thread pool multiply.  On a 24 atom", &
         "anthracene GFN2 single point at OMP_NUM_THREADS=8 the fastest setting", &
         "measured was OPENBLAS_NUM_THREADS=1 (12 ms), ahead of --blas-threads 1", &
         "(99 ms) and of the untouched pool (583 ms); the same order held at 501", &
         "atoms.  OPENBLAS_NUM_THREADS is read before main, so the shell must", &
         "export it, while --blas-threads resizes the existing pool from inside", &
         "the run.  The -P and --parallel flags set the OpenMP thread count;", &
         "where both reach the same knob the flag that appears last wins.", &
         "See man 1 oxtb for the full table.", &
         "", &
         "Colour:", &
         "", &
         "this page is coloured on a terminal.  CLICOLOR_FORCE set to a value", &
         "other than 0 forces colour even into a pipe and wins over everything", &
         "else; otherwise NO_COLOR set to a non-empty value disables it; otherwise", &
         "colour follows the output terminal, so a piped or redirected page", &
         "carries no escape sequence and is the plain page.", &
         "", &
         "Output Conventions:", &
         "", &
         "total energies are given in atomic units (Eh)", &
         "gaps/HL energies are given in eV", &
         "", &
         "More information can be obtained by `man 1 oxtb` and `man 7 xcontrol`", &
         "or at https://xtb-docs.readthedocs.io/en/latest/", &
         ""
   end subroutine write_tail

   !> The error of an unknown topic: it names the categories and the usage.
   subroutine write_unknown(iunit, topic, use_color)
      integer, intent(in) :: iunit
      character(len=*), intent(in) :: topic
      logical, intent(in) :: use_color
      integer :: icat
      write(iunit,'(a)') &
         "", &
         "Unknown help topic '" // trim(topic) // "'.", &
         "", &
         styled("Usage:", sgr_usage, use_color) // " oxtb --help [TOPIC]", &
         "  " // styled("oxtb --help", sgr_usage, use_color) // &
            "              grouped index of every option", &
         "  " // styled("oxtb --help", sgr_usage, use_color) // &
            " <OPTION>     full text of one option, e.g. oxtb --help --method", &
         "  " // styled("oxtb --help", sgr_usage, use_color) // &
            " <CATEGORY>   one category of the index", &
         "", &
         "Categories:"
      do icat = 1, n_cat
         write(iunit,'(a)') "  " // trim(cat_key(icat)) // "  " // &
            trim(cat_title(icat))
      end do
   end subroutine write_unknown

end subroutine help

subroutine definebanner
   use, intrinsic :: iso_fortran_env, only : id => output_unit
   write(id,"(""           _"")")
   write(id,"(""          | |"")")
   write(id,"(""  _ __ ___| |_ _   _ _ __ _ __        __"")")
   write(id,"("" | '__/ _ \ __| | | | '__| '_ \      / _|"")")
   write(id,"("" | | |  __/ |_| |_| | |  | | | |___ | |_"")")
   write(id,"("" |_|  \___|\__|\__,_|_|  |_| |_/ _ \|  _|"")")
   write(id,"(""         _       __ _         | (_) | |"")")
   write(id,"(""        | |     / _(_)         \___/|_|"")")
   write(id,"(""      __| | ___| |_ _ _ __   ___"")")
   write(id,"(""     / _` |/ _ \  _| | '_ \ / _ \"")")
   write(id,"(""    | (_| |  __/ | | | | | |  __/"")")
   write(id,"(""     \__,_|\___|_| |_|_| |_|\___|"")")
end subroutine definebanner
