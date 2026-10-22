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

!> Unit tests of the landscape settings and the four-criterion convergence test.
module test_landscape_settings
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, &
      & ieee_positive_inf
   use testdrive, only : new_unittest, unittest_type, error_type, check
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : aatoau
   use xtb_landscape_settings
   implicit none
   private

   public :: collect_landscape_settings

contains

!> Collect all exported unit tests
subroutine collect_landscape_settings(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("defaults-common", test_defaults_common), &
      new_unittest("defaults-opt", test_defaults_opt), &
      new_unittest("defaults-sp", test_defaults_sp), &
      new_unittest("convergence-boundary", test_convergence_boundary), &
      new_unittest("convergence-measures", test_convergence_measures), &
      new_unittest("convergence-levels", test_convergence_levels), &
      new_unittest("level-step-budget", test_level_step_budget), &
      new_unittest("identifier-names", test_identifier_names), &
      new_unittest("algorithm-lookup", test_algorithm_lookup), &
      new_unittest("level-lookup", test_level_lookup), &
      new_unittest("enum-case", test_enum_case), &
      new_unittest("level-list", test_level_list), &
      new_unittest("number-strings", test_number_strings), &
      new_unittest("restart-names", test_restart_names) &
      ]

end subroutine collect_landscape_settings

!> Thresholds inherited by every landscape task
subroutine test_defaults_common(error)
   type(error_type), allocatable, intent(out) :: error
   type(landscape_opt_settings) :: set

   call check(error, set%max_steps, 256)
   if (allocated(error)) return
   call check(error, set%e_tol, 5.0e-6_wp)
   if (allocated(error)) return
   call check(error, set%g_max, 3.0e-4_wp)
   if (allocated(error)) return
   call check(error, set%g_rms, 1.0e-4_wp)
   if (allocated(error)) return
   call check(error, set%d_max, 4.0e-3_wp)
   if (allocated(error)) return
   call check(error, set%d_rms, 2.0e-3_wp)
   if (allocated(error)) return
   call check(error, set%printlevel, 1)

end subroutine test_defaults_common

!> Defaults of the geometry minimisation settings
subroutine test_defaults_opt(error)
   type(error_type), allocatable, intent(out) :: error
   type(landscape_opt_settings) :: set

   call check(error, p_landscape_opt_sd, 1)
   if (allocated(error)) return
   call check(error, p_landscape_opt_lbfgs, 2)
   if (allocated(error)) return
   call check(error, p_landscape_opt_bfgs, 3)
   if (allocated(error)) return
   call check(error, p_landscape_opt_diis, 4)
   if (allocated(error)) return
   call check(error, p_landscape_opt_rfo, 5)
   if (allocated(error)) return
   call check(error, p_landscape_opt_cg, 6)
   if (allocated(error)) return
   call check(error, p_landscape_hess_update_none, 0)
   if (allocated(error)) return
   call check(error, p_landscape_hess_update_bfgs, 1)
   if (allocated(error)) return
   call check(error, p_landscape_hess_update_sr1, 2)
   if (allocated(error)) return
   call check(error, p_landscape_hess_update_bofill, 3)
   if (allocated(error)) return
   call check(error, p_landscape_hess_update_powell, 4)
   if (allocated(error)) return

   call check(error, set%algorithm, p_landscape_opt_lbfgs)
   if (allocated(error)) return
   call check(error, set%lbfgs_memory, 5)
   if (allocated(error)) return
   call check(error, set%hessian_update, p_landscape_hess_update_bofill)
   if (allocated(error)) return
   call check(error, set%max_displ, 3.7795e-1_wp)
   if (allocated(error)) return
   call check(error, set%trust_radius, 3.7795e-1_wp)
   if (allocated(error)) return
   call check(error, set%trust_min, 1.0e-3_wp*aatoau)
   if (allocated(error)) return
   call check(error, set%trust_max, 1.0_wp*aatoau)
   if (allocated(error)) return
   call check(error, set%curvature, 70.0_wp)
   if (allocated(error)) return
   call check(error, set%hlow, 1.0e-2_wp)
   if (allocated(error)) return
   call check(error, set%line_search)
   if (allocated(error)) return
   call check(error, set%evals_eps, 1.0e-10_wp)
   if (allocated(error)) return
   call check(error, set%mu_margin, 1.0e-8_wp)
   if (allocated(error)) return
   call check(error, set%max_bisect_it, 60)
   if (allocated(error)) return
   ! recalc = 5: a numerical Hessian costs 6N single points, so recomputing it
   ! on every RFO step (`recalc = 1`) would be prohibitive
   call check(error, set%recalc, 5)
   if (allocated(error)) return
   call check(error, set%eta_shrink, 0.75_wp)
   if (allocated(error)) return
   call check(error, set%eta_expand, 1.75_wp)
   if (allocated(error)) return
   call check(error, set%armijo_c1, 1.0e-4_wp)
   if (allocated(error)) return
   call check(error, set%max_backtrack, 10)
   if (allocated(error)) return
   call check(error, set%diis_memory, 6)
   if (allocated(error)) return
   call check(error, set%diis_store_every, 5)
   if (allocated(error)) return
   call check(error, set%diis_min_snapshots, 3)
   if (allocated(error)) return
   call check(error, set%diis_regularization, 1.0e-10_wp)
   if (allocated(error)) return
   call check(error, set%write_files)
   if (allocated(error)) return
   call check(error, set%output_prefix, "xtb")
   if (allocated(error)) return
   call check(error, set%hessian_source, "numerical")
   if (allocated(error)) return
   call check(error, set%task(), "opt")
   if (allocated(error)) return

   call check(error, set%engine%etemp, 300.0_wp)
   if (allocated(error)) return
   call check(error, set%engine%maxiter, 250)
   if (allocated(error)) return
   call check(error, set%engine%acc, 1.0_wp)
   if (allocated(error)) return
   call check(error, set%engine%restart_mode, p_landscape_restart_off)
   if (allocated(error)) return
   call check(error, set%engine%hess_step, 5.0e-3_wp)
   if (allocated(error)) return
   call check(error, set%engine%hess_acc, 0.3_wp)
   if (allocated(error)) return
   call check(error, set%engine%odlr, .false.)
   if (allocated(error)) return

end subroutine test_defaults_opt

!> Defaults of the single point settings
subroutine test_defaults_sp(error)
   type(error_type), allocatable, intent(out) :: error
   type(landscape_sp_settings) :: set

   call check(error, set%write_gradient, .false.)
   if (allocated(error)) return
   call check(error, set%grad_threshold, 0.0_wp)
   if (allocated(error)) return
   call check(error, set%write_hessian, .false.)
   if (allocated(error)) return
   call check(error, set%max_steps, 256)
   if (allocated(error)) return
   call check(error, set%hessian_source, "numerical")
   if (allocated(error)) return
   call check(error, set%task(), "sp")

end subroutine test_defaults_sp

!> Boundary behaviour of the four-criterion convergence test.
!>
!> A hundred component vector is used so that a single component can sit
!> exactly on the max threshold without violating the RMS threshold.  The RMS
!> boundaries are probed by setting the threshold to the measured value, which
!> makes the `<=` comparison bit-exact.
subroutine test_convergence_boundary(error)
   type(error_type), allocatable, intent(out) :: error
   type(landscape_opt_settings) :: set
   real(wp) :: g(100), d(100)
   real(wp) :: measured

   g = 0.0_wp
   d = 0.0_wp
   g(1) = 3.0e-4_wp
   d(1) = 4.0e-3_wp

   ! exactly on the max thresholds and no energy change
   call check(error, set%is_converged(0.0_wp, g, d))
   if (allocated(error)) return

   ! g_max slightly above the threshold
   g(1) = 3.0001e-4_wp
   call check(error, .not. set%is_converged(0.0_wp, g, d))
   if (allocated(error)) return
   g(1) = 3.0e-4_wp

   ! d_max slightly above the threshold
   d(1) = 4.0001e-3_wp
   call check(error, .not. set%is_converged(0.0_wp, g, d))
   if (allocated(error)) return
   d(1) = 4.0e-3_wp

   ! g_rms exactly on a threshold that equals the measured RMS
   d = 0.0_wp
   g = 0.0_wp
   g(1) = 1.0e-4_wp
   measured = sqrt(sum(g**2)/real(size(g), wp))
   set%g_rms = measured
   call check(error, set%is_converged(0.0_wp, g, d))
   if (allocated(error)) return
   set%g_rms = measured - 1.0e-12_wp*measured
   call check(error, .not. set%is_converged(0.0_wp, g, d))
   if (allocated(error)) return

   ! d_rms exactly on a threshold that equals the measured RMS
   g = 0.0_wp
   d = 0.0_wp
   d(1) = 2.0e-3_wp
   measured = sqrt(sum(d**2)/real(size(d), wp))
   set%d_rms = measured
   call check(error, set%is_converged(0.0_wp, g, d))
   if (allocated(error)) return
   set%d_rms = measured - 1.0e-12_wp*measured
   call check(error, .not. set%is_converged(0.0_wp, g, d))
   if (allocated(error)) return

   ! energy guard
   d = 0.0_wp
   call check(error, set%is_converged(5.0e-6_wp, g, d))
   if (allocated(error)) return
   call check(error, .not. set%is_converged(5.0001e-6_wp, g, d))

end subroutine test_convergence_boundary

!> The RMS values are normalised by the number of vector components (3N)
subroutine test_convergence_measures(error)
   type(error_type), allocatable, intent(out) :: error
   real(wp) :: g(100), d(100)
   real(wp) :: g_max, g_rms, d_max, d_rms

   g = 0.0_wp
   d = 0.0_wp
   g(1) = 3.0e-4_wp
   d(1) = 4.0e-3_wp

   call landscape_convergence_measures(g, d, g_max, g_rms, d_max, d_rms)

   ! relative tolerances: an absolute 1e-18 threshold is below the rounding
   ! error of sqrt() on a 1e-5 quantity (~1e-21) only by chance.
   call check(error, g_max, 3.0e-4_wp, thr=1.0e-13_wp*3.0e-4_wp)
   if (allocated(error)) return
   call check(error, g_rms, 3.0e-5_wp, thr=1.0e-13_wp*3.0e-5_wp)
   if (allocated(error)) return
   call check(error, d_max, 4.0e-3_wp, thr=1.0e-13_wp*4.0e-3_wp)
   if (allocated(error)) return
   call check(error, d_rms, 4.0e-4_wp, thr=1.0e-13_wp*4.0e-4_wp)

end subroutine test_convergence_measures

!> Named convergence levels overwrite the four displacement/force thresholds
!> A named convergence level drives four axes (CLI redesign section 3.5, R5):
!> the four thresholds, `e_tol`, the SCC accuracy of the engine and the
!> multiplier of the step budget.  The values are the ones of the matching
!> upstream level (`get_optthr`, `src/optimizer.f90:48-119`), so the bridge of
!> `--opt <LEVEL>` is an identity on the numbers.
subroutine test_convergence_levels(error)
   type(error_type), allocatable, intent(out) :: error
   type(landscape_opt_settings) :: set

   call set%set_conv_level("loose")
   call check(error, set%g_max, 2.0e-3_wp)
   if (allocated(error)) return
   call check(error, set%g_rms, 5.0e-4_wp)
   if (allocated(error)) return
   call check(error, set%d_max, 1.0e-2_wp)
   if (allocated(error)) return
   call check(error, set%d_rms, 7.0e-3_wp)
   if (allocated(error)) return
   call check(error, set%e_tol, 5.0e-5_wp)
   if (allocated(error)) return
   call check(error, set%engine%acc, 2.0_wp)
   if (allocated(error)) return
   call check(error, set%max_steps_multiplier, 2)
   if (allocated(error)) return

   call set%set_conv_level("tight")
   call check(error, set%g_max, 1.0e-4_wp)
   if (allocated(error)) return
   call check(error, set%g_rms, 3.0e-5_wp)
   if (allocated(error)) return
   call check(error, set%d_max, 1.0e-3_wp)
   if (allocated(error)) return
   call check(error, set%d_rms, 6.0e-4_wp)
   if (allocated(error)) return
   call check(error, set%e_tol, 1.0e-6_wp)
   if (allocated(error)) return
   call check(error, set%engine%acc, 0.2_wp)
   if (allocated(error)) return
   call check(error, set%max_steps_multiplier, 5)
   if (allocated(error)) return

   call set%set_conv_level("normal")
   call check(error, set%g_max, 3.0e-4_wp)
   if (allocated(error)) return
   call check(error, set%g_rms, 1.0e-4_wp)
   if (allocated(error)) return
   call check(error, set%d_max, 4.0e-3_wp)
   if (allocated(error)) return
   call check(error, set%d_rms, 2.0e-3_wp)
   if (allocated(error)) return
   call check(error, set%e_tol, 5.0e-6_wp)
   if (allocated(error)) return
   call check(error, set%engine%acc, 1.0_wp)
   if (allocated(error)) return
   call check(error, set%max_steps_multiplier, 3)
   if (allocated(error)) return

   call set%set_conv_level("superloose")
   call check(error, set%g_max, 4.6e-3_wp)
   if (allocated(error)) return
   call check(error, set%e_tol, 5.0e-4_wp)
   if (allocated(error)) return
   call check(error, set%engine%acc, 3.0_wp)
   if (allocated(error)) return
   call check(error, set%max_steps_multiplier, 1)
   if (allocated(error)) return

   call set%set_conv_level("vtight")
   call check(error, set%g_max, 3.0e-5_wp)
   if (allocated(error)) return
   call check(error, set%e_tol, 1.0e-7_wp)
   if (allocated(error)) return
   call check(error, set%engine%acc, 0.05_wp)
   if (allocated(error)) return
   call check(error, set%max_steps_multiplier, 20)

end subroutine test_convergence_levels

!> The step budget a level derives from the atom count (CLI redesign section
!> 3.5, R5).
!>
!> The rule is `clamp(n_atoms * multiplier, 200, 10000)` of the upstream
!> `get_optthr`; the default level is `normal` (multiplier 3), so a water
!> molecule gets 200 steps - **fewer** than the fixed default of 256 this code
!> used before the level became the source of the axis.  An explicit budget
!> (`max_steps_source`, written by `--l*-max-steps` and by the `--cycles`
!> bridge) wins over the derived one.
subroutine test_level_step_budget(error)
   type(error_type), allocatable, intent(out) :: error
   type(landscape_opt_settings) :: set

   ! default level: `normal`, multiplier 3
   call set%level_step_budget(3)
   call check(error, set%max_steps, 200)
   if (allocated(error)) return
   ! a 500 atom system grows with the atom count
   call set%level_step_budget(500)
   call check(error, set%max_steps, 1500)
   if (allocated(error)) return
   ! and the clamp of the upstream rule holds on both ends
   call set%level_step_budget(20000)
   call check(error, set%max_steps, 10000)
   if (allocated(error)) return

   call set%set_conv_level("vtight")
   call set%level_step_budget(10)
   call check(error, set%max_steps, 200)
   if (allocated(error)) return
   call set%level_step_budget(2000)
   call check(error, set%max_steps, 10000)
   if (allocated(error)) return

   ! an explicit budget is not touched by the derivation, and it keeps its
   ! value when the level is named afterwards (`--lopt-max-steps 8
   ! --lopt-level tight` and the reverse order are the same run)
   set%max_steps = 8
   set%max_steps_source = '--lopt-max-steps'
   call set%set_conv_level("tight")
   call set%level_step_budget(1000)
   call check(error, set%max_steps, 8)

end subroutine test_level_step_budget

!> Human readable names of the algorithm and update identifiers
subroutine test_identifier_names(error)
   type(error_type), allocatable, intent(out) :: error

   call check(error, landscape_opt_algorithm_name(p_landscape_opt_sd), "sd")
   if (allocated(error)) return
   call check(error, landscape_opt_algorithm_name(p_landscape_opt_cg), "cg")
   if (allocated(error)) return
   call check(error, landscape_opt_algorithm_name(p_landscape_opt_lbfgs), "lbfgs")
   if (allocated(error)) return
   call check(error, landscape_opt_algorithm_name(p_landscape_opt_rfo), "rfo")
   if (allocated(error)) return
   call check(error, landscape_hessian_update_name(p_landscape_hess_update_bofill), &
      & "bofill")
   if (allocated(error)) return
   call check(error, landscape_hessian_update_name(p_landscape_hess_update_none), "none")

end subroutine test_identifier_names

!> Every algorithm name must map back to its identifier, unknown names to 0
subroutine test_algorithm_lookup(error)
   type(error_type), allocatable, intent(out) :: error
   integer :: algorithm

   do algorithm = 1, p_landscape_opt_cg
      call check(error, landscape_opt_algorithm_from_name( &
         & landscape_opt_algorithm_name(algorithm)), algorithm)
      if (allocated(error)) return
   end do
   call check(error, landscape_opt_algorithm_from_name("unknown"), 0)
   if (allocated(error)) return
   ! names are matched case insensitively, see `enum-case`
   call check(error, landscape_opt_algorithm_from_name("LBFGS"), &
      & p_landscape_opt_lbfgs)

end subroutine test_algorithm_lookup

!> Only the five ORCA levels are recognised
subroutine test_level_lookup(error)
   type(error_type), allocatable, intent(out) :: error

   call check(error, landscape_conv_level_known("superloose"))
   if (allocated(error)) return
   call check(error, landscape_conv_level_known("loose"))
   if (allocated(error)) return
   call check(error, landscape_conv_level_known("normal"))
   if (allocated(error)) return
   call check(error, landscape_conv_level_known("tight"))
   if (allocated(error)) return
   call check(error, landscape_conv_level_known("vtight"))
   if (allocated(error)) return
   call check(error, .not. landscape_conv_level_known("medium"))
   if (allocated(error)) return
   call check(error, .not. landscape_conv_level_known(""))

end subroutine test_level_lookup

!> Every enumeration name is matched case insensitively (CLI redesign
!> section 7.5).  The suffix `-lookup` cases above cover the lower case
!> spelling; this case pins the upper case one for every entry point the
!> command line reaches, including the two that *store* a name rather than an
!> identifier.
subroutine test_enum_case(error)
   type(error_type), allocatable, intent(out) :: error
   type(landscape_opt_settings) :: upper
   type(landscape_opt_settings) :: lower

   call check(error, landscape_opt_algorithm_from_name("BFGS"), &
      & p_landscape_opt_bfgs)
   if (allocated(error)) return
   call check(error, landscape_irc_method_from_name("LQA"), p_landscape_irc_lqa)
   if (allocated(error)) return
   call check(error, landscape_ts_method_from_name("PRFO"), p_landscape_ts_prfo)
   if (allocated(error)) return
   call check(error, landscape_scan_coord_from_name("DIHEDRAL"), &
      & p_landscape_scan_coord_dihedral)
   if (allocated(error)) return
   call check(error, landscape_conv_level_known("Vtight"))
   if (allocated(error)) return
   call check(error, landscape_scan_mode_known("RIGID"))
   if (allocated(error)) return
   call check(error, landscape_ts_init_mode_known("FORCE"))
   if (allocated(error)) return
   call check(error, landscape_ts_ninit_known("Random"))
   if (allocated(error)) return
   ! the level setter has to normalise the name as well, otherwise a name that
   ! passed `landscape_conv_level_known` would silently leave the thresholds
   ! at their default
   call upper%set_conv_level("TIGHT")
   call lower%set_conv_level("tight")
   call check(error, upper%g_max, lower%g_max)
   if (allocated(error)) return
   call check(error, upper%g_rms, lower%g_rms)
   if (allocated(error)) return
   call check(error, upper%d_max, lower%d_max)
   if (allocated(error)) return
   call check(error, upper%d_rms, lower%d_rms)
   if (allocated(error)) return
   call check(error, upper%e_tol, lower%e_tol)
   if (allocated(error)) return
   call check(error, upper%engine%acc, lower%engine%acc)
   if (allocated(error)) return
   call check(error, upper%max_steps_multiplier, lower%max_steps_multiplier)

end subroutine test_enum_case

!> The accepted convergence levels exist in one list.
!>
!> The set used to be spelled three times (the `select case` of
!> `landscape_settings_set_conv_level`, `landscape_conv_level_known` and the
!> `v_level` table of `src/xhelp.f90`).  The first two read
!> `landscape_conv_level_names` now and the threshold tables of the setter are
!> declared with `landscape_conv_level_count`, so a level added to the list
!> without a threshold row does not compile; `src/xhelp.f90` is a pre-existing
!> file that may not `use` the suite and stays on its own `v_level` table,
!> which `xtb/landscape-help-values` ties to the binary's acceptance.  This case
!> pins the list itself: the order documented in
!> `landscape_conv_level_index`, the count and the lookup.
subroutine test_level_list(error)
   type(error_type), allocatable, intent(out) :: error
   type(landscape_opt_settings) :: level
   integer :: i

   call check(error, landscape_conv_level_count, 5)
   if (allocated(error)) return
   call check(error, trim(landscape_conv_level_names(1)), "superloose")
   if (allocated(error)) return
   call check(error, trim(landscape_conv_level_names(3)), "normal")
   if (allocated(error)) return
   call check(error, trim(landscape_conv_level_names(5)), "vtight")
   if (allocated(error)) return

   ! every listed name is known, resolves to its own index and sets thresholds
   do i = 1, landscape_conv_level_count
      call check(error, landscape_conv_level_known(landscape_conv_level_names(i)))
      if (allocated(error)) return
      call check(error, landscape_conv_level_index(landscape_conv_level_names(i)), i)
      if (allocated(error)) return
      ! the setter really overwrites every axis: a sentinel value that no
      ! level carries stays only if the name missed its threshold row
      level%g_max = -1.0_wp
      level%g_rms = -1.0_wp
      level%d_max = -1.0_wp
      level%d_rms = -1.0_wp
      level%e_tol = -1.0_wp
      level%engine%acc = -1.0_wp
      level%max_steps_multiplier = -1
      call level%set_conv_level(landscape_conv_level_names(i))
      call check(error, min(level%g_max, level%g_rms, level%d_max, &
         & level%d_rms, level%e_tol, level%engine%acc, &
         & real(level%max_steps_multiplier, wp)) > 0.0_wp)
      if (allocated(error)) return
   end do

   ! and nothing else
   call check(error, landscape_conv_level_index("medium"), 0)
   if (allocated(error)) return
   call check(error, .not. landscape_conv_level_known("extreme"))

end subroutine test_level_list

!> The two shared number formatters, byte for byte.
!>
!> `landscape_real_to_string` takes the edit descriptor as an argument because
!> the five copies it replaced used three different ones: `es15.6` in the three
!> scan modules, `es16.8` in the Hessian reader and `f0.4` in the mode
!> follower.  The reference strings are written here with the bare format, so a
!> changed default or a changed `trim` shows up.
subroutine test_number_strings(error)
   type(error_type), allocatable, intent(out) :: error
   real(wp), parameter :: x = 1.23456789e-3_wp
   character(len=32) :: buffer
   real(wp) :: nan, inf

   call check(error, landscape_int_to_string(0), "0")
   if (allocated(error)) return
   call check(error, landscape_int_to_string(42), "42")
   if (allocated(error)) return
   call check(error, landscape_int_to_string(-7), "-7")
   if (allocated(error)) return

   write(buffer, '(es15.6)') x
   call check(error, landscape_real_to_string(x, '(es15.6)'), trim(buffer))
   if (allocated(error)) return
   write(buffer, '(es16.8)') x
   call check(error, landscape_real_to_string(x, '(es16.8)'), trim(buffer))
   if (allocated(error)) return
   write(buffer, '(f0.4)') x
   call check(error, landscape_real_to_string(x, '(f0.4)'), trim(buffer))
   if (allocated(error)) return

   ! `abs(x) <= huge` would raise `invalid` on a NaN in a checked build;
   ! `ieee_is_finite` does not.
   nan = ieee_value(0.0_wp, ieee_quiet_nan)
   inf = ieee_value(0.0_wp, ieee_positive_inf)
   call check(error, landscape_real_is_finite(x))
   if (allocated(error)) return
   call check(error, .not. landscape_real_is_finite(nan))
   if (allocated(error)) return
   call check(error, .not. landscape_real_is_finite(inf))
   if (allocated(error)) return
   call check(error, all(landscape_real_is_finite([1.0_wp, -2.0_wp])))
   if (allocated(error)) return
   call check(error, .not. any(landscape_real_is_finite([inf, nan])))

end subroutine test_number_strings

!> Human readable names of the engine restart policies
subroutine test_restart_names(error)
   type(error_type), allocatable, intent(out) :: error

   call check(error, landscape_restart_name(p_landscape_restart_off), "off")
   if (allocated(error)) return
   call check(error, landscape_restart_name(p_landscape_restart_pin), "pin")
   if (allocated(error)) return
   call check(error, landscape_restart_name(p_landscape_restart_warm), "warm")
   if (allocated(error)) return
   call check(error, landscape_restart_name(99), "unknown")

end subroutine test_restart_names

end module test_landscape_settings
