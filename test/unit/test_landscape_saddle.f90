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

!> Tests of the transition state search of the landscape suite (M4, stage 1).
!>
!> The coverage is analytic wherever possible, so that the algorithms are
!> checked against exact references and not against a second implementation of
!> the same idea:
!>
!>  * P-RFO shifts: the analytic `lambda_p = (b_k + sqrt(b_k^2 + 4 F_k^2))/2`
!>    of Eq. 4.18a, the sign and the secular residual
!>    `sum F_i^2/(lambda_n - b_i) - lambda_n` of Eq. 4.18b.
!>  * Mueller-Brown surface: P-RFO, B-P-RFO and the dimer converge from
!>    different starts to the same saddle point `(-0.8220016, 0.6243128)`
!>    (PLAN section 9.5); the dimer never calls `hessian`.
!>  * mode following: `argmax|V^T v|` picks the tracked mode through an
!>    eigenvalue crossing, the sign is realigned, a lost mode is reselected and
!>    the phase convention makes the largest component positive.
!>  * rigid body projection: the six basis vectors are orthonormal, a projected
!>    displacement is orthogonal to all of them, a raw P-RFO step is not, and a
!>    linear geometry yields five independent directions.
!>  * safeguards: the curvature lift with more than one negative eigenvalue, the
!>    trust radius clip and the invalid input error paths.
!>  * dimer: the curvature of Eq. 2 against an analytic Hessian, both branches
!>    of the translation force of Eq. 21, Eq. 13 on an exact Fourier model (with
!>    the difference to the earlier Eq. 6 recorded), the torque reduction of the
!>    rotation, the VASP rotation diagnostic, the quick-min update and the
!>    reproducibility of the random initial axis.
module test_landscape_saddle
   use testdrive, only : new_unittest, unittest_type, error_type, check
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use xtb_type_environment
   use xtb_type_molecule
   use xtb_type_restart, only : TRestart
   use xtb_landscape_settings, only : landscape_ts_settings, &
      & p_landscape_ts_prfo, p_landscape_ts_bprfo, p_landscape_ts_dimer, &
      & landscape_ts_method_name, landscape_ts_method_from_name, &
      & landscape_ts_method_known, landscape_ts_init_mode_known
   use xtb_landscape_results, only : landscape_result
   use xtb_landscape_engine, only : landscape_energy_model, landscape_engine, &
      & new_landscape_engine
   use xtb_landscape_saddle_eigenfollow, only : landscape_transition_mode, &
      & landscape_spectrum, landscape_select_transition_mode, &
      & landscape_follow_transition_mode, landscape_rigid_modes, &
      & landscape_project_rigid
   use xtb_landscape_tasks_ts, only : landscape_ts_curvature_band_hit, &
      & landscape_ts_curvature_band_cm1, landscape_ts_can_tighten, &
      & confirm_saddle
   use xtb_landscape_saddle_prfo, only : landscape_prfo_step, &
      & landscape_saddle_optimize
   use xtb_landscape_saddle_bprfo, only : landscape_bprfo_step
   use xtb_landscape_saddle_dimer, only : landscape_dimer_optimize, &
      & landscape_dimer_curvature, landscape_dimer_translation_force, &
      & landscape_dimer_rotation_direction, landscape_dimer_rotation_angle, &
      & landscape_dimer_velocity_update, landscape_dimer_initial_axis, &
      & landscape_dimer_mobility, &
      & landscape_dimer_rotate
   use xtb_landscape_testpot, only : landscape_muller_brown_model, &
      & new_muller_brown_model
   implicit none
   private

   !> Exact Mueller-Brown saddle point (dimensionless surface coordinates)
   real(wp), parameter :: mb_saddle(2) = [-0.82200156_wp, 0.62431280_wp]

   !> Mueller-Brown surface whose Hessian is not available.
   !>
   !> The dimer has to find the saddle point without it; a `hessian` call that
   !> returns `fail = .true.` aborts the search, so a converged run proves that
   !> the dimer never asked for one (PLAN section 9.5).
   type, extends(landscape_muller_brown_model) :: no_hessian_model
      !> Number of Hessian requests (must stay zero)
      integer :: n_hess = 0
   contains
      !> Hessian, deliberately unavailable
      procedure :: hessian => no_hessian_hessian
   end type no_hessian_model

   !> A nine dimensional quadratic with one internal and one rigid force.
   !>
   !> `E = 1/2 (x-x0)^T H (x-x0) + c^T (x-x0)`; `H` has one negative internal
   !> eigenvalue (a saddle), two positive internal ones and a small positive
   !> curvature on the rigid body subspace, while `c` carries a component along
   !> a rigid body direction.  A search that projects its step onto the internal
   !> space cannot follow that component and stays on the internal manifold; a
   !> search that does not drifts along it.  This is the falsifiable coverage of
   !> the driver level rigid body projection (M4 of
   !> `29-m4-independent-review.md`).
   type, extends(landscape_energy_model) :: rigid_drift_model
      !> Hessian (9, 9)
      real(wp), allocatable :: h(:, :)
      !> Linear gradient term of the rigid direction
      real(wp), allocatable :: c(:)
      !> Expansion point (bohr)
      real(wp), allocatable :: x0(:)
   contains
      !> Energy and gradient
      procedure :: evaluate => drift_evaluate
      !> Exact (constant) Hessian
      procedure :: hessian => drift_hessian
   end type rigid_drift_model

   !> Synthetic Hessian with an internal mode inside the curvature band.
   !>
   !> `confirm_saddle` tightens its confirmation Hessian only when the spectrum
   !> holds a mode inside `0 +/- 3 cm^-1`, and the `hessian_source` spelling
   !> decides whether that re-decision may happen.  A real GFN2 surface whose
   !> confirmation Hessian has an internal mode inside the band is not
   !> constructible (`docs/_landscape-recon/30-m4-fix-report.md`), so the wiring
   !> is driven with the isotropic constant Hessian `H = curvature * I`: after
   !> the mass weighting the projected internal modes of water sit at
   !> `curvature = 1e-10` hartree/bohr^2 at `0.0484`, `0.0491` and
   !> `0.0512 cm^-1` (measured), all inside the band, while the six projected
   !> rigid modes are exactly zero.  The matrix does not depend on the finite
   !> difference step, so the default and the tightened resolutions return the
   !> same curvature.
   type, extends(landscape_energy_model) :: band_hessian_model
      !> Isotropic curvature (hartree/bohr^2)
      real(wp) :: curvature = 1.0e-10_wp
   contains
      !> Energy and gradient (the confirmation never asks for them)
      procedure :: evaluate => band_evaluate
      !> Constant isotropic Hessian
      procedure :: hessian => band_hessian
   end type band_hessian_model

   !> Water geometry in bohr (a nonlinear three atom molecule)
   integer, parameter :: nat_water = 3
   !> Element symbols of the water fixture
   character(len=*), parameter :: sym_water(nat_water) = ["O", "H", "H"]
   real(wp), parameter :: xyz_water(3, nat_water) = reshape([ &
      &  0.0_wp,  0.0_wp,  0.0_wp, &
      &  0.0_wp,  1.8_wp,  0.0_wp, &
      &  0.0_wp, -0.9_wp,  1.5_wp], [3, nat_water])

   !> Linear CO2 in bohr (only two independent rotations)
   integer, parameter :: nat_co2 = 3
   real(wp), parameter :: xyz_co2(3, nat_co2) = reshape([ &
      &  0.0_wp,  0.0_wp, -2.2_wp, &
      &  0.0_wp,  0.0_wp,  0.0_wp, &
      &  0.0_wp,  0.0_wp,  2.2_wp], [3, nat_co2])

   public :: collect_landscape_saddle

contains

!> Collect all exported unit tests
subroutine collect_landscape_saddle(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("method-names", test_method_names), &
      new_unittest("spectrum-ascending", test_spectrum_ascending), &
      new_unittest("mode-select-negative", test_mode_select_negative), &
      new_unittest("mode-select-no-negative", test_mode_select_no_negative), &
      new_unittest("mode-select-lowest-internal", test_mode_select_lowest_internal), &
      new_unittest("mode-follow-argmax", test_mode_follow_argmax), &
      new_unittest("mode-follow-sign", test_mode_follow_sign), &
      new_unittest("mode-follow-low-overlap", test_mode_follow_low_overlap), &
      new_unittest("mode-phase", test_mode_phase), &
      new_unittest("rigid-basis-orthonormal", test_rigid_basis), &
      new_unittest("rigid-linear-basis", test_rigid_linear), &
      new_unittest("rigid-projection", test_rigid_projection), &
      new_unittest("prfo-shifts", test_prfo_shifts), &
      new_unittest("prfo-lambda-n", test_prfo_lambda_n), &
      new_unittest("prfo-lambda-n-no-root", test_prfo_lambda_n_no_root), &
      new_unittest("prfo-trust-clip", test_prfo_trust_clip), &
      new_unittest("prfo-hlow-lift", test_prfo_hlow_lift), &
      new_unittest("prfo-noise-negative", test_prfo_noise_negative), &
      new_unittest("prfo-negative-downhill", test_prfo_negative_downhill), &
      new_unittest("prfo-init-mode-case", test_prfo_init_mode_case), &
      new_unittest("prfo-invalid", test_prfo_invalid), &
      new_unittest("prfo-rigid-projection", test_prfo_rigid), &
      new_unittest("bprfo-radius-split", test_bprfo_radius_split), &
      new_unittest("mb-prfo-saddle", test_mb_prfo), &
      new_unittest("mb-bprfo-saddle", test_mb_bprfo), &
      new_unittest("mb-methods-agree", test_mb_agree), &
      new_unittest("ts-curvature-band", test_ts_curvature_band), &
      new_unittest("ts-can-tighten-case", test_ts_can_tighten_case), &
      new_unittest("ts-confirm-tighten-case", test_ts_confirm_tighten_case), &
      new_unittest("driver-rigid-projection", test_driver_rigid_projection), &
      new_unittest("dimer-curvature", test_dimer_curvature), &
      new_unittest("dimer-translation-force", test_dimer_translation_force), &
      new_unittest("dimer-rotation-angle", test_dimer_rotation_angle), &
      new_unittest("dimer-rotation-alignment", test_dimer_rotation_alignment), &
      new_unittest("dimer-rotation-warning", test_dimer_rotation_warning), &
      new_unittest("dimer-velocity-update", test_dimer_velocity), &
      new_unittest("dimer-mobility", test_dimer_mobility), &
      new_unittest("dimer-initial-axis", test_dimer_initial_axis), &
      new_unittest("dimer-axis-case", test_dimer_axis_case), &
      new_unittest("dimer-invalid", test_dimer_invalid), &
      new_unittest("dimer-mb-saddle", test_dimer_mb), &
      new_unittest("dimer-no-hessian", test_dimer_no_hessian), &
      new_unittest("dimer-reproducible", test_dimer_reproducible) &
      ]

end subroutine collect_landscape_saddle

!> Method identifiers, names and lookups of the transition state settings.
subroutine test_method_names(error)
   type(error_type), allocatable, intent(out) :: error

   type(landscape_ts_settings) :: set

   call check(error, trim(landscape_ts_method_name(p_landscape_ts_prfo)) == 'prfo')
   if (allocated(error)) return
   call check(error, trim(landscape_ts_method_name(p_landscape_ts_bprfo)) == 'bprfo')
   if (allocated(error)) return
   call check(error, trim(landscape_ts_method_name(p_landscape_ts_dimer)) == 'dimer')
   if (allocated(error)) return
   call check(error, landscape_ts_method_from_name('bprfo') == p_landscape_ts_bprfo)
   if (allocated(error)) return
   call check(error, .not. landscape_ts_method_known('whatever'))
   if (allocated(error)) return
   call check(error, landscape_ts_init_mode_known('force'))
   if (allocated(error)) return
   call check(error, .not. landscape_ts_init_mode_known('h1'))
   if (allocated(error)) return
   ! the default method of the settings is P-RFO
   call check(error, set%method == p_landscape_ts_prfo)

end subroutine test_method_names

!> The spectrum helper returns ascending eigenvalues and orthonormal vectors.
subroutine test_spectrum_ascending(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   real(wp), allocatable :: hess(:, :), evals(:), evecs(:, :)
   logical :: fail

   call init(env)
   allocate(hess(3, 3), source=0.0_wp)
   hess(1, 1) = 2.0_wp
   hess(2, 2) = -1.0_wp
   hess(3, 3) = 0.5_wp
   hess(1, 2) = 0.2_wp
   hess(2, 1) = 0.2_wp

   call landscape_spectrum(env, hess, evals, evecs, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, evals(1) <= evals(2) .and. evals(2) <= evals(3))
   if (allocated(error)) return
   call check(error, abs(dot_product(evecs(:, 1), evecs(:, 2))) < 1.0e-12_wp)
   if (allocated(error)) return
   call check(error, abs(sum(evecs(:, 1)**2) - 1.0_wp) < 1.0e-12_wp)

end subroutine test_spectrum_ascending

!> A Hessian with negative eigenvalues selects the most negative mode.
subroutine test_mode_select_negative(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   real(wp) :: evals(3), evecs(3, 3)
   type(landscape_transition_mode) :: mode
   logical :: fail

   call init(env)
   evals = [-1.0_wp, 0.5_wp, 2.0_wp]
   call identity(evecs)
   call landscape_select_transition_mode(env, evals, evecs, 1, mode, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, mode%index == 1)
   if (allocated(error)) return
   call check(error, abs(abs(mode%vector(1)) - 1.0_wp) < 1.0e-12_wp)
   if (allocated(error)) return
   call check(error, .not. mode%followed)

end subroutine test_mode_select_negative

!> A Hessian without a negative eigenvalue still starts a search: the lowest
!> mode is followed and a warning is emitted.
subroutine test_mode_select_no_negative(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   real(wp) :: evals(3), evecs(3, 3)
   type(landscape_transition_mode) :: mode
   logical :: fail
   integer :: nlog

   call init(env)
   evals = [0.4_wp, 0.9_wp, 2.0_wp]
   call identity(evecs)
   nlog = env%nLog
   call landscape_select_transition_mode(env, evals, evecs, 1, mode, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, mode%index == 1)
   if (allocated(error)) return
   call check(error, env%nLog > nlog)

end subroutine test_mode_select_no_negative

!> A spectrum whose lowest entries are rigid body zero modes selects the lowest
!> **internal** mode (MEDIUM-1 of `31-m4-review-round2.md`).
!>
!> The ascending spectrum of a molecule starts with the six (five for a linear
!> molecule) rigid body directions, which the projected numerical Hessian leaves
!> at `O(eps ||H||)` with a random sign.  Selecting index 1 on a minimum
!> therefore followed a rigid body direction whose ascent is removed again by
!> the step projection, and the search degenerated into a minimisation.  The
!> synthetic spectrum below is exactly that situation: two entries inside the
!> curvature noise band, then the two internal modes.  A fixed `index = 1`
!> assignment returns the `-1e-16` mode and fails this test.
subroutine test_mode_select_lowest_internal(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   real(wp) :: evals(4), evecs(4, 4)
   type(landscape_transition_mode) :: mode
   logical :: fail
   integer :: nlog

   call init(env)
   evals = [-1.0e-16_wp, -1.0e-17_wp, 0.3_wp, 0.5_wp]
   call identity(evecs)
   nlog = env%nLog
   call landscape_select_transition_mode(env, evals, evecs, 1, mode, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, mode%index == 3,&
      & "a rigid body zero mode was followed instead of the lowest internal mode")
   if (allocated(error)) return
   call check(error, abs(mode%eigenvalue - 0.3_wp) < 1.0e-14_wp)
   if (allocated(error)) return
   call check(error, abs(abs(mode%vector(3)) - 1.0_wp) < 1.0e-12_wp)
   if (allocated(error)) return
   ! the mismatch between "no negative eigenvalue" and the request is reported
   call check(error, env%nLog > nlog)
   if (allocated(error)) return

   ! with one negative eigenvalue but a higher target the most negative mode is
   ! still the one to follow (it is the first entry of the ascending spectrum)
   evals = [-0.5_wp, -1.0e-16_wp, 0.3_wp, 0.5_wp]
   call identity(evecs)
   call landscape_select_transition_mode(env, evals, evecs, 2, mode, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, mode%index == 1,&
      & "an unreachable target mode did not fall back to the most negative mode")
   if (allocated(error)) return
   call check(error, abs(mode%eigenvalue + 0.5_wp) < 1.0e-14_wp)

end subroutine test_mode_select_lowest_internal

!> `argmax|V^T v|` follows the tracked mode through an eigenvalue crossing.
!>
!> The new spectrum has the eigenvector aligned with the tracked mode at
!> position 3; a fixed index assignment (position 1 here) would lose it.
subroutine test_mode_follow_argmax(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   real(wp) :: evals(3), evecs(3, 3)
   type(landscape_transition_mode) :: mode
   logical :: fail

   call init(env)
   evals = [-1.0_wp, 0.5_wp, 2.0_wp]
   call identity(evecs)
   call landscape_select_transition_mode(env, evals, evecs, 1, mode, fail)
   if (allocated(error)) return

   ! crossed spectrum: the mode that was index 1 is now index 3
   evals = [-1.2_wp, 0.4_wp, 1.9_wp]
   call identity(evecs)
   call swap_columns(evecs, 1, 3)
   call landscape_follow_transition_mode(env, evals, evecs, mode, 0.3_wp, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, mode%followed)
   if (allocated(error)) return
   call check(error, mode%index == 3,&
      & "the followed mode is not the argmax overlap candidate")
   if (allocated(error)) return
   call check(error, mode%overlap > 0.99_wp)

end subroutine test_mode_follow_argmax

!> The new eigenvector is realigned to keep the sign of the tracked vector.
!>
!> The tracked vector is deliberately set to `-e1` before the update: after the
!> following step it has to stay in that direction (`dot(v_new, v_old) > 0`),
!> which is only true with the `sign(overlap)` realignment.  The phase
!> convention is *not* applied here - it belongs to the initial selection and
!> would silently repair a missing realignment whenever the largest component
!> happens to be positive.
subroutine test_mode_follow_sign(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   real(wp) :: evals(3), evecs(3, 3), v_old(3)
   type(landscape_transition_mode) :: mode
   logical :: fail

   call init(env)
   evals = [-1.0_wp, 0.5_wp, 2.0_wp]
   call identity(evecs)
   call landscape_select_transition_mode(env, evals, evecs, 1, mode, fail)
   if (allocated(error)) return
   mode%vector = -mode%vector
   v_old = mode%vector

   ! the same eigenvector with the opposite sign
   evals = [-1.1_wp, 0.5_wp, 2.0_wp]
   call identity(evecs)
   call landscape_follow_transition_mode(env, evals, evecs, mode, 0.3_wp, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, mode%followed)
   if (allocated(error)) return
   call check(error, dot_product(mode%vector, v_old) > 0.99_wp,&
      & "the sign realignment lost the direction of the tracked vector")

end subroutine test_mode_follow_sign

!> An overlap below the guard reselects the mode instead of following it.
subroutine test_mode_follow_low_overlap(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   real(wp) :: evals(3), evecs(3, 3)
   type(landscape_transition_mode) :: mode
   logical :: fail
   integer :: nlog

   call init(env)
   evals = [-1.0_wp, 0.5_wp, 2.0_wp]
   call identity(evecs)
   call landscape_select_transition_mode(env, evals, evecs, 1, mode, fail)
   if (allocated(error)) return

   ! a new frame whose eigenvectors are all orthogonal to the tracked e1
   evals = [-0.8_wp, 0.6_wp, 2.1_wp]
   evecs = 0.0_wp
   evecs(2, 1) = 1.0_wp
   evecs(3, 2) = 1.0_wp
   evecs(2, 3) = sqrt(0.5_wp)
   evecs(3, 3) = sqrt(0.5_wp)
   nlog = env%nLog
   call landscape_follow_transition_mode(env, evals, evecs, mode, 0.3_wp, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, .not. mode%followed)
   if (allocated(error)) return
   call check(error, mode%overlap < 0.3_wp)
   if (allocated(error)) return
   call check(error, env%nLog > nlog)

end subroutine test_mode_follow_low_overlap

!> The phase convention makes the largest component of the vector positive.
subroutine test_mode_phase(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   real(wp) :: evals(3), evecs(3, 3)
   type(landscape_transition_mode) :: mode
   real(wp) :: vec(3)
   logical :: fail

   call init(env)
   evals = [-1.0_wp, 0.5_wp, 2.0_wp]
   vec = [0.1_wp, -0.9_wp, 0.2_wp]
   evecs = 0.0_wp
   evecs(:, 1) = vec/sqrt(sum(vec**2))
   evecs(:, 2) = [1.0_wp, 0.0_wp, 0.0_wp]
   evecs(:, 3) = [0.0_wp, 0.0_wp, 1.0_wp]
   call landscape_select_transition_mode(env, evals, evecs, 1, mode, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, mode%phase_index == 2)
   if (allocated(error)) return
   call check(error, mode%phase_value > 0.0_wp)
   if (allocated(error)) return
   call check(error, mode%vector(2) > 0.0_wp)

end subroutine test_mode_phase

!> The rigid body basis of a nonlinear molecule has six orthonormal columns.
subroutine test_rigid_basis(error)
   type(error_type), allocatable, intent(out) :: error
   real(wp), allocatable :: basis(:, :), x(:), v(:), nonrigid(:)
   integer :: n_rigid, i, j
   real(wp) :: rigid_change, nonrigid_change

   x = reshape(xyz_water, [3*nat_water])
   call landscape_rigid_modes(x, basis, n_rigid)
   call check(error, n_rigid == 6)
   if (allocated(error)) return
   do i = 1, n_rigid
      do j = i, n_rigid
         if (i == j) then
            call check(error, abs(dot_product(basis(:, i), basis(:, j)) - 1.0_wp) < 1.0e-12_wp)
         else
            call check(error, abs(dot_product(basis(:, i), basis(:, j))) < 1.0e-12_wp)
         end if
         if (allocated(error)) return
      end do
   end do

   ! A rigid displacement changes the internal coordinates only to second order;
   ! an internal displacement of the same size changes them to first order.  The
   ! comparison is what makes the basis meaningful: any set of six orthonormal
   ! vectors would pass the Gram matrix check above.
   v = x
   do i = 1, n_rigid
      v = v + 1.0e-3_wp*basis(:, i)
   end do
   nonrigid = x
   nonrigid(2) = nonrigid(2) + 1.0e-3_wp
   rigid_change = internal_distance_change(x, v)
   nonrigid_change = internal_distance_change(x, nonrigid)
   call check(error, nonrigid_change > 1.0e-4_wp)
   if (allocated(error)) return
   call check(error, rigid_change < 0.1_wp*nonrigid_change)

end subroutine test_rigid_basis

!> A linear molecule has only five independent rigid body directions.
subroutine test_rigid_linear(error)
   type(error_type), allocatable, intent(out) :: error
   real(wp), allocatable :: basis(:, :), x(:)
   integer :: n_rigid

   x = reshape(xyz_co2, [3*nat_co2])
   call landscape_rigid_modes(x, basis, n_rigid)
   call check(error, n_rigid == 5,&
      & "a linear molecule must have five rigid body directions")

end subroutine test_rigid_linear

!> Projection removes the component along every rigid body direction, while a
!> raw displacement keeps it (the M1 debt that M4 has to close).
subroutine test_rigid_projection(error)
   type(error_type), allocatable, intent(out) :: error
   real(wp), allocatable :: basis(:, :), x(:), v(:), raw(:)
   integer :: n_rigid, i
   real(wp) :: max_before, max_after

   x = reshape(xyz_water, [3*nat_water])
   call landscape_rigid_modes(x, basis, n_rigid)
   if (allocated(error)) return

   allocate(raw(3*nat_water))
   raw = 0.0_wp
   ! a displacement with an obvious rigid body component: a pure translation
   raw(2) = 0.2_wp
   raw(5) = 0.2_wp
   raw(8) = 0.2_wp
   v = raw
   call landscape_project_rigid(v, basis, n_rigid)

   max_before = 0.0_wp
   max_after = 0.0_wp
   do i = 1, n_rigid
      max_before = max(max_before, abs(dot_product(basis(:, i), raw)))
      max_after = max(max_after, abs(dot_product(basis(:, i), v)))
   end do
   call check(error, max_before > 1.0e-3_wp,&
      & "the test vector should have a rigid body component")
   if (allocated(error)) return
   call check(error, max_after < 1.0e-12_wp,&
      & "projection left a rigid body component")

end subroutine test_rigid_projection

!> The analytic `lambda_p` of Eq. 4.18a and the sign condition of `lambda_n`.
subroutine test_prfo_shifts(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   real(wp) :: hess(2, 2), gradient(2), trans(2), step(2)
   real(wp) :: lam_p, lam_n, predicted, ref_p
   logical :: fail, lifted

   call init(env)
   hess = 0.0_wp
   hess(1, 1) = -0.5_wp
   hess(2, 2) = 2.0_wp
   gradient = [0.1_wp, 0.3_wp]
   trans = [1.0_wp, 0.0_wp]

   call landscape_prfo_step(env, hess, gradient, trans, 10.0_wp, 1.0e-2_wp, &
      & step, predicted, fail, lambda_p=lam_p, lambda_n=lam_n, lifted=lifted)
   call check(error, .not. fail)
   if (allocated(error)) return
   ref_p = 0.5_wp*(-0.5_wp) + 0.5_wp*sqrt(0.25_wp + 4.0_wp*0.01_wp)
   call check(error, abs(lam_p - ref_p) < 1.0e-12_wp, "lambda_p does not match Eq. 4.18a")
   if (allocated(error)) return
   call check(error, lam_n < 0.0_wp, "lambda_n must be negative")
   if (allocated(error)) return
   call check(error, .not. lifted)

end subroutine test_prfo_shifts

!> `lambda_n` solves the secular equation `sum F_i^2/(lambda_n - b_i) = lambda_n`.
!>
!> A solution that skips the iteration and returns zero leaves a residual of
!> `-sum F_i^2/b_i`, which is far from zero for this case.
subroutine test_prfo_lambda_n(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   real(wp) :: hess(3, 3), gradient(3), trans(3), step(3)
   real(wp), allocatable :: evals(:), evecs(:, :)
   real(wp) :: lam_n, predicted, residual
   logical :: fail

   call init(env)
   hess = 0.0_wp
   hess(1, 1) = 0.1_wp
   hess(2, 2) = 1.0_wp
   hess(3, 3) = 4.0_wp
   gradient = [0.5_wp, 1.0_wp, 1.0_wp]
   trans = [1.0_wp, 0.0_wp, 0.0_wp]

   call landscape_prfo_step(env, hess, gradient, trans, 10.0_wp, 1.0e-2_wp, &
      & step, predicted, fail, lambda_n=lam_n)
   call check(error, .not. fail)
   if (allocated(error)) return
   call landscape_spectrum(env, hess, evals, evecs, fail)
   if (allocated(error)) return
   residual = -lam_n
   residual = residual + 1.0_wp**2/(lam_n - 1.0_wp)
   residual = residual + 1.0_wp**2/(lam_n - 4.0_wp)
   call check(error, abs(residual) < 1.0e-6_wp,&
      & "lambda_n does not solve the RFO secular equation")
   if (allocated(error)) return
   call check(error, lam_n < 0.0_wp)
   if (allocated(error)) return
   ! the step follows -F_i/(b_i - lambda_n) with that one shift, so
   ! dq3/dq2 = (b2 - lambda_n)/(b3 - lambda_n); a Newton step (lambda_n = 0)
   ! would give b2/b3 = 0.25 instead of about 0.37
   call check(error, abs(step(3)/step(2) - (1.0_wp - lam_n)/(4.0_wp - lam_n)) < 1.0e-6_wp)

end subroutine test_prfo_lambda_n

!> A `lambda_n` equation without a root in its bracket is reported and does not
!> return a non-root (MEDIUM-2 of `31-m4-review-round2.md`).
!>
!> `residual` is strictly decreasing on `(-inf, min_{i/=k} b_i)` and the bracket
!> upper end is `bmin - mu_margin`, so `residual(bmin - mu_margin) > 0` means
!> that the root lies above the bracket - the equation has no downhill branch at
!> all - and bisection would return the bracket end.  It takes a downhill
!> gradient of `|f_i| <~ 1e-8` for that (the inequality is documented in
!> `landscape_rfo_downhill_shift`), together with a downhill eigenvalue that is
!> negative but inside the curvature noise band, so this is the residual
!> verification of the equation and not of its input.
subroutine test_prfo_lambda_n_no_root(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   real(wp) :: hess(3, 3), gradient(3), trans(3), step(3)
   real(wp) :: lam_n, predicted, residual_hi, bmin, margin
   logical :: fail, lifted
   integer :: nlog

   call init(env)
   hess = 0.0_wp
   hess(1, 1) = -0.5_wp          ! the followed mode
   hess(2, 2) = -1.0e-16_wp      ! projected rigid body zero mode (noise band)
   hess(3, 3) = 0.5_wp
   ! the followed mode carries no gradient (`f_k = 0`): a finite `f_k` would
   ! swamp the `1e-22` downhill contributions of the sum `sum F_i^2` in double
   ! precision and turn the case into the documented all-zero-gradient branch
   gradient = [0.0_wp, 1.0e-11_wp, 1.0e-11_wp]
   trans = [1.0_wp, 0.0_wp, 0.0_wp]

   ! the residual at the upper bracket end, computed exactly as the solver does
   bmin = -1.0e-16_wp
   margin = 1.0e-8_wp
   residual_hi = gradient(2)**2/(bmin - margin - bmin) + &
      & gradient(3)**2/(bmin - margin - hess(3, 3)) - (bmin - margin)
   call check(error, residual_hi > 0.0_wp,&
      & "the bracket of this test does not miss the root")
   if (allocated(error)) return

   nlog = env%nLog
   call landscape_prfo_step(env, hess, gradient, trans, 1.0e-1_wp, 1.0e-2_wp, &
      & step, predicted, fail, lambda_n=lam_n, lifted=lifted)
   call check(error, .not. fail)
   if (allocated(error)) return
   ! the equation has no negative root: the shift degrades to the degenerate
   ! lambda_n = 0 (a plain Newton descent), it must not return the bracket end
   call check(error, lam_n == 0.0_wp, &
      & "a shift that does not solve the secular equation was returned")
   if (allocated(error)) return
   call check(error, env%nLog > nlog,&
      & "the missing root was not reported")
   if (allocated(error)) return
   call check(error, .not. lifted)
   if (allocated(error)) return
   call check(error, maxval(abs(step)) <= 1.0e-1_wp + 1.0e-12_wp)

end subroutine test_prfo_lambda_n_no_root

!> The step is clipped to the trust radius.
subroutine test_prfo_trust_clip(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   real(wp) :: hess(2, 2), gradient(2), trans(2), step(2), predicted
   real(wp) :: trust
   logical :: fail

   call init(env)
   hess = 0.0_wp
   hess(1, 1) = -0.5_wp
   hess(2, 2) = 0.5_wp
   gradient = [1.0_wp, 1.0_wp]
   trans = [1.0_wp, 0.0_wp]
   trust = 0.05_wp

   call landscape_prfo_step(env, hess, gradient, trans, trust, 1.0e-2_wp, &
      & step, predicted, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, sqrt(sum(step**2)) <= trust + 1.0e-12_wp,&
      & "the step exceeded the trust radius")

end subroutine test_prfo_trust_clip

!> More than one negative eigenvalue triggers the `hlow` lift.
subroutine test_prfo_hlow_lift(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   real(wp) :: hess(3, 3), gradient(3), trans(3), step(3), predicted, lam_n
   logical :: fail, lifted
   integer :: nlog

   call init(env)
   hess = 0.0_wp
   hess(1, 1) = -0.5_wp
   hess(2, 2) = -0.2_wp
   hess(3, 3) = 2.0_wp
   gradient = [0.1_wp, 0.2_wp, 0.3_wp]
   trans = [1.0_wp, 0.0_wp, 0.0_wp]

   nlog = env%nLog
   call landscape_prfo_step(env, hess, gradient, trans, 10.0_wp, 1.0e-2_wp, &
      & step, predicted, fail, lambda_n=lam_n, lifted=lifted)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, lifted, "the curvature lift was not reported")
   if (allocated(error)) return
   call check(error, env%nLog > nlog)
   if (allocated(error)) return
   ! after the lift the second downhill eigenvalue is hlow, so lambda_n stays
   ! below it
   call check(error, lam_n < 1.0e-2_wp)

end subroutine test_prfo_hlow_lift

!> A round off negative eigenvalue in the downhill subspace is not curvature.
!>
!> The rigid body directions of a numerical Hessian that were not projected
!> away are lifted to `~eps*||H||`; on a clean GFN2 molecule that is about
!> 1e-11 and made the `hlow` lift engage on *every* step (H4 of
!> `29-m4-independent-review.md`).  The noise threshold of
!> `landscape_curvature_noise_tol` must keep `lifted` false and emit no
!> warning for such an eigenvalue, while `prfo-negative-downhill` still pins
!> that a genuine negative curvature is lifted.
subroutine test_prfo_noise_negative(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   real(wp) :: hess(3, 3), gradient(3), trans(3), step(3), predicted
   logical :: fail, lifted
   integer :: nlog

   call init(env)
   hess = 0.0_wp
   hess(1, 1) = -0.5_wp
   hess(2, 2) = -1.0e-11_wp   ! rigid body round off, not curvature
   hess(3, 3) = 2.0_wp
   gradient = [0.1_wp, 0.2_wp, 0.3_wp]
   trans = [1.0_wp, 0.0_wp, 0.0_wp]

   nlog = env%nLog
   call landscape_prfo_step(env, hess, gradient, trans, 10.0_wp, 1.0e-2_wp, &
      & step, predicted, fail, lifted=lifted)
   call check(error, .not. fail, "the step could not be constructed")
   if (allocated(error)) return
   call check(error, .not. lifted, &
      & "a round off negative eigenvalue triggered the curvature lift")
   if (allocated(error)) return
   call check(error, env%nLog == nlog, &
      & "a spurious curvature lift warning was emitted")

end subroutine test_prfo_noise_negative

!> A single negative eigenvalue that is *not* the followed mode must be lifted.
!>
!> The followed mode (`trans = e1`) has positive curvature here, so the negative
!> mode at index 2 sits in the downhill subspace and would make the secular
!> equation unsolvable below zero (`residual(0) > 0`); the lift brings the
!> downhill subspace back under `hlow`.
subroutine test_prfo_negative_downhill(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   real(wp) :: hess(3, 3), gradient(3), trans(3), step(3), predicted, lam_n
   logical :: fail, lifted

   call init(env)
   hess = 0.0_wp
   hess(1, 1) = 0.2_wp
   hess(2, 2) = -0.5_wp
   hess(3, 3) = 2.0_wp
   gradient = [0.1_wp, 0.2_wp, 0.3_wp]
   trans = [1.0_wp, 0.0_wp, 0.0_wp]

   call landscape_prfo_step(env, hess, gradient, trans, 10.0_wp, 1.0e-2_wp, &
      & step, predicted, fail, lambda_n=lam_n, lifted=lifted)
   call check(error, .not. fail, "the step could not be constructed")
   if (allocated(error)) return
   call check(error, lifted, "the downhill negative curvature was not lifted")
   if (allocated(error)) return
   call check(error, lam_n < 0.0_wp)

end subroutine test_prfo_negative_downhill

!> Invalid inputs are reported instead of producing a step.
subroutine test_prfo_invalid(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   real(wp) :: hess(2, 2), gradient(2), trans(2), step(2), predicted
   logical :: fail

   call init(env)
   hess = 0.0_wp
   hess(1, 1) = -0.5_wp
   hess(2, 2) = 0.5_wp
   gradient = [1.0_wp, 1.0_wp]
   trans = [1.0_wp, 0.0_wp]

   call landscape_prfo_step(env, hess, gradient, trans, 0.0_wp, 1.0e-2_wp, &
      & step, predicted, fail)
   call check(error, fail, "a non-positive trust radius must be rejected")

end subroutine test_prfo_invalid

!> A raw P-RFO step of a rigid-mixed gradient carries a rigid body component;
!> the projection used by the driver removes it.
subroutine test_prfo_rigid(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   real(wp), allocatable :: basis(:, :), x(:), hess(:, :), gradient(:), trans(:)
   real(wp), allocatable :: step_null(:), step_proj(:)
   real(wp) :: predicted, max_rigid, max_proj
   real(wp) :: lam_p, lam_n
   logical :: fail, lifted
   integer :: n_rigid, i, n

   call init(env)
   x = reshape(xyz_water, [3*nat_water])
   n = size(x)
   call landscape_rigid_modes(x, basis, n_rigid)
   if (allocated(error)) return

   allocate(hess(n, n), gradient(n), trans(n), step_null(n), step_proj(n))
   hess = 0.0_wp
   do i = 1, n
      hess(i, i) = 0.8_wp
   end do
   hess(1, 1) = -0.5_wp
   gradient = 0.05_wp
   trans = 0.0_wp
   trans(1) = 1.0_wp

   call landscape_prfo_step(env, hess, gradient, trans, 10.0_wp, 1.0e-2_wp, &
      & step_null, predicted, fail, lambda_p=lam_p, lambda_n=lam_n, lifted=lifted)
   call check(error, .not. fail)
   if (allocated(error)) return
   max_rigid = 0.0_wp
   do i = 1, n_rigid
      max_rigid = max(max_rigid, abs(dot_product(basis(:, i), step_null)))
   end do
   call check(error, max_rigid > 1.0e-3_wp,&
      & "the raw step should carry a rigid body component")

   step_proj = step_null
   call landscape_project_rigid(step_proj, basis, n_rigid)
   max_proj = 0.0_wp
   do i = 1, n_rigid
      max_proj = max(max_proj, abs(dot_product(basis(:, i), step_proj)))
   end do
   call check(error, max_proj < 1.0e-12_wp,&
      & "the projected step is not orthogonal to the rigid body space")

end subroutine test_prfo_rigid

!> B-P-RFO partitions the trust radius between the two subspaces.
subroutine test_bprfo_radius_split(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   real(wp) :: hess(3, 3), gradient(3), trans(3), step(3), predicted
   real(wp) :: trust, r_p, r_n, uphill, downhill
   logical :: fail

   call init(env)
   hess = 0.0_wp
   hess(1, 1) = -0.5_wp
   hess(2, 2) = 0.4_wp
   hess(3, 3) = 1.6_wp
   gradient = [0.7_wp, 0.5_wp, 0.9_wp]
   trans = [1.0_wp, 0.0_wp, 0.0_wp]
   trust = 0.3_wp

   call landscape_bprfo_step(env, hess, gradient, trans, trust, 0.5_wp, &
      & 1.0e-2_wp, step, predicted, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   r_p = 0.5_wp*trust
   r_n = sqrt(0.75_wp)*trust
   ! the two subspaces are the followed mode (here e1) and its complement
   uphill = abs(step(1))
   downhill = sqrt(step(2)**2 + step(3)**2)
   call check(error, abs(uphill - r_p) < 1.0e-9_wp,&
      & "the uphill component is not at its sub-radius")
   if (allocated(error)) return
   call check(error, abs(downhill - r_n) < 1.0e-9_wp,&
      & "the downhill component is not at its sub-radius")
   if (allocated(error)) return
   call check(error, sqrt(sum(step**2)) <= trust + 1.0e-12_wp)

end subroutine test_bprfo_radius_split

!> P-RFO converges to the Mueller-Brown saddle point.
subroutine test_mb_prfo(error)
   type(error_type), allocatable, intent(out) :: error

   call run_mb(error, p_landscape_ts_prfo, [-0.6_wp, 0.4_wp])

end subroutine test_mb_prfo

!> B-P-RFO converges to the same Mueller-Brown saddle point.
subroutine test_mb_bprfo(error)
   type(error_type), allocatable, intent(out) :: error

   call run_mb(error, p_landscape_ts_bprfo, [-0.7_wp, 0.45_wp])

end subroutine test_mb_bprfo

!> P-RFO and B-P-RFO agree to better than 2e-5 bohr on Mueller-Brown.
!>
!> The delivered RMSD is 6.3022e-06 bohr (the M4 acceptance number of
!> `AGENTS.md` section 2), so the gate is 2e-5 bohr: 3x above the delivered
!> value and 50x tighter than the 1e-3 tripwire it replaces.  A change that
!> moves the two searches apart by more than 1e-5 bohr is a change of the
!> saddle machinery and is re-baselined in a report, not absorbed here.
!>
!> The quantity is a difference of two independently converged iterates, so a
!> different convergence path can move it by the last step length.  Measured
!> on both compilers the project builds with - gfortran 12.2.0 and the gating
!> CI's 14.4.0 - the value is **identical** (6.3022e-06 bohr printed, and the
!> two converged points agree to all 14 printed digits), so the datum is not a
!> compiler artefact of this machine.  The tolerance is nevertheless backed by
!> the two searches' own residuals below: a run that accepted a sloppy step
!> fails on the residual rather than silently moving the pair apart.
subroutine test_mb_agree(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(landscape_muller_brown_model) :: model
   real(wp), allocatable :: x_prfo(:), x_bprfo(:)
   real(wp) :: rmsd, gmax_prfo, gmax_bprfo
   logical :: ok

   call init(env)
   call new_muller_brown_model(model)
   call run_one(env, error, model, p_landscape_ts_prfo, [-0.6_wp, 0.4_wp], x_prfo, &
      & ok, gmax_prfo)
   if (allocated(error)) return
   call run_one(env, error, model, p_landscape_ts_bprfo, [-0.7_wp, 0.45_wp], &
      & x_bprfo, ok, gmax_bprfo)
   if (allocated(error)) return
   ! the searches' own residuals (measured 9.6478e-06 and 1.0395e-06
   ! hartree/bohr): the pair comparison is only meaningful while both are
   ! converged, so this is the gate a different LAPACK is allowed to move
   call check(error, gmax_prfo < 1.0e-4_wp, &
      & "the P-RFO search did not converge tightly enough for the comparison")
   if (allocated(error)) return
   call check(error, gmax_bprfo < 1.0e-4_wp, &
      & "the B-P-RFO search did not converge tightly enough for the comparison")
   if (allocated(error)) return
   rmsd = sqrt(sum((x_prfo - x_bprfo)**2)/real(size(x_prfo), wp))
   write(stdout, '(1x,a,es12.4)') "  prfo vs bprfo RMSD (bohr): ", rmsd
   call check(error, rmsd < 2.0e-5_wp,&
      & "the two methods did not reach the same saddle point")

end subroutine test_mb_agree

!> A mode inside `0 +/- 3 cm^-1` triggers the finite difference tightening.
!>
!> This is the gate of M4 hard constraint 2 (H1 of
!> `29-m4-independent-review.md`): `confirm_saddle` recomputes a Hessian whose
!> spectrum has an internal mode inside the band with `step = 2e-3` bohr and
!> `sccacc = 0.1` before it trusts the sign.  The predicate is what decides
!> that, so it is pinned here: a frequency vector with an internal mode inside
!> the band hits, one whose internal modes are all outside does not, and the
!> projected zero modes at the front are ignored.
subroutine test_ts_curvature_band(error)
   type(error_type), allocatable, intent(out) :: error
   real(wp) :: freq(8)

   ! five projected zero modes, then three internal modes
   freq(1:5) = 0.0_wp
   freq(6) = 0.5_wp          ! inside 0 +/- 3 cm^-1
   freq(7) = 100.0_wp
   freq(8) = 1000.0_wp
   call check(error, landscape_ts_curvature_band_hit(freq, 5, &
      & landscape_ts_curvature_band_cm1), &
      & "a mode inside the band was not detected")
   if (allocated(error)) return

   ! a negative mode inside the band (an almost flat saddle) hits too
   freq(6) = -2.9_wp
   call check(error, landscape_ts_curvature_band_hit(freq, 5, &
      & landscape_ts_curvature_band_cm1), &
      & "an imaginary mode inside the band was not detected")
   if (allocated(error)) return

   ! all internal modes outside the band: no tightening
   freq(6) = 4.0_wp
   call check(error, .not. landscape_ts_curvature_band_hit(freq, 5, &
      & landscape_ts_curvature_band_cm1), &
      & "a mode outside the band triggered the tightening")
   if (allocated(error)) return

   ! the band constant is the documented 3 cm^-1
   call check(error, abs(landscape_ts_curvature_band_cm1 - 3.0_wp) < 1.0e-12_wp)

end subroutine test_ts_curvature_band

!> The curvature band re-decision accepts every spelling of `numerical`.
!>
!> `landscape_ts_can_tighten` is the gate of the tightened recomputation of
!> `confirm_saddle`: only the numerical provider has the `step`/`sccacc`
!> parameters the re-decision sets.  The validator
!> (`landscape_hessian_source_known`) and the provider factory both match the
!> source name case insensitively, so the gate has to as well - a literal
!> comparison accepted `--lts-hessian NUMERICAL` and then reported the
!> resolution limited curvature sign instead of tightening (library entry
!> point R2 of `docs/_simplify-recon/21-stage1-fix-report.md`).
subroutine test_ts_can_tighten_case(error)
   type(error_type), allocatable, intent(out) :: error

   call check(error, landscape_ts_can_tighten('numerical'))
   if (allocated(error)) return
   call check(error, landscape_ts_can_tighten('NUMERICAL'), &
      & "the upper case numerical source cannot tighten")
   if (allocated(error)) return
   call check(error, landscape_ts_can_tighten('NuMeRiCaL'), &
      & "the mixed case numerical source cannot tighten")
   if (allocated(error)) return
   call check(error, .not. landscape_ts_can_tighten('model'), &
      & "the model Hessian claims to support tighter differences")
   if (allocated(error)) return
   call check(error, .not. landscape_ts_can_tighten('MODEL'), &
      & "the upper case model Hessian claims to support tighter differences")
   if (allocated(error)) return
   call check(error, .not. landscape_ts_can_tighten('file'))
   if (allocated(error)) return
   call check(error, .not. landscape_ts_can_tighten('nonsense'))

end subroutine test_ts_can_tighten_case

!> The curvature-band tightening is wired to `landscape_ts_can_tighten`.
!>
!> `confirm_saddle` decides through the public predicate whether a near-zero
!> curvature may be re-decided with tighter finite differences.  The CLI lower
!> cases `--lts-hessian` before it stores it (`main.F90:2966`), so no command
!> line run can tell the normalising predicate from a literal comparison; the
!> gate is therefore driven here through the production subroutine, with a
!> synthetic Hessian that puts a genuine internal mode inside the band.  The
!> two spellings `'NUMERICAL'` and `'numerical'` have to produce the same
!> observable - the `recomputing the Hessian` notice and the `tightened
!> re-decision` line - whereas a call site that compares literally would skip
!> the tightening for the upper case spelling and report the resolution
!> limited curvature (`cannot be re-decided`) instead.  Nothing else about the
!> run depends on the spelling, so the two confirmations are byte identical.
subroutine test_ts_confirm_tighten_case(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(landscape_engine) :: engine
   type(landscape_ts_settings) :: set
   type(landscape_result) :: result
   type(band_hessian_model) :: model
   integer :: unit_lower, unit_upper, nlog_lower, nlog_upper
   character(len=:), allocatable :: lower, upper

   call init(env)
   call init(mol, sym_water, xyz_water)

   set = landscape_ts_settings()
   set%printlevel = 1

   ! lower case spelling: the gate has to tighten
   call new_landscape_engine(engine)
   call engine%set_model(model)
   set%hessian_source = 'numerical'
   open(newunit=unit_lower, status='scratch')
   env%unit = unit_lower
   nlog_lower = env%nLog
   call confirm_saddle(env, mol, chk, engine, set, result)
   nlog_lower = env%nLog - nlog_lower
   lower = read_unit_text(unit_lower)
   close(unit_lower)

   ! upper case spelling: the same gate
   call new_landscape_engine(engine)
   call engine%set_model(model)
   set%hessian_source = 'NUMERICAL'
   open(newunit=unit_upper, status='scratch')
   env%unit = unit_upper
   nlog_upper = env%nLog
   call confirm_saddle(env, mol, chk, engine, set, result)
   nlog_upper = env%nLog - nlog_upper
   upper = read_unit_text(unit_upper)
   close(unit_upper)
   env%unit = stdout

   ! the synthetic Hessian really reaches the curvature band ...
   call check(error, index(lower, 'near-zero curvature:') > 0, &
      & "the synthetic Hessian did not trigger the curvature band")
   if (allocated(error)) return
   ! ... and the numerical source really takes the tightening branch
   call check(error, index(lower, 'tightened re-decision:') > 0, &
      & "the lower case numerical source did not tighten")
   if (allocated(error)) return
   call check(error, index(upper, 'tightened re-decision:') > 0, &
      & "the upper case numerical source did not tighten")
   if (allocated(error)) return
   ! a literal comparison logs the `cannot be re-decided` warning instead, so
   ! the two spellings must not only print the same text but also log the same
   ! number of messages
   call check(error, upper == lower, &
      & "the two spellings of numerical gave different confirmations")
   if (allocated(error)) return
   call check(error, nlog_upper == nlog_lower, &
      & "the two spellings of numerical logged a different number of messages")

end subroutine test_ts_confirm_tighten_case

!> Energy and gradient of the band-Hessian model (never evaluated).
subroutine band_evaluate(self, x, energy, gradient, fail)
   !> Model
   class(band_hessian_model), intent(inout) :: self
   !> Coordinates (bohr)
   real(wp), intent(in) :: x(:)
   !> Energy (hartree)
   real(wp), intent(out) :: energy
   !> Gradient (hartree/bohr)
   real(wp), intent(out) :: gradient(:)
   !> Failure flag
   logical, intent(out) :: fail

   fail = .false.
   energy = 0.0_wp
   gradient = 0.0_wp

end subroutine band_evaluate

!> Constant isotropic Hessian of the band-Hessian model.
subroutine band_hessian(self, x, hess, fail)
   !> Model
   class(band_hessian_model), intent(inout) :: self
   !> Coordinates (bohr)
   real(wp), intent(in) :: x(:)
   !> Hessian (hartree/bohr^2)
   real(wp), intent(out) :: hess(:, :)
   !> Failure flag
   logical, intent(out) :: fail

   integer :: i

   fail = .false.
   hess = 0.0_wp
   do i = 1, size(hess, 1)
      hess(i, i) = self%curvature
   end do

end subroutine band_hessian

!> The forced initial transition vector follows the spelling of `init_mode`.
!>
!> `landscape_ts_init_mode_known` accepts `FORCE`, so the gate that turns the
!> accepted name into the initial vector of P-RFO has to accept it too.  The
!> observable is the report of the selection: only the forced branch prints
!> `transition vector from the largest gradient component`, while the `h0`
!> branch reports the eigenvector it selected.  Two runs of the analytic
!> Mueller-Brown surface with the two spellings therefore have to print the
!> same selection and take the same number of steps (library entry point R2 of
!> `docs/_simplify-recon/21-stage1-fix-report.md`).
subroutine test_prfo_init_mode_case(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(landscape_muller_brown_model) :: model
   type(landscape_ts_settings) :: set
   type(landscape_result) :: result
   integer :: unit, n_steps
   character(len=:), allocatable :: lower, upper
   character(len=*), parameter :: marker = &
      & "transition vector from the largest gradient component:"

   call init(env)
   call new_muller_brown_model(model)

   set = landscape_ts_settings()
   set%method = p_landscape_ts_prfo
   set%printlevel = 1
   set%write_files = .false.
   set%verify_ts = .false.
   set%max_steps = 6
   ! the surface is dimensionless and its curvatures are of order 1 kcal/mol;
   ! the GFN2 oriented hlow would lift every downhill mode
   set%hlow = 1.0e-6_wp

   open(newunit=unit, status='scratch')
   env%unit = unit

   set%init_mode = 'force'
   call landscape_saddle_optimize(env, model, set, [-0.6_wp, 0.4_wp], result)
   n_steps = result%n_steps
   lower = read_unit_text(unit)

   rewind(unit)
   set%init_mode = 'FORCE'
   call landscape_saddle_optimize(env, model, set, [-0.6_wp, 0.4_wp], result)
   upper = read_unit_text(unit)

   close(unit)
   env%unit = stdout

   call check(error, index(lower, marker) > 0, &
      & "the lower case init_mode did not force the transition vector")
   if (allocated(error)) return
   call check(error, index(upper, marker) > 0, &
      & "the upper case init_mode did not force the transition vector")
   if (allocated(error)) return
   call check(error, result%n_steps, n_steps, &
      & "the upper case init_mode changed the search")

end subroutine test_prfo_init_mode_case

!> The dimer initial axis follows the spelling of `n_init`.
!>
!> `landscape_ts_ninit_known` accepts `FORCE`/`RANDOM`, so
!> `landscape_dimer_initial_axis` has to accept them as well: a literal
!> comparison sent the accepted name into the `case default` branch and
!> rejected the run with an error (library entry point R2 of
!> `docs/_simplify-recon/21-stage1-fix-report.md`).  Both spellings must give
!> the same axis - bitwise for the deterministic `force` source and for the
!> seeded `random` one - and the rejection still echoes the name as given.
subroutine test_dimer_axis_case(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   real(wp) :: gradient(4), forced(4), upper_forced(4)
   real(wp) :: random(4), upper_random(4)
   logical :: fail

   call init(env)
   gradient = [0.1_wp, -0.9_wp, 0.2_wp, 0.05_wp]

   call landscape_dimer_initial_axis(env, gradient, 'force', 42, forced, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call landscape_dimer_initial_axis(env, gradient, 'FORCE', 42, &
      & upper_forced, fail)
   call check(error, .not. fail, &
      & "the upper case force axis was rejected")
   if (allocated(error)) return
   call check(error, maxval(abs(upper_forced - forced)), 0.0_wp, thr=0.0_wp, &
      & message="the upper case force axis differs from the lower case one")
   if (allocated(error)) return

   call landscape_dimer_initial_axis(env, gradient, 'random', 42, random, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call landscape_dimer_initial_axis(env, gradient, 'RANDOM', 42, &
      & upper_random, fail)
   call check(error, .not. fail, &
      & "the upper case random axis was rejected")
   if (allocated(error)) return
   call check(error, maxval(abs(upper_random - random)), 0.0_wp, thr=0.0_wp, &
      & message="the upper case random axis differs from the lower case one")
   if (allocated(error)) return

   call landscape_dimer_initial_axis(env, gradient, 'BOGUS', 42, forced, fail)
   call check(error, fail, "an unknown upper case initial axis was accepted")
   if (allocated(error)) return
   call check(error, index(env%log(env%nLog)%message, 'BOGUS') > 0, &
      & "the rejected axis name is not echoed as it was given")

end subroutine test_dimer_axis_case

!> The driver removes the rigid body components from the spectrum and the step.
!>
!> Mueller-Brown is two dimensional (`n_rigid = 0`), so the projection in the
!> driver is a no-op for every other test of this file and removing its call
!> left the whole suite green (M4 of `29-m4-independent-review.md`).  This test
!> uses a three atom (nine dimensional) model with a deliberately rigid
!> gradient component, so the rigid body basis is non-empty:
!>
!>  * `remove_rigid = .true.` keeps the centre of mass and the inertia frame of
!>    every frame (the search only moves along internal coordinates);
!>  * `remove_rigid = .false.` lets the search drift along the rigid force, so
!>    the two trajectories and their final geometries differ.
subroutine test_driver_rigid_projection(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(rigid_drift_model) :: model
   type(landscape_ts_settings) :: set
   type(landscape_result) :: projected, raw_result
   real(wp), allocatable :: x0(:), x_proj(:), x_raw(:), traj_proj(:, :)
   real(wp), allocatable :: traj_raw(:, :)
   real(wp) :: drift_proj, drift_raw, difference

   call init(env)
   x0 = reshape(xyz_water, [3*nat_water])
   call new_rigid_drift_model(model, x0)

   set = landscape_ts_settings()
   set%method = p_landscape_ts_prfo
   set%printlevel = 0
   set%write_files = .false.
   set%verify_ts = .false.
   set%max_steps = 200
   set%hlow = 1.0e-6_wp

   set%remove_rigid = .true.
   call landscape_saddle_optimize(env, model, set, x0, projected, x_proj, &
      & traj_proj)
   if (allocated(error)) return
   set%remove_rigid = .false.
   call landscape_saddle_optimize(env, model, set, x0, raw_result, x_raw, &
      & traj_raw)
   if (allocated(error)) return

   ! the projection keeps the rigid coordinates of the start in every frame
   drift_proj = rigid_frame_drift(traj_proj)
   drift_raw = rigid_frame_drift(traj_raw)
   write(stdout, '(1x,a,es12.4,a,es12.4)') &
      & "  rigid drift (projected / raw): ", drift_proj, " / ", drift_raw
   call check(error, drift_proj < 1.0e-8_wp, &
      & "the projected trajectory still drifts along a rigid body mode")
   if (allocated(error)) return
   call check(error, drift_raw > 1.0e-4_wp, &
      & "the unprojected trajectory did not drift (the model has no rigid force)")
   if (allocated(error)) return

   difference = maxval(abs(x_proj - x_raw))
   call check(error, difference > 1.0e-3_wp, &
      & "the two trajectories reached the same geometry; the projection did "//&
      & "nothing")

end subroutine test_driver_rigid_projection

!> Construct the rigid-drift test model from a geometry.
subroutine new_rigid_drift_model(model, x0)
   !> New model
   type(rigid_drift_model), intent(out) :: model
   !> Expansion point (bohr, length 3N)
   real(wp), intent(in) :: x0(:)

   real(wp), allocatable :: basis(:, :), w(:, :), cand(:)
   real(wp) :: nrm
   integer :: n, n_rigid, n_int, i, k

   model%x0 = x0
   n = size(x0)
   call landscape_rigid_modes(x0, basis, n_rigid)

   ! an orthonormal basis of the internal (complement of the rigid) subspace
   allocate(w(n, 3), source=0.0_wp)
   allocate(cand(n), source=0.0_wp)
   n_int = 0
   do i = 1, n
      if (n_int >= 3) exit
      cand = 0.0_wp
      cand(i) = 1.0_wp
      do k = 1, n_rigid
         cand = cand - dot_product(basis(:, k), cand)*basis(:, k)
      end do
      do k = 1, n_int
         cand = cand - dot_product(w(:, k), cand)*w(:, k)
      end do
      nrm = sqrt(dot_product(cand, cand))
      if (nrm <= 1.0e-8_wp) cycle
      n_int = n_int + 1
      w(:, n_int) = cand/nrm
   end do

   allocate(model%h(n, n), source=0.0_wp)
   call add_outer(model%h, w(:, 1), -0.5_wp)   ! the followed (negative) mode
   call add_outer(model%h, w(:, 2), 0.8_wp)
   call add_outer(model%h, w(:, 3), 1.2_wp)
   do k = 1, n_rigid
      call add_outer(model%h, basis(:, k), 0.05_wp)
   end do

   allocate(model%c(n), source=0.0_wp)
   model%c = 0.02_wp*basis(:, 1) - 0.03_wp*w(:, 1)

end subroutine new_rigid_drift_model

!> Energy and gradient of the rigid-drift model.
subroutine drift_evaluate(self, x, energy, gradient, fail)
   !> Model
   class(rigid_drift_model), intent(inout) :: self
   !> Coordinates (bohr)
   real(wp), intent(in) :: x(:)
   !> Energy (hartree)
   real(wp), intent(out) :: energy
   !> Gradient (hartree/bohr)
   real(wp), intent(out) :: gradient(:)
   !> Failure flag
   logical, intent(out) :: fail

   real(wp), allocatable :: d(:)

   fail = .false.
   energy = 0.0_wp
   gradient = 0.0_wp
   if (.not. allocated(self%h) .or. size(x) /= size(self%x0)) then
      fail = .true.
      return
   end if
   d = x - self%x0
   energy = 0.5_wp*dot_product(d, matmul(self%h, d)) + dot_product(self%c, d)
   gradient = matmul(self%h, d) + self%c

end subroutine drift_evaluate

!> Constant Hessian of the rigid-drift model.
subroutine drift_hessian(self, x, hess, fail)
   !> Model
   class(rigid_drift_model), intent(inout) :: self
   !> Coordinates (bohr)
   real(wp), intent(in) :: x(:)
   !> Hessian (hartree/bohr^2)
   real(wp), intent(out) :: hess(:, :)
   !> Failure flag
   logical, intent(out) :: fail

   fail = .false.
   hess = 0.0_wp
   if (.not. allocated(self%h) .or. size(x) /= size(self%x0)) then
      fail = .true.
      return
   end if
   hess = self%h

end subroutine drift_hessian

!> `a <- a + s v v^T` for a symmetric rank one update.
subroutine add_outer(a, v, s)
   !> Matrix
   real(wp), intent(inout) :: a(:, :)
   !> Vector
   real(wp), intent(in) :: v(:)
   !> Scale
   real(wp), intent(in) :: s

   integer :: i

   do i = 1, size(v)
      a(:, i) = a(:, i) + s*v(i)*v
   end do

end subroutine add_outer

!> Maximum centre-of-mass drift of a trajectory relative to its first frame.
function rigid_frame_drift(traj) result(drift)
   !> Trajectory `(3N, nframes)`
   real(wp), intent(in) :: traj(:, :)
   !> Maximum centre-of-mass displacement (bohr)
   real(wp) :: drift

   real(wp) :: com(3), com0(3)
   integer :: nat, f, a, ic

   nat = size(traj, 1)/3
   com0 = 0.0_wp
   do a = 1, nat
      ic = 3*(a - 1)
      com0 = com0 + traj(ic + 1:ic + 3, 1)
   end do
   com0 = com0/real(nat, wp)

   drift = 0.0_wp
   do f = 2, size(traj, 2)
      com = 0.0_wp
      do a = 1, nat
         ic = 3*(a - 1)
         com = com + traj(ic + 1:ic + 3, f)
      end do
      com = com/real(nat, wp)
      drift = max(drift, sqrt(sum((com - com0)**2)))
   end do

end function rigid_frame_drift

!> Run one search on Mueller-Brown and check the saddle point.
subroutine run_mb(error, method, xstart)
   type(error_type), allocatable, intent(out) :: error
   !> Method identifier
   integer, intent(in) :: method
   !> Starting point
   real(wp), intent(in) :: xstart(:)

   type(TEnvironment) :: env
   type(landscape_muller_brown_model) :: model
   real(wp), allocatable :: xfinal(:)
   logical :: ok

   call init(env)
   call new_muller_brown_model(model)
   call run_one(env, error, model, method, xstart, xfinal, ok)
   if (allocated(error)) return
   call check(error, ok, "the search did not converge")
   if (allocated(error)) return
   call check(error, maxval(abs(xfinal - mb_saddle)) < 1.0e-3_wp,&
      & "the search did not reach the Mueller-Brown saddle point")

end subroutine run_mb

!> Run one transition state search with analytic Mueller-Brown data.
subroutine run_one(env, error, model, method, xstart, xfinal, converged, gmax)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Failure report
   type(error_type), allocatable, intent(out) :: error
   !> Potential
   class(landscape_energy_model), intent(inout) :: model
   !> Method identifier
   integer, intent(in) :: method
   !> Starting point
   real(wp), intent(in) :: xstart(:)
   !> Converged geometry
   real(wp), allocatable, intent(out) :: xfinal(:)
   !> Convergence flag
   logical, intent(out) :: converged
   !> Largest gradient component at the converged point (hartree/bohr)
   real(wp), intent(out), optional :: gmax

   type(landscape_ts_settings) :: set
   type(landscape_result) :: result

   set = landscape_ts_settings()
   set%method = method
   set%printlevel = 0
   set%write_files = .false.
   set%max_steps = 400
   set%verify_ts = .false.
   ! the surface is dimensionless and its curvatures are of order 1 kcal/mol;
   ! the GFN2 oriented hlow would lift every downhill mode
   set%hlow = 1.0e-6_wp

   call landscape_saddle_optimize(env, model, set, xstart, result, xfinal)
   converged = result%converged
   if (present(gmax)) gmax = result%gmax
   write(stdout, '(1x,a,a,a,i6,a,l1,a,es12.4,a,es12.4)') "  ", &
      & landscape_ts_method_name(method), " steps ", result%n_steps, &
      & " converged ", result%converged, ", gmax ", result%gmax, &
      & ", dist ", maxval(abs(xfinal - mb_saddle))
   ! full precision converged point: the PLAN experiment 5 cross-validation
   ! report compares it with the analytic Mueller-Brown saddle point
   write(stdout, '(1x,a,2(1x,f20.14))') "  converged point ", xfinal

end subroutine run_one

!> The curvature of Eq. 2 is the energy curvature along the dimer axis.
!>
!> The image labels matter: image 1 sits at `R + delta*N` and image 2 at
!> `R - delta*N`, which is what makes `C = N^T H N` negative for a
!> reaction-coordinate axis, i.e. what selects the branch of Eq. 21.  The
!> quantitative check is the force difference of the analytic Mueller-Brown
!> Hessian, which pins both the sign and the factor `2*delta`.
subroutine test_dimer_curvature(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(landscape_muller_brown_model) :: model
   real(wp), allocatable :: f_plus(:), f_minus(:), gradient(:), hess(:, :)
   real(wp) :: delta, energy, curvature, projected
   integer :: ireference
   real(wp) :: axis(2), x(2)
   logical :: fail

   call init(env)
   delta = 1.0e-2_wp
   allocate(f_plus(2), f_minus(2), gradient(2), hess(2, 2))

   ! H = diag(-2, 1) in force form: F = -H x at the images x = +/- delta e_i
   axis = [1.0_wp, 0.0_wp]
   f_plus = [2.0_wp*delta, 0.0_wp]
   f_minus = [-2.0_wp*delta, 0.0_wp]
   call check(error, abs(landscape_dimer_curvature(f_plus, f_minus, axis, delta) &
      & + 2.0_wp) < 1.0e-12_wp, "the curvature along the negative mode is not -2")
   if (allocated(error)) return
   axis = [0.0_wp, 1.0_wp]
   f_plus = [0.0_wp, -delta]
   f_minus = [0.0_wp, delta]
   call check(error, abs(landscape_dimer_curvature(f_plus, f_minus, axis, delta) &
      & - 1.0_wp) < 1.0e-12_wp, "the curvature along the positive mode is not 1")
   if (allocated(error)) return

   ! quantitative check against the analytic Hessian of Mueller-Brown
   call new_muller_brown_model(model)
   x = [-0.6_wp, 0.4_wp]
   do ireference = 1, 2
      if (ireference > 1) then
         axis = [0.0_wp, 1.0_wp]
      else
         axis = [1.0_wp, 0.0_wp]
      end if
      call model%evaluate(x + delta*axis, energy, gradient, fail)
      if (allocated(error)) return
      f_plus = -gradient
      call model%evaluate(x - delta*axis, energy, gradient, fail)
      f_minus = -gradient
      call model%hessian(x, hess, fail)
      curvature = landscape_dimer_curvature(f_plus, f_minus, axis, delta)
      projected = dot_product(axis, matmul(hess, axis))
      write(stdout, '(1x,a,es12.4,a,es12.4)') &
         & "  dimer curvature (force/analytic): ", curvature, " / ", projected
      ! the central difference carries an O(delta^2) truncation error of the
      ! fourth derivative; measured at the default delta it is 0.07 % of the
      ! curvature along e1 and 0.09 % along e2, so the tolerance is 0.2 %
      call check(error, abs(curvature - projected) < 2.0e-3_wp*abs(projected), &
         & "the finite difference curvature differs from N^T H N")
      if (allocated(error)) return
   end do

end subroutine test_dimer_curvature

!> Both branches of the translation force of Eq. 21.
!>
!> The parallel component of the midpoint force is *reversed* in both branches
!> (the axis is climbed); the branches differ in the perpendicular subspace,
!> which is relaxed only for a negative curvature.  A sign flip of the branch
!> condition therefore moves the perpendicular component.
subroutine test_dimer_translation_force(error)
   type(error_type), allocatable, intent(out) :: error
   real(wp) :: force(2), axis(2), f_dagger(2)
   real(wp) :: parallel, perpendicular

   force = [0.4_wp, 0.3_wp]
   axis = [1.0_wp, 0.0_wp]

   ! negative curvature (reaction coordinate): ascend along the axis, relax
   ! the perpendicular subspace
   call landscape_dimer_translation_force(force, axis, -1.0e-2_wp, f_dagger)
   parallel = dot_product(f_dagger, axis)
   perpendicular = f_dagger(2)
   call check(error, abs(parallel + 0.4_wp) < 1.0e-12_wp, &
      & "the parallel force is not reversed for C < 0")
   if (allocated(error)) return
   call check(error, abs(perpendicular - 0.3_wp) < 1.0e-12_wp, &
      & "the perpendicular force is not kept for C < 0")
   if (allocated(error)) return

   ! positive curvature: the perpendicular relaxation is dropped
   call landscape_dimer_translation_force(force, axis, 1.0e-2_wp, f_dagger)
   call check(error, abs(dot_product(f_dagger, axis) + 0.4_wp) < 1.0e-12_wp, &
      & "the parallel force is not reversed for C > 0")
   if (allocated(error)) return
   call check(error, abs(f_dagger(2)) < 1.0e-12_wp, &
      & "the perpendicular force is kept for C > 0")

end subroutine test_dimer_translation_force

!> Eq. 13 recovers the root of an exact two term Fourier rotational force.
!>
!> The rotational force is `sin(2(theta* - theta))` here, i.e. a pure two term
!> Fourier series (the model the paper assumes) with a **positive** value at
!> `theta = 0`, like the magnitude of the rotational force of a real dimer.  The
!> angle that nulls it is known exactly.  Eq. 13 finds it up to the `O(dtheta)`
!> error of the finite difference derivative; the printed reference values show
!> what the earlier Eq. 6 of the paper would return instead.
subroutine test_dimer_rotation_angle(error)
   type(error_type), allocatable, intent(out) :: error
   real(wp) :: delta_theta, f0, f0_prime, theta, theta6
   real(wp) :: target, angle, probe

   ! small misalignment: the linear regime, where Eq. 6 differs by dtheta/2
   target = 1.0e-3_wp
   delta_theta = 1.0e-4_wp
   f0 = fourier_rotation_force(0.0_wp, target)
   probe = fourier_rotation_force(delta_theta, target)
   f0_prime = (probe - f0)/delta_theta
   theta = landscape_dimer_rotation_angle(f0, f0_prime)
   theta6 = -(f0 + probe)/(2.0_wp*f0_prime)
   write(stdout, '(1x,a,es12.4,a,es12.4,a,es12.4)') &
      & "  dimer rotation angle (Eq. 13 / Eq. 6 / exact): ", theta, " / ", &
      & theta6, " / ", target
   call check(error, abs(theta - target) < 1.0e-7_wp, &
      & "Eq. 13 does not recover the small angle Fourier root")
   if (allocated(error)) return
   ! the linear interpolation of Eq. 6 cannot pass that tolerance
   call check(error, abs(theta6 - target) > 1.0e-6_wp, &
      & "the Eq. 6 reference has to differ measurably")
   if (allocated(error)) return

   ! large misalignment: the arctan of Eq. 13 saturates, Eq. 6 does not
   target = 0.3_wp
   f0 = fourier_rotation_force(0.0_wp, target)
   probe = fourier_rotation_force(delta_theta, target)
   f0_prime = (probe - f0)/delta_theta
   theta = landscape_dimer_rotation_angle(f0, f0_prime)
   theta6 = -(f0 + probe)/(2.0_wp*f0_prime)
   write(stdout, '(1x,a,es12.4,a,es12.4,a,es12.4)') &
      & "  dimer rotation angle (Eq. 13 / Eq. 6 / exact): ", theta, " / ", &
      & theta6, " / ", target
   call check(error, abs(theta - target) < 1.0e-3_wp, &
      & "Eq. 13 does not recover the large angle Fourier root")
   if (allocated(error)) return
   call check(error, abs(theta6 - target) > 1.0e-2_wp, &
      & "the Eq. 6 reference has to differ measurably")
   if (allocated(error)) return

   ! a vanishing rotational force needs no rotation
   call check(error, abs(landscape_dimer_rotation_angle(0.0_wp, 1.0_wp)) &
      & < 1.0e-15_wp)
   ! a large misalignment: the root has to be the *first positive* one, i.e.
   ! the minimum of the curvature along the probe direction (see the function
   ! documentation); the principal branch of the arctan would return the
   ! maximum, 90 degrees away
   call check(error, abs(landscape_dimer_rotation_angle(0.0136_wp, &
      & 0.0148_wp) - 1.0352_wp) < 1.0e-3_wp, &
      & "the rotation angle is not the minimising root")

end subroutine test_dimer_rotation_angle

!> The rotation lowers the torque of a misaligned dimer on Mueller-Brown.
subroutine test_dimer_rotation_alignment(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(landscape_muller_brown_model) :: model
   type(landscape_ts_settings) :: set
   type(landscape_transition_mode) :: mode
   real(wp), allocatable :: rigid(:, :), f_plus(:), f_minus(:), q(:)
   real(wp), allocatable :: energy(:), gradient(:)
   real(wp) :: x(2), torque_before, torque_after, magnitude, rms
   logical :: fail, converged
   integer :: n_rigid, iterations, n_eval

   call init(env)
   call new_muller_brown_model(model)
   x = [-0.6_wp, 0.4_wp]
   allocate(q(2), energy(2), gradient(2))
   call landscape_rigid_modes(x, rigid, n_rigid)

   ! the torque of the starting axis
   mode%vector = [1.0_wp, 0.0_wp]
   allocate(f_plus(2), f_minus(2))
   call dimer_images_in_test(model, x, mode%vector, 1.8897e-2_wp, f_plus, &
      & f_minus)
   call landscape_dimer_rotation_direction(f_plus, f_minus, mode%vector, q, &
      & magnitude, torque_before)

   set = landscape_ts_settings()
   set%printlevel = 0
   set%rot_alpha = 1.0_wp
   set%rot_max_iter = 5
   call landscape_dimer_rotate(env, model, x, set, mode, rigid, n_rigid, &
      & f_plus, f_minus, iterations, torque_after, magnitude, converged, &
      & n_eval, fail)
   write(stdout, '(1x,a,es12.4,a,es12.4,a,l1)') &
      & "  dimer torque rms before/after rotation: ", torque_before, " / ", &
      & torque_after, ", converged=", converged
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, torque_after < 0.2_wp*torque_before, &
      & "the rotation did not lower the rotational force")
   if (allocated(error)) return
   call check(error, abs(sum(mode%vector**2) - 1.0_wp) < 1.0e-12_wp, &
      & "the rotated axis is not normalised")

end subroutine test_dimer_rotation_alignment

!> Exhausting the rotation budget above the torque threshold warns.
!>
!> An unreachable threshold (`rot_f_rms_th = 0`) is the pathological case of
!> the VASP diagnostic (the torque does not fall because the forces are not
!> accurate enough): the rotation has to stop after `rot_max_iter` rotations,
!> report the failure and emit the documented warning.
subroutine test_dimer_rotation_warning(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(landscape_muller_brown_model) :: model
   type(landscape_ts_settings) :: set
   type(landscape_transition_mode) :: mode
   real(wp), allocatable :: rigid(:, :), f_plus(:), f_minus(:)
   real(wp) :: x(2), torque_rms, torque_max
   logical :: fail, converged
   integer :: n_rigid, iterations, n_eval, nlog

   call init(env)
   call new_muller_brown_model(model)
   x = [-0.6_wp, 0.4_wp]
   call landscape_rigid_modes(x, rigid, n_rigid)
   mode%vector = [1.0_wp, 0.0_wp]

   set = landscape_ts_settings()
   set%printlevel = 0
   set%rot_alpha = 1.0_wp
   set%rot_max_iter = 2
   set%rot_f_rms_th = -1.0_wp
   nlog = env%nLog
   call landscape_dimer_rotate(env, model, x, set, mode, rigid, n_rigid, &
      & f_plus, f_minus, iterations, torque_rms, torque_max, converged, &
      & n_eval, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, .not. converged)
   if (allocated(error)) return
   call check(error, iterations == 2, &
      & "the rotation did not run the full budget")
   if (allocated(error)) return
   call check(error, env%nLog > nlog, &
      & "the exhausted rotation did not warn")
   if (allocated(error)) return
   call check(error, index(env%log(env%nLog)%message, &
      & "dimer rotation did not converge") > 0, &
      & "the warning is not the documented VASP diagnostic")

end subroutine test_dimer_rotation_warning

!> The quick-min update: increment, momentum, cap and uphill reset.
subroutine test_dimer_velocity(error)
   type(error_type), allocatable, intent(out) :: error
   real(wp) :: force(2), velocity(2), step(2)

   force = [1.0_wp, 0.0_wp]

   ! a zero incoming velocity gives exactly the increment
   velocity = 0.0_wp
   call landscape_dimer_velocity_update(force, velocity, 0.1_wp, 0.25_wp, step)
   call check(error, abs(step(1) - 0.1_wp) < 1.0e-14_wp)
   if (allocated(error)) return

   ! a parallel velocity accumulates (the momentum of quick-min)
   velocity = [0.1_wp, 0.0_wp]
   call landscape_dimer_velocity_update(force, velocity, 0.1_wp, 0.25_wp, step)
   call check(error, abs(step(1) - 0.2_wp) < 1.0e-14_wp, &
      & "the parallel velocity did not accumulate")
   if (allocated(error)) return

   ! the momentum is capped at step_max
   velocity = [0.3_wp, 0.0_wp]
   call landscape_dimer_velocity_update(force, velocity, 0.1_wp, 0.25_wp, step)
   call check(error, abs(step(1) - 0.25_wp) < 1.0e-14_wp, &
      & "the displacement exceeds step_max")
   if (allocated(error)) return

   ! an uphill velocity is discarded, so the step never opposes the force
   velocity = [-0.2_wp, 0.0_wp]
   call landscape_dimer_velocity_update(force, velocity, 0.05_wp, 0.25_wp, step)
   call check(error, all(abs(step) < 1.0e-14_wp), &
      & "an uphill velocity was not discarded")
   if (allocated(error)) return
   call check(error, all(abs(velocity) < 1.0e-14_wp))
   if (allocated(error)) return

   ! a vanishing force gives no step
   velocity = [0.1_wp, 0.0_wp]
   call landscape_dimer_velocity_update([0.0_wp, 0.0_wp], velocity, 0.1_wp, &
      & 0.25_wp, step)
   call check(error, all(abs(step) < 1.0e-14_wp))

end subroutine test_dimer_velocity

!> The mobility limiter of the translation (the consumer of `dR`).
!>
!> A non-zero curvature caps the displacement at the harmonic length scale
!> `1/|kappa|` of the direction that is followed - for a **positive** curvature
!> (stiff direction) as well as for a **negative** one (a direction that is
!> being climbed, which would otherwise grow without limit, M8 of
!> `29-m4-independent-review.md`); a flat direction (below the floor) leaves
!> the propagation mobility alone.
subroutine test_dimer_mobility(error)
   type(error_type), allocatable, intent(out) :: error
   real(wp) :: mobility, tau

   mobility = 2.0_wp
   ! flat curvature: unchanged
   call check(error, abs(landscape_dimer_mobility(mobility, 1.0e-9_wp) &
      & - mobility) < 1.0e-15_wp)
   if (allocated(error)) return
   ! negative curvature: bounded by 1/|kappa| (M8)
   call check(error, abs(landscape_dimer_mobility(mobility, -2.0_wp) &
      & - 0.5_wp) < 1.0e-12_wp, &
      & "a negative curvature is not bounded by 1/|kappa|")
   if (allocated(error)) return
   ! stiff direction: the harmonic length scale binds (measured on GFN2 HCN
   ! the first translation reports mobility 1.2070 at kappa 0.8285 Eh/bohr^2)
   tau = landscape_dimer_mobility(mobility, 0.8285_wp)
   call check(error, abs(tau - 1.0_wp/0.8285_wp) < 1.0e-12_wp, &
      & "the harmonic length scale does not bound the mobility")
   if (allocated(error)) return
   ! soft direction: the propagation mobility is the smaller one
   call check(error, abs(landscape_dimer_mobility(mobility, 0.1_wp) &
      & - mobility) < 1.0e-15_wp)

end subroutine test_dimer_mobility

!> The initial axis: the deterministic `force` source and the seeded random one.
!>
!> PLAN N10 requires an explicit seed: the axis has to be a pure function of
!> `settings%seed`, so two draws with the same seed are identical and a
!> different seed gives a different axis.
subroutine test_dimer_initial_axis(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   real(wp) :: gradient(4), axis(4), first(4), second(4), other(4)
   real(wp) :: difference
   logical :: fail
   integer :: nlog

   call init(env)
   gradient = [0.1_wp, -0.9_wp, 0.2_wp, 0.05_wp]

   call landscape_dimer_initial_axis(env, gradient, 'force', 42, axis, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   ! the axis is a line, so the print/sign convention is irrelevant here
   call check(error, abs(abs(axis(2)) - 1.0_wp) < 1.0e-14_wp, &
      & "the deterministic axis is not the largest gradient component")
   if (allocated(error)) return
   call check(error, abs(sum(axis**2) - 1.0_wp) < 1.0e-14_wp)
   if (allocated(error)) return

   call landscape_dimer_initial_axis(env, gradient, 'random', 42, first, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call landscape_dimer_initial_axis(env, gradient, 'random', 42, second, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, all(abs(first - second) < 1.0e-15_wp), &
      & "the random axis is not reproducible for one seed")
   if (allocated(error)) return
   call check(error, abs(sum(first**2) - 1.0_wp) < 1.0e-14_wp)
   if (allocated(error)) return

   call landscape_dimer_initial_axis(env, gradient, 'random', 7, other, fail)
   difference = maxval(abs(first - other))
   write(stdout, '(1x,a,es12.4)') &
      & "  dimer random axis: seed 42 vs seed 7, max difference ", difference
   call check(error, difference > 1.0e-3_wp, &
      & "the axis does not depend on the seed")
   if (allocated(error)) return

   nlog = env%nLog
   call landscape_dimer_initial_axis(env, gradient, 'bogus', 42, axis, fail)
   call check(error, fail, "an unknown initial axis must be rejected")
   if (allocated(error)) return
   call check(error, env%nLog > nlog)

end subroutine test_dimer_initial_axis

!> Unusable dimer parameters abort the search instead of producing a step.
!>
!> Every length scale of the dimer (the separation, the finite difference step,
!> the probe angle and the two step lengths) enters a division, a normalisation
!> or an exponential, so a value that is not positive has to be rejected before
!> the first evaluation.
subroutine test_dimer_invalid(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(landscape_muller_brown_model) :: model
   type(landscape_ts_settings) :: set
   type(landscape_result) :: result
   real(wp), allocatable :: xfinal(:)
   integer :: nlog

   call init(env)
   call new_muller_brown_model(model)

   set = landscape_ts_settings()
   set%method = p_landscape_ts_dimer
   set%printlevel = 0
   set%write_files = .false.
   set%verify_ts = .false.
   set%delta = 0.0_wp
   nlog = env%nLog
   call landscape_dimer_optimize(env, model, set, [-0.75_wp, 0.5_wp], result, &
      & xfinal)
   call check(error, .not. result%converged, &
      & "a zero dimer separation must not start a search")
   if (allocated(error)) return
   call check(error, env%nLog > nlog)
   if (allocated(error)) return
   call check(error, index(env%log(env%nLog)%message, &
      & "must be positive") > 0, "the rejection is not reported")
   if (allocated(error)) return

   set%delta = 1.8897e-2_wp
   set%step_max = -1.0_wp
   call landscape_dimer_optimize(env, model, set, [-0.75_wp, 0.5_wp], result, &
      & xfinal)
   call check(error, .not. result%converged, &
      & "a negative step_max must not start a search")

end subroutine test_dimer_invalid

!> The dimer converges to the Mueller-Brown saddle point without a Hessian.
subroutine test_dimer_mb(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(landscape_muller_brown_model) :: model
   type(landscape_result) :: result
   real(wp), allocatable :: x_dimer(:), x_prfo(:)
   real(wp) :: rmsd, gmax_prfo
   logical :: ok

   call init(env)
   call new_muller_brown_model(model)
   call run_dimer(env, model, [-0.75_wp, 0.5_wp], x_dimer, result, 2, 'random')
   call run_one(env, error, model, p_landscape_ts_prfo, [-0.75_wp, 0.5_wp], &
      & x_prfo, ok, gmax_prfo)
   if (allocated(error)) return
   rmsd = sqrt(sum((x_dimer - x_prfo)**2)/real(size(x_dimer), wp))
   write(stdout, '(1x,a,i4,a,l1,a,es12.4,a,es12.4,a,es12.4)') &
      & "  dimer steps ", result%n_steps, " converged ", result%converged, &
      & ", gmax ", result%gmax, ", dist ", &
      & maxval(abs(x_dimer - mb_saddle)), ", rmsd(prfo) ", rmsd
   call check(error, result%converged, "the dimer search did not converge")
   if (allocated(error)) return
   call check(error, maxval(abs(x_dimer - mb_saddle)) < 1.0e-3_wp, &
      & "the dimer did not reach the Mueller-Brown saddle point")
   if (allocated(error)) return
   ! delivered 7.7397e-05 bohr (AGENTS.md section 2); 2e-4 is 2.6x above it
   ! and 5x tighter than the 1e-3 tripwire of the other searches.  Measured
   ! identically on gfortran 12.2.0 and on the gating CI's gfortran 14.4.0
   ! (both 7.7397e-05 bohr, the converged points equal to 14 digits), so the
   ! datum is not compiler specific; the two searches' own residuals below
   ! carry the robustness of the comparison.
   call check(error, gmax_prfo < 1.0e-4_wp, &
      & "the P-RFO reference did not converge tightly enough for the comparison")
   if (allocated(error)) return
   ! measured 1.0413e-04 hartree/bohr for the dimer
   call check(error, result%gmax < 1.0e-3_wp, &
      & "the dimer did not converge tightly enough for the comparison")
   if (allocated(error)) return
   call check(error, rmsd < 2.0e-4_wp, &
      & "the dimer and P-RFO results differ")

end subroutine test_dimer_mb

!> The dimer never asks for a Hessian.
!>
!> The potential is a Mueller-Brown surface whose `hessian` is unavailable; a
!> single call would abort the search, so the converged run itself is the
!> assertion (PLAN section 9.5: the dimer is the method of choice when a
!> Hessian is too expensive).
subroutine test_dimer_no_hessian(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(no_hessian_model) :: model
   type(landscape_result) :: result
   real(wp), allocatable :: x_dimer(:)

   call init(env)
   call run_dimer(env, model, [-0.75_wp, 0.5_wp], x_dimer, result, 2, 'random')
   write(stdout, '(1x,a,i4,a,i0,a,l1)') "  dimer Hessian requests ", &
      & model%n_hess, " (steps ", result%n_steps, ", converged ", &
      & result%converged, ")"
   call check(error, model%n_hess == 0, &
      & "the dimer requested a Hessian")
   if (allocated(error)) return
   call check(error, result%converged, &
      & "the Hessian-free search did not converge")

end subroutine test_dimer_no_hessian

!> Two runs with the same seed are bitwise identical.
subroutine test_dimer_reproducible(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(landscape_muller_brown_model) :: model
   type(landscape_result) :: first_result, second_result, other_result
   real(wp), allocatable :: first(:), second(:), other(:)
   real(wp) :: difference

   call init(env)
   call new_muller_brown_model(model)
   call run_dimer(env, model, [-0.75_wp, 0.5_wp], first, first_result, 2, &
      & 'random')
   call run_dimer(env, model, [-0.75_wp, 0.5_wp], second, second_result, 2, &
      & 'random')
   call check(error, all(first == second), &
      & "two runs with the same seed differ")
   if (allocated(error)) return
   call check(error, first_result%n_steps == second_result%n_steps)
   if (allocated(error)) return

   call run_dimer(env, model, [-0.75_wp, 0.5_wp], other, other_result, 7, &
      & 'random')
   difference = maxval(abs(first - other))
   write(stdout, '(1x,a,es12.4,a,l1,a,l1)') &
      & "  dimer seed 2 vs 7: geometry difference ", difference, &
      & ", converged both ", first_result%converged, &
      & " / ", other_result%converged
   call check(error, first_result%converged .and. other_result%converged)
   if (allocated(error)) return
   ! A different draw gives a different initial axis but, on this surface, the
   ! same saddle point: measured, the two displacements agree to 1.0e-6 bohr
   ! (the seeded draw itself is pinned by `dimer-initial-axis`).
   call check(error, difference < 1.0e-4_wp, &
      & "two seeds reached different stationary points")

end subroutine test_dimer_reproducible

!> Rotational force of the exact two term Fourier model used by the Eq. 13 test.
!>
!> `F(theta) = sin(2 (theta* - theta))` is the form the paper assumes, with the
!> root `theta*` known analytically.
pure function fourier_rotation_force(theta, target) result(force)
   !> Rotation angle (rad)
   real(wp), intent(in) :: theta
   !> Exact root of the rotational force (rad)
   real(wp), intent(in) :: target
   !> Rotational force (arbitrary unit)
   real(wp) :: force

   force = sin(2.0_wp*(target - theta))

end function fourier_rotation_force

!> Image forces of a dimer in a unit test (no failure handling needed).
subroutine dimer_images_in_test(model, x, axis, delta, f_plus, f_minus)
   !> Potential
   class(landscape_muller_brown_model), intent(inout) :: model
   !> Midpoint geometry (bohr)
   real(wp), intent(in) :: x(:)
   !> Unit dimer axis
   real(wp), intent(in) :: axis(:)
   !> Centre-to-image distance (bohr)
   real(wp), intent(in) :: delta
   !> Force at `R + delta*axis` (hartree/bohr)
   real(wp), intent(out) :: f_plus(:)
   !> Force at `R - delta*axis` (hartree/bohr)
   real(wp), intent(out) :: f_minus(:)

   real(wp), allocatable :: gradient(:)
   real(wp) :: energy
   logical :: fail

   allocate(gradient(size(x)))
   call model%evaluate(x + delta*axis, energy, gradient, fail)
   f_plus = -gradient
   call model%evaluate(x - delta*axis, energy, gradient, fail)
   f_minus = -gradient

end subroutine dimer_images_in_test

!> Run one dimer search of the analytic test surfaces.
!>
!> The settings are the defaults of `--lts` with the method, the seed and the
!> starting axis overridden; warnings are echoed once so that the diagnostic of
!> the rotation shows up in the test output.
subroutine run_dimer(env, model, xstart, xfinal, result, seed, ninit)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Potential
   class(landscape_energy_model), intent(inout) :: model
   !> Starting point (bohr)
   real(wp), intent(in) :: xstart(:)
   !> Converged geometry (bohr)
   real(wp), allocatable, intent(out) :: xfinal(:)
   !> Search result
   type(landscape_result), intent(out) :: result
   !> Seed of the random initial axis
   integer, intent(in) :: seed
   !> Source of the initial axis
   character(len=*), intent(in) :: ninit

   type(landscape_ts_settings) :: set
   integer :: nlog

   set = landscape_ts_settings()
   set%method = p_landscape_ts_dimer
   set%printlevel = 0
   set%write_files = .false.
   set%verify_ts = .false.
   set%max_steps = 400
   set%seed = seed
   set%n_init = ninit
   ! The Mueller-Brown surface is dimensionless and its energy scale is
   ! `kcaltoau`; the default step lengths (0.28 - 0.38 bohr, calibrated for a
   ! GFN2 geometry) are a third of the whole surface and make the translation
   ! jump past the saddle.  Both lengths are therefore calibrated here; the
   ! algorithm is untouched (measured: 23 steps to 1.4e-4 bohr of the saddle
   ! with these values, no convergence in 400 steps with the defaults).
   set%step0 = 5.0e-2_wp
   set%step_max = 5.0e-2_wp
   nlog = env%nLog
   call landscape_dimer_optimize(env, model, set, xstart, result, xfinal)
   if (env%nLog > nlog) then
      write(stdout, '(1x,a)') "  dimer warning: "//env%log(env%nLog)%message
   end if
   ! full precision converged point for the PLAN experiment 5 cross-validation
   if (allocated(xfinal)) write(stdout, '(1x,a,2(1x,f20.14))') &
      & "  dimer converged point ", xfinal

end subroutine run_dimer

!> Identity matrix.
pure subroutine identity(matrix)
   !> Matrix to fill
   real(wp), intent(out) :: matrix(:, :)

   integer :: i

   matrix = 0.0_wp
   do i = 1, min(size(matrix, 1), size(matrix, 2))
      matrix(i, i) = 1.0_wp
   end do

end subroutine identity

!> Hessian of the Hessian-free test surface, always unavailable.
subroutine no_hessian_hessian(self, x, hess, fail)
   !> Model
   class(no_hessian_model), intent(inout) :: self
   !> Coordinates (bohr)
   real(wp), intent(in) :: x(:)
   !> Hessian (hartree/bohr^2), never filled
   real(wp), intent(out) :: hess(:, :)
   !> Failure flag, always set
   logical, intent(out) :: fail

   self%n_hess = self%n_hess + 1
   hess = 0.0_wp
   fail = .true.

end subroutine no_hessian_hessian

!> Swap two columns of a matrix.
pure subroutine swap_columns(matrix, a, b)
   !> Matrix
   real(wp), intent(inout) :: matrix(:, :)
   !> Column indices
   integer, intent(in) :: a, b

   real(wp) :: buffer(size(matrix, 1))

   buffer = matrix(:, a)
   matrix(:, a) = matrix(:, b)
   matrix(:, b) = buffer

end subroutine swap_columns

!> Change of the two O-H distances of the water test geometry under a displaced
!> geometry; zero for a rigid body displacement.
pure function internal_distance_change(xref, x) result(change)
   !> Reference geometry (bohr)
   real(wp), intent(in) :: xref(:)
   !> Displaced geometry (bohr)
   real(wp), intent(in) :: x(:)
   !> Maximum absolute distance change (bohr)
   real(wp) :: change

   integer :: i

   change = 0.0_wp
   do i = 2, 3
      change = max(change, abs(norm2(x(3*i - 2:3*i) - x(1:3)) &
         & - norm2(xref(3*i - 2:3*i) - xref(1:3))))
   end do

end function internal_distance_change

!> Every record written to a scratch unit, separated by newlines.
!>
!> The startup banner of a task is a property of the printed output, so the
!> only way to assert it from a unit test is to point `env%unit` at a scratch
!> file and read it back.  The unit is rewound first, so the helper also works
!> after a previous read reached the endfile record (a rewrite then needs the
!> rewind anyway).
function read_unit_text(unit) result(text)
   !> Scratch unit written by the task under test
   integer, intent(in) :: unit
   !> Records of the unit, joined by newlines
   character(len=:), allocatable :: text

   character(len=256) :: line
   integer :: stat

   text = ''
   rewind(unit)
   do
      read(unit, '(a)', iostat=stat) line
      if (stat /= 0) exit
      text = text//trim(line)//new_line('a')
   end do

end function read_unit_text

end module test_landscape_saddle
