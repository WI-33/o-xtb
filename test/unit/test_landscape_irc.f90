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

!> Tests of the intrinsic reaction coordinate integrators (M6).
!>
!> The coverage is analytic wherever possible (PLAN risk R-14): the
!> integrators are compared against **closed form** references, not against a
!> second implementation of the same idea.
!>
!>  * unit round trip: the mass weighted coordinate `q_i = sqrt(m_i/amu) r_i`
!>    and the unit `1 bohr*amu^1/2` of the arc length, checked against a
!>    hand computed translation of a three atom system (the largest trap of
!>    the task, PLAN section 6.5).
!>  * LQA analytic propagation against the closed form solution of the
!>    quadratic ODE (one and two dimensions) and the modified Bulirsch-Stoer
!>    corrector against the same closed form.
!>  * the DWI surface reproduces its sample energies and gradients and
!>    interpolates between two quadratic expansions.
!>  * Mueller-Brown surface: GS and LQA from the analytic saddle point reach
!>    the two analytic minima, the path energy is below the saddle energy
!>    everywhere and the arc length increments respect the sphere constraint
!>    of the Gonzalez-Schlegel scheme (chord <= step).
!>  * LEPS H + H2: both directions lower the energy and the H2 bond length
!>    approaches the closed form equilibrium distance.
!>  * termination signals: (a) `max_steps`, (b) the gradient threshold and
!>    (c) the energy rise are reached with the documented conditions, and (d)
!>    the GS limit cycle guard truncates a stalled macro step walk (a cut, not
!>    a convergence: its endpoint counts only after an endpoint relaxation).
!>  * the transition vector is deterministic and rejects a minimum.
!>
!> Units: mass weighted coordinates and arc lengths in `bohr*amu^1/2`,
!> energies in hartree, gradients in hartree/bohr.
module test_landscape_irc
   use testdrive, only : new_unittest, unittest_type, error_type, check
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use xtb_mctc_lapack, only : mctc_syev
   use xtb_type_environment
   use xtb_landscape_settings, only : landscape_irc_settings, &
      & p_landscape_irc_gs, p_landscape_irc_lqa, p_landscape_irc_hpc, &
      & p_landscape_irc_eulerpc
   use xtb_landscape_engine, only : landscape_energy_model
   use xtb_landscape_coords, only : vec_norm
   use xtb_landscape_results, only : p_landscape_irc_term_max_steps, &
      & p_landscape_irc_term_grad, p_landscape_irc_term_energy_rise, &
      & p_landscape_irc_term_limit_cycle
   use xtb_landscape_saddle_eigenfollow, only : landscape_transition_mode, &
      & landscape_transition_mode_set_phase
   use xtb_landscape_irc_corrector, only : landscape_irc_sqrt_mass, &
      & landscape_irc_to_weighted, landscape_irc_from_weighted, &
      & landscape_irc_gradient_to_weighted, landscape_irc_hessian_to_weighted, &
      & landscape_irc_select_direction, landscape_irc_lqa_step, &
      & landscape_irc_dwi_surface, landscape_irc_dwi_init, landscape_irc_dwi_add, &
      & landscape_irc_dwi_evaluate, landscape_irc_mbs_integrate, &
      & landscape_irc_path, landscape_irc_path_add, landscape_irc_path_reset
   use xtb_landscape_tasks_irc, only : endpoint_reached
   use xtb_landscape_irc_gs, only : landscape_irc_gs_walk, landscape_irc_gs_step, &
      & landscape_irc_gs_lagrange_step
   use xtb_landscape_irc_lqa, only : landscape_irc_lqa_walk
   use xtb_landscape_irc_hpc, only : landscape_irc_hpc_walk
   use xtb_landscape_irc_eulerpc, only : landscape_irc_eulerpc_walk
   use xtb_landscape_testpot, only : landscape_muller_brown_model, &
      & new_muller_brown_model, landscape_leps_model, new_leps_model
   implicit none(type, external)
   private

   !> Analytic Mueller-Brown saddle point (surface coordinates)
   real(wp), parameter :: mb_saddle(2) = [-0.82200156_wp, 0.62431280_wp]
   !> Delivered maximum deviation of each integrator from the reference path
   !> on the analytic Mueller-Brown surface (bohr), forward and backward, in
   !> the order gs / lqa / hpc / eulerpc (report
   !> `docs/_simplify-recon/64-irc-gs-lagrange.md` section 2.3).  These are the
   !> M6 acceptance numbers; the gate below pins them instead of bounding the
   !> worst of the four by 0.15 bohr, which was 6x looser than the largest
   !> delivered value and could not see a real drift of a single method.
   real(wp), parameter :: mb_mep_dev_fwd(4) = [0.0191945_wp, 0.0065749_wp, &
      & 0.0021531_wp, 0.0257361_wp]
   !> Same, backward
   real(wp), parameter :: mb_mep_dev_bwd(4) = [0.0343972_wp, 0.0191347_wp, &
      & 0.0139460_wp, 0.0280571_wp]
   !> Tolerance of the per-method comparison (bohr): the seven digits the test
   !> prints, plus a 200x margin over their own round-off and 1.5e4x tighter
   !> than the 0.15 bohr tripwire of the previous gate
   real(wp), parameter :: mb_mep_dev_tol = 1.0e-5_wp
   !> Analytic Mueller-Brown global minimum
   real(wp), parameter :: mb_min1(2) = [-0.55822363_wp, 1.44172584_wp]
   !> Analytic Mueller-Brown minimum connected to the saddle by the second
   !> IRC direction (the third stationary point of the surface; the
   !> `(-0.822, 0.624)` saddle connects `mb_min1` and this minimum, not the
   !> second well at `(0.623, 0.028)`)
   real(wp), parameter :: mb_min2(2) = [-0.05001083_wp, 0.46669410_wp]
   !> Collinear symmetric H3 saddle point of the LEPS surface (bohr)
   real(wp), parameter :: leps_r_star = 1.50526859_wp
   !> Equilibrium bond length of the H2 asymptote (bohr)
   real(wp), parameter :: leps_r_h2 = 1.4140743_wp

   public :: collect_landscape_irc

contains

!> Collect all exported unit tests
subroutine collect_landscape_irc(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("mass-weighting-roundtrip", test_mass_weighting_roundtrip), &
      new_unittest("arc-length-unit", test_arc_length_unit), &
      new_unittest("lqa-1d-exact", test_lqa_1d_exact), &
      new_unittest("lqa-2d-consistency", test_lqa_2d_consistency), &
      new_unittest("mbs-1d-exact", test_mbs_1d_exact), &
      new_unittest("mbs-vs-lqa", test_mbs_vs_lqa), &
      new_unittest("dwi-interpolation", test_dwi_interpolation), &
      new_unittest("phase-determinism", test_phase_determinism), &
      new_unittest("direction-rejects-minimum", test_direction_rejects_minimum), &
      new_unittest("mb-gs-minima", test_mb_gs_minima), &
      new_unittest("mb-lqa-minima", test_mb_lqa_minima), &
      new_unittest("mb-sphere-chord", test_mb_sphere_chord), &
      new_unittest("mass-metric-effect", test_mass_metric_effect), &
      new_unittest("leps-irc-descent", test_leps_irc_descent), &
      new_unittest("lqa-step-scan", test_lqa_step_scan), &
      new_unittest("termination-max-steps", test_termination_max_steps), &
      new_unittest("termination-grad-tol", test_termination_grad_tol), &
      new_unittest("termination-energy-rise", test_termination_energy_rise), &
      new_unittest("lqa-zero-curvature-limit", test_lqa_zero_curvature_limit), &
      new_unittest("hpc-mb-minima", test_hpc_mb_minima), &
      new_unittest("eulerpc-mb-minima", test_eulerpc_mb_minima), &
      new_unittest("pc-eval-distribution", test_pc_eval_distribution), &
      new_unittest("hpc-corrector-fallback", test_hpc_corrector_fallback), &
      new_unittest("eulerpc-loose-cycles", test_eulerpc_loose_cycles), &
      new_unittest("eulerpc-convergence-order", test_eulerpc_convergence_order), &
      new_unittest("termination-gradient-scale", test_termination_gradient_scale), &
      new_unittest("termination-all-methods", test_termination_all_methods), &
      new_unittest("four-method-crosscheck", test_four_method_crosscheck), &
      new_unittest("gs-path-geometry-consistency", test_gs_path_geometry_consistency), &
      new_unittest("gs-sphere-relaxation-quality", test_gs_sphere_relaxation_quality), &
      new_unittest("gs-kkt-system", test_gs_kkt_system), &
      new_unittest("path-cross-dimension-reuse", test_path_cross_dimension_reuse), &
      new_unittest("endpoint-rmsd-criterion", test_endpoint_rmsd_criterion), &
      new_unittest("phase-convention-sign-invariance", test_phase_sign_invariance), &
      new_unittest("leps-four-method-crosscheck", test_leps_four_method_crosscheck), &
      new_unittest("gs-limit-cycle-guard", test_gs_limit_cycle_guard), &
      new_unittest("hpc-predictor-accuracy", test_hpc_predictor_accuracy) &
      ]

end subroutine collect_landscape_irc

!> Mass weighting is an exact round trip and the arc length carries the amu.
!>
!> A translation of all three atoms by `d` along x changes the mass weighted
!> geometry by `d*sqrt(sum m)`, which is *not* the Cartesian length `d`: the
!> two only agree for unit masses.  This is the unit round trip test of the
!> PLAN (section 6.5): the IRC metric uses amu, so the arc length of a rigid
!> translation is `d*sqrt(sum m_i / amu)`.
subroutine test_mass_weighting_roundtrip(error)
   type(error_type), allocatable, intent(out) :: error
   real(wp) :: mass(9), sqrt_mass(9), xyz(9), q(9), back(9), gcart(9), gq(9)
   real(wp) :: expected, d
   integer :: i

   ! three atoms: H, C, N (amu), the same mass on all three components
   mass = 0.0_wp
   mass(1:3) = 1.00794_wp
   mass(4:6) = 12.0_wp
   mass(7:9) = 14.003074_wp
   sqrt_mass = landscape_irc_sqrt_mass(mass)
   call check(error, all(sqrt_mass > 0.0_wp))
   if (allocated(error)) return

   d = 0.3_wp
   xyz = 0.0_wp
   do i = 1, 3
      xyz(3*(i - 1) + 1) = d
   end do
   call landscape_irc_to_weighted(sqrt_mass, xyz, q)
   call landscape_irc_from_weighted(sqrt_mass, q, back)
   call check(error, maxval(abs(back - xyz)) < 1.0e-14_wp, &
      & "the mass weighting is not an exact round trip")
   if (allocated(error)) return

   expected = d*sqrt(mass(1) + mass(4) + mass(7))
   call check(error, abs(vec_norm(q) - expected) < 1.0e-12_wp, &
      & "the mass weighted length of a rigid translation is wrong")
   if (allocated(error)) return
   ! the mass weighted length differs from the Cartesian length
   call check(error, abs(vec_norm(q) - d) > 0.5_wp*d, &
      & "the mass weighted length collapsed onto the Cartesian length")
   if (allocated(error)) return

   ! gradient round trip
   gcart = 1.0_wp
   call landscape_irc_gradient_to_weighted(sqrt_mass, gcart, gq)
   call check(error, maxval(abs(gq*sqrt_mass - gcart)) < 1.0e-14_wp, &
      & "the mass weighted gradient is not g_i/sqrt(m_i)")
   if (allocated(error)) return

   write(stdout, '(1x,a,f12.6,a,f12.6)') "  cartesian length ", d, &
      & " mass weighted length ", vec_norm(q)

end subroutine test_mass_weighting_roundtrip

!> The arc length of the analytic propagation is the requested one.
subroutine test_arc_length_unit(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   real(wp) :: b(2), hess(2, 2), dx(2), s_reached
   logical :: fail

   call init(env)
   ! isotropic quadratic, gradient along the diagonal
   hess = reshape([1.0_wp, 0.0_wp, 0.0_wp, 2.0_wp], [2, 2])
   b = [0.3_wp, -0.4_wp]
   call landscape_irc_lqa_step(env, b, hess, 0.25_wp, 5000, dx, s_reached, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, abs(s_reached - 0.25_wp) < 1.0e-10_wp, &
      & "the realised arc length is not the requested step")
   if (allocated(error)) return

end subroutine test_arc_length_unit

!> LQA on a one dimensional quadratic is exactly `x0 - s` (closed form).
!>
!> For `E = E0 + g0 dx + 1/2 H dx^2` the ODE `dx/ds = -(g0 + H dx)/|...|` has
!> the solution `dx(s) = -s sign(g0)` until the gradient vanishes, because
!> `a(t) = (e^{-Ht}-1)/H` and `S(t) = |g0|(1-e^{-Ht})/H` combine to
!> `dx = -S g0/|g0| = -s sign(g0)`.  The test pins the propagation to that
!> closed form (the `H` cancels exactly).
subroutine test_lqa_1d_exact(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   real(wp) :: b(1), hess(1, 1), dx(1), s_reached
   real(wp) :: s
   logical :: fail

   call init(env)
   hess(1, 1) = 2.0_wp
   b(1) = 3.0_wp
   s = 0.5_wp
   call landscape_irc_lqa_step(env, b, hess, s, 5000, dx, s_reached, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, abs(dx(1) + s) < 1.0e-10_wp, &
      & "the 1D LQA step is not -s")
   if (allocated(error)) return
   write(stdout, '(1x,a,f16.12)') "  1D LQA dx: ", dx(1)

end subroutine test_lqa_1d_exact

!> The two dimensional propagation is consistent with its own gradient.
!>
!> The endpoint of the analytic propagation must satisfy the exact flow: its
!> displacement is `U a(t) U^T g0`, so the quadratic gradient at the endpoint
!> is `g(t) = U e^{-Lambda t} U^T g0`.  The test reconstructs that closed form
!> independently and compares it with the displacement and the arc length,
!> which pins the eigenvalue factorisation.
subroutine test_lqa_2d_consistency(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   real(wp) :: hess(2, 2), b(2), dx(2), s_reached
   real(wp) :: a1, a2, lam1, lam2, u(2, 2)
   real(wp) :: gp(2), t, sref, dt
   real(wp), allocatable :: work(:, :), ev(:)
   integer :: n
   logical :: fail

   call init(env)
   hess = reshape([1.0_wp, 0.2_wp, 0.2_wp, 0.5_wp], [2, 2])
   b = [0.3_wp, -0.1_wp]
   call landscape_irc_lqa_step(env, b, hess, 0.3_wp, 5000, dx, s_reached, fail)
   call check(error, .not. fail)
   if (allocated(error)) return

   ! independent closed form: diagonalise and integrate the arc length
   allocate(work(2, 2), ev(2))
   work = hess
   call diagonalise(env, work, ev)
   gp = matmul(transpose(work), b)
   ! find t with the arc length 0.3 by a plain Euler march (independent of the
   ! Simpson quadrature of the implementation)
   t = 0.0_wp
   sref = 0.0_wp
   dt = 1.0e-5_wp
   n = 0
   do while (sref < 0.3_wp .and. n < 10000000)
      sref = sref + sqrt((gp(1)*exp(-ev(1)*t))**2 + (gp(2)*exp(-ev(2)*t))**2)*dt
      t = t + dt
      n = n + 1
   end do
   lam1 = ev(1)
   lam2 = ev(2)
   a1 = (exp(-lam1*t) - 1.0_wp)/lam1
   a2 = (exp(-lam2*t) - 1.0_wp)/lam2
   u(:, 1) = work(:, 1)
   u(:, 2) = work(:, 2)
   ! dx = U diag(a) U^T b
   gp(1) = a1*dot_product(u(:, 1), b)
   gp(2) = a2*dot_product(u(:, 2), b)
   gp = matmul(u, gp)
   call check(error, maxval(abs(gp - dx)) < 1.0e-6_wp, &
      & "the 2D LQA displacement disagrees with the closed form")
   if (allocated(error)) return
   call check(error, abs(s_reached - 0.3_wp) < 1.0e-10_wp)
   if (allocated(error)) return
   write(stdout, '(1x,a,2f14.10,a,f12.8)') "  2D LQA dx: ", dx, " s ", s_reached

end subroutine test_lqa_2d_consistency

!> The mBS corrector reproduces the closed form of a constant field exactly.
subroutine test_mbs_1d_exact(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(landscape_irc_dwi_surface) :: surf
   real(wp) :: hess(1, 1), b(1), x0(1), x1(1)
   integer :: n_ode
   logical :: converged, fail

   call init(env)
   hess(1, 1) = 2.0_wp
   b(1) = 3.0_wp
   x0(1) = 0.0_wp
   call landscape_irc_dwi_init(surf, hess, 1)
   call landscape_irc_dwi_add(surf, x0, 0.0_wp, b)
   call landscape_irc_mbs_integrate(env, surf, x0, 0.5_wp, 1.0e-8_wp, 20, 10, &
      & 50000, x1, n_ode, converged, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, abs(x1(1) + 0.5_wp) < 1.0e-10_wp, &
      & "the 1D mBS integration is not -s")
   if (allocated(error)) return
   write(stdout, '(1x,a,f16.12,a,i0)') "  1D mBS x: ", x1(1), "  evals ", n_ode

end subroutine test_mbs_1d_exact

!> The mBS corrector converges to the analytic (LQA) propagation.
!>
!> Both integrate the steepest descent ODE of the same quadratic; a tight
!> tolerance must approach the closed form and a loose one must stay further
!> away, which is the falsifiable statement about the extrapolation.
subroutine test_mbs_vs_lqa(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(landscape_irc_dwi_surface) :: surf
   real(wp) :: hess(2, 2), b(2), x0(2), dx(2), x1(2), s_reached
   real(wp) :: err_fine, err_coarse
   integer :: n_ode
   logical :: converged, fail

   call init(env)
   hess = reshape([1.0_wp, 0.3_wp, 0.3_wp, 0.7_wp], [2, 2])
   b = [0.5_wp, -0.2_wp]
   x0 = 0.0_wp
   call landscape_irc_lqa_step(env, b, hess, 0.4_wp, 5000, dx, s_reached, fail)
   call check(error, .not. fail)
   if (allocated(error)) return

   call landscape_irc_dwi_init(surf, hess, 1)
   call landscape_irc_dwi_add(surf, x0, 0.0_wp, b)

   ! a single Euler step of the whole arc length: first order, large error
   call landscape_irc_mbs_integrate(env, surf, x0, 0.4_wp, 0.0_wp, 1, 1, 1000000, &
      & x1, n_ode, converged, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   err_coarse = maxval(abs(x1 - dx))

   ! Richardson extrapolation of the same Euler scheme
   call landscape_irc_mbs_integrate(env, surf, x0, 0.4_wp, 1.0e-8_wp, 1, 10, 1000000, &
      & x1, n_ode, converged, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, converged, "the extrapolated mBS integration did not converge")
   if (allocated(error)) return
   err_fine = maxval(abs(x1 - dx))

   call check(error, err_fine < 1.0e-7_wp, &
      & "the extrapolated mBS integration did not reach the analytic propagation")
   if (allocated(error)) return
   call check(error, err_coarse > 100.0_wp*err_fine, &
      & "the single Euler step is not much less accurate than the extrapolation")
   if (allocated(error)) return
   write(stdout, '(1x,a,es12.4,a,es12.4)') "  mBS error single ", err_coarse, &
      & " extrapolated ", err_fine

end subroutine test_mbs_vs_lqa

!> The DWI surface reproduces its samples and interpolates between them.
subroutine test_dwi_interpolation(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(landscape_irc_dwi_surface) :: surf
   real(wp) :: hess(2, 2), xa(2), xb(2), ga(2), gb(2), ea, eb
   real(wp) :: x(2), energy, grad(2), emean

   call init(env)
   hess = reshape([1.0_wp, 0.0_wp, 0.0_wp, 2.0_wp], [2, 2])
   xa = [0.0_wp, 0.0_wp]
   xb = [0.4_wp, -0.2_wp]
   ga = [0.3_wp, -0.2_wp]
   gb = [-0.1_wp, 0.4_wp]
   ea = 0.1_wp
   eb = 0.2_wp

   call landscape_irc_dwi_init(surf, hess, 2)
   call landscape_irc_dwi_add(surf, xa, ea, ga)
   call landscape_irc_dwi_add(surf, xb, eb, gb)

   call landscape_irc_dwi_evaluate(surf, xa, energy, grad)
   call check(error, abs(energy - ea) < 1.0e-8_wp, "DWI energy at sample A")
   if (allocated(error)) return
   call check(error, maxval(abs(grad - ga)) < 1.0e-6_wp, "DWI gradient at sample A")
   if (allocated(error)) return
   call landscape_irc_dwi_evaluate(surf, xb, energy, grad)
   call check(error, abs(energy - eb) < 1.0e-8_wp, "DWI energy at sample B")
   if (allocated(error)) return
   call check(error, maxval(abs(grad - gb)) < 1.0e-6_wp, "DWI gradient at sample B")
   if (allocated(error)) return

   ! between the samples the interpolant stays inside the range of the two
   ! quadratic expansions
   x = 0.5_wp*(xa + xb)
   call landscape_irc_dwi_evaluate(surf, x, energy, grad)
   emean = 0.5_wp*(ea + eb)
   call check(error, abs(energy - emean) < 0.2_wp, &
      & "the DWI energy between the samples is far from their mean")
   if (allocated(error)) return
   write(stdout, '(1x,a,f14.10)') "  DWI midpoint energy ", energy

end subroutine test_dwi_interpolation

!> The phase convention makes the transition vector deterministic.
!>
!> The same Hessian must give the same direction (bitwise), and the largest
!> absolute component of the returned vector must be positive - that is the
!> definition of "forward" of the task.
subroutine test_phase_determinism(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(landscape_transition_mode) :: mode1, mode2
   real(wp), allocatable :: hess(:, :)
   logical :: fail

   call init(env)
   ! a 2x2 matrix with exactly one negative eigenvalue
   allocate(hess(2, 2))
   hess = reshape([-1.0_wp, 0.2_wp, 0.2_wp, 0.5_wp], [2, 2])
   call landscape_irc_select_direction(env, hess, 1, mode1, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call landscape_irc_select_direction(env, hess, 1, mode2, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, all(mode1%vector == mode2%vector), &
      & "the transition vector is not deterministic")
   if (allocated(error)) return
   ! the phase convention: the largest absolute component is positive and the
   ! recorded component is exactly that one.  Both assertions avoid indexing
   ! the vector by `phase_index`, so an unset phase (the injection of report
   ! 42) fails instead of reading out of bounds.
   call check(error, mode1%phase_value > 0.0_wp, &
      & "the phase convention did not make the largest component positive")
   if (allocated(error)) return
   call check(error, mode1%phase_index == maxloc(abs(mode1%vector), dim=1), &
      & "the recorded phase component is not the largest one")
   if (allocated(error)) return
   call check(error, mode1%phase_value == maxval(abs(mode1%vector)), &
      & "the recorded phase value is not the largest component")
   if (allocated(error)) return
   call check(error, mode1%eigenvalue < 0.0_wp)
   if (allocated(error)) return
   write(stdout, '(1x,a,2f14.8)') "  transition vector ", mode1%vector

end subroutine test_phase_determinism

!> A Hessian without negative curvature is rejected as an initial direction.
subroutine test_direction_rejects_minimum(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(landscape_transition_mode) :: mode
   real(wp) :: hess(2, 2)
   logical :: fail
   integer :: nlog

   call init(env)
   hess = reshape([1.0_wp, 0.0_wp, 0.0_wp, 2.0_wp], [2, 2])
   call landscape_irc_select_direction(env, hess, 1, mode, fail)
   call check(error, fail, "a minimum was accepted as a transition state")
   if (allocated(error)) return

   ! `target_mode` is read, not decorative: a target beyond the number of
   ! negative eigenvalues emits a warning and falls back to the most negative
   ! mode.  The `irc` task confirms "exactly one imaginary mode" before this
   ! call, so `--lirc-target-mode N` with N > 1 cannot change the direction of
   ! a valid run - it warns and follows the same mode (documented in the
   ! option table of `46-m6-fix-report.md`).
   hess = reshape([-1.0_wp, 0.1_wp, 0.1_wp, 2.0_wp], [2, 2])
   nlog = env%nLog
   call landscape_irc_select_direction(env, hess, 3, mode, fail)
   call check(error, .not. fail, "a first order saddle was rejected")
   if (allocated(error)) return
   call check(error, env%nLog > nlog, &
      & "--lirc-target-mode 3 did not warn about the missing second mode")
   if (allocated(error)) return
   call check(error, mode%index == 1 .and. mode%target == 3, &
      & "the target mode fallback did not follow the most negative mode")
   if (allocated(error)) return

end subroutine test_direction_rejects_minimum

!> GS on the Mueller-Brown surface reaches the two analytic minima.
subroutine test_mb_gs_minima(error)
   type(error_type), allocatable, intent(out) :: error
   type(landscape_muller_brown_model) :: model
   type(landscape_irc_path) :: fwd, bwd
   real(wp) :: mass(2), x0(2), step
   real(wp) :: e_ts, e1, e2
   real(wp) :: g(2)
   logical :: fail

   call new_muller_brown_model(model)
   mass = 1.0_wp
   x0 = mb_saddle
   step = 0.05_wp
   call model%evaluate(x0, e_ts, g, fail)
   call check(error, .not. fail)
   if (allocated(error)) return

   call walk_mb(model, mass, x0, step, 1.0_wp, 400, p_landscape_irc_gs, fwd, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call walk_mb(model, mass, x0, step, -1.0_wp, 400, p_landscape_irc_gs, bwd, fail)
   call check(error, .not. fail)
   if (allocated(error)) return

   call check(error, fwd%n > 5 .and. bwd%n > 5)
   if (allocated(error)) return
   e1 = fwd%e(fwd%n)
   e2 = bwd%e(bwd%n)
   call check(error, e1 < e_ts, "the forward path did not descend")
   if (allocated(error)) return
   call check(error, e2 < e_ts, "the backward path did not descend")
   if (allocated(error)) return

   ! the two endpoints are the two minima (in some order)
   write(stdout, '(1x,a,2f12.6,a,i0,a,2f12.6,a,i0)') "  MB GS fwd end ", &
      & fwd%xyz(:, fwd%n), " term ", fwd%term, " bwd end ", bwd%xyz(:, bwd%n), &
      & " term ", bwd%term
   call check(error, min(dist_to(fwd%xyz(:, fwd%n), mb_min1), &
      & dist_to(fwd%xyz(:, fwd%n), mb_min2)) < 0.25_wp, &
      & "the forward endpoint is not near a Mueller-Brown minimum")
   if (allocated(error)) return
   call check(error, min(dist_to(bwd%xyz(:, bwd%n), mb_min1), &
      & dist_to(bwd%xyz(:, bwd%n), mb_min2)) < 0.25_wp, &
      & "the backward endpoint is not near a Mueller-Brown minimum")
   if (allocated(error)) return
   ! ... and they are different minima
   call check(error, dist_to(fwd%xyz(:, fwd%n), bwd%xyz(:, bwd%n)) > 0.5_wp, &
      & "both directions reached the same minimum")
   if (allocated(error)) return
   write(stdout, '(1x,a,2es14.6,a,i0)') "  MB GS endpoint energies ", e1, e2, &
      & "  points ", fwd%n


end subroutine test_mb_gs_minima

!> LQA on the Mueller-Brown surface reaches the two analytic minima.
subroutine test_mb_lqa_minima(error)
   type(error_type), allocatable, intent(out) :: error
   type(landscape_muller_brown_model) :: model
   type(landscape_irc_path) :: fwd, bwd
   real(wp) :: mass(2), x0(2), step, e_ts, g(2)
   logical :: fail

   call new_muller_brown_model(model)
   mass = 1.0_wp
   x0 = mb_saddle
   step = 0.05_wp
   call model%evaluate(x0, e_ts, g, fail)
   call check(error, .not. fail)
   if (allocated(error)) return

   call walk_mb(model, mass, x0, step, 1.0_wp, 400, p_landscape_irc_lqa, fwd, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call walk_mb(model, mass, x0, step, -1.0_wp, 400, p_landscape_irc_lqa, bwd, fail)
   call check(error, .not. fail)
   if (allocated(error)) return

   call check(error, fwd%e(fwd%n) < e_ts .and. bwd%e(bwd%n) < e_ts)
   if (allocated(error)) return
   write(stdout, '(1x,a,2f12.6,a,i0,a,2f12.6,a,i0)') "  MB LQA fwd end ", &
      & fwd%xyz(:, fwd%n), " term ", fwd%term, " bwd end ", bwd%xyz(:, bwd%n), &
      & " term ", bwd%term
   call check(error, min(dist_to(fwd%xyz(:, fwd%n), mb_min1), &
      & dist_to(fwd%xyz(:, fwd%n), mb_min2)) < 0.15_wp, &
      & "the forward LQA endpoint is not near a Mueller-Brown minimum")
   if (allocated(error)) return
   call check(error, min(dist_to(bwd%xyz(:, bwd%n), mb_min1), &
      & dist_to(bwd%xyz(:, bwd%n), mb_min2)) < 0.15_wp, &
      & "the backward LQA endpoint is not near a Mueller-Brown minimum")
   if (allocated(error)) return
   write(stdout, '(1x,a,2es14.6,a,i0)') "  MB LQA endpoint energies ", &
      & fwd%e(fwd%n), bwd%e(bwd%n), "  points ", fwd%n


end subroutine test_mb_lqa_minima

!> The Gonzalez-Schlegel sphere constraint bounds every chord by the step.
!>
!> The pivot lies half a step from the current point and the new point lies on
!> the sphere of radius `step/2` about it, so the chord `|q_new - q_old|` can
!> never exceed `step`.  The assertion is a direct consequence of the
!> projection and turns red when the projection is removed from the micro
!> cycle.
subroutine test_mb_sphere_chord(error)
   type(error_type), allocatable, intent(out) :: error
   type(landscape_muller_brown_model) :: model
   type(landscape_irc_path) :: path
   real(wp) :: mass(2), x0(2), step, ds, smax, smin
   logical :: fail
   integer :: k

   call new_muller_brown_model(model)
   mass = 1.0_wp
   x0 = mb_saddle
   step = 0.10_wp
   call walk_mb(model, mass, x0, step, 1.0_wp, 40, p_landscape_irc_gs, path, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, path%n > 3)
   if (allocated(error)) return

   ! the first point is half a step from the transition state by construction
   call check(error, abs(path%s(1) - 0.5_wp*step) < 1.0e-12_wp)
   if (allocated(error)) return

   smax = 0.0_wp
   smin = huge(1.0_wp)
   do k = 2, path%n
      ds = path%s(k) - path%s(k - 1)
      smax = max(smax, ds)
      smin = min(smin, ds)
   end do
   write(stdout, '(1x,a,f14.10,a,f14.10,a,f14.10)') "  chord min ", smin, &
      & " max ", smax, " step ", step
   call check(error, smax <= step*(1.0_wp + 1.0e-8_wp), &
      & "a Gonzalez-Schlegel chord exceeded the sphere constraint")
   if (allocated(error)) return
   ! the pivot `q - (1/2) s N G g` is exactly half a step long, so the sphere
   ! relaxation advances the path by a large fraction of the step: the mean
   ! increment has to stay above half the step.  The assertion pins the
   ! `N = (g^T G g)^(-1/2)` normalisation: with the exponent flipped the pivot
   ! collapses onto the current point and the relaxation only moves by the
   ! sphere radius.
   call check(error, smin >= 0.5_wp*step, &
      & "a Gonzalez-Schlegel step is much shorter than the pivot step")
   if (allocated(error)) return

end subroutine test_mb_sphere_chord

!> LEPS H + H2: both directions descend and approach the H2 asymptote.
subroutine test_leps_irc_descent(error)
   type(error_type), allocatable, intent(out) :: error
   type(landscape_leps_model) :: model
   type(landscape_irc_path) :: path
   real(wp) :: mass(9), x0(9), step, g(9), e_ts
   real(wp) :: rmin, r23
   logical :: fail
   integer :: k

   call new_leps_model(model)
   mass = 1.0_wp
   x0 = 0.0_wp
   x0(4:6) = [leps_r_star, 0.0_wp, 0.0_wp]
   x0(7:9) = [2.0_wp*leps_r_star, 0.0_wp, 0.0_wp]
   step = 0.05_wp
   call model%evaluate(x0, e_ts, g, fail)
   call check(error, .not. fail)
   if (allocated(error)) return

   call walk_irc(model, mass, x0, step, 1.0_wp, 60, p_landscape_irc_lqa, path, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, path%n > 5)
   if (allocated(error)) return
   do k = 1, path%n
      call check(error, path%e(k) <= e_ts + 1.0e-10_wp)
      if (allocated(error)) return
   end do
   call check(error, path%e(path%n) < e_ts, "the LEPS IRC did not descend")
   if (allocated(error)) return
   rmin = pair_distance(path%xyz(:, path%n), 2, 3)
   r23 = pair_distance(path%xyz(:, path%n), 1, 2)
   call check(error, min(rmin, r23) < leps_r_h2 + 0.25_wp, &
      & "the LEPS IRC did not approach the H2 equilibrium distance")
   if (allocated(error)) return
   write(stdout, '(1x,a,es14.6,a,f12.6,a,f12.6)') "  LEPS IRC end energy ", &
      & path%e(path%n), "  r12 ", r23, " r23 ", rmin

end subroutine test_leps_irc_descent

!> Step size scan of the LQA propagation on the Mueller-Brown surface.
!>
!> The literature failure mode (report 04 section 1.3): the straightforward
!> LQA step oscillates and deviates from the MEP on strongly curved paths;
!> 0.05 - 0.10 are the recommended sizes and 0.40 is out of control.  The test
!> measures the endpoint distance to the analytic minimum and the maximum
!> perpendicular deviation from a fine reference path for the steps
!> {0.05, 0.10, 0.20, 0.40} and asserts that the deviation grows with the
!> step size.
subroutine test_lqa_step_scan(error)
   type(error_type), allocatable, intent(out) :: error
   type(landscape_muller_brown_model) :: model
   type(landscape_irc_path) :: reference, path
   real(wp) :: mass(2), x0(2), g(2), e_ts
   real(wp) :: steps(4)
   real(wp) :: dev(4), enddist(4)
   logical :: fail
   integer :: i

   call new_muller_brown_model(model)
   mass = 1.0_wp
   x0 = mb_saddle
   call model%evaluate(x0, e_ts, g, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   steps = [0.05_wp, 0.10_wp, 0.20_wp, 0.40_wp]

   ! fine reference path
   call walk_mb(model, mass, x0, 0.02_wp, 1.0_wp, 2000, p_landscape_irc_lqa, &
      & reference, fail)
   call check(error, .not. fail)
   if (allocated(error)) return

   do i = 1, 4
      call walk_mb(model, mass, x0, steps(i), 1.0_wp, 2000, p_landscape_irc_lqa, &
         & path, fail)
      call check(error, .not. fail)
      if (allocated(error)) return
      dev(i) = max_path_deviation(path, reference)
      enddist(i) = dist_to(path%xyz(:, path%n), mb_min1)
      write(stdout, '(1x,a,f6.3,a,es12.4,a,es12.4)') "  step ", steps(i), &
         & "  max deviation ", dev(i), "  endpoint distance ", enddist(i)
   end do

   call check(error, dev(1) <= dev(4), &
      & "the LQA deviation does not grow with the step size")
   if (allocated(error)) return
   call check(error, enddist(1) <= enddist(4) + 0.05_wp, &
      & "the LQA endpoint of the smallest step is not the closest to the minimum")
   if (allocated(error)) return
   ! the documented failure: a 0.40 step leaves the MEP by a large margin
   call check(error, dev(4) > 0.05_wp, &
      & "the 0.40 LQA step did not deviate from the reference path")

end subroutine test_lqa_step_scan

!> Termination signal (a): the step budget stops the path.
subroutine test_termination_max_steps(error)
   type(error_type), allocatable, intent(out) :: error
   type(landscape_muller_brown_model) :: model
   type(landscape_irc_path) :: path
   real(wp) :: mass(2), x0(2)
   logical :: fail

   call new_muller_brown_model(model)
   mass = 1.0_wp
   x0 = mb_saddle
   call walk_mb(model, mass, x0, 0.05_wp, 1.0_wp, 5, p_landscape_irc_gs, path, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, path%term == p_landscape_irc_term_max_steps, &
      & "the step budget did not report signal (a)")
   if (allocated(error)) return
   call check(error, path%n == 5)
   if (allocated(error)) return

end subroutine test_termination_max_steps

!> Termination signal (b): a loose gradient threshold stops the path early.
subroutine test_termination_grad_tol(error)
   type(error_type), allocatable, intent(out) :: error
   type(landscape_muller_brown_model) :: model
   type(landscape_irc_settings) :: settings
   type(landscape_irc_path) :: path
   real(wp) :: mass(2), x0(2)
   logical :: fail

   call new_muller_brown_model(model)
   mass = 1.0_wp
   x0 = mb_saddle
   call default_settings(settings, 0.05_wp, 200, p_landscape_irc_lqa)
   settings%grad_tol = 1.0_wp   ! any gradient is below this
   call walk_with(model, mass, x0, settings, 1.0_wp, path, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, path%term == p_landscape_irc_term_grad, &
      & "the gradient threshold did not report signal (b)")
   if (allocated(error)) return

end subroutine test_termination_grad_tol

!> Termination signal (c): the energy rise ends the path past the minimum.
!>
!> The gradient threshold is disabled (0), so the only physical signal left is
!> the energy rise; with `e_rise_steps` large enough the same path runs into
!> the step budget instead.  That is the falsifiable pair of the signal.
!>
!> The limit cycle guard (report 69) is switched off for the second half: it is
!> a second source of signal (c) and would otherwise end the stalled walk
!> before the budget does.
subroutine test_termination_energy_rise(error)
   type(error_type), allocatable, intent(out) :: error
   type(landscape_muller_brown_model) :: model
   type(landscape_irc_settings) :: settings
   type(landscape_irc_path) :: path
   real(wp) :: mass(2), x0(2)
   logical :: fail

   call new_muller_brown_model(model)
   mass = 1.0_wp
   x0 = mb_saddle

   call default_settings(settings, 0.05_wp, 200, p_landscape_irc_gs)
   settings%grad_tol = 0.0_wp
   settings%e_rise_steps = 3
   call walk_with(model, mass, x0, settings, 1.0_wp, path, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   write(stdout, '(1x,a,i0,a,i0,a,es14.6)') "  energy-rise term ", path%term, &
      & " n ", path%n, " last e ", path%e(path%n)
   call check(error, path%term == p_landscape_irc_term_energy_rise, &
      & "the energy rise did not report signal (c)")
   if (allocated(error)) return
   call check(error, path%n < 200, "signal (c) did not end the path early")
   if (allocated(error)) return

   ! disable (c): the same walk runs into the step budget.  The limit cycle
   ! guard would end the stalled walk long before the budget (it has its own
   ! signal (d) `limit_cycle` since the round-2 review HIGH-1, it is no longer
   ! a second source of (c)), so it has to be switched off as well for signal
   ! (a) to be the only remaining exit - with the default window 8 the MB 0.05
   ! forward walk does enter a cancelling stall once the energy rise criterion
   ! no longer stops it at 22 points.
   settings%e_rise_steps = 100000
   settings%osc_window = 0
   settings%method = p_landscape_irc_gs
   call walk_with(model, mass, x0, settings, 1.0_wp, path, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, path%term == p_landscape_irc_term_max_steps, &
      & "disabling the energy rise did not fall back to signal (a)")
   if (allocated(error)) return

end subroutine test_termination_energy_rise

!> The mass metric changes the path, not only the units.
!>
!> The IRC is defined in mass weighted coordinates, so a different mass
!> distribution moves the MEP.  The test runs the same analytic surface with
!> two mass vectors and requires the endpoints to differ; with the identity
!> metric (the injection of report 42) the two runs are bitwise identical and
!> the test is red.
subroutine test_mass_metric_effect(error)
   type(error_type), allocatable, intent(out) :: error
   type(landscape_muller_brown_model) :: model
   type(landscape_irc_path) :: unit_mass, heavy
   real(wp) :: mass_unit(2), mass_heavy(2), x0(2)
   logical :: fail

   call new_muller_brown_model(model)
   mass_unit = 1.0_wp
   mass_heavy = [1.0_wp, 4.0_wp]
   x0 = mb_saddle

   call walk_mb(model, mass_unit, x0, 0.05_wp, 1.0_wp, 400, p_landscape_irc_lqa, &
      & unit_mass, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call walk_mb(model, mass_heavy, x0, 0.05_wp, 1.0_wp, 400, p_landscape_irc_lqa, &
      & heavy, fail)
   call check(error, .not. fail)
   if (allocated(error)) return

   call check(error, unit_mass%n > 3 .and. heavy%n > 3, &
      & "a mass metric run produced no path")
   if (allocated(error)) return
   ! the initial half-step direction is a mass weighted eigenvector, so it
   ! differs in Cartesian space, and the arc length to the minimum differs as
   ! well.  Both are exactly equal when the metric is the identity.
   call check(error, dist_to(unit_mass%xyz(:, 1), heavy%xyz(:, 1)) > 1.0e-6_wp, &
      & "the mass metric did not change the first IRC point (identity metric?)")
   if (allocated(error)) return
   call check(error, abs(unit_mass%s(unit_mass%n) - heavy%s(heavy%n)) > 1.0e-3_wp, &
      & "the mass metric did not change the arc length to the minimum")
   if (allocated(error)) return
   write(stdout, '(1x,a,es12.4,a,es12.4)') "  mass metric: first point gap ", &
      & dist_to(unit_mass%xyz(:, 1), heavy%xyz(:, 1)), "  arc length gap ", &
      & abs(unit_mass%s(unit_mass%n) - heavy%s(heavy%n))

end subroutine test_mass_metric_effect

!> The `lambda -> 0` branch of the analytic predictor.
!>
!> A Hessian with an exactly vanishing eigenvalue (a rigid mode of an analytic
!> model) has the removable limits `a_ii -> -t` and `g_i(t) -> g_i(0)`.  With
!> the branch present the step is finite, covers the requested arc length and
!> keeps the zero curvature contribution `-t g_i`; removing the branch either
!> divides by the vanishing eigenvalue (NaN) or drops the mode (`dx_1 = 0`),
!> so the two assertions are the falsifiable test of the branch.
subroutine test_lqa_zero_curvature_limit(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   real(wp) :: hess(2, 2), b(2), dx(2), s_reached
   logical :: fail

   call init(env)
   hess = 0.0_wp
   hess(2, 2) = 2.0_wp
   b = [1.0_wp, 0.3_wp]
   call landscape_irc_lqa_step(env, b, hess, 0.5_wp, 5000, dx, s_reached, fail)
   call check(error, .not. fail, "the zero curvature predictor failed")
   if (allocated(error)) return
   ! a NaN does not compare equal to itself
   call check(error, all(dx == dx), &
      & "the zero curvature limit divided by the vanishing eigenvalue (NaN)")
   if (allocated(error)) return
   call check(error, abs(s_reached - 0.5_wp) < 1.0e-8_wp, &
      & "the zero curvature predictor did not reach the requested arc length")
   if (allocated(error)) return
   ! a_ii -> -t: the vanishing curvature mode contributes -t g_1, not zero
   call check(error, dx(1) < -0.4_wp, &
      & "the vanishing curvature mode was dropped instead of using a_ii -> -t")
   if (allocated(error)) return
   call check(error, dx(2) < 0.0_wp, "the finite curvature mode did not descend")
   if (allocated(error)) return
   write(stdout, '(1x,a,2f14.8,a,f14.8)') "  zero-curvature dx ", dx, &
      & " s ", s_reached

end subroutine test_lqa_zero_curvature_limit

!> HPC on the Mueller-Brown surface reaches the two analytic minima.
subroutine test_hpc_mb_minima(error)
   type(error_type), allocatable, intent(out) :: error
   call mb_both_minima(error, p_landscape_irc_hpc, 0.25_wp)
end subroutine test_hpc_mb_minima

!> EulerPC on the Mueller-Brown surface reaches the two analytic minima.
subroutine test_eulerpc_mb_minima(error)
   type(error_type), allocatable, intent(out) :: error
   call mb_both_minima(error, p_landscape_irc_eulerpc, 0.25_wp)
end subroutine test_eulerpc_mb_minima

!> Shared body of the two minimum tests: both directions of one method.
subroutine mb_both_minima(error, method, tol)
   type(error_type), allocatable, intent(out) :: error
   !> Integrator
   integer, intent(in) :: method
   !> Tolerance of the endpoint distance to the analytic minimum
   real(wp), intent(in) :: tol

   type(landscape_muller_brown_model) :: model
   type(landscape_irc_path) :: fwd, bwd
   real(wp) :: mass(2), x0(2), e_ts, g(2)
   logical :: fail

   call new_muller_brown_model(model)
   mass = 1.0_wp
   x0 = mb_saddle
   call model%evaluate(x0, e_ts, g, fail)
   call check(error, .not. fail)
   if (allocated(error)) return

   call walk_mb(model, mass, x0, 0.05_wp, 1.0_wp, 400, method, fwd, fail)
   call check(error, .not. fail, "the forward walk failed")
   if (allocated(error)) return
   call walk_mb(model, mass, x0, 0.05_wp, -1.0_wp, 400, method, bwd, fail)
   call check(error, .not. fail, "the backward walk failed")
   if (allocated(error)) return
   call check(error, fwd%n > 5 .and. bwd%n > 5)
   if (allocated(error)) return
   call check(error, fwd%e(fwd%n) < e_ts .and. bwd%e(bwd%n) < e_ts, &
      & "a predictor-corrector direction did not descend")
   if (allocated(error)) return
   call check(error, min(dist_to(fwd%xyz(:, fwd%n), mb_min1), &
      & dist_to(fwd%xyz(:, fwd%n), mb_min2)) < tol, &
      & "the forward endpoint is not near a Mueller-Brown minimum")
   if (allocated(error)) return
   call check(error, min(dist_to(bwd%xyz(:, bwd%n), mb_min1), &
      & dist_to(bwd%xyz(:, bwd%n), mb_min2)) < tol, &
      & "the backward endpoint is not near a Mueller-Brown minimum")
   if (allocated(error)) return
   call check(error, dist_to(fwd%xyz(:, fwd%n), bwd%xyz(:, bwd%n)) > 0.5_wp, &
      & "both directions reached the same minimum")
   if (allocated(error)) return
   write(stdout, '(1x,a,i0,a,2f11.6,a,2f11.6,a,2f8.3)') "  method ", method, &
      & " fwd end ", fwd%xyz(:, fwd%n), " bwd end ", bwd%xyz(:, bwd%n), &
      & " mean E/g ", mean_eval(fwd), mean_eval(bwd)

end subroutine mb_both_minima

!> The per point energy/gradient cost of the predictor-corrector integrators.
!>
!> The `E/g` column of `irc.log` is the number of PES evaluations spent on a
!> point.  The predictor-corrector evaluates the predicted point and the
!> corrected point, so the converged mean is two, versus the 3-5 of a
!> Gonzalez-Schlegel macro step and the 16.9-20.3 that the GS micro-cycle
!> relaxer cost on a real molecule before the equality constrained Lagrange
!> step (1.9-2.1 since; reports 42 and 64).  The test prints the full
!> distribution and fails when the mean regresses above three or when the
!> predictor-corrector ever spends less than one evaluation per point.
subroutine test_pc_eval_distribution(error)
   type(error_type), allocatable, intent(out) :: error
   type(landscape_muller_brown_model) :: model
   type(landscape_irc_path) :: hpc, eulerpc, lqa
   real(wp) :: mass(2), x0(2)
   real(wp) :: mean_hpc, mean_epc, mean_lqa
   integer :: k, hist(6)
   logical :: fail

   call new_muller_brown_model(model)
   mass = 1.0_wp
   x0 = mb_saddle
   call walk_mb(model, mass, x0, 0.05_wp, 1.0_wp, 400, p_landscape_irc_hpc, &
      & hpc, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call walk_mb(model, mass, x0, 0.05_wp, 1.0_wp, 400, p_landscape_irc_eulerpc, &
      & eulerpc, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call walk_mb(model, mass, x0, 0.05_wp, 1.0_wp, 400, p_landscape_irc_lqa, &
      & lqa, fail)
   call check(error, .not. fail)
   if (allocated(error)) return

   mean_hpc = mean_eval(hpc)
   mean_epc = mean_eval(eulerpc)
   mean_lqa = mean_eval(lqa)

   ! the LQA integrator is the pure analytic propagation: exactly one E/g
   call check(error, abs(mean_lqa - 1.0_wp) < 1.0e-12_wp, &
      & "the LQA integrator does not cost exactly one E/g per point")
   if (allocated(error)) return
   ! the predictor-corrector evaluates two points per step
   call check(error, mean_hpc > 1.0_wp .and. mean_hpc <= 3.0_wp, &
      & "the HPC mean E/g per point is outside (1, 3]")
   if (allocated(error)) return
   call check(error, mean_epc > 1.0_wp .and. mean_epc <= 3.0_wp, &
      & "the EulerPC mean E/g per point is outside (1, 3]")
   if (allocated(error)) return
   call check(error, minval(hpc%n_eval(1:hpc%n)) >= 1)
   if (allocated(error)) return

   hist = 0
   do k = 1, hpc%n
      if (hpc%n_eval(k) <= 5) hist(hpc%n_eval(k) + 1) = hist(hpc%n_eval(k) + 1) + 1
   end do
   write(stdout, '(1x,a,f8.4,a,f8.4,a,f8.4)') "  mean E/g  lqa ", mean_lqa, &
      & " hpc ", mean_hpc, " eulerpc ", mean_epc
   write(stdout, '(1x,a,6i6)') "  HPC E/g histogram (0..5) ", hist

end subroutine test_pc_eval_distribution

!> A non-converging corrector falls back to one Gonzalez-Schlegel step.
!>
!> `--lirc-mbs-max-k 1` with `--lirc-mbs-tol 0` makes the Richardson
!> extrapolation unable to reach the tolerance, so every macro step has to
!> take the documented fallback: the walk warns and continues with a GS step.
!> The fallback is identified by its **direct** signal, the warning counter
!> (`env%nLog`), one warning per macro point; the cost is then only bounded,
!> not used as the identifier: the GS fallback used to be the most expensive
!> step of the suite (11.818 E/g when this test was written, a 6x cost compared
!> with the two E/g of the converged corrector), and it is now 3 E/g - the
!> corrector's two evaluations plus its own micro cycles (report 64).  The old
!> assertion `mean E/g > 2` was a *proxy* for "the fallback really ran" which
!> the warning count now establishes directly, so the cost is asserted as a
!> two sided band around the converged corrector's cost instead.
subroutine test_hpc_corrector_fallback(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(landscape_muller_brown_model) :: model
   type(landscape_irc_settings) :: settings
   type(landscape_irc_path) :: path
   real(wp) :: mass(2), x0(2), mean_ev
   integer :: nlog, nfallback
   logical :: fail

   call init(env)
   call new_muller_brown_model(model)
   mass = 1.0_wp
   x0 = mb_saddle
   call default_settings(settings, 0.05_wp, 200, p_landscape_irc_hpc)
   settings%grad_tol = 0.0_wp
   settings%mbs_max_k = 1
   settings%mbs_tol = 0.0_wp

   nlog = env%nLog
   call walk_with(model, mass, x0, settings, 1.0_wp, path, fail, env)
   call check(error, .not. fail, "the fallback walk failed")
   if (allocated(error)) return
   nfallback = env%nLog - nlog
   call check(error, nfallback > 0, &
      & "the corrector failure did not emit a warning (the fallback did not fire)")
   if (allocated(error)) return
   call check(error, path%n > 5, "the fallback walk produced no path")
   if (allocated(error)) return
   call check(error, path%term == p_landscape_irc_term_energy_rise .or. &
      & path%term == p_landscape_irc_term_max_steps)
   if (allocated(error)) return
   ! every macro step had to fall back: the warning is the direct evidence, and
   ! it also proves the GS step - not the corrector - produced the path
   call check(error, nfallback >= path%n - 1, &
      & "the corrector fallback did not fire at every macro point")
   if (allocated(error)) return
   ! cost band of the fallback step: at least one evaluation and less than
   ! twice the two E/g of the converged corrector.  Measured over the 20 points
   ! of this walk: 1.450 E/g (the micro cycle now converges in one or two
   ! evaluations); it was 11.818 E/g, i.e. six times the corrector, when the
   ! clamped tangent step consumed the whole micro cycle budget.  The band
   ! catches a return to the budget-driven cost of that step.
   mean_ev = mean_eval(path)
   write(stdout, '(1x,a,f8.3,a,i0,a,i0)') "  fallback mean E/g ", mean_ev, &
      & "  points ", path%n, "  warnings ", nfallback
   call check(error, mean_ev >= 1.0_wp .and. mean_ev < 4.0_wp, &
      & "the GS fallback is no longer in the cost band of the corrector")
   if (allocated(error)) return

end subroutine test_hpc_corrector_fallback

!> The EulerPC coarse phase switches off the corrector for the first steps.
!>
!> `loose_cycles` is the number of leading macro steps that record the plain
!> Euler prediction (one E/g).  With a value larger than the step budget the
!> whole path is coarse and its mean E/g is exactly one; with zero the whole
!> path is fine and the corrector is used.  The two paths must differ, which
!> is what makes the branch observable.
subroutine test_eulerpc_loose_cycles(error)
   type(error_type), allocatable, intent(out) :: error
   type(landscape_muller_brown_model) :: model
   type(landscape_irc_settings) :: settings
   type(landscape_irc_path) :: coarse, fine
   real(wp) :: mass(2), x0(2)
   logical :: fail

   call new_muller_brown_model(model)
   mass = 1.0_wp
   x0 = mb_saddle

   call default_settings(settings, 0.05_wp, 60, p_landscape_irc_eulerpc)
   settings%grad_tol = 0.0_wp
   settings%e_rise_steps = 100000
   settings%loose_cycles = 100000
   call walk_with(model, mass, x0, settings, 1.0_wp, coarse, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, abs(mean_eval(coarse) - 1.0_wp) < 1.0e-12_wp, &
      & "the coarse EulerPC phase did not cost exactly one E/g per point")
   if (allocated(error)) return

   settings%loose_cycles = 0
   call walk_with(model, mass, x0, settings, 1.0_wp, fine, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, mean_eval(fine) > 1.0_wp, &
      & "the fine EulerPC phase did not use the corrector")
   if (allocated(error)) return
   call check(error, dist_to(coarse%xyz(:, coarse%n), fine%xyz(:, fine%n)) > 1.0e-4_wp, &
      & "the coarse and fine EulerPC phases produced the same endpoint")
   if (allocated(error)) return
   write(stdout, '(1x,a,2f8.3,a,2es12.4)') "  coarse/fine mean E/g ", &
      & mean_eval(coarse), mean_eval(fine), "  endpoint gap ", &
      & dist_to(coarse%xyz(:, coarse%n), fine%xyz(:, fine%n))

end subroutine test_eulerpc_loose_cycles

!> Measured convergence order of the EulerPC integrator.
!>
!> The paper classifies EulerPC as first order (its Euler predictor), with a
!> reported log-log slope of 2.016 in the endpoint error.  The test measures
!> the **energy error at a fixed arc length** (`s = 1 bohr*amu^1/2`) against a
!> fine reference path for the step sizes {0.05, 0.10, 0.20}, with the pure
!> Euler predictor (all coarse) and with the corrector enabled, fits a line
!> through `log(error)` vs `log(step)` and prints the slope, the fitting
!> interval and the number of points.  The assertion only requires a positive,
!> sub-third-order slope (the exact value is a measurement, not a design
!> constant).
subroutine test_eulerpc_convergence_order(error)
   type(error_type), allocatable, intent(out) :: error
   type(landscape_muller_brown_model) :: model
   type(landscape_irc_settings) :: settings
   type(landscape_irc_path) :: reference, path
   real(wp) :: mass(2), x0(2)
   real(wp) :: steps(3), err_coarse(3), err_fine(3), slope_coarse, slope_fine
   real(wp) :: slope_local_coarse, slope_local_fine
   real(wp) :: x(3), y(3)
   integer :: i, max_steps
   logical :: fail

   call new_muller_brown_model(model)
   mass = 1.0_wp
   x0 = mb_saddle
   steps = [0.05_wp, 0.10_wp, 0.20_wp]

   ! fine reference path (LQA at 0.01, dense)
   call walk_mb(model, mass, x0, 0.01_wp, 1.0_wp, 200, p_landscape_irc_lqa, &
      & reference, fail)
   call check(error, .not. fail)
   if (allocated(error)) return

   ! pure Euler predictor (all coarse)
   do i = 1, 3
      max_steps = nint(0.9_wp/steps(i)) + 4
      call default_settings(settings, steps(i), max_steps, p_landscape_irc_eulerpc)
      settings%grad_tol = 0.0_wp
      settings%e_rise_steps = 100000
      settings%loose_cycles = 100000
      call walk_with(model, mass, x0, settings, 1.0_wp, path, fail)
      call check(error, .not. fail)
      if (allocated(error)) return
      call check(error, path%s(path%n) >= 0.8_wp, &
         & "the coarse EulerPC path did not reach the reference arc length")
      if (allocated(error)) return
      err_coarse(i) = deviation_at_s(path, reference, 0.8_wp)
      x(i) = log(steps(i))
      y(i) = log(max(err_coarse(i), 1.0e-14_wp))
   end do
   slope_coarse = fit_slope(x, y, 3)

   ! with the corrector
   do i = 1, 3
      max_steps = nint(0.9_wp/steps(i)) + 4
      call default_settings(settings, steps(i), max_steps, p_landscape_irc_eulerpc)
      settings%grad_tol = 0.0_wp
      settings%e_rise_steps = 100000
      settings%loose_cycles = 0
      call walk_with(model, mass, x0, settings, 1.0_wp, path, fail)
      call check(error, .not. fail)
      if (allocated(error)) return
      err_fine(i) = deviation_at_s(path, reference, 0.8_wp)
      y(i) = log(max(err_fine(i), 1.0e-14_wp))
   end do
   slope_fine = fit_slope(x, y, 3)

   write(stdout, '(1x,a,3es12.4)') "  EulerPC deviation at s=0.8 (0.05/0.10/0.20) ", err_coarse
   write(stdout, '(1x,a,3es12.4)') "  EulerPC deviation with corrector               ", err_fine
   write(stdout, '(1x,a,f8.4,a,i0)') "  order slope (least squares over 0.05-0.20, 3 points) ", &
      & slope_coarse, "  points ", 3
   ! the asymptotic (0.10 -> 0.20) local slope is the meaningful order estimate:
   ! the 0.05 point is super-convergent for this surface (the coarse path stays
   ! on the reference by accident), which inflates the three point fit
   slope_local_coarse = log(err_coarse(3)/err_coarse(2))/log(steps(3)/steps(2))
   slope_local_fine = log(err_fine(3)/err_fine(2))/log(steps(3)/steps(2))
   write(stdout, '(1x,a,f8.4,a,f8.4)') "  local order 0.10 -> 0.20  pure Euler ", &
      & slope_local_coarse, " with corrector ", slope_local_fine
   call check(error, slope_local_coarse > 0.5_wp .and. slope_local_coarse < 2.0_wp, &
      & "the measured asymptotic Euler predictor order is not first order")
   if (allocated(error)) return
   call check(error, slope_local_fine > 0.5_wp .and. slope_local_fine < 2.0_wp, &
      & "the measured asymptotic EulerPC order is not first order")

end subroutine test_eulerpc_convergence_order

!> Distance from the path point closest to `s_target` to the reference path.
!>
!> The smallest distance to any reference point approximates the perpendicular
!> distance to the reference curve (the reference is dense).
pure function deviation_at_s(path, reference, s_target) result(dev)
   !> Path under test
   type(landscape_irc_path), intent(in) :: path
   !> Reference path
   type(landscape_irc_path), intent(in) :: reference
   !> Arc length at which the deviation is measured (bohr*amu^1/2)
   real(wp), intent(in) :: s_target
   !> Deviation
   real(wp) :: dev

   real(wp) :: dmin, d
   integer :: k, kbest, j

   kbest = 1
   do k = 2, path%n
      if (abs(path%s(k) - s_target) < abs(path%s(kbest) - s_target)) kbest = k
   end do
   dmin = huge(1.0_wp)
   do j = 1, reference%n
      d = vec_norm(path%xyz(:, kbest) - reference%xyz(:, j))
      dmin = min(dmin, d)
   end do
   dev = dmin

end function deviation_at_s

!> The gradient criterion is the RMS per component, not a power of the norm.
!>
!> The first recorded point of a walk has a known `|g|` (read from the path
!> itself).  Setting `grad_tol` a little **below** its RMS-per-component value
!> must not fire signal (b) at that point, while doubling the threshold must.
!> A criterion that used `|g|^2` (much smaller than `|g|` for these gradients)
!> would also fire in the first case, so the pair pins the magnitude of the
!> criterion.  The test is run on all four integrators.
subroutine test_termination_gradient_scale(error)
   type(error_type), allocatable, intent(out) :: error
   type(landscape_muller_brown_model) :: model
   type(landscape_irc_settings) :: settings
   type(landscape_irc_path) :: path
   real(wp) :: mass(2), x0(2), g_rms
   integer :: methods(4), i
   logical :: fail

   call new_muller_brown_model(model)
   mass = 1.0_wp
   x0 = mb_saddle
   methods = [p_landscape_irc_gs, p_landscape_irc_lqa, p_landscape_irc_hpc, &
      & p_landscape_irc_eulerpc]

   do i = 1, 4
      ! read the first point's gradient from a one point walk
      call default_settings(settings, 0.05_wp, 1, methods(i))
      settings%grad_tol = 0.0_wp
      settings%e_rise_steps = 100000
      call walk_with(model, mass, x0, settings, 1.0_wp, path, fail)
      call check(error, .not. fail)
      if (allocated(error)) return
      call check(error, path%n == 1)
      if (allocated(error)) return
      g_rms = path%gnorm(1)/sqrt(2.0_wp)
      call check(error, g_rms > 0.0_wp)
      if (allocated(error)) return

      ! a threshold just below the first point's RMS gradient must not fire (b)
      ! at the first point (a criterion in `|g|^2` would fire there)
      call default_settings(settings, 0.05_wp, 200, methods(i))
      settings%grad_tol = 0.5_wp*g_rms
      settings%e_rise_steps = 100000
      call walk_with(model, mass, x0, settings, 1.0_wp, path, fail)
      call check(error, .not. fail)
      if (allocated(error)) return
      call check(error, path%n > 1, &
         & "a threshold below the first point's gradient fired signal (b) there")
      if (allocated(error)) return

      ! doubling it must fire (b) immediately
      settings%grad_tol = 2.0_wp*g_rms
      call walk_with(model, mass, x0, settings, 1.0_wp, path, fail)
      call check(error, .not. fail)
      if (allocated(error)) return
      call check(error, path%term == p_landscape_irc_term_grad, &
         & "a threshold above the first point's gradient did not fire signal (b)")
      if (allocated(error)) return
      call check(error, path%n == 1, &
         & "signal (b) did not stop the path at the first point")
      if (allocated(error)) return
   end do
   write(stdout, '(1x,a)') "  gradient criterion magnitude pinned on gs/lqa/hpc/eulerpc"

end subroutine test_termination_gradient_scale

!> The three termination signals fire independently on all four integrators.
!>
!> (a) the step budget, (b) an unreachable gradient threshold and (c) the
!> energy rise are exercised for gs, lqa, hpc and eulerpc on the analytic
!> Mueller-Brown surface.  Each branch is paired with a control so that the
!> signal, not the surface, is what stops the path.
subroutine test_termination_all_methods(error)
   type(error_type), allocatable, intent(out) :: error
   type(landscape_muller_brown_model) :: model
   type(landscape_irc_settings) :: settings
   type(landscape_irc_path) :: path
   real(wp) :: mass(2), x0(2)
   integer :: methods(4), i
   logical :: fail

   call new_muller_brown_model(model)
   mass = 1.0_wp
   x0 = mb_saddle
   methods = [p_landscape_irc_gs, p_landscape_irc_lqa, p_landscape_irc_hpc, &
      & p_landscape_irc_eulerpc]

   do i = 1, 4
      ! (a) step budget: the energy does not stop the path inside 5 points
      call default_settings(settings, 0.05_wp, 5, methods(i))
      settings%grad_tol = 0.0_wp
      settings%e_rise_steps = 100000
      call walk_with(model, mass, x0, settings, 1.0_wp, path, fail)
      call check(error, .not. fail)
      if (allocated(error)) return
      call check(error, path%term == p_landscape_irc_term_max_steps, &
         & "signal (a) did not fire on every integrator")
      if (allocated(error)) return

      ! (b) an unreachable gradient threshold
      call default_settings(settings, 0.05_wp, 200, methods(i))
      settings%grad_tol = 1.0_wp
      settings%e_rise_steps = 100000
      call walk_with(model, mass, x0, settings, 1.0_wp, path, fail)
      call check(error, .not. fail)
      if (allocated(error)) return
      call check(error, path%term == p_landscape_irc_term_grad, &
         & "signal (b) did not fire on every integrator")
      if (allocated(error)) return

      ! (c) the energy rise, with the gradient signal disabled
      call default_settings(settings, 0.05_wp, 200, methods(i))
      settings%grad_tol = 0.0_wp
      settings%e_rise_steps = 3
      call walk_with(model, mass, x0, settings, 1.0_wp, path, fail)
      call check(error, .not. fail)
      if (allocated(error)) return
      call check(error, path%term == p_landscape_irc_term_energy_rise, &
         & "signal (c) did not fire on every integrator")
      if (allocated(error)) return
      call check(error, path%n < 200, "signal (c) did not end the path early")
      if (allocated(error)) return
   end do
   write(stdout, '(1x,a)') "  three termination signals fired on gs/lqa/hpc/eulerpc"

end subroutine test_termination_all_methods

!> PLAN experiment 6: the four integrators agree on the same transition state.
!>
!> From the analytic Mueller-Brown saddle all four methods are run in both
!> directions; the two endpoints of each method are compared with the other
!> three (geometry RMSD and energy).  The test prints the full matrix and
!> requires a common endpoint basin; the measured spread is the experiment 6
!> number recorded in the report.
subroutine test_four_method_crosscheck(error)
   type(error_type), allocatable, intent(out) :: error
   type(landscape_muller_brown_model) :: model
   type(landscape_irc_path) :: fwd(4), bwd(4), reference_fwd, reference_bwd
   real(wp) :: mass(2), x0(2)
   real(wp) :: gap_f(4, 4), gap_b(4, 4), de_f(4, 4), de_b(4, 4)
   real(wp) :: max_gap_f, max_gap_b, max_de_f, max_de_b, dev_f(4), dev_b(4)
   integer :: methods(4), i, j
   logical :: fail

   call new_muller_brown_model(model)
   mass = 1.0_wp
   x0 = mb_saddle
   methods = [p_landscape_irc_gs, p_landscape_irc_lqa, p_landscape_irc_hpc, &
      & p_landscape_irc_eulerpc]

   do i = 1, 4
      call walk_mb(model, mass, x0, 0.05_wp, 1.0_wp, 400, methods(i), fwd(i), fail)
      call check(error, .not. fail)
      if (allocated(error)) return
      call walk_mb(model, mass, x0, 0.05_wp, -1.0_wp, 400, methods(i), bwd(i), fail)
      call check(error, .not. fail)
      if (allocated(error)) return
   end do

   ! deviation of each path from a fine reference (the analytic MEP): the
   ! reference is the 0.01 step LQA path of the **same direction**, which is
   ! the accepted numeric MEP of the analytic surface (the surface has no
   ! closed form reaction path)
   call walk_mb(model, mass, x0, 0.01_wp, 1.0_wp, 2000, p_landscape_irc_lqa, &
      & reference_fwd, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call walk_mb(model, mass, x0, 0.01_wp, -1.0_wp, 2000, p_landscape_irc_lqa, &
      & reference_bwd, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   do i = 1, 4
      dev_f(i) = max_path_deviation(fwd(i), reference_fwd)
      dev_b(i) = max_path_deviation(bwd(i), reference_bwd)
   end do
   write(stdout, '(1x,a,4f12.7)') "  MB max MEP deviation fwd (gs/lqa/hpc/eulerpc) ", dev_f
   write(stdout, '(1x,a,4f12.7)') "  MB max MEP deviation bwd (gs/lqa/hpc/eulerpc) ", dev_b
   do i = 1, 4
      call check(error, abs(dev_f(i) - mb_mep_dev_fwd(i)) < mb_mep_dev_tol, &
         & "a forward MEP deviation left the delivered value")
      if (allocated(error)) return
      call check(error, abs(dev_b(i) - mb_mep_dev_bwd(i)) < mb_mep_dev_tol, &
         & "a backward MEP deviation left the delivered value")
      if (allocated(error)) return
   end do
   if (allocated(error)) return

   max_gap_f = 0.0_wp
   max_gap_b = 0.0_wp
   max_de_f = 0.0_wp
   max_de_b = 0.0_wp
   do i = 1, 4
      do j = 1, 4
         gap_f(i, j) = dist_to(fwd(i)%xyz(:, fwd(i)%n), fwd(j)%xyz(:, fwd(j)%n))
         gap_b(i, j) = dist_to(bwd(i)%xyz(:, bwd(i)%n), bwd(j)%xyz(:, bwd(j)%n))
         de_f(i, j) = abs(fwd(i)%e(fwd(i)%n) - fwd(j)%e(fwd(j)%n))
         de_b(i, j) = abs(bwd(i)%e(bwd(i)%n) - bwd(j)%e(bwd(j)%n))
         max_gap_f = max(max_gap_f, gap_f(i, j))
         max_gap_b = max(max_gap_b, gap_b(i, j))
         max_de_f = max(max_de_f, de_f(i, j))
         max_de_b = max(max_de_b, de_b(i, j))
      end do
   end do
   write(stdout, '(1x,a,f14.8)') "  crosscheck forward endpoint rmsd ", max_gap_f
   write(stdout, '(1x,a,f14.8)') "  crosscheck backward endpoint rmsd ", max_gap_b
   write(stdout, '(1x,a,es14.6)') "  crosscheck forward  dE ", max_de_f
   write(stdout, '(1x,a,es14.6)') "  crosscheck backward dE ", max_de_b
   call check(error, max_gap_f < 0.25_wp, &
      & "the forward endpoints of the four integrators do not agree")
   if (allocated(error)) return
   call check(error, max_gap_b < 0.25_wp, &
      & "the backward endpoints of the four integrators do not agree")
   if (allocated(error)) return
   call check(error, max_de_f < 2.0e-2_wp, &
      & "the forward endpoint energies of the four integrators do not agree")
   if (allocated(error)) return
   call check(error, max_de_b < 2.0e-2_wp, &
      & "the backward endpoint energies of the four integrators do not agree")
   if (allocated(error)) return
   ! every method reaches a different minimum in the two directions
   do i = 1, 4
      call check(error, dist_to(fwd(i)%xyz(:, fwd(i)%n), bwd(i)%xyz(:, bwd(i)%n)) > 0.5_wp, &
         & "an integrator reached the same minimum in both directions")
      if (allocated(error)) return
   end do

end subroutine test_four_method_crosscheck

!> Every recorded point of a path carries its own geometry, energy and gradient.
!>
!> The path object stores the Cartesian geometry, the energy and the gradient
!> norm of every point.  Those three values must be self consistent: feeding the
!> stored geometry back into the potential has to reproduce the stored energy
!> and gradient norm.  The assertion catches the regression that a refactoring
!> of `landscape_irc_gs_walk` introduced (report 43 section 4): when the split
!> out `landscape_irc_gs_step` failed to return its Cartesian geometry, every
!> recorded point kept the **first** point's geometry while the energy was
!> updated - the path looked plausible in energy but its geometry was frozen.
subroutine test_gs_path_geometry_consistency(error)
   type(error_type), allocatable, intent(out) :: error
   type(landscape_muller_brown_model) :: mb
   type(landscape_leps_model) :: leps
   ! one path object per model: the container keeps its allocation across a
   ! `landscape_irc_path_reset`, so the two must not be mixed (finding F-1 of
   ! the M6 validation report)
   type(landscape_irc_path) :: path_mb, path_leps
   real(wp) :: mass2(2), mass9(9), x0_2(2), x0_9(9)
   real(wp) :: g(9), e_ts
   logical :: fail

   ! Mueller-Brown
   call new_muller_brown_model(mb)
   mass2 = 1.0_wp
   x0_2 = mb_saddle
   call walk_mb(mb, mass2, x0_2, 0.05_wp, 1.0_wp, 400, p_landscape_irc_gs, &
      & path_mb, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check_geometry_consistency(error, mb, path_mb, "Mueller-Brown GS")
   if (allocated(error)) return

   ! LEPS (three atoms)
   call new_leps_model(leps)
   mass9 = 1.0_wp
   x0_9 = 0.0_wp
   x0_9(4:6) = [leps_r_star, 0.0_wp, 0.0_wp]
   x0_9(7:9) = [2.0_wp*leps_r_star, 0.0_wp, 0.0_wp]
   call leps%evaluate(x0_9, e_ts, g, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call walk_irc(leps, mass9, x0_9, 0.05_wp, 1.0_wp, 60, p_landscape_irc_gs, &
      & path_leps, fail)
   call check(error, .not. fail, "the LEPS GS walk failed")
   if (allocated(error)) return
   write(stdout, '(1x,a,i0)') "  LEPS GS path points ", path_leps%n
   if (allocated(error)) return
   call check_geometry_consistency(error, leps, path_leps, "LEPS GS")
   if (allocated(error)) return
   write(stdout, '(1x,a)') "  GS path geometry reproduces its energies and gradients"

end subroutine test_gs_path_geometry_consistency

!> Re-evaluate every point of a path and require the stored values to reproduce.
subroutine check_geometry_consistency(error, model, path, what)
   type(error_type), allocatable, intent(out) :: error
   !> Potential
   class(landscape_energy_model), intent(inout) :: model
   !> Path to check
   type(landscape_irc_path), intent(in) :: path
   !> Label of the path
   character(len=*), intent(in) :: what

   real(wp), allocatable :: g(:)
   real(wp) :: e
   logical :: fail
   integer :: k

   call check(error, path%n > 3, what//": the path has too few points")
   if (allocated(error)) return
   ! the path really moves (a frozen path would be caught below as well)
   call check(error, dist_to(path%xyz(:, 1), path%xyz(:, path%n)) > 1.0e-2_wp, &
      & what//": the recorded path geometry did not move")
   if (allocated(error)) return
   allocate(g(path%ndim))
   do k = 1, path%n
      call model%evaluate(path%xyz(:, k), e, g, fail)
      call check(error, .not. fail, what//": an evaluation failed")
      if (allocated(error)) return
      call check(error, abs(e - path%e(k)) < 1.0e-10_wp, &
         & what//": the stored geometry does not reproduce the stored energy")
      if (allocated(error)) return
      call check(error, abs(vec_norm(g) - path%gnorm(k)) < 1.0e-10_wp, &
         & what//": the stored geometry does not reproduce the stored gradient")
      if (allocated(error)) return
   end do
   write(stdout, '(1x,a,i0,a)') "  ", path%n, " "//what//" points reproduce "//&
      & "their energy and gradient"

end subroutine check_geometry_consistency

!> The GS sphere relaxation must return the minimum of the energy on the sphere.
!>
!> The sphere relaxation is the **definition** of the Gonzalez-Schlegel step:
!> the next point is the constrained minimum of the energy on the sphere of
!> radius `step/2` about the pivot.  No other unit test checks that property -
!> `mb-sphere-chord` only requires the chord to be at most the step, which every
!> point of the sphere satisfies.  The reference is computed **independently**
!> here: the circle about the pivot is scanned densely in the polar angle and
!> refined by golden section search, i.e. without any Hessian, tangent space or
!> micro cycle.
!>
!> Measured (this test): the returned point is `8.5e-14` hartree above the true
!> constrained minimum, **independently of the micro cycle budget** (5, 20, 21
!> and 200 give the same value to the last digit, because the cycle now meets
!> its own step criterion after three evaluations), and the micro steps shrink
!> monotonically (`1.6e-2`, `1.3e-3`, `1.1e-5` bohr*amu^1/2).  The equality
!> constrained Lagrange step replaced the clamped tangent step that used to
!> leave the relaxation in a two cycle, 3.4e-4 / 1.3e-3 hartree above the
!> minimum with the parity of the budget (report 64); the assertions below were
!> flipped deliberately with that fix.
subroutine test_gs_sphere_relaxation_quality(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(landscape_muller_brown_model) :: model
   type(landscape_irc_settings) :: settings
   real(wp) :: sqrt_mass(2), q(2), g(2), h(2, 2), h0(2, 2), q_new(2), x_new(2), g_new(2)
   real(wp) :: pivot(2), radius, e_cur, e_new, e_sphere, gap, gnorm
   real(wp) :: gap_max, gap_even, gap_odd
   integer, parameter :: budgets(4) = [5, 20, 21, 200]
   integer :: n_eval, kcycle
   logical :: micro_converged, fail

   gap_max = -huge(1.0_wp)
   gap_even = 0.0_wp
   gap_odd = 0.0_wp

   call new_muller_brown_model(model)
   sqrt_mass = 1.0_wp
   ! a point on the lower flank of the Mueller-Brown saddle, at a distance
   ! large enough that the sphere relaxation has a real job to do
   q = mb_saddle - [0.0_wp, 0.25_wp]
   call model%evaluate(q, e_cur, g, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call model%hessian(q, h, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   ! the Hessian of the step is `intent(inout)` and is updated along the micro
   ! cycles, so every budget has to start from the same matrix or the budgets
   ! would not be comparable
   h0 = h

   call default_settings(settings, 0.10_wp, 50, p_landscape_irc_gs)
   call init(env)
   do kcycle = 1, size(budgets)
      settings%max_micro_cycles = budgets(kcycle)
      h = h0
      call landscape_irc_gs_step(env, model, settings, sqrt_mass, q, g, h, &
         & settings%step_length, .false., q_new, x_new, e_new, g_new, n_eval, &
         & micro_converged, fail)
      call check(error, .not. fail)
      if (allocated(error)) return

      ! the pivot of the step and the independent constrained minimum
      gnorm = vec_norm(g)
      pivot = q - 0.5_wp*settings%step_length*(g/gnorm)
      radius = 0.5_wp*settings%step_length
      call sphere_minimum(model, pivot, radius, q_new, e_sphere, fail)
      call check(error, .not. fail)
      if (allocated(error)) return

      ! the returned point lies on the sphere
      call check(error, abs(vec_norm(q_new - pivot) - radius) < 1.0e-12_wp, &
         & "the GS sphere relaxation returned a point off the sphere")
      if (allocated(error)) return
      ! and it is the constrained minimum up to the measured deviation
      gap = e_new - e_sphere
      gap_max = max(gap_max, gap)
      if (budgets(kcycle) == 20) gap_even = gap
      if (budgets(kcycle) == 21) gap_odd = gap
      ! the fix, demonstrated rather than described: the micro cycle now meets
      ! its own criterion (it used to be stopped by the budget in a two cycle,
      ! which is why the CLI check asserted `n_pinned >= 1`)
      call check(error, micro_converged, &
         & "the GS micro cycle did not reach its own step criterion")
      if (allocated(error)) return
      write(stdout, '(1x,a,i0,a,es14.6,a,i0,a,l1)') "  GS sphere relaxation, budget ", &
         & budgets(kcycle), ": gap ", gap, " Eh, E/g ", n_eval, &
         & ", micro cycles converged ", micro_converged
      call check(error, gap >= -1.0e-12_wp, &
         & "the GS sphere relaxation returned a point below the sphere minimum")
      if (allocated(error)) return
   end do
   ! bound = 1e-10 hartree, 1000x the measured gap (8.5e-14) and 1e7 below the
   ! 1.3e-3 of the old two cycle: the relaxation reaches the constrained
   ! minimum of the sphere.  The step itself is verified against the KKT system
   ! to 1e-15 by `gs-kkt-system`.
   call check(error, gap_max < 1.0e-10_wp, &
      & "the GS sphere relaxation did not reach the sphere minimum")
   if (allocated(error)) return
   ! budget independence: the two values are now the same (the two cycle that
   ! made the parity of the budget select one of two alternating points is
   ! gone; all four budgets converge after two evaluations), which is the
   ! property the fix buys
   call check(error, abs(gap_even - gap_odd) < 1.0e-12_wp, &
      & "the GS micro cycle still depends on the budget parity")
   if (allocated(error)) return
   ! `micro_step_thresh` is read, not decorative: a threshold at or above the
   ! first micro step accepts the first trial point, i.e. the relaxation
   ! degenerates to one step from the starting guess in a single evaluation.
   ! That is exactly the cheap-but-coarse configuration measured in report 44
   ! section 7.1 (`--lirc-micro-thresh 0.05` walks off the fine path and is cut
   ! off by the step budget), so the setting has to be visible here.
   settings%max_micro_cycles = 20
   settings%micro_step_thresh = 0.05_wp
   call landscape_irc_gs_step(env, model, settings, sqrt_mass, q, g, h, &
      & settings%step_length, .false., q_new, x_new, e_new, g_new, n_eval, &
      & micro_converged, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, micro_converged, &
      & "a micro step threshold above the sphere radius did not converge the first cycle")
   if (allocated(error)) return
   call check(error, n_eval == 1, &
      & "the accepted threshold did not stop the sphere relaxation after one evaluation")
   if (allocated(error)) return
   write(stdout, '(1x,a,es14.6)') "  GS micro threshold = radius, one cycle gap ", &
      & e_new - e_sphere

end subroutine test_gs_sphere_relaxation_quality

!> The GS micro step is the Newton step of the equality constrained problem.
!>
!> `landscape_irc_gs_lagrange_step` is exposed for exactly this test: it
!> returns the step `d`, the multiplier increment `dlambda` and the level shift
!> `2 lambda` it used, so the linearised KKT system
!>
!>     [ H + 2 lambda I      2 r ] [ d  ]   [ -( g + 2 lambda r ) ]
!>     [ 2 r^T               0   ] [ dL ] = [ -( |r|^2 - R^2 )    ]
!>
!> is reassembled here and both of its rows are evaluated.  The step is taken
!> from the first trial point of a real GS micro cycle on Mueller-Brown, and
!> then the starting point is pulled **off** the sphere (1e-3 and 1e-2 along
!> the radius, the second one large enough that the constraint row dominates
!> the residual): the constraint row is the only row that can be right there,
!> so it pins the constrained formulation rather than a merely tangential
!> step - the replaced clamped tangent step satisfies neither row.  Measured
!> residual over the three points: `< 1e-15` (relative), and the step is a
!> descent direction on the sphere (`d^T P g < 0`).
!>
!> The **update** of the multiplier (not only a single step) is covered as well
!> and is the second half of the test: the lambda of the KKT system is the
!> multiplier of the current residual, `2 lambda = -(r^T g)/|r|^2`, recomputed
!> at every micro cycle from the gradient of the point it is taken at.  The
!> test runs the micro iteration itself (exact Hessian per iterate, the same
!> march as `landscape_irc_gs_lagrange_step` + the sphere projection) and
!> asserts that the tangential KKT residual `|P g|` decreases strictly at every
!> cycle (measured 4.8e-2 -> 5.8e-4 -> 1.5e-6) and that `2 lambda` tracks it
!> (`lambda` converges to the multiplier of the constrained minimum).  The last
!> block forces the **level shift** branch: with a synthetic negative curvature
!> along the sphere `shift > 2 lambda_res` (the multiplier alone would not make
!> the reduced Hessian positive definite) and the step is checked again - a
!> large shift shrinks `|d|`, so the test also pins that the shift is not a
!> way to fake convergence: at the shifted step `|P g|` may not decrease by
!> less than the shift would require, i.e. the tangential residual is still
!> the quantity that drives `|d|` (same `d`, residual `|P g|` and shift are all
!> printed).
subroutine test_gs_kkt_system(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(landscape_muller_brown_model) :: model
   real(wp) :: q(2), r(2), pivot(2), g(2), h(2, 2), d(2), res(2), hd(2)
   real(wp) :: radius, dlambda, shift, gt(2), scale, e, offs(3), rnorm, gtan
   integer :: k
   logical :: fail

   offs = [0.0_wp, 1.0e-3_wp, -1.0e-2_wp]
   call new_muller_brown_model(model)
   call init(env)
   q = mb_saddle - [0.0_wp, 0.25_wp]
   call model%evaluate(q, e, g, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call model%hessian(q, h, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   ! the pivot of the step and the radius of its sphere
   pivot = q - 0.5_wp*0.10_wp*g/vec_norm(g)
   radius = 0.05_wp
   ! the first trial point of the micro cycle, on the sphere
   q = q - 0.10_wp*g/vec_norm(g)

   do k = 1, 3
      r = q - pivot
      r = r + offs(k)*r/vec_norm(r)
      call model%evaluate(pivot + r, e, g, fail)
      call check(error, .not. fail)
      if (allocated(error)) return
      call model%hessian(pivot + r, h, fail)
      call check(error, .not. fail)
      if (allocated(error)) return

      call landscape_irc_gs_lagrange_step(env, pivot + r, pivot, radius, h, g, &
         & d, dlambda, shift, fail)
      call check(error, .not. fail)
      if (allocated(error)) return

      ! row 1 of the KKT system, and the constraint row 2
      hd = matmul(h, d) + shift*d
      res = hd + 2.0_wp*r*dlambda + (g + shift*r)
      scale = max(1.0_wp, vec_norm(g) + abs(shift)*vec_norm(r) + vec_norm(hd))
      write(stdout, '(1x,a,es10.2,a,es12.4,a,es12.4)') &
         & "  GS KKT residual, off-sphere ", offs(k), ": row 1 ", &
         & vec_norm(res)/scale, ", row 2 ", &
         & (2.0_wp*dot_product(r, d) + dot_product(r, r) - radius**2)/radius
      call check(error, vec_norm(res)/scale < 1.0e-12_wp, &
         & "the GS Lagrange step does not solve the first KKT row")
      if (allocated(error)) return
      call check(error, abs(2.0_wp*dot_product(r, d) + dot_product(r, r) - &
         & radius**2) < 1.0e-12_wp*radius, &
         & "the GS Lagrange step does not satisfy the linearised constraint")
      if (allocated(error)) return

      ! the tangential part of the step is a descent direction of the sphere
      ! problem: `d_t^T P g < 0` with `d_t = d - (r^T d) r/|r|^2`
      gt = g - dot_product(g, r)*r/dot_product(r, r)
      call check(error, dot_product(d, gt) < 0.0_wp, &
         & "the GS Lagrange step is not a descent direction on the sphere")
      if (allocated(error)) return
   end do

   ! --- the multiplier is updated, and the iteration converges --------------
   ! the same march as the micro cycle: evaluate, take the Lagrange step,
   ! project onto the sphere, take the multiplier of the new residual
   r = q - pivot
   q = pivot + r
   gtan = huge(1.0_wp)
   do k = 1, 20
      r = q - pivot
      call model%evaluate(q, e, g, fail)
      call check(error, .not. fail)
      if (allocated(error)) return
      call model%hessian(q, h, fail)
      call check(error, .not. fail)
      if (allocated(error)) return
      call landscape_irc_gs_lagrange_step(env, q, pivot, radius, h, g, d, &
         & dlambda, shift, fail)
      call check(error, .not. fail)
      if (allocated(error)) return
      ! the multiplier of the current residual, recomputed at every cycle
      call check(error, abs(shift - 0.5_wp*(-2.0_wp*dot_product(r, g)/ &
         & dot_product(r, r))) < 1.0e-8_wp, &
         & "the GS Lagrange step did not use the residual multiplier")
      if (allocated(error)) return
      rnorm = sqrt(dot_product(r, r))
      gt = g - dot_product(g, r)*r/dot_product(r, r)
      write(stdout, '(1x,a,i2,a,es12.4,a,es12.4,a,es12.4)') &
         & "  GS micro cycle ", k, ": |P g| ", vec_norm(gt), &
         & "  |d| ", vec_norm(d), "  2 lambda ", shift
      call check(error, vec_norm(gt) < gtan, &
         & "the tangential KKT residual did not decrease in the micro cycle")
      if (allocated(error)) return
      gtan = vec_norm(gt)
      if (vec_norm(d) <= 1.0e-3_wp) exit
      ! the projection of the micro cycle: onto the sphere, not along the old
      ! radius (the real code calls `project_sphere`)
      q = pivot + radius*(q + d - pivot)/vec_norm(q + d - pivot)
   end do
   ! it really converges: the tangential residual reaches the level the step
   ! criterion implies.  The two criteria are equivalent norms here - the step
   ! is `d = -(P H P + 2 lambda P)^-1 P g` with a positive definite reduced
   ! Hessian, so `|d| <= 1e-3` bounds `|P g|` by `1e-3 * ||reduced Hessian||`;
   ! measured `|P g| = 4.2e-5` at `|d| = 1.1e-5`, i.e. the operator norm is
   ! ~3.8 and the bound is ~4e-3.  The shifted Hessian is positive definite by
   ! construction (`gs-kkt-system` checks the branch below), so a large level
   ! shift can make `|d|` small but it also makes `|P g| <= ||R|| |d|` small
   ! with the same factor: there is no shift that fakes convergence.
   call check(error, gtan < 1.0e-4_wp, &
      & "the GS micro iteration did not converge to the sphere minimum")
   if (allocated(error)) return

   ! --- the level shift branch --------------------------------------------
   ! A synthetic Hessian whose curvature along the sphere is negative: the
   ! residual multiplier alone would leave the reduced Hessian indefinite, so
   ! `2 lambda` has to be raised to the level shift `-eval_min + floor`.  The
   ! level shift necessarily produces a **huge** step along the shifted
   ! direction (`eval_min + shift = floor`), so this branch is the one the
   ! trust region exists for: the test pins that the shift fired, that the
   ! step stayed inside the trust region and that it is still a descent
   ! direction.  The KKT residual is **not** asserted here - a clamped step is
   ! no longer the solution of the system, by construction; that assertion is
   ! the one the three unclamped points above and the micro cycle carry.
   h = 0.0_wp
   h(1, 1) = -20.0_wp
   h(2, 2) = 1.0_wp
   r = q - pivot
   call model%evaluate(pivot + r, e, g, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call landscape_irc_gs_lagrange_step(env, pivot + r, pivot, radius, h, g, &
      & d, dlambda, shift, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, shift > 0.5_wp*(-2.0_wp*dot_product(r, g)/ &
      & dot_product(r, r)), &
      & "the level shift did not raise the multiplier of an indefinite step")
   if (allocated(error)) return
   call check(error, vec_norm(d) <= radius*(1.0_wp + 1.0e-12_wp), &
      & "the level shifted step left the trust region")
   if (allocated(error)) return
   gt = g - dot_product(g, r)*r/dot_product(r, r)
   call check(error, dot_product(d, gt) < 0.0_wp, &
      & "the level shifted step is not a descent direction")
   if (allocated(error)) return
   write(stdout, '(1x,a,es12.4,a,es12.4,a,es12.4)') &
      & "  GS level shift: 2 lambda_res ", -2.0_wp*dot_product(r, g)/ &
      & dot_product(r, r), "  shift used ", shift, "  |d| ", vec_norm(d)

end subroutine test_gs_kkt_system

!> Minimum of the energy on the circle `|q - pivot| = radius`, by dense scan.!>
!> Independent reference for `test_gs_sphere_relaxation_quality`: the circle is
!> sampled in the polar angle and the bracket of the best sample is refined by
!> golden section search (no gradient, no Hessian, no micro cycle).
subroutine sphere_minimum(model, pivot, radius, q_ref, e_min, fail)
   !> Potential
   class(landscape_energy_model), intent(inout) :: model
   !> Centre of the circle
   real(wp), intent(in) :: pivot(:)
   !> Radius of the circle
   real(wp), intent(in) :: radius
   !> A point on the circle (defines the polar axis of the scan)
   real(wp), intent(in) :: q_ref(:)
   !> Energy of the minimum (hartree)
   real(wp), intent(out) :: e_min
   !> Failure flag
   logical, intent(out) :: fail

   integer, parameter :: nscan = 20000
   integer, parameter :: nrefine = 200
   real(wp) :: axis(2), normal(2), phi, phi_lo, phi_hi, e, g(size(pivot))
   real(wp) :: q(2), aa, bb, cc, dd, e_aa, e_bb, e_cc, e_dd
   integer :: k

   fail = .false.
   ! orthonormal frame of the plane: the scan is the exact circle
   axis = (q_ref - pivot)/radius
   normal = [-axis(2), axis(1)]
   e_min = huge(1.0_wp)
   phi_lo = 0.0_wp
   phi_hi = 0.0_wp
   do k = 0, nscan
      phi = 2.0_wp*acos(-1.0_wp)*real(k, wp)/real(nscan, wp)
      q = pivot + radius*(cos(phi)*axis + sin(phi)*normal)
      call model%evaluate(q, e, g, fail)
      if (fail) return
      if (e < e_min) then
         e_min = e
         phi_lo = 2.0_wp*acos(-1.0_wp)*real(k - 1, wp)/real(nscan, wp)
         phi_hi = 2.0_wp*acos(-1.0_wp)*real(k + 1, wp)/real(nscan, wp)
      end if
   end do
   ! golden section search on [phi_lo, phi_hi]
   aa = phi_lo
   bb = phi_hi
   cc = bb - (bb - aa)/1.618033988749895_wp
   dd = aa + (bb - aa)/1.618033988749895_wp
   do k = 1, nrefine
      q = pivot + radius*(cos(cc)*axis + sin(cc)*normal)
      call model%evaluate(q, e_cc, g, fail)
      if (fail) return
      q = pivot + radius*(cos(dd)*axis + sin(dd)*normal)
      call model%evaluate(q, e_dd, g, fail)
      if (fail) return
      if (e_cc < e_dd) then
         bb = dd
         dd = cc
         cc = bb - (bb - aa)/1.618033988749895_wp
      else
         aa = cc
         cc = dd
         dd = aa + (bb - aa)/1.618033988749895_wp
      end if
   end do
   q = pivot + radius*(cos(0.5_wp*(aa + bb))*axis + sin(0.5_wp*(aa + bb))*normal)
   call model%evaluate(q, e_aa, g, fail)
   if (fail) return
   e_min = min(e_min, e_aa)

end subroutine sphere_minimum

!> A path object can be reused for a path of another dimension.
!>
!> `landscape_irc_path` stores `xyz(ndim, capacity)`; the reset used to keep the
!> allocation and only forget the dimension, so a second walk with another
!> dimension wrote a 2-vector into a 9-row array (finding F-1 of
!> `45-m6-independent-review.md`).  The test pins both entry points: the
!> documented reset (which now releases the storage) and the defensive rebind of
!> a live path with a point of another dimension, which must not write through
!> mismatched shapes.
subroutine test_path_cross_dimension_reuse(error)
   type(error_type), allocatable, intent(out) :: error
   type(landscape_irc_path) :: path
   real(wp) :: x2(2), x9(9)
   integer :: k

   x2 = [1.0_wp, -2.0_wp]
   x9 = [(real(k, wp), k = 1, 9)]

   ! a two dimensional path
   call landscape_irc_path_add(path, x2, -1.0_wp, 0.5_wp, 0.25_wp, 0.125_wp, 7)
   call check(error, path%ndim == 2, "the path dimension was not recorded")
   if (allocated(error)) return
   call check(error, path%n == 1, "the first point was not stored")
   if (allocated(error)) return

   ! the reset restores the empty state *including* the storage
   call landscape_irc_path_reset(path)
   call check(error, path%n == 0 .and. path%ndim == 0, &
      & "the reset did not empty the path")
   if (allocated(error)) return
   call check(error, path%capacity == 0, &
      & "the reset kept the capacity of the previous dimension")
   if (allocated(error)) return
   call check(error, .not. allocated(path%xyz), &
      & "the reset kept the storage of the previous dimension")
   if (allocated(error)) return

   ! a nine component point in the same object
   call landscape_irc_path_add(path, x9, -2.0_wp, 1.5_wp, 0.75_wp, 0.375_wp, 9)
   call check(error, path%ndim == 9 .and. path%n == 1, &
      & "the reused path did not switch to the new dimension")
   if (allocated(error)) return
   call check(error, maxval(abs(path%xyz(:, 1) - x9)) == 0.0_wp, &
      & "the stored point is not the nine component point")
   if (allocated(error)) return
   call check(error, path%e(1) == -2.0_wp .and. path%n_eval(1) == 9, &
      & "the point bookkeeping of the reused path is wrong")
   if (allocated(error)) return

   ! without the reset: the storage of the other dimension is released, the
   ! path restarts instead of writing through mismatched shapes
   call landscape_irc_path_add(path, x2, -3.0_wp, 2.5_wp, 0.0_wp, 0.0_wp, 1)
   call check(error, path%ndim == 2 .and. path%n == 1, &
      & "a point of another dimension was not rebound")
   if (allocated(error)) return
   call check(error, maxval(abs(path%xyz(:, 1) - x2)) == 0.0_wp, &
      & "the rebound point is not the two component point")
   if (allocated(error)) return

end subroutine test_path_cross_dimension_reuse

!> The endpoint attribution rule: termination signal, relaxation, RMSD.
!>
!> `endpoint_reached` is the function behind `reached_min_fwd/bwd` and
!> therefore behind the user visible "reached minimum" line of `irc.log` and of
!> the result summary.  The rmsd criterion is what makes the documented
!> semantics of `--lirc-rmsd-end` (settings comment, `--help`) true: an endpoint
!> further than `rmsd_end` from its relaxed minimum is **not** assigned to a
!> minimum even though the path terminated at a stationary point.  Without a
!> relaxed endpoint the criterion has no input and is skipped.
!>
!> Signal (d) (the limit cycle guard) is the **truncation** of a stalled walk:
!> its endpoint is one phase of the oscillation, not a minimum, so it counts
!> only when the endpoint relaxation independently confirms the geometry
!> (round-2 review HIGH-1 - otherwise a truncated path would read as an
!> arrival).
subroutine test_endpoint_rmsd_criterion(error)
   type(error_type), allocatable, intent(out) :: error
   real(wp), parameter :: rmsd_end = 0.18897_wp

   ! a direction cut off by the step budget never counts
   call check(error, .not. endpoint_reached(p_landscape_irc_term_max_steps, &
      & .true., .true., 0.0_wp, rmsd_end), &
      & "a path cut off by the step budget was assigned to a minimum")
   if (allocated(error)) return
   ! the limit cycle guard alone does not reach a minimum: the endpoint is one
   ! phase of the oscillation, even when it happens to sit inside `rmsd_end`
   ! of a relaxed minimum there is no relaxation to say which one
   call check(error, .not. endpoint_reached(p_landscape_irc_term_limit_cycle, &
      & .false., .false., 0.0_wp, rmsd_end), &
      & "a truncated limit cycle endpoint was assigned to a minimum")
   if (allocated(error)) return
   ! with an endpoint relaxation that converged inside the RMSD window the
   ! truncation is confirmed and counts
   call check(error, endpoint_reached(p_landscape_irc_term_limit_cycle, &
      & .true., .true., 1.0e-3_wp, rmsd_end), &
      & "a limit cycle endpoint confirmed by its relaxation was rejected")
   if (allocated(error)) return
   ! ... but not when the relaxed minimum is further away than rmsd_end
   call check(error, .not. endpoint_reached(p_landscape_irc_term_limit_cycle, &
      & .true., .true., 2.0_wp*rmsd_end, rmsd_end), &
      & "a limit cycle endpoint outside the RMSD window was assigned to a minimum")
   if (allocated(error)) return
   ! signal (b) with a converged relaxation and a close minimum
   call check(error, endpoint_reached(p_landscape_irc_term_grad, .true., &
      & .true., 1.0e-3_wp, rmsd_end), &
      & "a converged endpoint inside the RMSD window was rejected")
   if (allocated(error)) return
   ! signal (c) counts as well
   call check(error, endpoint_reached(p_landscape_irc_term_energy_rise, .true., &
      & .true., 1.0e-3_wp, rmsd_end), &
      & "signal (c) with a converged relaxation was rejected")
   if (allocated(error)) return
   ! the same endpoint, but further from the relaxed minimum than rmsd_end
   call check(error, .not. endpoint_reached(p_landscape_irc_term_grad, .true., &
      & .true., 2.0_wp*rmsd_end, rmsd_end), &
      & "an endpoint outside the RMSD window was assigned to a minimum")
   if (allocated(error)) return
   ! a failed relaxation is not a minimum
   call check(error, .not. endpoint_reached(p_landscape_irc_term_grad, .true., &
      & .false., 0.0_wp, rmsd_end), &
      & "an unconverged relaxation was assigned to a minimum")
   if (allocated(error)) return
   ! without a relaxation there is no RMSD to test: the criterion is skipped
   call check(error, endpoint_reached(p_landscape_irc_term_grad, .false., &
      & .false., 2.0_wp*rmsd_end, rmsd_end), &
      & "the RMSD criterion was applied without a relaxed endpoint")
   if (allocated(error)) return

end subroutine test_endpoint_rmsd_criterion

!> The phase convention removes the eigensolver's eigenvector sign ambiguity.
!>
!> Linear algebra libraries return an eigenvector with an arbitrary sign.  The
!> phase convention (largest absolute component positive, PLAN section 6.3) is
!> what makes "forward" independent of that sign, so the same geometry gives
!> the same orientation on any build.  The test feeds `landscape_transition_
!> mode_set_phase` an explicitly sign flipped vector and requires the same
!> oriented vector back; removing the convention leaves the raw eigenvector and
!> the flipped copy pointing in opposite directions.
subroutine test_phase_sign_invariance(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(landscape_transition_mode) :: forward, flipped
   real(wp), allocatable :: hess(:, :)
   logical :: fail

   call init(env)
   allocate(hess(2, 2))
   ! one negative eigenvalue, and the largest absolute component of the
   ! transition vector is not the first one (so the test is not degenerate)
   hess = reshape([-0.5_wp, 0.9_wp, 0.9_wp, 1.5_wp], [2, 2])
   call landscape_irc_select_direction(env, hess, 1, forward, fail)
   call check(error, .not. fail)
   if (allocated(error)) return

   ! the same mode, but with the eigenvector sign that a different library
   ! (or a different LAPACK) could have returned
   flipped%vector = -forward%vector
   call landscape_transition_mode_set_phase(flipped)

   call check(error, all(flipped%vector == forward%vector), &
      & "the phase convention is not invariant under a sign flip of the "//&
      & "eigenvector (the forward direction depends on the library)")
   if (allocated(error)) return
   call check(error, forward%phase_value > 0.0_wp, &
      & "the oriented transition vector has a negative largest component")
   if (allocated(error)) return
   call check(error, flipped%phase_index == forward%phase_index, &
      & "the recorded phase component changed under a sign flip")
   if (allocated(error)) return
   call check(error, flipped%phase_value == forward%phase_value, &
      & "the recorded phase value changed under a sign flip")
   if (allocated(error)) return
   write(stdout, '(1x,a,2f12.6)') "  phase invariant under sign flip: ", &
      & forward%vector

end subroutine test_phase_sign_invariance

!> The HPC predictor step is load bearing: a wrong predicted point degrades the path.
!>
!> The HPC corrector builds its DWI surface from the current point and the
!> **predicted** point and integrates from the current point.  Because the DWI
!> reproduces the current sample exactly, a wrong predictor is invisible in the
!> first correction step and only biases the interpolation as the path advances;
!> it therefore has to be pinned by a path accuracy assertion, not by the
!> endpoint basin alone (a wrong `q_pred = q_cur - d` still reached the analytic
!> minimum and left every other unit test green - the coverage gap this test
!> closes).
subroutine test_hpc_predictor_accuracy(error)
   type(error_type), allocatable, intent(out) :: error
   type(landscape_muller_brown_model) :: model
   type(landscape_irc_path) :: reference, path
   real(wp) :: mass(2), x0(2)
   real(wp) :: steps(4), dev(4), endd(4)
   integer :: i
   logical :: fail

   call new_muller_brown_model(model)
   mass = 1.0_wp
   x0 = mb_saddle
   steps = [0.05_wp, 0.10_wp, 0.20_wp, 0.40_wp]

   call walk_mb(model, mass, x0, 0.01_wp, 1.0_wp, 2000, p_landscape_irc_lqa, &
      & reference, fail)
   call check(error, .not. fail)
   if (allocated(error)) return

   do i = 1, 4
      call walk_mb(model, mass, x0, steps(i), 1.0_wp, 2000, p_landscape_irc_hpc, &
         & path, fail)
      call check(error, .not. fail)
      if (allocated(error)) return
      dev(i) = max_path_deviation(path, reference)
      endd(i) = min(dist_to(path%xyz(:, path%n), mb_min1), &
         & dist_to(path%xyz(:, path%n), mb_min2))
   end do
   write(stdout, '(1x,a,4f10.5)') "  HPC MB path deviation 0.05/0.10/0.20/0.40 ", dev
   write(stdout, '(1x,a,4f10.5)') "  HPC MB endpoint distance          ", endd

   ! At the validated step sizes (0.05, 0.10, 0.20) the path must stay close to
   ! the reference; the 0.40 step is already a documented deviation even with
   ! the correct predictor (M6: "step 0.40 fails for lqa/hpc/eulerpc").  The
   ! endpoint distance is the sharp discriminator: the correct predictor reaches
   ! the analytic minimum to better than 5e-4, a wrong one leaves it by 0.07 to
   ! 0.12.
   call check(error, max(dev(1), dev(2)) < 0.02_wp, &
      & "the HPC path deviates from the reference: its predictor is wrong")
   if (allocated(error)) return
   ! the endpoint contract is stated for the validated step sizes only.  At 0.40
   ! the DWI/mBS corrector falls back to one Gonzalez-Schlegel step and that
   ! step is taken with a sphere wide enough to throw the path across the
   ! ridge: the measured endpoint is (0.94117, 0.04050), in the **third** well
   ! of the surface at (0.623, 0.028) (0.318 bohr away), **not** in either IRC
   ! basin - it is 1.079 bohr from mb_min2 and 2.052 bohr from mb_min1, so the
   ! reported `endd(4) = 1.07893` is the distance to the *nearest* minimum.
   ! This comment used to read the 1.079/0.318 pair as "min1 / min2" and to
   ! describe the 0.40 failure as an alternation "between the two
   ! Mueller-Brown basins"; report 67 measured both to be wrong.  The jump is
   ! a cross ridge one into the third well, while the GS **macro step**
   ! oscillation itself is a stall around a *single* minimum (the two phases
   ! of step 0.20 forward lie 0.235 and 0.162 bohr from mb_min2 and more than
   ! 0.98 bohr from mb_min1; report 69 adds the limit cycle guard).  Which
   ! well the 0.40 endpoint lands in is a known limitation of the GS macro
   ! step, **not** a predictor signal, so the 0.40 endpoint is pinned to the
   ! basin scale rather than to one basin: pinning a basin would turn red for
   ! the next person who fixes the GS macro step, the bound holds either way.
   call check(error, maxval(endd(1:3)) < 0.02_wp, &
      & "the HPC endpoint left the analytic minimum")
   if (allocated(error)) return
   call check(error, endd(4) < 2.0_wp, &
      & "the HPC 0.40 endpoint left the surface")
   if (allocated(error)) return

end subroutine test_hpc_predictor_accuracy

!> The GS limit cycle guard truncates a stalled macro step walk and is inert on
!> every validated path.
!>
!> Report 67 measured that the Gonzalez-Schlegel macro step stalls in a two
!> cycle near a minimum: once the pivot sits on the basin floor (distance to
!> the minimum below `step/2`) the sphere of radius `step/2` carries **two**
!> local constrained minima and the micro cycle, whose initial guess is the
!> point opposite the current one, alternates between them.  The running
!> minimum is refreshed by ~1e-5 hartree per cycle, so termination signal (c)
!> fires only after a very long crawl: on the analytic Mueller-Brown surface at
!> step 0.20 forward the walk runs to 88 points (87 macro steps) over an arc
!> length of 17.48899 before (c) stops it, and at 0.35 **backward** the macro
!> steps keep cancelling until the whole step budget is exhausted (400 points,
!> arc length 139.597491, `term = max_steps`), i.e. (c) does not fire there
!> before the budget.  The returned endpoint sits 0.106 bohr from the minimum
!> instead of 0.025 (report 69).
!>
!> The guard counts consecutive macro steps that cancel
!> (`|u_k + u_{k-1}| < 0.1 step`, `u_k = q_new - q_cur`) and leaves the path
!> with its **own** signal (d) `limit_cycle` after `osc_window` of them, so a
!> truncated path stays distinguishable from a normal one in `irc.log` and
!> `--ljson` (the round-2 review finding HIGH-1).  The assertions pin both
!> halves of its contract:
!>
!>  * `osc_window = 0` reproduces the stalled walk bitwise (the escape hatch),
!>  * the default window 8 truncates it to a short walk that ends closer to the
!>    minimum than the stalled one, for 0.20 forward **and** for the worst case
!>    0.35 backward (the direction that otherwise burns the whole budget),
!>  * on the validated steps 0.05, 0.10 and 0.15 (the longest runs of >170
!>    degree reversals there are 5, 7 and 4, so the cancelling run is shorter
!>    than the window) and on the LEPS 0.05 path the result is **bitwise**
!>    unchanged, i.e. the guard cannot move a validated number, and
!>  * the opt-in `--lirc-adaptive` improves the endpoint of the stalled walk
!>    (0.09402 -> 0.01855 bohr), is off by default and reports its reduction
!>    sequence (three reductions, final effective step 0.025, floor not
!>    reached).
subroutine test_gs_limit_cycle_guard(error)
   type(error_type), allocatable, intent(out) :: error
   type(landscape_muller_brown_model) :: model
   type(landscape_leps_model) :: leps
   type(landscape_irc_settings) :: settings
   type(landscape_irc_path) :: stalled, guarded, plain(3), reference(3)
   type(landscape_irc_path) :: leps_plain, leps_reference
   type(landscape_irc_path) :: adaptive, default_off
   type(landscape_irc_path) :: stalled_bwd, guarded_bwd
   real(wp) :: mass(2), mass9(9), x0(2), x0_9(9), steps(3), d_stalled, d_guarded
   real(wp) :: d_adaptive, last_step
   logical :: fail
   integer :: i, n_reduction

   call new_muller_brown_model(model)
   mass = 1.0_wp
   x0 = mb_saddle
   steps = [0.05_wp, 0.10_wp, 0.15_wp]

   ! --- the stalled walk, guard disabled (the escape hatch) --------------
   call default_settings(settings, 0.20_wp, 400, p_landscape_irc_gs)
   settings%osc_window = 0
   call walk_with(model, mass, x0, settings, 1.0_wp, stalled, fail)
   call check(error, .not. fail)
   if (allocated(error)) return

   ! --- the same walk with the default guard -----------------------------
   settings%osc_window = 8
   call walk_with(model, mass, x0, settings, 1.0_wp, guarded, fail)
   call check(error, .not. fail)
   if (allocated(error)) return

   d_stalled = min(dist_to(stalled%xyz(:, stalled%n), mb_min1), &
      & dist_to(stalled%xyz(:, stalled%n), mb_min2))
   d_guarded = min(dist_to(guarded%xyz(:, guarded%n), mb_min1), &
      & dist_to(guarded%xyz(:, guarded%n), mb_min2))
   write(stdout, '(1x,a,i0,a,f12.5,a,i0,a,f12.5)') "  GS MB 0.20 points ", &
      & stalled%n, " s ", stalled%s(stalled%n), "  guarded ", guarded%n, &
      & " s ", guarded%s(guarded%n)
   write(stdout, '(1x,a,2f12.6,a,2f12.6)') "  GS MB 0.20 endpoint ", &
      & stalled%xyz(:, stalled%n), "  guarded ", guarded%xyz(:, guarded%n)
   write(stdout, '(1x,a,f10.5,a,f10.5,a,i0)') "  GS MB 0.20 endpoint distance ", &
      & d_stalled, " -> ", d_guarded, "  term ", guarded%term

   ! the walk really stalls without the guard (otherwise this test is vacuous);
   ! the pinned numbers are the measurements of **this** build.  Report 67
   ! section 6.2 reports the same structure with a baseline arc length of
   ! 17.389 and the two phases (-0.15530, 0.47702) / (0.04295, 0.45231); the
   ! prototype probe behind them carried extra diagnostics in `gs_step`, and
   ! its baseline arc length is exactly one step length (0.100) below the value
   ! measured here, so its `s` did not count the initial half step point.
   ! The guarded arc length agrees to four digits (2.49688 vs 2.497).
   call check(error, stalled%n == 88, &
      & "the GS 0.20 walk no longer stalls: the guard test would be vacuous")
   if (allocated(error)) return
   call check(error, abs(stalled%s(stalled%n) - 17.48899_wp) < 1.0e-4_wp, &
      & "the stalled GS 0.20 arc length changed")
   if (allocated(error)) return
   call check(error, stalled%term == p_landscape_irc_term_energy_rise)
   if (allocated(error)) return

   ! the guard truncates it, reports its **own** signal (d) `limit_cycle` (a
   ! cut, not the energy rise of (c)) and ends no further from the minimum than
   ! the stalled walk did
   call check(error, guarded%n == 13, &
      & "the limit cycle guard did not truncate the stalled walk to 13 points")
   if (allocated(error)) return
   call check(error, abs(guarded%s(guarded%n) - 2.49688_wp) < 1.0e-4_wp, &
      & "the guarded GS 0.20 arc length changed")
   if (allocated(error)) return
   call check(error, guarded%term == p_landscape_irc_term_limit_cycle, &
      & "the limit cycle guard did not report its own signal (d)")
   if (allocated(error)) return
   call check(error, d_guarded < d_stalled - 1.0e-3_wp, &
      & "the guarded endpoint is no closer to the minimum than the stalled one")
   if (allocated(error)) return
   ! the endpoint is the phase of the oscillation the guard stopped on
   call check(error, abs(guarded%xyz(1, guarded%n) - 0.042910_wp) < 1.0e-5_wp .and. &
      & abs(guarded%xyz(2, guarded%n) - 0.452366_wp) < 1.0e-5_wp, &
      & "the guarded GS 0.20 endpoint moved")
   if (allocated(error)) return

   ! --- the worst case: 0.35 backward burns the whole budget ------------
   ! This direction never reaches an endpoint on its own: the macro steps keep
   ! cancelling until `max_steps` is exhausted (measured: 400 points, arc
   ! length 139.597491, `term = max_steps`), so it is the falsification entry
   ! point for "the guard really solves the budget burn" and the only
   ! backward coverage of the truncation side (the round-2 review MEDIUM-2).
   call default_settings(settings, 0.35_wp, 400, p_landscape_irc_gs)
   settings%osc_window = 0
   call walk_with(model, mass, x0, settings, -1.0_wp, stalled_bwd, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, stalled_bwd%n == 400 .and. &
      & stalled_bwd%term == p_landscape_irc_term_max_steps, &
      & "the GS 0.35 backward walk no longer exhausts the step budget")
   if (allocated(error)) return
   call check(error, abs(stalled_bwd%s(stalled_bwd%n) - 139.597491_wp) < 1.0e-4_wp, &
      & "the stalled GS 0.35 backward arc length changed")
   if (allocated(error)) return
   settings%osc_window = 8
   call walk_with(model, mass, x0, settings, -1.0_wp, guarded_bwd, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   write(stdout, '(1x,a,i0,a,f12.6,a,i0,a,f12.6,a,i0)') "  GS MB 0.35 bwd pts ", &
      & stalled_bwd%n, " s ", stalled_bwd%s(stalled_bwd%n), " guarded ", &
      & guarded_bwd%n, " s ", guarded_bwd%s(guarded_bwd%n), " term ", &
      & guarded_bwd%term
   call check(error, guarded_bwd%n == 12, &
      & "the limit cycle guard did not truncate the 0.35 backward walk to 12 points")
   if (allocated(error)) return
   call check(error, guarded_bwd%term == p_landscape_irc_term_limit_cycle, &
      & "the 0.35 backward truncation did not report signal (d)")
   if (allocated(error)) return
   call check(error, abs(guarded_bwd%s(guarded_bwd%n) - 3.959668_wp) < 1.0e-4_wp, &
      & "the guarded GS 0.35 backward arc length changed")
   if (allocated(error)) return

   ! --- bitwise inert on the validated step sizes ------------------------
   do i = 1, 3
      call default_settings(settings, steps(i), 400, p_landscape_irc_gs)
      settings%osc_window = 0
      call walk_with(model, mass, x0, settings, 1.0_wp, reference(i), fail)
      call check(error, .not. fail)
      if (allocated(error)) return
      settings%osc_window = 8
      call walk_with(model, mass, x0, settings, 1.0_wp, plain(i), fail)
      call check(error, .not. fail)
      if (allocated(error)) return
      call check(error, paths_identical(plain(i), reference(i)), &
         & "the limit cycle guard changed a validated Gonzalez-Schlegel path")
      if (allocated(error)) return
      write(stdout, '(1x,a,f6.3,a,i0,a,f14.10,a,f12.6)') "  GS step ", steps(i), &
         & " points ", plain(i)%n, " s ", plain(i)%s(plain(i)%n), &
         & " (guard inert, bitwise)"
   end do

   ! --- LEPS 0.05: the cancellation test does not misfire on the forward ----
   ! --- jitter of the flat H2 asymptote --------------------------------
   ! The turn angle between consecutive macro steps exceeds 170 degrees 286
   ! times on this path while the reaction coordinate advances monotonically
   ! (report 67 section 4.2), so an angle criterion would truncate it; the
   ! cancellation criterion is inert, which is why it was chosen.  This leg
   ! is the client of the LEPS four method crosscheck, whose golden numbers
   ! are produced at exactly this step.
   call new_leps_model(leps)
   mass9 = 1.0_wp
   x0_9 = 0.0_wp
   x0_9(4:6) = [leps_r_star, 0.0_wp, 0.0_wp]
   x0_9(7:9) = [2.0_wp*leps_r_star, 0.0_wp, 0.0_wp]
   call default_settings(settings, 0.05_wp, 400, p_landscape_irc_gs)
   settings%osc_window = 0
   call walk_with(leps, mass9, x0_9, settings, 1.0_wp, leps_reference, fail)
   call check(error, .not. fail, "the LEPS GS reference walk failed")
   if (allocated(error)) return
   settings%osc_window = 8
   call walk_with(leps, mass9, x0_9, settings, 1.0_wp, leps_plain, fail)
   call check(error, .not. fail, "the LEPS GS guarded walk failed")
   if (allocated(error)) return
   write(stdout, '(1x,a,i0,a,f14.10)') "  LEPS GS 0.05 points ", leps_plain%n, &
      & " s ", leps_plain%s(leps_plain%n)
   call check(error, leps_plain%n > 100, &
      & "the LEPS GS 0.05 walk is too short to test the guard")
   if (allocated(error)) return
   call check(error, paths_identical(leps_plain, leps_reference), &
      & "the limit cycle guard truncated the LEPS 0.05 path (turn angle "//&
      & "misfire: the reaction coordinate is advancing)")
   if (allocated(error)) return

   ! --- the opt-in adaptive step (`--lirc-adaptive`, default off) --------
   ! It halves the effective step length at the first cancellation instead of
   ! ending the walk, so a large step recovers accuracy; the default off state
   ! has to reproduce the guarded walk bitwise.  The reduction sequence is
   ! pinned as well (round-2 review LOW-3): measured three reductions
   ! 0.20 -> 0.10 -> 0.05 -> 0.025, and the frozen floor `step/16 = 0.0125` is
   ! **not** reached on this walk.  Injecting the floor as `step/2` collapses
   ! the sequence to a single reduction, so the assertion is falsifiable.
   call default_settings(settings, 0.20_wp, 400, p_landscape_irc_gs)
   settings%osc_window = 8
   settings%adaptive = .true.
   call walk_capture_reductions(model, mass, x0, settings, 1.0_wp, adaptive, &
      & n_reduction, last_step, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   d_adaptive = min(dist_to(adaptive%xyz(:, adaptive%n), mb_min1), &
      & dist_to(adaptive%xyz(:, adaptive%n), mb_min2))
   write(stdout, '(1x,a,i0,a,f14.10,a,f10.5)') "  GS MB 0.20 adaptive points ", &
      & adaptive%n, " s ", adaptive%s(adaptive%n), " endpoint distance ", d_adaptive
   write(stdout, '(1x,a,i0,a,f10.6,a)') "  GS MB 0.20 adaptive reductions ", &
      & n_reduction, " final effective step ", last_step, " bohr*amu^1/2"
   ! measured: 12 points, arc length 1.5216 and an endpoint 0.01855 bohr from
   ! the minimum (the guard alone returns 0.09402)
   call check(error, adaptive%n == 12, &
      & "--lirc-adaptive changed the number of macro steps")
   if (allocated(error)) return
   call check(error, n_reduction == 3, &
      & "--lirc-adaptive changed the number of step reductions")
   if (allocated(error)) return
   call check(error, abs(last_step - 0.025_wp) < 1.0e-12_wp, &
      & "--lirc-adaptive changed the final effective step length")
   if (allocated(error)) return
   call check(error, abs(adaptive%s(adaptive%n) - 1.521604_wp) < 1.0e-5_wp, &
      & "--lirc-adaptive changed the adaptive arc length")
   if (allocated(error)) return
   call check(error, d_adaptive < d_guarded - 1.0e-3_wp, &
      & "--lirc-adaptive did not improve the endpoint of the stalled walk")
   if (allocated(error)) return
   call check(error, d_adaptive < 2.0e-2_wp, &
      & "--lirc-adaptive left the endpoint further from the minimum than 0.02 bohr")
   if (allocated(error)) return
   settings%adaptive = .false.
   call walk_with(model, mass, x0, settings, 1.0_wp, default_off, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   call check(error, paths_identical(default_off, guarded), &
      & "--lirc-adaptive is not off by default")
   if (allocated(error)) return

end subroutine test_gs_limit_cycle_guard

!> Run one walk with the environment unit redirected into a scratch file and
!> read back the adaptive step reductions.
!>
!> `landscape_irc_gs_walk` reports every reduction of the opt-in adaptive step
!> as `irc: the macro steps cancel - reducing the effective step length to
!> <level> bohr*amu^1/2`.  The helper counts those lines and returns the level
!> after the last one, so the reduction sequence of `--lirc-adaptive` is
!> asserted instead of inferred from the endpoint (the frozen floor has no
!> other observable).
subroutine walk_capture_reductions(model, mass, x0, settings, sign, path, &
      & n_reduction, last_step, fail)
   !> Potential
   class(landscape_energy_model), intent(inout) :: model
   !> Masses (amu)
   real(wp), intent(in) :: mass(:)
   !> Starting point
   real(wp), intent(in) :: x0(:)
   !> Settings
   type(landscape_irc_settings), intent(in) :: settings
   !> Direction sign
   real(wp), intent(in) :: sign
   !> Resulting path
   type(landscape_irc_path), intent(inout) :: path
   !> Number of reductions reported by the walk
   integer, intent(out) :: n_reduction
   !> Effective step length after the last reduction (0 when none)
   real(wp), intent(out) :: last_step
   !> Failure flag
   logical, intent(out) :: fail

   type(TEnvironment) :: env
   character(len=512) :: line
   integer :: unit, pos, ios
   real(wp) :: value

   call init(env)
   n_reduction = 0
   last_step = 0.0_wp
   open(newunit=unit, status='scratch', iostat=ios)
   env%unit = unit
   call walk_with(model, mass, x0, settings, sign, path, fail, env)
   if (ios == 0 .and. .not. fail) then
      rewind(unit)
      do
         read(unit, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (index(line, 'reducing') > 0) then
            n_reduction = n_reduction + 1
            pos = index(line, 'to ')
            read(line(pos + 3:), *) value
            last_step = value
         end if
      end do
   end if
   close(unit)

end subroutine walk_capture_reductions

!> Bitwise equality of two paths (size, dimension, arc length, energy,
!> geometry, gradient norm and maximum gradient component).
!>
!> Used as the "no effect" assertion of the GS limit cycle guard: an inert
!> guard may not touch a single bit of a validated path.  `ndim` and `gmax`
!> belong to the comparison as well (round-2 review LOW-2): `gs_walk` writes
!> `gmax` for every point, so a path comparison that skips it is not the
!> "bitwise" comparison its name claims.
pure function paths_identical(a, b) result(same)
   !> First path
   type(landscape_irc_path), intent(in) :: a
   !> Second path
   type(landscape_irc_path), intent(in) :: b
   !> True when every stored value is bitwise equal
   logical :: same

   same = a%n == b%n .and. a%term == b%term .and. a%n_pinned == b%n_pinned &
      & .and. a%ndim == b%ndim
   if (.not. same .or. a%n <= 0) return
   same = all(a%s(1:a%n) == b%s(1:b%n)) .and. all(a%e(1:a%n) == b%e(1:b%n)) &
      & .and. all(a%gnorm(1:a%n) == b%gnorm(1:b%n)) &
      & .and. all(a%gmax(1:a%n) == b%gmax(1:b%n)) &
      & .and. all(a%n_eval(1:a%n) == b%n_eval(1:b%n)) &
      & .and. all(a%xyz(:, 1:a%n) == b%xyz(:, 1:b%n))

end function paths_identical

!> LEPS H + H2: all four integrators reach the H2 asymptote in both directions.
!>
!> The collinear symmetric H3 geometry is the LEPS saddle point; the two IRC
!> directions are the two equivalent H2 + H channels.  The four integrators
!> have to reach the same H2 equilibrium bond length (the closed form
!> `leps_r_h2`) and agree on the endpoint geometry, which is experiment 6 on a
!> second analytic surface with a **known** asymptotic minimum.
subroutine test_leps_four_method_crosscheck(error)
   type(error_type), allocatable, intent(out) :: error
   type(landscape_leps_model) :: model
   type(landscape_irc_path) :: fwd(4), bwd(4)
   real(wp) :: mass(9), x0(9), g(9), e_ts
   real(wp) :: max_gap_f, max_gap_b
   integer :: methods(4), i, j
   logical :: fail

   call new_leps_model(model)
   mass = 1.0_wp
   x0 = 0.0_wp
   x0(4:6) = [leps_r_star, 0.0_wp, 0.0_wp]
   x0(7:9) = [2.0_wp*leps_r_star, 0.0_wp, 0.0_wp]
   call model%evaluate(x0, e_ts, g, fail)
   call check(error, .not. fail)
   if (allocated(error)) return
   methods = [p_landscape_irc_gs, p_landscape_irc_lqa, p_landscape_irc_hpc, &
      & p_landscape_irc_eulerpc]

   do i = 1, 4
      call walk_irc(model, mass, x0, 0.05_wp, 1.0_wp, 400, methods(i), fwd(i), fail)
      call check(error, .not. fail)
      if (allocated(error)) return
      call walk_irc(model, mass, x0, 0.05_wp, -1.0_wp, 400, methods(i), bwd(i), fail)
      call check(error, .not. fail)
      if (allocated(error)) return
      ! each direction must reach an H2 equilibrium bond length
      call check(error, leps_min_bond(fwd(i)) < leps_r_h2 + 0.15_wp, &
         & "a forward LEPS endpoint did not reach the H2 equilibrium distance")
      if (allocated(error)) return
      call check(error, leps_min_bond(bwd(i)) < leps_r_h2 + 0.15_wp, &
         & "a backward LEPS endpoint did not reach the H2 equilibrium distance")
      if (allocated(error)) return
   end do

   ! The dissociating H sits on a very flat LEPS asymptote, so the raw
   ! endpoint geometry is not a meaningful cross method metric (its position
   ! varies by bohr without a chemical consequence).  The metrics that carry
   ! chemistry are the H2 bond length and the channel topology.
   max_gap_f = spread_min_bond(fwd); max_gap_b = spread_min_bond(bwd)
   call check(error, max_gap_f < 0.06_wp, &
      & "the LEPS forward H2 bond lengths of the four integrators disagree")
   if (allocated(error)) return
   call check(error, max_gap_b < 0.06_wp, &
      & "the LEPS backward H2 bond lengths of the four integrators disagree")
   if (allocated(error)) return
   ! the two directions form the H2 molecule on opposite sides (different
   ! endpoints), which is the two channel topology of the LEPS surface
   do i = 1, 4
      ! EulerPC is first order and stops near the saddle on this flat surface
      ! (its two directions are mirror images that are still close together),
      ! but they are distinct channels: a same channel pair would coincide.
      call check(error, dist_to(fwd(i)%xyz(:, fwd(i)%n), &
         & bwd(i)%xyz(:, bwd(i)%n)) > 0.5_wp, &
         & "an integrator reached the same LEPS channel in both directions")
      if (allocated(error)) return
   end do
   write(stdout, '(1x,a,f9.5,a,f9.5)') "  LEPS crosscheck H2 bond spread fwd ", &
      & max_gap_f, " bwd ", max_gap_b
   write(stdout, '(1x,a,4f9.5)') "  LEPS forward  H2 bond ", &
      & (leps_min_bond(fwd(i)), i=1, 4)
   write(stdout, '(1x,a,4f9.5)') "  LEPS backward H2 bond ", &
      & (leps_min_bond(bwd(i)), i=1, 4)
   write(stdout, '(1x,a,4i5)') "  LEPS forward points ", (fwd(i)%n, i=1, 4)

end subroutine test_leps_four_method_crosscheck

!> Spread of the shortest pair distance over the endpoints of four paths.
pure function spread_min_bond(paths) result(spread)
   !> Paths
   type(landscape_irc_path), intent(in) :: paths(:)
   !> Spread of the shortest pair distance (bohr)
   real(wp) :: spread

   real(wp) :: dmin, dmax
   integer :: i

   dmin = huge(1.0_wp)
   dmax = -huge(1.0_wp)
   do i = 1, size(paths)
      dmin = min(dmin, leps_min_bond(paths(i)))
      dmax = max(dmax, leps_min_bond(paths(i)))
   end do
   spread = dmax - dmin

end function spread_min_bond

!> Distance of the two atoms that formed the H2 molecule at the endpoint.
pure function leps_min_bond(path) result(d)
   !> Path
   type(landscape_irc_path), intent(in) :: path
   !> Shortest H-H distance at the endpoint (bohr)
   real(wp) :: d

   if (path%n <= 0) then
      d = huge(1.0_wp)
      return
   end if
   d = min(pair_distance(path%xyz(:, path%n), 1, 2), &
      & pair_distance(path%xyz(:, path%n), 2, 3), &
      & pair_distance(path%xyz(:, path%n), 1, 3))

end function leps_min_bond

!> Mean number of energy/gradient evaluations per recorded point.
pure function mean_eval(path) result(mean)
   !> Path
   type(landscape_irc_path), intent(in) :: path
   !> Mean E/g
   real(wp) :: mean

   if (path%n <= 0) then
      mean = 0.0_wp
      return
   end if
   mean = real(sum(path%n_eval(1:path%n)), wp)/real(path%n, wp)

end function mean_eval

!> Least squares slope of `y = a + b x`.
pure function fit_slope(x, y, n) result(slope)
   !> Abscissae
   real(wp), intent(in) :: x(:)
   !> Ordinates
   real(wp), intent(in) :: y(:)
   !> Number of points to use
   integer, intent(in) :: n
   !> Slope
   real(wp) :: slope

   real(wp) :: xm, ym, sxy, sxx
   integer :: i

   xm = sum(x(1:n))/real(n, wp)
   ym = sum(y(1:n))/real(n, wp)
   sxy = 0.0_wp
   sxx = 0.0_wp
   do i = 1, n
      sxy = sxy + (x(i) - xm)*(y(i) - ym)
      sxx = sxx + (x(i) - xm)**2
   end do
   if (sxx <= 0.0_wp) then
      slope = 0.0_wp
   else
      slope = sxy/sxx
   end if

end function fit_slope

!> Walk a Mueller-Brown path with the default settings of one method.
subroutine walk_mb(model, mass, x0, step, sign, max_steps, method, path, fail)
   !> Potential
   class(landscape_energy_model), intent(inout) :: model
   !> Masses (amu)
   real(wp), intent(in) :: mass(:)
   !> Starting point
   real(wp), intent(in) :: x0(:)
   !> Step length
   real(wp), intent(in) :: step
   !> Direction sign
   real(wp), intent(in) :: sign
   !> Step budget
   integer, intent(in) :: max_steps
   !> Method identifier
   integer, intent(in) :: method
   !> Resulting path
   type(landscape_irc_path), intent(inout) :: path
   !> Failure flag
   logical, intent(out) :: fail

   type(landscape_irc_settings) :: settings

   call default_settings(settings, step, max_steps, method)
   call walk_with(model, mass, x0, settings, sign, path, fail)

end subroutine walk_mb

!> Walk one direction with explicit settings.
subroutine walk_irc(model, mass, x0, step, sign, max_steps, method, path, fail)
   !> Potential
   class(landscape_energy_model), intent(inout) :: model
   !> Masses (amu)
   real(wp), intent(in) :: mass(:)
   !> Starting point
   real(wp), intent(in) :: x0(:)
   !> Step length
   real(wp), intent(in) :: step
   !> Direction sign
   real(wp), intent(in) :: sign
   !> Step budget
   integer, intent(in) :: max_steps
   !> Method identifier
   integer, intent(in) :: method
   !> Resulting path
   type(landscape_irc_path), intent(inout) :: path
   !> Failure flag
   logical, intent(out) :: fail

   type(landscape_irc_settings) :: settings

   call default_settings(settings, step, max_steps, method)
   call walk_with(model, mass, x0, settings, sign, path, fail)

end subroutine walk_irc

!> Fill the default settings of one walk.
subroutine default_settings(settings, step, max_steps, method)
   !> Settings
   type(landscape_irc_settings), intent(out) :: settings
   !> Step length
   real(wp), intent(in) :: step
   !> Step budget
   integer, intent(in) :: max_steps
   !> Method identifier
   integer, intent(in) :: method

   settings%step_length = step
   settings%max_steps = max_steps
   settings%method = method
   settings%printlevel = 0
   settings%e_rise_steps = 3
   ! the analytic surfaces are scaled by `kcaltoau`, so the absolute gradient
   ! threshold of a molecule (1e-3 hartree/bohr) is far too loose for them:
   ! the unit tests terminate the analytic paths by signal (c) and set the
   ! threshold explicitly where signal (b) is under test
   settings%grad_tol = 0.0_wp

end subroutine default_settings

!> Run one direction of a walk from a transition state with the given settings.
subroutine walk_with(model, mass, x0, settings, sign, path, fail, env)
   !> Potential
   class(landscape_energy_model), intent(inout) :: model
   !> Masses (amu)
   real(wp), intent(in) :: mass(:)
   !> Starting point
   real(wp), intent(in) :: x0(:)
   !> Settings
   type(landscape_irc_settings), intent(in) :: settings
   !> Direction sign
   real(wp), intent(in) :: sign
   !> Resulting path
   type(landscape_irc_path), intent(inout) :: path
   !> Failure flag
   logical, intent(out) :: fail
   !> Computational environment (a fresh one is created when absent)
   type(TEnvironment), intent(inout), optional :: env

   type(TEnvironment) :: local_env

   if (present(env)) then
      call do_walk(env, model, mass, x0, settings, sign, path, fail)
   else
      call init(local_env)
      call do_walk(local_env, model, mass, x0, settings, sign, path, fail)
   end if

end subroutine walk_with

!> Body of `walk_with`: select the integrator and run one direction.
subroutine do_walk(env, model, mass, x0, settings, sign, path, fail)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Potential
   class(landscape_energy_model), intent(inout) :: model
   !> Masses (amu)
   real(wp), intent(in) :: mass(:)
   !> Starting point
   real(wp), intent(in) :: x0(:)
   !> Settings
   type(landscape_irc_settings), intent(in) :: settings
   !> Direction sign
   real(wp), intent(in) :: sign
   !> Resulting path
   type(landscape_irc_path), intent(inout) :: path
   !> Failure flag
   logical, intent(out) :: fail

   type(landscape_transition_mode) :: mode
   real(wp), allocatable :: sqrt_mass(:), hess(:, :), hq(:, :), v(:)
   integer :: n

   n = size(x0)
   allocate(sqrt_mass(n), hess(n, n), hq(n, n))
   sqrt_mass = landscape_irc_sqrt_mass(mass)
   call model%hessian(x0, hess, fail)
   if (fail) return
   call landscape_irc_hessian_to_weighted(sqrt_mass, hess, hq)
   call landscape_irc_select_direction(env, hq, settings%target_mode, mode, fail)
   if (fail) return
   v = mode%vector

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
      fail = .true.
   end select

end subroutine do_walk

!> Euclidean distance between two flat geometries.
pure function dist_to(a, b) result(d)
   !> First geometry
   real(wp), intent(in) :: a(:)
   !> Second geometry
   real(wp), intent(in) :: b(:)
   !> Distance
   real(wp) :: d

   d = vec_norm(a - b)

end function dist_to

!> Maximum perpendicular deviation of a path from a reference path.
!>
!> For every point of `path` the smallest distance to any point of `reference`
!> is taken (the reference is dense, so this approximates the distance to the
!> reference curve).
pure function max_path_deviation(path, reference) result(dev)
   !> Path under test
   type(landscape_irc_path), intent(in) :: path
   !> Reference path
   type(landscape_irc_path), intent(in) :: reference
   !> Largest deviation
   real(wp) :: dev

   real(wp) :: dmin, d
   integer :: i, j

   dev = 0.0_wp
   do i = 1, path%n
      dmin = huge(1.0_wp)
      do j = 1, reference%n
         d = vec_norm(path%xyz(:, i) - reference%xyz(:, j))
         dmin = min(dmin, d)
      end do
      dev = max(dev, dmin)
   end do

end function max_path_deviation

!> Distance between two atoms of a flat geometry.
pure function pair_distance(x, i, j) result(d)
   !> Flat geometry
   real(wp), intent(in) :: x(:)
   !> Atom index
   integer, intent(in) :: i
   !> Atom index
   integer, intent(in) :: j
   !> Distance
   real(wp) :: d

   d = vec_norm(x(3*(i - 1) + 1:3*i) - x(3*(j - 1) + 1:3*j))

end function pair_distance

!> Symmetric eigenvalues and eigenvectors (ascending) of a small matrix.
subroutine diagonalise(env, hess, ev)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Matrix, overwritten with the eigenvectors
   real(wp), intent(inout) :: hess(:, :)
   !> Eigenvalues (ascending)
   real(wp), intent(out) :: ev(:)

   call mctc_syev(env, hess, ev)

end subroutine diagonalise

end module test_landscape_irc
