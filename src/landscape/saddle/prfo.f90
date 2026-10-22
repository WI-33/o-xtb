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

!> Partitioned rational function optimisation (P-RFO) for transition state
!> searches, and the macro driver that turns the step into a search
!> (PLAN sections 6.3 and 9.5).
!>
!> The step is the **double shift** Newton step of Baker (J. Comput. Chem. 7,
!> 385 (1986); ORCA 6.1 section 4.3.1):
!>
!> ```text
!> E    = E0 + g^T dq + 1/2 dq^T H dq                        (4.16)
!> dq_k = -(F_k V_k)/(b_k - lambda_p),  F_i = V_i^T g         (4.18a, 4.19)
!> lambda_p = 1/2 b_k + 1/2 sqrt(b_k^2 + 4 F_k^2)             (4.18a)
!> sum_{i/=k} F_i^2/(lambda_n - b_i) = lambda_n,  lambda_n < 0 (4.18b)
!> dq_i = -(F_i V_i)/(b_i - lambda_n),  i /= k               (4.19)
!> dq   = sum_j dq_j                                          (4.20)
!> ```
!>
!> The two shifts split the spectrum: the followed mode `k` is *maximised*
!> (`lambda_p` gives an ascent along `V_k`, so the search walks uphill towards
!> the saddle) while the orthogonal complement is *minimised* (`lambda_n < 0`
!> keeps the rest of the step a descent step).  The eigenvalues `b_i` and
!> vectors `V_i` come from a symmetric diagonalisation in **ascending** order
!> (`xtb_mctc_lapack::mctc_syev`); the mode `k` is identified by the largest
!> overlap with the transition vector of `xtb_landscape_saddle_eigenfollow`.
!>
!> Safeguards (PLAN section 6.3):
!>
!>  * More than one negative eigenvalue: the downhill modes are lifted,
!>    `b_i -> max(b_i, hlow)`, so the only direction with negative curvature is
!>    the one that is deliberately followed.  The lift is reported through
!>    `lifted`.
!>  * The `lambda_n` secular equation is solved by bisection on a strictly
!>    decreasing function; it has exactly one root in `(-inf, min_i b_i)` and
!>    that root is negative, which is asserted.
!>  * `evals_eps` guards every denominator, and the step is clipped to the
!>    trust radius.
!>  * The Cartesian step is projected onto the complement of the rigid body
!>    directions by the caller, so the trust radius is spent on internal
!>    motion only.
!>
!> Units: Hessian in hartree/bohr^2, gradient in hartree/bohr, displacements
!> and trust radii in bohr, energies in hartree.
module xtb_landscape_saddle_prfo
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment, only : TEnvironment
   use xtb_landscape_settings, only : landscape_ts_settings, p_landscape_ts_prfo, &
      & p_landscape_ts_bprfo, p_landscape_hess_update_none, &
      & landscape_convergence_measures, landscape_lowercase, landscape_real_is_finite
   use xtb_landscape_results, only : landscape_result, landscape_step, &
      & landscape_result_add_step, landscape_result_trim, landscape_report_add, &
      & landscape_print_step
   use xtb_landscape_coords, only : vec_norm, vec_maxabs
   use xtb_landscape_constraints, only : landscape_project_fixed_atoms
   use xtb_landscape_engine, only : landscape_energy_model
   use xtb_landscape_hessian_update, only : update_by_enum
   use xtb_landscape_saddle_eigenfollow, only : landscape_transition_mode, &
      & landscape_spectrum, landscape_select_transition_mode, &
      & landscape_follow_transition_mode, landscape_rigid_modes, &
      & landscape_project_rigid, landscape_project_hessian_rigid, &
      & landscape_curvature_noise_tol, landscape_rfo_downhill_shift, &
      & landscape_transition_mode_report, landscape_transition_mode_set_phase
   use xtb_landscape_saddle_bprfo, only : landscape_bprfo_step
   implicit none(type, external)
   private

   public :: landscape_prfo_step, landscape_prfo_step_spectrum
   public :: landscape_saddle_optimize

   !> Number of consecutive failed steps before the search is terminated
   integer, parameter :: stagnate_patience = 8

   !> Model energy change below which the initial step counts as an ascent
   !> (hartree, negative = an ascent).
   !>
   !> The threshold only has to separate a genuine ascent from the numerically
   !> zero change of a stationary point: measured model changes of the two
   !> regimes are `-4.83e-2` hartree (Mueller-Brown, a climbing start) and
   !> `O(1e-20)` (a converged minimum), so `1e-12` leaves eight orders of margin
   !> on both sides.
   real(wp), parameter :: escape_ascent_tol = 1.0e-12_wp

   !> Length of the escape walk out of a minimum (bohr).
   !>
   !> The push loop of `escape_from_minimum` is bounded by this distance, so a
   !> mode whose curvature never turns negative (a very soft coordinate, or a
   !> molecule that falls apart under the displacement) cannot walk away
   !> indefinitely.  Three bohr is the order of magnitude on which the bend of
   !> GFN2 HCN turns negative (1.4 bohr measured), with a safety factor of two.
   real(wp), parameter :: escape_max_walk = 3.0_wp

contains

!> P-RFO step from a Hessian and a gradient.
!>
!> Thin wrapper that diagonalises the (symmetrised) Hessian and delegates to
!> `landscape_prfo_step_spectrum`.  The optional outputs expose the shifts and
!> the curvature lift for reporting and tests.
subroutine landscape_prfo_step(env, hess, gradient, trans, trust, hlow, step, &
      & predicted, fail, evals_eps, mu_margin, max_bisect_it, lambda_p, &
      & lambda_n, lifted, lifted_evals, warn)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Hessian (hartree/bohr^2)
   real(wp), intent(in) :: hess(:, :)
   !> Gradient (hartree/bohr)
   real(wp), intent(in) :: gradient(:)
   !> Transition vector (unit, Cartesian)
   real(wp), intent(in) :: trans(:)
   !> Trust radius (bohr)
   real(wp), intent(in) :: trust
   !> Curvature floor of the downhill modes (hartree/bohr^2)
   real(wp), intent(in) :: hlow
   !> P-RFO step (bohr)
   real(wp), intent(out) :: step(:)
   !> Predicted energy change of the quadratic model (hartree)
   real(wp), intent(out) :: predicted
   !> Failure flag
   logical, intent(out) :: fail
   !> Relative denominator guard, defaults to 1e-10
   real(wp), intent(in), optional :: evals_eps
   !> Margin of the `lambda_n` bracket, defaults to 1e-8
   real(wp), intent(in), optional :: mu_margin
   !> Maximum number of bisections, defaults to 60
   integer, intent(in), optional :: max_bisect_it
   !> Shift of the followed mode
   real(wp), intent(out), optional :: lambda_p
   !> Shift of the downhill subspace
   real(wp), intent(out), optional :: lambda_n
   !> True if at least one downhill eigenvalue was lifted to `hlow`
   logical, intent(out), optional :: lifted
   !> Original eigenvalues that were lifted to `hlow` (hartree/bohr^2)
   real(wp), allocatable, intent(out), optional :: lifted_evals(:)
   !> Emit the curvature lift warning (defaults to `.true.`); the macro driver
   !> disables all but the first one so a repeated lift does not flood the log
   logical, intent(in), optional :: warn

   real(wp), allocatable :: evals(:), evecs(:, :)

   fail = .false.
   call landscape_spectrum(env, hess, evals, evecs, fail)
   if (fail) return
   call landscape_prfo_step_spectrum(env, evals, evecs, gradient, trans, trust, &
      & hlow, step, predicted, fail, evals_eps=evals_eps, mu_margin=mu_margin, &
      & max_bisect_it=max_bisect_it, lambda_p=lambda_p, lambda_n=lambda_n, &
      & lifted=lifted, lifted_evals=lifted_evals, warn=warn)

end subroutine landscape_prfo_step

!> P-RFO step from an explicit Hessian spectrum.
!>
!> `evals` must be ascending and `evecs` must hold the matching orthonormal
!> eigenvectors as columns.  Keeping the spectrum as an argument lets the
!> macro driver reuse the diagonalisation it already performed for mode
!> following and lets tests drive the kernel with synthetic spectra.
subroutine landscape_prfo_step_spectrum(env, evals, evecs, gradient, trans, &
      & trust, hlow, step, predicted, fail, evals_eps, mu_margin, &
      & max_bisect_it, lambda_p, lambda_n, lifted, lifted_evals, warn)
   character(len=*), parameter :: source = 'landscape_prfo_step'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Ascending eigenvalues (hartree/bohr^2)
   real(wp), intent(in) :: evals(:)
   !> Eigenvectors as columns
   real(wp), intent(in) :: evecs(:, :)
   !> Gradient (hartree/bohr)
   real(wp), intent(in) :: gradient(:)
   !> Transition vector (unit, Cartesian)
   real(wp), intent(in) :: trans(:)
   !> Trust radius (bohr)
   real(wp), intent(in) :: trust
   !> Curvature floor of the downhill modes (hartree/bohr^2)
   real(wp), intent(in) :: hlow
   !> P-RFO step (bohr)
   real(wp), intent(out) :: step(:)
   !> Predicted energy change of the quadratic model (hartree)
   real(wp), intent(out) :: predicted
   !> Failure flag
   logical, intent(out) :: fail
   !> Relative denominator guard, defaults to 1e-10
   real(wp), intent(in), optional :: evals_eps
   !> Margin of the `lambda_n` bracket, defaults to 1e-8
   real(wp), intent(in), optional :: mu_margin
   !> Maximum number of bisections, defaults to 60
   integer, intent(in), optional :: max_bisect_it
   !> Shift of the followed mode
   real(wp), intent(out), optional :: lambda_p
   !> Shift of the downhill subspace
   real(wp), intent(out), optional :: lambda_n
   !> True if at least one downhill eigenvalue was lifted to `hlow`
   logical, intent(out), optional :: lifted
   !> Original eigenvalues that were lifted to `hlow` (hartree/bohr^2)
   real(wp), allocatable, intent(out), optional :: lifted_evals(:)
   !> Emit the curvature lift warning, defaults to `.true.`
   logical, intent(in), optional :: warn

   real(wp), allocatable :: f(:), b(:), coeff(:), lifted_list(:)
   real(wp) :: eps_guard, margin, lam_p, lam_n, denom, neg_tol
   real(wp) :: nrm, scale, shs
   character(len=160) :: msg
   logical :: was_lifted, solved, emit_warning
   integer :: n, i, k, nneg, n_lifted, max_bisect

   fail = .false.
   step = 0.0_wp
   predicted = 0.0_wp
   if (present(lambda_p)) lambda_p = 0.0_wp
   if (present(lambda_n)) lambda_n = 0.0_wp
   if (present(lifted)) lifted = .false.

   eps_guard = 1.0e-10_wp
   if (present(evals_eps)) eps_guard = evals_eps
   margin = 1.0e-8_wp
   if (present(mu_margin)) margin = mu_margin
   max_bisect = 60
   if (present(max_bisect_it)) max_bisect = max(1, max_bisect_it)

   n = size(evals)
   if (n == 0 .or. size(evecs, 1) /= n .or. size(gradient) /= n) then
      call env%error("the P-RFO step requires a matching Hessian and gradient", &
         & source)
      fail = .true.
      return
   end if
   if (trust <= 0.0_wp) then
      call env%error("the P-RFO trust radius must be positive", source)
      fail = .true.
      return
   end if

   allocate(f(n), b(n), coeff(n))
   do i = 1, n
      f(i) = dot_product(evecs(:, i), gradient)
   end do
   b = evals

   ! the followed mode is the eigenvector with the largest overlap with the
   ! tracked transition vector
   k = 1
   if (size(trans) == n) then
      k = maxloc(abs(matmul(transpose(evecs), trans)), dim=1)
   end if

   ! A negative eigenvalue in the downhill subspace makes the secular equation
   ! of Eq. 4.18b unsolvable below zero (its residual is already positive at
   ! `lambda = 0`), so the downhill modes are lifted to the `hlow` floor.  This
   ! covers both the common case of a genuine negative eigenvalue in the
   ! downhill subspace and the case of a single negative mode that was not
   ! selected as the followed one (for example `--lts-init-mode force`).
   !
   ! Only eigenvalues below `-neg_tol` gate the lift: the rigid body
   ! directions of a numerical Hessian whose spectrum was not projected are
   ! lifted by round off (~eps*||H||) and must not be mistaken for negative
   ! curvature (H4 of `29-m4-independent-review.md`).  Once a significant
   ! negative eigenvalue exists the whole downhill subspace is raised to the
   ! `hlow` floor, exactly as before, so a soft or near-zero downhill mode
   ! cannot blow up the step.
   neg_tol = landscape_curvature_noise_tol(evals)
   nneg = 0
   do i = 1, n
      if (i == k) cycle
      if (evals(i) < -neg_tol) nneg = nneg + 1
   end do
   was_lifted = .false.
   n_lifted = 0
   allocate(lifted_list(n))
   if (nneg > 0) then
      do i = 1, n
         if (i == k) cycle
         if (b(i) < hlow) then
            n_lifted = n_lifted + 1
            lifted_list(n_lifted) = b(i)
            b(i) = hlow
         end if
      end do
      was_lifted = n_lifted > 0
   end if
   emit_warning = .true.
   if (present(warn)) emit_warning = warn
   if (was_lifted .and. emit_warning) then
      write(msg, '(a,i0,a,es12.4,a)') "the downhill subspace has ", n_lifted, &
         & " negative curvature eigenvalue(s) (most negative ", &
         & minval(lifted_list(1:n_lifted)), &
         & "); the downhill modes were lifted to hlow"
      call env%warning(trim(msg), source)
   end if
   if (present(lifted_evals)) then
      if (n_lifted > 0) then
         lifted_evals = lifted_list(1:n_lifted)
      else
         allocate(lifted_evals(0))
      end if
   end if

   ! lambda_p: ascend along the followed mode (Eq. 4.18a)
   lam_p = 0.5_wp*b(k) + 0.5_wp*sqrt(b(k)**2 + 4.0_wp*f(k)**2)

   denom = b(k) - lam_p
   if (abs(denom) <= eps_guard) then
      coeff(k) = 0.0_wp
   else
      coeff(k) = -f(k)/denom
   end if

   ! lambda_n: the downhill shift (Eq. 4.18b)
   call landscape_rfo_downhill_shift(env, f, b, k, margin, max_bisect, lam_n, &
      & solved, fail, source)
   if (fail) then
      call env%error("the P-RFO secular equation could not be solved", source)
      return
   end if
   if (solved .and. lam_n >= 0.0_wp) then
      ! Cannot happen because the secular function is negative at zero for a
      ! non-empty downhill subspace; asserted so that a refactor cannot
      ! silently lose the sign condition of Eq. 4.18b.
      call env%error("the P-RFO downhill shift is not negative", source)
      fail = .true.
      return
   end if
   do i = 1, n
      if (i == k) cycle
      denom = b(i) - lam_n
      if (abs(denom) <= eps_guard) then
         coeff(i) = 0.0_wp
      else
         coeff(i) = -f(i)/denom
      end if
   end do

   do i = 1, n
      step = step + coeff(i)*evecs(:, i)
   end do

   ! trust region: the two shifts already regulate the step, the radius is the
   ! hard safeguard of PLAN section 6.3
   nrm = vec_norm(step)
   if (nrm > trust) then
      scale = trust/nrm
      step = scale*step
   end if

   ! predicted change of the (lifted) quadratic model along the final step
   shs = 0.0_wp
   do i = 1, n
      shs = shs + b(i)*dot_product(evecs(:, i), step)**2
   end do
   predicted = -(dot_product(gradient, step) + 0.5_wp*shs)

   if (present(lambda_p)) lambda_p = lam_p
   if (present(lambda_n)) lambda_n = lam_n
   if (present(lifted)) lifted = was_lifted

end subroutine landscape_prfo_step_spectrum

!> Macro driver of a P-RFO / B-P-RFO transition state search.
!>
!> The loop mirrors `xtb_landscape_optimize::landscape_optimize_model`, but the
!> target is a stationary point of index one rather than a minimum:
!>
!>  * the transition vector is selected from the initial Hessian and then
!>    followed across the macro iterations (`landscape_follow_transition_mode`);
!>  * every step is projected onto the complement of the rigid body directions
!>    before it is taken, which closes the M1 debt of
!>    `xtb_landscape_rfo::landscape_rfo_step` (PLAN section 9.5);
!>  * the step is always taken when the trial evaluation succeeds, because for
!>    a saddle the energy may rise; the trust radius, not an energy decrease,
!>    regulates the magnitude.  A failed evaluation shrinks the radius and
!>    repeats the geometry.
!>
!> `x_final` returns the converged geometry and `traj` the recorded
!> trajectory (`(ndim, nframes)`, first frame is the start).
subroutine landscape_saddle_optimize(env, objective, settings, x0, result, &
      & x_final, traj)
   character(len=*), parameter :: source = 'landscape_saddle_optimize'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Objective function
   class(landscape_energy_model), intent(inout) :: objective
   !> Transition state settings
   type(landscape_ts_settings), intent(in) :: settings
   !> Starting coordinates (bohr)
   real(wp), intent(in) :: x0(:)
   !> Search result
   type(landscape_result), intent(inout) :: result
   !> Converged coordinates (bohr)
   real(wp), allocatable, intent(out), optional :: x_final(:)
   !> Recorded trajectory (bohr)
   real(wp), allocatable, intent(out), optional :: traj(:, :)

   type(landscape_transition_mode) :: mode
   type(landscape_step) :: record
   real(wp), allocatable :: x(:), g(:), step(:), x_trial(:), g_trial(:)
   real(wp), allocatable :: hmat(:, :), hnew(:, :), frames(:, :)
   real(wp), allocatable :: evals(:), evecs(:, :)
   real(wp), allocatable :: rigid(:, :), lifted_evals(:)
   real(wp) :: energy, e_trial, de, trust, rho, predicted, lam_p, lam_n
   real(wp) :: g_max, g_rms, d_max, d_rms
   logical :: fail, ok, accepted, lifted, check_displ, warned_lift
   integer :: istep, nreject, n_evals, ndim, nframes, n_rigid, method

   if (size(x0) == 0) return
   call reset_result(result)
   ndim = size(x0)
   method = settings%method

   allocate(x(ndim), g(ndim), step(ndim), x_trial(ndim), g_trial(ndim))
   allocate(hmat(ndim, ndim), hnew(ndim, ndim), source=0.0_wp)
   x = x0
   g = 0.0_wp
   step = 0.0_wp
   e_trial = 0.0_wp
   g_trial = 0.0_wp
   de = 0.0_wp
   nreject = 0
   n_evals = 0
   warned_lift = .false.

   if (.not. all(landscape_real_is_finite(x))) then
      result%message = "initial coordinates are not finite"
      call finish()
      return
   end if

   call objective%evaluate(x, energy, g, fail)
   n_evals = n_evals + 1
   if (fail .or. .not. landscape_real_is_finite(energy) .or. &
         & .not. all(landscape_real_is_finite(g))) then
      result%message = "initial energy or gradient is not finite"
      call finish()
      return
   end if

   call objective%hessian(x, hmat, fail)
   if (fail) then
      result%message = "initial Hessian unavailable"
      call finish()
      return
   end if
   hmat = 0.5_wp*(hmat + transpose(hmat))

   ! The rigid body directions are removed from the Hessian *spectrum* as well
   ! as from the step: the numerical Hessian lifts its null directions by round
   ! off, and those tiny eigenvalues must not be counted as negative curvature
   ! (H4 of `29-m4-independent-review.md`).  The basis is rebuilt from the
   ! current geometry at the top of every macro iteration.
   n_rigid = 0
   if (settings%remove_rigid) then
      call landscape_rigid_modes(x, rigid, n_rigid)
      call landscape_project_hessian_rigid(hmat, rigid, n_rigid)
   end if

   call landscape_spectrum(env, hmat, evals, evecs, fail)
   if (fail) then
      result%message = "initial Hessian could not be diagonalised"
      call finish()
      return
   end if
   call select_initial_mode(env, settings, evals, evecs, g, mode, fail)
   if (fail) then
      result%message = "initial transition vector unavailable"
      call finish()
      return
   end if
   call report_initial_mode_diagnostics()

   trust = min(settings%trust, settings%trust_max)
   if (trust <= 0.0_wp) then
      call env%error("the transition state trust radius must be positive", source)
      result%message = "invalid trust radius"
      call finish()
      return
   end if

   ! Escape from a minimum: P-RFO can only start uphill along the followed mode
   ! when the initial model step ascends.  On a stationary minimum the ascent
   ! vanishes with the gradient and the search re-minimises, so the start has to
   ! be displaced along the lowest internal mode (MEDIUM-1 of
   ! `31-m4-review-round2.md`).
   if (.not. initial_step_ascends()) call escape_from_minimum()
   if (fail) then
      call finish()
      return
   end if

   allocate(frames(ndim, 32))
   nframes = 1
   frames(:, 1) = x

   do istep = 1, settings%max_steps
      step = 0.0_wp
      predicted = 0.0_wp
      accepted = .false.
      lifted = .false.

      ! rebuild the rigid body basis from the current geometry before the
      ! spectrum is refreshed, so the Hessian of this step is projected too
      if (settings%remove_rigid) then
         call landscape_rigid_modes(x, rigid, n_rigid)
      end if

      if (settings%recalc >= 1 .and. istep > 1 .and. &
            & mod(istep - 1, settings%recalc) == 0) then
         call objective%hessian(x, hnew, fail)
         if (.not. fail) then
            hmat = 0.5_wp*(hnew + transpose(hnew))
            if (n_rigid > 0) call landscape_project_hessian_rigid(hmat, rigid, n_rigid)
            call landscape_spectrum(env, hmat, evals, evecs, fail)
            if (.not. fail) then
               if (settings%printlevel >= 2) then
                  call landscape_follow_transition_mode(env, evals, evecs, mode, &
                     & settings%mode_overlap_min, fail, unit=env%unit)
               else
                  call landscape_follow_transition_mode(env, evals, evecs, mode, &
                     & settings%mode_overlap_min, fail)
               end if
            end if
         end if
         if (fail) then
            call env%warning("Hessian recalculation failed; keeping the "//&
               & "reduced Hessian", source)
            fail = .false.
         end if
      end if

      select case (method)
      case (p_landscape_ts_prfo)
         call landscape_prfo_step(env, hmat, g, mode%vector, trust, &
            & settings%hlow, step, predicted, fail, &
            & evals_eps=settings%evals_eps, mu_margin=settings%mu_margin, &
            & max_bisect_it=settings%max_bisect_it, lambda_p=lam_p, &
            & lambda_n=lam_n, lifted=lifted, lifted_evals=lifted_evals, &
            & warn=.not. warned_lift)
      case (p_landscape_ts_bprfo)
         call landscape_bprfo_step(env, hmat, g, mode%vector, trust, &
            & settings%alpha, settings%hlow, step, predicted, fail, &
            & evals_eps=settings%evals_eps, mu_margin=settings%mu_margin, &
            & max_micro_cycles=settings%max_micro_cycles, &
            & lambda_p=lam_p, lambda_n=lam_n, lifted=lifted, &
            & lifted_evals=lifted_evals, warn=.not. warned_lift)
      case default
         call env%error("unknown transition state method", source)
         result%message = "unknown transition state method"
         call finish()
         return
      end select

      ! The curvature lift warning is emitted once per run (H4 of
      ! `29-m4-independent-review.md`): the per-step record still carries
      ! `lifted=T`, but the log is not flooded by the same message every step.
      if (lifted .and. .not. warned_lift) then
         warned_lift = .true.
         if (settings%printlevel > 0 .and. allocated(lifted_evals) .and. &
               & size(lifted_evals) > 0) then
            write(env%unit, '(1x,a,i0,a,es15.8,a)') "  curvature lift: ", &
               & size(lifted_evals), " downhill eigenvalue(s) lifted to "//&
               & "hlow (most negative ", minval(lifted_evals), " Eh/bohr^2)"
         end if
      end if

      if (fail) then
         trust = max(settings%trust_min, settings%eta_shrink*trust)
         nreject = nreject + 1
         if (nreject >= stagnate_patience) then
            result%message = "no step could be constructed"
            exit
         end if
         call record_step(result, env, settings, istep, energy, 0.0_wp, g, step, &
            & trust, method, .true.)
         cycle
      end if

      if (n_rigid > 0) then
         call landscape_project_rigid(step, rigid, n_rigid)
      end if
      ! `$fix`: the accepted step must not move a fixed atom.  The gradient is
      ! already zero there (the calculator does that), but the P-RFO step is
      ! expanded in the eigenvectors of the Hessian, which mix the frozen and
      ! the free components - the same defect the minimiser had.  Projecting
      ! the step (and not the Hessian) is deliberate: unlike a rigid body mode,
      ! a `$fix` direction is **not** a null direction of the Hessian, so
      ! removing it from the matrix would create three artificial zero
      ! eigenvalues per fixed atom and distort the eigenvalue selection and the
      ! `hlow` lift.  The projections of the rigid directions are unaffected -
      ! both are applied to the same step, one after the other.
      call landscape_project_fixed_atoms(step)

      x_trial = x + step
      call objective%evaluate(x_trial, e_trial, g_trial, fail)
      n_evals = n_evals + 1
      if (fail .or. .not. landscape_real_is_finite(e_trial) .or. &
            & .not. all(landscape_real_is_finite(g_trial))) then
         trust = max(settings%trust_min, settings%eta_shrink*trust)
         nreject = nreject + 1
         if (nreject >= stagnate_patience) then
            result%message = "trial evaluations keep failing"
            exit
         end if
         call record_step(result, env, settings, istep, energy, 0.0_wp, g, step, &
            & trust, method, .true.)
         cycle
      end if

      ! the ratio only regulates the trust radius; for a saddle the step has to
      ! be taken even when the energy rises
      rho = 0.0_wp
      if (predicted > 0.0_wp) rho = (energy - e_trial)/predicted
      if (rho < 0.25_wp) then
         trust = max(settings%trust_min, settings%eta_shrink*trust)
      else if (rho > 0.75_wp) then
         trust = min(settings%trust_max, settings%eta_expand*trust)
      end if

      if (settings%hessian_update /= p_landscape_hess_update_none) then
         call update_by_enum(settings%hessian_update, step, g_trial - g, &
            & hmat, ok)
      end if

      de = e_trial - energy
      x = x_trial
      energy = e_trial
      g = g_trial
      accepted = .true.
      nreject = 0
      nframes = nframes + 1
      call append_frame(frames, nframes, x)

      call record_step(result, env, settings, istep, energy, de, g, step, &
         & trust, method, .false., lam_p, lam_n, lifted, mode%followed)

      check_displ = .true.
      if (settings%is_converged(de, g, step, check_displ=check_displ)) then
         result%converged = .true.
         result%message = "converged"
         exit
      end if
   end do

   if (.not. allocated(result%message)) then
      if (settings%max_steps >= 1) then
         result%message = "maximum number of steps reached"
      else
         result%message = "no step performed"
      end if
   end if

   result%energy = energy
   result%gnorm = vec_norm(g)
   result%gmax = vec_maxabs(g)
   result%n_eval = n_evals
   call landscape_convergence_measures(g, step, g_max, g_rms, d_max, d_rms)
   result%conv%converged = [g_max <= settings%g_max, g_rms <= settings%g_rms, &
      & d_max <= settings%d_max, d_rms <= settings%d_rms]

   ! Report the *final* transition vector at the converged geometry with the
   ! same phase convention the initial selection used, so the printed direction
   ! is well defined and matches the dimer (which also phases its final axis).
   ! M6 (IRC) consumes this direction to choose the forward branch.
   if (allocated(mode%vector)) then
      call landscape_transition_mode_set_phase(mode)
      if (settings%printlevel > 0) then
         write(env%unit, '(1x,a)') "final transition vector"
         call landscape_transition_mode_report(mode, env%unit)
      end if
   end if

   call finish()

contains

   !> Report how the initial transition vector relates to the rigid body space.
   !>
   !> `||R v||` with `R = sum_a b_a b_a^T` is the weight of the selected mode in
   !> the translation/rotation subspace of the current geometry: ~1 for a rigid
   !> body direction and ~0 for an internal mode.  It is printed because a run
   !> that follows a rigid body zero mode (the MEDIUM-1 failure of
   !> `31-m4-review-round2.md`) shows up here before it shows up in the outcome.
   subroutine report_initial_mode_diagnostics()
      real(wp), allocatable :: probe(:)
      integer :: a

      if (settings%printlevel <= 0) return
      if (n_rigid <= 0 .or. .not. allocated(mode%vector)) return
      allocate(probe(ndim), source=0.0_wp)
      do a = 1, n_rigid
         probe = probe + dot_product(rigid(:, a), mode%vector)*rigid(:, a)
      end do
      write(env%unit, '(1x,a,f12.6)') "selected mode, rigid body overlap: ", &
         & vec_norm(probe)
   end subroutine report_initial_mode_diagnostics

   !> True if the quadratic model of the initial geometry predicts an ascent.
   !>
   !> The probe is the first step of the search, evaluated for its model change
   !> only.  `predicted` is the *decrease* of the (lifted) quadratic model, so
   !> `predicted < -escape_ascent_tol` means the step walks uphill - the search
   !> has a direction to ascend along and must not be displaced (measured on
   !> Mueller-Brown started at `(-0.6, 0.4)`: `predicted = -4.83e-2` hartree,
   !> the search climbs to the saddle in six steps).  A descent, or a change
   !> below the tolerance, is what a (near) stationary minimum gives.
   logical function initial_step_ascends()
      real(wp) :: step_probe(ndim), predicted_probe
      logical :: probe_fail

      initial_step_ascends = .false.
      select case (method)
      case (p_landscape_ts_prfo)
         call landscape_prfo_step(env, hmat, g, mode%vector, trust, &
            & settings%hlow, step_probe, predicted_probe, probe_fail, &
            & evals_eps=settings%evals_eps, mu_margin=settings%mu_margin, &
            & max_bisect_it=settings%max_bisect_it, warn=.false.)
      case (p_landscape_ts_bprfo)
         call landscape_bprfo_step(env, hmat, g, mode%vector, trust, &
            & settings%alpha, settings%hlow, step_probe, predicted_probe, &
            & probe_fail, evals_eps=settings%evals_eps, &
            & mu_margin=settings%mu_margin, &
            & max_micro_cycles=settings%max_micro_cycles, warn=.false.)
      case default
         return
      end select
      if (probe_fail) then
         initial_step_ascends = .true.
         return
      end if
      initial_step_ascends = predicted_probe < -escape_ascent_tol
   end function initial_step_ascends

   !> Walk the start geometry uphill along the lowest internal mode on a minimum.
   !>
   !> P-RFO maximises along the followed mode, so it can only start uphill when
   !> the quadratic model of the initial geometry ascends.  On a stationary
   !> minimum the ascent component `-F_k/(b_k - lambda_p)` vanishes like the
   !> gradient and the search descends back into the minimum: measured on GFN2
   !> HCN started on its linear minimum, the corrected mode selection alone
   !> converges back to that minimum in three steps (0 imaginary modes, 0.0486
   !> bohr first step, `-2.03 kcal/mol`).  The caller only enters this walk when
   !> the initial model step is *not* an ascent, and the walk itself only runs
   !> while the spectrum has no negative curvature.  It pushes the geometry
   !> along the lowest internal mode in steps of `settings%escape_step` bohr,
   !> follows the mode across the pushes and stops as soon as the followed mode
   !> has turned into a genuine negative curvature, which is what P-RFO needs to
   !> continue on its own.
   !>
   !> A single push is not enough in general: on the linear HCN minimum the
   !> lowest mode (the HCN bend) carries `+0.0824 Eh/bohr^2` at the minimum and
   !> is still `+0.0960 Eh/bohr^2` after 0.2 bohr, so its curvature only becomes
   !> negative after 1.4 bohr of walking (measured; the walk stops at
   !> `-0.0053 Eh/bohr^2`).  The walk is bounded by `escape_max_walk` bohr (one
   !> Hessian per push, at most `ceiling(escape_max_walk/escape_step)` pushes)
   !> and reports the number of pushes, the total displacement and the final
   !> curvature of the mode.
   !>
   !> Every failure (a rejected evaluation, an unavailable Hessian, a zero
   !> displacement) only warns and leaves the last valid geometry as the start
   !> of the search, so the user keeps control.  `settings%escape_step <= 0`
   !> disables the walk; a walk that runs out of pushes also leaves the last
   !> geometry in place - the search then starts on a point that still has no
   !> negative curvature, which the confirmation reports as a minimum (H2).
   subroutine escape_from_minimum()
      real(wp), allocatable :: dvec(:), x_esc(:), g_esc(:), h_esc(:, :)
      real(wp), allocatable :: evals_esc(:), evecs_esc(:, :)
      real(wp) :: e_esc, tol
      integer :: itry, npush, ntries
      logical :: local_fail

      if (settings%escape_step <= 0.0_wp) return
      if (method /= p_landscape_ts_prfo .and. method /= p_landscape_ts_bprfo) return
      if (.not. allocated(mode%vector)) return

      tol = landscape_curvature_noise_tol(evals)
      if (any(evals < -tol)) return
      ntries = max(1, ceiling(escape_max_walk/settings%escape_step))

      allocate(dvec(ndim), x_esc(ndim), g_esc(ndim), h_esc(ndim, ndim))
      npush = 0
      do itry = 1, ntries
         dvec = settings%escape_step*mode%vector
         if (n_rigid > 0) call landscape_project_rigid(dvec, rigid, n_rigid)
         ! the escape walk is a geometry displacement as well, so it may not
         ! move a fixed atom either; a mode that lives entirely on the fixed
         ! coordinates therefore stops the walk instead of being followed
         call landscape_project_fixed_atoms(dvec)
         if (vec_norm(dvec) <= tiny(1.0_wp)) then
            call env%warning("the lowest mode has no internal component; the "//&
               & "escape walk is stopped", source)
            exit
         end if

         x_esc = x + dvec
         call objective%evaluate(x_esc, e_esc, g_esc, local_fail)
         n_evals = n_evals + 1
         if (local_fail .or. .not. landscape_real_is_finite(e_esc) .or. &
               & .not. all(landscape_real_is_finite(g_esc))) then
            call env%warning("the escape displacement could not be evaluated; "//&
               & "the last valid geometry is kept as the start", source)
            exit
         end if
         call objective%hessian(x_esc, h_esc, local_fail)
         if (local_fail) then
            call env%warning("the escape displacement has no Hessian; the "//&
               & "last valid geometry is kept as the start", source)
            exit
         end if
         h_esc = 0.5_wp*(h_esc + transpose(h_esc))
         n_rigid = 0
         if (settings%remove_rigid) then
            call landscape_rigid_modes(x_esc, rigid, n_rigid)
            call landscape_project_hessian_rigid(h_esc, rigid, n_rigid)
         end if
         call landscape_spectrum(env, h_esc, evals_esc, evecs_esc, local_fail)
         if (local_fail) then
            call env%warning("the escaped Hessian could not be diagonalised; "//&
               & "the last valid geometry is kept as the start", source)
            exit
         end if

         ! commit the push
         x = x_esc
         energy = e_esc
         g = g_esc
         hmat = h_esc
         call move_alloc(evals_esc, evals)
         call move_alloc(evecs_esc, evecs)
         npush = npush + 1
         call select_initial_mode(env, settings, evals, evecs, g, mode, local_fail)
         if (local_fail) then
            result%message = "initial transition vector unavailable"
            fail = .true.
            return
         end if
         mode%followed = .false.
         tol = landscape_curvature_noise_tol(evals)
         if (settings%printlevel >= 2) then
            write(env%unit, '(1x,a,i0,a,f10.6,a,es14.6,a)') &
               & "escape push ", npush, ", total displacement ", &
               & real(npush, wp)*settings%escape_step, " bohr, lowest mode ", &
               & mode%eigenvalue, " Eh/bohr^2"
         end if
         if (any(evals < -tol)) exit
      end do

      if (npush > 0 .and. settings%printlevel > 0) then
         write(env%unit, '(1x,a,i0,a,f10.6,a,es14.6,a)') &
            & "escape from a minimum: ", npush, " push(es), total "//&
            & "displacement ", real(npush, wp)*settings%escape_step, &
            & " bohr, lowest mode ", mode%eigenvalue, " Eh/bohr^2"
      end if
      call report_initial_mode_diagnostics()
   end subroutine escape_from_minimum

   !> Record the final coordinates and trajectory.
   subroutine finish()
      call landscape_result_trim(result)
      if (present(x_final)) x_final = x
      if (present(traj)) then
         if (allocated(frames)) then
            traj = frames(:, 1:nframes)
         else
            allocate(traj(size(x), 1))
            traj(:, 1) = x
         end if
      end if
   end subroutine finish

   !> Append one macro step to the result and print it.
   subroutine record_step(res, env_, settings_, istep, energy_, de_, g_, step_, &
         & trust_, method_, rejected, lam_p_, lam_n_, lifted_, followed_)
      !> Result
      type(landscape_result), intent(inout) :: res
      !> Environment
      type(TEnvironment), intent(inout) :: env_
      !> Settings
      type(landscape_ts_settings), intent(in) :: settings_
      !> Iteration number
      integer, intent(in) :: istep
      !> Energy (hartree)
      real(wp), intent(in) :: energy_
      !> Energy change (hartree)
      real(wp), intent(in) :: de_
      !> Gradient
      real(wp), intent(in) :: g_(:)
      !> Step (bohr)
      real(wp), intent(in) :: step_(:)
      !> Trust radius (bohr)
      real(wp), intent(in) :: trust_
      !> Method identifier
      integer, intent(in) :: method_
      !> True for a rejected step
      logical, intent(in) :: rejected
      !> Shift of the followed mode
      real(wp), intent(in), optional :: lam_p_
      !> Downhill shift
      real(wp), intent(in), optional :: lam_n_
      !> Curvature lift flag
      logical, intent(in), optional :: lifted_
      !> Mode following flag
      logical, intent(in), optional :: followed_

      record%istep = istep
      record%energy = energy_
      record%delta_e = de_
      record%gnorm = vec_norm(g_)
      record%gmax = vec_maxabs(g_)
      record%dispnorm = vec_norm(step_)
      record%disp_max = vec_maxabs(step_)
      record%method = method_
      record%fallback = rejected
      record%trust_radius = trust_
      call landscape_result_add_step(res, record)
      call landscape_report_add(res%conv, record)
      if (settings_%printlevel >= 1) then
         call landscape_print_step(record, env_%unit)
         if (present(lam_p_) .and. present(lam_n_)) then
            write(env_%unit, '(1x,a,es14.6,a,es14.6,a,l1,a,l1)') &
               & "  shifts: lambda_p=", lam_p_, " lambda_n=", lam_n_, &
               & " lifted=", lifted_, " followed=", followed_
         end if
      end if
   end subroutine record_step

   !> Select the initial transition vector from the settings.
   subroutine select_initial_mode(env_, settings_, evals_, evecs_, g_, mode_, fail_)
      !> Environment
      type(TEnvironment), intent(inout) :: env_
      !> Settings
      type(landscape_ts_settings), intent(in) :: settings_
      !> Spectrum
      real(wp), intent(in) :: evals_(:)
      !> Eigenvectors
      real(wp), intent(in) :: evecs_(:, :)
      !> Gradient
      real(wp), intent(in) :: g_(:)
      !> Transition mode
      type(landscape_transition_mode), intent(inout) :: mode_
      !> Failure flag
      logical, intent(out) :: fail_

      real(wp) :: direction(size(g_))
      integer :: imax

      fail_ = .false.
      ! `landscape_ts_init_mode_known` accepts any spelling, so the gate that
      ! turns the accepted name into a transition vector has to as well.
      if (landscape_lowercase(settings_%init_mode) == 'force') then
         imax = maxloc(abs(g_), dim=1)
         direction = 0.0_wp
         direction(imax) = 1.0_wp
         if (allocated(mode_%vector)) deallocate(mode_%vector)
         allocate(mode_%vector(size(g_)))
         mode_%vector = direction
         mode_%index = 0
         mode_%eigenvalue = 0.0_wp
         mode_%overlap = 1.0_wp
         mode_%followed = .false.
         mode_%target = max(1, settings_%target_mode)
         mode_%phase_index = imax
         mode_%phase_value = 1.0_wp
         if (settings_%printlevel > 0) then
            write(env_%unit, '(1x,a,f12.6,a,i0)') &
               & "transition vector from the largest gradient component: ", &
               & g_(imax), " (component ", imax, ")"
         end if
      else if (settings_%printlevel > 0) then
         call landscape_select_transition_mode(env_, evals_, evecs_, &
            & settings_%target_mode, mode_, fail_, unit=env_%unit)
      else
         call landscape_select_transition_mode(env_, evals_, evecs_, &
            & settings_%target_mode, mode_, fail_)
      end if
   end subroutine select_initial_mode

end subroutine landscape_saddle_optimize

!> Reset a result before a new search.
subroutine reset_result(result)
   !> Result
   type(landscape_result), intent(inout) :: result

   result%converged = .false.
   result%n_steps = 0
   result%energy = 0.0_wp
   result%gnorm = 0.0_wp
   result%gmax = 0.0_wp
   result%n_eval = 0
   result%n_cache_miss = 0
   result%n_fallback = 0
   result%n_imag = -1
   result%imag_freq = 0.0_wp
   result%is_saddle = .false.
   if (allocated(result%message)) deallocate(result%message)
   if (allocated(result%steps)) deallocate(result%steps)
   if (allocated(result%conv%history_energy)) deallocate(result%conv%history_energy)
   if (allocated(result%conv%history_gnorm)) deallocate(result%conv%history_gnorm)
   if (allocated(result%conv%history_gmax)) deallocate(result%conv%history_gmax)
   if (allocated(result%conv%history_dispnorm)) deallocate(result%conv%history_dispnorm)
   if (allocated(result%conv%history_dispmax)) deallocate(result%conv%history_dispmax)
   if (allocated(result%conv%history_trust)) deallocate(result%conv%history_trust)
   if (allocated(result%conv%history_method)) deallocate(result%conv%history_method)
   result%conv%n_iter = 0
   result%conv%converged = .false.

end subroutine reset_result

!> Append a geometry to a trajectory buffer, growing it geometrically.
subroutine append_frame(frames, nframes, x)
   !> Trajectory buffer (bohr)
   real(wp), allocatable, intent(inout) :: frames(:, :)
   !> Number of recorded frames
   integer, intent(inout) :: nframes
   !> Geometry to append (bohr)
   real(wp), intent(in) :: x(:)

   real(wp), allocatable :: buffer(:, :)
   integer :: capacity

   if (nframes > size(frames, 2)) then
      capacity = 2*size(frames, 2)
      allocate(buffer(size(x), capacity))
      buffer(:, 1:nframes - 1) = frames(:, 1:nframes - 1)
      call move_alloc(buffer, frames)
   end if
   frames(:, nframes) = x

end subroutine append_frame

end module xtb_landscape_saddle_prfo
