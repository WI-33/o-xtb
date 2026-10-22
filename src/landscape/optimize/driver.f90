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

!> Unified geometry minimisation core of the landscape suite.
!>
!> `landscape_optimize_model` runs the minimisation for an arbitrary
!> `landscape_energy_model`, i.e. for analytic test potentials as well as for
!> the engine backed model built by `xtb_landscape_tasks_opt`.  It contains
!> only the pure mathematical kernel: algorithm dispatch, the macro iteration,
!> the convergence test and the trajectory buffer.  The task level wrapper and
!> the file output live in `tasks/opt.f90` and `types/log.f90` (PLAN section
!> 4.2), so this layer never touches the xtb private types.
!>
!> Algorithm dispatch follows `landscape_opt_settings%algorithm`
!> (sd, cg, lbfgs, bfgs, diis, rfo).  Every macro iteration records a
!> `landscape_step`, the four convergence criteria are tested on the last
!> displacement and the gradient, and a stagnation counter terminates runs that
!> can no longer make progress.  Any step produced by the steepest descent
!> backstop is flagged and counted (`landscape_result%n_fallback`,
!> `landscape_step%fallback`), so a quasi-Newton run that silently degraded to
!> SD is visible to the caller.
!>
!> Units: coordinates in bohr, energies in hartree, gradients in hartree/bohr.
module xtb_landscape_optimize
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment, only : TEnvironment
   use xtb_landscape_settings, only : landscape_value_present, landscape_opt_settings, &
      & p_landscape_opt_sd, p_landscape_opt_cg, p_landscape_opt_lbfgs, &
      & p_landscape_opt_bfgs, p_landscape_opt_diis, p_landscape_opt_rfo, &
      & landscape_opt_algorithm_name, landscape_convergence_measures, &
      & landscape_lowercase, landscape_real_is_finite
   use xtb_landscape_results, only : landscape_result, landscape_step, &
      & landscape_result_add_step, landscape_result_trim, landscape_report_add, &
      & landscape_print_step
   use xtb_landscape_coords, only : vec_norm, vec_maxabs
   use xtb_landscape_constraints, only : landscape_project_frozen
   use xtb_landscape_engine, only : landscape_energy_model
   use xtb_landscape_line_search, only : landscape_line_search
   use xtb_landscape_sd, only : landscape_sd_direction
   use xtb_landscape_cg, only : landscape_cg_type, new_landscape_cg
   use xtb_landscape_lbfgs, only : landscape_lbfgs_type, new_landscape_lbfgs
   use xtb_landscape_bfgs, only : landscape_bfgs_type, new_landscape_bfgs
   use xtb_landscape_diis, only : landscape_diis_type, new_landscape_diis
   use xtb_landscape_rfo, only : landscape_rfo_step, landscape_rfo_update_trust
   use xtb_landscape_hessian_update, only : update_by_enum
   implicit none(type, external)
   private

   public :: landscape_optimize_model

   !> Number of consecutive stagnated iterations before a run is terminated
   integer, parameter :: stagnate_patience = 10
   !> Number of consecutive rejected RFO steps before an SD fallback step
   integer, parameter :: rfo_reject_patience = 3

contains

!> Core minimisation for an arbitrary potential.
!>
!> `x_final` returns the optimised coordinates and `traj` the recorded
!> trajectory (`(ndim, nframes)`, the first frame is the starting geometry).
subroutine landscape_optimize_model(env, objective, settings, x0, result, &
      & x_final, traj, frozen)
   character(len=*), parameter :: source = 'landscape_optimize_model'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Objective function
   class(landscape_energy_model), intent(inout) :: objective
   !> Optimisation settings
   type(landscape_opt_settings), intent(in) :: settings
   !> Starting coordinates (bohr)
   real(wp), intent(in) :: x0(:)
   !> Optimisation result
   type(landscape_result), intent(inout) :: result
   !> Optimised coordinates (bohr)
   real(wp), allocatable, intent(out), optional :: x_final(:)
   !> Recorded trajectory (bohr)
   real(wp), allocatable, intent(out), optional :: traj(:, :)
   !> Mask of the Cartesian components held fixed by `$fix` (length `3N`,
   !> `xtb_landscape_constraints::landscape_frozen_mask`).  Every search
   !> direction is projected onto the subspace this mask leaves free, which is
   !> what makes the RFO step respect `$fix` (see the module documentation of
   !> `xtb_landscape_constraints`).  Omitted, or all-false, means unconstrained
   !> and reproduces the previous behaviour bitwise.
   logical, intent(in), optional :: frozen(:)
   type(landscape_lbfgs_type) :: lbfgs
   type(landscape_cg_type) :: cg
   type(landscape_bfgs_type) :: bfgs
   type(landscape_diis_type) :: diis
   type(landscape_step) :: record
   real(wp), allocatable :: x(:), g(:), d(:), step(:), x_trial(:), g_trial(:)
   real(wp), allocatable :: hmat(:, :), hnew(:, :), hinit(:, :), frames(:, :)
   logical, allocatable :: fix_mask(:)
   real(wp) :: energy, e_trial, de, trust, trust_hi, rho, predicted
   real(wp) :: gmax, grms, dmax, drms
   logical :: fail, ok, ls_fail, steepest, accepted, use_trust, check_displ
   logical :: fallback, fallback_step, constrained
   integer :: istep, nstagnate, nreject, n_evals, ls_evals, ndim, nframes
   integer :: method

   if (size(x0) == 0) return
   call reset_result(result)
   ndim = size(x0)

   ! `$fix`: a second order step must not move a fixed atom even though its
   ! gradient is zero, so every direction is projected onto the subspace the
   ! mask leaves free (see `xtb_landscape_constraints`).  With an absent or
   ! all-false mask `constrained` stays false and nothing below changes.
   allocate(fix_mask(ndim), source=.false.)
   if (present(frozen)) then
      if (size(frozen) /= ndim) then
         result%message = "the $fix mask does not match the coordinates"
         call env%error("the frozen component mask has the wrong length", source)
         return
      end if
      fix_mask = frozen
   end if
   constrained = any(fix_mask)

   allocate(x(ndim), g(ndim), d(ndim), step(ndim), x_trial(ndim), g_trial(ndim))
   x = x0
   g = 0.0_wp
   d = 0.0_wp
   step = 0.0_wp
   e_trial = 0.0_wp
   g_trial = 0.0_wp
   de = 0.0_wp
   nstagnate = 0
   nreject = 0
   n_evals = 0
   ! The trajectory is allocated **before** the first evaluation: every early
   ! return below (a non-finite start, a failed initial evaluation) calls
   ! `finish`, which copies the trajectory, and a trajectory that was never
   ! allocated would make that copy undefined.  The first frame is the
   ! starting geometry on every path, exactly as it was before (the allocation
   ! used to happen after the first evaluation, which left `frames` and
   ! `nframes` undefined on the failure path; the stack layout of the release
   ! build decided whether that was fatal).
   allocate(frames(ndim, 32))
   nframes = 1
   frames(:, 1) = x

   ! A non-finite starting geometry poisons every downstream comparison (the
   ! line search would silently reject every trial); reject it explicitly.
   if (.not. all(landscape_real_is_finite(x))) then
      result%message = "initial coordinates are not finite"
      return
   end if

   call objective%evaluate(x, energy, g, fail)
   n_evals = n_evals + 1
   if (fail) then
      result%message = "initial energy evaluation failed"
      call finish()
      return
   end if
   if (.not. landscape_real_is_finite(energy) .or. .not. all(landscape_real_is_finite(g))) then
      result%message = "initial energy or gradient is not finite"
      call finish()
      return
   end if

   ! `max_displ` is the global single-step limit; the RFO trust radius is
   ! additionally clamped by it so that no algorithm exceeds the documentated
   ! maximum displacement.
   trust = min(settings%trust_radius, settings%max_displ)
   trust_hi = min(settings%trust_max, settings%max_displ)
   use_trust = settings%algorithm == p_landscape_opt_rfo
   if (use_trust) then
      allocate(hmat(ndim, ndim), hnew(ndim, ndim), source=0.0_wp)
      call objective%hessian(x, hmat, fail)
      if (fail) then
         hmat = 0.0_wp
         do istep = 1, ndim
            hmat(istep, istep) = settings%hlow
         end do
         call env%warning("initial Hessian unavailable, using hlow*I", source)
      end if
      hmat = 0.5_wp*(hmat + transpose(hmat))
   else if (settings%algorithm == p_landscape_opt_bfgs) then
      call new_landscape_bfgs(bfgs, ndim, settings%hessian_update, settings%hlow)
      ! The model Hessian (`hessian_source = 'model'`) is a cheap initial guess
      ! and seeds the quasi-Newton approximation (PLAN section 5.2, D3).  The
      ! expensive numerical Hessian is never built just to seed BFGS, so the
      ! seed is taken from the objective only for a model source.  The name is
      ! matched case insensitively, like `landscape_hessian_source_known` and
      ! the provider factory (review finding M2 of `20-stage1-review.md`).
      if (landscape_lowercase(settings%hessian_source) == 'model') then
         allocate(hinit(ndim, ndim), source=0.0_wp)
         call objective%hessian(x, hinit, fail)
         if (fail) then
            call env%warning("initial model Hessian unavailable, using identity", &
               & source)
         else
            call bfgs%set_hessian(hinit)
         end if
         deallocate(hinit)
      end if
   else if (settings%algorithm == p_landscape_opt_cg) then
      call new_landscape_cg(cg, ndim)
   else if (settings%algorithm == p_landscape_opt_lbfgs) then
      call new_landscape_lbfgs(lbfgs, ndim, settings%lbfgs_memory, &
         & curvature=settings%curvature)
   else if (settings%algorithm == p_landscape_opt_diis) then
      call new_landscape_diis(diis, ndim, settings%diis_memory, &
         & settings%diis_regularization, settings%diis_min_snapshots)
      call diis%push(x, g)
   end if

   do istep = 1, settings%max_steps
      d = 0.0_wp
      predicted = 0.0_wp
      accepted = .false.
      fallback = .false.

      if (use_trust .and. settings%recalc >= 1 .and. istep > 1) then
         if (mod(istep - 1, settings%recalc) == 0) then
            call objective%hessian(x, hnew, fail)
            if (.not. fail) hmat = 0.5_wp*(hnew + transpose(hnew))
         end if
      end if

      ! `method` records which algorithm actually produced the search
      ! direction and `fallback_step` is set whenever the configured algorithm
      ! had to be replaced or completed by steepest descent.  Both are exposed
      ! through the result so that a silent degradation of a quasi-Newton run
      ! to SD can be detected (PLAN section 5.4).
      method = settings%algorithm
      fallback_step = .false.

      select case (settings%algorithm)
      case (p_landscape_opt_sd)
         call landscape_sd_direction(g, d)
      case (p_landscape_opt_cg)
         call cg%direction(g, d)
      case (p_landscape_opt_lbfgs)
         call lbfgs%direction(g, d, ok, steepest)
      case (p_landscape_opt_bfgs)
         call bfgs%direction(env, g, d, ok)
         if (.not. ok) then
            call landscape_sd_direction(g, d)
            method = p_landscape_opt_sd
            fallback_step = .true.
         end if
      case (p_landscape_opt_diis)
         call diis%direction(env, x, d, ok)
         if (.not. ok) then
            call landscape_sd_direction(g, d)
            method = p_landscape_opt_sd
            fallback_step = .true.
         end if
      case (p_landscape_opt_rfo)
         call landscape_rfo_step(env, hmat, g, trust, settings%hlow, d, predicted, &
            & fail, evals_eps=settings%evals_eps, mu_margin=settings%mu_margin, &
            & max_bisect_it=settings%max_bisect_it)
         if (fail) then
            call landscape_sd_direction(g, d)
            method = p_landscape_opt_sd
            fallback_step = .true.
         end if
      case default
         call env%error("unknown minimisation algorithm", source)
         result%message = "unknown minimisation algorithm"
         call finish()
         return
      end select

      ! `$fix`: drop the frozen components of the direction.  For SD/CG/L-BFGS/
      ! GDIIS this is a provable no-op (they are built from the zeroed gradient
      ! and from history vectors that inherit its zeros); BFGS only when its
      ! seed is decoupled (the default `hlow*I` is, a Lindh-D2 seed from
      ! `--lopt-hessian model` is not), and for RFO it is the projection onto
      ! the subspace `$fix` leaves free, which is the whole fix (see
      ! `xtb_landscape_constraints`).  The gradient
      ! vanishes on the frozen components, so the projected direction is still
      ! a descent direction and `predicted` stays an upper bound of the model
      ! decrease of the step that is actually taken.
      if (constrained) call landscape_project_frozen(d, fix_mask)

      if (.not. use_trust) then
         call clip_step(d, settings%max_displ)
      end if

      if (use_trust .and. predicted > 0.0_wp) then
         step = d
         x_trial = x + d
         call objective%evaluate(x_trial, e_trial, g_trial, fail)
         n_evals = n_evals + 1
         if (.not. fail) accepted = e_trial <= energy
         if (accepted) then
            rho = (energy - e_trial)/predicted
            call landscape_rfo_update_trust(trust, rho, settings%eta_shrink, &
               & settings%eta_expand, settings%trust_min, trust_hi)
            nreject = 0
         else
            step = 0.0_wp
            nreject = nreject + 1
            call landscape_rfo_update_trust(trust, 0.0_wp, settings%eta_shrink, &
               & settings%eta_expand, settings%trust_min, trust_hi)
            if (nreject >= rfo_reject_patience) then
               call steepest_fallback(objective, x, energy, g, settings, step, &
                  & e_trial, g_trial, ls_evals, ls_fail, accepted)
               n_evals = n_evals + ls_evals
               fallback_step = .true.
               if (accepted) then
                  nreject = 0
                  fallback = .true.
                  method = p_landscape_opt_sd
               end if
            end if
         end if
      else
         call take_descent_step(objective, x, energy, g, d, settings, step, &
            & e_trial, g_trial, ls_evals, ls_fail, accepted)
         n_evals = n_evals + ls_evals
         if (.not. accepted) then
            call steepest_fallback(objective, x, energy, g, settings, step, &
               & e_trial, g_trial, ls_evals, ls_fail, accepted)
            n_evals = n_evals + ls_evals
            fallback = accepted
            fallback_step = .true.
            if (accepted) method = p_landscape_opt_sd
         end if
      end if

      if (.not. accepted) then
         step = 0.0_wp
         de = 0.0_wp
         ! no progress along the search direction: the stored correction pair
         ! of the conjugate gradient is no longer valid
         if (settings%algorithm == p_landscape_opt_cg) call cg%reset()
      else
         de = e_trial - energy
         select case (settings%algorithm)
         case (p_landscape_opt_cg)
            ! after a steepest descent fallback the accepted step was not taken
            ! along `d`, so the pair is dropped instead of being stored
            if (fallback) then
               call cg%reset()
            else
               call cg%push(g, d)
            end if
         case (p_landscape_opt_lbfgs)
            call lbfgs%push(step, g_trial - g)
         case (p_landscape_opt_bfgs)
            call bfgs%update_hessian(step, g_trial - g, ok)
         case (p_landscape_opt_rfo)
            call update_by_enum(settings%hessian_update, step, g_trial - g, hmat, ok)
         end select
         x = x + step
         energy = e_trial
         g = g_trial
         if (settings%algorithm == p_landscape_opt_diis) then
            if (settings%diis_store_every <= 1 .or. &
                  & mod(istep, settings%diis_store_every) == 0) then
               call diis%push(x, g)
            end if
         end if
         nframes = nframes + 1
         call append_frame(frames, nframes, x)
      end if

      ! A fallback is only a degradation when its step is actually taken; a
      ! rejected fallback (for example once the gradient has already vanished)
      ! does not mean the configured algorithm was replaced by SD.
      if (fallback_step .and. accepted) result%n_fallback = result%n_fallback + 1

      record%istep = istep
      record%energy = energy
      record%delta_e = de
      record%gnorm = vec_norm(g)
      record%gmax = vec_maxabs(g)
      record%dispnorm = vec_norm(step)
      record%disp_max = vec_maxabs(step)
      record%method = method
      record%fallback = fallback_step
      if (use_trust) then
         record%trust_radius = trust
      else
         record%trust_radius = settings%max_displ
      end if
      call landscape_result_add_step(result, record)
      call landscape_report_add(result%conv, record)
      if (settings%printlevel >= 1) call landscape_print_step(record, env%unit)

      ! A rejected step has a zero displacement; its displacement criteria are
      ! meaningless, so only the gradient and energy criteria are tested.
      check_displ = accepted
      if (settings%is_converged(de, g, step, check_displ=check_displ)) then
         result%converged = .true.
         result%message = "converged"
         exit
      end if

      if (.not. accepted) then
         nstagnate = nstagnate + 1
      else
         nstagnate = 0
      end if
      if (nstagnate >= stagnate_patience) then
         result%message = "stagnated: no acceptable step found"
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

   call finish()

contains

   !> Clip a direction to the maximum displacement.
   pure subroutine clip_step(direction, max_displ)
      !> Direction (bohr), scaled in place
      real(wp), intent(inout) :: direction(:)
      !> Maximum displacement (bohr)
      real(wp), intent(in) :: max_displ

      real(wp) :: dnorm

      dnorm = vec_norm(direction)
      if (dnorm > max_displ) direction = direction*(max_displ/dnorm)

   end subroutine clip_step

   !> Take one descent step with or without the Armijo line search.
   !>
   !> When `settings%line_search` is false the full direction is evaluated once
   !> and only accepted if it does not increase the energy.
   subroutine take_descent_step(objective, x, energy, g, d, settings, step, &
         & e_trial, g_trial, n_eval, fail, accepted)
      !> Objective function
      class(landscape_energy_model), intent(inout) :: objective
      !> Current coordinates (bohr)
      real(wp), intent(in) :: x(:)
      !> Current energy (hartree)
      real(wp), intent(in) :: energy
      !> Current gradient (hartree/bohr)
      real(wp), intent(in) :: g(:)
      !> Search direction (bohr)
      real(wp), intent(in) :: d(:)
      !> Optimisation settings
      type(landscape_opt_settings), intent(in) :: settings
      !> Accepted displacement (bohr)
      real(wp), intent(out) :: step(:)
      !> Energy at the accepted point (hartree)
      real(wp), intent(out) :: e_trial
      !> Gradient at the accepted point (hartree/bohr)
      real(wp), intent(out) :: g_trial(:)
      !> Number of objective evaluations
      integer, intent(out) :: n_eval
      !> Failure flag
      logical, intent(out) :: fail
      !> True if a step was accepted
      logical, intent(out) :: accepted

      real(wp) :: alpha

      step = 0.0_wp
      if (settings%line_search) then
         call landscape_line_search(objective, x, energy, g, d, 1.0_wp, step, &
            & alpha, e_trial, g_trial, n_eval, fail, &
            & c1=settings%armijo_c1, max_backtrack=settings%max_backtrack)
         accepted = .not. fail
      else
         n_eval = 1
         call objective%evaluate(x + d, e_trial, g_trial, fail)
         if (.not. fail .and. e_trial <= energy) then
            step = d
            accepted = .true.
         else
            fail = .true.
            accepted = .false.
         end if
      end if

   end subroutine take_descent_step

   !> Fall back to a steepest descent step after a rejected quasi-Newton step.
   subroutine steepest_fallback(objective, x, energy, g, settings, step, &
         & e_trial, g_trial, n_eval, fail, accepted)
      !> Objective function
      class(landscape_energy_model), intent(inout) :: objective
      !> Current coordinates (bohr)
      real(wp), intent(in) :: x(:)
      !> Current energy (hartree)
      real(wp), intent(in) :: energy
      !> Current gradient (hartree/bohr)
      real(wp), intent(in) :: g(:)
      !> Optimisation settings
      type(landscape_opt_settings), intent(in) :: settings
      !> Accepted displacement (bohr)
      real(wp), intent(out) :: step(:)
      !> Energy at the accepted point (hartree)
      real(wp), intent(out) :: e_trial
      !> Gradient at the accepted point (hartree/bohr)
      real(wp), intent(out) :: g_trial(:)
      !> Number of objective evaluations
      integer, intent(out) :: n_eval
      !> Failure flag
      logical, intent(out) :: fail
      !> True if a step was accepted
      logical, intent(out) :: accepted

      real(wp), allocatable :: dir(:)

      allocate(dir(size(g)))
      call landscape_sd_direction(g, dir)
      call clip_step(dir, settings%max_displ)
      call take_descent_step(objective, x, energy, g, dir, settings, step, &
         & e_trial, g_trial, n_eval, fail, accepted)

   end subroutine steepest_fallback

   !> Finalise the result: convergence measures, coordinates and trajectory.
   subroutine finish()
      real(wp), allocatable :: trimmed(:, :)

      result%energy = energy
      result%gnorm = vec_norm(g)
      result%gmax = vec_maxabs(g)
      result%n_eval = n_evals
      call landscape_result_trim(result)
      call landscape_convergence_measures(g, step, gmax, grms, dmax, drms)
      ! A non-finite measure (the caller has already classified such a gradient)
      ! must not be compared: `<` raises the invalid exception on a NaN, so
      ! such a run is reported as not converged instead.
      if (all(landscape_value_present([gmax, grms, dmax, drms]))) then
         result%conv%converged = [gmax <= settings%g_max, grms <= settings%g_rms, &
            & dmax <= settings%d_max, drms <= settings%d_rms]
      else
         result%conv%converged = .false.
      end if
      if (present(x_final)) x_final = x
      if (present(traj)) then
         allocate(trimmed(ndim, nframes), source=frames(:, :nframes))
         call move_alloc(trimmed, traj)
      end if
   end subroutine finish

end subroutine landscape_optimize_model

!> Copy x into the trajectory buffer as frame `n` (1-based).
subroutine append_frame(frames, n, x)
   !> Trajectory buffer `(ndim, capacity)`
   real(wp), allocatable, intent(inout) :: frames(:, :)
   !> 1-based index of the new frame
   integer, intent(in) :: n
   !> New frame
   real(wp), intent(in) :: x(:)

   real(wp), allocatable :: buffer(:, :)
   integer :: capacity

   if (.not. allocated(frames)) allocate(frames(size(x), 32))
   if (n > size(frames, 2)) then
      capacity = max(2*size(frames, 2), n)
      ! `source=` requires the source shape to equal the shape of the
      ! allocate-object; the grown buffer is larger, so the recorded frames are
      ! copied into its front explicitly (the new columns are written before
      ! they are ever read).
      allocate(buffer(size(frames, 1), capacity))
      buffer(:, 1:size(frames, 2)) = frames
      call move_alloc(buffer, frames)
   end if
   frames(:, n) = x

end subroutine append_frame

!> Reset all mutable parts of a result object.
subroutine reset_result(result)
   !> Result to reset
   type(landscape_result), intent(inout) :: result

   result%converged = .false.
   result%n_steps = 0
   result%energy = 0.0_wp
   result%gnorm = 0.0_wp
   result%gmax = 0.0_wp
   result%n_eval = 0
   result%n_cache_miss = 0
   result%n_fallback = 0
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

end module xtb_landscape_optimize
