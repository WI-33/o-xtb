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

!> Dimer transition state search of the landscape suite
!> (PLAN sections 6.3 and 9.5).
!>
!> The dimer method of Henkelman and Jonsson (J. Chem. Phys. **111**, 7010
!> (1999)) locates a first order saddle point using **only first
!> derivatives**.  Two images of the system, separated by `2*delta` around a
!> common midpoint, measure the curvature of the PES along the dimer axis; the
!> axis is rotated into the direction of lowest curvature and the midpoint is
!> translated uphill along that direction while the perpendicular subspace is
!> relaxed.  No Hessian is ever constructed, which makes the dimer the natural
!> complement of the eigenvector following search of the same milestone
!> (`xtb_landscape_saddle_prfo`, which builds a numerical Hessian on every
!> macro iteration).
!>
!> Equations (numbering of Henkelman & Jonsson 1999):
!>
!> ```text
!> C        = (F2 - F1).N/(2 delta) = N^T H N              (Eq. 2)
!> F'       = F1' - F2',   Fi' = Fi - (Fi.N) N            (Eq. 5)
!> dtheta   = -1/2 arctan(2 F0/F0')                       (Eq. 13)
!> F^dagger = F_R - 2 F_par                      if C < 0 (Eq. 21)
!> F^dagger = -F_par                             if C > 0 (Eq. 21)
!> ```
!>
!> **Curvature.**  The paper writes `C = (F2 - F1).N/(2 delta)`, with `F1` and
!> `F2` the forces at the two images and `N` the unit dimer axis, and quotes
!> the equivalent estimate `C = (E2 - E0)/delta^2` from the image energy `E2`
!> and the midpoint energy `E0` of its Eq. 3.  The first form only has a sign
!> once the two images are *named*: with image 1 at `R + delta*N` and image 2
!> at `R - delta*N` (the convention of this module) the force difference is
!> `F2 - F1 = 2 delta (H N)`, so `C = N^T H N` is the **energy curvature**
!> along the axis and the C<0 branch of Eq. 21 becomes the physical one
!> (ascend along the axis, descend perpendicularly).  Swapping the image
!> labels flips the sign of `C` and would invert Eq. 21; the assignment is
!> therefore part of the convention and is pinned by
!> `landscape_dimer_translation_force`.
!>
!> The two forms of Eq. 2 agree at a stationary point along the axis:
!> `(E2 - E0)/delta^2 = C/2 - (g.N)/delta`, so the energy form carries a
!> factor 1/2 relative to the force form (the images sit at `+/-delta` from the
!> midpoint, not at `+/-delta/2`) plus a gradient term that only vanishes for
!> `g.N = 0`.  A further factor 2 appears if the denominator is written with
!> the **image separation** `2*delta`: that is exactly the ambiguity of PLAN
!> section 5.3.6, where VASP's `DdR` denotes twice the image-to-midpoint
!> distance while `delta` here denotes the centre-to-image distance
!> (pysisyphus convention, `1.0e-2_wp*aatoau` bohr = `1e-2` Angstrom).  This
!> module uses the force form with the centre-to-image `delta`, so the same
!> `delta` also sets the denominator `2*delta` and the displacement `+/-`
!> `delta*N` of the two images.
!>
!> **Sign resolution** (M4 hard constraint 2).  `C` is a difference of two
!> forces divided by `2*delta`: with the force noise of a cold SCC start
!> (`sigma_F ~ 1e-9` hartree/bohr, the cold-start SCC gradient noise) and the
!> default `delta = 1.0e-2_wp*aatoau` bohr the statistical error of the curvature
!> is `sqrt(2) sigma_F/(2 delta) ~ 4e-8` hartree/bohr^2.  That is five orders of
!> magnitude below the 3 cm^-1 band (about `1e-4` hartree/bohr^2) quoted by
!> M4 hard constraint 2, so the branch decision is not noise limited for any
!> mode outside that band; a mode *inside* the band is below the resolution of
!> any finite difference curvature and needs the remedy of that constraint
!> (tighter SCC accuracy) rather than a different `delta`.
!>
!> **Translation.**  The midpoint is propagated with quick-min along the
!> effective force `F^dagger` of Eq. 21.  The displacement is *proportional to
!> the force*: `step = tau F^dagger`, with the mobility `tau` fixed once, from
!> the initial effective force, so that `min(step0, step_max)` is the **upper
!> bound** of the first translation and later ones follow the force (and shrink
!> with it, which is what makes the displacement criteria reachable at the
!> saddle).  The bound is only an upper limit: the finite difference `dR` of the
!> paper supplies the curvature along the effective force, which limits the
!> mobility to `1/|kappa|` when that is smaller.  Measured on GFN2 HCN the first
!> translation is 0.1032 bohr (mobility 0.8930 at `kappa = 1.1198` Eh/bohr^2),
!> i.e. the curvature limiter, not `min(step0, step_max) = 0.28346`, sets the
!> length.  A stiff direction therefore cannot overshoot its harmonic length
!> scale, and a negative-curvature direction is bounded by the same `1/|kappa|`
!> instead of growing without limit.
!>
!> *Measured deviation from the paper's step length rule.*  Taking the step
!> length itself as the Newton length `|F^dagger|/kappa` (the literal reading of
!> the paper's one dimensional minimisation along `F^dagger` with step `dR`)
!> makes the search **ping-pong** around the saddle: the effective force is not
!> a gradient field - its parallel part ascends while its perpendicular part
!> descends - so `kappa = u^T H u` changes sign near the transition state and
!> `|F^dagger|/kappa` stops shrinking.  Measured on Mueller-Brown (delta =
!> `1.0e-2_wp*aatoau` bohr): from `(-0.6, 0.4)` the Newton form keeps taking
!> the maximal step near the saddle and orbits it forever (200 steps, distance
!> the saddle alternating between 0.036 and 0.015, `kappa` between -1.20 and
!> -0.96), while the force proportional form below converges in 23 steps to
!> 1.4e-4 bohr of the saddle.  The mobility limiter `1/kappa` did not bind in
!> that run (the largest `kappa` on the path was 0.9, against `1/tau = 2.7`);
!> it is a guard, not the step control.
!>
!> Conventions:
!>
!>  * flattened geometries are `xyz(3, n)` (`x1,y1,z1,x2,...`) and displacement
!>    vectors have length `3n` in the same order;
!>  * the dimer axis is a *line*: `N` and `-N` describe the same dimer, so no
!>    phase convention is needed inside the algorithm.  The reported
!>    transition vector is nevertheless phase fixed (largest component made
!>    positive) by `xtb_landscape_saddle_eigenfollow`, so that the same
!>    geometry always prints the same direction;
!>  * the random initial axis is the only place in the suite where a random
!>    number is drawn, and it is seeded explicitly from `settings%seed`
!>    (PLAN N10);
!>  * the Cartesian translation step and the axis are projected onto the
!>    complement of the translations and rotations, so neither the displacement
!>    bound (`step_max`) nor the rotation is spent on rigid body motion.
!>
!> Units: coordinates, `delta`, `dR`, `step0`, `step_max` and displacements in
!> bohr, energies in hartree, forces and gradients in hartree/bohr, curvatures
!> in hartree/bohr^2, rotation angles (`rot_dtheta`, Eq. 13) in rad.
module xtb_landscape_saddle_dimer
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment, only : TEnvironment
   use xtb_landscape_settings, only : landscape_ts_settings, &
      & landscape_convergence_measures, landscape_value_present, &
      & landscape_lowercase
   use xtb_landscape_results, only : landscape_result, landscape_step, &
      & landscape_result_add_step, landscape_result_trim, landscape_report_add, &
      & landscape_print_step, p_landscape_bound_step_max
   use xtb_landscape_coords, only : vec_norm, vec_maxabs, vec_rms
   use xtb_landscape_constraints, only : landscape_project_fixed_atoms
   use xtb_landscape_engine, only : landscape_energy_model
   use xtb_landscape_saddle_eigenfollow, only : landscape_rigid_modes, &
      & landscape_project_rigid, landscape_transition_mode, &
      & landscape_transition_mode_report, landscape_transition_mode_set_phase
   implicit none(type, external)
   private

   public :: landscape_dimer_optimize
   public :: landscape_dimer_curvature
   public :: landscape_dimer_translation_force
   public :: landscape_dimer_rotation_direction
   public :: landscape_dimer_rotation_angle
   public :: landscape_dimer_velocity_update
   public :: landscape_dimer_initial_axis
   public :: landscape_dimer_rotate
   public :: landscape_dimer_mobility

   !> Curvature floor of the translation step (hartree/bohr^2).  Below it the
   !> Newton length `|F^dagger|/kappa` of the quick-min update is not usable (a
   !> flat or negative curvature along the effective force) and the step falls
   !> back to `step0`.
   real(wp), parameter :: dimer_kappa_floor = 1.0e-6_wp

   !> Norm below which a vector counts as zero
   real(wp), parameter :: dimer_tiny = 1.0e-12_wp

   !> Pi, only needed by the degenerate limit of Eq. 13
   real(wp), parameter :: dimer_pi = 3.14159265358979323846_wp

contains

!> Curvature of the PES along the dimer axis (H&J Eq. 2).
!>
!> `f_plus` and `f_minus` are the **forces** (`-gradient`) at the images
!> `R + delta*N` and `R - delta*N`.  With those image labels the result is the
!> energy curvature `N^T H N`: negative when the axis points along the
!> reaction coordinate, positive otherwise.  `delta > 0` is required; the
!> caller validates it.
pure function landscape_dimer_curvature(f_plus, f_minus, axis, delta) result(curvature)
   !> Force at the image `R + delta*N` (hartree/bohr)
   real(wp), intent(in) :: f_plus(:)
   !> Force at the image `R - delta*N` (hartree/bohr)
   real(wp), intent(in) :: f_minus(:)
   !> Unit dimer axis
   real(wp), intent(in) :: axis(:)
   !> Centre-to-image distance (bohr, positive)
   real(wp), intent(in) :: delta
   !> Energy curvature along the axis (hartree/bohr^2)
   real(wp) :: curvature

   curvature = dot_product(f_minus - f_plus, axis)/(2.0_wp*delta)

end function landscape_dimer_curvature

!> Effective translation force of the dimer (H&J Eq. 21).
!>
!> `F_par = (F_R.N) N` is the component of the midpoint force along the axis.
!> With a negative curvature the axis is the reaction coordinate: the
!> perpendicular component is kept (the perpendicular subspace is relaxed)
!> and the parallel component is reversed (the axis is climbed).  With a
!> positive curvature the lowest curvature mode is still the direction to
!> ascend, but the perpendicular relaxation is deliberately dropped, as in the
!> paper: near a minimum it would only pull the dimer back into the basin it
!> tries to leave.
pure subroutine landscape_dimer_translation_force(force, axis, curvature, &
      & force_dagger)
   !> Force at the midpoint, `F_R = -g` (hartree/bohr)
   real(wp), intent(in) :: force(:)
   !> Unit dimer axis
   real(wp), intent(in) :: axis(:)
   !> Curvature along the axis (hartree/bohr^2)
   real(wp), intent(in) :: curvature
   !> Effective translation force (hartree/bohr)
   real(wp), intent(out) :: force_dagger(:)

   real(wp) :: f_par(size(force))

   f_par = dot_product(force, axis)*axis
   if (curvature < 0.0_wp) then
      force_dagger = force - 2.0_wp*f_par
   else
      force_dagger = -f_par
   end if

end subroutine landscape_dimer_translation_force

!> Direction and magnitude of the rotational force (H&J Eq. 5).
!>
!> `F' = F1' - F2'` is the difference of the image forces with the parallel
!> component removed.  Because `F1 - F2 = -2 delta (H N)`, `F'` is proportional
!> to minus the projection of `H N` on the unit sphere, i.e. it points into the
!> direction in which the axis has to be rotated to **lower** the curvature.
!> `rms` normalises `magnitude` by the number of Cartesian components and is
!> the quantity compared against `rot_f_rms_th`; `direction` is normalised and
!> left at zero when the rotational force vanishes (the axis is already an
!> eigenvector of the Hessian).  `force_rot` optionally returns the
!> perpendicular force difference itself, so the caller can report its largest
!> component.
pure subroutine landscape_dimer_rotation_direction(f_plus, f_minus, axis, &
      & direction, magnitude, rms, force_rot)
   !> Force at the image `R + delta*N` (hartree/bohr)
   real(wp), intent(in) :: f_plus(:)
   !> Force at the image `R - delta*N` (hartree/bohr)
   real(wp), intent(in) :: f_minus(:)
   !> Unit dimer axis
   real(wp), intent(in) :: axis(:)
   !> Unit direction of the rotational force
   real(wp), intent(out) :: direction(:)
   !> Norm of the perpendicular force difference (hartree/bohr)
   real(wp), intent(out) :: magnitude
   !> RMS of the perpendicular force difference (hartree/bohr)
   real(wp), intent(out) :: rms
   !> Perpendicular force difference (hartree/bohr), optional
   real(wp), intent(out), optional :: force_rot(:)

   real(wp) :: f_rot(size(f_plus))
   real(wp) :: parallel

   f_rot = f_plus - f_minus
   parallel = dot_product(f_rot, axis)
   f_rot = f_rot - parallel*axis
   magnitude = vec_norm(f_rot)
   rms = vec_rms(f_rot)
   if (magnitude > dimer_tiny) then
      direction = f_rot/magnitude
   else
      direction = 0.0_wp
   end if
   if (present(force_rot)) force_rot = f_rot

end subroutine landscape_dimer_rotation_direction

!> Optimal rotation angle of the dimer, Eq. 13 of Henkelman & Jonsson.
!>
!> The rotational force along the probe direction, `F(theta)`, is a two term
!> Fourier series in `2*theta` (`A sin 2theta + B cos 2theta`; the period is
!> `pi` because the dimer axis is a line).  With `F0 = F(0) = B` and
!> `F0' = dF/dtheta|_0 = 2A`, the two roots of the model are `pi/2` apart and
!> the rotation that **lowers** the curvature is the first positive one,
!> `theta = 0.5 atan2(-2 F0, F0')` (plus `pi/2` when that is negative).  The
!> earlier Eq. 6 of the paper, the linear interpolation
!> `(F + F*)/(-2F')`, differs from Eq. 13 by `O(dtheta)` and is not
!> implemented.
!>
!> Taking the principal branch of the paper's form `-1/2 arctan(2 F0/F0')`
!> instead is not equivalent: it selects the root **nearest to zero**, which is
!> the maximum of the curvature along the probe direction whenever the model
!> needs a rotation of more than 45 degrees, and the search then rotates
!> *uphill*.  Measured on the GFN2 HCN start geometry (an axis chosen by the
!> largest gradient component): the principal branch returns -0.54 rad where
!> the correct root is +1.04 rad, and the search follows the maximising
!> direction into a stiff mode (curvature +2.6 Eh/bohr^2, energy rising to
!> +1e2 Eh); with the branch above the same start keeps the axis on the
!> negative mode (curvature -0.126 Eh/bohr^2).
!>
!> A vanishing force means a converged rotation and yields a zero angle; a
!> vanishing derivative degenerates to the limit of the expression, `pi/4`.
pure function landscape_dimer_rotation_angle(f0, f0_prime) result(theta)
   !> Rotational force at the current orientation (hartree/bohr, positive)
   real(wp), intent(in) :: f0
   !> Finite difference derivative of the rotational force with respect to the
   !> rotation angle (hartree/bohr per rad)
   real(wp), intent(in) :: f0_prime
   !> Optimal rotation angle (rad), positive towards the rotational force
   real(wp) :: theta

   theta = 0.0_wp
   if (f0 <= 0.0_wp) return
   theta = 0.5_wp*atan2(-2.0_wp*f0, f0_prime)
   if (theta < 0.0_wp) theta = theta + 0.5_wp*dimer_pi

end function landscape_dimer_rotation_angle

!> One quick-min displacement update of the dimer translation.
!>
!> Quick-min (H&J section II D; Sheppard, Terrell & Henkelman, J. Chem. Phys.
!> **128**, 134106 (2008)) propagates a velocity along the effective force and
!> keeps its own momentum:
!>
!> ```text
!> v    <- v + (increment/|F^dagger|) F^dagger
!> v    <- 0                                if v.F^dagger <= 0
!> v    <- v * min(1, step_max/|v|)         (bounded momentum)
!> step = v
!> ```
!>
!> The velocity of an uphill move (one that opposes the force that produced it)
!> is discarded, so a step never continues against that force; the momentum is
!> capped at `step_max`, which therefore bounds both the displacement of one
!> update and the momentum the next update can add to.  With a zero incoming
!> velocity the displacement has exactly the length `increment`, which the
!> caller derives from the mobility of the translation and the magnitude of the
!> effective force.
pure subroutine landscape_dimer_velocity_update(force, velocity, increment, &
      & step_max, step)
   !> Effective translation force (hartree/bohr)
   real(wp), intent(in) :: force(:)
   !> Quick-min velocity (bohr), updated in place
   real(wp), intent(inout) :: velocity(:)
   !> Length added along the force by this update (bohr)
   real(wp), intent(in) :: increment
   !> Upper bound of the displacement and of the momentum (bohr)
   real(wp), intent(in) :: step_max
   !> Displacement of this update (bohr)
   real(wp), intent(out) :: step(:)

   real(wp) :: f_norm, v_norm

   step = 0.0_wp
   f_norm = vec_norm(force)
   if (f_norm <= dimer_tiny) then
      velocity = 0.0_wp
      return
   end if
   velocity = velocity + (increment/f_norm)*force
   if (dot_product(velocity, force) <= 0.0_wp) then
      velocity = 0.0_wp
      return
   end if
   v_norm = vec_norm(velocity)
   if (v_norm > step_max) velocity = velocity*(step_max/v_norm)
   step = velocity

end subroutine landscape_dimer_velocity_update

!> Effective mobility of one quick-min translation.
!>
!> `mobility` is the displacement produced by a unit effective force
!> (bohr^2/hartree) and `kappa` the curvature of the PES along that force,
!> measured with the paper's finite difference step `dR`.  A non-zero
!> curvature bounds the displacement by the harmonic length scale `1/|kappa|`
!> of the followed direction, so neither a stiff direction (`kappa > 0`) nor a
!> direction that is being climbed (`kappa < 0`) can overshoot; without that
!> bound the quick-min velocity of a negative-curvature direction grows without
!> limit (measured on a one dimensional model, report 27 section 4).  A flat
!> direction (`|kappa| <= dimer_kappa_floor`) leaves the propagation mobility
!> unchanged rather than reversing the step.
pure function landscape_dimer_mobility(mobility, kappa) result(tau)
   !> Displacement per unit force of the quick-min propagation (bohr^2/hartree)
   real(wp), intent(in) :: mobility
   !> Curvature along the effective force (hartree/bohr^2)
   real(wp), intent(in) :: kappa
   !> Effective mobility (bohr^2/hartree)
   real(wp) :: tau

   tau = mobility
   if (abs(kappa) > dimer_kappa_floor) tau = min(mobility, 1.0_wp/abs(kappa))

end function landscape_dimer_mobility

!> Initial dimer axis from the settings.
!>
!> `n_init = 'random'` draws a uniform direction and is the only place in the
!> suite where a random number is used.  The generator state is set from
!> `seed` before the draw (PLAN N10), so the axis is a pure function of the
!> seed: two runs with the same seed are bitwise identical, which the unit
!> tests assert.  `n_init = 'force'` is the deterministic substitute: the unit
!> vector of the largest Cartesian gradient component, i.e. the Cartesian
!> direction in which the starting geometry is displaced most.
subroutine landscape_dimer_initial_axis(env, gradient, n_init, seed, axis, fail)
   character(len=*), parameter :: source = 'landscape_dimer_initial_axis'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Gradient at the starting geometry (hartree/bohr)
   real(wp), intent(in) :: gradient(:)
   !> Source of the axis, `random` or `force`
   character(len=*), intent(in) :: n_init
   !> Seed of the generator (`random` only)
   integer, intent(in) :: seed
   !> Unit initial axis
   real(wp), intent(out) :: axis(:)
   !> Failure flag
   logical, intent(out) :: fail

   integer, allocatable :: seed_array(:)
   integer :: nseed, imax
   real(wp) :: nrm

   fail = .false.
   axis = 0.0_wp
   if (size(axis) /= size(gradient) .or. size(axis) == 0) then
      call env%error("the dimer axis and the gradient must have the same "//&
         & "non-zero length", source)
      fail = .true.
      return
   end if

   ! `landscape_ts_ninit_known` accepts any spelling, so the selection has to
   ! as well; the rejection message below echoes the name as it was given.
   select case (landscape_lowercase(n_init))
   case ('force')
      imax = maxloc(abs(gradient), dim=1)
      axis(imax) = 1.0_wp
   case ('random')
      call random_seed(size=nseed)
      allocate(seed_array(nseed))
      seed_array = seed
      call random_seed(put=seed_array)
      call random_number(axis)
      axis = axis - 0.5_wp
   case default
      call env%error("unknown dimer initial axis '"//trim(n_init)//"'", source)
      fail = .true.
      return
   end select

   nrm = vec_norm(axis)
   if (nrm <= dimer_tiny) then
      call env%error("the dimer initial axis is zero", source)
      fail = .true.
      return
   end if
   axis = axis/nrm

end subroutine landscape_dimer_initial_axis

!> Rotate the dimer axis into the direction of lowest curvature.
!>
!> Each iteration evaluates the two image forces at the current axis (Eq. 5),
!> stops as soon as the RMS of the perpendicular force difference is below
!> `rot_f_rms_th`, and otherwise probes the rotational force one probe angle
!> `rot_dtheta` further along the rotational force to build the finite
!> difference derivative used by Eq. 13.  The applied rotation is
!> `rot_alpha*(Eq. 13)` (`rot_alpha = 1` takes the full angle of the two term
!> Fourier model), and the axis is re-projected onto the complement of the
!> rigid body directions and renormalised afterwards.
!>
!> Exhausting `rot_max_iter` while the torque is still above the threshold is
!> the failure mode the VASP manual warns about (a torque that does not fall
!> because the forces are not accurate enough): it is reported with the
!> diagnostic warning "... dimer rotation did not converge ...".  That warning
!> is not fatal - a partially rotated axis still gives a usable search
!> direction - so the caller continues with the best axis found.
!>
!> On return `f_plus`/`f_minus` are the image forces **of the returned axis**
!> (one extra evaluation pair is spent when the last action was a rotation), so
!> the caller reads the curvature without paying for another pair.  `fail` is
!> reserved for a failed *image* evaluation, i.e. for a state without any axial
!> forces; a failed probe only ends the rotation.
subroutine landscape_dimer_rotate(env, objective, x, settings, mode, rigid, &
      & n_rigid, f_plus, f_minus, iterations, torque_rms, torque_max, &
      & converged, n_eval, fail)
   character(len=*), parameter :: source = 'landscape_dimer_rotate'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Objective function
   class(landscape_energy_model), intent(inout) :: objective
   !> Midpoint geometry (bohr)
   real(wp), intent(in) :: x(:)
   !> Transition state settings (`delta`, `rot_dtheta`, `rot_max_iter`,
   !> `rot_alpha`, `rot_f_rms_th`)
   type(landscape_ts_settings), intent(in) :: settings
   !> Transition mode carrying the axis, rotated in place
   type(landscape_transition_mode), intent(inout) :: mode
   !> Orthonormal rigid body basis of the midpoint geometry
   real(wp), intent(in) :: rigid(:, :)
   !> Number of valid rigid body directions
   integer, intent(in) :: n_rigid
   !> Force at the image `R + delta*N` of the returned axis (hartree/bohr)
   real(wp), allocatable, intent(out) :: f_plus(:)
   !> Force at the image `R - delta*N` of the returned axis (hartree/bohr)
   real(wp), allocatable, intent(out) :: f_minus(:)
   !> Number of rotations performed
   integer, intent(out) :: iterations
   !> RMS of the perpendicular force difference at the returned axis
   real(wp), intent(out) :: torque_rms
   !> Maximum absolute component of that difference (hartree/bohr)
   real(wp), intent(out) :: torque_max
   !> True if the torque fell below the threshold
   logical, intent(out) :: converged
   !> Number of objective evaluations spent
   integer, intent(out) :: n_eval
   !> Failure flag (a failed image evaluation)
   logical, intent(out) :: fail

   real(wp), allocatable :: q(:), q_probe(:), axis_probe(:), f_rot(:)
   real(wp), allocatable :: f_plus_probe(:), f_minus_probe(:)
   real(wp) :: torque, torque_probe, f0_prime, theta, c, s, dummy
   logical :: have_forces
   integer :: iter, max_iter

   fail = .false.
   converged = .false.
   n_eval = 0
   iterations = 0
   torque_rms = 0.0_wp
   torque_max = 0.0_wp
   allocate(q(size(x)), q_probe(size(x)), axis_probe(size(x)), f_rot(size(x)))
   allocate(f_plus_probe(size(x)), f_minus_probe(size(x)))
   allocate(f_plus(size(x)), f_minus(size(x)))
   max_iter = max(1, settings%rot_max_iter)
   have_forces = .false.

   do iter = 1, max_iter
      iterations = iter
      call dimer_image_forces(objective, x, mode%vector, settings%delta, &
         & f_plus, f_minus, fail)
      n_eval = n_eval + 2
      if (fail) return
      have_forces = .true.
      call landscape_dimer_rotation_direction(f_plus, f_minus, mode%vector, &
         & q, torque, torque_rms, force_rot=f_rot)
      torque_max = vec_maxabs(f_rot)
      call keep_internal(q, rigid, n_rigid, mode%vector)
      if (torque_rms <= settings%rot_f_rms_th) then
         converged = .true.
         exit
      end if

      ! probe the rotational force one probe angle further along `q`
      c = cos(settings%rot_dtheta)
      s = sin(settings%rot_dtheta)
      axis_probe = c*mode%vector + s*q
      call dimer_image_forces(objective, x, axis_probe, settings%delta, &
         & f_plus_probe, f_minus_probe, fail)
      n_eval = n_eval + 2
      if (fail) then
         call env%warning("the dimer rotation probe failed; keeping the "//&
            & "current axis", source)
         fail = .false.
         exit
      end if
      have_forces = .false.
      call landscape_dimer_rotation_direction(f_plus_probe, f_minus_probe, &
         & axis_probe, q_probe, torque_probe, dummy)
      f0_prime = (torque_probe - torque)/settings%rot_dtheta
      theta = settings%rot_alpha*landscape_dimer_rotation_angle(torque, &
         & f0_prime)
      mode%vector = cos(theta)*mode%vector + sin(theta)*q
      call normalise_axis(mode%vector)
      if (settings%printlevel >= 2) then
         write(env%unit, '(1x,a,i0,a,es14.6,a,f12.8,a,es14.6,a)') &
            & "  dimer rotation ", iter, ": torque rms ", torque, &
            & " Eh/bohr, angle ", theta, " rad, probe ", torque_probe, &
            & " Eh/bohr"
      end if
   end do

   ! The image forces have to belong to the axis that is returned; after a
   ! rotation (and after an exhausted loop) they do not, so they are rebuilt.
   if (.not. have_forces) then
      call dimer_image_forces(objective, x, mode%vector, settings%delta, &
         & f_plus, f_minus, fail)
      n_eval = n_eval + 2
      if (fail) return
      call landscape_dimer_rotation_direction(f_plus, f_minus, mode%vector, &
         & q, torque, torque_rms, force_rot=f_rot)
      torque_max = vec_maxabs(f_rot)
   end if

   if (.not. converged) then
      ! VASP diagnostic (PLAN section 6.3): a torque that does not fall while
      ! the rotation is iterated means the forces are not accurate enough, so
      ! either a larger separation or a tighter SCC is needed.
      call env%warning("dimer rotation did not converge; consider increasing "//&
         & "delta or tightening SCF", source)
   end if

end subroutine landscape_dimer_rotate

!> Restrict a rotation direction to the internal (non-rigid) space.
!>
!> The rotations and translations of a molecule are - in exact arithmetic -
!> null vectors of the Hessian, so they carry the lowest curvature of the
!> surface.  The rotational force of a *numerical* Hessian therefore has a
!> spurious component along them (measured on GFN2 HCN it dominates), and a
!> dimer that follows it rotates its axis out of the internal space: the
!> projection that follows removes the component again but wastes every
!> rotation on it.  Removing the rigid directions from the *direction* keeps
!> the whole rotation inside the internal space, in which the axis already
!> lives.  The component along the axis itself is removed as well, so that a
!> rotation cannot turn into a (meaningless) rescaling of the axis.
pure subroutine keep_internal(direction, rigid, n_rigid, axis)
   !> Rotation direction, modified in place
   real(wp), intent(inout) :: direction(:)
   !> Orthonormal rigid body basis
   real(wp), intent(in) :: rigid(:, :)
   !> Number of valid rigid body directions
   integer, intent(in) :: n_rigid
   !> Current dimer axis
   real(wp), intent(in) :: axis(:)

   real(wp) :: parallel

   if (n_rigid > 0) call landscape_project_rigid(direction, rigid, n_rigid)
   parallel = dot_product(direction, axis)
   direction = direction - parallel*axis
   call normalise_axis(direction)

end subroutine keep_internal

!> Macro driver of the dimer transition state search.
!>
!> One macro iteration is the rotation of the axis followed by one (or, with
!> `translate_max_iter > 1`, several) translation of the midpoint, as in
!> Henkelman & Jonsson:
!>
!>  * the axis is rotated into the lowest curvature direction
!>    (`landscape_dimer_rotate`);
!>  * the curvature along the axis follows from the two image forces of that
!>    rotation (Eq. 2) and gives the effective translation force of Eq. 21;
!>  * the translation is a quick-min propagation along the effective force; its
!>    mobility is fixed by the initial force (`min(step0, step_max)` bounds the
!>    first translation; the `dR` curvature limiter can shorten it further),
!>    length of the first translation), limited by the curvature that the `dR`
!>    probe measures along that force, and bounded by `step_max`.
!>
!> The search has converged when the inherited four criteria are met **and**
!> the curvature along the final axis is negative (the VASP convergence rule):
!> a vanishing gradient with a positive curvature along the axis belongs to a
!> minimum, not to the requested saddle point, and is reported as a warning
!> instead of a result.
!>
!> `x_final` returns the final geometry and `traj` the recorded trajectory
!> (`(ndim, nframes)`, first frame is the start).
subroutine landscape_dimer_optimize(env, objective, settings, x0, result, &
      & x_final, traj)
   character(len=*), parameter :: source = 'landscape_dimer_optimize'
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
   !> Final coordinates (bohr)
   real(wp), allocatable, intent(out), optional :: x_final(:)
   !> Recorded trajectory (bohr)
   real(wp), allocatable, intent(out), optional :: traj(:, :)

   type(landscape_transition_mode) :: mode
   type(landscape_step) :: record
   real(wp), allocatable :: x(:), g(:), g_probe(:), unit(:), force(:)
   real(wp), allocatable :: force_dagger(:), velocity(:), step(:), x_trial(:)
   real(wp), allocatable :: rigid(:, :), f_plus(:), f_minus(:), frames(:, :)
   real(wp), allocatable :: axis(:)
   real(wp) :: energy, e_trial, de, curvature, kappa, h, step_norm
   real(wp) :: torque_rms, torque_max, g_max, g_rms, d_max, d_rms, mobility
   real(wp) :: tau_eff
   logical :: fail, rot_converged, reached, warned_minimum, unusable
   integer :: istep, itrans, ndim, nframes, n_rigid, n_rot, n_eval, n_eval_rot

   if (size(x0) == 0) return
   call reset_result(result)
   ndim = size(x0)

   allocate(x(ndim), g(ndim), g_probe(ndim), unit(ndim), force(ndim))
   allocate(force_dagger(ndim), velocity(ndim), step(ndim), x_trial(ndim))
   allocate(f_plus(ndim), f_minus(ndim), axis(ndim))
   x = x0
   g = 0.0_wp
   step = 0.0_wp
   velocity = 0.0_wp
   energy = 0.0_wp
   de = 0.0_wp
   curvature = 0.0_wp
   torque_rms = 0.0_wp
   torque_max = 0.0_wp
   mobility = 0.0_wp
   tau_eff = 0.0_wp
   n_eval = 0
   warned_minimum = .false.

   if (.not. all(landscape_value_present(x))) then
      result%message = "initial coordinates are not finite"
      call finish()
      return
   end if
   if (.not. dimer_settings_valid(settings)) then
      call env%error("the dimer parameters delta, dR, rot_dtheta, step0 and "//&
         & "step_max must be positive", source)
      result%message = "invalid dimer parameters"
      call finish()
      return
   end if

   ! The rigid body basis is always built so that an allocated array can be
   ! handed to the rotation; `n_rigid = 0` (below) disables the projection.
   call landscape_rigid_modes(x, rigid, n_rigid)
   if (.not. settings%remove_rigid) n_rigid = 0

   call objective%evaluate(x, energy, g, fail)
   n_eval = n_eval + 1
   if (fail .or. .not. landscape_value_present(energy) .or. &
         & .not. all(landscape_value_present(g))) then
      result%message = "initial energy or gradient is not finite"
      call finish()
      return
   end if

   call landscape_dimer_initial_axis(env, g, settings%n_init, settings%seed, &
      & axis, fail)
   if (fail) then
      result%message = "initial dimer axis unavailable"
      call finish()
      return
   end if
   if (n_rigid > 0) then
      call landscape_project_rigid(axis, rigid, n_rigid)
      call normalise_axis(axis)
      if (vec_norm(axis) <= dimer_tiny) then
         call env%error("the dimer initial axis is entirely a rigid body "//&
            & "direction", source)
         result%message = "initial dimer axis unavailable"
         call finish()
         return
      end if
   end if
   mode%vector = axis
   mode%index = 0
   mode%overlap = 1.0_wp
   mode%followed = .true.
   mode%target = 1
   if (settings%printlevel > 0) then
      write(env%unit, '(1x,a,a,a,i0,a)') "dimer initial axis: ", &
         & trim(settings%n_init), " direction (seed ", settings%seed, ")"
      write(env%unit, '(1x,a,es13.6,a,es13.6,a)') "dimer separation:  ", &
         & settings%delta, " bohr, finite difference step ", settings%dR, &
         & " bohr"
   end if

   allocate(frames(ndim, 32))
   nframes = 1
   frames(:, 1) = x

   do istep = 1, settings%max_steps
      step = 0.0_wp
      de = 0.0_wp
      kappa = 0.0_wp
      tau_eff = mobility
      reached = .false.
      if (settings%remove_rigid) call landscape_rigid_modes(x, rigid, n_rigid)

      ! 1. rotate the axis into the direction of lowest curvature; the image
      !    forces of the returned axis give the curvature without extra work
      call landscape_dimer_rotate(env, objective, x, settings, mode, rigid, &
         & n_rigid, f_plus, f_minus, n_rot, torque_rms, torque_max, &
         & rot_converged, n_eval_rot, fail)
      n_eval = n_eval + n_eval_rot
      if (fail) then
         result%message = "a dimer image evaluation failed"
         exit
      end if
      curvature = landscape_dimer_curvature(f_plus, f_minus, mode%vector, &
         & settings%delta)
      mode%eigenvalue = curvature

      ! 2. effective translation force of Eq. 21
      force = -g
      call landscape_dimer_translation_force(force, mode%vector, curvature, &
         & force_dagger)

      ! 3. quick-min translation, bounded by step_max.  The mobility is set
      !    once, from the first effective force, so that min(step0, step_max)
      !    bounds the first displacement; the dR curvature probe can shorten it.
      step_norm = vec_norm(force_dagger)
      if (mobility <= 0.0_wp .and. step_norm > dimer_tiny) then
         mobility = min(settings%step0, settings%step_max)/step_norm
      end if
      do itrans = 1, max(1, settings%translate_max_iter)
         step_norm = vec_norm(force_dagger)
         if (step_norm <= dimer_tiny) exit
         unit = force_dagger/step_norm
         kappa = 0.0_wp
         unusable = .false.
         call objective%evaluate(x + settings%dR*unit, e_trial, g_probe, fail)
         n_eval = n_eval + 1
         if (fail .or. .not. landscape_value_present(e_trial) .or. &
               & .not. all(landscape_value_present(g_probe))) then
            fail = .false.
            unusable = .true.
         else
            ! kappa = u^T H u from the force difference over the probe step
            kappa = -dot_product((-g_probe) - force, unit)/settings%dR
         end if
         if (unusable .and. settings%printlevel >= 2) then
            write(env%unit, '(1x,a)') "  dimer translation: the dR probe "//&
               & "failed, keeping the propagation mobility"
         end if
         tau_eff = landscape_dimer_mobility(mobility, kappa)
         h = tau_eff*step_norm
         call landscape_dimer_velocity_update(force_dagger, velocity, &
            & min(h, settings%step_max), settings%step_max, step)
         if (vec_norm(step) <= dimer_tiny) then
            ! quick-min discards a velocity that opposes the effective force,
            ! which would leave this macro iteration idle (a zero displacement
            ! still costs the rotation and the probe); fall back to a fresh
            ! step along the effective force so the iteration makes progress.
            step = min(tau_eff*step_norm, settings%step_max)*unit
         end if
         if (n_rigid > 0) call landscape_project_rigid(step, rigid, n_rigid)
         ! `$fix`: the translation step must not move a fixed atom.  The dimer
         ! does not need a Hessian, but the mobility `tau` scales the *whole*
         ! force vector, and the finite difference force of a fixed atom is
         ! zero only up to the SCF gradient zeroing of the calculator, so the
         ! projection is applied here exactly like the minimiser does it.
         call landscape_project_fixed_atoms(step)
         if (vec_norm(step) <= dimer_tiny) exit

         x_trial = x + step
         call objective%evaluate(x_trial, e_trial, g_probe, fail)
         n_eval = n_eval + 1
         if (fail .or. .not. landscape_value_present(e_trial) .or. &
               & .not. all(landscape_value_present(g_probe))) then
            call env%warning("the dimer translation failed; keeping the "//&
               & "previous geometry", source)
            fail = .false.
            velocity = 0.0_wp
            exit
         end if
         de = e_trial - energy
         x = x_trial
         g = g_probe
         energy = e_trial

         if (settings%is_converged(de, g, step, check_displ=.true.)) then
            reached = .true.
            exit
         end if
      end do

      nframes = nframes + 1
      call append_frame(frames, nframes, x)
      call record_step(settings, istep, energy, de, g, step, curvature, n_rot, &
         & torque_rms, torque_max, rot_converged, vec_norm(force_dagger), &
         & tau_eff, kappa)

      if (reached) then
         if (curvature < 0.0_wp) then
            result%converged = .true.
            result%message = "converged"
            exit
         end if
         if (.not. warned_minimum) then
            call env%warning("the convergence criteria are met but the "//&
               & "curvature along the dimer axis is not negative: the point "//&
               & "is a minimum rather than a saddle point", source)
            warned_minimum = .true.
         end if
      end if
   end do

   if (.not. allocated(result%message)) then
      if (settings%max_steps >= 1) then
         result%message = "maximum number of steps reached"
      else
         result%message = "no step performed"
      end if
   end if

   call landscape_transition_mode_set_phase(mode)
   if (settings%printlevel > 0) then
      write(env%unit, '(1x,a,es15.6,a)') "dimer final curvature: ", curvature, &
         & " Eh/bohr^2"
      call landscape_transition_mode_report(mode, env%unit)
   end if

   result%energy = energy
   result%gnorm = vec_norm(g)
   result%gmax = vec_maxabs(g)
   result%n_eval = n_eval
   call landscape_convergence_measures(g, step, g_max, g_rms, d_max, d_rms)
   result%conv%converged = [g_max <= settings%g_max, g_rms <= settings%g_rms, &
      & d_max <= settings%d_max, d_rms <= settings%d_rms]
   call finish()

contains

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

   !> Append one macro step to the result and print it together with the dimer
   !> diagnostics (rotation count and torque, curvature, displacement).
   subroutine record_step(settings_, istep, energy_, de_, g_, step_, &
         & curvature_, n_rot_, torque_rms_, torque_max_, rot_converged_, &
         & force_norm_, mobility_, kappa_)
      !> Settings
      type(landscape_ts_settings), intent(in) :: settings_
      !> Iteration number
      integer, intent(in) :: istep
      !> Energy (hartree)
      real(wp), intent(in) :: energy_
      !> Energy change (hartree)
      real(wp), intent(in) :: de_
      !> Gradient (hartree/bohr)
      real(wp), intent(in) :: g_(:)
      !> Displacement of the macro step (bohr)
      real(wp), intent(in) :: step_(:)
      !> Curvature along the axis (hartree/bohr^2)
      real(wp), intent(in) :: curvature_
      !> Number of rotations
      integer, intent(in) :: n_rot_
      !> RMS torque (hartree/bohr)
      real(wp), intent(in) :: torque_rms_
      !> Maximum torque component (hartree/bohr)
      real(wp), intent(in) :: torque_max_
      !> Rotation convergence flag
      logical, intent(in) :: rot_converged_
      !> Norm of the effective translation force (hartree/bohr)
      real(wp), intent(in) :: force_norm_
      !> Effective mobility of the translation (bohr^2/hartree)
      real(wp), intent(in) :: mobility_
      !> Curvature along the effective force (hartree/bohr^2)
      real(wp), intent(in) :: kappa_

      record%istep = istep
      record%energy = energy_
      record%delta_e = de_
      record%gnorm = vec_norm(g_)
      record%gmax = vec_maxabs(g_)
      record%dispnorm = vec_norm(step_)
      record%disp_max = vec_maxabs(step_)
      record%method = settings_%method
      record%fallback = .false.
      ! The dimer has no trust radius; the recorded value is the active bound
      ! of the translation (`step_max`), and `bound_kind` makes the step table
      ! label that column `stepmax/bohr` instead of `trust/bohr`.
      record%trust_radius = settings_%step_max
      record%bound_kind = p_landscape_bound_step_max
      call landscape_result_add_step(result, record)
      call landscape_report_add(result%conv, record)
      if (settings_%printlevel >= 1) then
         call landscape_print_step(record, env%unit)
         write(env%unit, '(1x,a,i0,a,es13.6,a,es13.6,a,l1)') &
            & "  dimer rotation: iterations ", n_rot_, ", torque rms ", &
            & torque_rms_, " max ", torque_max_, " Eh/bohr, converged=", &
            & rot_converged_
         write(env%unit, '(1x,a,es15.6,a)') "  dimer curvature: ", curvature_, &
            & " Eh/bohr^2 (Eq. 2)"
         write(env%unit, '(1x,a,es13.6,a,es13.6,a,es13.6,a,es13.6,a)') &
            & "  dimer translation: |F_dagger| ", force_norm_, &
            & " Eh/bohr, step ", vec_norm(step_), " bohr, mobility ", &
            & mobility_, " (kappa ", kappa_, " Eh/bohr^2)"
      end if
   end subroutine record_step

end subroutine landscape_dimer_optimize

!> Forces (`-gradient`) at the two images of a dimer, i.e. at `R +/- delta*N`.
subroutine dimer_image_forces(objective, x, axis, delta, f_plus, f_minus, fail)
   !> Objective function
   class(landscape_energy_model), intent(inout) :: objective
   !> Midpoint geometry (bohr)
   real(wp), intent(in) :: x(:)
   !> Unit dimer axis
   real(wp), intent(in) :: axis(:)
   !> Centre-to-image distance (bohr)
   real(wp), intent(in) :: delta
   !> Force at `R + delta*N` (hartree/bohr)
   real(wp), intent(out) :: f_plus(:)
   !> Force at `R - delta*N` (hartree/bohr)
   real(wp), intent(out) :: f_minus(:)
   !> Failure flag
   logical, intent(out) :: fail

   real(wp), allocatable :: gradient(:)
   real(wp) :: e_tmp
   integer :: n

   fail = .false.
   n = size(x)
   allocate(gradient(n))

   call objective%evaluate(x + delta*axis, e_tmp, gradient, fail)
   if (fail .or. .not. landscape_value_present(e_tmp) .or. &
         & .not. all(landscape_value_present(gradient))) then
      fail = .true.
      return
   end if
   f_plus = -gradient

   call objective%evaluate(x - delta*axis, e_tmp, gradient, fail)
   if (fail .or. .not. landscape_value_present(e_tmp) .or. &
         & .not. all(landscape_value_present(gradient))) then
      fail = .true.
      return
   end if
   f_minus = -gradient

end subroutine dimer_image_forces

!> True if every dimer length scale is usable.
pure function dimer_settings_valid(settings) result(valid)
   !> Transition state settings
   type(landscape_ts_settings), intent(in) :: settings
   !> True if the parameters are positive
   logical :: valid

   valid = settings%delta > 0.0_wp .and. settings%dR > 0.0_wp .and. &
      & settings%rot_dtheta > 0.0_wp .and. settings%step0 > 0.0_wp .and. &
      & settings%step_max > 0.0_wp

end function dimer_settings_valid

!> Normalise an axis in place; a zero vector is left untouched.
pure subroutine normalise_axis(axis)
   !> Axis
   real(wp), intent(inout) :: axis(:)

   real(wp) :: nrm

   nrm = vec_norm(axis)
   if (nrm > dimer_tiny) axis = axis/nrm

end subroutine normalise_axis

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

end module xtb_landscape_saddle_dimer
