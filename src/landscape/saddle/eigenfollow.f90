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

!> Eigenvector following and Cartesian coordinate hygiene for transition
!> state searches (PLAN section 6.3).
!>
!> Two independent pieces live here because both are pure linear algebra on a
!> Cartesian displacement vector and are shared by the P-RFO and B-P-RFO
!> kernels:
!>
!>  * **Mode selection and following.**  The transition vector is the
!>    eigenvector of the Hessian along which the PES is *maximised*.  It is
!>    selected once (the most negative mode, or a target-corrected one) and
!>    then followed across the macro iterations.  Ordinary eigenvalue sorting
!>    reorders the eigenvectors whenever two eigenvalues cross, so the
!>    following step interposes an explicit overlap assignment
!>    `idx = argmax_i |V_i^T v|` together with a sign realignment
!>    `v_new = sign(V_idx^T v) V_idx` (pysisyphus `PRFO.py:634-672`,
!>    PLAN section 6.3).  An overlap below `mode_overlap_min` is treated as a
!>    lost transition vector: the mode is reselected and the caller is warned.
!>
!>  * **Rigid body projection.**  The Cartesian Hessian of a nonlinear
!>    molecule has six exactly (or, for a numerical Hessian, numerically)
!>    vanishing directions - the three translations and the three rotations -
!>    that are independent of the atomic masses.  A trust region step that
!>    keeps a component along them spends part of its radius on a motion that
!>    cannot change the energy; the M1 RFO kernel did exactly that, because
!>    the lifted eigenvalues (`max(b_i, hlow)`) leave a small residual force
!>    along the rigid directions.  `landscape_rigid_modes` builds an
!>    orthonormal basis of that six dimensional space and
!>    `landscape_project_rigid` removes the component of a step along it.
!>
!> Conventions:
!>
!>  * `mctc_syev` returns eigenvalues in **ascending** order and the
!>    eigenvectors as the columns of the overwritten matrix
!>    (`xtb_mctc_lapack`).  "Most negative mode" therefore means index 1.
!>  * the overall sign of an eigenvector is arbitrary.  The phase convention
!>    of a freshly selected transition vector fixes it: the component with the
!>    largest absolute value is made positive and reported together with its
!>    atom number (PLAN section 6.3, R04 section 35).  Following keeps the
!>    sign continuous, it does not re-apply the phase convention.
!>  * flattened geometries are `xyz(3, n)` (`x1,y1,z1,x2,...`) and displacement
!>    vectors have length `3n` in the same order.
!>
!> Units: coordinates in bohr, Hessian elements in hartree/bohr^2,
!> eigenvalues in hartree/bohr^2.
module xtb_landscape_saddle_eigenfollow
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_lapack, only : mctc_syev
   use xtb_type_environment, only : TEnvironment
   use xtb_landscape_coords, only : vec_norm
   use xtb_landscape_settings, only : landscape_log_failed, landscape_real_to_string
   implicit none(type, external)
   private

   public :: landscape_transition_mode
   public :: landscape_spectrum
   public :: landscape_select_transition_mode
   public :: landscape_follow_transition_mode
   public :: landscape_transition_mode_report
   public :: landscape_transition_mode_set_phase
   public :: landscape_rigid_modes
   public :: landscape_project_rigid
   public :: landscape_project_hessian_rigid
   public :: landscape_curvature_noise_tol
   public :: landscape_rfo_downhill_shift, landscape_rfo_downhill_residual

   !> Norm below which a candidate rigid body vector is considered linearly
   !> dependent on the ones already accepted (bohr, relative to unit vectors)
   real(wp), parameter :: rigid_dependent_tol = 1.0e-8_wp

   !> Relative threshold below which a Hessian eigenvalue is numerical noise
   !> rather than a curvature: the negative eigenvalues of a numerical Hessian
   !> whose rigid body directions were not projected away are of the order
   !> `eps * ||H||` (measured ~1e-11 for a GFN2 molecule), whereas a physical
   !> curvature is never that small relative to the spectrum.  An eigenvalue
   !> above `-tol` therefore does not count as negative curvature and does not
   !> trigger the `hlow` lift.
   real(wp), parameter :: curvature_noise_relative = 1.0e-8_wp

   !> The followed transition vector and its provenance.
   !>
   !> `vector` is always a unit vector in Cartesian displacement space
   !> (length 3N).  `index` is its position in the ascending eigenvalue list of
   !> the Hessian it was last assigned from, `eigenvalue` that eigenvalue and
   !> `overlap` the absolute overlap the assignment was based on
   !> (`overlap = 1` for a fresh selection).  `followed` is false whenever the
   !> vector was (re)selected instead of tracked, so a run that lost the
   !> transition vector is visible in the result.
   type :: landscape_transition_mode
      !> Unit transition vector (Cartesian, length 3N)
      real(wp), allocatable :: vector(:)
      !> 1-based index in the ascending eigenvalue list
      integer :: index = 0
      !> Eigenvalue of the mode (hartree/bohr^2)
      real(wp) :: eigenvalue = 0.0_wp
      !> Absolute overlap of the last following step
      real(wp) :: overlap = 1.0_wp
      !> True if the last update tracked the vector, false if it reselected it
      logical :: followed = .false.
      !> Target negative mode requested by the settings (1-based)
      integer :: target = 1
      !> Flat index of the largest absolute component
      integer :: phase_index = 0
      !> Value of the largest absolute component after the phase convention
      real(wp) :: phase_value = 0.0_wp
   end type landscape_transition_mode

contains

!> Symmetric eigenvalue decomposition of a Hessian.
!>
!> The matrix is symmetrised on the way in (the reduced Hessian of an update
!> formula may carry round off asymmetry) and passed to `mctc_syev`, which
!> returns the eigenvalues in ascending order and the eigenvectors as columns.
!> A failure of the solver is reported through `fail`; the environment is
!> probed only for the messages appended by this call.
subroutine landscape_spectrum(env, hess, evals, evecs, fail)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Hessian (hartree/bohr^2), only read
   real(wp), intent(in) :: hess(:, :)
   !> Ascending eigenvalues (hartree/bohr^2)
   real(wp), allocatable, intent(out) :: evals(:)
   !> Eigenvectors as columns (dimensionless)
   real(wp), allocatable, intent(out) :: evecs(:, :)
   !> Failure flag
   logical, intent(out) :: fail

   real(wp), allocatable :: work(:, :)
   integer :: n, nlog

   fail = .false.
   n = size(hess, 1)
   if (n == 0 .or. size(hess, 2) /= n) then
      call env%error("the Hessian of the spectrum must be a square matrix", &
         & 'landscape_spectrum')
      fail = .true.
      return
   end if

   allocate(work(n, n))
   work = 0.5_wp*(hess + transpose(hess))
   allocate(evals(n))
   nlog = env%nLog
   call mctc_syev(env, work, evals)
   if (landscape_log_failed(env, nlog)) then
      fail = .true.
      return
   end if
   call move_alloc(work, evecs)

end subroutine landscape_spectrum

!> Select the initial transition vector from a Hessian spectrum.
!>
!> `target` selects which negative mode is followed (1 = the most negative,
!> matching the `tsroot` semantics).  A Hessian without that many negative
!> eigenvalues is not an error - at a minimum the search has to start
!> uphill anyway - but it is reported: the lowest **internal** mode is
!> followed instead and a warning names the mismatch.
!>
!> Only eigenvalues below `-landscape_curvature_noise_tol` count as negative
!> curvature, and the fallback skips every eigenvalue inside `+/-tol`.  This is
!> the MEDIUM-1 repair of `31-m4-review-round2.md`: the rigid body directions
!> of a molecule are lifted to `eps * ||H||` by the numerical Hessian (and to a
!> few `1e-17` hartree/bohr^2 by the `P H P` projection), so an ascending
!> spectrum starts with the six (five for a linear molecule) numerically zero
!> modes.  Selecting index 1 of a **minimum** therefore picked a rigid body
!> zero mode whose ascent direction is removed again by the step projection,
!> and the search degenerated into a minimisation - the documented "the lowest
!> mode is followed" did not hold.  With the tolerance the fallback selects the
!> lowest mode that is neither rigid nor numerically unresolved.
subroutine landscape_select_transition_mode(env, evals, evecs, target, mode, fail, unit)
   character(len=*), parameter :: source = 'landscape_select_transition_mode'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Ascending eigenvalues (hartree/bohr^2)
   real(wp), intent(in) :: evals(:)
   !> Eigenvectors as columns
   real(wp), intent(in) :: evecs(:, :)
   !> Requested negative mode (1 = most negative)
   integer, intent(in) :: target
   !> Transition mode, overwritten
   type(landscape_transition_mode), intent(inout) :: mode
   !> Failure flag
   logical, intent(out) :: fail
   !> Output unit for the phase report, optional
   integer, intent(in), optional :: unit

   real(wp) :: tol
   integer :: n, nneg, index

   fail = .false.
   n = size(evals)
   if (n == 0 .or. size(evecs, 1) /= n) then
      call env%error("the spectrum of the transition mode is empty", source)
      fail = .true.
      return
   end if

   tol = landscape_curvature_noise_tol(evals)
   nneg = count(evals < -tol)
   if (target >= 1 .and. target <= nneg) then
      ! ascending order puts the significant negative eigenvalues first, so the
      ! target-th most negative mode is at index `target`
      index = target
   else if (nneg >= 1) then
      ! fewer negative eigenvalues than requested: follow the most negative one
      ! (still the first entry, because every near-zero value is larger than a
      ! significant negative eigenvalue in an ascending spectrum)
      index = 1
      call env%warning("the Hessian has fewer negative eigenvalues than "//&
         & "the requested target mode; the most negative mode is followed", &
         & source)
   else
      ! No negative eigenvalue at all: follow the lowest mode that is neither a
      ! rigid body zero mode nor numerically unresolved.
      index = landscape_lowest_internal_mode(evals, tol)
      call env%warning("the Hessian has no negative eigenvalue; the "//&
         & "lowest internal mode is followed as the initial transition vector", &
         & source)
   end if

   if (allocated(mode%vector)) deallocate(mode%vector)
   allocate(mode%vector(n))
   mode%vector = evecs(:, index)
   call normalise(mode%vector)
   mode%index = index
   mode%eigenvalue = evals(index)
   mode%overlap = 1.0_wp
   mode%followed = .false.
   mode%target = max(1, target)
   call landscape_transition_mode_set_phase(mode)
   if (present(unit)) call landscape_transition_mode_report(mode, unit)

end subroutine landscape_select_transition_mode

!> Follow the transition vector through a new Hessian spectrum.
!>
!> The new vector is the eigenvector with the largest absolute overlap with the
!> tracked one, realigned so that the overlap is positive.  This is the
!> standard cure for the arbitrary reordering (and sign) of a fresh
!> diagonalisation.  An absolute overlap below `overlap_min` means the vector
!> was lost between the two Hessians: the mode is reselected through
!> `landscape_select_transition_mode` and a warning is emitted.  In that case
!> `mode%followed` is false and `mode%overlap` holds the failing overlap.
subroutine landscape_follow_transition_mode(env, evals, evecs, mode, overlap_min, &
      & fail, unit)
   character(len=*), parameter :: source = 'landscape_follow_transition_mode'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Ascending eigenvalues (hartree/bohr^2)
   real(wp), intent(in) :: evals(:)
   !> Eigenvectors as columns
   real(wp), intent(in) :: evecs(:, :)
   !> Transition mode, updated in place
   type(landscape_transition_mode), intent(inout) :: mode
   !> Overlap below which the mode is reselected
   real(wp), intent(in) :: overlap_min
   !> Failure flag
   logical, intent(out) :: fail
   !> Output unit for the phase report, optional
   integer, intent(in), optional :: unit

   real(wp), allocatable :: overlaps(:)
   real(wp) :: overlap
   integer :: i, n, index, target

   fail = .false.
   n = size(evals)
   if (n == 0 .or. size(evecs, 1) /= n) then
      call env%error("the spectrum of the transition mode is empty", source)
      fail = .true.
      return
   end if
   if (.not. allocated(mode%vector) .or. size(mode%vector) /= n) then
      call env%error("the transition mode has not been selected yet", source)
      fail = .true.
      return
   end if

   allocate(overlaps(n))
   do i = 1, n
      overlaps(i) = dot_product(evecs(:, i), mode%vector)
   end do
   index = maxloc(abs(overlaps), dim=1)
   overlap = abs(overlaps(index))

   if (overlap < overlap_min) then
      call env%warning("mode following lost the transition vector (overlap "//&
         & trim(landscape_real_to_string(overlap, '(f0.4)'))//" below the threshold "//&
         & trim(landscape_real_to_string(overlap_min, '(f0.4)'))//"); reselecting the mode", source)
      target = mode%target
      call landscape_select_transition_mode(env, evals, evecs, target, mode, fail)
      mode%overlap = overlap
      mode%followed = .false.
      if (present(unit) .and. .not. fail) call landscape_transition_mode_report(mode, unit)
      return
   end if

   ! sign realignment: the new eigenvector is oriented so that its dot product
   ! with the tracked vector stays positive (`v_new = sign(overlap) V_index`),
   ! so a global sign flip of the diagonalisation cannot invert the ascent
   ! direction.
   if (overlaps(index) < 0.0_wp) then
      mode%vector = -evecs(:, index)
   else
      mode%vector = evecs(:, index)
   end if
   mode%index = index
   mode%eigenvalue = evals(index)
   mode%overlap = overlap
   mode%followed = .true.

   ! The phase convention is deliberately **not** re-applied here: it belongs
   ! to the initial selection (it defines the forward direction of the report),
   ! and re-applying it would undo the sign realignment above whenever the
   ! largest component of the new eigenvector happens to be negative.
   if (present(unit)) call landscape_transition_mode_report(mode, unit)

end subroutine landscape_follow_transition_mode

!> Print the followed mode: index, eigenvalue, largest component and overlap.
!>
!> The printed component is the current largest component, not the phase marker
!> of the initial selection; the forward direction and its atom number are
!> defined (and printed) once by `landscape_select_transition_mode`.
subroutine landscape_transition_mode_report(mode, unit)
   !> Transition mode
   type(landscape_transition_mode), intent(in) :: mode
   !> Output unit
   integer, intent(in) :: unit

   integer :: iatom, icart, imax

   write(unit, '(1x,a,i0,a,f16.8,a)') "transition mode: index ", mode%index, &
      & ", eigenvalue ", mode%eigenvalue, " Eh/bohr^2"
   if (allocated(mode%vector)) then
      imax = maxloc(abs(mode%vector), dim=1)
      iatom = (imax - 1)/3 + 1
      icart = mod(imax - 1, 3) + 1
      if (mod(size(mode%vector), 3) == 0) then
         write(unit, '(1x,a,f12.6,a,i0,a,i0,a)') "transition vector component: ", &
            & mode%vector(imax), " (cartesian ", icart, " of atom ", iatom, ")"
      else
         write(unit, '(1x,a,f12.6,a,i0,a)') "transition vector component: ", &
            & mode%vector(imax), " (component ", imax, ")"
      end if
   end if
   write(unit, '(1x,a,f12.6)') "transition vector overlap: ", mode%overlap

end subroutine landscape_transition_mode_report

!> Apply the phase convention: the largest absolute component becomes positive.
!>
!> The component index and the (positive) value are recorded, so the caller can
!> report the direction that defines "forward" (PLAN section 6.3).  The
!> convention is shared by the eigenvector following methods (applied once, at
!> selection) and by the dimer (applied to the final axis), so that the same
!> geometry always produces the same printed transition vector.
subroutine landscape_transition_mode_set_phase(mode)
   !> Transition mode, vector updated in place
   type(landscape_transition_mode), intent(inout) :: mode

   integer :: imax

   if (.not. allocated(mode%vector)) return
   imax = maxloc(abs(mode%vector), dim=1)
   if (mode%vector(imax) < 0.0_wp) mode%vector = -mode%vector
   mode%phase_index = imax
   mode%phase_value = mode%vector(imax)

end subroutine landscape_transition_mode_set_phase

!> Solve the RFO secular equation `sum_{i/=k} F_i^2/(lambda - b_i) = lambda`.
!>
!> This is Eq. 4.18b of PLAN section 6.3, i.e. the shift `lambda_n < 0` that
!> keeps the complement of the followed mode a descent subspace.  It is shared
!> by the double shift P-RFO and the restricted step B-P-RFO (the latter falls
!> back to it whenever the downhill step fits its sub-radius).
!>
!> `residual(lambda) = sum F_i^2/(lambda - b_i) - lambda` is strictly
!> decreasing (its derivative is `-sum F_i^2/(lambda - b_i)^2 - 1`), tends to
!> `+inf` for `lambda -> -inf` and to `-inf` when `lambda` approaches
!> `min_{i/=k} b_i` from below, so exactly one root exists below that
!> eigenvalue, and that root is negative because `residual(0) < 0` for positive
!> downhill eigenvalues.  The bracket is opened by doubling and then bisected.
subroutine landscape_rfo_downhill_shift(env, f, b, k, margin, max_bisect, lambda, &
      & solved, fail, source)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Gradient in the eigenbasis `F_i = V_i^T g`
   real(wp), intent(in) :: f(:)
   !> (Lifted) eigenvalues (hartree/bohr^2)
   real(wp), intent(in) :: b(:)
   !> Index of the followed mode
   integer, intent(in) :: k
   !> Margin below the smallest downhill eigenvalue
   real(wp), intent(in) :: margin
   !> Maximum number of bisections
   integer, intent(in) :: max_bisect
   !> Root of the secular equation
   real(wp), intent(out) :: lambda
   !> True if the secular equation had a non-empty downhill subspace
   logical, intent(out) :: solved
   !> Failure flag
   logical, intent(out) :: fail
   !> Name of the calling routine for the diagnostics
   character(len=*), intent(in) :: source

   integer, parameter :: max_bracket = 200
   real(wp) :: bmin, lo, hi, mid, gap, fl
   integer :: i, n, nbr, ibis

   fail = .false.
   solved = .false.
   lambda = 0.0_wp
   n = size(b)

   bmin = huge(1.0_wp)
   do i = 1, n
      if (i == k) cycle
      bmin = min(bmin, b(i))
   end do
   if (bmin >= huge(1.0_wp)) return

   ! an all-zero downhill gradient has no downhill step at all
   if (sum(f**2) - f(k)**2 <= tiny(1.0_wp)) then
      lambda = 0.0_wp
      return
   end if
   solved = .true.

   hi = bmin - max(margin, 1.0e-14_wp*max(1.0_wp, abs(bmin)))
   gap = max(1.0_wp, abs(hi))
   lo = min(hi, 0.0_wp) - gap
   nbr = 0
   do while (landscape_rfo_downhill_residual(lo, f, b, k) <= 0.0_wp .and. &
         & nbr < max_bracket)
      lo = lo - gap
      nbr = nbr + 1
   end do
   fl = landscape_rfo_downhill_residual(lo, f, b, k)
   if (fl <= 0.0_wp) then
      call env%warning("the lambda_n bracket did not enclose the root", source)
      fail = .true.
      return
   end if

   ! The upper end of the bracket has to be on the other side of the root as
   ! well: `residual` is strictly decreasing, so `residual(hi) > 0` means that
   ! its zero lies *above* `hi` - outside the branch the shift is defined on -
   ! and bisection would silently return `hi` instead of a root (MEDIUM-2 of
   ! `31-m4-review-round2.md`).  The root is missing when every downhill
   ! gradient component is tiny,
   ! `sum_{i/=k} f_i^2 <= margin * (margin - bmin)`, i.e. `|f_i| <~ 1e-8` at the
   ! default `mu_margin = 1e-8`, while a downhill eigenvalue is negative but
   ! below the curvature noise threshold (so the `hlow` lift did not engage and
   ! `bmin` itself is that near-zero value).  Measured on the M4 searches
   ! (GFN2 HCN from the bent guess and from the escaped minimum, both methods)
   ! the residual at `hi` is `-5e-4 ... -2e-2` and the residual of the returned
   ! root is `|1e-16|`, i.e. eleven orders of magnitude away from the failure
   ! mode; the numerical `bmin ~ -2e-16` of the projected Hessian is harmless.
   ! If the inequality ever breaks, the equation has no downhill branch at all
   ! and the honest answer is the degenerate `lambda_n = 0` (a plain Newton
   ! descent of the downhill subspace), reported as such instead of returning a
   ! shift that is not a root.
   if (landscape_rfo_downhill_residual(hi, f, b, k) > 0.0_wp) then
      call env%warning("the lambda_n secular equation has no root below the "//&
         & "smallest downhill eigenvalue; the downhill step uses lambda_n = 0", &
         & source)
      lambda = 0.0_wp
      solved = .false.
      return
   end if

   do ibis = 1, max_bisect
      mid = 0.5_wp*(lo + hi)
      if (landscape_rfo_downhill_residual(mid, f, b, k) > 0.0_wp) then
         lo = mid
      else
         hi = mid
      end if
   end do
   lambda = 0.5_wp*(lo + hi)

end subroutine landscape_rfo_downhill_shift

!> Left hand side of the RFO secular equation minus its right hand side.
pure function landscape_rfo_downhill_residual(lambda, f, b, k) result(value)
   !> Trial shift
   real(wp), intent(in) :: lambda
   !> Gradient in the eigenbasis
   real(wp), intent(in) :: f(:)
   !> (Lifted) eigenvalues
   real(wp), intent(in) :: b(:)
   !> Index of the followed mode
   integer, intent(in) :: k
   !> Residual
   real(wp) :: value

   integer :: i

   value = -lambda
   do i = 1, size(b)
      if (i == k) cycle
      value = value + f(i)**2/(lambda - b(i))
   end do

end function landscape_rfo_downhill_residual

!> Build an orthonormal basis of the six rigid body directions of a geometry.
!>
!> The basis is mass independent: the translations are the unit vectors of the
!> three Cartesian axes repeated on every atom, the rotations are
!> `e_a x (x_i - x_c)` about the geometric centre `x_c`.  Both are proper null
!> vectors of the Cartesian Hessian for any geometry, so removing them from a
!> step removes motion that cannot change the energy.  The rotations are
!> Gram-Schmidt orthogonalised against the translations and each other; for a
!> linear molecule one of them is dependent (only two rotational directions
!> exist) and is dropped, leaving `n_rigid = 5`.
!>
!> `basis` is allocated with six columns, the first `n_rigid` of which are
!> orthonormal.  A geometry whose length is not a multiple of three has no
!> rigid body interpretation and yields `n_rigid = 0`.
subroutine landscape_rigid_modes(x, basis, n_rigid)
   !> Flattened coordinates (bohr, length 3N)
   real(wp), intent(in) :: x(:)
   !> Orthonormal rigid body basis (length 3N, six columns)
   real(wp), allocatable, intent(out) :: basis(:, :)
   !> Number of valid columns (0, 3, 5 or 6)
   integer, intent(out) :: n_rigid

   real(wp), allocatable :: raw(:, :), col(:)
   real(wp) :: centre(3), rvec(3), proj, nrm
   integer :: nat, dim, i, c, a, k

   n_rigid = 0
   dim = size(x)
   if (mod(dim, 3) /= 0) then
      allocate(basis(dim, 0))
      return
   end if
   nat = dim/3
   allocate(basis(dim, 6), source=0.0_wp)
   if (nat <= 1) return
   ! Two translations are the only rigid motions of a single atom; even that is
   ! degenerate, so an isolated atom gets an empty basis.

   centre = 0.0_wp
   do i = 1, nat
      centre = centre + x(3*(i - 1) + 1:3*(i - 1) + 3)
   end do
   centre = centre/real(nat, wp)

   allocate(raw(dim, 6), source=0.0_wp)
   do c = 1, 3
      do i = 1, nat
         raw(3*(i - 1) + c, c) = 1.0_wp
      end do
   end do
   ! rotations: column 3 + c holds (e_c x (x_i - centre))
   do i = 1, nat
      rvec = x(3*(i - 1) + 1:3*(i - 1) + 3) - centre
      raw(3*(i - 1) + 1, 4) = 0.0_wp
      raw(3*(i - 1) + 2, 4) = -rvec(3)
      raw(3*(i - 1) + 3, 4) = rvec(2)
      raw(3*(i - 1) + 1, 5) = rvec(3)
      raw(3*(i - 1) + 2, 5) = 0.0_wp
      raw(3*(i - 1) + 3, 5) = -rvec(1)
      raw(3*(i - 1) + 1, 6) = -rvec(2)
      raw(3*(i - 1) + 2, 6) = rvec(1)
      raw(3*(i - 1) + 3, 6) = 0.0_wp
   end do

   allocate(col(dim))
   do k = 1, 6
      col = raw(:, k)
      do a = 1, n_rigid
         proj = dot_product(basis(:, a), col)
         col = col - proj*basis(:, a)
      end do
      nrm = vec_norm(col)
      if (nrm <= rigid_dependent_tol) cycle
      n_rigid = n_rigid + 1
      basis(:, n_rigid) = col/nrm
   end do

end subroutine landscape_rigid_modes

!> Remove the component of `v` along the (normalised) rigid body basis.
!>
!> `v` is modified in place; the result is orthogonal to every basis vector up
!> to round off.  No scaling is applied, so a projected trust region step keeps
!> its original radius as an upper bound.
pure subroutine landscape_project_rigid(v, basis, n_rigid)
   !> Displacement vector (bohr), projected in place
   real(wp), intent(inout) :: v(:)
   !> Orthonormal rigid body basis
   real(wp), intent(in) :: basis(:, :)
   !> Number of valid basis columns
   integer, intent(in) :: n_rigid

   real(wp) :: proj
   integer :: a

   do a = 1, n_rigid
      proj = dot_product(basis(:, a), v)
      v = v - proj*basis(:, a)
   end do

end subroutine landscape_project_rigid

!> Remove the rigid body components from a Hessian, `H <- P H P` with the
!> orthogonal projector `P = I - sum_a b_a b_a^T` of `landscape_rigid_modes`.
!>
!> The step built by the double shift methods is projected onto the complement
!> of the rigid body space, but the *spectrum* the step is built from was not:
!> the numerical Hessian of a molecule lifts its six null directions to
!> `eps * ||H||`, and those tiny negative eigenvalues were counted as downhill
!> negative curvature, so the `hlow` lift engaged on every step of a clean
!> run.  Projecting the Hessian as well removes the rigid body block before the
!> eigenvalues are counted, so `lifted` reflects genuine negative curvature.
!> The projection is a similarity transform of the symmetric matrix and keeps
!> it symmetric (up to round off).
subroutine landscape_project_hessian_rigid(hess, basis, n_rigid)
   !> Hessian (hartree/bohr^2), projected in place
   real(wp), intent(inout) :: hess(:, :)
   !> Orthonormal rigid body basis (from `landscape_rigid_modes`)
   real(wp), intent(in) :: basis(:, :)
   !> Number of valid basis columns
   integer, intent(in) :: n_rigid

   integer :: j

   if (n_rigid <= 0) return
   ! P H : project every column
   do j = 1, size(hess, 2)
      call landscape_project_rigid(hess(:, j), basis, n_rigid)
   end do
   ! (P H) P^T = P H P : project every row (P is symmetric)
   do j = 1, size(hess, 1)
      call landscape_project_rigid(hess(j, :), basis, n_rigid)
   end do

end subroutine landscape_project_hessian_rigid

!> Index of the lowest Hessian mode that is not numerically zero.
!>
!> The ascending spectrum of a molecular Hessian starts with the rigid body
!> directions, which a projected Hessian leaves at `O(eps ||H||)` and an
!> unprojected numerical Hessian at `O(eps ||H||)` as well.  Those modes carry
!> no internal motion, so they must never be followed as a transition vector
!> (MEDIUM-1 of `31-m4-review-round2.md`).  The function returns the first
!> index with `evals(i) > tol`; when every eigenvalue is inside `+/-tol` the
!> numerically most significant mode is returned instead, so the caller always
!> gets a usable index.
pure function landscape_lowest_internal_mode(evals, tol) result(index)
   !> Ascending Hessian eigenvalues (hartree/bohr^2)
   real(wp), intent(in) :: evals(:)
   !> Noise threshold (hartree/bohr^2)
   real(wp), intent(in) :: tol
   !> 1-based index of the selected mode
   integer :: index

   integer :: i

   index = 0
   do i = 1, size(evals)
      if (evals(i) > tol) then
         index = i
         return
      end if
   end do
   if (size(evals) > 0) index = maxloc(abs(evals), dim=1)

end function landscape_lowest_internal_mode

!> Magnitude below which a Hessian eigenvalue counts as noise (hartree/bohr^2).
!>
!> The threshold is relative to the largest absolute eigenvalue of the
!> spectrum, so it scales with the system exactly like the round off of the
!> numerical Hessian it is meant to filter.  Eigenvalues above `-tol` are not
!> negative curvature.
pure function landscape_curvature_noise_tol(evals) result(tol)
   !> Ascending Hessian eigenvalues (hartree/bohr^2)
   real(wp), intent(in) :: evals(:)
   !> Noise threshold (hartree/bohr^2)
   real(wp) :: tol

   real(wp) :: scale

   scale = 0.0_wp
   if (size(evals) > 0) scale = maxval(abs(evals))
   tol = curvature_noise_relative*max(scale, tiny(1.0_wp))

end function landscape_curvature_noise_tol

!> Normalise a vector in place; a zero vector is left untouched.
pure subroutine normalise(v)
   !> Vector to normalise
   real(wp), intent(inout) :: v(:)

   real(wp) :: nrm

   nrm = vec_norm(v)
   if (nrm > tiny(1.0_wp)) v = v/nrm

end subroutine normalise

end module xtb_landscape_saddle_eigenfollow
