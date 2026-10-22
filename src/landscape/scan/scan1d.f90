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

!> One-dimensional scan driver of the landscape suite.
!>
!> The scan is anchored at the coordinate value `q0` of the input geometry.
!> A **forward** sweep visits `q0 + k*step` for `k = 0..nstep`; with
!> `bidirectional` a **reverse** sweep follows immediately, visiting the same
!> values in decreasing order (`k = nstep-1, ..., 0`) without returning to the
!> starting geometry in between.  The two branches therefore form a closed
!> cycle over the same grid, and the point counts are
!>
!>     nstep + 1                    (unidirectional)
!>     2*nstep + 1                  (bidirectional, turnaround point shared)
!>
!> The difference `E_forward(q) - E_reverse(q)` on the shared grid is the
!> hysteresis of the scan.  It is non-zero whenever the relaxed structure is
!> history dependent (two basins separated by a barrier along the scan
!> coordinate), which is exactly the dominant failure mode that R-04 documents
!> for relaxed torsional scans
!> (`docs/_landscape-recon/04-web-research-ts-irc-state-of-art.md` section 4.6)
!> and what experiment 11 of PLAN section 8.2 asks to quantify.  The report
!> contains the maximum difference together with the grid point at which it
!> occurs and the largest energy jump between consecutive points (the
!> discontinuity of the profile).
!>
!> `warm_start` passes the geometry of the previous point as the starting
!> geometry of the next one.  Only the **geometry** is carried over; the
!> wavefunction policy stays entirely with the engine restart mode
!> (`landscape_engine_control%restart_mode`), because mixing warm and cold
!> restarts inside one task would introduce a SCF-history dependence of about
!> 4e-7 Eh/bohr (`docs/_landscape-recon/07-*-findings.md`).
!>
!> Units: coordinates in bohr, energies in hartree, distances in bohr, angles
!> in radian, gradients in hartree/bohr (or hartree/rad for an angle).
module xtb_landscape_scan1d
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : pi
   use xtb_mctc_convert, only : autoaa
   use xtb_type_environment, only : TEnvironment
   use xtb_landscape_settings, only : landscape_scan_settings, &
      & landscape_scan_coord_settings, p_landscape_scan_coord_distance, &
      & landscape_missing_value, landscape_lowercase, landscape_real_to_string
   use xtb_landscape_engine, only : landscape_energy_model
   use xtb_landscape_coords, only : geometry_unpack
   use xtb_landscape_scan_relax, only : landscape_scan_point, &
      & landscape_scan_relax_point, landscape_scan_rigid_point, &
      & landscape_scan_coord_value, landscape_scan_coord_wrap_delta
   use xtb_landscape_results, only : landscape_result
   implicit none(type, external)
   private

   public :: landscape_scan_hysteresis, landscape_scan_profile
   public :: landscape_scan1d, landscape_scan_solve_point
   public :: landscape_scan_profile_write, scan_jump_warning

   !> Energy jump between two consecutive scan points above which the driver
   !> reports a warning (hartree).  A jump of 0.5 hartree is far beyond any
   !> smooth change of a GFN2 surface per scan step and indicates a basin
   !> switch or a SCF discontinuity (PLAN section 6.2).
   real(wp), parameter :: scan_jump_warning = 0.5_wp

   !> Hysteresis report of a bidirectional scan.
   !>
   !> Points without an energy (NaN) take no part in the comparison: `NaN`
   !> compares false against every bound, so a pair that involves a missing
   !> energy is skipped instead of poisoning the maximum (the profile is
   !> broken there anyway, which `n_fail` reports).
   type :: landscape_scan_hysteresis
      !> Number of grid points compared between the two branches
      integer :: n_compared = 0
      !> Largest `|E_forward - E_reverse|` (hartree)
      real(wp) :: max_delta = 0.0_wp
      !> Coordinate value at which `max_delta` occurs (bohr or radian)
      real(wp) :: at_value = 0.0_wp
      !> Largest energy jump between consecutive points of the profile
      !> (hartree), i.e. the discontinuity of the scan curve
      real(wp) :: max_jump = 0.0_wp
      !> Coordinate value before the largest jump (bohr or radian)
      real(wp) :: jump_from = 0.0_wp
      !> Coordinate value after the largest jump (bohr or radian)
      real(wp) :: jump_to = 0.0_wp
   end type landscape_scan_hysteresis

   !> Result of a one-dimensional scan.
   type :: landscape_scan_profile
      !> Number of stored scan points
      integer :: n_points = 0
      !> Number of points of the forward branch
      integer :: n_forward = 0
      !> True if the profile contains a reverse branch
      logical :: bidirectional = .false.
      !> Scanned coordinate value per point (bohr or radian)
      real(wp), allocatable :: values(:)
      !> Coordinate value reached per point (bohr or radian); a relaxed point
      !> stays offset from the requested value by the finite stiffness of the
      !> constraint bias (PLAN R-04)
      real(wp), allocatable :: reached(:)
      !> De-biased energy per point (hartree), NaN if the point has none
      real(wp), allocatable :: energy(:)
      !> True if the point carries a usable energy; a frame of `xtbscan.log` is
      !> written for these points only
      logical, allocatable :: energy_valid(:)
      !> Biased energy per point (hartree)
      real(wp), allocatable :: total(:)
      !> Constraint residual per point (hartree/bohr or hartree/rad)
      real(wp), allocatable :: residual(:)
      !> Inner minimisation steps per point
      integer, allocatable :: n_steps(:)
      !> Inner convergence flag per point
      logical, allocatable :: converged(:)
      !> Geometry per point (bohr), stored in visit order
      real(wp), allocatable :: xyz(:, :)
      !> Hysteresis report of the scan
      type(landscape_scan_hysteresis) :: hysteresis
      !> Number of points whose inner minimisation failed
      integer :: n_fail = 0
   end type landscape_scan_profile

contains

!> Solve one scan point, relaxed or rigid.
!>
!> This is the single point solver shared by the 1D driver below and the 2D
!> grid (`xtb_landscape_scan2d`).  It constrains `ncoord` coordinates
!> simultaneously, so a 2D grid point satisfies both scanned coordinates at
!> once.
subroutine landscape_scan_solve_point(env, base, at, settings, coords, targets, &
      & ncoord, x0, xout, point, fail)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Potential
   class(landscape_energy_model), intent(inout) :: base
   !> Atom numbers (1-based)
   integer, intent(in) :: at(:)
   !> Scan settings
   type(landscape_scan_settings), intent(in) :: settings
   !> Scanned coordinates
   type(landscape_scan_coord_settings), intent(in) :: coords(:)
   !> Target values of the coordinates (bohr or radian)
   real(wp), intent(in) :: targets(:)
   !> Number of constraints
   integer, intent(in) :: ncoord
   !> Starting geometry (bohr)
   real(wp), intent(in) :: x0(:)
   !> Geometry of the point (bohr)
   real(wp), allocatable, intent(out) :: xout(:)
   !> Scan point record
   type(landscape_scan_point), intent(out) :: point
   !> Failure flag
   logical, intent(out) :: fail

   type(landscape_result) :: result

   if (landscape_lowercase(settings%mode) == 'rigid') then
      call landscape_scan_rigid_point(env, base, coords, targets, ncoord, at, &
         & x0, xout, point, fail)
   else
      call landscape_scan_relax_point(env, base, settings, coords, targets, &
         & ncoord, at, x0, xout, point, result, fail)
   end if

end subroutine landscape_scan_solve_point

!> Run a one-dimensional scan.
!>
!> On return `profile` holds the scan curve in visit order together with the
!> geometries, the constraint residuals and the hysteresis report.  A failing
!> point (failed geometry construction or single point, unconverged inner
!> relaxation) is recorded, counted in `profile%n_fail` and reported with a
!> NaN energy (`profile%energy_valid` is false); the scan continues so that the
!> driver can report the whole curve, and `fail` is set once at least one
!> point failed.  An unconverged but evaluated point keeps its energy: only a
!> discarded evaluation leaves the energy missing.
subroutine landscape_scan1d(env, base, at, settings, x0, profile, fail)
   character(len=*), parameter :: source = 'landscape_scan1d'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Potential
   class(landscape_energy_model), intent(inout) :: base
   !> Atom numbers (1-based)
   integer, intent(in) :: at(:)
   !> Scan settings
   type(landscape_scan_settings), intent(in) :: settings
   !> Starting geometry (bohr)
   real(wp), intent(in) :: x0(:)
   !> Scan result
   type(landscape_scan_profile), intent(out) :: profile
   !> Failure flag
   logical, intent(out) :: fail

   type(landscape_scan_point) :: point
   type(landscape_scan_coord_settings) :: coords(1)
   real(wp), allocatable :: xstart(:), xout(:), xyz(:, :)
   real(wp) :: q0, step, value, targets(1)
   integer :: nstep, k, index, npoints
   logical :: point_fail

   fail = .false.
   profile%bidirectional = settings%bidirectional
   nstep = settings%coord(1)%nstep
   step = settings%coord(1)%step
   coords(1) = settings%coord(1)

   allocate(xyz(3, size(at)))
   call geometry_unpack(x0, xyz)
   q0 = landscape_scan_coord_value(settings%coord(1), xyz)

   profile%n_forward = nstep + 1
   if (settings%bidirectional) then
      npoints = 2*nstep + 1
   else
      npoints = nstep + 1
   end if
   profile%n_points = npoints
   allocate(profile%values(npoints), profile%reached(npoints), &
      & profile%energy(npoints), profile%total(npoints), &
      & profile%residual(npoints))
   allocate(profile%n_steps(npoints), profile%converged(npoints), &
      & profile%energy_valid(npoints))
   allocate(profile%xyz(3*size(at), npoints))
   profile%values = 0.0_wp
   profile%reached = 0.0_wp
   profile%energy = landscape_missing_value()
   profile%total = 0.0_wp
   profile%residual = 0.0_wp
   profile%xyz = 0.0_wp
   profile%n_steps = 0
   profile%converged = .false.
   profile%energy_valid = .false.

   allocate(xstart(size(x0)), xout(size(x0)))
   xstart = x0

   ! Forward branch: k = 0, 1, ..., nstep
   do k = 0, nstep
      index = k + 1
      value = q0 + real(k, wp)*step
      targets(1) = value
      call landscape_scan_solve_point(env, base, at, settings, coords, targets, &
         & 1, xstart, xout, point, point_fail)
      call store_point(profile, index, value, xout, point, point_fail, fail)
      if (settings%warm_start .and. allocated(xout)) xstart = xout
   end do

   ! Reverse branch: k = nstep-1, ..., 0 (the turnaround point belongs to the
   ! forward branch and is not repeated)
   if (settings%bidirectional) then
      do k = nstep - 1, 0, -1
         index = profile%n_forward + (nstep - k)
         value = q0 + real(k, wp)*step
         targets(1) = value
         call landscape_scan_solve_point(env, base, at, settings, coords, &
            & targets, 1, xstart, xout, point, point_fail)
         call store_point(profile, index, value, xout, point, point_fail, fail)
         if (settings%warm_start .and. allocated(xout)) xstart = xout
      end do
   end if

   call hysteresis_report(profile)
   call jump_report(profile%values, profile%energy, profile%n_points, &
      & profile%hysteresis)

   if (profile%hysteresis%max_jump > scan_jump_warning) then
      call env%warning("the scan profile jumps by "//&
         & trim(landscape_real_to_string(profile%hysteresis%max_jump, '(es15.6)'))//&
         & " hartree between two consecutive points; this indicates a basin "//&
         & "switch or a SCF discontinuity", source)
   end if

   if (profile%n_fail > 0) then
      fail = .true.
      if (settings%printlevel > 0) then
         write(env%unit, '(1x,a,i0,a)') "landscape scan: ", profile%n_fail, &
            & " point(s) carry no energy (their single point failed)"
      end if
   end if

end subroutine landscape_scan1d

!> Store one solved point in the profile and update the failure counter.
subroutine store_point(profile, index, value, xout, point, point_fail, fail)
   !> Profile to fill
   type(landscape_scan_profile), intent(inout) :: profile
   !> Index of the point
   integer, intent(in) :: index
   !> Requested coordinate value (bohr or radian)
   real(wp), intent(in) :: value
   !> Geometry of the point (bohr), may be unallocated on failure
   real(wp), allocatable, intent(in) :: xout(:)
   !> Scan point record
   type(landscape_scan_point), intent(in) :: point
   !> True if this point failed
   logical, intent(in) :: point_fail
   !> Overall failure flag
   logical, intent(inout) :: fail

   profile%values(index) = value
   profile%reached(index) = point%values(1)
   profile%total(index) = point%total
   profile%residual(index) = point%residuals(1)
   profile%n_steps(index) = point%n_steps
   profile%converged(index) = point%converged
   profile%energy_valid(index) = point%energy_valid
   if (point%energy_valid) then
      profile%energy(index) = point%energy
   else
      ! Never store the `0.0` of an initialised but never evaluated point: a
      ! failed point has *no* energy (see `landscape_missing_value`).
      profile%energy(index) = landscape_missing_value()
   end if
   if (allocated(xout)) profile%xyz(:, index) = xout
   if (point_fail) then
      profile%n_fail = profile%n_fail + 1
      fail = .true.
   end if

end subroutine store_point

!> Difference between the two branches of a bidirectional scan.
!>
!> The forward branch occupies the indices `1..n_forward` (values `k = 0..n`),
!> the reverse branch the indices `n_forward+1..n_points` (values
!> `k = n-1..0`, i.e. value `k` sits at `n_forward + (n - k)`).  The turnaround
!> point `k = n` is shared by construction and therefore not compared, so a
!> bidirectional run compares exactly `nstep` grid points.
subroutine hysteresis_report(profile)
   !> Scan profile, the hysteresis field is filled
   type(landscape_scan_profile), intent(inout) :: profile

   integer :: k, ifwd, irev

   profile%hysteresis%n_compared = 0
   profile%hysteresis%max_delta = 0.0_wp
   profile%hysteresis%at_value = 0.0_wp
   if (.not. profile%bidirectional) return

   do k = 0, profile%n_forward - 2
      ifwd = k + 1
      irev = profile%n_forward + (profile%n_forward - 1 - k)
      if (irev < 1 .or. irev > profile%n_points) cycle
      ! A pair with a missing energy is not compared (NaN compares false
      ! against everything), so it must not be *counted* either - otherwise an
      ! all-NaN profile reports `n_compared` points with a maximum difference
      ! of exactly 0.0, which reads like a perfectly reproducible scan.
      if (.not. profile%energy_valid(ifwd)) cycle
      if (.not. profile%energy_valid(irev)) cycle
      profile%hysteresis%n_compared = profile%hysteresis%n_compared + 1
      if (abs(profile%energy(ifwd) - profile%energy(irev)) > &
            & profile%hysteresis%max_delta) then
         profile%hysteresis%max_delta = &
            & abs(profile%energy(ifwd) - profile%energy(irev))
         profile%hysteresis%at_value = profile%values(ifwd)
      end if
   end do

end subroutine hysteresis_report

!> Largest energy jump between consecutive points of a profile.
subroutine jump_report(values, energy, npoints, hysteresis)
   !> Coordinate values of the points
   real(wp), intent(in) :: values(:)
   !> Energies of the points (hartree)
   real(wp), intent(in) :: energy(:)
   !> Number of points
   integer, intent(in) :: npoints
   !> Hysteresis report whose jump fields are filled
   type(landscape_scan_hysteresis), intent(inout) :: hysteresis

   real(wp) :: delta
   integer :: i

   hysteresis%max_jump = 0.0_wp
   hysteresis%jump_from = 0.0_wp
   hysteresis%jump_to = 0.0_wp
   do i = 2, npoints
      delta = abs(energy(i) - energy(i - 1))
      ! A failed point carries the NaN sentinel; comparing it would raise the
      ! invalid exception, so such a pair contributes no jump at all.
      if (.not. ieee_is_finite(delta)) cycle
      if (delta > hysteresis%max_jump) then
         hysteresis%max_jump = delta
         hysteresis%jump_from = values(i - 1)
         hysteresis%jump_to = values(i)
      end if
   end do

end subroutine jump_report

!> Print the scan profile and the hysteresis report.
!>
!> `coord_type` selects the user facing unit of the first column: angstrom for
!> a distance, degree for an angle or a dihedral.
subroutine landscape_scan_profile_write(profile, unit, coord_type)
   !> Scan profile
   type(landscape_scan_profile), intent(in) :: profile
   !> Output unit
   integer, intent(in) :: unit
   !> Coordinate type of the scanned coordinate
   integer, intent(in) :: coord_type

   real(wp) :: scale
   character(len=3) :: cunit
   integer :: i

   if (coord_type == p_landscape_scan_coord_distance) then
      scale = autoaa
      cunit = "Ang"
   else
      scale = 180.0_wp/pi
      cunit = "deg"
   end if

   ! The reached value and the drift are printed next to the requested value:
   ! the constraint is a finite spring, so a relaxed point sits at `q_target -
   ! F_phys/(2 fc)` rather than at `q_target`, and the size of that offset is
   ! the only per-point evidence of the systematic error of a scan
   ! (`docs/_landscape-recon/14-m2-validation-experiments.md` section 4.2).
   ! The drift column is the **wrapped** difference
   ! (`landscape_scan_coord_wrap_delta`), so a scan that walks through
   ! `+-180` degrees reports the real drift instead of a full turn
   ! (`docs/_landscape-recon/17-m2-independent-review.md` section 3.3).
   write(unit, '(1x,a,a,2(1x,a),1x,a)') "   point      value/", cunit, &
      & "reached/"//cunit, "drift/"//cunit, &
      & "        energy/Eh      residual      steps  conv"
   do i = 1, profile%n_points
      write(unit, '(1x,i8,1x,f16.8,1x,f16.8,1x,f16.8,1x,f20.12,1x,es12.4,1x,i6,1x,l1)') &
         & i, profile%values(i)*scale, profile%reached(i)*scale, &
         & landscape_scan_coord_wrap_delta(coord_type, &
            & profile%reached(i) - profile%values(i))*scale, &
         & profile%energy(i), profile%residual(i), profile%n_steps(i), &
         & profile%converged(i)
   end do
   if (profile%bidirectional) then
      write(unit, '(1x,a,i0)') "hysteresis: points compared ", &
         & profile%hysteresis%n_compared
      write(unit, '(1x,a,es15.6,a)') "hysteresis: max |E_fwd - E_rev| ", &
         & profile%hysteresis%max_delta, " Eh"
      write(unit, '(1x,a,f16.8,1x,a)') "hysteresis: at coordinate value ", &
         & profile%hysteresis%at_value*scale, trim(cunit)
   end if
   write(unit, '(1x,a,es15.6,a)') "largest jump between points: ", &
      & profile%hysteresis%max_jump, " Eh"

end subroutine landscape_scan_profile_write

end module xtb_landscape_scan1d
