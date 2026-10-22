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

!> Two-dimensional scan driver of the landscape suite.
!>
!> The grid is the Cartesian product of the two coordinate grids,
!> `values1(i) x values2(j)` for `i = 1..nstep1+1` and `j = 1..nstep2+1`,
!> anchored at the coordinate values of the input geometry (PLAN section 6.2).
!> It contains exactly `n1*n2` distinct points, each visited once.
!>
!> The visit order is *serpentine* (boustrophedon): rows are traversed in
!> alternating directions, so that consecutive points are always nearest
!> neighbours on the grid.  Together with `warm_start` this keeps the starting
!> geometry of every point as close as possible to the previous one, which is
!> the mitigation PLAN section 6.2 prescribes against hopping between basins
!> across the grid.
!>
!> The single scan point itself comes from
!> `xtb_landscape_scan1d::landscape_scan_solve_point`, so a 2D scan supports
!> the same `relaxed` and `rigid` modes as the 1D scan.
!>
!> Units: coordinates in bohr, energies in hartree, distances in bohr, angles
!> in radian.
module xtb_landscape_scan2d
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : pi
   use xtb_mctc_convert, only : autoaa
   use xtb_type_environment, only : TEnvironment
   use xtb_landscape_settings, only : landscape_scan_settings, &
      & p_landscape_scan_coord_distance, landscape_missing_value, &
      & landscape_value_present, landscape_real_to_string
   use xtb_landscape_engine, only : landscape_energy_model
   use xtb_landscape_coords, only : geometry_unpack
   use xtb_landscape_scan_relax, only : landscape_scan_point, &
      & landscape_scan_coord_value
   use xtb_landscape_scan1d, only : landscape_scan_solve_point, &
      & scan_jump_warning
   implicit none(type, external)
   private

   public :: landscape_scan_grid, landscape_scan2d, landscape_scan_grid_write

   !> Result of a two-dimensional scan.
   !>
   !> The energy, residual and convergence arrays are indexed by the grid
   !> coordinates `(i, j)`; `xyz`, `visit_i` and `visit_j` record the visit
   !> order of the serpentine traversal.
   type :: landscape_scan_grid
      !> Number of grid points along coordinate 1
      integer :: n1 = 0
      !> Number of grid points along coordinate 2
      integer :: n2 = 0
      !> Number of visited points, `n1*n2`
      integer :: n_points = 0
      !> Grid values of coordinate 1 (bohr or radian)
      real(wp), allocatable :: values1(:)
      !> Grid values of coordinate 2 (bohr or radian)
      real(wp), allocatable :: values2(:)
      !> De-biased energy (hartree), `(n1, n2)`, NaN if the point has none
      real(wp), allocatable :: energy(:, :)
      !> True if the grid point carries a usable energy; a frame of
      !> `xtbscan.log` is written for these points only, `(n1, n2)`
      logical, allocatable :: energy_valid(:, :)
      !> Biased energy (hartree), `(n1, n2)`
      real(wp), allocatable :: total(:, :)
      !> Constraint residual, `(n1, n2)`
      real(wp), allocatable :: residual(:, :)
      !> Inner convergence flag, `(n1, n2)`
      logical, allocatable :: converged(:, :)
      !> Geometry of every visited point (bohr), `(3N, n_points)`
      real(wp), allocatable :: xyz(:, :)
      !> Grid index of coordinate 1 in visit order
      integer, allocatable :: visit_i(:)
      !> Grid index of coordinate 2 in visit order
      integer, allocatable :: visit_j(:)
      !> Number of points whose inner minimisation failed
      integer :: n_fail = 0
      !> Largest energy jump between consecutive visited points (hartree)
      real(wp) :: max_jump = 0.0_wp
      !> Coordinate 1 value before the largest jump
      real(wp) :: jump_from1 = 0.0_wp
      !> Coordinate 2 value before the largest jump
      real(wp) :: jump_from2 = 0.0_wp
      !> Coordinate 1 value after the largest jump
      real(wp) :: jump_to1 = 0.0_wp
      !> Coordinate 2 value after the largest jump
      real(wp) :: jump_to2 = 0.0_wp
      !> Lowest de-biased energy of the grid (hartree)
      real(wp) :: min_energy = 0.0_wp
      !> Grid index of coordinate 1 at the minimum
      integer :: min_i = 0
      !> Grid index of coordinate 2 at the minimum
      integer :: min_j = 0
      !> Index in visit order at the minimum
      integer :: min_visit = 0
   end type landscape_scan_grid

contains

!> Run a two-dimensional Cartesian product scan.
!>
!> `bidirectional` governs the 1D scan only; the 2D grid is defined by its
!> two value lists and therefore covers the whole product by construction.
subroutine landscape_scan2d(env, base, at, settings, x0, grid, fail)
   character(len=*), parameter :: source = 'landscape_scan2d'
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
   type(landscape_scan_grid), intent(out) :: grid
   !> Failure flag
   logical, intent(out) :: fail

   type(landscape_scan_point) :: point
   real(wp), allocatable :: xstart(:), xout(:), xyz(:, :)
   real(wp) :: q1, q2, delta
   integer :: n1, n2, i, j, ipoint, istep, jstep
   logical :: point_fail

   fail = .false.
   n1 = settings%coord(1)%nstep + 1
   n2 = settings%coord(2)%nstep + 1

   allocate(xyz(3, size(at)))
   call geometry_unpack(x0, xyz)
   q1 = landscape_scan_coord_value(settings%coord(1), xyz)
   q2 = landscape_scan_coord_value(settings%coord(2), xyz)

   grid%n1 = n1
   grid%n2 = n2
   grid%n_points = n1*n2
   allocate(grid%values1(n1), grid%values2(n2))
   do i = 1, n1
      grid%values1(i) = q1 + real(i - 1, wp)*settings%coord(1)%step
   end do
   do j = 1, n2
      grid%values2(j) = q2 + real(j - 1, wp)*settings%coord(2)%step
   end do

   allocate(grid%energy(n1, n2), grid%total(n1, n2), grid%residual(n1, n2), &
      & grid%converged(n1, n2), grid%energy_valid(n1, n2))
   allocate(grid%xyz(3*size(at), n1*n2), grid%visit_i(n1*n2), &
      & grid%visit_j(n1*n2))
   grid%energy = landscape_missing_value()
   grid%total = 0.0_wp
   grid%residual = 0.0_wp
   grid%converged = .false.
   grid%energy_valid = .false.
   grid%xyz = 0.0_wp
   grid%visit_i = 0
   grid%visit_j = 0

   allocate(xstart(size(x0)), xout(size(x0)))
   xstart = x0

   ipoint = 0
   do jstep = 1, n2
      if (mod(jstep, 2) == 1) then
         do istep = 1, n1
            call solve_grid_point(istep, jstep)
         end do
      else
         do istep = n1, 1, -1
            call solve_grid_point(istep, jstep)
         end do
      end if
   end do

   ! The lowest energy point is searched over the *converged* points that
   ! carry an energy only: a failed point has no meaningful energy and must
   ! not be reported as the minimum of the grid (`min_energy` stays NaN when
   ! no point qualifies).
   grid%min_energy = landscape_missing_value()
   grid%min_i = 1
   grid%min_j = 1
   grid%min_visit = 1
   do ipoint = 1, grid%n_points
      i = grid%visit_i(ipoint)
      j = grid%visit_j(ipoint)
      if (.not. grid%converged(i, j)) cycle
      if (.not. grid%energy_valid(i, j)) cycle
      ! Block if, not one `.or.` condition: the first valid point has to set
      ! the baseline, and `min_energy` starts from the NaN sentinel, which
      ! Fortran does not require the right operand to short circuit around.
      if (.not. landscape_value_present(grid%min_energy)) then
         grid%min_energy = grid%energy(i, j)
         grid%min_i = i
         grid%min_j = j
         grid%min_visit = ipoint
      else if (grid%energy(i, j) < grid%min_energy) then
         grid%min_energy = grid%energy(i, j)
         grid%min_i = i
         grid%min_j = j
         grid%min_visit = ipoint
      end if
   end do

   do ipoint = 2, grid%n_points
      delta = abs(grid%energy(grid%visit_i(ipoint), grid%visit_j(ipoint)) - &
         & grid%energy(grid%visit_i(ipoint - 1), grid%visit_j(ipoint - 1)))
      ! A failed point carries the NaN sentinel; comparing it would raise the
      ! invalid exception, so such a pair contributes no jump at all.
      if (.not. ieee_is_finite(delta)) cycle
      if (delta > grid%max_jump) then
         grid%max_jump = delta
         grid%jump_from1 = grid%values1(grid%visit_i(ipoint - 1))
         grid%jump_from2 = grid%values2(grid%visit_j(ipoint - 1))
         grid%jump_to1 = grid%values1(grid%visit_i(ipoint))
         grid%jump_to2 = grid%values2(grid%visit_j(ipoint))
      end if
   end do

   if (grid%max_jump > scan_jump_warning) then
      call env%warning("the 2D scan grid jumps by "//&
         & trim(landscape_real_to_string(grid%max_jump, '(es15.6)'))//" hartree between two "//&
         & "consecutive points; this indicates a basin switch or a SCF "//&
         & "discontinuity", source)
   end if

   if (grid%n_fail > 0) then
      fail = .true.
      if (settings%printlevel > 0) then
         write(env%unit, '(1x,a,i0,a)') "landscape scan: ", grid%n_fail, &
            & " grid point(s) carry no energy (their single point failed)"
      end if
   end if

contains

   !> Solve one grid point, store it and record the visit order.
   subroutine solve_grid_point(istep, jstep)
      !> Grid index of coordinate 1
      integer, intent(in) :: istep
      !> Grid index of coordinate 2
      integer, intent(in) :: jstep

      real(wp) :: targets(2)

      targets(1) = grid%values1(istep)
      targets(2) = grid%values2(jstep)
      call landscape_scan_solve_point(env, base, at, settings, &
         & settings%coord(1:2), targets, 2, xstart, xout, point, point_fail)
      if (point_fail) then
         grid%n_fail = grid%n_fail + 1
         fail = .true.
      end if
      if (settings%warm_start .and. allocated(xout)) xstart = xout

      ipoint = ipoint + 1
      grid%visit_i(ipoint) = istep
      grid%visit_j(ipoint) = jstep
      grid%total(istep, jstep) = point%total
      grid%residual(istep, jstep) = point%residuals(1)
      grid%converged(istep, jstep) = point%converged
      grid%energy_valid(istep, jstep) = point%energy_valid
      if (point%energy_valid) then
         grid%energy(istep, jstep) = point%energy
      else
         ! A failed grid point has no energy; never store the `0.0` of an
         ! initialised but never evaluated point.
         grid%energy(istep, jstep) = landscape_missing_value()
      end if
      if (allocated(xout)) grid%xyz(:, ipoint) = xout

   end subroutine solve_grid_point

end subroutine landscape_scan2d

!> Print the grid in a compact table (one row per grid point of coordinate 2).
subroutine landscape_scan_grid_write(grid, unit, coord_type1, coord_type2)
   !> Scan grid
   type(landscape_scan_grid), intent(in) :: grid
   !> Output unit
   integer, intent(in) :: unit
   !> Coordinate type of the first coordinate
   integer, intent(in) :: coord_type1
   !> Coordinate type of the second coordinate
   integer, intent(in) :: coord_type2

   real(wp) :: scale1, scale2
   character(len=3) :: unit1, unit2
   integer :: i, j

   scale1 = user_scale(coord_type1)
   scale2 = user_scale(coord_type2)
   unit1 = user_unit(coord_type1)
   unit2 = user_unit(coord_type2)

   write(unit, '(1x,a,i0,1x,a,i0,a)') "2D scan grid: ", grid%n1, "x", &
      & grid%n2, " points"
   write(unit, '(1x,a,2(1x,a,a),1x,a)') "   coord2/", unit2, "   coord1/", &
      & unit1, "        energy/Eh      residual"
   do j = 1, grid%n2
      do i = 1, grid%n1
         write(unit, '(2(1x,f14.8),1x,f20.12,1x,es12.4)') &
            & grid%values2(j)*scale2, grid%values1(i)*scale1, &
            & grid%energy(i, j), grid%residual(i, j)
      end do
   end do
   write(unit, '(1x,a,es15.6,a)') "lowest grid energy: ", grid%min_energy, " Eh"
   write(unit, '(1x,a,2(1x,f14.8))') "at grid point:      ", &
      & grid%values1(grid%min_i)*scale1, grid%values2(grid%min_j)*scale2
   write(unit, '(1x,a,es15.6,a)') "largest jump between consecutive points: ", &
      & grid%max_jump, " Eh"

end subroutine landscape_scan_grid_write

!> User facing conversion factor of a coordinate type.
pure function user_scale(coord_type) result(scale)
   !> Coordinate type
   integer, intent(in) :: coord_type
   !> Conversion factor (bohr -> angstrom, radian -> degree)
   real(wp) :: scale

   if (coord_type == p_landscape_scan_coord_distance) then
      scale = autoaa
   else
      scale = 180.0_wp/pi
   end if

end function user_scale

!> User facing unit string of a coordinate type.
pure function user_unit(coord_type) result(unit)
   !> Coordinate type
   integer, intent(in) :: coord_type
   !> Unit string
   character(len=3) :: unit

   if (coord_type == p_landscape_scan_coord_distance) then
      unit = "Ang"
   else
      unit = "deg"
   end if

end function user_unit

end module xtb_landscape_scan2d
