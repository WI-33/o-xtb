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

!> Climbing image NEB (CI-NEB) force of the chain-of-states searches.
!>
!> Henkelman, Uberuaga & Jonsson, J. Chem. Phys. **113**, 9901 (2000).  One
!> image of the band - the one with the highest energy - is turned into a
!> "climbing image": it feels **no spring force** and the component of the
!> physical force along the band tangent is **inverted**, so it moves uphill
!> along the path and downhill perpendicular to it, towards the saddle point
!> on the minimum energy path.  With `F = -grad E` and the normalized band
!> tangent `tau`.
!>
!>     `F_ci = F - 2 (F . tau) tau`                                  (Eq. 5)
!>
!> The other images keep the plain NEB force (`F_perp + F_spring`) of
!> `xtb_landscape_path_neb`.  The two-stage protocol (the band is propagated
!> as a plain NEB until its largest perpendicular force falls below
!> `climb_switch_rms`, and only then is the climbing image switched on) is
!> kept as a safety margin, **not** as a measured accuracy requirement:
!> on the analytic LEPS surface and on the GFN2 ammonia inversion, igniting
!> at iteration 0 (the threshold raised until the initial band already
!> passes it) gives a bitwise identical barrier, saddle geometry and
!> iteration count, so the expected failure mode "the early climbing image
!> distorts the energy profile" was **not reproduced** (report 37 section
!> 3.1(e), unit test `climb-switch-calibration` and CLI tests
!> `cli-cineb-t1`/`cli-cineb-t10` in `m5-path-data/`).  The mechanism of Eq.
!> 5 is what makes the ignition moment unimportant here: the force always
!> pulls the image towards the saddle point along the band, whatever
!> geometry it starts from.
!>
!> The guard is retained because the climbing force treats the **band
!> tangent** as the reaction coordinate: while the band is still a straight
!> line or an IDPP path between distant endpoints (HCN -> HNC), that tangent
!> is the interpolation direction and not the direction of the path at the
!> highest image, so the reversed parallel component can push the image away
!> from the saddle instead of towards it.  Igniting late also avoids
!> freezing a climbing image index (it is frozen after ignition) chosen from
!> an image distribution that has not settled yet.  The cost of the guard is
!> the plain NEB iterations the band needs anyway.  The module holds no state: it
!> provides the force of Eq. 5, the deterministic choice of the image to
!> ignite and the ignition predicate, and the driver in
!> `xtb_landscape_path_neb` sequences them.
!>
!> Units: coordinates in bohr, energies in hartree, forces and gradients in
!> hartree/bohr; the tangent is a unit vector in the flattened Cartesian
!> space of one image.
module xtb_landscape_path_cineb
   use xtb_mctc_accuracy, only : wp
   implicit none(type, external)
   private

   public :: landscape_cineb_force
   public :: landscape_cineb_ignition_ready
   public :: landscape_cineb_highest_image

contains

!> Climbing image force of Eq. 5: `F - 2 (F . tau) tau`.
!>
!> `F` is the physical force (`-grad E`) on the climbing image and `tau` the
!> normalized band tangent at that image.  The projection removes the
!> component that would move the image off the path and reverses the sign of
!> the component along the path, turning the saddle into an attractor.
pure subroutine landscape_cineb_force(force, tau, f_ci)
   !> Physical force on the climbing image (hartree/bohr)
   real(wp), intent(in) :: force(:)
   !> Normalized band tangent at the climbing image (dimensionless)
   real(wp), intent(in) :: tau(:)
   !> Climbing image force (hartree/bohr)
   real(wp), intent(out) :: f_ci(:)

   real(wp) :: parallel

   parallel = dot_product(force, tau)
   f_ci = force - 2.0_wp*parallel*tau

end subroutine landscape_cineb_force

!> True if the climbing image may be switched on.
!>
!> The band is first a plain NEB; ignition waits until the largest
!> per-image perpendicular force RMS is **below** `climb_switch_rms`
!> (strictly below, so the strict inequality of the guard is the documented
!> behaviour).  `landscape_path_settings%climb_switch_rms` documents why the
!> measured gradient noise (`<= 1.2e-15` hartree/bohr random, `3.9e-7` for a
!> warm restart) can never trigger this threshold spuriously: the threshold
!> is 13 respectively 4 orders of magnitude larger, so it is the state of
!> the band that decides, not the numerical noise.
pure function landscape_cineb_ignition_ready(f_rms_max, climb_switch_rms) result(ready)
   !> Largest per-image perpendicular force RMS of the band (hartree/bohr)
   real(wp), intent(in) :: f_rms_max
   !> Ignition threshold (hartree/bohr)
   real(wp), intent(in) :: climb_switch_rms
   !> True if the climbing image may be switched on
   logical :: ready

   ready = f_rms_max < climb_switch_rms

end function landscape_cineb_ignition_ready

!> Index of the highest energy image among the internal images of a band.
!>
!> The endpoints are never candidates (they are fixed).  Ties are broken by
!> the lowest index, so the choice is deterministic.  Returns 0 if the band
!> has no internal image.
pure function landscape_cineb_highest_image(energies) result(image)
   !> Image energies (hartree), the endpoints included
   real(wp), intent(in) :: energies(:)
   !> 1-based index of the highest internal image, 0 if there is none
   integer :: image

   integer :: i

   image = 0
   do i = 2, size(energies) - 1
      ! A block if, not `image == 0 .or. ...`: Fortran does not require the
      ! right operand to be short circuited, so the first iteration would read
      ! `energies(0)`, which a bounds checked build rejects.
      if (image == 0) then
         image = i
      else if (energies(i) > energies(image)) then
         image = i
      end if
   end do

end function landscape_cineb_highest_image

end module xtb_landscape_path_cineb
