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

!> Constraint relaxation kernel of the landscape scan task.
!>
!> A relaxed scan point is a constrained minimisation of the electronic energy
!> `E_phys(x)` under the constraint `q(x) = q_target`.  Following the xtb
!> convention (`src/scan_driver.f90:69`, `src/constrain_pot.f90:313`) the
!> constraint is realised as an additive bias potential
!>
!>     E_bias(x) = fc * (q(x) - q_target)^expo        (expo = `spring_expo`)
!>
!> so that the inner minimiser works on `E_total = E_phys + E_bias`.  The
!> reported scan energy is the de-biased value `E_phys = E_total - E_bias` and
!> the residual
!>
!>     r = |dE_total/dq| = |grad_x E_total . grad_x q| / |grad_x q|^2
!>
!> measures how well the point satisfies `dE_phys/dq = -dE_bias/dq = lambda`,
!> the Lagrange multiplier of the holonomic constraint.  The residual is the
!> mitigation of risk R-04 (torsional scans of the reference implementation
!> show hysteresis and catastrophe points, see
!> `docs/_landscape-recon/04-web-research-ts-irc-state-of-art.md` section 4.6):
!> it is the per-point evidence that a silently non-converged inner relaxation
!> is visible in the scan output.
!>
!> All constraint machinery comes from `xtb_constrainpot::constrain_pot` and
!> `xtb_type_setvar::constr_setvar`, but the parameter object is built here
!> with explicit `fc`/`expo` arguments.  The module level globals of
!> `xtb_scanparam` (`fcconstr`, `springexpo`, `potset`, ...) are deliberately
!> never touched: they are hidden state, and PLAN N14 forbids relying on them.
!> `constraint_fc` is in hartree/bohr^2 for the documented default
!> `spring_expo = 2`.
!>
!> Units: coordinates in bohr, energies in hartree, gradients in hartree/bohr.
!> Distances are in bohr, angles and dihedrals in radian, all atom indices are
!> 1-based.  Flattened vectors follow the `x1, y1, z1, x2, ...` convention of
!> `xtb_landscape_coords`.
module xtb_landscape_scan_relax
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : pi
   use xtb_mctc_param, only : covalent_radius_2009
   use xtb_type_environment, only : TEnvironment
   use xtb_type_setvar, only : constr_setvar
   use xtb_constrainpot, only : constrain_pot
   use xtb_basic_geo, only : bangl, valijkl
   use xtb_landscape_settings, only : landscape_scan_settings, &
      & landscape_scan_coord_settings, landscape_opt_settings, &
      & p_landscape_scan_coord_distance, p_landscape_scan_coord_angle, &
      & p_landscape_scan_coord_dihedral, landscape_scan_coord_known, &
      & landscape_scan_coord_natoms, landscape_scan_coord_name, &
      & landscape_scan_mode_known, landscape_opt_algorithm_from_name, &
      & landscape_lowercase, landscape_missing_value, landscape_value_present, &
      & landscape_real_to_string
   use xtb_landscape_engine, only : landscape_energy_model
   use xtb_landscape_coords, only : geometry_pack, geometry_unpack
   use xtb_landscape_optimize, only : landscape_optimize_model
   use xtb_landscape_constraints, only : landscape_frozen_mask
   use xtb_landscape_results, only : landscape_result
   implicit none(type, external)
   private

   public :: landscape_scan_point
   public :: landscape_scan_relax_point, landscape_scan_rigid_point
   public :: landscape_scan_rigid_geometry
   public :: landscape_scan_coord_value, landscape_scan_coord_gradient
   public :: landscape_scan_coord_wrap_delta
   public :: landscape_scan_check, landscape_scan_inner_settings
   public :: landscape_scan_fragments
   public :: landscape_scan_bias_model, new_landscape_scan_bias_model

   !> Relative bond tolerance of the connectivity used by the rigid driver.
   !> A pair is bonded if its distance is below `bond_factor` times the sum of
   !> the covalent radii; 1.2 is the factor of
   !> `xtb_scanparam::constrain_all_bonds` (`src/scanparam.f90:322`).
   real(wp), parameter :: bond_factor = 1.2_wp
   !> Step of the central differences used for the Hessian of the constraint
   !> bias (bohr).  The bias is a smooth analytic function of the geometry, so
   !> this step keeps the truncation error far below the accuracy of the
   !> numerical Hessian of the electronic energy itself.
   real(wp), parameter :: bias_hess_step = 1.0e-4_wp
   !> Guard for divisions by distances and coordinate gradients
   real(wp), parameter :: geom_guard = 1.0e-10_wp
   !> Tolerated drift of a rigid scan point (bohr or radian).
   !>
   !> A rigid point is a pure geometry construction, so the reached coordinate
   !> has to reproduce its target up to round off; anything above a few times
   !> the double precision resolution of a Cartesian rotation is a
   !> meaningless construction (the moving fragment did not separate the two
   !> sides), not a numerical error.  `1e-8` is about `1e6` times the drift
   !> observed on 1,2-dichloroethane (`< 5e-9` degrees,
   !> `docs/_landscape-recon/16-m2-fix-report.md`) and still far below the
   !> `1e-3` radian resolution of the finest scan step the CLI offers.
   real(wp), parameter :: rigid_drift_tol = 1.0e-8_wp

   !> Result of one scan point.
   !>
   !> `total` is the biased energy at the returned geometry and
   !> `energy = total - bias` the de-biased electronic energy written to the
   !> scan log (the `Etot - Efix` of `src/scan_driver.f90:98`).  A point of an
   !> `ndim = 2` scan carries the requested, reached and residual values of
   !> both coordinates.
   !>
   !> Invariants of a point whose single point was discarded (unconverged SCC,
   !> terminated calculator): `energy_valid` is false, `energy` and `residuals`
   !> are quiet NaNs and `converged` is false - reporting the `0.0` an
   !> initialised but never evaluated point carries, or a convergence flag
   !> derived from the discarded gradient, would fabricate data
   !> (`docs/_landscape-recon/14-m2-validation-experiments.md` sections 5.3,
   !> 5.4).  `values` keeps the reached coordinate whenever the geometry exists
   !> (that value is a property of the construction, not of the discarded
   !> evaluation) and is NaN only when the geometry itself could not be built.
   type :: landscape_scan_point
      !> Number of scanned coordinates of the point (1 or 2)
      integer :: ncoord = 1
      !> Requested values of the scanned coordinates (bohr or radian)
      real(wp) :: targets(2) = 0.0_wp
      !> Values actually reached (bohr or radian)
      real(wp) :: values(2) = 0.0_wp
      !> Constraint residuals `|dE_total/dq_i|` (hartree/bohr or hartree/rad)
      real(wp) :: residuals(2) = 0.0_wp
      !> De-biased electronic energy (hartree), NaN if there is none
      real(wp) :: energy = 0.0_wp
      !> True if `energy` is a usable value
      logical :: energy_valid = .false.
      !> Constraint bias energy (hartree)
      real(wp) :: bias = 0.0_wp
      !> Biased energy (hartree)
      real(wp) :: total = 0.0_wp
      !> True if the inner minimisation converged (`rigid` points are exact)
      logical :: converged = .false.
      !> Macro steps of the inner minimisation (0 for a rigid point)
      integer :: n_steps = 0
      !> True if the point is a rigid geometry construction
      logical :: rigid = .false.
   end type landscape_scan_point

   !> Potential `E_phys + E_bias` on which a relaxed scan point minimises.
   type, extends(landscape_energy_model) :: landscape_scan_bias_model
      !> Underlying de-biased potential
      class(landscape_energy_model), allocatable :: base
      !> Single-constraint bias potential, built with explicit `fc`/`expo`
      type(constr_setvar) :: fix
      !> Number of atoms
      integer :: nat = 0
      !> Atom numbers as expected by `constrain_pot` (1-based, `nat` entries)
      integer, allocatable :: at(:)
   contains
      !> Energy and gradient of `E_phys + E_bias`
      procedure :: evaluate => scan_bias_model_evaluate
      !> Hessian of `E_phys + E_bias`
      procedure :: hessian => scan_bias_model_hessian
   end type landscape_scan_bias_model

contains

!> Value of an internal coordinate.
!>
!> Distances are returned in bohr, angles and dihedrals in radian.  The angle
!> and dihedral conventions are those of `xtb_basic_geo` (`bangl`, `valijkl`),
!> i.e. exactly the conventions of the constraint code in `xtb_constrainpot`;
!> the dihedral is in `(-pi, pi]`.
function landscape_scan_coord_value(coord, xyz) result(value)
   !> Coordinate definition
   type(landscape_scan_coord_settings), intent(in) :: coord
   !> Cartesian geometry (bohr)
   real(wp), intent(in) :: xyz(:, :)
   !> Coordinate value (bohr or radian)
   real(wp) :: value

   integer :: i, j, k, l, nat

   nat = size(xyz, 2)
   select case (coord%ctype)
   case (p_landscape_scan_coord_distance)
      i = coord%atoms(1)
      j = coord%atoms(2)
      value = norm2(xyz(:, j) - xyz(:, i))
   case (p_landscape_scan_coord_angle)
      i = coord%atoms(1)
      j = coord%atoms(2)
      k = coord%atoms(3)
      call bangl(xyz, i, j, k, value)
   case (p_landscape_scan_coord_dihedral)
      i = coord%atoms(1)
      j = coord%atoms(2)
      k = coord%atoms(3)
      l = coord%atoms(4)
      value = valijkl(nat, xyz, i, j, k, l)
   case default
      value = 0.0_wp
   end select

end function landscape_scan_coord_value

!> Cartesian gradient `grad_x q` of an internal coordinate.
!>
!> The gradient is obtained from the constraint code itself: `constrain_pot`
!> is evaluated with a *linearised* single constraint for which
!> `E_bias(x) = q(x)` up to an irrelevant constant, so the returned gradient
!> is exactly `grad_x q`.  Reusing the constraint code (instead of a second,
!> independently written derivative) keeps the residual of a scan point and
!> the bias potential that generated it mathematically consistent.
subroutine landscape_scan_coord_gradient(coord, at, xyz, grad)
   !> Coordinate definition
   type(landscape_scan_coord_settings), intent(in) :: coord
   !> Atom numbers (1-based, `size(xyz, 2)` entries)
   integer, intent(in) :: at(:)
   !> Cartesian geometry (bohr)
   real(wp), intent(in) :: xyz(:, :)
   !> Cartesian gradient of the coordinate
   real(wp), intent(out) :: grad(:, :)

   type(constr_setvar) :: fix
   type(landscape_scan_coord_settings) :: single(1)
   real(wp) :: e, target(1)

   grad = 0.0_wp
   e = 0.0_wp
   single(1) = coord
   select case (coord%ctype)
   case (p_landscape_scan_coord_distance)
      ! E = 1 * d^1 = r - r0  ->  grad = (x_j - x_i)/r = grad r
      target(1) = 0.0_wp
      call landscape_scan_fix_setup(fix, size(xyz, 2), single, target, 1, &
         & 1.0_wp, 1.0_wp)
   case (p_landscape_scan_coord_angle)
      ! E = (theta - theta0)^2 with theta0 = theta - 1/2  ->  dE/dtheta = 1
      target(1) = landscape_scan_coord_value(coord, xyz) - 0.5_wp
      call landscape_scan_fix_setup(fix, size(xyz, 2), single, target, 1, &
         & 1.0_wp, 2.0_wp)
   case (p_landscape_scan_coord_dihedral)
      ! E = 1 + cos(phi0 - phi + pi) with phi0 = phi - pi/2  ->  dE/dphi = 1
      target(1) = landscape_scan_coord_value(coord, xyz) - 0.5_wp*pi
      call landscape_scan_fix_setup(fix, size(xyz, 2), single, target, 1, &
         & 1.0_wp, 2.0_wp)
   case default
      return
   end select
   call constrain_pot(fix, size(xyz, 2), at, xyz, grad, e)

end subroutine landscape_scan_coord_gradient

!> Reduce a coordinate difference to a single period.
!>
!> A distance is not periodic, so its difference is returned unchanged.
!> Angles and dihedrals live on a circle: the difference of two values of
!> `(-pi, pi]` spans `(-2 pi, 2 pi]`, and a target the scan reaches across
!> the branch cut makes the raw difference look like a full turn.  The
!> wrapped difference lies in `[-pi, pi)`, so a scan that walks through
!> `+-180` degrees reports the real drift instead of `360` degrees minus it
!> (`docs/_landscape-recon/17-m2-independent-review.md` section 3.3).
!>
!> The function is the shared definition of the drift of a scan point: the
!> rigid drift criterion and the `drift` column of the 1D scan table both
!> use it.
pure function landscape_scan_coord_wrap_delta(ctype, delta) result(wrapped)
   !> Coordinate type, one of the `p_landscape_scan_coord_*` constants
   integer, intent(in) :: ctype
   !> Raw difference `reached - target` (bohr or radian)
   real(wp), intent(in) :: delta
   !> Difference reduced to one period (bohr or radian)
   real(wp) :: wrapped

   if (ctype == p_landscape_scan_coord_distance) then
      wrapped = delta
   else
      wrapped = modulo(delta + pi, 2.0_wp*pi) - pi
   end if

end function landscape_scan_coord_wrap_delta

!> Build the constraint parameter object of `ncoord` scanned coordinates.
!>
!> The force constant and the exponent are explicit arguments; the resulting
!> bias is `E = fc * sum_i (q_i - target_i)^expo`.  A mixed set such as one
!> distance and one angle is supported because `constr_setvar` keeps a
!> separate list per coordinate type.  The dihedral branch of `constrain_pot`
!> uses the periodic form `E = fc * (1 + cos(phi0 - phi + pi))`, which agrees
!> with the harmonic well to leading order and is bounded by `2 fc`.
subroutine landscape_scan_fix_setup(fix, nat, coords, targets, ncoord, fc, expo)
   !> Constraint set, allocated on return
   type(constr_setvar), intent(out) :: fix
   !> Number of atoms
   integer, intent(in) :: nat
   !> Scanned coordinates
   type(landscape_scan_coord_settings), intent(in) :: coords(:)
   !> Target values of the coordinates (bohr or radian)
   real(wp), intent(in) :: targets(:)
   !> Number of coordinates to set up
   integer, intent(in) :: ncoord
   !> Force constant (hartree/bohr^2 for the default exponent two)
   real(wp), intent(in) :: fc
   !> Exponent of the distance bias
   real(wp), intent(in) :: expo

   integer :: i, n

   call fix%allocate(nat, ncoord, fc=fc, expo=expo)
   fix%dist%n = 0
   fix%angle%n = 0
   fix%dihedral%n = 0

   do i = 1, ncoord
      select case (coords(i)%ctype)
      case (p_landscape_scan_coord_distance)
         n = fix%dist%n + 1
         fix%dist%n = n
         fix%dist%atoms(2*n - 1:2*n) = coords(i)%atoms(1:2)
         fix%dist%val(n) = targets(i)
         fix%dist%expo(n) = expo
      case (p_landscape_scan_coord_angle)
         n = fix%angle%n + 1
         fix%angle%n = n
         fix%angle%atoms(3*n - 2:3*n) = coords(i)%atoms(1:3)
         fix%angle%val(n) = targets(i)
      case (p_landscape_scan_coord_dihedral)
         n = fix%dihedral%n + 1
         fix%dihedral%n = n
         fix%dihedral%atoms(4*n - 3:4*n) = coords(i)%atoms(1:4)
         fix%dihedral%val(n) = targets(i)
      end select
   end do

end subroutine landscape_scan_fix_setup

!> Bias energy and Cartesian gradient of a single-constraint set.
subroutine landscape_scan_bias_potential(fix, at, xyz, energy, gradient)
   !> Constraint set
   type(constr_setvar), intent(in) :: fix
   !> Atom numbers (1-based, `size(xyz, 2)` entries)
   integer, intent(in) :: at(:)
   !> Cartesian geometry (bohr)
   real(wp), intent(in) :: xyz(:, :)
   !> Bias energy (hartree)
   real(wp), intent(out) :: energy
   !> Cartesian bias gradient (hartree/bohr)
   real(wp), intent(out) :: gradient(:, :)

   energy = 0.0_wp
   gradient = 0.0_wp
   call constrain_pot(fix, size(xyz, 2), at, xyz, gradient, energy)

end subroutine landscape_scan_bias_potential

!> Construct the biased potential of one scan point.
subroutine new_landscape_scan_bias_model(self, base, coords, targets, ncoord, &
      & at, fc, expo)
   !> New model
   type(landscape_scan_bias_model), intent(out) :: self
   !> Underlying de-biased potential
   class(landscape_energy_model), intent(in) :: base
   !> Scanned coordinates
   type(landscape_scan_coord_settings), intent(in) :: coords(:)
   !> Target values of the coordinates (bohr or radian)
   real(wp), intent(in) :: targets(:)
   !> Number of coordinates to constrain
   integer, intent(in) :: ncoord
   !> Atom numbers (1-based)
   integer, intent(in) :: at(:)
   !> Force constant of the bias
   real(wp), intent(in) :: fc
   !> Exponent of the distance bias
   real(wp), intent(in) :: expo

   allocate(self%base, source=base)
   self%nat = size(at)
   self%at = at
   call landscape_scan_fix_setup(self%fix, self%nat, coords, targets, ncoord, &
      & fc, expo)

end subroutine new_landscape_scan_bias_model

!> Energy and gradient of the biased potential at the flattened coordinates.
subroutine scan_bias_model_evaluate(self, x, energy, gradient, fail)
   !> Biased potential
   class(landscape_scan_bias_model), intent(inout) :: self
   !> Flattened coordinates (bohr)
   real(wp), intent(in) :: x(:)
   !> Biased energy (hartree)
   real(wp), intent(out) :: energy
   !> Flattened biased gradient (hartree/bohr)
   real(wp), intent(out) :: gradient(:)
   !> Failure flag
   logical, intent(out) :: fail

   real(wp), allocatable :: xyz(:, :), gbias(:, :)
   real(wp) :: ebias

   call self%base%evaluate(x, energy, gradient, fail)
   if (fail) return

   allocate(xyz(3, self%nat), gbias(3, self%nat))
   call geometry_unpack(x, xyz)
   call landscape_scan_bias_potential(self%fix, self%at, xyz, ebias, gbias)
   energy = energy + ebias
   gradient = gradient + geometry_pack(gbias)

end subroutine scan_bias_model_evaluate

!> Hessian of the biased potential.
!>
!> The Hessian of the bias is obtained by central differences of its analytic
!> gradient (no electronic structure evaluation), the Hessian of the base
!> potential is asked from the base model itself.  The biased Hessian is only
!> needed when a relaxed scan uses the `rfo` inner minimiser.
subroutine scan_bias_model_hessian(self, x, hess, fail)
   !> Biased potential
   class(landscape_scan_bias_model), intent(inout) :: self
   !> Flattened coordinates (bohr)
   real(wp), intent(in) :: x(:)
   !> Dense Hessian (hartree/bohr^2)
   real(wp), intent(out) :: hess(:, :)
   !> Failure flag
   logical, intent(out) :: fail

   real(wp), allocatable :: xp(:), xyz(:, :), gplus(:, :), gminus(:, :)
   real(wp) :: e
   integer :: i, j, ic, n3

   call self%base%hessian(x, hess, fail)
   if (fail) return

   n3 = 3*self%nat
   allocate(xp(n3), xyz(3, self%nat), gplus(3, self%nat), gminus(3, self%nat))
   do j = 1, n3
      xp = x
      xp(j) = xp(j) + bias_hess_step
      call geometry_unpack(xp, xyz)
      call landscape_scan_bias_potential(self%fix, self%at, xyz, e, gplus)
      xp(j) = x(j) - bias_hess_step
      call geometry_unpack(xp, xyz)
      call landscape_scan_bias_potential(self%fix, self%at, xyz, e, gminus)
      do i = 1, self%nat
         ic = 3*(i - 1)
         hess(ic + 1:ic + 3, j) = hess(ic + 1:ic + 3, j) &
            & + (gplus(:, i) - gminus(:, i))/(2.0_wp*bias_hess_step)
      end do
   end do
   hess = 0.5_wp*(hess + transpose(hess))

end subroutine scan_bias_model_hessian

!> Validate the scan settings against a molecule with `nat` atoms.
!>
!> Every precondition the scan driver relies on is checked here, so the
!> failure paths are reachable without running a single SCF: an unknown
!> coordinate type, an atom index outside `1..nat`, a repeated atom inside one
!> coordinate, `nstep < 1`, a zero step size, an unknown scan mode or an
!> unknown inner algorithm.  All problems are reported through `env`.
function landscape_scan_check(env, settings, nat) result(ok)
   character(len=*), parameter :: source = 'landscape_scan_check'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Scan settings
   type(landscape_scan_settings), intent(in) :: settings
   !> Number of atoms of the molecule
   integer, intent(in) :: nat
   !> True if all preconditions hold
   logical :: ok

   integer :: i, j, natoms, iatom
   character(len=160) :: buffer
   character(len=8) :: step_unit
   character(len=:), allocatable :: cname, flag_suffix, atom_list

   ok = .true.

   if (settings%ndim < 1 .or. settings%ndim > 2) then
      write(buffer, '(i0)') settings%ndim
      call env%error("a landscape scan supports ndim = 1 or 2, not "//&
         & trim(buffer), source)
      ok = .false.
      return
   end if
   if (nat < 2) then
      call env%error("a landscape scan needs at least two atoms", source)
      ok = .false.
      return
   end if
   if (.not. landscape_scan_mode_known(settings%mode)) then
      call env%error("unknown scan mode '"//trim(settings%mode)//&
         & "', expected 'relaxed' or 'rigid'", source)
      ok = .false.
      return
   end if
   if (landscape_lowercase(settings%mode) == 'relaxed') then
      if (landscape_opt_algorithm_from_name(settings%method) == 0) then
         call env%error("unknown inner minimiser '"//trim(settings%method)//&
            & "' for a relaxed scan", source)
         ok = .false.
         return
      end if
      ! A non-finite value reaches this comparison unchecked (`nan <= 0` is
      ! false), and it then poisons every biased energy of the scan: the CLI
      ! parser accepts `nan`/`inf` for `--lscan-fc`, so the guard has to name
      ! the value instead of letting the run produce a NaN profile.
      if (.not. ieee_is_finite(settings%constraint_fc)) then
         write(buffer, '(a,es15.6,a)') "--lscan-fc ", settings%constraint_fc, &
            & " is not a finite force constant"
         call env%error(trim(buffer), source)
         ok = .false.
         return
      end if
      if (settings%constraint_fc <= 0.0_wp) then
         call env%error("the constraint force constant must be positive", source)
         ok = .false.
         return
      end if
      if (settings%spring_expo <= 0.0_wp) then
         call env%error("the constraint potential exponent must be positive", source)
         ok = .false.
         return
      end if
   end if

   do i = 1, settings%ndim
      if (.not. landscape_scan_coord_known(settings%coord(i)%ctype)) then
         write(buffer, '(i0)') settings%coord(i)%ctype
         call env%error("unknown internal coordinate type "//trim(buffer)//&
            & " for coordinate 1-based index "//one_digit(i), source)
         ok = .false.
         return
      end if
      cname = landscape_scan_coord_name(settings%coord(i)%ctype)
      natoms = landscape_scan_coord_natoms(settings%coord(i)%ctype)
      ! The messages point at the flag the user really has to pass (CLI
      ! redesign section 8.2).  The second coordinate carries the `2` suffixed
      ! flags (`--lscan-atoms2`, `--lscan-step2`, `--lscan-nstep2`), and the
      ! unit of the step is angstrom for a distance but degree for an angle or
      ! a dihedral (the command line converts both to atomic units).
      flag_suffix = ""
      if (i == 2) flag_suffix = "2"
      atom_list = "1"
      do j = 2, natoms
         atom_list = atom_list//","//nat_to_string(j)
      end do
      if (settings%coord(i)%ctype == p_landscape_scan_coord_distance) then
         step_unit = "angstrom"
      else
         step_unit = "degree"
      end if
      do j = 1, natoms
         iatom = settings%coord(i)%atoms(j)
         if (iatom < 1 .or. iatom > nat) then
            write(buffer, '(i0)') iatom
            call env%error("the "//cname//" coordinate needs --lscan-atoms"//&
               & flag_suffix//" "//atom_list//" (atom indices are 1-based, "//&
               & "1.."//nat_to_string(nat)//"); index "//trim(buffer)//&
               & " is invalid", source)
            ok = .false.
            return
         end if
      end do
      do j = 1, natoms
         if (count(settings%coord(i)%atoms(1:natoms) == &
               & settings%coord(i)%atoms(j)) > 1) then
            call env%error("the "//cname//" coordinate uses atom "//&
               & nat_to_string(settings%coord(i)%atoms(j))//&
               & " more than once", source)
            ok = .false.
            return
         end if
      end do
      if (settings%coord(i)%nstep < 1) then
         call env%error("coordinate "//one_digit(i)//" has nstep = "//&
            & nat_to_string(settings%coord(i)%nstep)//&
            & "; set it with --lscan-nstep"//flag_suffix//" <n>", source)
         ok = .false.
         return
      end if
      if (settings%coord(i)%step == 0.0_wp) then
         call env%error("coordinate "//one_digit(i)//" has a zero step "//&
            & "size; set it with --lscan-step"//flag_suffix//" <value> ("//&
            & trim(step_unit)//")", source)
         ok = .false.
         return
      end if
   end do

end function landscape_scan_check

!> Map the scan settings onto the inner minimisation settings.
!>
!> The inner minimiser shares the convergence thresholds of the scan and never
!> writes files (the scan driver owns all output).  Its algorithm is the one
!> named by `landscape_scan_settings%method`.
function landscape_scan_inner_settings(settings) result(opt)
   !> Scan settings
   type(landscape_scan_settings), intent(in) :: settings
   !> Settings of the inner minimisation
   type(landscape_opt_settings) :: opt

   opt = landscape_opt_settings()
   opt%algorithm = landscape_opt_algorithm_from_name(settings%method)
   opt%printlevel = 0
   opt%write_files = .false.
   opt%max_steps = settings%max_steps
   opt%e_tol = settings%e_tol
   opt%g_max = settings%g_max
   opt%g_rms = settings%g_rms
   opt%d_max = settings%d_max
   opt%d_rms = settings%d_rms
   opt%engine = settings%engine

end function landscape_scan_inner_settings

!> Relax one scan point under the constraints `q_i = target_i`.
!>
!> The point is the minimisation of `E_phys + E_bias` started from `x0`; the
!> returned geometry `xout` and the de-biased energy `point%energy` are the
!> scan output.  `result` is the inner minimisation result and `fail` reports
!> an inner failure such as an unconverged SCC or a missing geometry.
subroutine landscape_scan_relax_point(env, base, settings, coords, targets, &
      & ncoord, at, x0, xout, point, result, fail)
   character(len=*), parameter :: source = 'landscape_scan_relax_point'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> De-biased potential
   class(landscape_energy_model), intent(inout) :: base
   !> Scan settings
   type(landscape_scan_settings), intent(in) :: settings
   !> Scanned coordinates
   type(landscape_scan_coord_settings), intent(in) :: coords(:)
   !> Target values of the coordinates (bohr or radian)
   real(wp), intent(in) :: targets(:)
   !> Number of constraints
   integer, intent(in) :: ncoord
   !> Atom numbers (1-based)
   integer, intent(in) :: at(:)
   !> Starting geometry (bohr)
   real(wp), intent(in) :: x0(:)
   !> Relaxed geometry (bohr)
   real(wp), allocatable, intent(out) :: xout(:)
   !> Scan point record
   type(landscape_scan_point), intent(out) :: point
   !> Inner minimisation result
   type(landscape_result), intent(out) :: result
   !> Failure flag
   logical, intent(out) :: fail

   type(landscape_scan_bias_model) :: model
   type(landscape_opt_settings) :: opt
   real(wp), allocatable :: xyz(:, :), gbias(:, :), gcoord(:, :)
   real(wp), allocatable :: gphys(:), gflat(:)
   integer :: i

   fail = .false.
   point%ncoord = ncoord
   point%targets = 0.0_wp
   point%targets(1:ncoord) = targets(1:ncoord)
   point%rigid = .false.
   point%converged = .false.
   point%energy = landscape_missing_value()
   point%energy_valid = .false.
   point%values = landscape_missing_value()
   point%residuals = landscape_missing_value()
   ! A point without an energy has no biased energy either; the `0.0` of an
   ! initialised but never evaluated point is not a value.
   point%total = landscape_missing_value()

   call new_landscape_scan_bias_model(model, base, coords, targets, ncoord, &
      & at, settings%constraint_fc, settings%spring_expo)
   opt = landscape_scan_inner_settings(settings)

   call landscape_optimize_model(env, model, opt, x0, result, xout, &
      & frozen=landscape_frozen_mask(size(at)))
   if (.not. allocated(xout)) then
      call env%error("the inner minimisation returned no geometry", source)
      fail = .true.
      return
   end if

   ! De-bias the point: evaluate the electronic energy and gradient at the
   ! relaxed geometry (a cache hit inside the engine) and the bias at the same
   ! geometry, then form the residual on every scanned coordinate.
   allocate(xyz(3, size(at)), gbias(3, size(at)), gcoord(3, size(at)))
   allocate(gphys(size(xout)), gflat(size(xout)))
   call geometry_unpack(xout, xyz)
   ! The reached coordinates are a property of the geometry, so they are
   ! reported even when the energy of that geometry is discarded below.
   do i = 1, ncoord
      point%values(i) = landscape_scan_coord_value(coords(i), xyz)
   end do
   call base%evaluate(xout, point%energy, gphys, fail)
   if (fail) then
      ! The point has no usable energy: keep the NaN of the missing value and
      ! the `.false.` of the convergence flag instead of reporting the `0.0`
      ! left behind by the failed evaluation.
      point%energy = landscape_missing_value()
      point%energy_valid = .false.
      point%converged = .false.
      call env%error("the de-biased energy of a scan point failed", source)
      return
   end if
   call landscape_scan_bias_potential(model%fix, model%at, xyz, point%bias, gbias)

   gflat = gphys + geometry_pack(gbias)
   do i = 1, ncoord
      call landscape_scan_coord_gradient(coords(i), at, xyz, gcoord)
      point%residuals(i) = coordinate_residual(gflat, gcoord)
   end do
   point%total = point%energy + point%bias
   point%converged = result%converged
   point%n_steps = result%n_steps
   if (.not. landscape_value_present(point%energy)) then
      ! A non-finite de-biased energy is not a value either; it must reach the
      ! caller as "missing", never as a number.
      point%energy = landscape_missing_value()
      point%energy_valid = .false.
      point%converged = .false.
      fail = .true.
      call env%error("the de-biased energy of a scan point is not finite", &
         & source)
      return
   end if
   point%energy_valid = .true.

end subroutine landscape_scan_relax_point

!> Build one rigid scan point: move the coordinates to their targets and
!> evaluate the resulting geometry once.
!>
!> No relaxation is performed, so the point is a plain single point of the
!> constructed geometry.  Its residuals are the forces along the scanned
!> coordinates, `|dE_phys/dq_i|`, which are generally **not** zero - for a
!> rigid scan the residual is exactly the quantity that quantifies how far the
!> rigid geometry is from the relaxed one.
subroutine landscape_scan_rigid_point(env, base, coords, targets, ncoord, at, &
      & x0, xout, point, fail)
   character(len=*), parameter :: source = 'landscape_scan_rigid_point'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Potential
   class(landscape_energy_model), intent(inout) :: base
   !> Scanned coordinates
   type(landscape_scan_coord_settings), intent(in) :: coords(:)
   !> Target values of the coordinates (bohr or radian)
   real(wp), intent(in) :: targets(:)
   !> Number of constraints
   integer, intent(in) :: ncoord
   !> Atom numbers (1-based)
   integer, intent(in) :: at(:)
   !> Starting geometry (bohr)
   real(wp), intent(in) :: x0(:)
   !> Constructed geometry (bohr)
   real(wp), allocatable, intent(out) :: xout(:)
   !> Scan point record
   type(landscape_scan_point), intent(out) :: point
   !> Failure flag
   logical, intent(out) :: fail

   real(wp), allocatable :: xyz(:, :), gcoord(:, :), xwork(:, :)
   real(wp), allocatable :: gphys(:)
   real(wp) :: drift
   integer :: i

   fail = .false.
   point%ncoord = ncoord
   point%targets = 0.0_wp
   point%targets(1:ncoord) = targets(1:ncoord)
   point%rigid = .true.
   ! A rigid point is exact by construction, so it converges as soon as its
   ! single point succeeds; a discarded single point leaves it unconverged
   ! with no energy instead (see the invariant of `landscape_scan_point`).
   point%converged = .false.
   point%n_steps = 0
   point%bias = 0.0_wp
   point%energy = landscape_missing_value()
   point%energy_valid = .false.
   point%values = landscape_missing_value()
   point%residuals = landscape_missing_value()
   ! ... see `landscape_scan_relax_point`: no energy means no biased energy
   point%total = landscape_missing_value()

   ! Apply the rigid operations one after another on the same geometry.  The
   ! second operation can disturb the first coordinate weakly (the fragment
   ! rotated for the second coordinate may contain atoms of the first one);
   ! this is inherent to a rigid scan and is reported through the residuals.
   allocate(xwork(3, size(at)))
   call geometry_unpack(x0, xwork)
   do i = 1, ncoord
      call landscape_scan_rigid_geometry(env, coords(i), at, xwork, &
         & targets(i), fail)
      if (fail) return
   end do
   xout = geometry_pack(xwork)

   allocate(xyz(3, size(at)), gcoord(3, size(at)))
   allocate(gphys(size(xout)))
   call geometry_unpack(xout, xyz)
   ! the constructed coordinate values are reported independently of the
   ! single point that follows
   do i = 1, ncoord
      point%values(i) = landscape_scan_coord_value(coords(i), xyz)
   end do
   ! A rigid point is a pure geometry construction, so the reached value has
   ! to reproduce the target up to round off.  A larger (wrapped) difference
   ! means the moving fragment did not separate the two sides and the scanned
   ! coordinate did not move; without this check the point would be reported
   ! as a perfectly converged one
   ! (`docs/_landscape-recon/17-m2-independent-review.md` section 3.2).  For
   ! `ndim = 2` the second operation can disturb the first coordinate weakly
   ! (see above); that is inherent to a coupled rigid scan, so there the
   ! drift is only a warning.
   do i = 1, ncoord
      drift = landscape_scan_coord_wrap_delta(coords(i)%ctype, &
         & point%values(i) - targets(i))
      if (abs(drift) <= rigid_drift_tol) cycle
      if (ncoord > 1) then
         call env%warning("rigid scan coordinate "//one_digit(i)//&
            & " missed its target by "//trim(landscape_real_to_string(drift, '(es15.6)'))//&
            & "; the second rigid operation disturbed the first one", &
            & source)
         cycle
      end if
      call env%error("a rigid scan point missed its target by "//&
         & trim(landscape_real_to_string(drift, '(es15.6)'))//" (bohr or radian): the moving "//&
         & "fragment did not separate the two sides of the molecule", &
         & source)
      fail = .true.
      return
   end do
   call base%evaluate(xout, point%energy, gphys, fail)
   if (fail .or. .not. landscape_value_present(point%energy)) then
      point%energy = landscape_missing_value()
      point%energy_valid = .false.
      point%converged = .false.
      fail = .true.
      call env%error("the single point of a rigid scan point failed", source)
      return
   end if
   point%energy_valid = .true.
   point%converged = .true.

   point%total = point%energy
   do i = 1, ncoord
      call landscape_scan_coord_gradient(coords(i), at, xyz, gcoord)
      point%residuals(i) = coordinate_residual(gphys, gcoord)
   end do

end subroutine landscape_scan_rigid_point

!> Residual force along an internal coordinate, `|dE/dq|` (hartree/bohr or /rad).
!>
!> `|dE/dq| = |g . g_q| / |g_q|^2` with `g` the Cartesian gradient and
!> `g_q = grad_x q`.  A vanishing `g_q` (a near-linear angle, a collapsed
!> distance) is reported as a zero residual: the coordinate is not
!> differentiable there and no meaningful projection exists.
pure function coordinate_residual(g, gq) result(residual)
   !> Flat Cartesian gradient (hartree/bohr), length 3N
   real(wp), intent(in) :: g(:)
   !> Cartesian gradient of the coordinate (3, N)
   real(wp), intent(in) :: gq(:, :)
   !> `|dE/dq|`
   real(wp) :: residual

   real(wp), allocatable :: gflat(:)
   real(wp) :: gnorm2

   allocate(gflat(size(g)))
   gflat = geometry_pack(gq)
   gnorm2 = sum(gflat**2)
   if (gnorm2 > geom_guard) then
      residual = abs(sum(g*gflat))/gnorm2
   else
      residual = 0.0_wp
   end if

end function coordinate_residual

!> Rigidly move the internal coordinate to `target`.
!>
!> The geometry is modified in place by a rigid body operation of the
!> molecular fragments defined by the scanned coordinate (PLAN section 6.2):
!>
!>  * distance - the two fragments are translated along the bond axis by
!>    `(target - r)/2` each, so exactly the scanned distance changes while
!>    every intra-fragment distance is preserved;
!>  * angle - the fragment of the first atom is rotated about the axis through
!>    the central atom perpendicular to the plane, which changes exactly the
!>    scanned angle;
!>  * dihedral - the far side of the central bond `x_j - x_k` is rotated about
!>    it, i.e. the fragment that contains the third atom `k` and therefore the
!>    last one `l`; the whole side turns as a rigid body, which changes exactly
!>    the scanned dihedral (the `valijkl` convention increases the dihedral
!>    under this rotation).  The two central atoms lie on the rotation axis and
!>    do not move, so rotating that side and rotating only `l` change the
!>    scanned dihedral by the same amount; they are *not* interchangeable,
!>    however, because a rigid scan has to move every substituent of the far
!>    side (see the `rigid-dihedral-fragment` unit test).
!>
!> The operation on a two atom fragment is exact and hand checkable, which is
!> what the `rigid-scan-*` unit tests assert.  The cut has to separate the
!> molecule into two sides; a ring, a non bonded cut and a bridge that the
!> distance based connectivity swallows are refused by
!> `landscape_scan_check_cut` instead of silently moving both sides (which
!> would freeze the scanned coordinate, the failure mode of the wrong
!> dihedral seeding, `docs/_landscape-recon/17-m2-independent-review.md`
!> section 3.2).
subroutine landscape_scan_rigid_geometry(env, coord, at, xyz, target, fail)
   character(len=*), parameter :: source = 'landscape_scan_rigid_geometry'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Scanned coordinate
   type(landscape_scan_coord_settings), intent(in) :: coord
   !> Atom numbers (1-based)
   integer, intent(in) :: at(:)
   !> Cartesian geometry (bohr), moved in place
   real(wp), intent(inout) :: xyz(:, :)
   !> Target value of the coordinate (bohr or radian)
   real(wp), intent(in) :: target
   !> Failure flag
   logical, intent(out) :: fail

   real(wp) :: axis(3), delta, alpha
   logical, allocatable :: moving(:), moving2(:)
   integer :: nat, i, j, k, l

   fail = .false.
   nat = size(xyz, 2)

   select case (coord%ctype)
   case (p_landscape_scan_coord_distance)
      i = coord%atoms(1)
      j = coord%atoms(2)
      axis = xyz(:, j) - xyz(:, i)
      delta = norm2(axis)
      if (delta <= geom_guard) then
         call env%error("cannot scan a distance between coincident atoms", source)
         fail = .true.
         return
      end if
      axis = axis/delta
      alpha = target - delta
      call landscape_scan_fragments(at, xyz, i, j, moving)
      call landscape_scan_fragments(at, xyz, j, i, moving2)
      call landscape_scan_check_cut(env, moving, j, fail)
      if (fail) return
      call landscape_scan_check_cut(env, moving2, i, fail)
      if (fail) return
      do k = 1, nat
         if (moving(k)) xyz(:, k) = xyz(:, k) - 0.5_wp*alpha*axis
         if (moving2(k)) xyz(:, k) = xyz(:, k) + 0.5_wp*alpha*axis
      end do
   case (p_landscape_scan_coord_angle)
      i = coord%atoms(1)
      j = coord%atoms(2)
      k = coord%atoms(3)
      alpha = landscape_scan_coord_value(coord, xyz) - target
      axis = cross(xyz(:, i) - xyz(:, j), xyz(:, k) - xyz(:, j))
      delta = norm2(axis)
      if (delta <= geom_guard) then
         call env%error("cannot scan a linear angle", source)
         fail = .true.
         return
      end if
      call landscape_scan_fragments(at, xyz, i, j, moving)
      call landscape_scan_check_cut(env, moving, j, fail)
      if (fail) return
      call rotate_fragment(xyz, moving, xyz(:, j), axis/delta, alpha)
   case (p_landscape_scan_coord_dihedral)
      i = coord%atoms(1)
      j = coord%atoms(2)
      k = coord%atoms(3)
      l = coord%atoms(4)
      alpha = target - landscape_scan_coord_value(coord, xyz)
      axis = xyz(:, k) - xyz(:, j)
      delta = norm2(axis)
      if (delta <= geom_guard) then
         call env%error("cannot scan a dihedral with coincident central atoms", &
            & source)
         fail = .true.
         return
      end if
      ! The moving side is the one behind the central bond `j-k`, i.e. the
      ! fragment seeded at `k` (which contains `l` and every substituent of
      ! `k`).  Seeding at `l` instead cuts the wrong edge and moves a single
      ! atom: for `l` bonded to `k` the fragment collapses to `{l}`, the
      ! substituents of `k` stay behind and, as soon as the connectivity
      ! happens to bridge the cut, the fragment swallows the `i` side and the
      ! scanned dihedral stops changing.
      call landscape_scan_fragments(at, xyz, k, j, moving)
      call landscape_scan_check_cut(env, moving, j, fail)
      if (fail) return
      call rotate_fragment(xyz, moving, xyz(:, j), axis/delta, alpha)
   case default
      call env%error("unknown internal coordinate type", source)
      fail = .true.
      return
   end select

end subroutine landscape_scan_rigid_geometry

!> Refuse a cut that does not separate the molecule into two sides.
!>
!> `moving` is the fragment returned by `landscape_scan_fragments(seed,
!> blocked)` and `blocked` is the atom that must stay behind.  Whenever the
!> blocked atom is part of the fragment, the cut did not separate the
!> system: the two sides of a ring are still connected, the two atoms are
!> not bonded at all, or a third atom bridges the cut.  The rigid operation
!> would then move both sides (or the whole molecule), the scanned
!> coordinate would not change and the point would silently be wrong, so the
!> scan is stopped with an error instead.
subroutine landscape_scan_check_cut(env, moving, blocked, fail)
   character(len=*), parameter :: source = 'landscape_scan_check_cut'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Mask of the fragment
   logical, intent(in) :: moving(:)
   !> Atom that has to stay outside the fragment
   integer, intent(in) :: blocked
   !> Failure flag, set on return
   logical, intent(inout) :: fail

   if (.not. moving(blocked)) return
   call env%error("the bond cut of a rigid scan does not separate the two "//&
      & "sides of the molecule: atom "//nat_to_string(blocked)//" is part "//&
      & "of the moving fragment, so the scanned coordinate would not "//&
      & "change; a ring or a bridged system needs a relaxed scan", source)
   fail = .true.

end subroutine landscape_scan_check_cut

!> Fragment of the molecular graph that contains `seed` once the bond
!> `seed-blocked` is cut.
!>
!> The adjacency is built from the covalent radii of the atoms with the
!> `bond_factor` tolerance of `xtb_scanparam::constrain_all_bonds`.  The
!> returned mask has `size(xyz, 2)` entries; `blocked` itself is never part of
!> the fragment.
pure subroutine landscape_scan_fragments(at, xyz, seed, blocked, moving)
   !> Atom numbers (1-based)
   integer, intent(in) :: at(:)
   !> Cartesian geometry (bohr)
   real(wp), intent(in) :: xyz(:, :)
   !> Atom the fragment must contain
   integer, intent(in) :: seed
   !> Neighbour that is not crossed
   integer, intent(in) :: blocked
   !> Mask of the fragment
   logical, allocatable, intent(out) :: moving(:)

   integer, allocatable :: stack(:)
   integer :: nat, nstack, i, j
   real(wp) :: cutoff

   nat = size(xyz, 2)
   allocate(moving(nat), source=.false.)
   allocate(stack(nat), source=0)
   moving(seed) = .true.
   nstack = 1
   stack(1) = seed
   do while (nstack > 0)
      i = stack(nstack)
      nstack = nstack - 1
      do j = 1, nat
         if (moving(j)) cycle
         if (i == seed .and. j == blocked) cycle
         cutoff = bond_factor*(covalent_radius_2009(at(i)) + &
            & covalent_radius_2009(at(j)))
         if (norm2(xyz(:, j) - xyz(:, i)) > cutoff) cycle
         moving(j) = .true.
         nstack = nstack + 1
         stack(nstack) = j
      end do
   end do

end subroutine landscape_scan_fragments

!> Rotate the atoms of `moving` about the axis `axis` through `pivot`.
pure subroutine rotate_fragment(xyz, moving, pivot, axis, angle)
   !> Cartesian geometry (bohr), modified in place
   real(wp), intent(inout) :: xyz(:, :)
   !> Mask of the atoms to rotate
   logical, intent(in) :: moving(:)
   !> Point on the rotation axis (bohr)
   real(wp), intent(in) :: pivot(3)
   !> Unit rotation axis
   real(wp), intent(in) :: axis(3)
   !> Rotation angle (radian)
   real(wp), intent(in) :: angle

   real(wp) :: ca, sa, rel(3), rot(3)
   integer :: i

   ca = cos(angle)
   sa = sin(angle)
   do i = 1, size(xyz, 2)
      if (.not. moving(i)) cycle
      rel = xyz(:, i) - pivot
      rot = rel*ca + cross(axis, rel)*sa + axis*sum(axis*rel)*(1.0_wp - ca)
      xyz(:, i) = pivot + rot
   end do

end subroutine rotate_fragment

!> Cross product of two three-vectors.
pure function cross(a, b) result(c)
   !> First vector
   real(wp), intent(in) :: a(3)
   !> Second vector
   real(wp), intent(in) :: b(3)
   !> `a x b`
   real(wp) :: c(3)

   c(1) = a(2)*b(3) - a(3)*b(2)
   c(2) = a(3)*b(1) - a(1)*b(3)
   c(3) = a(1)*b(2) - a(2)*b(1)

end function cross

!> Decimal representation of a single digit index (1 or 2).
pure function one_digit(value) result(string)
   !> Index (at most two)
   integer, intent(in) :: value
   !> Decimal string
   character(len=:), allocatable :: string

   select case (value)
   case (1)
      string = "1"
   case (2)
      string = "2"
   case default
      string = "?"
   end select

end function one_digit

!> Decimal representation of an atom count for error messages.
function nat_to_string(value) result(string)
   !> Value to convert
   integer, intent(in) :: value
   !> Decimal string
   character(len=:), allocatable :: string

   character(len=16) :: buffer

   write(buffer, '(i0)') value
   string = trim(buffer)

end function nat_to_string

end module xtb_landscape_scan_relax
