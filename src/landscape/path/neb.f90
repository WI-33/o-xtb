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

!> Nudged elastic band (NEB) of the chain-of-states searches (PLAN section
!> 6.4).  This module is the pure mathematical kernel: it never touches a
!> `TMolecule` or a calculator, it evaluates a `class(landscape_energy_model)`
!> at the geometries of the band and returns the optimised band as plain
!> arrays.  The task layer (`xtb_landscape_tasks_ts`) supplies the model and
!> writes the files.
!>
!> Algorithm (Henkelman & Jonsson, J. Chem. Phys. **113**, 9978 (2000)):
!>
!>  * **Improved tangent** (Eqs. 8-11).  With `dm = R_i - R_{i-1}`,
!>    `dp = R_{i+1} - R_i` the tangent of image `i` is
!>    - the unit vector of `dp` if `E_{i-1} < E_i < E_{i+1}` (a
!>      monotonically rising stretch, the tangent points uphill),
!>    - the unit vector of `dm` if `E_{i-1} > E_i > E_{i+1}` (the mirrored
!>      case on a falling stretch),
!>    - otherwise the weighted average
!>      `tau ~ dps*max(|E_{i+1}-E_i|,|E_{i-1}-E_i|) + other*min(...)`, where
!>      the larger weight is attached to the vector of the **higher energy**
!>      neighbour (Eq. 11).  This is the energy based upwind rule: it keeps
!>      the tangent continuous around an extremum, where a naive `dp + dm`
!>      average would cancel the two contributions and produce a kink.
!>    The implementation falls back to `dp + dm` and finally to `e_1` when
!>    the weighted average degenerates, so the tangent is always a unit
!>    vector and never a NaN.
!>  * **NEB force** `F_i = F_i^perp + k (|R_{i+1}-R_i| - |R_i-R_{i-1}|) tau_i`
!>    (Eq. 12), i.e. the perpendicular part of the physical force plus a
!>    spring force **along the tangent only**, built from the **scalar**
!>    difference of the two image distances.  The original elastic band
!>    spring of the 1998 formulation adds the whole vector
!>    `k (R_{i+1} - 2 R_i + R_{i-1})`, which pulls the images off the path
!>    (its component perpendicular to the tangent is not removed); it is
!>    deliberately not implemented (PLAN section 6.4).
!>  * **Dynamic springs** (ORCA style):
!>    `k_i = k_max - (k_max - k_min) exp(-decay dE_i / dE_max)` with
!>    `dE_i = max(E_i - max(E_{i-1},E_{i+1}), 0)` and `dE_max = max_i dE_i`;
!>    the two endpoints get `k_max` (their springs are never used).  The
!>    literal formula puts the **lower** end of the window at a vanishing
!>    rise (`dE_i = 0` gives `exp(0) = 1`, i.e. `k_i = k_min`) and lifts the
!>    spring of the image that rises above **both** of its neighbours
!>    towards `k_max` (its rise is `dE_max`, so `k_i = k_max -
!>    (k_max - k_min) exp(-decay)`, e.g. `4.5e-2` hartree/bohr^2 for the
!>    calibrated `k_min = 1e-2`, `k_max = 1e-1`, `decay = 0.5`).  The
!>    ordering is therefore the opposite of the intuitive one ("the barrier
!>    top gets the weak spring because its position is set by the physics
!>    and a flat stretch gets the stiff one"): the formula is the contract,
!>    it is what the reference implementation's expression says, and the
!>    unit test `dynamic-springs` pins it value by value (report 35 section
!>    2.3 records the contradiction between that expression and the comment
!>    of the reference implementation).
!>
!>    Measured on the analytic LEPS surface of H + H2 with 9 images, 400
!>    iteration budget (report 37 sections 2.2 and 3.2, raw output in
!>    `m5-path-data/`): the dynamic schedule of the calibrated defaults
!>    converges in **48** iterations and the PLAN initial values
!>    (`k_min = 2.083e-3`, `k_max = 2.083e-2`) in **30**, both with a barrier
!>    error below `1e-4` kcal/mol against the closed form `8.9985` kcal/mol.
!>    **Static** springs in between do *not* converge within the budget
!>    (`1e-2`: barrier error `1.1e6` kcal/mol; `1e-1`: `2.3e4`), and a static
!>    `1e2` does not either (`5.8e5` kcal/mol, the VASP warning about an
!>    oversized spring constant); static `1e-5 ... 1e-3`, `1` and `10` do
!>    converge (31, 30, 31, 37 and 43 iterations).  On the molecular surface
!>    (GFN2 ammonia inversion) both static ends converge as well (11 and 10
!>    iterations), so the static failure window is a property of the steep
!>    analytic surface: the schedule, not a single constant, is what makes
!>    the defaults work.
!>  * **Convergence** uses geomeTRIC's dual criterion: the average and the
!>    largest of the per-image perpendicular force RMS must both be below
!>    their thresholds.  A single criterion converges "on average" while one
!>    image is still badly placed; the two together cannot.
!>  * **Alignment** (Kabsch, `remove_rotation_and_translation`): every
!>    internal image is rotated onto the reactant before each evaluation.
!>    The alignment only fixes the frame in which the geometries are
!>    compared, interpolated and printed; the two endpoints are fixed and
!>    are never moved by it.
!>  * **Fixed endpoints** by default (`relax_ends = .false.`): the first and
!>    the last image do not take part in the optimisation.  The optional
!>    endpoint relaxation is performed by the task with the M1 minimiser
!>    before the band is built.  Because the endpoints never move, their
!>    energy and gradient are evaluated **once** and reused by every later
!>    iteration (`landscape_band_endpoint_cache`); re-evaluating them would
!>    recompute the same pair of values, since the engine's `off` restart
!>    policy restores the recorded wavefunction before every evaluation.
!>    That equality is exact for `off` and `pin`; with `--lts-restart warm`
!>    it is not, see the note on the cache.
!>  * **Propagation and its measured failure window.**  The band is propagated
!>    by a truncated L-BFGS step on the band force (steps are scaled down, not
!>    projected, when the largest component exceeds `max_displ`) and the step
!>    is taken as it is: there is **no line search**, because the two natural
!>    merit functions are not descent functions of the propagation.  Measured
!>    on the GFN2 ammonia inversion (report 39 section 6.1): at the third
!>    iteration the step raises the band force norm at *every* halting of the
!>    step, for the L-BFGS direction and for the steepest descent direction
!>    alike, although the band is still a factor 2.8 above the force
!>    threshold; the tangential springs are what moves the images there.  On
!>    the analytic LEPS surface the norm grows by a factor 10 for several
!>    iterations in a row and a step that lowers the *energy* of the images
!>    can raise the force norm, so an energy based Armijo test rejects the
!>    steps that spread the band.  The consequence is that convergence is
!>    **not monotone in the number of images**: measured with the IDPP path
!>    and a 512 iteration budget, 3, 5, 9 and 11 images converge with the
!>    default bound of `0.37795` bohr while 7, 13 and 15 do not, and with
!>    `max_displ = 1e-1` bohr 3, 5, 7, 9 and 11 converge while 13 and 15 do
!>    not - the initial band force does not depend on the number of images
!>    (`f_rms_max = 1.199e-1` with IDPP, `5.687e-2` with the linear
!>    interpolation for 3 to 15 images).  Two families of failures are
!>    distinguished in the result: a band that **diverges** (the force norm
!>    grows beyond `band_divergence_factor` times its initial value) is
!>    aborted by the safety net with `diverged = .true.`, and a band that
!>    merely fails to meet the criterion within the budget is reported as
!>    unconverged by the iteration limit.  A third, weaker case sits between
!>    them: a band whose force norm rises in `band_descent_warn_streak`
!>    consecutive iterations is *reported* (once, naming `--lts-max-displ`)
!>    but neither moved nor aborted, because it may still recover - only the
!>    10^3 net is an abort.  All of these are failures of the search,
!>    never silent: a band that returns `converged = .true.` reproduces the
!>    analytic barrier (checked by the unit test `band-initial-path-scan` on
!>    every configuration of the scan).  AutoNEB inherits the window.
!>
!> Units: coordinates in bohr, energies in hartree, gradients and forces in
!> hartree/bohr, spring constants in hartree/bohr^2, band distances in bohr.
module xtb_landscape_path_neb
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : autokcal
   use xtb_type_environment, only : TEnvironment
   use xtb_landscape_settings, only : landscape_path_settings, landscape_value_present
   use xtb_landscape_engine, only : landscape_energy_model
   use xtb_landscape_coords, only : geometry_pack, geometry_unpack, align_geometry
   use xtb_landscape_lbfgs, only : landscape_lbfgs_type, new_landscape_lbfgs
   use xtb_landscape_path_idpp, only : landscape_idpp_path
   use xtb_landscape_path_cineb, only : landscape_cineb_force, &
      & landscape_cineb_ignition_ready, landscape_cineb_highest_image
   implicit none(type, external)
   private

   public :: landscape_band_type, landscape_band_endpoint_cache
   public :: landscape_band_optimize
   public :: landscape_band_build_path, landscape_band_energies, landscape_band_forces
   public :: landscape_band_align, landscape_band_measures, landscape_band_converged
   public :: landscape_neb_tangent, landscape_neb_springs, landscape_neb_force
   public :: landscape_band_pack, landscape_band_unpack
   public :: landscape_image_force_rms, landscape_band_distances

   !> Abort the propagation once the band force norm has grown beyond this
   !> factor of its initial value (the divergence safety net)
   real(wp), parameter :: band_divergence_factor = 1.0e3_wp

   !> Consecutive iterations of a band whose force norm rises before the
   !> propagation is reported as non-descending (H3 of report 220 section 4).
   !> The divergence safety net above only fires at 10^3 times the initial
   !> norm, so a *slowly* non-descending band used to spend the whole step
   !> budget and then report a bare "maximum number of steps".  Ten iterations
   !> is above the longest rise the calibrated runs show (the LEPS band with
   !> the climbing image grows by a factor 10 for a few iterations in a row)
   !> and short enough to be visible inside a default budget.
   integer, parameter :: band_descent_warn_streak = 10

   !> Result of a band optimisation.
   type :: landscape_band_type
      !> True if the dual convergence criterion was met
      logical :: converged = .false.
      !> Number of band iterations performed
      integer :: n_iter = 0
      !> Number of images (endpoints included)
      integer :: n_images = 0
      !> Number of energy/gradient evaluations performed by the kernel: the
      !> initial band (`n_images`), one sweep over the internal images per
      !> iteration (`n_images - 2`), because the two endpoints are evaluated
      !> once and cached
      integer :: n_eval = 0
      !> True if the propagation was aborted by the divergence safety net
      logical :: diverged = .false.
      !> Ignited climbing image (0 if the climbing image was never switched
      !> on, i.e. for a plain NEB run)
      integer :: climbing_image = 0
      !> Band iteration at which the climbing image was ignited (-1 if never)
      integer :: ignition_iteration = -1
      !> Highest energy internal image of the final band
      integer :: highest_image = 0
      !> Number of iterations after ignition in which the highest energy
      !> image was **not** the climbing image.  The climbing image is frozen
      !> at ignition (the reference behaviour), so a non-zero value is
      !> reported as a warning: it means the path topology moved under the
      !> climbing image.
      integer :: n_climb_switch = 0
      !> Average perpendicular force RMS over the non climbing internal
      !> images (hartree/bohr)
      real(wp) :: f_rms_avg = 0.0_wp
      !> Largest per-image perpendicular force RMS (hartree/bohr)
      real(wp) :: f_rms_max = 0.0_wp
      !> RMS of the **unprojected** force at the climbing image
      !> (hartree/bohr), zero without a climbing image
      real(wp) :: ci_f_rms = 0.0_wp
      !> Energy of the transition state candidate (hartree): the climbing
      !> image if one was ignited, otherwise the highest internal image
      real(wp) :: ts_energy = 0.0_wp
      !> Barrier of the final band relative to the reactant (hartree)
      real(wp) :: barrier = 0.0_wp
      !> Energy of every image of the final band (hartree)
      real(wp), allocatable :: energies(:)
      !> Distance between consecutive images (bohr), length `n_images - 1`
      real(wp), allocatable :: distances(:)
      !> Human readable termination message
      character(len=:), allocatable :: message
   end type landscape_band_type

!> Cached energy and gradient of the two fixed endpoint images of a band.
!>
!> `landscape_band_pack` moves only the internal images `2 .. n_images-1`, so
!> the two endpoints keep their geometry for the whole run; evaluating them
!> again in every iteration recomputes the same pair of values.  This type
!> holds them, and `landscape_band_energies` fills it on the first sweep and
!> reuses it afterwards.
!>
!> **The reuse changes no number for the `off` and `pin` restart policies.**
!> The engine's default restart policy `off` records the incoming
!> wavefunction and restores it before every evaluation (`xtb_landscape_engine`),
!> so a repeated geometry returns a bitwise identical energy and gradient; the
!> cached values are therefore the values the re-evaluation would have
!> produced (measured byte for byte, see report 223).  With
!> `--lts-restart warm` the engine keeps the incoming wavefunction instead
!> (`engine.f90`, `p_landscape_restart_warm`), so the second evaluation of the
!> same geometry can land on a different SCF solution - `AGENTS.md` section 6
!> measures a 3.9e-07 Eh/bohr gradient jump for that policy - and the cache
!> freezes the first sweep's endpoint values where a re-evaluation would have
!> returned slightly different ones.  A warm-restart band is therefore not
!> bitwise equal to the pre-cache run; the band kernel cannot ask the model
!> for its restart mode, so the difference is documented here rather than
!> avoided (report 224, finding 5).
!>
!> The cache belongs to **one** band: fill it from the first evaluation of a
!> band whose endpoints do not move.  `landscape_band_energies` enforces that
!> by keeping copies of the two endpoint geometries it filled the cache from
!> and dropping the cache when the band hands over different ones (or a band
!> with a different atom count), so a cache reused by a second band cannot
!> return the first band's values.  The two consumers are
!> `landscape_band_optimize` and the string/GSM optimisation, whose endpoints
!> are fixed as well, and both create it right before the run.
type :: landscape_band_endpoint_cache
   !> True once both endpoint evaluations are stored
   logical :: valid = .false.
   !> Endpoint energies (hartree), `(1)` reactant and `(2)` product
   real(wp) :: energy(2) = 0.0_wp
   !> Endpoint gradients (hartree/bohr), `(3, nat, 2)`
   real(wp), allocatable :: gradient(:, :, :)
   !> Geometries the cache was filled from (bohr), `(3, nat, 2)`; the cache is
   !> dropped when the band presents different endpoints
   real(wp), allocatable :: geometry(:, :, :)
end type landscape_band_endpoint_cache

contains

!> Build the initial path of a band.
!>
!> The interior images are placed by linear interpolation between the two
!> endpoints (the product is first Kabsch aligned to the reactant when
!> `settings%align` is set, so the interpolation happens in a common frame)
!> and, with `settings%ifidpp`, replaced by the IDPP path of
!> `xtb_landscape_path_idpp`.  The endpoints themselves are copied verbatim:
!> they are fixed images.
subroutine landscape_band_build_path(env, settings, reactant, product, images, fail)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Band settings
   type(landscape_path_settings), intent(in) :: settings
   !> Reactant geometry (bohr), shape `(3, nat)`
   real(wp), intent(in) :: reactant(:, :)
   !> Product geometry (bohr), shape `(3, nat)`
   real(wp), intent(in) :: product(:, :)
   !> Band images `(3, nat, n_images)`, allocated on return
   real(wp), allocatable, intent(out) :: images(:, :, :)
   !> Failure flag
   logical, intent(out) :: fail

   real(wp), allocatable :: ref(:, :), prod(:, :)
   real(wp) :: lam
   integer :: n_images, nat, k

   fail = .true.
   nat = size(reactant, 2)
   if (nat < 2) return
   if (size(reactant, 1) /= 3 .or. size(product, 1) /= 3) return
   if (size(product, 2) /= nat) return
   n_images = settings%n_images
   if (n_images < 3) return

   allocate(ref(3, nat), source=reactant)
   allocate(prod(3, nat), source=product)
   if (settings%align) call align_geometry(env, prod, ref)

   allocate(images(3, nat, n_images))
   do k = 1, n_images
      lam = real(k - 1, wp)/real(n_images - 1, wp)
      images(:, :, k) = (1.0_wp - lam)*ref + lam*prod
   end do
   images(:, :, 1) = ref
   images(:, :, n_images) = prod

   if (settings%ifidpp) then
      call landscape_idpp_path(env, images, fail)
      if (fail) return
   end if
   if (settings%align) call landscape_band_align(env, images)
   fail = .false.

end subroutine landscape_band_build_path

!> Rotate every internal image of a band onto the reactant.
!>
!> The reactant (image 1) is the reference frame and is never modified; the
!> endpoint images are not touched either.  Alignment is a geometry
!> operation only (a proper rotation about the arithmetic centroid), so it
!> changes no energy and no distance.
subroutine landscape_band_align(env, images)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Band images `(3, nat, n_images)`, internal images overwritten in place
   real(wp), intent(inout) :: images(:, :, :)

   real(wp), allocatable :: ref(:, :)
   integer :: k, n_images

   n_images = size(images, 3)
   if (n_images < 3) return
   allocate(ref(3, size(images, 2)), source=images(:, :, 1))
   do k = 2, n_images - 1
      call align_geometry(env, images(:, :, k), ref)
   end do

end subroutine landscape_band_align

!> Evaluate all images of a band: energy and Cartesian gradient.
!>
!> With the optional `endpoints` cache the two fixed endpoint images are
!> evaluated only on the first sweep and read from the cache afterwards, which
!> removes `2/n_images` of the evaluations of a band run (40 % at five images,
!> 22 % at nine) without changing any value: see the note on
!> `landscape_band_endpoint_cache`.  The cache is dropped on entry when it was
!> filled from other endpoint geometries, so a cache handed to a second band
!> cannot answer with the first band's values.  `n_evaluated` reports how many
!> images this call actually sent to the potential, so a caller can keep an
!> honest cost counter instead of assuming `n_images`.
subroutine landscape_band_energies(env, model, images, energies, gradients, fail, &
      & endpoints, n_evaluated)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Potential
   class(landscape_energy_model), intent(inout) :: model
   !> Band images `(3, nat, n_images)`
   real(wp), intent(in) :: images(:, :, :)
   !> Energy of every image (hartree), allocated on return
   real(wp), allocatable, intent(out) :: energies(:)
   !> Cartesian gradient of every image (hartree/bohr), `(3, nat, n)`
   real(wp), allocatable, intent(out) :: gradients(:, :, :)
   !> Failure flag
   logical, intent(out) :: fail
   !> Cached endpoint values of this band (optional)
   type(landscape_band_endpoint_cache), intent(inout), optional :: endpoints
   !> Number of images evaluated on this call (optional)
   integer, intent(out), optional :: n_evaluated

   real(wp), allocatable :: x(:), g(:)
   integer :: k, nat, n_images, endpoint
   logical :: filled(2), cached

   fail = .false.
   nat = size(images, 2)
   n_images = size(images, 3)
   allocate(energies(n_images), gradients(3, nat, n_images))
   allocate(x(3*nat), g(3*nat))
   if (present(n_evaluated)) n_evaluated = 0
   filled = .false.
   if (present(endpoints)) then
      ! the cache is only usable for the band it was filled from: a different
      ! atom count or different endpoint geometries mean the stored values
      ! belong to another band (report 224, finding 4)
      if (allocated(endpoints%gradient)) then
         if (size(endpoints%gradient, 1) /= 3 .or. &
               & size(endpoints%gradient, 2) /= nat .or. &
               & size(endpoints%gradient, 3) /= 2) then
            deallocate(endpoints%gradient)
            endpoints%valid = .false.
         end if
      end if
      if (endpoints%valid) then
         if (n_images < 2 .or. .not. allocated(endpoints%geometry)) then
            endpoints%valid = .false.
         else if (size(endpoints%geometry, 1) /= 3 .or. &
               & size(endpoints%geometry, 2) /= nat .or. &
               & size(endpoints%geometry, 3) /= 2) then
            endpoints%valid = .false.
         else if (any(endpoints%geometry(:, :, 1) /= images(:, :, 1)) .or. &
               & any(endpoints%geometry(:, :, 2) /= images(:, :, n_images))) then
            endpoints%valid = .false.
         end if
         if (.not. endpoints%valid) then
            if (allocated(endpoints%gradient)) deallocate(endpoints%gradient)
            if (allocated(endpoints%geometry)) deallocate(endpoints%geometry)
         end if
      end if
   end if
   do k = 1, n_images
      endpoint = 0
      if (k == 1) endpoint = 1
      if (k == n_images) endpoint = 2
      cached = .false.
      if (endpoint > 0 .and. present(endpoints)) cached = endpoints%valid
      if (cached) then
         energies(k) = endpoints%energy(endpoint)
         gradients(:, :, k) = endpoints%gradient(:, :, endpoint)
         cycle
      end if
      x = geometry_pack(images(:, :, k))
      call model%evaluate(x, energies(k), g, fail)
      if (fail) then
         call env%error("a band image could not be evaluated", &
            & 'landscape_band_energies')
         return
      end if
      ! A diverged image (two atoms on top of each other, an overflowing
      ! exponential) must not be propagated as a usable energy: the failure
      ! is reported and the band is abandoned, exactly like the NaN energy of
      ! an unconverged scan point.
      if (.not. landscape_value_present(energies(k)) .or. &
            & .not. all(landscape_value_present(g))) then
         call env%error("a band image returned a non-finite energy or "//&
            & "gradient", 'landscape_band_energies')
         fail = .true.
         return
      end if
      gradients(:, :, k) = reshape(g, [3, nat])
      if (present(n_evaluated)) n_evaluated = n_evaluated + 1
      if (endpoint > 0 .and. present(endpoints)) then
         if (.not. allocated(endpoints%gradient)) &
            & allocate(endpoints%gradient(3, nat, 2))
         if (.not. allocated(endpoints%geometry)) &
            & allocate(endpoints%geometry(3, nat, 2))
         endpoints%energy(endpoint) = energies(k)
         endpoints%gradient(:, :, endpoint) = gradients(:, :, k)
         endpoints%geometry(:, :, endpoint) = images(:, :, k)
         filled(endpoint) = .true.
      end if
   end do
   ! `valid` means both endpoints are stored; a call that found one of them
   ! cached keeps the flag it found
   if (present(endpoints)) endpoints%valid = endpoints%valid .or. all(filled)

end subroutine landscape_band_energies

!> Improved tangent of image `i` (Henkelman & Jonsson Eqs. 8-11).
!>
!> `xm`, `x` and `xp` are the flattened coordinates of the images
!> `i-1, i, i+1` and `em`, `e`, `ep` their energies.  The result `tau` is a
!> unit vector of the same length; it is normalized with the fallbacks
!> documented in the module header, so a degenerate input still yields a
!> usable direction instead of a NaN.
pure subroutine landscape_neb_tangent(xm, x, xp, em, e, ep, tau)
   !> Flattened coordinates of image `i-1` (bohr)
   real(wp), intent(in) :: xm(:)
   !> Flattened coordinates of image `i` (bohr)
   real(wp), intent(in) :: x(:)
   !> Flattened coordinates of image `i+1` (bohr)
   real(wp), intent(in) :: xp(:)
   !> Energy of image `i-1` (hartree)
   real(wp), intent(in) :: em
   !> Energy of image `i` (hartree)
   real(wp), intent(in) :: e
   !> Energy of image `i+1` (hartree)
   real(wp), intent(in) :: ep
   !> Unit tangent (dimensionless)
   real(wp), intent(out) :: tau(:)

   real(wp), allocatable :: dm(:), dp(:), um(:), up(:)
   real(wp) :: norm, rise_up, rise_down

   allocate(dm(size(x)), dp(size(x)), um(size(x)), up(size(x)))
   dm = x - xm
   dp = xp - x
   um = dm
   up = dp
   call unit_vector(um)
   call unit_vector(up)

   if (ep > e .and. e > em) then
      ! monotonically rising: upwind is the forward difference
      tau = up
   else if (ep < e .and. e < em) then
      ! monotonically falling: upwind is the backward difference
      tau = um
   else
      ! energy extremum (Eq. 11): the neighbour with the higher energy gets
      ! the larger of the two energy differences as its weight
      rise_up = abs(ep - e)
      rise_down = abs(em - e)
      if (ep >= em) then
         tau = max(rise_up, rise_down)*up + min(rise_up, rise_down)*um
      else
         tau = min(rise_up, rise_down)*up + max(rise_up, rise_down)*um
      end if
   end if

   norm = sqrt(sum(tau**2))
   if (norm < 1.0e-16_wp) then
      tau = dp + dm
      norm = sqrt(sum(tau**2))
   end if
   if (norm < 1.0e-16_wp) then
      tau = 0.0_wp
      if (size(tau) > 0) tau(1) = 1.0_wp
      norm = 1.0_wp
   end if
   tau = tau/norm

end subroutine landscape_neb_tangent

!> Normalize a vector in place; a zero vector is left as the first unit
!> vector, so the caller never divides by zero.
pure subroutine unit_vector(v)
   !> Vector, normalized in place
   real(wp), intent(inout) :: v(:)

   real(wp) :: norm

   norm = sqrt(sum(v**2))
   if (norm > 1.0e-16_wp) then
      v = v/norm
   else
      v = 0.0_wp
      if (size(v) > 0) v(1) = 1.0_wp
   end if

end subroutine unit_vector

!> Spring constants of a band, one per image (endpoints included).
!>
!> With `use_dynamic` the ORCA style formula
!> `k_i = k_max - (k_max - k_min) exp(-decay dE_i/dE_max)` is used, where
!> `dE_i` is the rise of image `i` above the higher of its two neighbours and
!> `dE_max` the largest such rise; the two endpoints get `k_max` (they are
!> not used).  Without `use_dynamic` every spring is `k_max`.  A flat band
!> (`dE_max` at the level of the round off) falls back to `k_max` everywhere,
!> which is the well conditioned limit.
pure subroutine landscape_neb_springs(energies, k_min, k_max, decay, &
      & use_dynamic, springs)
   !> Image energies (hartree), endpoints included
   real(wp), intent(in) :: energies(:)
   !> Lower spring constant (hartree/bohr^2)
   real(wp), intent(in) :: k_min
   !> Upper spring constant (hartree/bohr^2)
   real(wp), intent(in) :: k_max
   !> Decay factor of the dynamic springs (dimensionless)
   real(wp), intent(in) :: decay
   !> Use the dynamic formula
   logical, intent(in) :: use_dynamic
   !> Spring constant per image (hartree/bohr^2)
   real(wp), intent(out) :: springs(:)

   real(wp), allocatable :: rise(:)
   real(wp) :: rise_max
   integer :: i, n

   n = size(energies)
   springs = k_max
   if (.not. use_dynamic .or. n < 3) return

   allocate(rise(n))
   rise = 0.0_wp
   do i = 2, n - 1
      rise(i) = max(energies(i) - max(energies(i - 1), energies(i + 1)), 0.0_wp)
   end do
   rise_max = maxval(rise)
   if (rise_max < 1.0e-12_wp) return

   do i = 2, n - 1
      springs(i) = k_max - (k_max - k_min)*exp(-decay*rise(i)/rise_max)
   end do

end subroutine landscape_neb_springs

!> NEB force of one image: perpendicular physical force plus the spring term
!> of Eq. 12.
pure subroutine landscape_neb_force(force_perp, tau, dm, dp, k, force)
   !> Perpendicular physical force (hartree/bohr)
   real(wp), intent(in) :: force_perp(:)
   !> Unit tangent (dimensionless)
   real(wp), intent(in) :: tau(:)
   !> Displacement to the previous image `R_i - R_{i-1}` (bohr)
   real(wp), intent(in) :: dm(:)
   !> Displacement to the next image `R_{i+1} - R_i` (bohr)
   real(wp), intent(in) :: dp(:)
   !> Spring constant (hartree/bohr^2)
   real(wp), intent(in) :: k
   !> NEB force (hartree/bohr)
   real(wp), intent(out) :: force(:)

   real(wp) :: d_next, d_prev

   d_next = sqrt(sum(dp**2))
   d_prev = sqrt(sum(dm**2))
   force = force_perp + k*(d_next - d_prev)*tau

end subroutine landscape_neb_force

!> Band forces of every image.
!>
!> Internal images get the NEB force; with `climbing_image` set, that image
!> gets the climbing image force of Eq. 5 instead (no spring, reversed
!> parallel component).  The two endpoints never move; their entry holds the
!> **physical** force, so the printed path table shows the residual gradient
!> of the endpoints, and `landscape_band_measures` deliberately skips them.
!>
!> The frozen spring constants come from `landscape_neb_springs` and are
!> recomputed from the current energies, i.e. the springs are dynamic as soon
!> as `settings%use_dynamic_k` is set.
subroutine landscape_band_forces(env, images, energies, gradients, settings, &
      & climbing_image, forces)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Band images `(3, nat, n_images)`
   real(wp), intent(in) :: images(:, :, :)
   !> Image energies (hartree)
   real(wp), intent(in) :: energies(:)
   !> Image gradients (hartree/bohr)
   real(wp), intent(in) :: gradients(:, :, :)
   !> Band settings
   type(landscape_path_settings), intent(in) :: settings
   !> Climbing image (0 for a plain NEB)
   integer, intent(in) :: climbing_image
   !> Band force of every image (hartree/bohr), allocated on return
   real(wp), allocatable, intent(out) :: forces(:, :, :)

   real(wp), allocatable :: springs(:), tau(:), xm(:), x(:), xp(:), f(:), fp(:), fci(:)
   real(wp), allocatable :: dm(:), dp(:)
   real(wp) :: parallel
   integer :: k, n_images, nat, ndim

   nat = size(images, 2)
   n_images = size(images, 3)
   ndim = 3*nat
   allocate(forces(3, nat, n_images), source=0.0_wp)

   allocate(springs(n_images), tau(ndim), xm(ndim), x(ndim), xp(ndim))
   allocate(f(ndim), fp(ndim), fci(ndim), dm(ndim), dp(ndim))
   call landscape_neb_springs(energies, settings%k_min, settings%k_max, &
      & settings%k_decay, settings%use_dynamic_k, springs)

   do k = 1, n_images
      if (k == 1 .or. k == n_images) then
         ! fixed endpoint: report the physical force, never move it
         forces(:, :, k) = -gradients(:, :, k)
         cycle
      end if
      xm = geometry_pack(images(:, :, k - 1))
      x = geometry_pack(images(:, :, k))
      xp = geometry_pack(images(:, :, k + 1))
      call landscape_neb_tangent(xm, x, xp, energies(k - 1), energies(k), &
         & energies(k + 1), tau)
      f = -reshape(gradients(:, :, k), [ndim])

      if (k == climbing_image) then
         call landscape_cineb_force(f, tau, fci)
         forces(:, :, k) = reshape(fci, [3, nat])
         cycle
      end if

      parallel = dot_product(f, tau)
      fp = f - parallel*tau
      dm = x - xm
      dp = xp - x
      call landscape_neb_force(fp, tau, dm, dp, springs(k), fci)
      forces(:, :, k) = reshape(fci, [3, nat])
   end do

end subroutine landscape_band_forces

!> Per-image perpendicular force RMS and its band average and maximum.
!>
!> Only the **internal** images contribute (the endpoints are fixed) and the
!> climbing image is excluded: its convergence is measured with the
!> unprojected force instead, because the climbing image force of Eq. 5 has
!> no spring term and therefore no convergent perpendicular part in the same
!> sense.  The per-image value is `sqrt(sum(F^2)/(3 nat))`, so the measure is
!> a component RMS like the gradient criteria of the minimiser.
pure subroutine landscape_band_measures(forces, climbing_image, f_rms_avg, f_rms_max)
   !> Band forces (hartree/bohr)
   real(wp), intent(in) :: forces(:, :, :)
   !> Climbing image, excluded from the measures (0 for a plain NEB)
   integer, intent(in) :: climbing_image
   !> Average of the per-image force RMS (hartree/bohr)
   real(wp), intent(out) :: f_rms_avg
   !> Largest per-image force RMS (hartree/bohr)
   real(wp), intent(out) :: f_rms_max

   real(wp) :: rms
   integer :: k, n_images, ncomp, ncount

   n_images = size(forces, 3)
   ncomp = 3*size(forces, 2)
   f_rms_avg = 0.0_wp
   f_rms_max = 0.0_wp
   ncount = 0
   do k = 2, n_images - 1
      if (k == climbing_image) cycle
      rms = sqrt(sum(forces(:, :, k)**2)/real(ncomp, wp))
      f_rms_avg = f_rms_avg + rms
      f_rms_max = max(f_rms_max, rms)
      ncount = ncount + 1
   end do
   if (ncount > 0) f_rms_avg = f_rms_avg/real(ncount, wp)

end subroutine landscape_band_measures

!> Dual band convergence criterion (geomeTRIC style).
!>
!> Both the average and the largest per-image force RMS must be below their
!> thresholds; with a climbing image the **unprojected** force RMS at that
!> image must meet the maximum threshold as well.  The `AND` is the point:
!> an average criterion alone declares a band converged while a single badly
!> placed image is still far from the path.
pure function landscape_band_converged(f_rms_avg, f_rms_max, ci_f_rms, &
      & has_climbing, settings) result(converged)
   !> Average per-image force RMS (hartree/bohr)
   real(wp), intent(in) :: f_rms_avg
   !> Largest per-image force RMS (hartree/bohr)
   real(wp), intent(in) :: f_rms_max
   !> Unprojected force RMS at the climbing image (hartree/bohr)
   real(wp), intent(in) :: ci_f_rms
   !> True if a climbing image is active
   logical, intent(in) :: has_climbing
   !> Band settings
   type(landscape_path_settings), intent(in) :: settings
   !> True if every criterion is fulfilled
   logical :: converged

   converged = f_rms_avg <= settings%f_rms_avg_tol &
      & .and. f_rms_max <= settings%f_rms_max_tol
   if (has_climbing) converged = converged .and. ci_f_rms <= settings%f_rms_max_tol

end function landscape_band_converged

!> RMS of the unprojected force at one image (hartree/bohr).
pure function landscape_image_force_rms(gradient) result(rms)
   !> Cartesian gradient of the image (hartree/bohr), `(3, nat)`
   real(wp), intent(in) :: gradient(:, :)
   !> Component RMS (hartree/bohr)
   real(wp) :: rms

   rms = sqrt(sum(gradient**2)/real(3*size(gradient, 2), wp))

end function landscape_image_force_rms

!> Flatten the internal images of a band into one vector.
pure subroutine landscape_band_pack(images, x)
   !> Band images `(3, nat, n_images)`
   real(wp), intent(in) :: images(:, :, :)
   !> Packed internal images (bohr), length `3 nat (n_images - 2)`
   real(wp), intent(out) :: x(:)

   real(wp), allocatable :: one(:)
   integer :: k, nat, offset

   nat = size(images, 2)
   allocate(one(3*nat))
   offset = 0
   do k = 2, size(images, 3) - 1
      one = geometry_pack(images(:, :, k))
      x(offset + 1:offset + 3*nat) = one
      offset = offset + 3*nat
   end do

end subroutine landscape_band_pack

!> Unflatten the internal images of a band from one vector.
pure subroutine landscape_band_unpack(x, images, fail)
   !> Packed internal images (bohr)
   real(wp), intent(in) :: x(:)
   !> Band images `(3, nat, n_images)`, internal images overwritten
   real(wp), intent(inout) :: images(:, :, :)
   !> Failure flag
   logical, intent(out) :: fail

   integer :: k, nat, offset

   fail = .false.
   nat = size(images, 2)
   if (size(x) /= 3*nat*(size(images, 3) - 2)) then
      fail = .true.
      return
   end if
   offset = 0
   do k = 2, size(images, 3) - 1
      call geometry_unpack(x(offset + 1:offset + 3*nat), images(:, :, k))
      offset = offset + 3*nat
   end do

end subroutine landscape_band_unpack

!> Optimise a band with NEB and, optionally, the climbing image refinement.
!>
!> On entry `reactant` and `product` are the two endpoint geometries (same
!> atoms in the same order); on return `images` holds the converged band,
!> `energies` its energies and `forces` its band forces.  The two endpoints
!> are copied verbatim and never modified, so their energy and gradient are
!> evaluated once and reused (see `landscape_band_endpoint_cache`);
!> `band%n_eval` counts the evaluations that actually happened,
!> `n_images + n_iter*(n_images - 2)` for a run that reaches `n_iter`
!> iterations.  `band%converged` reports the
!> dual criterion; an exhausted iteration budget or a failed image evaluation
!> is reported through `band%message` and `band%converged = .false.`, so the
!> task layer can escalate it to a non-zero exit status.
!>
!> `initial_images` replaces the constructed path: the given band is relaxed
!> instead.  It is the entry point of the incremental strategies - the
!> stages of AutoNEB and the optimisation phase of the string method - and of
!> the restart reuse of an earlier path (the endpoints and the atom count of
!> the supplied band are validated, the path itself is not rebuilt).  The
!> images are Kabsch aligned exactly like a freshly built path.
subroutine landscape_band_optimize(env, model, settings, reactant, product, &
      & images, energies, forces, band, initial_images)
   character(len=*), parameter :: source = 'landscape_band_optimize'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Potential
   class(landscape_energy_model), intent(inout) :: model
   !> Band settings
   type(landscape_path_settings), intent(in) :: settings
   !> Reactant geometry (bohr), shape `(3, nat)`
   real(wp), intent(in) :: reactant(:, :)
   !> Product geometry (bohr), shape `(3, nat)`
   real(wp), intent(in) :: product(:, :)
   !> Band images `(3, nat, n_images)` of the converged band
   real(wp), allocatable, intent(out) :: images(:, :, :)
   !> Energy of every image (hartree), allocated once the band was evaluated
   real(wp), allocatable, intent(out) :: energies(:)
   !> Band force of every image (hartree/bohr)
   real(wp), allocatable, intent(out) :: forces(:, :, :)
   !> Band result
   type(landscape_band_type), intent(out) :: band
   !> Path to relax instead of building one (bohr), shape `(3, nat, n)`
   real(wp), intent(in), optional :: initial_images(:, :, :)

   type(landscape_lbfgs_type) :: lbfgs
   type(landscape_band_endpoint_cache) :: endpoints
   real(wp), allocatable :: trial(:, :, :), x(:), x_new(:), g(:), g_new(:)
   real(wp), allocatable :: d(:), en_new(:), grad_new(:, :, :), f_new(:, :, :)
   real(wp) :: ci_rms, merit_initial
   !> Force norm of the previous iteration and of the current one (hartree/bohr)
   real(wp) :: prev_force_norm, force_norm
   !> Number of consecutive iterations with a rising force norm
   integer :: n_rise
   !> True once the non-descending propagation was reported
   logical :: descent_warned
   !> Formatting buffer of the descent report (one line, ~300 characters)
   character(len=400) :: buffer
   integer :: iter, n_images, nat, ndim, highest, n_evaluated
   logical :: fail, ok, steepest, converged, has_climbing, aborted

   band%converged = .false.
   call band_settings_check(env, settings, fail)
   if (fail) then
      band%message = "invalid band settings"
      return
   end if
   if (present(initial_images)) then
      call adopt_initial_images(env, settings, reactant, product, initial_images, &
         & images, fail)
      if (fail) then
         band%message = "the supplied path is not a band between the endpoints"
         return
      end if
   else
      call landscape_band_build_path(env, settings, reactant, product, images, fail)
      if (fail) then
         call env%error("the band path could not be built", source)
         band%message = "band path construction failed"
         return
      end if
   end if

   nat = size(images, 2)
   n_images = size(images, 3)
   ndim = 3*nat*(n_images - 2)
   band%n_images = n_images

   allocate(energies(n_images))
   allocate(trial(3, nat, n_images))
   allocate(x(ndim), x_new(ndim), g(ndim), g_new(ndim), d(ndim))
   call landscape_band_energies(env, model, images, energies, grad_new, fail, &
      & endpoints, n_evaluated)
   if (fail) then
      band%message = "a band image could not be evaluated"
      return
   end if
   band%n_eval = n_evaluated

   call landscape_band_forces(env, images, energies, grad_new, settings, 0, forces)
   call landscape_band_measures(forces, 0, band%f_rms_avg, band%f_rms_max)
   merit_initial = band_force_norm(forces)
   call landscape_band_pack(images, x)
   call new_landscape_lbfgs(lbfgs, ndim, settings%lbfgs_memory, settings%curvature)

   if (settings%printlevel > 0) then
      write(env%unit, '(1x,a,i0,a,i0,a,i0,a)') "band: ", n_images, &
         & " images (", n_images - 2, " internal, ", 2, " fixed endpoints)"
      write(env%unit, '(1x,a,es10.3,a,es10.3,a)') "band: k_min ", settings%k_min, &
         & " Eh/bohr^2, k_max ", settings%k_max, " Eh/bohr^2"
      if (settings%use_dynamic_k) then
         write(env%unit, '(1x,a,f8.3)') "band: dynamic springs, decay ", &
            & settings%k_decay
      else
         write(env%unit, '(1x,a)') "band: static springs (k_max everywhere)"
      end if
      if (settings%climb) then
         write(env%unit, '(1x,a,es10.3,a)') "band: climbing image switch at ", &
            & settings%climb_switch_rms, " Eh/bohr"
      else
         write(env%unit, '(1x,a)') "band: plain NEB (climbing image disabled)"
      end if
      write(env%unit, '(1x,a)') "   iter   f_rms_avg     f_rms_max     dE(hei)/kcal"
      ! the iter 0 row is the state the ignition switch of the climbing image
      ! compares against, so it is part of the table (report 37 F3)
      highest = landscape_cineb_highest_image(energies)
      write(env%unit, '(1x,i6,1x,es12.5,1x,es12.5,1x,f14.4)') 0, &
         & band%f_rms_avg, band%f_rms_max, (energies(max(highest, 1)) - energies(1))*autokcal
   end if

   iter = 0
   converged = .false.
   aborted = .false.
   ! H3 of report 220 section 4: the descent counter starts at the norm of the
   ! initial band, so iteration 1 is compared against the state the run began
   ! in.  A pure counter: it changes no step, no criterion and no exit.
   prev_force_norm = merit_initial
   n_rise = 0
   descent_warned = .false.
   do
      has_climbing = band%climbing_image > 0
      ci_rms = 0.0_wp
      if (has_climbing) then
         ci_rms = landscape_image_force_rms(grad_new(:, :, band%climbing_image))
      end if

      ! ignition: plain NEB first, the climbing image only once the band has
      ! relaxed below the switch threshold
      if (settings%climb .and. .not. has_climbing) then
         if (landscape_cineb_ignition_ready(band%f_rms_max, settings%climb_switch_rms)) then
            band%climbing_image = landscape_cineb_highest_image(energies)
            band%ignition_iteration = iter
            call landscape_band_forces(env, images, energies, grad_new, settings, &
               & band%climbing_image, forces)
            call landscape_band_measures(forces, band%climbing_image, &
               & band%f_rms_avg, band%f_rms_max)
            has_climbing = .true.
            ci_rms = landscape_image_force_rms(grad_new(:, :, band%climbing_image))
            if (settings%printlevel > 0) then
               write(env%unit, '(1x,a,i0,a,i0)') "band: climbing image ", &
                  & band%climbing_image, " ignited at iteration ", iter
            end if
         end if
      end if

      converged = landscape_band_converged(band%f_rms_avg, band%f_rms_max, &
         & ci_rms, has_climbing, settings)
      if (converged) exit
      if (iter >= settings%max_steps) exit

      g = -reshape(pack_internal_vector(forces), [ndim])
      call lbfgs%direction(g, d, ok, steepest)
      if (.not. ok) exit
      call clip_band_step(d, settings%max_displ)
      x_new = x + d
      trial = images
      call landscape_band_unpack(x_new, trial, fail)
      if (fail) then
         band%message = "the band step could not be applied"
         aborted = .true.
         exit
      end if
      if (settings%align) call landscape_band_align(env, trial)
      call landscape_band_pack(trial, x_new)

      call landscape_band_energies(env, model, trial, en_new, grad_new, fail, &
         & endpoints, n_evaluated)
      if (fail) then
         call env%error("a band image could not be evaluated", source)
         band%message = "a band image could not be evaluated"
         aborted = .true.
         exit
      end if
      band%n_eval = band%n_eval + n_evaluated
      call landscape_band_forces(env, trial, en_new, grad_new, settings, &
         & band%climbing_image, f_new)
      call landscape_band_measures(f_new, band%climbing_image, &
         & band%f_rms_avg, band%f_rms_max)

      g_new = -reshape(pack_internal_vector(f_new), [ndim])
      call lbfgs%push(x_new - x, g_new - g)

      images = trial
      energies = en_new
      forces = f_new
      x = x_new
      iter = iter + 1
      band%n_iter = iter

      ! Divergence safety net of the propagation.  The step is a truncated
      ! L-BFGS step without a line search (the merit functions that a line
      ! search could use are not monotone along the propagation, see the
      ! module documentation), so an overshooting step sends the band away
      ! from the surface and the force norm explodes: measured, the band force
      ! norm of the 11 image LEPS band grew to 1.3e4 times its initial value
      ! and the barrier to 6.7e5 kcal/mol while the criterion was never met.
      ! Running out the whole iteration budget on such a band is useless, so
      ! the run is aborted as soon as the norm exceeds
      ! `band_divergence_factor` times its initial value.  The threshold is
      ! two orders of magnitude above the largest transient of the calibrated
      ! runs (measured: 4.4e1 on the LEPS band with the climbing image,
      ! report 37 section 3.1(d)).
      force_norm = band_force_norm(forces)

      ! H3 of report 220 section 4: the net above only fires at 10^3 times the
      ! initial norm, so a *slowly* non-descending band burns the whole step
      ! budget and reports a bare "maximum number of steps" - the user cannot
      ! tell a propagation that was still improving from one that never
      ! descended.  The counter reports the latter once.  It is deliberately a
      ! **warning and not an abort**: the calibrated bands show rises of a few
      ! iterations in a row before recovering (report 37 section 3.1(d)), so an
      ! abort here would kill a band that would have converged - which is worse
      ! than the loud slow run, because a silent wrong answer is worse than
      ! both.  The counter is state only: no step, criterion or exit changes.
      if (force_norm > prev_force_norm) then
         n_rise = n_rise + 1
      else
         n_rise = 0
      end if
      prev_force_norm = force_norm
      if (n_rise >= band_descent_warn_streak .and. .not. descent_warned) then
         descent_warned = .true.
         write(buffer, '(a,i0,a,es9.3,a)') "the band force norm rose in each of "//&
            & "the last ", n_rise, " iterations (--lts-max-displ ", &
            & settings%max_displ, " bohr): the propagation is not descending "//&
            & "and can spend the whole step budget without converging; the "//&
            & "knob to look at is --lts-max-displ (the run continues)"
         call env%warning(trim(buffer), source)
      end if

      if (force_norm > band_divergence_factor*merit_initial) then
         band%diverged = .true.
         aborted = .true.
         band%message = "the band diverged from the path: the band force norm "//&
            & "grew beyond 1e3 times its initial value (reduce --lts-max-displ)"
         call env%error("the band diverged from the path; change --lts-max-displ "//&
            & "(the measured failure window of the propagation is documented "//&
            & "on xtb_landscape_path_neb)", source)
         exit
      end if

      if (settings%hessian_recalc > 0) then
         if (mod(iter, settings%hessian_recalc) == 0) call lbfgs%reset()
      end if

      highest = landscape_cineb_highest_image(energies)
      if (band%climbing_image > 0 .and. highest /= band%climbing_image) then
         band%n_climb_switch = band%n_climb_switch + 1
         if (band%n_climb_switch == 1) then
            call env%warning("the highest energy image moved away from the "//&
               & "climbing image after ignition; the climbing image stays "//&
               & "frozen at its ignition index", source)
         end if
      end if

      if (settings%printlevel > 0) then
         write(env%unit, '(1x,i6,1x,es12.5,1x,es12.5,1x,f14.4)') iter, &
            & band%f_rms_avg, band%f_rms_max, &
            & (energies(max(highest, 1)) - energies(1))*autokcal
      end if
   end do

   ! single exit: every field of the result is filled from the state of the
   ! band, whatever way the loop was left, so a caller can rely on the arrays
   ! and on the indices even for an aborted run
   band%converged = converged
   call band_finalise(band, images, energies, grad_new)
   if (converged) then
      band%message = "band converged"
   else if (aborted) then
      ! message already set
   else if (iter >= settings%max_steps) then
      band%message = "band reached the iteration limit without converging"
   else
      band%message = "band optimisation stopped"
   end if

   if (settings%printlevel > 0) then
      write(env%unit, '(1x,a,i0,a,es12.5,a,es12.5)') "band: iterations ", band%n_iter, &
         & ", f_rms_avg ", band%f_rms_avg, ", f_rms_max ", band%f_rms_max
      ! the step bound is repeated here: a run that did not converge has to
      ! name the knob that governs its robustness without the user having to
      ! scroll back to the header (report 37 section 8.1)
      write(env%unit, '(1x,a,i0,a,es12.4,a,i0,a,l1)') "band: ", band%n_images, &
         & " images, max_displ ", settings%max_displ, " bohr, max_steps ", &
         & settings%max_steps, ", diverged ", band%diverged
      write(env%unit, '(1x,a,i0,a,f20.12,a,f12.4,a)') "band: highest image ", &
         & band%highest_image, ", energy ", band%ts_energy, " Eh, barrier ", &
         & band%barrier*autokcal, " kcal/mol"
      if (band%climbing_image > 0) then
         write(env%unit, '(1x,a,i0,a,es12.5,a)') "band: climbing image ", &
            & band%climbing_image, ", force RMS ", band%ci_f_rms, " Eh/bohr"
      end if
   end if

end subroutine landscape_band_optimize

!> Validate the settings the propagation depends on.
!>
!> The dynamic spring formula stays inside the documented window
!> `k_min <= k_i <= k_max` only for a non negative decay: a negative decay
!> turns `exp(-decay dE/dE_max)` into a growth factor and can push the spring
!> above `k_max` (report 38 L-1).  An inverted window or a non positive
!> displacement bound would leave the band without a usable step, and a
!> non positive spring constant would let the images collapse onto each
!> other.  The check is performed here and not in the CLI layer because the
!> kernel is public and is called directly by the unit tests and by AutoNEB.
subroutine band_settings_check(env, settings, fail)
   character(len=*), parameter :: source = 'landscape_band_optimize'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Band settings
   type(landscape_path_settings), intent(in) :: settings
   !> Failure flag
   logical, intent(out) :: fail

   fail = .true.
   if (settings%max_displ <= 0.0_wp) then
      call env%error("the band displacement bound max_displ (--lts-max-displ) "//&
         & "must be positive", source)
      return
   end if
   if (settings%k_min <= 0.0_wp .or. settings%k_max <= 0.0_wp) then
      call env%error("the band spring constants k_min and k_max must be "//&
         & "positive", source)
      return
   end if
   if (settings%k_min > settings%k_max) then
      call env%error("the band spring window is inverted: k_min must not "//&
         & "exceed k_max", source)
      return
   end if
   if (settings%use_dynamic_k .and. settings%k_decay < 0.0_wp) then
      call env%error("the decay of the dynamic springs k_decay must not be "//&
         & "negative; a negative decay exceeds the documented window "//&
         & "k_min <= k_i <= k_max", source)
      return
   end if
   fail = .false.

end subroutine band_settings_check

!> Validate and adopt a supplied band path.
!>
!> The supplied path is the initial configuration of an incremental strategy
!> (an AutoNEB stage, the restart reuse of an AutoNEB run, the optimisation
!> phase of the string) or of a stage of the string method.  It is accepted
!> only if it is a genuine band between exactly the two endpoints: 3
!> dimensions, at least three images, a matching number of atoms and endpoints
!> that agree with the requested ones.
!>
!> The endpoint comparison is **rigid body invariant**: the pairwise distances
!> of the supplied endpoint are compared with those of the requested one
!> instead of the coordinates.  A coordinate comparison would reject the
!> path the driver itself produced, because with `settings%align` the
!> constructed band stores the product **rotated into the frame of the
!> reactant** (`landscape_band_build_path`); a Kabsch based comparison would
!> work as well but is unreliable for geometries with a degenerate rotation
!> (the two atom embedding of the unit tests), so the distance measure is the
!> right invariant here.
!>
!> The endpoints of the adopted path are set to the requested geometries:
!> image 1 becomes the reactant verbatim and image `n` the product, rotated
!> into the frame of the reactant when `settings%align` is set, exactly like a
!> freshly built path.  With the alignment off the frames are identical by
!> construction, so the endpoints are bitwise stable across the stages of an
!> incremental search.
subroutine adopt_initial_images(env, settings, reactant, product, initial_images, &
      & images, fail)
   character(len=*), parameter :: source = 'landscape_band_optimize'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Band settings (the alignment switch)
   type(landscape_path_settings), intent(in) :: settings
   !> Reactant geometry (bohr), shape `(3, nat)`
   real(wp), intent(in) :: reactant(:, :)
   !> Product geometry (bohr), shape `(3, nat)`
   real(wp), intent(in) :: product(:, :)
   !> Supplied path `(3, nat, n_images)`
   real(wp), intent(in) :: initial_images(:, :, :)
   !> Band images, allocated on return
   real(wp), allocatable, intent(out) :: images(:, :, :)
   !> Failure flag
   logical, intent(out) :: fail

   real(wp), allocatable :: dist_supplied(:), dist_wanted(:)
   real(wp), allocatable :: frame(:, :)
   integer :: nat, n_images

   fail = .true.
   if (size(reactant, 1) /= 3 .or. size(product, 1) /= 3) return
   nat = size(reactant, 2)
   if (size(product, 2) /= nat) return
   if (size(initial_images, 1) /= 3) return
   if (size(initial_images, 2) /= nat) return
   n_images = size(initial_images, 3)
   if (n_images < 3) return

   call pair_distance_vector(initial_images(:, :, 1), dist_supplied)
   call pair_distance_vector(reactant, dist_wanted)
   if (maxval(abs(dist_supplied - dist_wanted)) > 1.0e-4_wp) then
      call env%error("the supplied path starts at a geometry that is not the "//&
         & "reactant", source)
      return
   end if
   call pair_distance_vector(initial_images(:, :, n_images), dist_supplied)
   call pair_distance_vector(product, dist_wanted)
   if (maxval(abs(dist_supplied - dist_wanted)) > 1.0e-4_wp) then
      call env%error("the supplied path ends at a geometry that is not the "//&
         & "product", source)
      return
   end if

   allocate(images(3, nat, n_images), source=initial_images)
   images(:, :, 1) = reactant
   frame = product
   if (settings%align) then
      ! the same frame convention as `landscape_band_build_path`
      call align_geometry(env, frame, images(:, :, 1))
      call landscape_band_align(env, images)
   end if
   images(:, :, n_images) = frame
   fail = .false.

end subroutine adopt_initial_images

!> Pairwise distance vector of a geometry (bohr), upper triangle `a < b`.
subroutine pair_distance_vector(xyz, dist)
   !> Geometry (bohr), shape `(3, nat)`
   real(wp), intent(in) :: xyz(:, :)
   !> Pair distances (bohr), allocated on return
   real(wp), allocatable, intent(out) :: dist(:)

   integer :: a, b, nat, p

   nat = size(xyz, 2)
   allocate(dist(nat*(nat - 1)/2))
   p = 0
   do a = 1, nat
      do b = a + 1, nat
         p = p + 1
         dist(p) = sqrt(sum((xyz(:, a) - xyz(:, b))**2))
      end do
   end do

end subroutine pair_distance_vector

!> Fill the summary fields of a band result from the state of the band.
!>
!> Called on the single exit of `landscape_band_optimize`, so an aborted or
!> non converged band still carries its energies, distances, highest image
!> and climbing image; only a band that never produced an energy (a failed
!> first evaluation) leaves `band%energies` unallocated, and the callers test
!> for that.
subroutine band_finalise(band, images, energies, gradients)
   !> Band result, summary fields updated
   type(landscape_band_type), intent(inout) :: band
   !> Band images
   real(wp), intent(in) :: images(:, :, :)
   !> Image energies (hartree)
   real(wp), intent(in) :: energies(:)
   !> Image gradients (hartree/bohr)
   real(wp), intent(in) :: gradients(:, :, :)

   integer :: k, n_images

   n_images = size(images, 3)
   band%n_images = n_images
   if (allocated(band%energies)) deallocate(band%energies)
   allocate(band%energies(n_images), source=energies)
   call landscape_band_distances(images, band%distances)
   band%highest_image = landscape_cineb_highest_image(energies)
   k = band%highest_image
   if (band%climbing_image > 0) k = band%climbing_image
   if (k > 0) then
      band%ts_energy = energies(k)
      band%barrier = energies(k) - energies(1)
   end if
   band%ci_f_rms = 0.0_wp
   if (band%climbing_image > 0) then
      band%ci_f_rms = landscape_image_force_rms(gradients(:, :, band%climbing_image))
   end if

end subroutine band_finalise

!> Pack the internal band forces into one vector.
function pack_internal_vector(forces) result(vector)
   !> Band forces `(3, nat, n_images)`
   real(wp), intent(in) :: forces(:, :, :)
   !> Packed internal forces (hartree/bohr)
   real(wp), allocatable :: vector(:)

   integer :: k, nat, offset, ncomp

   nat = size(forces, 2)
   ncomp = 3*nat
   allocate(vector(ncomp*(size(forces, 3) - 2)))
   offset = 0
   do k = 2, size(forces, 3) - 1
      vector(offset + 1:offset + ncomp) = reshape(forces(:, :, k), [ncomp])
      offset = offset + ncomp
   end do

end function pack_internal_vector

!> Norm of the band force over the internal images (hartree/bohr).
!>
!> The divergence safety net of the propagation compares this norm with its
!> initial value: the step is a truncated L-BFGS step without a line search,
!> and the quantities a line search could use are not monotone along the
!> propagation (see the module documentation), so the only cheap global guard
!> is the force itself.  The climbing image is included: its force is the
!> climbing force of Eq. 5, which vanishes at the saddle point as well.
pure function band_force_norm(forces) result(norm)
   !> Band forces `(3, nat, n_images)`
   real(wp), intent(in) :: forces(:, :, :)
   !> Norm of the internal image forces (hartree/bohr)
   real(wp) :: norm

   integer :: k

   norm = 0.0_wp
   do k = 2, size(forces, 3) - 1
      norm = norm + sum(forces(:, :, k)**2)
   end do
   norm = sqrt(norm)

end function band_force_norm

!> Clip a band step so that no component exceeds `max_displ` (in place).
pure subroutine clip_band_step(step, max_displ)
   !> Step (bohr), overwritten in place
   real(wp), intent(inout) :: step(:)
   !> Upper bound of one component (bohr)
   real(wp), intent(in) :: max_displ

   real(wp) :: largest

   largest = maxval(abs(step))
   if (largest > max_displ .and. largest > 0.0_wp) then
      step = step*max_displ/largest
   end if

end subroutine clip_band_step

!> Distance between consecutive images of a band (bohr).
subroutine landscape_band_distances(images, distances)
   !> Band images `(3, nat, n_images)`
   real(wp), intent(in) :: images(:, :, :)
   !> Distances (bohr), allocated on return
   real(wp), allocatable, intent(out) :: distances(:)

   integer :: k, n_images

   n_images = size(images, 3)
   allocate(distances(max(n_images - 1, 0)))
   do k = 1, n_images - 1
      distances(k) = sqrt(sum((images(:, :, k + 1) - images(:, :, k))**2))
   end do

end subroutine landscape_band_distances

end module xtb_landscape_path_neb
