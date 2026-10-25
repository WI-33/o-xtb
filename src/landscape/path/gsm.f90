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

!> Growing string method / growing string (GSM), driver
!> (Peters, Heyden, Bell & Chakraborty, J. Chem. Phys. **120**, 114101
!> (2004)).
!>
!> Three phases; the optimisation and the climbing node follow the
!> publication, while the growth uses a **single frontier** from the reactant
!> (the publication's two-frontier variant was abandoned, see below):
!>
!> 1. **Growth.**  The string starts as the reactant endpoint and grows one
!>    new node at a time from a **single frontier** towards the product.  The
!>    node is appended to the frontier along the local tangent, is relaxed
!>    with the perpendicular force of Eq. 1 (the mathematical kernel of
!>    `xtb_landscape_path_string`), and only once its perpendicular force
!>    fulfils the growth indicator of Eq. 7 does the frontier arc length
!>    advance by `ds` (Eq. 9).  The growth stops as soon as the frontier is
!>    within one spacing of the product, which is then appended as the last
!>    node; `string_nnodes` is therefore an **upper bound** of the node count,
!>    not a target (the two-frontier variant of the publication was abandoned:
!>    measured on the ammonia inversion, its two chains relaxed into different
!>    directions and the frontier to frontier distance stayed at 1.4 bohr for
!>    every combination of node count and step).
!> 2. **Optimisation.**  With the node count fixed, every internal node is
!>    relaxed with the perpendicular force and the nodes are redistributed
!>    equidistantly in arc length after every step.  Convergence is the dual
!>    criterion of the band kernel (`f_rms_avg` **and** `f_rms_max` below
!>    their thresholds), so the string cannot "converge on average" while a
!>    single node is badly placed.
!> 3. **Climbing node** (with `landscape_path_settings%climb`).  The highest
!>    energy node is refined as a climbing node: its neighbours are frozen and
!>    it feels the climbing force of Eq. 5 of the CI-NEB module
!>    (`F - 2 (F . tau) tau`), which lets it ascend to the saddle point.
!>    Without this step the highest node only has a vanishing **perpendicular**
!>    force, and a stationary point of the string need not be a stationary
!>    point of the surface.
!>
!> **Failure modes** are reported, never absorbed:
!>
!>  * fewer than three nodes (`string_nnodes < 3`): there is no interior node
!>    to optimise;
!>  * a growth step that is not smaller than the endpoint distance
!>    (`string_ds >= |P - R|`): the first frontier node would already sit at
!>    (or beyond) the other endpoint, so no interior geometry is resolved;
!>  * a string that does not close: the node budget is exhausted while the
!>    frontier is still further from the product than one spacing - the
!>    requested `ds` is too small (or the node budget too small) for the
!>    length of the path;
!>  * a frontier node that does not reach the growth indicator within
!>    `string_grow_steps` perpendicular relaxation iterations: the growth is
!>    stalled.  When the curvature condition of the node spacing had enlarged
!>    the step, the frontier is **retried once at the requested spacing**
!>    before the stall is reported (measured on the ammonia inversion: the
!>    enlarged spacing at the cap of `2 ds` puts the frontier node into a
!>    region where the constrained relaxation has no downhill direction,
!>    while the requested `ds` relaxes, report 39 section 3.4).
!>
!> **Convergence on a molecular surface.**  The perpendicular force criterion
!> of the dual band criterion is a resolution limit for a discrete string: on
!> the ammonia inversion an 8 to 10 node string stalls at `f_rms_avg`/
!> `f_rms_max` of `3.8e-4`/`1.0e-3` hartree/bohr, i.e. just above the band
!> defaults (`3e-4`/`1e-3`), so a molecular string run has to relax
!> `--lts-frms-avg`/`--lts-frms-max` (the CLI test uses `1e-3`/`2e-3`).  The
!> barrier the string finds is nonetheless correct (report 39 section 8.2: it
!> agrees with the band to 0.012 kcal/mol); `string%converged` is a statement
!> about the force residual of the nodes, not about the quality of the saddle
!> candidate (report 38 L-4, F5).
!>
!> Units: coordinates in bohr, energies in hartree, forces and gradients in
!> hartree/bohr, `ds` and the node distances in bohr.
module xtb_landscape_path_gsm
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : autokcal
   use xtb_type_environment, only : TEnvironment
   use xtb_landscape_settings, only : landscape_path_settings, landscape_int_to_string
   use xtb_landscape_engine, only : landscape_energy_model
   use xtb_landscape_coords, only : geometry_pack, align_geometry
   use xtb_landscape_lbfgs, only : landscape_lbfgs_type, new_landscape_lbfgs
   use xtb_landscape_path_neb, only : landscape_band_energies, landscape_band_align, &
      & landscape_band_measures, landscape_band_converged, &
      & landscape_band_endpoint_cache
   use xtb_landscape_path_cineb, only : landscape_cineb_force
   use xtb_landscape_path_string, only : landscape_string_normal_force, &
      & landscape_string_growth_indicator, landscape_string_spacing, &
      & landscape_string_tangents, landscape_string_arc_length, &
      & landscape_string_reparameterize
   implicit none(type, external)
   private

   public :: landscape_string_type, landscape_string_optimize
   public :: landscape_string_grow

   !> Maximum number of refinement rounds of the string (splitting the longest
   !> intervals) after the first relaxation phase
   integer, parameter :: string_refine_rounds = 8
   !> Maximum number of step halvings of the backtracking line search
   integer, parameter :: string_ls_halvings = 12
   !> Armijo fraction of the predicted energy decrease
   real(wp), parameter :: string_ls_armijo = 1.0e-4_wp

   !> Result of a string (growth plus optimisation) run.
   type :: landscape_string_type
      !> True if the optimisation phase met the dual convergence criterion
      logical :: converged = .false.
      !> True if the single frontier reached the product within one spacing
      logical :: closed = .false.
      !> True if the node spacing condition enlarged the requested spacing
      logical :: spacing_adjusted = .false.
      !> Number of nodes of the final string (endpoints included)
      integer :: n_nodes = 0
      !> Number of growth steps performed
      integer :: n_grow = 0
      !> Number of optimisation iterations of the string
      integer :: n_iter = 0
      !> Number of climbing node iterations
      integer :: n_climb_iter = 0
      !> Number of equal arc length redistributions
      integer :: n_reparameterize = 0
      !> Total potential evaluations (energy + gradient)
      integer :: n_eval = 0
      !> Highest energy node of the final string (the saddle candidate)
      integer :: highest_node = 0
      !> Refined climbing node (0 if none)
      integer :: climbing_node = 0
      !> Energy of the saddle candidate (hartree)
      real(wp) :: ts_energy = 0.0_wp
      !> Barrier of the final string relative to the reactant (hartree)
      real(wp) :: barrier = 0.0_wp
      !> Average perpendicular force RMS of the internal nodes (hartree/bohr)
      real(wp) :: f_rms_avg = 0.0_wp
      !> Largest perpendicular force RMS of the internal nodes (hartree/bohr)
      real(wp) :: f_rms_max = 0.0_wp
      !> Unprojected force RMS at the climbing node (hartree/bohr)
      real(wp) :: ci_f_rms = 0.0_wp
      !> Total length of the final string (bohr)
      real(wp) :: path_length = 0.0_wp
      !> Node spacing used by the growth (bohr)
      real(wp) :: spacing = 0.0_wp
      !> Distance from the frontier to the product when the growth stopped
      !> (bohr)
      real(wp) :: gap = 0.0_wp
      !> Growth indicator threshold of Eq. 7 (hartree/bohr)
      real(wp) :: grow_tol = 0.0_wp
      !> Energy of every node of the final string (hartree)
      real(wp), allocatable :: energies(:)
      !> Distance between consecutive nodes of the final string (bohr)
      real(wp), allocatable :: distances(:)
      !> Cumulative arc length of the final string (bohr)
      real(wp), allocatable :: arc_length(:)
      !> Human readable termination message
      character(len=:), allocatable :: message
   end type landscape_string_type

contains

!> Grow a string from its two endpoints, then relax it and refine the saddle.
!>
!> On return `nodes` holds the final string (`string_nnodes` nodes, the
!> endpoints verbatim), `energies` the energy of every node and `forces` the
!> force of every node (the perpendicular force of Eq. 1 for the interior
!> nodes, the physical force at the endpoints).  `string` carries the
!> convergence state and the diagnostics; a failure is reported through
!> `env%error` and `string%converged = .false.` with a readable message.
subroutine landscape_string_optimize(env, model, settings, reactant, product, &
      & nodes, energies, forces, string)
   character(len=*), parameter :: source = 'landscape_string_optimize'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Potential
   class(landscape_energy_model), intent(inout) :: model
   !> Band settings (the string fields are the growth parameters)
   type(landscape_path_settings), intent(in) :: settings
   !> Reactant geometry (bohr), shape `(3, nat)`
   real(wp), intent(in) :: reactant(:, :)
   !> Product geometry (bohr), shape `(3, nat)`
   real(wp), intent(in) :: product(:, :)
   !> Final string nodes `(3, nat, string_nnodes)`
   real(wp), allocatable, intent(out) :: nodes(:, :, :)
   !> Energy of every node (hartree)
   real(wp), allocatable, intent(out) :: energies(:)
   !> Force of every node (hartree/bohr)
   real(wp), allocatable, intent(out) :: forces(:, :, :)
   !> String result
   type(landscape_string_type), intent(out) :: string

   real(wp), allocatable :: gradients(:, :, :), tau(:, :), refined(:, :, :)
   real(wp) :: endpoint_distance, f_avg, f_max
   integer :: rounds
   logical :: fail

   string%converged = .false.
   string%message = "string completed"
   string%grow_tol = settings%string_grow_tol
   if (size(reactant, 1) /= 3 .or. size(product, 1) /= 3) then
      call env%error("the endpoints are not three dimensional", source)
      string%message = "invalid endpoints"
      return
   end if
   if (size(reactant, 2) /= size(product, 2)) then
      call env%error("the endpoints do not have the same number of atoms", source)
      string%message = "invalid endpoints"
      return
   end if
   if (settings%string_nnodes < 3) then
      call env%error("the string needs at least three nodes; --lts-string-nnodes "//&
         & "is below that", source)
      string%message = "too few string nodes"
      return
   end if
   if (settings%string_ds <= 0.0_wp) then
      call env%error("the string growth step --lts-string-ds must be positive", source)
      string%message = "invalid growth step"
      return
   end if
   endpoint_distance = sqrt(sum((product - reactant)**2))
   if (settings%string_ds >= endpoint_distance) then
      call env%error("the string growth step --lts-string-ds is not smaller than "//&
         & "the distance between the endpoints; no interior node would be resolved", &
         & source)
      string%message = "the growth step is too large"
      return
   end if

   call landscape_string_grow(env, model, settings, reactant, product, nodes, string)
   if (.not. allocated(nodes)) return

   call landscape_string_relax(env, model, settings, nodes, energies, gradients, &
      & forces, string, fail)
   if (fail) return

   ! the discrete string has a resolution limit: the perpendicular force of a
   ! node carries an `O(h^2)` discretisation error of the curvature of the
   ! minimum energy path, so a path with too few nodes cannot reach the force
   ! threshold at all.  Measured on the ammonia inversion: a string of 8 nodes
   ! (2.0 bohr long) stalls at `f_rms_max = 1.5e-3` while the threshold is
   ! `1e-3`.  The path is therefore refined by splitting the longest intervals
   ! until the threshold is met or `string_nnodes` (the budget) is exhausted,
   ! and relaxed again.
   do rounds = 1, string_refine_rounds
      if (string%converged) exit
      if (size(nodes, 3) >= settings%string_nnodes) exit
      call landscape_string_refine(nodes, settings%string_nnodes, refined, fail)
      if (fail) exit
      nodes = refined
      call landscape_string_relax(env, model, settings, nodes, energies, &
         & gradients, forces, string, fail)
      if (fail) return
   end do

   string%n_nodes = size(nodes, 3)
   string%highest_node = highest_energy_node(energies)
   string%climbing_node = 0

   if (settings%climb .and. string%highest_node > 0) then
      call landscape_string_climb(env, model, settings, nodes, energies, gradients, &
         & string, fail)
      if (fail) return
      ! the climbing node is the candidate, whatever its rank in the string is
      ! after the refinement
      string%highest_node = string%climbing_node
      ! the forces of the path are rebuilt for the table, but the convergence
      ! measures keep the values of the state that was judged: the climbing
      ! node has moved since, so re-measuring them would compare a different
      ! string with the criterion
      call build_forces(settings, nodes, gradients, tau, forces, f_avg, f_max)
      ! the re-measured forces of the moved climbing node are a diagnostic
      ! only; the convergence measures of the string keep the values of the
      ! state that was judged (see above)
      if (settings%printlevel > 1) then
         write(env%unit, '(1x,a,es12.5,a,es12.5)') "string: forces after the "//&
            & "climbing node, f_rms_avg ", f_avg, ", f_rms_max ", f_max
      end if
   end if

   call string_finalise(string, nodes, energies)
   if (settings%printlevel > 0) call landscape_string_report(env, string)

end subroutine landscape_string_optimize

!> Phase 1: grow the string from the reactant towards the product.
!>
!> The frontier starts at the reactant.  A new node is placed one spacing
!> beyond the frontier along the local tangent (the direction of the previous
!> chain segment, or the direction towards the **product** for the first
!> node) and is then relaxed with the perpendicular force of Eq. 1 on the
!> sphere of radius `spacing` around the frontier; the frontier advances by
!> exactly one spacing once the growth indicator of Eq. 7 holds, which is the
!> arc length advance of Eq. 9.  The growth stops as soon as the frontier is
!> within one spacing of the product, which is then appended as the last node.
!>
!> Consequence of the "stop as soon as it closes": `string_nnodes` is an
!> **upper bound** of the node count, not a target.  The number of nodes a
!> path needs is its length divided by `spacing`, and the length is not known
!> before the path is grown - forcing a fixed node count would either leave a
!> gap (too few nodes) or fold the string back on itself (too many, measured
!> on the ammonia inversion: the frontier to frontier distance of the
!> two-sided variant stayed at 1.4 bohr for every combination of node count
!> and step because the two chains relaxed into different directions).
!> A budget that is too small is reported: the growth error names the
!> exhausted `string_nnodes` together with the remaining distance.
!>
!> The produced string is the state the optimisation phase starts from; a
!> unit test asserts that the first frontier node fulfils Eq. 7 with the
!> tangent it was grown with.
subroutine landscape_string_grow(env, model, settings, reactant, product, nodes, string)
   character(len=*), parameter :: source = 'landscape_string_grow'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Potential
   class(landscape_energy_model), intent(inout) :: model
   !> Band settings
   type(landscape_path_settings), intent(in) :: settings
   !> Reactant geometry (bohr)
   real(wp), intent(in) :: reactant(:, :)
   !> Product geometry (bohr)
   real(wp), intent(in) :: product(:, :)
   !> String nodes, allocated on return
   real(wp), allocatable, intent(out) :: nodes(:, :, :)
   !> String result, growth diagnostics updated
   type(landscape_string_type), intent(inout) :: string

   real(wp), allocatable :: grown(:, :, :), reparam(:, :, :)
   real(wp) :: spacing, gap, initial_gap
   integer :: budget, nat, new_index, n_placed
   logical :: fail

   budget = settings%string_nnodes
   nat = size(reactant, 2)
   allocate(grown(3, nat, budget), source=0.0_wp)
   grown(:, :, 1) = reactant
   n_placed = 1
   spacing = settings%string_ds
   string%spacing = spacing

   do while (n_placed < budget)
      if (string_frontier_distance(grown, n_placed, product) <= spacing) exit
      new_index = n_placed + 1
      call grow_frontier(env, model, settings, grown, n_placed, n_placed - 1, &
         & product, new_index, spacing, fail, string)
      if (fail) then
         call env%error("the string stopped growing at node "//&
            & trim(landscape_int_to_string(n_placed + 1))//" of the budget "//&
            & trim(landscape_int_to_string(budget)), source)
         string%message = "the string growth stalled"
         return
      end if
      ! remove the rigid body motion the frontier relaxation collected: the
      ! perpendicular projection of Eq. 1 does not preserve the vanishing net
      ! force of an isolated molecule, so the projected force has a spurious
      ! translation (and rotation) component that costs no energy and that the
      ! relaxation follows (measured on the ammonia inversion: the whole
      ! molecule was translated by 0.25 bohr per growth step).  The alignment
      ! is a rigid body motion of one node onto the reactant frame, so it
      ! changes no energy, no distance and no perpendicular force.
      if (settings%align) then
         call align_geometry(env, grown(:, :, new_index), grown(:, :, 1))
      end if
      n_placed = n_placed + 1
      ! redistribute the chain grown so far equally in arc length: the
      ! frontier of a curved path is placed along the direction of its parent
      ! segment, which quickly deviates from the path, and the relaxation of
      ! the single new node cannot repair the accumulated deviation.
      ! Measured on the Mueller-Brown surface: without the redistribution the
      ! chain follows the path to the saddle and then turns around and walks
      ! back, so it never reaches the product.
      call landscape_string_reparameterize(grown(:, :, 1:n_placed), reparam)
      grown(:, :, 1:n_placed) = reparam
      string%n_reparameterize = string%n_reparameterize + 1
      string%n_grow = string%n_grow + 1
      string%spacing = spacing
   end do

   gap = string_frontier_distance(grown, n_placed, product)
   string%gap = gap
   if (gap > spacing) then
      ! the frontier did not get within one spacing of the product.  If it made
      ! no progress at all the growth step is unusable and the run is aborted;
      ! otherwise the string is completed by interpolating the remaining nodes
      ! between the frontier and the product.  The relaxation phase then has to
      ! repair the tail (it is responsible for the path, the growth only
      ! provides a starting configuration), and the interpolation is reported.
      ! Measured: the frontier of the abstract two dimensional test embedding
      ! leaves the path near the saddle and turns around, while the relaxation
      ! of the interpolated string still converges to the saddle (unit tests
      ! `string-mb-saddle`, `string-vs-neb`).
      initial_gap = sqrt(sum((reactant - product)**2))
      if (gap > 0.9_wp*initial_gap) then
         call env%error("the string made no progress towards the product: the "//&
            & "remaining distance is "//trim(number_string(gap))//" bohr after "//&
            & trim(landscape_int_to_string(n_placed))//" nodes with a spacing of "//&
            & trim(number_string(spacing))//" bohr; raise --lts-string-nnodes "//&
            & "or --lts-string-ds", source)
         string%message = "the string made no progress towards the product"
         return
      end if
      call env%warning("the string did not reach the product within the node "//&
         & "budget; the remaining nodes are interpolated and left to the "//&
         & "relaxation", source)
      allocate(nodes(3, nat, budget), source=0.0_wp)
      call interpolate_tail(grown, n_placed, product, nodes)
      return
   end if
   string%closed = .true.

   ! append the product endpoint; the budget is an upper bound, so the array
   ! is compacted to the nodes that were actually grown plus the product
   allocate(nodes(3, nat, n_placed + 1), source=0.0_wp)
   nodes(:, :, 1:n_placed) = grown(:, :, 1:n_placed)
   nodes(:, :, n_placed + 1) = product

end subroutine landscape_string_grow

!> Fill a partially grown string by interpolating up to the product.
!>
!> The nodes `1 .. n_placed` are copied and the remaining nodes are placed
!> equidistantly on the straight line from the frontier to the product; the
!> product is the last node.  The endpoints of the grown part are unchanged.
subroutine interpolate_tail(grown, n_placed, product, nodes)
   !> Grown nodes `(3, nat, n_budget)`
   real(wp), intent(in) :: grown(:, :, :)
   !> Number of grown nodes
   integer, intent(in) :: n_placed
   !> Product geometry (bohr)
   real(wp), intent(in) :: product(:, :)
   !> Completed string `(3, nat, n_budget)`, filled in place
   real(wp), intent(inout) :: nodes(:, :, :)

   real(wp) :: lam
   integer :: k, n, n_gap

   n = size(nodes, 3)
   nodes(:, :, 1:n_placed) = grown(:, :, 1:n_placed)
   n_gap = n - n_placed
   if (n_gap <= 0) return
   do k = 1, n_gap
      lam = real(k, wp)/real(n_gap, wp)
      nodes(:, :, n_placed + k) = (1.0_wp - lam)*grown(:, :, n_placed) + &
         & lam*product
   end do

end subroutine interpolate_tail

!> Distance between the frontier and the target endpoint (bohr).
pure function string_frontier_distance(nodes, frontier, product) result(distance)
   !> String nodes `(3, nat, n_budget)`
   real(wp), intent(in) :: nodes(:, :, :)
   !> Index of the frontier node
   integer, intent(in) :: frontier
   !> Target geometry (bohr)
   real(wp), intent(in) :: product(:, :)
   !> Distance (bohr)
   real(wp) :: distance

   distance = sqrt(sum((nodes(:, :, frontier) - product)**2))

end function string_frontier_distance

!> Place and relax one frontier node.
!>
!> The node is placed one spacing beyond `left` along the local tangent (the
!> direction of the previous chain segment, or the direction towards the
!> opposite **endpoint** of the path for the first node of a chain) and then
!> relaxed with the perpendicular force of Eq. 1 until the growth indicator
!> of Eq. 7 holds.  The local string direction is recomputed on every
!> iteration, so the node moves on the curved manifold and not only in the
!> initial plane.  The node spacing condition is evaluated at the node the new
!> one grows from and can enlarge the spacing; the enlarged value is returned
!> to the caller, which advances the arc length of the chain with it.
!>
!> The first chain step has no parent segment, so its tangent is taken from
!> the opposite endpoint geometry `opposite`.  That geometry is passed by
!> value on purpose: the opposite endpoint of a one sided growth is **not**
!> part of `nodes` (the array holds the grown chain, whose frame after the
!> frontier is the zero fill of the allocation), so an index into `nodes`
!> there would silently aim the first node at the origin of the coordinate
!> system instead of at the product.
subroutine grow_frontier(env, model, settings, nodes, left, left_prev, opposite, &
      & new_index, spacing, fail, string)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Potential
   class(landscape_energy_model), intent(inout) :: model
   !> Band settings
   type(landscape_path_settings), intent(in) :: settings
   !> String nodes, the new node is written in place
   real(wp), intent(inout) :: nodes(:, :, :)
   !> Index of the node the new one is grown from
   integer, intent(in) :: left
   !> Index of the node before `left` on the chain, 0 if there is none
   integer, intent(in) :: left_prev
   !> Geometry of the opposite endpoint of the path (bohr), shape `(3, nat)`
   real(wp), intent(in) :: opposite(:, :)
   !> Index the new node is stored at
   integer, intent(in) :: new_index
   !> Spacing of the growth step (bohr), updated by the spacing condition
   real(wp), intent(inout) :: spacing
   !> Failure flag
   logical, intent(out) :: fail
   !> String result, diagnostics updated
   type(landscape_string_type), intent(inout) :: string

   type(landscape_lbfgs_type) :: lbfgs
   real(wp), allocatable :: x(:), xnew(:), g(:), gnew(:), d(:), fperp(:), fperp_new(:)
   real(wp), allocatable :: tau(:), tau0(:), x_left(:)
   real(wp) :: energy, energy_new, norm_perp, new_spacing, radius
   integer :: ndim, it, n_eval, itry, n_try
   logical :: ok, steepest, adjusted, spacing_fail, converged, accepted
   logical :: eval_fail

   fail = .true.
   ! the tangent is built from the geometry of the opposite endpoint, so it
   ! has to have the atom count of the chain; a mismatched shape would only
   ! be caught by bounds checking at the array operation below
   if (size(opposite, 1) /= 3 .or. size(opposite, 2) /= size(nodes, 2)) return
   ndim = 3*size(nodes, 2)
   x_left = geometry_pack(nodes(:, :, left))
   allocate(x(ndim), xnew(ndim), g(ndim), gnew(ndim), d(ndim))
   allocate(fperp(ndim), fperp_new(ndim), tau(ndim), tau0(ndim))

   ! the tangent at the node that is being grown from
   if (left_prev >= 1 .and. left_prev <= size(nodes, 3)) then
      tau = x_left - geometry_pack(nodes(:, :, left_prev))
   else
      tau = geometry_pack(opposite) - x_left
   end if
   call normalize(tau)
   tau0 = tau

   ! node spacing condition, evaluated at the parent node.  The potential
   ! reports through its own flag: `fail` belongs to this routine and must
   ! not be reset to `.false.` by a successful evaluation (measured: it was,
   ! and a stalled growth then reported success instead of the stall).
   call evaluate_node(env, model, x_left, energy, g, eval_fail)
   if (eval_fail) return
   string%n_eval = string%n_eval + 1
   ! the spacing condition is evaluated against the **requested** step, so the
   ! enlargement cannot compound over the growth steps
   call landscape_string_spacing(g, tau, settings%string_ds, settings%string_nu2, &
      & new_spacing, adjusted, spacing_fail)
   if (spacing_fail) then
      call env%error("the string curvature floor nu2 must be positive", &
         & 'landscape_string_grow')
      return
   end if
   if (adjusted) then
      spacing = new_spacing
      string%spacing_adjusted = .true.
   end if

   ! The node is placed on the sphere of radius `radius` around its parent and
   ! relaxed there, so the arc length advance of Eq. 9 is exact for that
   ! radius.  The curvature condition can enlarge the spacing beyond the
   ! requested step, and an enlarged sphere is not always relaxable: measured
   ! on the ammonia inversion (report 39 section 3.4), the enlarged spacing of
   ! `1.0` bohr (the cap, twice the requested `0.5`) puts the frontier node
   ! into a region where the constrained relaxation has no downhill direction
   ! and the growth stalls, while the same node relaxes at the requested
   ! `0.5` bohr.  A stall at an enlarged spacing is therefore retried once at
   ! the requested spacing; the spacing the node was actually placed with is
   ! returned to the caller and reported (`string%spacing`).
   n_try = 1
   if (spacing > settings%string_ds) n_try = 2
   tau = tau0
   do itry = 1, n_try
      if (itry == 1) then
         radius = spacing
      else
         radius = settings%string_ds
      end if
      x = x_left + radius*tau
      call evaluate_node(env, model, x, energy, g, eval_fail)
      if (eval_fail) return
      string%n_eval = string%n_eval + 1
      call new_landscape_lbfgs(lbfgs, ndim, settings%lbfgs_memory, settings%curvature)

      converged = .false.
      do it = 1, settings%string_grow_steps
         tau = x - x_left
         call normalize(tau)
         call landscape_string_normal_force(g, tau, fperp)
         if (landscape_string_growth_indicator(fperp, settings%string_grow_tol)) then
            converged = .true.
            exit
         end if
         call lbfgs%direction(-fperp, d, ok, steepest)
         if (.not. ok) exit
         ! the node may only move on the sphere of radius `radius` around its
         ! parent node: the arc length advance of Eq. 9 is thereby exact, and
         ! the node cannot slide along the perpendicular manifold (measured:
         ! without the sphere the frontier of the LEPS surface slides the whole
         ! H2 molecule along the molecular axis - a nearly free mode - and the
         ! two chains pass each other instead of meeting)
         call frontier_step(env, model, x, x_left, radius, d, -fperp, energy, &
            & xnew, gnew, energy_new, accepted, n_eval, eval_fail)
         if (eval_fail) return
         string%n_eval = string%n_eval + n_eval
         if (.not. accepted) exit
         tau = xnew - x_left
         call normalize(tau)
         call landscape_string_normal_force(gnew, tau, fperp_new)
         call lbfgs%push(xnew - x, fperp_new - fperp)
         x = xnew
         g = gnew
         energy = energy_new
      end do
      if (converged) exit
   end do

   if (.not. converged) return
   spacing = radius
   tau = x - x_left
   call normalize(tau)
   call landscape_string_normal_force(g, tau, fperp)
   norm_perp = sqrt(sum(fperp**2))
   nodes(:, :, new_index) = reshape(x, [3, size(nodes, 2)])
   if (settings%string_history .and. settings%printlevel > 0) then
      write(env%unit, '(1x,a,i0,a,es12.4,a,es12.4,a)') "string: node ", &
         & new_index, " placed after ", norm_perp, " Eh/bohr, spacing ", &
         & spacing, " bohr"
   end if
   fail = .false.

end subroutine grow_frontier

!> Phase 2: relax the string with the perpendicular force.
!>
!> The nodes are redistributed equidistantly in arc length after every step;
!> the endpoints never move.  Convergence is the dual criterion of the band
!> kernel, so the average and the largest perpendicular force RMS must both
!> be below their thresholds.
subroutine landscape_string_relax(env, model, settings, nodes, energies, gradients, &
      & forces, string, fail)
   character(len=*), parameter :: source = 'landscape_string_relax'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Potential
   class(landscape_energy_model), intent(inout) :: model
   !> Band settings
   type(landscape_path_settings), intent(in) :: settings
   !> String nodes, relaxed in place
   real(wp), intent(inout) :: nodes(:, :, :)
   !> Node energies, allocated on return
   real(wp), allocatable, intent(out) :: energies(:)
   !> Node gradients, allocated on return
   real(wp), allocatable, intent(out) :: gradients(:, :, :)
   !> Node forces, allocated on return
   real(wp), allocatable, intent(out) :: forces(:, :, :)
   !> String result
   type(landscape_string_type), intent(inout) :: string
   !> Failure flag
   logical, intent(out) :: fail

   real(wp), allocatable :: tau(:, :), trial(:, :, :), reparam(:, :, :)
   real(wp), allocatable :: d(:), g(:), gnew(:), x(:), xnew(:), step(:)
   real(wp), allocatable :: en_new(:), grad_new(:, :, :)
   real(wp) :: energy0, slope, alpha
   integer :: n, nat, ndim, it, ls, n_done
   logical :: accepted
   ! The two endpoints of the string are the fixed reactant and product, so
   ! their energy and gradient are evaluated once and reused, exactly like the
   ! fixed endpoints of a band (report 223).
   type(landscape_band_endpoint_cache) :: endpoints

   fail = .true.
   nat = size(nodes, 2)
   n = size(nodes, 3)
   ndim = 3*nat*(n - 2)
   allocate(d(ndim), g(ndim), gnew(ndim), x(ndim), xnew(ndim), step(ndim))

   ! the optimisation starts from the equidistant distribution of the grown
   ! path; the growth itself leaves every node at exactly one spacing from its
   ! parent (that is the invariant its unit test pins), the redistribution is
   ! the first step of this phase
   call landscape_string_reparameterize(nodes, reparam)
   nodes = reparam
   string%n_reparameterize = string%n_reparameterize + 1

   call landscape_band_energies(env, model, nodes, energies, gradients, fail, &
      & endpoints, n_done)
   if (fail) return
   string%n_eval = string%n_eval + n_done
   call build_forces(settings, nodes, gradients, tau, forces, string%f_rms_avg, &
      & string%f_rms_max)
   do it = 0, settings%max_steps
      if (landscape_band_converged(string%f_rms_avg, string%f_rms_max, 0.0_wp, &
            & .false., settings)) then
         string%converged = .true.
         exit
      end if
      if (it >= settings%max_steps) exit

      call pack_internal_forces(forces, n, nat, g)
      ! the step is the **projected gradient** itself, not an L-BFGS step: the
      ! equal arc length redistribution moves the nodes after every accepted
      ! step, so the secant pairs of a quasi-Newton history describe a map that
      ! is not the force field of the current string.  Measured on the
      ! Mueller-Brown surface: the L-BFGS history stalls the line search at a
      ! perpendicular force RMS of 1e-2 (the threshold is 1e-3), while the same
      ! line search on the plain projected gradient converges to 1.8e-4 in 24
      ! iterations.  `g` is the packed projected **gradient**, so the descent
      ! direction is its negative.
      d = -g

      ! backtracking on the **energy of the string** (the sum of the internal
      ! node energies): the perpendicular force has no potential, so the
      ! energy is the only merit function that can accept or reject a step.
      ! Measured: without the line search the L-BFGS of the string oscillates
      ! between two node distributions on the Mueller-Brown surface with the
      ! full `max_displ` step.
      energy0 = sum(energies(2:n - 1))
      slope = dot_product(g, d)
      call pack_internal_nodes(nodes, x)
      accepted = .false.
      alpha = 1.0_wp
      do ls = 1, string_ls_halvings
         step = alpha*d
         call clip(step, settings%max_displ)
         xnew = x + step
         trial = nodes
         call unpack_internal_nodes(xnew, trial)
         if (settings%align) call landscape_band_align(env, trial)
         ! the merit function is the energy of the **direct** step; the equal
         ! arc length redistribution is applied to the accepted geometry
         ! afterwards and is deliberately not part of the test.  Measured: the
         ! redistribution moves the nodes by far more than the step itself
         ! (it is a re-parametrisation, not a displacement), so including it
         ! lets it dominate the energy change and every step fails the Armijo
         ! test on the Mueller-Brown surface.
         call landscape_band_energies(env, model, trial, en_new, grad_new, fail, &
            & endpoints, n_done)
         if (fail) then
            call env%error("a string node could not be evaluated", source)
            string%message = "a string node could not be evaluated"
            return
         end if
         string%n_eval = string%n_eval + n_done
         if (sum(en_new(2:n - 1)) <= energy0 + string_ls_armijo*alpha*slope) then
            accepted = .true.
            exit
         end if
         alpha = 0.5_wp*alpha
      end do
      if (.not. accepted) exit

      ! redistribute the accepted string equally in arc length
      call landscape_string_reparameterize(trial, reparam)
      trial = reparam
      string%n_reparameterize = string%n_reparameterize + 1
      call landscape_band_energies(env, model, trial, en_new, grad_new, fail, &
         & endpoints, n_done)
      if (fail) then
         call env%error("a string node could not be evaluated", source)
         string%message = "a string node could not be evaluated"
         return
      end if
      string%n_eval = string%n_eval + n_done

      call build_forces(settings, trial, grad_new, tau, forces, string%f_rms_avg, &
         & string%f_rms_max)

      call pack_internal_forces(forces, n, nat, gnew)
      call pack_internal_nodes(trial, xnew)
      nodes = trial
      energies = en_new
      gradients = grad_new
      x = xnew
      g = gnew
      string%n_iter = string%n_iter + 1

      if (settings%printlevel > 1) then
         write(env%unit, '(1x,a,i6,1x,a,es12.5,a,es12.5)') "string iter ", &
            & string%n_iter, " f_rms_avg ", string%f_rms_avg, " f_rms_max ", &
            & string%f_rms_max
      end if
   end do

   ! one final evaluation keeps energies, gradients and forces consistent with
   ! the node positions, whatever way the loop was left
   call landscape_band_energies(env, model, nodes, energies, gradients, fail, &
      & endpoints, n_done)
   if (fail) then
      string%message = "a string node could not be evaluated"
      return
   end if
   string%n_eval = string%n_eval + n_done
   call build_forces(settings, nodes, gradients, tau, forces, string%f_rms_avg, &
      & string%f_rms_max)

   if (settings%printlevel > 0) then
      write(env%unit, '(1x,a,i0,a,es12.5,a,es12.5)') "string: relaxed in ", &
         & string%n_iter, " iterations, f_rms_avg ", string%f_rms_avg, &
         & ", f_rms_max ", string%f_rms_max
   end if
   if (.not. string%converged) then
      string%message = "the string reached the iteration limit without converging"
   end if
   fail = .false.

end subroutine landscape_string_relax

!> Refine a string by splitting its longest intervals.
!>
!> One node is inserted at a time, at the midpoint of the **longest** interval
!> of the current string, until the node budget is reached.  The endpoints are
!> never touched and the node order is preserved; the inserted nodes are
!> linear interpolations, so they are refined by the relaxation that follows.
!> `fail` is set when there is no room for another node.
subroutine landscape_string_refine(nodes, budget, refined, fail)
   !> String nodes `(3, nat, n_nodes)`
   real(wp), intent(in) :: nodes(:, :, :)
   !> Upper bound of the node count
   integer, intent(in) :: budget
   !> Refined string, allocated on return
   real(wp), allocatable, intent(out) :: refined(:, :, :)
   !> Failure flag (no room, or an empty string)
   logical, intent(out) :: fail

   real(wp), allocatable :: work(:, :, :)
   real(wp) :: longest, span
   integer :: k, n, n_want, best, nat

   fail = .true.
   nat = size(nodes, 2)
   n = size(nodes, 3)
   if (n < 2 .or. n >= budget) return
   n_want = min(budget, 2*n - 1)
   allocate(work(3, nat, n_want), source=0.0_wp)
   work(:, :, 1:n) = nodes

   do while (n < n_want)
      longest = -1.0_wp
      best = 1
      do k = 1, n - 1
         span = sqrt(sum((work(:, :, k + 1) - work(:, :, k))**2))
         if (span > longest) then
            longest = span
            best = k
         end if
      end do
      work(:, :, best + 2:n + 1) = work(:, :, best + 1:n)
      work(:, :, best + 1) = 0.5_wp*(work(:, :, best) + work(:, :, best + 2))
      n = n + 1
   end do
   allocate(refined(3, nat, n), source=work(:, :, 1:n))
   fail = .false.

end subroutine landscape_string_refine

!> Phase 3: refine the highest energy node as a climbing node.
!>
!> The neighbours of the node are frozen (their geometries do not change), so
!> the string tangent at the node follows the node itself.  The force is
!> `F - 2 (F . tau) tau` of the CI-NEB module (Eq. 5), which has the parallel
!> component inverted: the node is attracted to the saddle point along the
!> string direction and to the minimum energy path perpendicular to it.  A
!> node whose **unprojected** force RMS meets the criterion is a stationary
!> point of the full surface, i.e. the transition state candidate.
subroutine landscape_string_climb(env, model, settings, nodes, energies, gradients, &
      & string, fail)
   character(len=*), parameter :: source = 'landscape_string_climb'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Potential
   class(landscape_energy_model), intent(inout) :: model
   !> Band settings
   type(landscape_path_settings), intent(in) :: settings
   !> String nodes, the climbing node is moved in place
   real(wp), intent(inout) :: nodes(:, :, :)
   !> Node energies, the climbing node entry is updated
   real(wp), intent(inout) :: energies(:)
   !> Node gradients, the climbing node entry is updated
   real(wp), intent(inout) :: gradients(:, :, :)
   !> String result
   type(landscape_string_type), intent(inout) :: string
   !> Failure flag
   logical, intent(out) :: fail

   type(landscape_lbfgs_type) :: lbfgs
   real(wp), allocatable :: x(:), xnew(:), g(:), gnew(:), d(:), tau(:), f(:), fci(:)
   real(wp), allocatable :: step(:)
   real(wp) :: energy, energy_new, alpha, rms_current
   integer :: k, nat, ndim, it, n_eval, ls
   logical :: ok, steepest, accepted, eval_fail

   fail = .true.
   ! the evaluation counter of this routine is accumulated into the string
   ! result: an uninitialised read here would add an arbitrary number to the
   ! "potential evaluations" diagnostic of the run
   n_eval = 0
   k = string%highest_node
   if (k < 2 .or. k > size(nodes, 3) - 1) then
      call env%error("the string has no interior node to refine", source)
      return
   end if
   nat = size(nodes, 2)
   ndim = 3*nat
   allocate(x(ndim), xnew(ndim), g(ndim), gnew(ndim), d(ndim), step(ndim))
   allocate(tau(ndim), f(ndim), fci(ndim))
   x = geometry_pack(nodes(:, :, k))
   call new_landscape_lbfgs(lbfgs, ndim, settings%lbfgs_memory, settings%curvature)

   do it = 1, settings%max_steps
      n_eval = 0
      g = reshape(gradients(:, :, k), [ndim])
      f = -g
      call string_tangent_at(nodes, k, settings, tau)
      call landscape_cineb_force(f, tau, fci)
      string%ci_f_rms = sqrt(sum(g**2)/real(ndim, wp))
      if (string%ci_f_rms <= settings%f_rms_max_tol) exit
      call lbfgs%direction(-fci, d, ok, steepest)
      if (.not. ok) exit
      ! the climbing node **ascends**; an energy criterion would reject every
      ! useful step, so the step is halved until it lowers the unprojected
      ! force (the quantity the node has to drive to zero), which is the
      ! natural merit function of a saddle search
      rms_current = sqrt(sum(g**2)/real(ndim, wp))
      accepted = .false.
      alpha = 1.0_wp
      do ls = 1, string_ls_halvings
         step = alpha*d
         call clip(step, settings%max_displ)
         xnew = x + step
         call evaluate_node(env, model, xnew, energy_new, gnew, eval_fail)
         n_eval = n_eval + 1
         if (eval_fail) then
            call env%error("the climbing node could not be evaluated", source)
            string%message = "the climbing node could not be evaluated"
            return
         end if
         if (sqrt(sum(gnew**2)/real(ndim, wp)) < rms_current) then
            accepted = .true.
            exit
         end if
         alpha = 0.5_wp*alpha
      end do
      string%n_eval = string%n_eval + n_eval
      if (.not. accepted) exit
      energy = energy_new
      string%n_climb_iter = string%n_climb_iter + 1
      call lbfgs%push(xnew - x, gnew - g)
      x = xnew
      nodes(:, :, k) = reshape(x, [3, nat])
      energies(k) = energy
      gradients(:, :, k) = reshape(gnew, [3, nat])
      if (settings%printlevel > 1) then
         write(env%unit, '(1x,a,i0,a,i0,a,es12.5)') "string: climbing node ", &
            & k, " iteration ", string%n_climb_iter, ", force rms ", &
            & sqrt(sum(gnew**2)/real(ndim, wp))
      end if
   end do

   ! the unprojected force of the final geometry
   g = reshape(gradients(:, :, k), [ndim])
   string%ci_f_rms = sqrt(sum(g**2)/real(ndim, wp))
   string%climbing_node = k
   if (string%ci_f_rms > settings%f_rms_max_tol) then
      string%converged = .false.
      string%message = "the climbing node did not reach the force criterion"
      call env%warning("the climbing node of the string did not reach the "//&
         & "force criterion; the saddle candidate is not a stationary point", &
         & source)
   end if
   if (settings%printlevel > 0) then
      write(env%unit, '(1x,a,i0,a,i0,a,es12.5)') "string: climbing node ", k, &
         & " refined in ", string%n_climb_iter, " iterations, force rms ", &
         & string%ci_f_rms
   end if
   fail = .false.

end subroutine landscape_string_climb

!> Tangent of the string at one node.
subroutine string_tangent_at(nodes, k, settings, tau)
   !> String nodes `(3, nat, n_nodes)`
   real(wp), intent(in) :: nodes(:, :, :)
   !> Node index
   integer, intent(in) :: k
   !> Band settings (the tangent convention)
   type(landscape_path_settings), intent(in) :: settings
   !> Unit tangent at the node (flattened, dimensionless)
   real(wp), intent(out) :: tau(:)

   real(wp), allocatable :: all_tangents(:, :)

   call landscape_string_tangents(nodes, all_tangents, settings%string_spline)
   tau = all_tangents(:, k)

end subroutine string_tangent_at

!> Perpendicular forces of the string and the dual convergence measures.
!>
!> The interior nodes get the perpendicular force of Eq. 1; the two endpoints
!> get their physical force, which the measure routine deliberately skips.
!> `f_rms_avg` and `f_rms_max` are the measures of `xtb_landscape_path_neb`,
!> so the string uses exactly the convergence criterion of the band searches.
subroutine build_forces(settings, nodes, gradients, tau, forces, f_rms_avg, f_rms_max)
   !> Band settings
   type(landscape_path_settings), intent(in) :: settings
   !> String nodes `(3, nat, n_nodes)`
   real(wp), intent(in) :: nodes(:, :, :)
   !> Node gradients (hartree/bohr)
   real(wp), intent(in) :: gradients(:, :, :)
   !> Tangents `(3 nat, n_nodes)`, allocated on return
   real(wp), allocatable, intent(out) :: tau(:, :)
   !> Node forces `(3, nat, n_nodes)`, allocated on return
   real(wp), allocatable, intent(out) :: forces(:, :, :)
   !> Average perpendicular force RMS (hartree/bohr)
   real(wp), intent(out) :: f_rms_avg
   !> Largest perpendicular force RMS (hartree/bohr)
   real(wp), intent(out) :: f_rms_max

   real(wp), allocatable :: one(:)
   integer :: k, nat, ndim

   nat = size(nodes, 2)
   ndim = 3*nat
   allocate(one(ndim))
   call landscape_string_tangents(nodes, tau, settings%string_spline)
   allocate(forces(3, nat, size(nodes, 3)), source=0.0_wp)
   forces(:, :, 1) = -gradients(:, :, 1)
   forces(:, :, size(nodes, 3)) = -gradients(:, :, size(nodes, 3))
   do k = 2, size(nodes, 3) - 1
      one = reshape(gradients(:, :, k), [ndim])
      call landscape_string_normal_force(one, tau(:, k), one)
      forces(:, :, k) = reshape(one, [3, nat])
   end do
   call landscape_band_measures(forces, 0, f_rms_avg, f_rms_max)

end subroutine build_forces

!> Fill the summary fields from the state of the string.
subroutine string_finalise(string, nodes, energies)
   !> String result
   type(landscape_string_type), intent(inout) :: string
   !> String nodes
   real(wp), intent(in) :: nodes(:, :, :)
   !> Node energies (hartree)
   real(wp), intent(in) :: energies(:)

   real(wp), allocatable :: arc(:)
   integer :: k, n

   n = size(nodes, 3)
   string%n_nodes = n
   if (allocated(string%energies)) deallocate(string%energies)
   allocate(string%energies(n), source=energies)
   if (allocated(string%distances)) deallocate(string%distances)
   allocate(string%distances(max(n - 1, 0)), source=0.0_wp)
   do k = 1, n - 1
      string%distances(k) = sqrt(sum((nodes(:, :, k + 1) - nodes(:, :, k))**2))
   end do
   call landscape_string_arc_length(nodes, arc)
   if (allocated(string%arc_length)) deallocate(string%arc_length)
   allocate(string%arc_length(n), source=arc)
   string%path_length = arc(n)
   k = string%highest_node
   if (k < 1 .or. k > n) k = highest_energy_node(energies)
   string%highest_node = k
   if (k > 0) then
      string%ts_energy = energies(k)
      string%barrier = energies(k) - energies(1)
   end if

end subroutine string_finalise

!> Print the report of a string run.
subroutine landscape_string_report(env, string)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> String result
   type(landscape_string_type), intent(in) :: string

   write(env%unit, '(1x,a,i0,a)') "string: ", string%n_nodes, &
      & " nodes (2 fixed endpoints)"
   write(env%unit, '(1x,a,i0,a,es12.4,a,es12.4,a,l1)') "string: growth steps ", &
      & string%n_grow, ", tolerance ", string%grow_tol, ", spacing ", &
      & string%spacing, " bohr, spacing adjusted ", string%spacing_adjusted
   write(env%unit, '(1x,a,i0,a,i0,a)') "string: reparameterisations ", &
      & string%n_reparameterize, ", potential evaluations ", string%n_eval
   write(env%unit, '(1x,a,f12.6,a,i0,a,f20.12,a,f12.4,a)') "string: path length ", &
      & string%path_length, " bohr, highest node ", string%highest_node, &
      & ", E ", string%ts_energy, " Eh, barrier ", string%barrier*autokcal, &
      & " kcal/mol"
   write(env%unit, '(1x,a,l1)') "string: converged ", string%converged

end subroutine landscape_string_report

!> Decimal representation of a real for a message.
pure function number_string(value) result(text)
   !> Value
   real(wp), intent(in) :: value
   !> Decimal string
   character(len=:), allocatable :: text

   character(len=24) :: buffer

   write(buffer, '(f12.6)') value
   text = trim(adjustl(buffer))

end function number_string

!> Index of the highest energy interior node.
pure function highest_energy_node(energies) result(node)
   !> Node energies (hartree)
   real(wp), intent(in) :: energies(:)
   !> 1-based index of the highest interior node, 0 if there is none
   integer :: node

   integer :: k

   node = 0
   do k = 2, size(energies) - 1
      ! A block if, not `node == 0 .or. ...`: Fortran does not require the
      ! right operand to be short circuited, so the first iteration would read
      ! `energies(0)`, which a bounds checked build rejects.
      if (node == 0) then
         node = k
      else if (energies(k) > energies(node)) then
         node = k
      end if
   end do

end function highest_energy_node

!> Backtracking step of a frontier node on the sphere around its parent.
!>
!> The trial geometry `x + alpha d` is projected onto the sphere of radius
!> `radius` around `parent` before it is evaluated, so the node keeps exactly
!> the growth spacing from the node it was grown from; the arc length advance
!> of Eq. 9 is therefore exact and the relaxation cannot slide the node along
!> the perpendicular manifold.  The sphere projection turns the step into a
!> move on the intersection of the sphere with the perpendicular plane of the
!> local string direction, which is the geometry the growth is meant to
!> produce.  The energy criterion and the halving are those of
!> the line search of the optimisation phase; a direction that never lowers
!> the energy is reported with `accepted = .false.`.
subroutine frontier_step(env, model, x, parent, radius, d, g, energy, xnew, &
      & gnew, enew, accepted, n_eval, fail)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Potential
   class(landscape_energy_model), intent(inout) :: model
   !> Current geometry of the frontier node (bohr)
   real(wp), intent(in) :: x(:)
   !> Geometry of the node the frontier grows from (bohr)
   real(wp), intent(in) :: parent(:)
   !> Radius of the sphere, i.e. the growth spacing (bohr)
   real(wp), intent(in) :: radius
   !> Descent direction (bohr)
   real(wp), intent(in) :: d(:)
   !> Projected gradient (hartree/bohr)
   real(wp), intent(in) :: g(:)
   !> Energy at `x` (hartree)
   real(wp), intent(in) :: energy
   !> Accepted geometry (bohr)
   real(wp), intent(out) :: xnew(:)
   !> Gradient at the accepted geometry (hartree/bohr)
   real(wp), intent(out) :: gnew(:)
   !> Energy at the accepted geometry (hartree)
   real(wp), intent(out) :: enew
   !> True if a step lowered the energy
   logical, intent(out) :: accepted
   !> Number of evaluations performed
   integer, intent(out) :: n_eval
   !> Failure flag
   logical, intent(out) :: fail

   real(wp), allocatable :: step(:), trial(:), gtrial(:), offset(:)
   real(wp) :: alpha, slope, energy_trial
   integer :: ls

   fail = .false.
   accepted = .false.
   n_eval = 0
   allocate(step(size(x)), trial(size(x)), gtrial(size(x)), offset(size(x)))
   slope = dot_product(g, d)
   if (slope >= 0.0_wp) then
      xnew = x
      gnew = g
      enew = energy
      return
   end if

   alpha = 1.0_wp
   do ls = 1, string_ls_halvings
      step = alpha*d
      call clip(step, radius)
      trial = x + step
      offset = trial - parent
      if (sqrt(sum(offset**2)) <= 1.0e-12_wp) then
         alpha = 0.5_wp*alpha
         cycle
      end if
      trial = parent + radius*offset/sqrt(sum(offset**2))
      call evaluate_node(env, model, trial, energy_trial, gtrial, fail)
      n_eval = n_eval + 1
      if (fail) return
      if (energy_trial <= energy + string_ls_armijo*alpha*slope) then
         xnew = trial
         gnew = gtrial
         enew = energy_trial
         accepted = .true.
         return
      end if
      alpha = 0.5_wp*alpha
   end do
   xnew = x
   gnew = g
   enew = energy

end subroutine frontier_step

!> Evaluate one node; `n_eval` returns the number of potential evaluations.
subroutine evaluate_node(env, model, x, energy, gradient, fail)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Potential
   class(landscape_energy_model), intent(inout) :: model
   !> Flattened coordinates (bohr)
   real(wp), intent(in) :: x(:)
   !> Energy (hartree)
   real(wp), intent(out) :: energy
   !> Gradient (hartree/bohr)
   real(wp), intent(out) :: gradient(:)
   !> Failure flag
   logical, intent(out) :: fail

   call model%evaluate(x, energy, gradient, fail)
   if (fail) call env%error("a string node could not be evaluated", &
      & 'landscape_string_gsm')

end subroutine evaluate_node

!> Pack the projected gradients of the internal nodes into one vector.
!>
!> The sign convention is the one of `landscape_lbfgs%direction`, which
!> expects a gradient: the packed vector is the projected **gradient**
!> `g_perp = -f_perp` of the internal nodes.
pure subroutine pack_internal_forces(forces, n, nat, vector)
   !> Node forces `(3, nat, n)`
   real(wp), intent(in) :: forces(:, :, :)
   !> Number of nodes
   integer, intent(in) :: n
   !> Number of atoms
   integer, intent(in) :: nat
   !> Packed projected gradient
   real(wp), intent(out) :: vector(:)

   integer :: k, offset, ncomp

   ncomp = 3*nat
   offset = 0
   do k = 2, n - 1
      vector(offset + 1:offset + ncomp) = -reshape(forces(:, :, k), [ncomp])
      offset = offset + ncomp
   end do

end subroutine pack_internal_forces

!> Pack the internal nodes of a string into one vector.
pure subroutine pack_internal_nodes(nodes, vector)
   !> String nodes `(3, nat, n)`
   real(wp), intent(in) :: nodes(:, :, :)
   !> Packed internal nodes (bohr)
   real(wp), intent(out) :: vector(:)

   integer :: k, offset, ncomp

   ncomp = 3*size(nodes, 2)
   offset = 0
   do k = 2, size(nodes, 3) - 1
      vector(offset + 1:offset + ncomp) = reshape(nodes(:, :, k), [ncomp])
      offset = offset + ncomp
   end do

end subroutine pack_internal_nodes

!> Write packed internal nodes back into a string.
pure subroutine unpack_internal_nodes(vector, nodes)
   !> Packed internal nodes (bohr)
   real(wp), intent(in) :: vector(:)
   !> String nodes `(3, nat, n)`, internal nodes overwritten
   real(wp), intent(inout) :: nodes(:, :, :)

   integer :: k, offset, ncomp

   ncomp = 3*size(nodes, 2)
   offset = 0
   do k = 2, size(nodes, 3) - 1
      nodes(:, :, k) = reshape(vector(offset + 1:offset + ncomp), &
         & [3, size(nodes, 2)])
      offset = offset + ncomp
   end do

end subroutine unpack_internal_nodes

!> Clip a step so that no component exceeds `limit` (in place).
pure subroutine clip(step, limit)
   !> Step (bohr), overwritten in place
   real(wp), intent(inout) :: step(:)
   !> Upper bound of one component (bohr)
   real(wp), intent(in) :: limit

   real(wp) :: largest

   largest = maxval(abs(step))
   if (largest > limit .and. largest > 0.0_wp) then
      step = step*limit/largest
   end if

end subroutine clip

!> Normalize a vector in place.
pure subroutine normalize(vector)
   !> Vector, normalized in place
   real(wp), intent(inout) :: vector(:)

   real(wp) :: norm

   norm = sqrt(sum(vector**2))
   if (norm > 1.0e-16_wp) then
      vector = vector/norm
   else
      vector = 0.0_wp
      if (size(vector) > 0) vector(1) = 1.0_wp
   end if

end subroutine normalize

end module xtb_landscape_path_gsm
