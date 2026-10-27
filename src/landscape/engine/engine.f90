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

!> Engine boundary of the xtb landscape task suite.
!>
!> The engine is the only place where the landscape algorithms touch the
!> electronic structure: it evaluates energy + gradient, the (dense,
!> symmetrised) Hessian and Hessian-vector products through the abstract
!> interfaces below.  Algorithms therefore never call `singlepoint` or
!> `calc%hessian` directly.
!>
!> Two backends can be attached to an engine:
!>
!>  * a `class(TCalculator)` instance (GFN0/GFN1/GFN2/FF/...), evaluated with
!>    the xtb single point driver `xtb_single::singlepoint` (`src/single.f90:37`)
!>  * a `class(landscape_energy_model)` instance, i.e. a pure analytic
!>    potential (harmonic oscillator, Rosenbrock function, Mueller-Brown
!>    surface).  This backend exists so that the optimisers can be unit tested
!>    without running a single SCF, and so that the Hessian provider can be
!>    replaced by a future analytic implementation without touching any
!>    algorithm (PLAN D3).
!>
!> The Hessian is produced by a `class(landscape_hessian_provider)` attached to
!> the engine (`xtb_landscape_hessian`); the default provider is the self-built
!> central difference implementation.  A provider failure (unknown source,
!> unsupported approximation) is turned into `fail = .true.` by probing only the
!> environment messages appended by the provider call.
!>
!> Units: coordinates in bohr, energies in hartree, gradients in hartree/bohr
!> and Hessian elements in hartree/bohr^2.
!>
!> Caching / reproducibility contract:
!>
!>  * `energy_gradient` returns the cached values if and only if the current
!>    geometry is *bitwise* identical to the geometry of the previous call
!>    (`all(mol%xyz == cache_xyz)`) **and** the atom list, electronic
!>    temperature and SCC accuracy are unchanged.  A cache hit does not
!>    increase `n_eval`.
!>  * `set_xyz` and `invalidate` unconditionally drop the cache.
!>  * `n_eval`, `n_cache_hit`, `n_hess`, `n_hvp` and `n_scf_fail` are exposed
!>    for cost accounting and determinism tests.
!>
!> Restart policy (settings `restart_mode`):
!>
!>  * `p_landscape_restart_off` (default): the wavefunction is reset to a fixed
!>    reference before every single point, so the same geometry always gives
!>    bitwise identical energy and gradient (see `07-*-findings.md`, F3).
!>  * `p_landscape_restart_pin`: the wavefunction of the first evaluation is
!>    recorded and reused, `engine_pin_state` may overwrite the record.
!>  * `p_landscape_restart_warm`: the incoming wavefunction is reused (fast,
!>    but reintroduces a SCF-history dependency, measured at ~4e-7 Eh/bohr).
!>
!>  The policy is applied to every energy/gradient evaluation, including the
!>  central differences of the numerical Hessian and of the HVP.  Bitwise
!>  reproducibility of `off` additionally requires a fixed `OMP_NUM_THREADS`:
!>  OpenMP reduction ordering changes the last bits of the SCF (measured at
!>  ~1e-16 Eh/bohr in-process, `docs/_landscape-recon/07-*.md`).
module xtb_landscape_engine
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   use xtb_type_restart, only : TRestart
   use xtb_type_calculator, only : TCalculator
   use xtb_xtb_calculator, only : TxTBCalculator
   use xtb_type_data, only : scc_results
   use xtb_single, only : singlepoint
   use xtb_landscape_settings, only : landscape_engine_control, &
      & p_landscape_restart_off, p_landscape_restart_pin, p_landscape_restart_warm, &
      & landscape_restart_name, landscape_log_failed, landscape_value_present
   use xtb_landscape_coords, only : geometry_pack
   use xtb_landscape_hessian, only : landscape_hessian_provider, &
      & landscape_numerical_hessian_provider, landscape_model_hessian_provider, &
      & new_landscape_hessian_provider
   implicit none(type, external)
   private

   public :: landscape_energy_model
   public :: landscape_hessian_provider, landscape_numerical_hessian_provider
   public :: landscape_model_hessian_provider, new_landscape_hessian_provider
   public :: landscape_engine, new_landscape_engine

   !> Abstract potential energy surface used interchangeably as engine backend
   !> and as objective function of the optimisers.
   !>
   !> `evaluate` returns the flattened gradient (length 3N, ordered
   !> `gx1, gy1, gz1, gx2, ...`), `hessian` the dense, symmetric Hessian
   !> `(3N, 3N)`.  Both routines report failure through the `fail` flag; they
   !> must not terminate the program.
   type, abstract :: landscape_energy_model
   contains
      !> Energy and gradient at the flattened coordinates `x`
      procedure(model_evaluate), deferred :: evaluate
      !> Dense symmetric Hessian at the flattened coordinates `x`
      procedure(model_hessian), deferred :: hessian
   end type landscape_energy_model

   abstract interface
      subroutine model_evaluate(self, x, energy, gradient, fail)
         import :: landscape_energy_model, wp
         !> Potential instance
         class(landscape_energy_model), intent(inout) :: self
         !> Flattened coordinates (bohr)
         real(wp), intent(in) :: x(:)
         !> Total energy (hartree)
         real(wp), intent(out) :: energy
         !> Flattened gradient (hartree/bohr)
         real(wp), intent(out) :: gradient(:)
         !> Failure flag
         logical, intent(out) :: fail
      end subroutine model_evaluate

      subroutine model_hessian(self, x, hess, fail)
         import :: landscape_energy_model, wp
         !> Potential instance
         class(landscape_energy_model), intent(inout) :: self
         !> Flattened coordinates (bohr)
         real(wp), intent(in) :: x(:)
         !> Dense Hessian (hartree/bohr^2)
         real(wp), intent(out) :: hess(:, :)
         !> Failure flag
         logical, intent(out) :: fail
      end subroutine model_hessian
   end interface

   !> Abstract Hessian provider (PLAN D3).
   !>
   !> The provider itself lives in `xtb_landscape_hessian`, together with the
   !> factory `new_landscape_hessian_provider` that turns the `hessian_source`
   !> setting into a concrete implementation.  It is re-exported here because
   !> the engine is the only place that calls it.

   !> Landscape engine holding the calculation backend, the Hessian provider,
   !> the evaluation cache and the cost counters.
   type :: landscape_engine
      !> Polymorphic single point calculator (may be unallocated)
      class(TCalculator), allocatable :: calc
      !> Injectable analytic potential (may be unallocated)
      class(landscape_energy_model), allocatable :: model
      !> Hessian provider
      class(landscape_hessian_provider), allocatable :: hp
      !> Run-time control forwarded to the single point driver
      type(landscape_engine_control) :: ctrl
      !> Deterministic reference wavefunction of the restart policy
      type(TRestart) :: wfn_ref
      !> True once `wfn_ref` has been recorded
      logical :: ref_valid = .false.
      !> Geometry of the reference state (bohr)
      real(wp), allocatable :: ref_xyz(:, :)
      !> Current geometry set through `set_xyz` (bohr)
      real(wp), allocatable :: current_xyz(:, :)
      !> Cache validity flag
      logical :: cache_valid = .false.
      !> Geometry of the cached evaluation (bohr)
      real(wp), allocatable :: cache_xyz(:, :)
      !> Atom numbers of the cached evaluation
      integer, allocatable :: cache_at(:)
      !> Electronic temperature of the cached evaluation (kelvin)
      real(wp) :: cache_etemp = 0.0_wp
      !> SCC accuracy of the cached evaluation
      real(wp) :: cache_acc = 0.0_wp
      !> Cached energy (hartree)
      real(wp) :: cache_energy = 0.0_wp
      !> Cached gradient (hartree/bohr)
      real(wp), allocatable :: cache_gradient(:, :)
      !> Cached strain derivatives (hartree)
      real(wp) :: cache_sigma(3, 3) = 0.0_wp
      !> Cached HOMO-LUMO gap (hartree)
      real(wp) :: cache_hlgap = 0.0_wp
      !> Cached detailed results
      type(scc_results) :: cache_results
      !> Number of energy/gradient evaluations (cache misses)
      integer :: n_eval = 0
      !> Number of cache hits
      integer :: n_cache_hit = 0
      !> Number of Hessian evaluations
      integer :: n_hess = 0
      !> Number of Hessian evaluations that reached the numerical provider,
      !> that is the builds that consumed `ctrl%hess_step`.  `n_hess` counts
      !> the model and file providers too, which document the finite
      !> difference step as unused, so the two counters are not the same
      !> number.  This one is the single source of the "(not used: ...)"
      !> notice of `report_unused_hess_step`.
      integer :: n_hess_fd = 0
      !> Number of Hessian-vector products
      integer :: n_hvp = 0
      !> Number of single point calls with unconverged SCC
      integer :: n_scf_fail = 0
   contains
      !> Attach a calculator (deep copy)
      procedure :: set_calculator => engine_set_calculator
      !> Retrieve a copy of the attached calculator
      procedure :: get_calculator => engine_get_calculator
      !> Attach an analytic potential
      procedure :: set_model => engine_set_model
      !> Attach a Hessian provider
      procedure :: set_hessian_provider => engine_set_hessian_provider
      !> Overwrite the current geometry and drop the cache
      procedure :: set_xyz => engine_set_xyz
      !> Read the current geometry
      procedure :: get_xyz => engine_get_xyz
      !> Read the geometry of the last cached evaluation
      procedure :: get_cached_xyz => engine_get_cached_xyz
      !> Record the current geometry and wavefunction as restart reference
      procedure :: pin_state => engine_pin_state
      !> Reset the wavefunction according to the restart policy
      procedure :: apply_restart => engine_apply_restart
      !> Drop the cached evaluation
      procedure :: invalidate => engine_invalidate
      !> Energy and gradient at the current geometry
      procedure :: energy_gradient => engine_energy_gradient
      !> Dense symmetric Hessian at the current geometry
      procedure :: hessian => engine_hessian
      !> Hessian-vector product at the current geometry
      procedure :: hvp => engine_hvp
      !> Print the cost counters
      procedure :: report => engine_report
      !> Report a named finite difference step no build of this run consumed
      procedure :: report_unused_hess_step => engine_report_unused_hess_step
      !> Effective electronic temperature read back from the calculator
      procedure :: effective_etemp => engine_effective_etemp
      !> Release all allocatable components
      procedure :: destroy => engine_destroy
   end type landscape_engine

contains

!> Construct a landscape engine with the numerical Hessian provider.
subroutine new_landscape_engine(self, ctrl)
   !> New engine
   type(landscape_engine), intent(out) :: self
   !> Run-time control
   type(landscape_engine_control), intent(in), optional :: ctrl

   if (present(ctrl)) self%ctrl = ctrl
   allocate(landscape_numerical_hessian_provider :: self%hp)
   self%hp%name = 'numerical'
   self%hp%exact = .false.
   call set_numerical_accuracy(self)

end subroutine new_landscape_engine

!> Copy the engine Hessian accuracy into the numerical provider.
subroutine set_numerical_accuracy(self)
   !> Engine
   class(landscape_engine), intent(inout) :: self

   select type (hp => self%hp)
   type is (landscape_numerical_hessian_provider)
      hp%sccacc = self%ctrl%hess_acc
   end select

end subroutine set_numerical_accuracy

!> Release all allocatable components of an engine.
subroutine engine_destroy(self)
   !> Engine to clean up
   class(landscape_engine), intent(inout) :: self

   if (allocated(self%calc)) deallocate(self%calc)
   if (allocated(self%model)) deallocate(self%model)
   if (allocated(self%hp)) deallocate(self%hp)
   if (allocated(self%cache_xyz)) deallocate(self%cache_xyz)
   if (allocated(self%cache_gradient)) deallocate(self%cache_gradient)
   if (allocated(self%cache_at)) deallocate(self%cache_at)
   if (allocated(self%ref_xyz)) deallocate(self%ref_xyz)
   if (allocated(self%current_xyz)) deallocate(self%current_xyz)
   self%cache_valid = .false.
   self%ref_valid = .false.

end subroutine engine_destroy

!> Attach a single point calculator to the engine (deep copy of the dynamic
!> type) and invalidate the cache.
subroutine engine_set_calculator(self, calc)
   !> Engine
   class(landscape_engine), intent(inout) :: self
   !> Calculator to copy
   class(TCalculator), intent(in) :: calc

   if (allocated(self%calc)) deallocate(self%calc)
   allocate(self%calc, source=calc)
   ! The engine control is the single source of truth for the SCC accuracy and
   ! the electronic temperature: `xtb_single::singlepoint` receives `acc` and
   ! `et` as arguments but drops both (`acc` never reaches the calculator,
   ! `et` is ignored, `src/single.f90:52`), so the values have to be pushed
   ! into the calculator itself here (`landscape_engine_control%acc` and
   ! `%etemp`).
   self%calc%accuracy = self%ctrl%acc
   call push_calculator_etemp(self)
   self%ref_valid = .false.
   call self%invalidate()

end subroutine engine_set_calculator

!> Push the electronic temperature of the engine control into the calculator.
!>
!> An xTB calculator uses its **own** `etemp` field (`scf(...)` is called with
!> `self%etemp`, `src/xtb/calculator.f90:258`), which the constructor copies
!> once from the global `set%etemp` (`src/xtb/calculator.f90:144`) because the
!> single point driver of the suite cannot forward it.  Without this push the
!> control field would be dead and a caller that sets
!> `landscape_engine_control%etemp` would be silently ignored (measured in
!> `docs/_landscape-recon/21-m3-validation-experiments.md`, F1).
!>
!> Only `TxTBCalculator` carries an electronic temperature; the other
!> calculators (GFN-FF, external codes) have no such field and keep their own
!> settings untouched.  The value is verified from the outside by
!> `engine_report`, which prints it next to the SCC accuracy.
subroutine push_calculator_etemp(self)
   !> Engine
   class(landscape_engine), intent(inout) :: self

   select type (calc => self%calc)
   type is (TxTBCalculator)
      calc%etemp = self%ctrl%etemp
   end select

end subroutine push_calculator_etemp

!> Return a copy of the attached calculator, if any.
subroutine engine_get_calculator(self, calc)
   !> Engine
   class(landscape_engine), intent(in) :: self
   !> Copy of the attached calculator
   class(TCalculator), allocatable, intent(out) :: calc

   if (allocated(self%calc)) allocate(calc, source=self%calc)

end subroutine engine_get_calculator

!> Attach an analytic potential and invalidate the cache.
subroutine engine_set_model(self, model)
   !> Engine
   class(landscape_engine), intent(inout) :: self
   !> Potential to copy
   class(landscape_energy_model), intent(in) :: model

   if (allocated(self%model)) deallocate(self%model)
   allocate(self%model, source=model)
   self%ref_valid = .false.
   call self%invalidate()

end subroutine engine_set_model

!> Attach a Hessian provider and invalidate the cache.
subroutine engine_set_hessian_provider(self, hp)
   !> Engine
   class(landscape_engine), intent(inout) :: self
   !> Provider to copy
   class(landscape_hessian_provider), intent(in) :: hp

   if (allocated(self%hp)) deallocate(self%hp)
   allocate(self%hp, source=hp)
   call self%invalidate()

end subroutine engine_set_hessian_provider

!> Overwrite the current geometry.  The caller is responsible for keeping the
!> shape of `xyz` unchanged; the cache is dropped unconditionally.
subroutine engine_set_xyz(self, xyz)
   !> Engine
   class(landscape_engine), intent(inout) :: self
   !> New geometry (bohr)
   real(wp), intent(in) :: xyz(:, :)

   if (.not. allocated(self%current_xyz)) then
      allocate(self%current_xyz(3, size(xyz, 2)))
   else if (size(self%current_xyz, 2) /= size(xyz, 2)) then
      deallocate(self%current_xyz)
      allocate(self%current_xyz(3, size(xyz, 2)))
   end if
   self%current_xyz = xyz
   self%cache_valid = .false.

end subroutine engine_set_xyz

!> Read the current geometry previously stored with `set_xyz`.
subroutine engine_get_xyz(self, xyz)
   !> Engine
   class(landscape_engine), intent(in) :: self
   !> Current geometry (bohr)
   real(wp), allocatable, intent(out) :: xyz(:, :)

   if (allocated(self%current_xyz)) then
      allocate(xyz, source=self%current_xyz)
   else
      allocate(xyz(3, 0))
   end if

end subroutine engine_get_xyz

!> Read the geometry of the last cached evaluation, if available.
subroutine engine_get_cached_xyz(self, xyz)
   !> Engine
   class(landscape_engine), intent(in) :: self
   !> Geometry of the cached evaluation (bohr)
   real(wp), allocatable, intent(out) :: xyz(:, :)

   if (allocated(self%cache_xyz)) then
      allocate(xyz, source=self%cache_xyz)
   else
      allocate(xyz(3, 0))
   end if

end subroutine engine_get_cached_xyz

!> Record the current geometry and wavefunction as the restart reference.
!>
!> Used by the `p_landscape_restart_pin` policy and to make the `off` policy
!> reproduce a well defined state.  The next energy/gradient evaluation will
!> start from this wavefunction.
subroutine engine_pin_state(self, mol, chk)
   !> Engine
   class(landscape_engine), intent(inout) :: self
   !> Molecular structure
   type(TMolecule), intent(in) :: mol
   !> Wavefunction to record
   type(TRestart), intent(in) :: chk

   call self%wfn_ref%copy(chk)
   if (allocated(self%ref_xyz)) deallocate(self%ref_xyz)
   allocate(self%ref_xyz(3, mol%n), source=mol%xyz)
   self%ref_valid = .true.

end subroutine engine_pin_state

!> Drop the cached energy/gradient evaluation.
subroutine engine_invalidate(self)
   !> Engine
   class(landscape_engine), intent(inout) :: self

   self%cache_valid = .false.

end subroutine engine_invalidate

!> Evaluate energy and gradient at the current geometry.
!>
!> A cache hit (bitwise identical geometry, atom numbers, electronic
!> temperature and SCC accuracy) returns the stored values without calling the
!> calculator.  The cache holds **one** evaluation - it is keyed on the *last*
!> geometry, not on a history - so it can only help when the very same point is
!> evaluated twice in a row.  No kernel of the suite does that: an
!> optimisation moves every geometry and a band or an IRC walk interleaves its
!> images/points, so the measured hit count is **0 in all six tasks** (report
!> 218).  What a band does need - its two fixed endpoints, which *are*
!> re-evaluated at an unchanged geometry - is cached by
!> `xtb_landscape_path_neb::landscape_band_endpoint_cache`, not here.  On a
!> cache miss the backend is evaluated exactly once and the result is stored.
!>
!> An unconverged SCC is never returned as a usable energy/gradient: the
!> failure is counted (`n_scf_fail`), reported with `env%warning` and turned
!> into `fail = .true.` so the caller has to decide what to do.  A warning (and
!> not an error) is used at this library boundary because a single displaced
!> geometry may fail to converge while the surrounding minimisation is still
!> meaningful; the task layer escalates any `n_scf_fail > 0` to an error
!> (`tasks/opt.f90`).  The fresh result is deliberately **not** cached, so a
!> later cache hit cannot silently resurrect it.
!>
!> A second, independent failure mode is checked as well: `singlepoint`
!> declares its `results` argument `intent(out)` but returns early without
!> assigning it whenever `env%check` reports a previous error ("Electronic
!> structure method terminated").  `results%converged` is then undefined and a
!> stale `.true.` from an earlier call passed as a successful evaluation -
!> measured: a scan point with a failed SCC reported `energy = 0.0` together
!> with `converged = .true.` (`docs/_landscape-recon/14-m2-validation-`
!> `experiments.md` section 5.3).  An error logged during the call is the
!> reliable signal, exactly like `TEnvironment%check` uses it, so it is tested
!> first and `results%converged` only afterwards (when `results` is defined).
subroutine engine_energy_gradient(self, env, mol, chk, energy, gradient, fail)
   character(len=*), parameter :: source = 'landscape_engine_energy_gradient'
   !> Engine
   class(landscape_engine), intent(inout) :: self
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure (geometry is not modified)
   type(TMolecule), intent(inout) :: mol
   !> Wavefunction data
   type(TRestart), intent(inout) :: chk
   !> Total energy (hartree)
   real(wp), intent(out) :: energy
   !> Cartesian gradient (hartree/bohr)
   real(wp), intent(out) :: gradient(3, mol%n)
   !> Failure flag
   logical, intent(out) :: fail

   real(wp), allocatable :: x(:), gflat(:)
   real(wp) :: egap
   integer :: nlog
   logical :: hit, scf_fail

   fail = .false.

   hit = self%cache_valid
   if (hit) then
      hit = allocated(self%cache_xyz) .and. allocated(self%cache_gradient)
   end if
   if (hit) then
      hit = size(self%cache_xyz, 1) == 3 .and. size(self%cache_xyz, 2) == mol%n
   end if
   if (hit) hit = all(mol%xyz == self%cache_xyz)
   if (hit) then
      hit = allocated(self%cache_at) .and. size(self%cache_at) == mol%n
      if (hit) hit = all(self%cache_at == mol%at)
   end if
   if (hit) hit = self%cache_etemp == self%ctrl%etemp
   if (hit) hit = self%cache_acc == self%ctrl%acc

   if (hit) then
      energy = self%cache_energy
      gradient = self%cache_gradient
      self%n_cache_hit = self%n_cache_hit + 1
      return
   end if

   if (allocated(self%model)) then
      x = geometry_pack(mol%xyz)
      allocate(gflat(size(x)))
      call self%model%evaluate(x, energy, gflat, fail)
      if (fail) then
         call env%error("analytic energy model failed", source)
         return
      end if
      gradient = reshape(gflat, [3, mol%n])
      self%cache_sigma = 0.0_wp
      self%cache_hlgap = 0.0_wp
   else if (allocated(self%calc)) then
      call self%apply_restart(env, mol, chk)
      ! `egap` is intent(inout) and is used as an initial guess by the
      ! calculator; passing the previous gap would leak state across calls and
      ! break the bitwise reproducibility of the `off` policy.
      egap = 0.0_wp
      nlog = env%nLog
      call singlepoint(env, mol, chk, self%calc, egap, self%ctrl%etemp, &
         & self%ctrl%maxiter, self%ctrl%printlevel, .true., self%ctrl%lgrad, &
         & self%ctrl%acc, energy, gradient, self%cache_sigma, self%cache_results)
      self%cache_hlgap = egap
      if (landscape_log_failed(env, nlog)) then
         ! The calculator terminated, so `results` was never assigned and the
         ! `energy` still holds the `0.0` `singlepoint` wrote before its SCF.
         scf_fail = .true.
      else
         scf_fail = .not. self%cache_results%converged
      end if
      if (scf_fail) then
         self%n_scf_fail = self%n_scf_fail + 1
         self%n_eval = self%n_eval + 1
         call env%warning("SCC did not converge at the current geometry; the "//&
            & "energy and gradient are discarded", source)
         fail = .true.
         return
      end if
      if (.not. landscape_value_present(energy) .or. &
            & .not. all(landscape_value_present(gradient))) then
         self%n_scf_fail = self%n_scf_fail + 1
         self%n_eval = self%n_eval + 1
         call env%warning("the single point returned a non-finite energy or "//&
            & "gradient; the result is discarded", source)
         fail = .true.
         return
      end if
      if (self%ctrl%restart_mode == p_landscape_restart_pin .and. .not. self%ref_valid) then
         call self%pin_state(mol, chk)
      end if
   else
      call env%error("landscape engine has neither a calculator nor an energy model", &
         & source)
      fail = .true.
      return
   end if

   self%n_eval = self%n_eval + 1

   if (allocated(self%cache_xyz)) deallocate(self%cache_xyz)
   allocate(self%cache_xyz(3, mol%n), source=mol%xyz)
   if (allocated(self%cache_gradient)) deallocate(self%cache_gradient)
   allocate(self%cache_gradient(3, mol%n), source=gradient)
   if (allocated(self%cache_at)) deallocate(self%cache_at)
   allocate(self%cache_at(mol%n), source=mol%at)
   self%cache_etemp = self%ctrl%etemp
   self%cache_acc = self%ctrl%acc
   if (allocated(self%current_xyz)) deallocate(self%current_xyz)
   allocate(self%current_xyz(3, mol%n), source=mol%xyz)
   self%cache_energy = energy
   self%cache_valid = .true.

end subroutine engine_energy_gradient

!> Apply the configured restart policy before a single point evaluation.
!>
!> `off` and `pin` reset the incoming wavefunction from the stored reference;
!> `warm` keeps it.  For `off` the reference is recorded on the first call,
!> for `pin` it is recorded after the first evaluation (or by
!> `engine_pin_state`).  The policy therefore also governs the evaluation
!> sequence of the numerical Hessian and of the HVP.
subroutine engine_apply_restart(self, env, mol, chk)
   character(len=*), parameter :: source = 'landscape_engine_apply_restart'
   !> Engine
   class(landscape_engine), intent(inout) :: self
   !> Computational environment for the unknown-mode report
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure (only used to record the reference)
   type(TMolecule), intent(in) :: mol
   !> Wavefunction data
   type(TRestart), intent(inout) :: chk

   select case (self%ctrl%restart_mode)
   case (p_landscape_restart_off)
      if (.not. self%ref_valid) call self%pin_state(mol, chk)
      call chk%copy(self%wfn_ref)
   case (p_landscape_restart_pin)
      if (self%ref_valid) call chk%copy(self%wfn_ref)
   case (p_landscape_restart_warm)
      ! keep the incoming wavefunction
   case default
      ! An unknown mode must not silently degrade to `warm`; the caller sees
      ! the error through the environment and aborts the run.
      call env%error("unknown landscape restart mode: "//&
         & trim(landscape_restart_name(self%ctrl%restart_mode)), source)
   end select

end subroutine engine_apply_restart

!> Evaluate the dense, symmetric Hessian at the current geometry.
!>
!> An attached analytic potential takes precedence; otherwise the Hessian
!> provider (default: central finite differences) is used.  The result is
!> always symmetrised before it is returned.
subroutine engine_hessian(self, env, mol, chk, hess, fail, step, final_err)
   character(len=*), parameter :: source = 'landscape_engine_hessian'
   !> Engine
   class(landscape_engine), intent(inout) :: self
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure (geometry is not modified)
   type(TMolecule), intent(inout) :: mol
   !> Wavefunction data
   type(TRestart), intent(inout) :: chk
   !> Dense Hessian (3N, 3N) in hartree/bohr^2
   real(wp), intent(out) :: hess(3*mol%n, 3*mol%n)
   !> Failure flag
   logical, intent(out) :: fail
   !> Finite difference step (bohr), defaults to `ctrl%hess_step`
   real(wp), intent(in), optional :: step
   !> Residual error estimate of the numerical Hessian
   real(wp), intent(out), optional :: final_err

   real(wp), allocatable :: x(:)
   real(wp) :: step_used, err
   integer :: nlog

   fail = .false.
   err = 0.0_wp
   step_used = self%ctrl%hess_step
   if (present(step)) step_used = step

   if (allocated(self%model)) then
      x = geometry_pack(mol%xyz)
      call self%model%hessian(x, hess, fail)
      if (fail) then
         call env%error("analytic Hessian model failed", source)
         return
      end if
   else if (allocated(self%calc)) then
      if (.not. allocated(self%hp)) then
         call env%error("landscape engine has no Hessian provider attached", source)
         fail = .true.
         return
      end if
      call set_numerical_accuracy(self)
      ! Start the differentiation sequence from the configured restart state so
      ! that the whole Hessian is evaluated under a single, consistent policy.
      call self%apply_restart(env, mol, chk)
      ! Providers add to the provided matrix, so start from zero
      hess = 0.0_wp
      ! Probe only the messages appended by the provider: a provider that
      ! refuses (unknown source, unsupported approximation) must not be
      ! mistaken for a successful call.
      nlog = env%nLog
      call self%hp%compute(env, mol, chk, self%calc, hess, step_used, &
         & self%ctrl%odlr, err)
      if (landscape_log_failed(env, nlog)) then
         fail = .true.
         return
      end if
      ! The provider evaluates its displaced geometries with the calculator
      ! directly (bypassing the engine cache), so its single points have to be
      ! added to the cost counters here.  Without this the reported "energy
      ! evaluations" would grossly underestimate a numerical Hessian (6N).
      self%n_eval = self%n_eval + self%hp%n_eval
      ! The finite difference step reached a provider exactly when the numeric
      ! provider computed: the `model` and `file` providers document their
      ! `step` argument as unused (`types/hessian.f90`).  The counter is the
      ! same fact as the `Hessians` line of `report`, so the unused notice
      ! cannot contradict the printed count.
      if (allocated(self%hp%name)) then
         if (self%hp%name == 'numerical') self%n_hess_fd = self%n_hess_fd + 1
      end if
   else
      call env%error("landscape engine has neither a calculator nor an energy model", &
         & source)
      fail = .true.
      return
   end if

   hess = 0.5_wp*(hess + transpose(hess))
   if (present(final_err)) final_err = err
   self%n_hess = self%n_hess + 1

end subroutine engine_hessian

!> Finite differences Hessian-vector product `H v`.
!>
!> `hv = (g(x + delta*v) - g(x - delta*v)) / (2 delta)`; the original geometry
!> is restored and the cache dropped before returning.  Each call costs two
!> energy/gradient evaluations.
subroutine engine_hvp(self, env, mol, chk, vector, hvec, fail, delta)
   !> Engine
   class(landscape_engine), intent(inout) :: self
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure, restored on return
   type(TMolecule), intent(inout) :: mol
   !> Wavefunction data
   type(TRestart), intent(inout) :: chk
   !> Direction to probe (dimensionless, length 3N)
   real(wp), intent(in) :: vector(:)
   !> Hessian-vector product (hartree/bohr^2)
   real(wp), intent(out) :: hvec(:)
   !> Failure flag
   logical, intent(out) :: fail
   !> Finite difference step (bohr), defaults to `ctrl%hvp_step`
   real(wp), intent(in), optional :: delta

   real(wp), allocatable :: x0(:, :), gplus(:, :), gminus(:, :)
   real(wp) :: step, eplus, eminus
   integer :: i, ic

   fail = .false.
   step = self%ctrl%hvp_step
   if (present(delta)) step = delta

   allocate(x0(3, mol%n), gplus(3, mol%n), gminus(3, mol%n))
   x0 = mol%xyz

   do i = 1, mol%n
      ic = 3*(i - 1)
      mol%xyz(:, i) = x0(:, i) + step*vector(ic + 1:ic + 3)
   end do
   call self%invalidate()
   call self%energy_gradient(env, mol, chk, eplus, gplus, fail)
   if (fail) then
      mol%xyz = x0
      call self%invalidate()
      if (allocated(self%current_xyz)) self%current_xyz = x0
      return
   end if

   do i = 1, mol%n
      ic = 3*(i - 1)
      mol%xyz(:, i) = x0(:, i) - step*vector(ic + 1:ic + 3)
   end do
   call self%invalidate()
   call self%energy_gradient(env, mol, chk, eminus, gminus, fail)
   if (fail) then
      mol%xyz = x0
      call self%invalidate()
      if (allocated(self%current_xyz)) self%current_xyz = x0
      return
   end if

   mol%xyz = x0
   call self%invalidate()
   ! Restore the reported current geometry as well: the two displaced
   ! evaluations above left `current_xyz` at `x0 - step*vector`.
   if (allocated(self%current_xyz)) self%current_xyz = x0
   hvec = reshape((gplus - gminus)/(2.0_wp*step), shape(hvec))
   self%n_hvp = self%n_hvp + 1

end subroutine engine_hvp

!> Effective electronic temperature of the engine (kelvin).
!>
!> The value is read back from the xTB calculator whenever one is attached, so
!> it is the temperature the SCF actually used (the one pushed by
!> `push_calculator_etemp`); without a calculator the control value is
!> returned.  The task layer records this value in its result object instead
!> of echoing the setting, so a missing push is visible from the result as well
!> (M4 hard constraint 1).
function engine_effective_etemp(self) result(etemp)
   !> Engine
   class(landscape_engine), intent(in) :: self
   !> Effective temperature (kelvin)
   real(wp) :: etemp

   etemp = self%ctrl%etemp
   if (allocated(self%calc)) then
      select type (calc => self%calc)
      type is (TxTBCalculator)
         etemp = calc%etemp
      end select
   end if

end function engine_effective_etemp

!> Print the cost counters of the engine.
subroutine engine_report(self, unit)
   !> Engine
   class(landscape_engine), intent(in) :: self
   !> Output unit
   integer, intent(in) :: unit

   write(unit, '(1x,a,i0)') "landscape engine: energy evaluations ", self%n_eval
   write(unit, '(1x,a,i0)') "landscape engine: cache hits         ", self%n_cache_hit
   write(unit, '(1x,a,i0)') "landscape engine: Hessians           ", self%n_hess
   write(unit, '(1x,a,i0)') "landscape engine: Hessian-vector     ", self%n_hvp
   write(unit, '(1x,a,i0)') "landscape engine: unconverged SCC    ", self%n_scf_fail
   write(unit, '(1x,a,a)') "landscape engine: restart policy     ", &
      & landscape_restart_name(self%ctrl%restart_mode)
   ! The effective SCC accuracy actually pushed into the calculator; printing it
   ! makes a missing CLI sync visible from the outside.
   write(unit, '(1x,a,f15.6)') "landscape engine: SCC accuracy       ", &
      & self%ctrl%acc
   ! The electronic temperature actually in effect: read back from the
   ! calculator (`push_calculator_etemp`), so a control value that did not
   ! reach it shows up here.
   if (allocated(self%calc)) then
      select type (calc => self%calc)
      type is (TxTBCalculator)
         write(unit, '(1x,a,f12.3,a)') "landscape engine: etemp              ", &
            & calc%etemp, " K"
      end select
   end if
   if (allocated(self%hp)) then
      write(unit, '(1x,a,a)') "landscape engine: hessian source     ", &
         & self%hp%name
   end if

end subroutine engine_report

!> Report a named finite difference knob that no Hessian build consumed.
!>
!> The notice comes from the counter `n_hess_fd` (the builds that actually
!> reached the numerical provider) and from the provider the engine used, so
!> it states the same fact as the `Hessians` line of `report` and cannot
!> contradict it: a non-zero `n_hess_fd` silences it, and a zero one only
!> occurs when the step was handed to no finite difference build at all.  The
!> claim is about the observed run, not about the selected method, so it stays
!> true when the run stops before its first Hessian (budget exhaustion): a
!> method that would have built one later simply reports that this run built
!> none.  Nothing is printed when the user did not name a flag of the knobs,
!> because the notice is about a knob the user asked for.
!>
!> The step (`hess_step`), the SCC accuracy (`hess_acc`) and the worker count
!> (`hess_jobs`) are one family: the numerical provider consumes all three and
!> the `model`/`file` providers
!> ignore all three, so a single notice covers them and the gate has to accept a
!> flag of **any** of them (`hess_step_source` **or** `hess_acc_source` **or**
!> `hess_jobs_source`).
!> Looking only at `hess_step_source` left `--lts-sccacc` / `--lvib-sccacc`
!> silent when they were the only knob given (`38-stage2b-fix2-review.md`
!> LOW-1).
subroutine engine_report_unused_hess_step(self, unit)
   !> Engine
   class(landscape_engine), intent(in) :: self
   !> Output unit
   integer, intent(in) :: unit

   character(len=:), allocatable :: reason

   reason = ""
   if (self%n_hess_fd > 0) then
      return
   else if (.not. allocated(self%ctrl%hess_step_source) .and. &
         &  .not. allocated(self%ctrl%hess_acc_source) .and. &
         &  .not. allocated(self%ctrl%hess_jobs_source)) then
      return
   else if (allocated(self%hp)) then
      if (allocated(self%hp%name)) then
         if (self%hp%name /= 'numerical') then
            reason = "the Hessian source '"//trim(self%hp%name)//&
               & "' has no finite differences"
         end if
      end if
   end if
   if (len(reason) == 0) then
      reason = "no numerical Hessian was built by this run"
   end if
   write(unit, '(1x,a,a,a)') "(not used: ", reason, ")"

end subroutine engine_report_unused_hess_step

end module xtb_landscape_engine
