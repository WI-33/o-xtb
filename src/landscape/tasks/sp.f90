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

!> Single point task of the landscape suite (energy, gradient, optional Hessian).
!>
!> This is the task level entry point of the `sp` task: it builds the engine,
!> evaluates energy and gradient through the engine cache and, if requested,
!> the Hessian through the configured provider.  Like `tasks/opt.f90` it is one
!> of the two places of the landscape suite that hold the xtb private types
!> (`TMolecule`, `TRestart`, `TCalculator`), keeping the lower layers free of
!> them (PLAN section 4.2).
!>
!> Units: coordinates in bohr, energies in hartree, gradients in hartree/bohr
!> and Hessian elements in hartree/bohr^2.
module xtb_landscape_sp
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : autokcal
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   use xtb_type_restart, only : TRestart
   use xtb_type_calculator, only : TCalculator
   use xtb_landscape_settings, only : landscape_sp_settings, landscape_restart_name
   use xtb_landscape_results, only : landscape_result
   use xtb_landscape_json, only : landscape_json_document, &
      & landscape_json_document_init, landscape_json_add_result, &
      & landscape_json_add_core, landscape_json_add_value, &
      & landscape_json_write, &
      & landscape_json_status_ok, landscape_json_status_error
   use xtb_landscape_coords, only : vec_norm, vec_maxabs, geometry_pack
   use xtb_landscape_engine, only : landscape_engine, new_landscape_engine, &
      & landscape_hessian_provider, new_landscape_hessian_provider
   implicit none(type, external)
   private

   public :: landscape_single_point

contains

!> Task level entry point: single point of `mol` with the calculator `calc`.
!>
!> The engine owns a copy of the calculator and caches the evaluation.  On
!> return `result` holds the energy, the gradient statistics and the cost
!> counters; the Hessian is reported but not stored in the result (the
!> vibrational analysis that consumes it is a follow-up milestone).
subroutine landscape_single_point(env, mol, chk, calc, settings, result)
   character(len=*), parameter :: source = 'landscape_single_point'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure (geometry is not modified)
   type(TMolecule), intent(inout) :: mol
   !> Wavefunction data
   type(TRestart), intent(inout) :: chk
   !> Single point calculator
   class(TCalculator), intent(in) :: calc
   !> Single point settings
   type(landscape_sp_settings), intent(in) :: settings
   !> Single point result
   type(landscape_result), intent(inout) :: result

   type(landscape_engine) :: engine
   class(landscape_hessian_provider), allocatable :: hp
   real(wp), allocatable :: gradient(:, :), hess(:, :)
   real(wp) :: energy, trace, asym
   logical :: fail
   integer :: i, nat

   nat = mol%n
   allocate(gradient(3, nat))

   call new_landscape_engine(engine, settings%engine)
   call new_landscape_hessian_provider(env, settings%hessian_source, hp, &
      & sccacc=settings%engine%hess_acc, filename=settings%hessian_file, &
      & printlevel=settings%printlevel)
   if (.not. allocated(hp)) return
   call engine%set_hessian_provider(hp)
   call engine%set_calculator(calc)

   if (settings%printlevel > 0) then
      write(env%unit, '(1x,a)') "landscape: single point"
      write(env%unit, '(1x,a,a)') "restart policy:           ", &
         & landscape_restart_name(settings%engine%restart_mode)
      write(env%unit, '(1x,a,a)') "Hessian source:           ", &
         & trim(settings%hessian_source)
   end if

   call engine%energy_gradient(env, mol, chk, energy, gradient, fail)
   if (fail) then
      call env%error("landscape single point failed", source)
      call engine%destroy()
      return
   end if

   result%energy = energy
   result%gnorm = vec_norm(geometry_pack(gradient))
   result%gmax = vec_maxabs(geometry_pack(gradient))
   result%n_eval = engine%n_eval
   result%n_cache_miss = engine%n_eval
   result%converged = engine%n_scf_fail == 0
   if (allocated(result%message)) deallocate(result%message)
   if (result%converged) then
      result%message = "single point converged"
   else
      result%message = "single point with unconverged SCC"
   end if

   if (settings%printlevel > 0) then
      write(env%unit, '(1x,a,f20.12,a)') "energy:          ", energy, " Eh"
      write(env%unit, '(1x,a,f20.12,a)') "energy:          ", &
         & energy*autokcal, " kcal/mol"
      write(env%unit, '(1x,a,es15.6,a)') "gradient norm:   ", result%gnorm, &
         & " Eh/bohr"
      write(env%unit, '(1x,a,es15.6,a)') "max gradient:    ", result%gmax, &
         & " Eh/bohr"
   end if

   if (settings%write_gradient) call write_gradient(env%unit, gradient, &
      & settings%grad_threshold)

   if (settings%write_hessian) then
      allocate(hess(3*nat, 3*nat), source=0.0_wp)
      call engine%hessian(env, mol, chk, hess, fail)
      if (fail) then
         call env%error("landscape single point Hessian failed", source)
         call engine%destroy()
         return
      end if
      trace = 0.0_wp
      do i = 1, 3*nat
         trace = trace + hess(i, i)
      end do
      asym = maxval(abs(hess - transpose(hess)))
      if (settings%printlevel > 0) then
         write(env%unit, '(1x,a,i0,a,i0,a)') "Hessian:         (", 3*nat, &
            & ", ", 3*nat, ")"
         write(env%unit, '(1x,a,es15.6,a)') "Hessian trace:   ", trace, &
            & " Eh/bohr^2"
         write(env%unit, '(1x,a,es15.6,a)') "asymmetry:       ", asym, &
            & " Eh/bohr^2"
      end if
   end if

   if (settings%printlevel > 0) call engine%report(env%unit)
   call engine%destroy()

   if (settings%ljson) call write_json(env, mol, settings, result)

end subroutine landscape_single_point

!> Write the structured result `landscape_sp.json`.
!>
!> The document carries the shared result fields plus the single point
!> specific request flags, so a regression diff shows whether a Hessian was
!> asked for as well.
subroutine write_json(env, mol, settings, result)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecule of the run (the core block's geometry)
   type(TMolecule), intent(in) :: mol
   !> Single point settings
   type(landscape_sp_settings), intent(in) :: settings
   !> Single point result
   type(landscape_result), intent(in) :: result

   type(landscape_json_document) :: doc = landscape_json_document()
   logical :: ok
   integer :: status

   status = landscape_json_status_ok
   ! A single point whose SCC did not converge is an error: it has no usable
   ! result, and `landscape_single_point` already reported it through `env`.
   if (.not. result%converged) status = landscape_json_status_error
   call landscape_json_document_init(doc, 'sp')
   call landscape_json_add_result(doc, result, status)
   call landscape_json_add_core(doc, mol, result%energy, settings%input_file)
   call landscape_json_add_value(doc, 'kcal_per_mol', result%energy*autokcal)
   call landscape_json_write(env, doc, ok)

end subroutine write_json

!> Print the gradient components above `threshold`.
!>
!> A threshold of zero prints every component (the MAPLE `SPParams.verbose=1`
!> behaviour, `dispatcher/sp/sp.py`).
subroutine write_gradient(unit, gradient, threshold)
   !> Output unit
   integer, intent(in) :: unit
   !> Cartesian gradient (hartree/bohr), shape (3, natoms)
   real(wp), intent(in) :: gradient(:, :)
   !> Print threshold (hartree/bohr)
   real(wp), intent(in) :: threshold

   integer :: iat, ic

   write(unit, '(1x,a)') "gradient (Eh/bohr), component order x, y, z"
   do iat = 1, size(gradient, 2)
      do ic = 1, 3
         if (abs(gradient(ic, iat)) >= threshold) then
            write(unit, '(1x,a,i0,a,i0,a,es15.6)') "atom ", iat, " component ", &
               & ic, ": ", gradient(ic, iat)
         end if
      end do
   end do

end subroutine write_gradient

end module xtb_landscape_sp
