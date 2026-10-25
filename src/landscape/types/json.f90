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

!> Structured result output of the landscape suite (`--ljson`).
!>
!> `--ljson` is the **single machine readable entry point for every task**
!> (report 117): the six suite tasks and the native run types (`--sp`,
!> `--opt`, `--hess`, `--ohess`, `--md`, ...) all write one
!> `landscape_<task>.json` document through this module.  The name keeps the
!> historical `landscape_` prefix for the native tasks too, so a single parser
!> finds every task under one naming convention.
!>
!> Every document carries a **uniform core block** whose keys and types are
!> the same across tasks, so one parser can read all of them:
!>
!>  * `format` (= `2`, the schema version a parser branches on);
!>  * `task`, `status`, `status_name`, `converged`, `message`;
!>  * `n_atoms` (integer) and `energy` (hartree, the canonical total energy);
!>  * `geometry` (an object with `unit` = `bohr`, `symbols` and
!>    `coordinates`, identical in shape for every task);
!>  * `provenance` (an object with `oxtb_version`, `engine_version`, `commit`,
!>    `command_line`, `method`, `charge`, `spin` and `input_file`).
!>
!> The task specific data stays in `counters`/`values`/`flags`/`history`, so a
!> parser can ignore what it does not know without failing.  The `format`
!> bump from 1 to 2 and the four new top level keys (`n_atoms`, `energy`,
!> `geometry`, `provenance`) are **additive**: every pre-existing key keeps its
!> name, type and meaning.
!>
!> The format is otherwise deliberately frozen:
!>
!>  * the top level keys are always the same and always in the same order
!>    (`format`, `task`, `status`, `status_name`, `converged`, `message`,
!>    `n_atoms`, `energy`, `geometry`, `provenance`, `counters`, `values`,
!>    `flags`, `history`);
!>  * all counters and values are written in the order the task added them, so
!>    a diff of two files of the same task is stable;
!>  * real numbers are written with 17 significant digits (`es25.17e3`), enough
!>    for a bitwise double round trip, and the conversion is deterministic for
!>    a fixed compiler, so two runs of the same input **at the same
!>    `OMP_NUM_THREADS`** produce byte identical files (different thread counts
!>    differ by up to 2.4e-10, PLAN section 12.16).  **No wall time, no path and
!>    no random quantity is written** -
!>    every one of them would break that property.  The command line in
!>    `provenance` is not a counterexample: it is the caller's own input, not a
!>    property of this process.
!>
!> Missing values are written as JSON `null`, never as a number: JSON has no
!> `NaN`/`Infinity` literal and a fabricated `0.0` would survive unnoticed in
!> a numerical diff (the rule PLAN section 5.5 established for the scan table).
!>
!> The counters are the *same* numbers the task prints on stdout, they are not
!> a second bookkeeping: each task adds them from its own result object, and
!> the CLI validation asserts the two agree
!> (`test/unit/landscape_e2e_check.py`, mode `ljson`).
!>
!> Implementation note: the document holds only intrinsic arrays.  An array
!> of a derived type with allocatable components would be the natural
!> container for the history series, but gfortran 12.2.0 has a stride defect
!> in the intrinsic assignment of exactly that shape
!> (`docs/_landscape-recon/10-m0-m1-fix-round-report.md` section 4.1); the
!> series are therefore stored as one flat real pool plus an offset/length
!> table per series, which is both safe and deterministic.
!>
!> Units: whatever the task documents; the file carries no unit strings, the
!> key names are the ones of the result types (`energy` in hartree, `gnorm` in
!> hartree/bohr, `etemp` in kelvin, ...).
module xtb_landscape_json
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   use xtb_mctc_version, only : version, engine_version, oxtb_version
   use xtb_setmod, only : p_gfn_method_ff
   use xtb_setparam, only : get_namespace, set
   use xtb_landscape_results, only : landscape_result, &
      & landscape_convergence_report
   implicit none(type, external)
   private

   public :: landscape_json_document
   public :: landscape_json_document_init
   public :: landscape_json_add_counter, landscape_json_add_value
   public :: landscape_json_add_flag, landscape_json_add_series
   public :: landscape_json_add_history, landscape_json_add_result
   public :: landscape_json_add_core
   public :: landscape_json_write, landscape_json_write_file
   public :: landscape_json_write_run
   public :: landscape_json_status_ok, landscape_json_status_incomplete
   public :: landscape_json_status_error, landscape_json_status_name

   !> Status of a landscape task in the JSON document: the task finished and
   !> its result is usable.
   integer, parameter :: landscape_json_status_ok = 0
   !> Status: the task finished but did not meet every criterion (an
   !> unconverged minimisation, a scan with unconverged points).  The result
   !> is still written, the caller decides whether it is usable.
   integer, parameter :: landscape_json_status_incomplete = 1
   !> Status: the task did not produce a usable result (an error was reported
   !> through `env`; the process exits non-zero).
   integer, parameter :: landscape_json_status_error = -1

   !> Maximum key length of a counter/value/flag/series name.
   integer, parameter :: landscape_json_key_len = 40
   !> Length of the task name field.
   integer, parameter :: landscape_json_task_len = 32
   !> Length of the message field; a longer termination message is truncated.
   !> Every message of the suite is far below this bound.
   integer, parameter :: landscape_json_message_len = 256

   !> Provenance of one run: everything a reader needs to trace a document
   !> back to the build, the invocation and the electronic structure setup.
   !>
   !> The two version strings are the ones the banner of `61fb57c` introduced
   !> (`oxtb_version`, the product version, and `engine_version`, the xTB
   !> engine it is built from); the commit is read back out of the `version`
   !> constant (the upstream `6.7.1 (<commit>)` spelling), never re-derived.
   !> `command_line` is the verbatim invocation (`get_command`), so a document
   !> is reproducible from the file alone.  `charge` and `spin` are the setup
   !> of the run (`mol%chrg`, `mol%uhf`).  The fields are all written, so the
   !> `provenance` block has the same shape in every task document.
   type :: landscape_provenance
      !> Product version of the program (`oxtb_version`, 0.3.2)
      character(len=:), allocatable :: oxtb_version
      !> Version of the xTB engine (`engine_version`, 6.7.1)
      character(len=:), allocatable :: engine_version
      !> Commit hash the binary was built from (from `version`)
      character(len=:), allocatable :: commit
      !> Verbatim command line of the run
      character(len=:), allocatable :: command_line
      !> Electronic structure method, `GFN0-xTB`/`GFN1-xTB`/`GFN2-xTB`/`GFN-FF`
      character(len=:), allocatable :: method
      !> Positional geometry input file (empty when the run has none)
      character(len=:), allocatable :: input_file
      !> Total charge of the run
      integer :: charge = 0
      !> Number of unpaired electrons (spin) of the run
      integer :: spin = 0
   end type landscape_provenance

   !> A task result as an ordered set of named scalars, series and a status.
   !>
   !> The document is a simple builder: the task appends counters, values,
   !> flags and series in a fixed order and `landscape_json_write` serialises
   !> them.  No module level state is involved, so two documents never
   !> interfere.
   !>
   !> The history series live in a single real pool (`series_value`); series
   !> `k` occupies `series_value(offset(k) : offset(k) + length(k) - 1)`.
   type :: landscape_json_document
      !> Task name, part of the file name (`landscape_<task>.json`)
      character(len=landscape_json_task_len) :: task = ''
      !> Status code, one of `landscape_json_status_*`
      integer :: status = landscape_json_status_error
      !> True if the task met every convergence criterion
      logical :: converged = .false.
      !> Human readable termination message
      character(len=landscape_json_message_len) :: message = ''
      !> Integer counters (names and values, insertion order)
      character(len=landscape_json_key_len), allocatable :: counter_key(:)
      integer, allocatable :: counter_value(:)
      !> Number of counters added
      integer :: n_counter = 0
      !> Real values (names and values, insertion order)
      character(len=landscape_json_key_len), allocatable :: value_key(:)
      real(wp), allocatable :: value_value(:)
      !> Number of values added
      integer :: n_value = 0
      !> Logical flags (names and values, insertion order)
      character(len=landscape_json_key_len), allocatable :: flag_key(:)
      logical, allocatable :: flag_value(:)
      !> Number of flags added
      integer :: n_flag = 0
      !> Names of the history series
      character(len=landscape_json_key_len), allocatable :: series_key(:)
      !> Start index of each series in `series_value`
      integer, allocatable :: series_offset(:)
      !> Number of values of each series
      integer, allocatable :: series_length(:)
      !> Flat pool holding every series back to back
      real(wp), allocatable :: series_value(:)
      !> Number of series added
      integer :: n_series = 0
      !> Number of atoms of the geometry block
      integer :: n_atoms = 0
      !> Canonical total energy of the run (hartree).  It duplicates the
      !> task specific energy inside `values` on purpose: a parser must not
      !> have to know a task specific key (`values%energy`, `values%energy_ts`)
      !> to find the energy of a document, so the core key is the one to read.
      real(wp) :: energy = 0.0_wp
      !> Atom symbols of the geometry block
      character(len=4), allocatable :: geometry_symbol(:)
      !> Cartesian coordinates (bohr, shape `(3, n_atoms)`)
      real(wp), allocatable :: geometry_xyz(:, :)
      !> Provenance of the run (version, command line, method, charge, spin)
      type(landscape_provenance) :: provenance
   end type landscape_json_document

contains

!> Initialise a document for one task.
subroutine landscape_json_document_init(self, task)
   !> Document to initialise
   type(landscape_json_document), intent(out) :: self
   !> Task name (`sp`, `opt`, `scan`, `vib`, `ts`, `irc`)
   character(len=*), intent(in) :: task

   self%task = task
   self%status = landscape_json_status_error
   self%converged = .false.
   self%message = ''
   allocate(self%counter_key(32))
   self%counter_key = ''

   allocate(self%counter_value(32), source=0)
   allocate(self%value_key(32))
   self%value_key = ''

   allocate(self%value_value(32), source=0.0_wp)
   allocate(self%flag_key(16))
   self%flag_key = ''

   allocate(self%flag_value(16), source=.false.)
   allocate(self%series_key(16))
   self%series_key = ''

   allocate(self%series_offset(16), source=0)
   allocate(self%series_length(16), source=0)
   allocate(self%series_value(0), source=0.0_wp)

   ! The provenance strings are allocated (to the empty string) so the writer
   ! can always emit the block with the same shape, even for a document that
   ! was built without `landscape_json_add_core` (the unit test drives the
   ! writer directly).
   self%provenance%oxtb_version = ''
   self%provenance%engine_version = ''
   self%provenance%commit = ''
   self%provenance%command_line = ''
   self%provenance%method = ''
   self%provenance%input_file = ''

end subroutine landscape_json_document_init

!> Append an integer counter.
subroutine landscape_json_add_counter(self, key, value)
   !> Document to extend
   type(landscape_json_document), intent(inout) :: self
   !> Counter name
   character(len=*), intent(in) :: key
   !> Counter value
   integer, intent(in) :: value

   integer :: capacity
   character(len=landscape_json_key_len), allocatable :: cbuffer(:)
   integer, allocatable :: ibuffer(:)

   self%n_counter = self%n_counter + 1
   if (self%n_counter > size(self%counter_key)) then
      capacity = 2*size(self%counter_key)
      allocate(cbuffer(capacity))
      cbuffer = ''
      cbuffer(:size(self%counter_key)) = self%counter_key
      call move_alloc(cbuffer, self%counter_key)
      allocate(ibuffer(capacity), source=0)
      ibuffer(:size(self%counter_value)) = self%counter_value
      call move_alloc(ibuffer, self%counter_value)
   end if
   self%counter_key(self%n_counter) = key
   self%counter_value(self%n_counter) = value

end subroutine landscape_json_add_counter

!> Append a real value.
subroutine landscape_json_add_value(self, key, value)
   !> Document to extend
   type(landscape_json_document), intent(inout) :: self
   !> Value name
   character(len=*), intent(in) :: key
   !> Value (in the unit the task documents)
   real(wp), intent(in) :: value

   integer :: capacity
   character(len=landscape_json_key_len), allocatable :: cbuffer(:)
   real(wp), allocatable :: rbuffer(:)

   self%n_value = self%n_value + 1
   if (self%n_value > size(self%value_key)) then
      capacity = 2*size(self%value_key)
      allocate(cbuffer(capacity))
      cbuffer = ''
      cbuffer(:size(self%value_key)) = self%value_key
      call move_alloc(cbuffer, self%value_key)
      allocate(rbuffer(capacity), source=0.0_wp)
      rbuffer(:size(self%value_value)) = self%value_value
      call move_alloc(rbuffer, self%value_value)
   end if
   self%value_key(self%n_value) = key
   self%value_value(self%n_value) = value

end subroutine landscape_json_add_value

!> Append a logical flag.
subroutine landscape_json_add_flag(self, key, value)
   !> Document to extend
   type(landscape_json_document), intent(inout) :: self
   !> Flag name
   character(len=*), intent(in) :: key
   !> Flag value
   logical, intent(in) :: value

   integer :: capacity
   character(len=landscape_json_key_len), allocatable :: cbuffer(:)
   logical, allocatable :: lbuffer(:)

   self%n_flag = self%n_flag + 1
   if (self%n_flag > size(self%flag_key)) then
      capacity = 2*size(self%flag_key)
      allocate(cbuffer(capacity))
      cbuffer = ''
      cbuffer(:size(self%flag_key)) = self%flag_key
      call move_alloc(cbuffer, self%flag_key)
      allocate(lbuffer(capacity), source=.false.)
      lbuffer(:size(self%flag_value)) = self%flag_value
      call move_alloc(lbuffer, self%flag_value)
   end if
   self%flag_key(self%n_flag) = key
   self%flag_value(self%n_flag) = value

end subroutine landscape_json_add_flag

!> Append a real series to the convergence history.
subroutine landscape_json_add_series(self, key, values)
   !> Document to extend
   type(landscape_json_document), intent(inout) :: self
   !> Series name
   character(len=*), intent(in) :: key
   !> Series values in iteration order
   real(wp), intent(in) :: values(:)

   integer :: capacity, offset
   character(len=landscape_json_key_len), allocatable :: cbuffer(:)
   integer, allocatable :: ibuffer(:)
   real(wp), allocatable :: rbuffer(:)

   if (self%n_series > 0) then
      offset = self%series_offset(self%n_series) + &
         & self%series_length(self%n_series)
   else
      offset = 1
   end if

   self%n_series = self%n_series + 1
   if (self%n_series > size(self%series_key)) then
      capacity = 2*size(self%series_key)
      allocate(cbuffer(capacity))
      cbuffer = ''
      cbuffer(:size(self%series_key)) = self%series_key
      call move_alloc(cbuffer, self%series_key)
      allocate(ibuffer(capacity), source=0)
      ibuffer(:size(self%series_offset)) = self%series_offset
      call move_alloc(ibuffer, self%series_offset)
      allocate(ibuffer(capacity), source=0)
      ibuffer(:size(self%series_length)) = self%series_length
      call move_alloc(ibuffer, self%series_length)
   end if

   allocate(rbuffer(offset - 1 + size(values)), source=0.0_wp)
   if (offset > 1) rbuffer(:offset - 1) = self%series_value
   rbuffer(offset:) = values
   call move_alloc(rbuffer, self%series_value)

   self%series_key(self%n_series) = key
   self%series_offset(self%n_series) = offset
   self%series_length(self%n_series) = size(values)

end subroutine landscape_json_add_series

!> Append the convergence history of a result.
!>
!> The series are the ones of `landscape_convergence_report`, in a fixed
!> order; the iteration count is added as the counter `n_iter`, so a reader
!> can detect a truncated series without trusting the array lengths.
subroutine landscape_json_add_history(self, conv)
   !> Document to extend
   type(landscape_json_document), intent(inout) :: self
   !> Convergence report (may hold no history at all)
   type(landscape_convergence_report), intent(in) :: conv

   call landscape_json_add_counter(self, 'n_iter', conv%n_iter)
   if (.not. allocated(conv%history_energy)) return
   call landscape_json_add_series(self, 'energy', conv%history_energy)
   call landscape_json_add_series(self, 'gnorm', conv%history_gnorm)
   call landscape_json_add_series(self, 'gmax', conv%history_gmax)
   call landscape_json_add_series(self, 'dispnorm', conv%history_dispnorm)
   call landscape_json_add_series(self, 'dispmax', conv%history_dispmax)
   call landscape_json_add_series(self, 'trust', conv%history_trust)
   if (allocated(conv%history_method)) then
      call landscape_json_add_series(self, 'method', &
         & real(conv%history_method, wp))
   end if

end subroutine landscape_json_add_history

!> Append the fields every landscape task shares.
!>
!> `status` is the task's own verdict, not a second convergence test: the
!> caller passes it in so that a task that reports a warning (a scan with
!> unconverged points, a frequency analysis with imaginary modes) can
!> distinguish "usable with caveats" from "converged".
subroutine landscape_json_add_result(self, result, status)
   !> Document to extend
   type(landscape_json_document), intent(inout) :: self
   !> Generic task result
   type(landscape_result), intent(in) :: result
   !> Task verdict, one of `landscape_json_status_*`
   integer, intent(in) :: status

   self%status = status
   self%converged = result%converged
   if (allocated(result%message)) then
      self%message = result%message
   else
      self%message = ''
   end if
   call landscape_json_add_counter(self, 'n_eval', result%n_eval)
   call landscape_json_add_counter(self, 'n_cache_miss', result%n_cache_miss)
   call landscape_json_add_counter(self, 'n_steps', result%n_steps)
   call landscape_json_add_value(self, 'energy', result%energy)
   call landscape_json_add_value(self, 'gnorm', result%gnorm)
   call landscape_json_add_value(self, 'gmax', result%gmax)
   call landscape_json_add_value(self, 'etemp', result%etemp)
   call landscape_json_add_counter(self, 'n_imag', result%n_imag)
   call landscape_json_add_value(self, 'imag_freq', result%imag_freq)
   call landscape_json_add_flag(self, 'is_saddle', result%is_saddle)
   call landscape_json_add_history(self, result%conv)

end subroutine landscape_json_add_result

!> Append the uniform core block of the document.
!>
!> This is the one routine every task calls, suite or native, so the core keys
!> (`n_atoms`, `energy`, `geometry`, `provenance`) have the same names, types
!> and shapes in every `landscape_<task>.json` and a single parser can read
!> them all.  `energy` is the canonical total energy in hartree (for the `irc`
!> task the transition state energy); `geometry` is written in bohr with the
!> unit stated inside the block.
!>
!> `provenance` is assembled here from the single sources of truth: the
!> `oxtb_version`/`engine_version` constants, the commit inside `version`, the
!> verbatim command line, the GFN method of the global settings, and the
!> charge/spin of the molecule.  The input file name is passed in because no
!> global carries it; `input_file` is an unallocated allocatable for the API
!> path (`sp` through the test driver), which is written as an empty string.
subroutine landscape_json_add_core(self, mol, energy, input_file)
   !> Document to extend
   type(landscape_json_document), intent(inout) :: self
   !> Molecule the document describes (charge, spin, symbols, coordinates)
   type(TMolecule), intent(in) :: mol
   !> Canonical total energy of the run (hartree)
   real(wp), intent(in) :: energy
   !> Positional input file name, unallocated when the run has none
   character(len=:), allocatable, intent(in) :: input_file

   integer :: i, length

   self%n_atoms = mol%n
   self%energy = energy

   if (allocated(self%geometry_symbol)) deallocate(self%geometry_symbol)
   allocate(self%geometry_symbol(mol%n))
   self%geometry_symbol = ''
   do i = 1, mol%n
      self%geometry_symbol(i) = mol%sym(i)
   end do
   if (allocated(self%geometry_xyz)) deallocate(self%geometry_xyz)
   allocate(self%geometry_xyz(3, mol%n))
   self%geometry_xyz = mol%xyz

   self%provenance%oxtb_version = oxtb_version
   self%provenance%engine_version = engine_version
   self%provenance%commit = commit_from_version(version)
   self%provenance%method = gfn_method_name(set%gfn_method)
   self%provenance%charge = nint(mol%chrg)
   self%provenance%spin = mol%uhf
   if (allocated(input_file)) then
      self%provenance%input_file = input_file
   else
      self%provenance%input_file = ''
   end if
   call get_command(length=length)
   if (length > 0) then
      ! `command_line` was allocated (to the empty string) by
      ! `landscape_json_document_init`; a deferred-length assignment
      ! reallocates it to the command line length, which is what gfortran
      ! supports here (an explicit `allocate` of the already allocated
      ! component aborts at run time).
      self%provenance%command_line = command_line_text(length)
   else
      self%provenance%command_line = ''
   end if

end subroutine landscape_json_add_core

!> Verbatim command line as a freshly allocated string of `length` characters.
!>
!> A helper keeps the `get_command` call (which writes into a fixed length
!> buffer) out of the assignment above: the deferred-length component then
!> only sees a plain assignment.  It is not `pure`: `get_command` reads the
!> process invocation, which is not a pure operation.
function command_line_text(length) result(text)
   !> Length of the command line in characters
   integer, intent(in) :: length
   !> Command line
   character(len=:), allocatable :: text

   allocate(character(len=length) :: text)
   call get_command(text)

end function command_line_text

!> Human readable name of a GFN method number.
!>
!> The spelling is the one `xtbout.json` uses for the numbered methods
!> (`GFN<level>-xTB`, `src/main/json.F90`) plus `GFN-FF` for the force field;
!> it is defined here once so the `method` provenance key has a single
!> spelling.  `p_gfn_method_ff` is the force field number, an unknown number
!> (the settings default `-1`) is reported as `unknown`.
pure function gfn_method_name(gfn_method) result(name)
   !> GFN method number (`set%gfn_method`)
   integer, intent(in) :: gfn_method
   !> Method name
   character(len=:), allocatable :: name

   select case (gfn_method)
   case (0)
      name = 'GFN0-xTB'
   case (1)
      name = 'GFN1-xTB'
   case (2)
      name = 'GFN2-xTB'
   case (p_gfn_method_ff)
      name = 'GFN-FF'
   case default
      name = 'unknown'
   end select

end function gfn_method_name

!> Extract the commit hash out of the `version` string.
!>
!> `version` has the upstream spelling `6.7.1 (<40-hex-commit>)`
!> (`assets/templates/version.f90`); the commit is the text between the
!> parentheses.  A version string without parentheses yields an empty commit
!> instead of an error, so an unusual build cannot break the writer.
pure function commit_from_version(v) result(commit)
   !> Version string
   character(len=*), intent(in) :: v
   !> Commit hash, empty when `v` does not carry one
   character(len=:), allocatable :: commit

   integer :: open, close

   open = index(v, '(')
   close = index(v, ')')
   if (open > 0 .and. close > open) then
      commit = v(open + 1:close - 1)
   else
      commit = ''
   end if

end function commit_from_version

!> Human readable name of a status code.
pure function landscape_json_status_name(status) result(name)
   !> Status code
   integer, intent(in) :: status
   !> Lower case status name
   character(len=:), allocatable :: name

   select case (status)
   case (landscape_json_status_ok)
      name = 'ok'
   case (landscape_json_status_incomplete)
      name = 'incomplete'
   case default
      name = 'error'
   end select

end function landscape_json_status_name

!> Serialise a document to `landscape_<task>.json`.
!>
!> The name goes through the xtb file namespace (`get_namespace`), exactly
!> like every other file of the suite; the name is printed once on `env%unit`
!> as `wrote <file>` (PLAN section 5.5).
subroutine landscape_json_write(env, self, ok)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Document to write
   type(landscape_json_document), intent(in) :: self
   !> True if the file was written
   logical, intent(out) :: ok

   character(len=:), allocatable :: fname

   fname = get_namespace('landscape_'//trim(self%task)//'.json')
   call landscape_json_write_file(env, self, fname, ok)

end subroutine landscape_json_write

!> Write the unified document of a native run type.
!>
!> The native run types (`--sp`, `--opt`, `--hess`, `--ohess`, `--md`, ...)
!> have no `landscape_result`; this routine fills the same core block from the
!> final molecule, the total energy and the SCC convergence flag, adds the
!> gradient norm as the one generic value, and writes the document under the
!> same `landscape_<task>.json` name the suite uses.  It is the native half of
!> the "one parser reads every task" contract of report 117.
!>
!> `task` is the short run type name (`sp`, `opt`, `hess`, `ohess`, `md`, ...),
!> `input_file` is the positional geometry file of the run.  The counters,
!> flags and history of the document are empty for a native run; `values`
!> carries only `gnorm` (hartree/bohr), so the task specific block stays
!> ignorable by a parser that does not know the native path.
subroutine landscape_json_write_run(env, mol, energy, gnorm, converged, task, &
      & input_file, ok)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Final molecule of the run
   type(TMolecule), intent(in) :: mol
   !> Total energy of the run (hartree)
   real(wp), intent(in) :: energy
   !> Norm of the final gradient (hartree/bohr)
   real(wp), intent(in) :: gnorm
   !> SCC convergence flag of the run
   logical, intent(in) :: converged
   !> Short run type name (`sp`, `opt`, `hess`, ...)
   character(len=*), intent(in) :: task
   !> Positional geometry input file of the run
   character(len=:), allocatable, intent(in) :: input_file
   !> True if the document was written
   logical, intent(out) :: ok

   type(landscape_json_document) :: doc

   call landscape_json_document_init(doc, task)
   if (converged) then
      doc%status = landscape_json_status_ok
   else
      doc%status = landscape_json_status_error
   end if
   doc%converged = converged
   call landscape_json_add_core(doc, mol, energy, input_file)
   call landscape_json_add_value(doc, 'gnorm', gnorm)
   call landscape_json_write(env, doc, ok)

end subroutine landscape_json_write_run

!> Serialise a document to an explicit file name (no namespace applied).
!>
!> This is the routine the determinism test drives directly
!> (`test_landscape_e2e`, `ljson-determinism`); the tasks use the namespaced
!> wrapper above.
subroutine landscape_json_write_file(env, self, filename, ok)
   character(len=*), parameter :: source = 'landscape_json_write_file'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Document to write
   type(landscape_json_document), intent(in) :: self
   !> File to write
   character(len=*), intent(in) :: filename
   !> True if the file was written
   logical, intent(out) :: ok

   character(len=:), allocatable :: fname
   integer :: unit, ierr, i

   ok = .false.
   fname = filename
   open(newunit=unit, file=fname, action='write', iostat=ierr)
   if (ierr /= 0) then
      call env%error("cannot write "//fname, source)
      return
   end if

   write(unit, '(a)') '{'
   write(unit, '(2x,a,i0,a)') '"format": ', 2, ','
   write(unit, '(2x,a,a,a)') '"task": "', trim(self%task), '",'
   write(unit, '(2x,a,i0,a)') '"status": ', self%status, ','
   write(unit, '(2x,a,a,a)') '"status_name": "', &
      & trim(landscape_json_status_name(self%status)), '",'
   write(unit, '(2x,a)', advance='no') '"converged": '
   call write_logical(unit, self%converged)
   write(unit, '(a)') ','
   write(unit, '(2x,a,a,a)') '"message": "', &
      & json_escape(trim(self%message)), '",'
   write(unit, '(2x,a,i0,a)') '"n_atoms": ', self%n_atoms, ','
   write(unit, '(2x,a)', advance='no') '"energy": '
   call write_real(unit, self%energy)
   write(unit, '(a)') ','
   call write_geometry(unit, self)
   call write_provenance(unit, self)

   write(unit, '(2x,a)') '"counters": {'
   do i = 1, self%n_counter
      if (i < self%n_counter) then
         write(unit, '(4x,a,a,a,i0,a)') '"', trim(self%counter_key(i)), '": ', &
            & self%counter_value(i), ','
      else
         write(unit, '(4x,a,a,a,i0)') '"', trim(self%counter_key(i)), '": ', &
            & self%counter_value(i)
      end if
   end do
   write(unit, '(2x,a)') '},'

   write(unit, '(2x,a)') '"values": {'
   do i = 1, self%n_value
      if (i < self%n_value) then
         write(unit, '(4x,a,a,a)', advance='no') '"', trim(self%value_key(i)), '": '
         call write_real(unit, self%value_value(i))
         write(unit, '(a)') ','
      else
         write(unit, '(4x,a,a,a)', advance='no') '"', trim(self%value_key(i)), '": '
         call write_real(unit, self%value_value(i))
         write(unit, *)
      end if
   end do
   write(unit, '(2x,a)') '},'

   write(unit, '(2x,a)') '"flags": {'
   do i = 1, self%n_flag
      if (i < self%n_flag) then
         write(unit, '(4x,a,a,a)', advance='no') '"', trim(self%flag_key(i)), '": '
         call write_logical(unit, self%flag_value(i))
         write(unit, '(a)') ','
      else
         write(unit, '(4x,a,a,a)', advance='no') '"', trim(self%flag_key(i)), '": '
         call write_logical(unit, self%flag_value(i))
         write(unit, *)
      end if
   end do
   write(unit, '(2x,a)') '},'

   write(unit, '(2x,a)') '"history": {'
   do i = 1, self%n_series
      write(unit, '(4x,a,a,a)', advance='no') '"', trim(self%series_key(i)), '": ['
      call write_real_range(unit, self%series_value, self%series_offset(i), &
         & self%series_length(i))
      if (i < self%n_series) then
         write(unit, '(a)') '],'
      else
         write(unit, '(a)') ']'
      end if
   end do
   write(unit, '(2x,a)') '}'

   write(unit, '(a)') '}'
   close(unit)

   write(env%unit, '(1x,a,a)') "wrote ", fname
   ok = .true.

end subroutine landscape_json_write_file

!> Serialise the geometry block.
!>
!> The shape is identical for every task: an object with the unit stated
!> explicitly (`bohr`, the internal unit of the suite), a `symbols` array and
!> a `coordinates` array of `[x, y, z]` triples, so a parser never has to
!> branch on the task.  A document built without `landscape_json_add_core`
!> has `n_atoms = 0` and writes two empty arrays.
subroutine write_geometry(unit, self)
   !> Output unit
   integer, intent(in) :: unit
   !> Document to serialise
   type(landscape_json_document), intent(in) :: self

   integer :: i

   write(unit, '(2x,a)') '"geometry": {'
   write(unit, '(4x,a)') '"unit": "bohr",'
   write(unit, '(4x,a)', advance='no') '"symbols": ['
   do i = 1, self%n_atoms
      if (i > 1) write(unit, '(a)', advance='no') ', '
      write(unit, '(a,a,a)', advance='no') '"', &
         & trim(self%geometry_symbol(i)), '"'
   end do
   write(unit, '(a)') '],'
   write(unit, '(4x,a)') '"coordinates": ['
   do i = 1, self%n_atoms
      write(unit, '(6x,a)', advance='no') '['
      call write_real(unit, self%geometry_xyz(1, i))
      write(unit, '(a)', advance='no') ', '
      call write_real(unit, self%geometry_xyz(2, i))
      write(unit, '(a)', advance='no') ', '
      call write_real(unit, self%geometry_xyz(3, i))
      if (i < self%n_atoms) then
         write(unit, '(a)') '],'
      else
         write(unit, '(a)') ']'
      end if
   end do
   write(unit, '(4x,a)') ']'
   write(unit, '(2x,a)') '},'

end subroutine write_geometry

!> Serialise the provenance block.
!>
!> All eight keys are always written (the strings may be empty), so the block
!> has the same shape in every task document.  The values come from
!> `landscape_json_add_core`; a document built without it writes empty strings
!> and zero charge/spin.
subroutine write_provenance(unit, self)
   !> Output unit
   integer, intent(in) :: unit
   !> Document to serialise
   type(landscape_json_document), intent(in) :: self

   write(unit, '(2x,a)') '"provenance": {'
   write(unit, '(4x,a,a,a)') '"oxtb_version": "', &
      & json_escape(self%provenance%oxtb_version), '",'
   write(unit, '(4x,a,a,a)') '"engine_version": "', &
      & json_escape(self%provenance%engine_version), '",'
   write(unit, '(4x,a,a,a)') '"commit": "', &
      & json_escape(self%provenance%commit), '",'
   write(unit, '(4x,a,a,a)') '"command_line": "', &
      & json_escape(self%provenance%command_line), '",'
   write(unit, '(4x,a,a,a)') '"method": "', &
      & json_escape(self%provenance%method), '",'
   write(unit, '(4x,a,i0,a)') '"charge": ', self%provenance%charge, ','
   write(unit, '(4x,a,i0,a)') '"spin": ', self%provenance%spin, ','
   write(unit, '(4x,a,a,a)') '"input_file": "', &
      & json_escape(self%provenance%input_file), '"'
   write(unit, '(2x,a)') '},'

end subroutine write_provenance

!> Write a logical in JSON syntax (`true` / `false`).
!>
!> The Fortran `L` edit descriptor prints `T`/`F`, which JSON does not accept;
!> the conversion is explicit so a reader can parse the file unmodified.
subroutine write_logical(unit, value)
   !> Output unit
   integer, intent(in) :: unit
   !> Value to write
   logical, intent(in) :: value

   if (value) then
      write(unit, '(a)', advance='no') 'true'
   else
      write(unit, '(a)', advance='no') 'false'
   end if

end subroutine write_logical

!> Write one real number in JSON syntax; a non finite value becomes `null`.
subroutine write_real(unit, value)
   !> Output unit
   integer, intent(in) :: unit
   !> Value to write
   real(wp), intent(in) :: value

   character(len=32) :: buffer

   if (.not. ieee_is_finite(value)) then
      ! NaN and infinity have no JSON literal; `null` is an honest placeholder
      ! that cannot be mistaken for a measured number.  `ieee_is_finite` is
      ! used because the obvious `value /= value` comparison would itself raise
      ! the invalid exception under a checked build.
      write(unit, '(a)', advance='no') 'null'
      return
   end if
   write(buffer, '(es25.17e3)') value
   write(unit, '(a)', advance='no') trim(adjustl(buffer))

end subroutine write_real

!> Write a contiguous range of the series pool in JSON syntax (one line).
subroutine write_real_range(unit, pool, offset, length)
   !> Output unit
   integer, intent(in) :: unit
   !> Series pool
   real(wp), intent(in) :: pool(:)
   !> Start index of the range in `pool`
   integer, intent(in) :: offset
   !> Number of values of the range (may be zero)
   integer, intent(in) :: length

   integer :: i

   do i = 0, length - 1
      if (i > 0) write(unit, '(a)', advance='no') ', '
      call write_real(unit, pool(offset + i))
   end do

end subroutine write_real_range

!> Escape the characters JSON cannot carry inside a string.
pure function json_escape(text) result(escaped)
   !> Input text
   character(len=*), intent(in) :: text
   !> Escaped text
   character(len=:), allocatable :: escaped

   character(len=1) :: c
   integer :: i

   escaped = ''
   do i = 1, len_trim(text)
      c = text(i:i)
      select case (c)
      case ('"')
         escaped = escaped//'\"'
      case ('\')
         escaped = escaped//'\\'
      case (achar(10))
         escaped = escaped//'\n'
      case (achar(13))
      case (achar(9))
         escaped = escaped//'\t'
      case default
         escaped = escaped//c
      end select
   end do

end function json_escape

end module xtb_landscape_json
