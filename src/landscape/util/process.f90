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

!> Worker processes of the concurrent finite-difference Hessian.
!>
!> The numerical Hessian provider (`xtb_landscape_hessian`) splits its 6N
!> displaced single points over `jobs` process workers when `jobs > 1`.  Each
!> worker is a `fork()` of the running program, so it starts from the parent's
!> calculator, geometry and reference wavefunction instead of rebuilding them
!> (a re-`exec` of the binary would pay the whole process start-up once per
!> point, which is several times the cost of a displaced single point at the
!> 30-atom scale); it writes its two Cartesian gradient columns into a pipe and
!> leaves with `_exit()`, so no atexit handler, no inherited output buffer and
!> no library shutdown path of the parent is executed in the worker.  The
!> parent reads the columns and performs the arithmetic of the differentiation
!> itself, in one process, which is what keeps the resulting Hessian
!> independent of the worker count.
!>
!> `jobs = 1` never enters this module: the provider keeps its serial loop, so
!> the default path has no dependency on process creation at all.
!>
!> The C implementations live in `util/process.c`; on a platform without
!> `fork()`/`pipe()` (Windows) `landscape_process_available` returns false and
!> the provider refuses `jobs > 1` with `env%error` rather than falling back to
!> a serial run the user did not ask for.
!>
!> Units: buffers are transferred byte-exactly in their own kinds, no
!> conversion happens here.  The payloads are `integer(c_int)` status words and
!> `real(c_double)` gradient arrays, and the provider relies on
!> `xtb_mctc_accuracy:wp` being `c_double` - a different working precision
!> would be rejected at compile time by the explicit interface rather than
!> silently truncated.
module xtb_landscape_process
   use, intrinsic :: iso_c_binding, only : c_int, c_long, c_double
   implicit none(type, external)
   private

   public :: landscape_process_available
   public :: landscape_process_pipe
   public :: landscape_process_fork
   public :: landscape_process_close
   public :: landscape_process_read_int, landscape_process_write_int
   public :: landscape_process_read_real, landscape_process_write_real
   public :: landscape_process_wait
   public :: landscape_process_status_text
   public :: landscape_process_worker_exit

   interface
      !> True if this build can create workers (POSIX `fork`/`pipe`)
      function landscape_process_available() result(available) bind(c, &
         & name='xtb_landscape_process_available')
         import :: c_int
         !> 1 if process creation is available, 0 otherwise
         integer(c_int) :: available
      end function landscape_process_available

      !> Create the result pipe of one worker
      function landscape_process_pipe(fds) result(status) bind(c, &
         & name='xtb_landscape_pipe')
         import :: c_int
         !> `fds(1)` is the read end, `fds(2)` the write end
         integer(c_int), intent(out) :: fds(2)
         !> 0 on success, -1 on error
         integer(c_int) :: status
      end function landscape_process_pipe

      !> Fork one worker
      function landscape_process_fork() result(pid) bind(c, &
         & name='xtb_landscape_fork')
         import :: c_int
         !> -1 on error, 0 in the worker, the worker's pid in the parent
         integer(c_int) :: pid
      end function landscape_process_fork

      !> Close one file descriptor
      function landscape_process_close(fd) result(status) bind(c, &
         & name='xtb_landscape_close')
         import :: c_int
         !> File descriptor
         integer(c_int), value :: fd
         !> 0 on success, -1 on error
         integer(c_int) :: status
      end function landscape_process_close

      !> Read `len` bytes, or as many as the writer produced
      function landscape_process_read_int(fd, buffer, len) result(got) bind(c, &
         & name='xtb_landscape_read')
         import :: c_int, c_long
         !> File descriptor
         integer(c_int), value :: fd
         !> Destination buffer
         integer(c_int), intent(out) :: buffer(*)
         !> Bytes requested
         integer(c_long), value :: len
         !> Bytes received, -1 on error
         integer(c_long) :: got
      end function landscape_process_read_int

      !> Write `len` bytes of integer data
      function landscape_process_write_int(fd, buffer, len) result(put) bind(c, &
         & name='xtb_landscape_write')
         import :: c_int, c_long
         !> File descriptor
         integer(c_int), value :: fd
         !> Source buffer
         integer(c_int), intent(in) :: buffer(*)
         !> Bytes to write
         integer(c_long), value :: len
         !> Bytes written, -1 on error
         integer(c_long) :: put
      end function landscape_process_write_int

      !> Read `len` bytes, or as many as the writer produced
      function landscape_process_read_real(fd, buffer, len) result(got) bind(c, &
         & name='xtb_landscape_read')
         import :: c_int, c_long, c_double
         !> File descriptor
         integer(c_int), value :: fd
         !> Destination buffer
         real(c_double), intent(out) :: buffer(*)
         !> Bytes requested
         integer(c_long), value :: len
         !> Bytes received, -1 on error
         integer(c_long) :: got
      end function landscape_process_read_real

      !> Write `len` bytes of real data
      function landscape_process_write_real(fd, buffer, len) result(put) bind(c, &
         & name='xtb_landscape_write')
         import :: c_int, c_long, c_double
         !> File descriptor
         integer(c_int), value :: fd
         !> Source buffer
         real(c_double), intent(in) :: buffer(*)
         !> Bytes to write
         integer(c_long), value :: len
         !> Bytes written, -1 on error
         integer(c_long) :: put
      end function landscape_process_write_real

      !> Wait for one worker
      function landscape_process_wait(pid, status) result(rc) bind(c, &
         & name='xtb_landscape_wait')
         import :: c_int
         !> Worker pid returned by `landscape_process_fork`
         integer(c_int), value :: pid
         !> Raw wait status; 0 means "exited normally with status 0"
         integer(c_int), intent(out) :: status
         !> 0 on success, -1 on error
         integer(c_int) :: rc
      end function landscape_process_wait

      !> Leave the current worker without running the parent's atexit handlers
      subroutine landscape_process_worker_exit(status) bind(c, &
         & name='xtb_landscape_worker_exit')
         import :: c_int
         !> Exit status of the worker
         integer(c_int), value :: status
      end subroutine landscape_process_worker_exit
   end interface

contains

!> Describe a raw worker wait status in words.
!>
!> `landscape_process_wait` stores the status exactly as `waitpid(2)` returns
!> it, i.e. as a packed word: the low seven bits are the terminating signal and
!> the next byte the exit status of a normal exit (`WIFSIGNALED`,
!> `WIFEXITED`/`WEXITSTATUS`).  Printed raw, a worker that left with `exit(1)`
!> claims "status 256" and one killed by `SIGSEGV` claims "status 11" - numbers
!> that read like corruption.  The two cases are decoded here so that a report
!> can name what happened; the decoding is POSIX specific, which is where
!> workers exist at all (elsewhere `landscape_process_available` is 0 and no
!> status is ever produced).
function landscape_process_status_text(status) result(text)
   !> Raw wait status of a reaped worker
   integer(c_int), intent(in) :: status
   !> One line, e.g. "exit status 1" or "terminated by signal 11"
   character(len=:), allocatable :: text

   character(len=64) :: buffer
   integer :: raw

   raw = int(status)
   if (iand(raw, 127) == 0) then
      write(buffer, '(a,i0)') "exit status ", ishft(raw, -8)
   else
      write(buffer, '(a,i0)') "terminated by signal ", iand(raw, 127)
   end if
   text = trim(buffer)

end function landscape_process_status_text

end module xtb_landscape_process
