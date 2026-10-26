/*
 * This file is part of xtb.
 *
 * Copyright (C) 2026 xtb contributors
 *
 * xtb is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * xtb is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with xtb.  If not, see <https://www.gnu.org/licenses/>.
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 *
 * Process primitives of the concurrent finite-difference Hessian.
 *
 * `landscape_numerical_hessian_provider` (`src/landscape/types/hessian.f90`)
 * evaluates the 6N displaced single points of a Hessian in `jobs` concurrent
 * worker *processes*.  Processes and not threads, because the calculator, the
 * global `set`, the error log and the library state of xtb are process-global
 * mutable state that a thread would share; a worker is created with `fork()`
 * from the running program, so it inherits the calculator, the geometry and
 * the reference wavefunction instead of having to rebuild them, and it writes
 * the two gradients of its Cartesian column into a pipe before it exits.
 *
 * The window between `fork()` and `_exit()` is short and contains nothing but
 * the worker's own single points, which is what makes the inherited state safe
 * to use: every displaced point starts from the *same* copied wavefunction
 * (the provider resets `chk` before each point), so no worker can observe a
 * state that another worker produced, and the parent's copies are never
 * written.  The worker leaves with `_exit()` and not with `exit()`: `exit()`
 * would run the atexit handlers and flush the *inherited* standard units,
 * which would duplicate the parent's buffered output and reach for library
 * state that does not belong to this process (the static-exit defect of
 * `docs/_simplify-recon/96-static-exit-sigsegv.md` is exactly a weak
 * `pthread` reference of that exit path).
 *
 * `fork()` and `pipe()` do not exist on Windows, so this file is the
 * portability boundary: on such a platform `xtb_landscape_process_available`
 * reports 0 and the provider refuses `jobs > 1` with `env%error` instead of
 * silently falling back.  `jobs = 1` never calls into this file, so the
 * default single-process path is complete without it.
 */

#if defined(_WIN32)

int
xtb_landscape_process_available(void)
{
   return 0;
}

int
xtb_landscape_pipe(int *fds)
{
   (void) fds;
   return -1;
}

int
xtb_landscape_fork(void)
{
   return -1;
}

int
xtb_landscape_close(int fd)
{
   (void) fd;
   return -1;
}

long
xtb_landscape_read(int fd, void *buf, long len)
{
   (void) fd;
   (void) buf;
   (void) len;
   return -1;
}

long
xtb_landscape_write(int fd, const void *buf, long len)
{
   (void) fd;
   (void) buf;
   (void) len;
   return -1;
}

int
xtb_landscape_wait(int pid, int *status)
{
   (void) pid;
   (void) status;
   return -1;
}

void
xtb_landscape_worker_exit(int status)
{
   (void) status;
}

#else

#include <errno.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

int
xtb_landscape_process_available(void)
{
   return 1;
}

/* Create the result pipe of one worker; `fds[0]` is read, `fds[1]` is write. */
int
xtb_landscape_pipe(int *fds)
{
   return pipe(fds);
}

/* One `fork()`: -1 on failure, 0 in the worker, the worker's pid otherwise. */
int
xtb_landscape_fork(void)
{
   return (int) fork();
}

int
xtb_landscape_close(int fd)
{
   return close(fd);
}

/*
 * Read until `len` bytes arrived, end of file or error; the return value is
 * the number of bytes actually transferred, or -1.  A short count means that
 * the worker died before it finished writing.
 */
long
xtb_landscape_read(int fd, void *buf, long len)
{
   char *pos = (char *) buf;
   long done = 0;
   while (done < len) {
      ssize_t got = read(fd, pos + done, (size_t) (len - done));
      if (got < 0) {
         if (errno == EINTR)
            continue;
         return -1;
      }
      if (got == 0)
         break;
      done += (long) got;
   }
   return done;
}

/* Write all of `len` bytes; the return value is `len`, or -1 on error. */
long
xtb_landscape_write(int fd, const void *buf, long len)
{
   const char *pos = (const char *) buf;
   long done = 0;
   while (done < len) {
      ssize_t put = write(fd, pos + done, (size_t) (len - done));
      if (put < 0) {
         if (errno == EINTR)
            continue;
         return -1;
      }
      done += (long) put;
   }
   return done;
}

/* Wait for `pid` and store its raw wait status; 0 on success, -1 on error. */
int
xtb_landscape_wait(int pid, int *status)
{
   int raw = 0;
   pid_t got;
   do {
      got = waitpid((pid_t) pid, &raw, 0);
   } while (got < 0 && errno == EINTR);
   if (got < 0)
      return -1;
   *status = raw;
   return 0;
}

/* Leave the worker without running the inherited atexit handlers. */
void
xtb_landscape_worker_exit(int status)
{
   _exit(status);
}

#endif
