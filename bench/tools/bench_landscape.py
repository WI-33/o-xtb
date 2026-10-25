#!/usr/bin/env python3
"""Wall-time benchmark of the six o-xtb PES tasks on in-repository fixtures.

The suite has grown tests that certify *numbers*, but nothing that measures
*time*: every performance claim in this tree was made by hand, on one machine,
with a command line written down in a report.  This tool is that command line.
It runs one fixed set of tasks against fixed input files and prints a table of
wall-clock times; the table is the baseline a later patch is compared against.

What it is not: it prints numbers, not verdicts.  Timings depend on the
machine, the load, the compiler and the BLAS, so **no gate uses this tool** and
it has no ``--check`` mode.  Compare two runs of *the same* machine and the
*same* binary against each other, never a run against a stored number.

The machine model is pinned the way the suite's own measurements are (AGENTS.md
section 6.12): ``OMP_NUM_THREADS=1`` and ``OPENBLAS_NUM_THREADS=1`` are forced
into the child environment, one process per timed run, a **fresh rundir per
repetition** and ``--norestart`` (a ``landscape3.xtbrestart`` left behind by an
earlier repetition changes both the iteration count and the energy).  ``XTBPATH``
is inherited, not set: the build tree finds its own parameter directory.

Usage::

    bench_landscape.py <oxtb> [--runs N] [--only NAME[,NAME]] [--list]
                       [--work-dir DIR] [--keep]

The default repetition count of a case is 3 for the fast tasks and 1 for the two
that cost of the order of 10 s per run (``lopt-rfo`` converges the 33-atom
system and ``hess-taxol`` is the 113-atom finite difference Hessian);
``--runs`` overrides every case, ``--only`` restricts the selection, ``--list``
prints the case table without running anything.
"""

import argparse
import os
import platform
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

#: Repository root: the grandparent of ``bench/tools/``.
ROOT = Path(__file__).resolve().parents[2]

#: One entry per benchmarked task.  ``input`` is the positional geometry of the
#: run, ``product`` the ``--lts-product`` endpoint where the task needs one, and
#: ``runs`` the default number of repetitions (median over them).
CASES = [
    {
        "name": "sp-gfn2",
        "atoms": 30,
        "input": "examples/sp/sp-01-gfn2-single-point-json/input.xyz",
        "args": ["--gfn", "2", "--sp"],
        "runs": 3,
    },
    {
        "name": "lopt-lbfgs",
        "atoms": 33,
        "input": "examples/lopt/lopt-01-lbfgs-from-distorted/input.xyz",
        "args": ["--gfn", "2", "--lopt", "--lopt-engine", "lbfgs",
                 "--lopt-level", "normal"],
        "runs": 3,
    },
    {
        "name": "lopt-rfo",
        "atoms": 33,
        "input": "examples/lopt/lopt-01-lbfgs-from-distorted/input.xyz",
        "args": ["--gfn", "2", "--lopt", "--lopt-engine", "rfo",
                 "--lopt-level", "normal"],
        "runs": 1,
    },
    {
        "name": "lvib",
        "atoms": 30,
        "input": "examples/hess/hess-01-numerical/input.xyz",
        "args": ["--gfn", "2", "--lvib"],
        "runs": 3,
    },
    {
        "name": "lts-band-neb",
        "atoms": 30,
        "input": "examples/structures/benzoic_acid_dimer/reactant.xyz",
        "product": "examples/structures/benzoic_acid_dimer/ts_guess.xyz",
        "args": ["--gfn", "2", "--lts", "--lts-method", "neb",
                 "--lts-images", "5", "--lts-max-steps", "5"],
        "runs": 3,
    },
    {
        "name": "lirc-hpc",
        "atoms": 30,
        "input": "examples/irc/irc-01-four-methods/input.xyz",
        "args": ["--gfn", "2", "--lirc", "--lirc-method", "hpc",
                 "--lirc-max-steps", "40"],
        "runs": 3,
    },
    {
        "name": "hess-taxol",
        "atoms": 113,
        "input": "assets/inputs/xyz/taxol.xyz",
        "args": ["--gfn", "2", "--hess"],
        "runs": 1,
    },
]

#: Forced into every child process: one thread, one BLAS thread, so the numbers
#: are comparable between runs of one binary (AGENTS.md section 6.12).
THREAD_ENV = {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}


def fail(message):
    sys.stderr.write("bench_landscape.py: error: %s\n" % message)
    raise SystemExit(2)


def case_arguments(case):
    """The complete argument list of one case, fixture paths resolved."""
    argv = [str(ROOT / case["input"])]
    argv += case["args"]
    if "product" in case:
        argv += ["--lts-product", str(ROOT / case["product"])]
    argv += ["--norestart"]
    return argv


def reported_count(output, pattern_prefix):
    """The last ``<prefix> <n>`` counter of the run's stdout, or ``-1``."""
    value = -1
    for line in output.splitlines():
        if pattern_prefix in line:
            token = line.split(pattern_prefix)[-1].split()[-1]
            try:
                value = int(token)
            except ValueError:
                pass
    return value


def time_once(binary, argv, rundir):
    """Run one case once and return (wall seconds, engine evaluations, exit)."""
    if rundir.exists():
        shutil.rmtree(rundir)
    rundir.mkdir(parents=True)
    env = dict(os.environ)
    env.update(THREAD_ENV)
    start = time.perf_counter()
    proc = subprocess.run(
        [str(binary), *argv],
        cwd=rundir,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=7200,
    )
    wall = time.perf_counter() - start
    evaluations = reported_count(proc.stdout, "landscape engine: energy evaluations")
    return wall, evaluations, proc.returncode


def binary_header(binary):
    """The version block of the binary, as two ``key: value`` lines."""
    proc = subprocess.run([str(binary), "--version"], stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True, timeout=300)
    lines = [line.strip() for line in proc.stdout.splitlines() if line.strip()]
    version = next((line for line in lines if "version" in line.lower()), "?")
    compiled = next((line for line in lines if "ompiled" in line), "?")
    return version, compiled


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Wall-time benchmark of the six o-xtb PES tasks "
                    "(prints numbers, not verdicts).")
    parser.add_argument("binary", type=Path, nargs="?",
                        help="the oxtb executable to time (required unless "
                             "--list is given)")
    parser.add_argument("--runs", type=int, default=0,
                        help="repetitions of every case (default: per case, 3 "
                             "for the fast tasks and 1 for the two expensive ones)")
    parser.add_argument("--only", default="",
                        help="comma separated case names to run (see --list)")
    parser.add_argument("--list", action="store_true",
                        help="print the case table and exit")
    parser.add_argument("--work-dir", type=Path, default=None,
                        help="where the per-run directories are created "
                             "(default: a fresh temporary directory)")
    parser.add_argument("--keep", action="store_true",
                        help="keep the per-run directories (for a post mortem)")
    args = parser.parse_args(argv)

    if args.list:
        for case in CASES:
            print("%-14s %3d atoms  runs %d  %s"
                  % (case["name"], case["atoms"], case["runs"],
                     " ".join(case["args"])))
        return 0

    binary = None if args.binary is None else args.binary.resolve()
    if binary is None:
        fail("no executable given (unless --list is used, one is required)")
    if not binary.is_file():
        fail("no executable at %s" % binary)

    selected = CASES
    if args.only:
        wanted = [name.strip() for name in args.only.split(",") if name.strip()]
        known = {case["name"] for case in CASES}
        unknown = [name for name in wanted if name not in known]
        if unknown:
            fail("unknown case(s) %s; choose from %s"
                 % (", ".join(unknown), ", ".join(sorted(known))))
        selected = [case for case in CASES if case["name"] in wanted]

    if args.work_dir is not None:
        work = args.work_dir.resolve()
        work.mkdir(parents=True, exist_ok=True)
        temporary = False
    else:
        work = Path(tempfile.mkdtemp(prefix="oxtb-bench-"))
        temporary = True

    thunks = []  # (case, wall, evaluations, exit)
    try:
        print("oxtb benchmark")
        version, compiled = binary_header(binary)
        print("  binary   %s" % binary)
        print("  version  %s" % version)
        print("  built    %s" % compiled)
        print("  machine  %s, %s, %d cpu" % (platform.system(),
                                              platform.machine(),
                                              os.cpu_count() or 0))
        print("  threads  OMP_NUM_THREADS=%s OPENBLAS_NUM_THREADS=%s"
              % (THREAD_ENV["OMP_NUM_THREADS"],
                 THREAD_ENV["OPENBLAS_NUM_THREADS"]))
        print()
        for case in selected:
            runs = args.runs if args.runs > 0 else case["runs"]
            walls = []
            evaluations = -1
            status = 0
            for repetition in range(runs):
                rundir = work / ("%s-%d" % (case["name"], repetition + 1))
                wall, evaluations, status = time_once(
                    binary, case_arguments(case), rundir)
                walls.append(wall)
                print("  %-14s run %d/%d  %8.3f s  exit %d"
                      % (case["name"], repetition + 1, runs, wall, status))
                sys.stdout.flush()
            thunks.append((case, walls, evaluations, status))
    finally:
        if temporary and not args.keep:
            shutil.rmtree(work, ignore_errors=True)
        elif temporary:
            print("\nwork directories kept in %s" % work)

    print()
    print("%-14s %5s %3s %11s %11s %11s %9s %4s"
          % ("case", "atoms", "n", "median/s", "min/s", "max/s", "evals", "exit"))
    print("-" * 78)
    for case, walls, evaluations, status in thunks:
        count = str(evaluations) if evaluations >= 0 else "-"
        print("%-14s %5d %3d %11.3f %11.3f %11.3f %9s %4d"
              % (case["name"], case["atoms"], len(walls),
                 statistics.median(walls), min(walls), max(walls),
                 count, status))
    total = sum(statistics.median(walls) for _, walls, _, _ in thunks)
    print("-" * 78)
    print("%-14s %5s %3s %11.3f" % ("total (medians)", "", "", total))
    print()
    print("Times are wall clock of one process at a pinned single thread; they "
          "are comparable between runs of this machine only, and a non-zero "
          "exit status of the band case is its step budget, not a failure.  "
          "The 'evals' column is the landscape engine's energy-evaluation "
          "counter, which the native --sp/--hess routes do not print (shown "
          "as '-').")
    return 0


if __name__ == "__main__":
    sys.exit(main())
