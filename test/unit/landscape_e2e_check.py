#!/usr/bin/env python3
#
# This file is part of xtb.
#
# Copyright (C) 2026 xtb contributors
#
# xtb is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# xtb is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with xtb.  If not, see <https://www.gnu.org/licenses/>.
"""End-to-end validation of the landscape suite (M7, PLAN section 9.8).

This script is the reproducible half of the M7 deliverable: it drives the
three test layers of PLAN section 8.1 through the command line and asserts the
cross-layer contracts that a unit test cannot see.

Layers
------
  L1  analytic potentials   the pure mathematics of the optimisers, saddle
                            searches, paths and IRC (the `landscape-*` unit
                            suites that build no SCF)
  L2  real GFN2 engine      the engine boundary: determinism, cost accounting
                            (`xtb/landscape-optimize-xtb`, `xtb/landscape-vib`)
  L3  end-to-end chain      `opt -> hess -> freq -> thermo` and `ts -> irc`

Modes (the first four mirror a CMake test of `test/unit/CMakeLists.txt`)
----------------------------------------------------------------------
  chain         <exe> <input> <rundir>
        Runs `--lopt` on `input`, then `--lvib` and `--hess` on the optimised
        `xtbopt.xyz`, and compares the frequencies, the zero point energy and
        the thermodynamic functions H(T)/G(T).  The comparison is **not**
        bitwise: the two routes project the same raw numerical Hessian with
        different conventions (Eckart on the mass weighted matrix for
        `--lvib`, the legacy `trproj` on the unweighted one for `--hess`).
        Measured on the converged water minimum this check drives (tight,
        |g| = 7.5e-6 Eh/bohr; report 49 section 1.1, re-measured in report
        51): max frequency difference **0.0100 cm^-1** and max |dZPE| = |dH|
        = |dG| = **1.43e-8 hartree**.  `FREQ_TOL` / `THERMO_TOL` below are
        20x and 70x those measurements; their comments state why that margin
        is needed instead of a tighter bound.
  chain-failure <exe> <input> <rundir>
        The failure path: step 1 (`--lopt`) succeeds, step 2 (`--lvib
        --lvib-hessian file` pointing at a Hessian file that does not exist) is
        misconfigured, so the runner must abort and must not write the
        artefacts of the failed step.  This is the check a "silently continue"
        chain would fail.
  ljson         <exe> <input> <rundir> <task>
        Runs one landscape task with `--ljson` and (a) parses the document
        with the JSON parser of the standard library, (b) asserts that every
        counter it reads is the number the task printed on stdout, (c) runs
        the same job a second time in a fresh directory and requires the two
        files to be byte identical, and (d) checks that `--namespace` prefixes
        the file name.  `task` is one of opt, scan, vib, ts, irc.
  ljson-methods <exe> <input> <rundir> <product>
        Runs all eight `ts` methods and all four `irc` methods with `--ljson`
        and asserts that each writes a parseable document whose `method`
        counter is the identifier of the method (the four band searches take
        a different code path in the task and return early).
  ljson-core <exe> <input> <rundir> <product>
        The single-parser contract of report 117: runs the six suite tasks,
        all eight `ts` methods, all four `irc` methods and the native run
        types (`--sp`, `--grad`, `--opt`, `--hess`, `--ohess`, `--md`) and
        feeds every `landscape_<task>.json` to one parser (`parse_document`),
        which asserts the uniform core block (keys, types, geometry shape,
        unit, `format: 2`) for all of them.  Removing the core block from any
        one task makes it red.
  scan-counters <exe> <input> <rundir>
        The point bookkeeping of a scan: runs a converged scan and a scan
        whose inner relaxation never converges (`--lscan-max-steps 1`) and
        asserts the document's `n_points`/`n_valid`/`n_fail`/`n_unconverged`
        against the frames and omitted points the task prints, plus the
        documented "no converged point => status error, non-zero exit".
  ts-counters <exe> <input> <rundir>
        The sampling point of the `ts` evaluation counter (PLAN section
        12.11).  Runs the same transition state search with and without the
        frequency confirmation (`--lts-no-verify`) and asserts that (a) the
        counter, the printed `cache misses` and the engine's `energy
        evaluations` of *one* run are the same number and (b) the confirmed
        run is more expensive by exactly whole 6N finite difference Hessians
        (one 6N confirmation Hessian for N atoms).  A counter sampled before
        the confirmation makes (a) and (b) fail.
  ts-failure <exe> <input> <rundir>
        A transition state search that runs out of steps (report 221, findings
        S2 and H1): neither the single geometry P-RFO nor the 9-image band may
        pay for the confirmation Hessian, and neither may lose its own failure
        message to the confirmation's verdict.  Asserts, for both, that the
        run exits non-zero, keeps its own message ("maximum number of steps
        reached" / "band reached the iteration limit"), reports `n_imag` = -1
        and `status` incomplete (not error, which is the "converged on
        something that is not a saddle" verdict), never prints the saddle
        confirmation, and costs exactly the same number of evaluations with and
        without `--lts-no-verify` - the search-only cost, pinned per fixture.
  ts-status <exe> <input> <rundir>
        The two failure statuses of the `ts` document, side by side (report
        224, finding 1): a search that ran out of steps reports `status` 1
        (`incomplete`) with its own message and `n_imag` -1, while a search
        that converged onto a *minimum* reports `status` -1 (`error`) with the
        confirmation's message and `n_imag` 0.  The second branch is the H2
        verdict, whose message the task raises through `env`; the status
        predicate must not be written on `converged`, which the confirmation
        clears, or the document is downgraded to `incomplete`.
  layers        <xtb_exe> <tester_exe> <input> <rundir> [<ts_input>]
        Runs L1 (every landscape unit suite that needs no SCF), L2 (the engine
        boundary suites) and L3 (the `opt -> hess -> freq -> thermo` chain,
        plus the `ts -> irc` chain when a transition state input is given) and
        writes a JSON summary to `rundir/layers.json`.
  maple         <xtb_exe> <input> <rundir>
        Records the MAPLE comparison methodology and the topology facts that
        can be compared.  MAPLE's fixtures use machine learning potentials
        (PLAN R-14), so they are **not** a GFN2 numerical baseline; the mode
        therefore compares topology only (number of minima and saddle points,
        endpoint assignment) and states explicitly whether MAPLE was executed.
        Set `MAPLE_DIR` to the MAPLE checkout; the mode reports "not executed"
        when the checkout or its dependencies are unavailable instead of
        pretending it ran.

The environment (XTBPATH, OMP_NUM_THREADS) is inherited from the test runner,
exactly like `landscape_cli_check.py`.
"""

import json
import os
import re
import shutil
import subprocess
import sys

# --- tolerances of the `--lvib` vs `--hess` comparison ---------------------
# Frequency (cm^-1): 20x the **measured** 0.0100 cm^-1 of the water minimum
# this check drives (tight level, |g| = 7.5e-6 Eh/bohr; report 49 section 1.1,
# re-measured in report 51).  The factor is not decoration: the compared
# quantity is the residual of two different rigid body projectors on the same
# raw finite difference matrix, and it grows with |g| - measured on the same
# molecule 0.02 bohr off the minimum it already reaches 3.3-5.2 cm^-1, and a
# *default* level `--lopt` minimum (|g| ~ 4e-4) gives 0.89 cm^-1 for methanol
# (report 49 section 1.2).  The margin therefore covers "the same check on a
# less converged minimum", not solver noise of this input; a caller that
# wants a tight assertion passes the per-molecule `freq_tol` argument.
FREQ_TOL = 0.2
# ZPE/H/G (hartree): 70x the **measured** 1.43e-8 hartree of the same water
# minimum.  It is a floor; the effective tolerance of a run is the larger of
# it and the worst case propagation of the frequency tolerance over the
# internal modes (`0.5 * n_modes * freq_tol / autorcm`), because a molecule
# with more modes accumulates more of the projection difference in the zero
# point energy (measured 2.98e-7 hartree for CH3OH at tight, 12 modes;
# report 49 section 1.1).  It cannot hide a parameter error: a 1 % frequency
# scaling moves G by ~1e-4 hartree (injection iv of the validation report).
THERMO_TOL = 1.0e-6
#> `autorcm` of xtb (`src/mctc/convert.f90:40`), cm^-1 per hartree
AUTORCM = 219474.63067

# Command line of every landscape task driven by the `ljson` mode.
TASK_ARGS = {
    "opt": ["--lopt", "--lopt-level", "tight"],
    "scan": ["--lscan", "--lscan-coord", "angle", "--lscan-atoms", "2,1,3",
             "--lscan-step", "5", "--lscan-nstep", "1", "--lscan-level", "loose"],
    "scan2d": ["--lscan", "--lscan-ndim", "2",
               "--lscan-coord", "distance", "--lscan-atoms", "1,2",
               "--lscan-step", "2", "--lscan-nstep", "1",
               "--lscan-coord2", "angle", "--lscan-atoms2", "2,1,3",
               "--lscan-step2", "5", "--lscan-nstep2", "1",
               "--lscan-level", "loose"],
    "vib": ["--lvib"],
    "ts": ["--lts", "--lts-method", "prfo"],
    "irc": ["--lirc", "--lirc-method", "hpc"],
}

# The `--ljson` document of a check is named after the *task*, not after the
# check: `scan2d` is a second command line of the `scan` task and writes the
# same `landscape_scan.json`.
DOCUMENT_TASK = {"scan2d": "scan"}


def run(exe, inp, args, rundir, extra_env=None):
    """Run xtb in `rundir` and return (status, combined output).

    `exe`, `inp` and `rundir` are resolved against the caller's working
    directory *before* the child starts.  Without that a relative `rundir`
    makes a relative `inp` unresolvable, because the child's `cwd` is
    `rundir` (M7 fix round finding F5: the failure surfaced as the easily
    misread `No input file given, so there is nothing to do`).
    """
    rundir = os.path.abspath(rundir)
    os.makedirs(rundir, exist_ok=True)
    env = dict(os.environ)
    if extra_env:
        env.update(extra_env)
    proc = subprocess.run(
        [os.path.abspath(exe), os.path.abspath(inp)] + args,
        cwd=rundir,
        capture_output=True,
        text=True,
        # xtb can emit Fortran record markers in its formatted output (a
        # pre-existing upstream defect, `28-stage2-fix-report.md` T2).  A
        # strict decode would make this chain test red for a reason that has
        # nothing to do with its assertions.
        errors="replace",
        env=env,
    )
    return proc.returncode, proc.stdout + proc.stderr


def fail(message, output=""):
    print("FAIL: " + message)
    if output:
        # keep the diagnosis readable: the last lines carry the failure
        lines = output.rstrip().splitlines()
        print("\n".join(lines[-40:]))
    sys.exit(1)


def read(path):
    with open(path) as handle:
        return handle.read()


def parse_float(pattern, text, what, output=""):
    match = re.search(pattern, text, re.M)
    if not match:
        fail("cannot read %s" % what, output or text)
    return float(match.group(1))


def parse_int(pattern, text, what, output=""):
    match = re.search(pattern, text, re.M)
    if not match:
        fail("cannot read %s" % what, output or text)
    return int(match.group(1))


def parse_pair(pattern, text, what, output=""):
    """Two integers of one printed line (the scan frame report)."""
    match = re.search(pattern, text, re.M)
    if not match:
        fail("cannot read %s" % what, output or text)
    return int(match.group(1)), int(match.group(2))


def parse_frequencies(text):
    """Internal vibrational frequencies (cm^-1) of a printed frequency table.

    Both routes print one or more `eigval :` lines; the projected zero modes
    are the leading zeros and are dropped.
    """
    # `--hess` prints the same table twice (the Hessian block and the
    # frequency printout), so only the first contiguous block is used.
    block = re.search(r"(?:^eigval :.*\n)+", text, re.M)
    if not block:
        fail("no frequency table in the output", text)
    values = [float(v) for v in re.findall(r"-?\d+\.\d+", block.group(0))]
    internal = [v for v in values if abs(v) > 1.0e-2]
    if not internal:
        fail("no internal frequency in the printed table", text)
    return internal


def parse_thermo(text):
    """ZPE, H(T) and G(T) (hartree) of a thermochemistry printout."""
    zpe = parse_float(r"zero point energy\s+(-?[\d.]+) Eh", text,
                      "the zero point energy")
    gibbs = parse_float(r"total free energy\s+(-?[\d.]+) Eh", text,
                        "the total free energy")
    match = re.search(r"TOTAL ENTHALPY\s+(-?[\d.]+) Eh", text)
    if match:
        enthalpy = float(match.group(1))
    else:
        # `--lvib` prints the enthalpy in its own summary instead of the
        # property printout of `--hess`.
        enthalpy = parse_float(r"total enthalpy H\(T\):\s+(-?[\d.]+) Eh", text,
                               "the enthalpy")
    return zpe, enthalpy, gibbs


def check_chain(exe, inp, rundir, freq_tol=None):
    """L3: `opt -> hess -> freq -> thermo`, compared against `--hess`.

    `freq_tol` overrides `FREQ_TOL`; the larger molecules of the validation
    document need it because the projection spread grows with the number of
    modes and with |g| (measured at a *tight* minimum: 0.0100 cm^-1 for H2O,
    0.0600 cm^-1 for CH3OH; at the **default** `--lopt` level CH3OH gives
    0.8900 cm^-1, which is what the 1.00 cm^-1 tolerance of the validation
    document has to cover - report 49 section 1.2, finding F2).
    """
    tol = FREQ_TOL if freq_tol is None else float(freq_tol)
    status, output = run(exe, inp, ["--gfn", "2", "--lopt", "--lopt-level",
                                    "tight", "--norestart"], rundir)
    if status != 0:
        fail("the chain optimisation exited %d" % status, output)
    if "landscape: converged" not in output:
        fail("the chain optimisation did not converge", output)
    geometry = os.path.join(rundir, "xtbopt.xyz")
    if not os.path.exists(geometry):
        fail("xtbopt.xyz was not written", output)
    # the chain consumes the optimised geometry, not the input file
    status, output = run(exe, geometry, ["--gfn", "2", "--lvib", "--norestart"],
                         rundir)
    if status != 0:
        fail("the chain frequency step exited %d" % status, output)
    if "landscape task terminated: frequency analysis completed" not in output:
        fail("the chain frequency step did not complete", output)
    vib_zpe, vib_h, vib_g = parse_thermo(output)
    vib_freq = parse_frequencies(output)

    status, output = run(exe, geometry, ["--gfn", "2", "--hess", "--norestart"],
                         rundir)
    if status != 0:
        fail("the --hess reference exited %d" % status, output)
    hess_zpe, hess_h, hess_g = parse_thermo(output)
    hess_freq = parse_frequencies(output)

    if len(vib_freq) != len(hess_freq):
        fail("the two routes report %d and %d internal frequencies"
             % (len(vib_freq), len(hess_freq)), output)
    max_dfreq = max(abs(a - b) for a, b in zip(vib_freq, hess_freq))
    thermo_tol = max(THERMO_TOL, 0.5*len(vib_freq)*tol/AUTORCM)
    if max_dfreq > tol:
        fail("the chain frequencies differ from --hess by %.4f cm^-1 "
             "(tolerance %.2f)" % (max_dfreq, tol), output)
    for name, ours, reference in (("ZPE", vib_zpe, hess_zpe),
                                  ("H(T)", vib_h, hess_h),
                                  ("G(T)", vib_g, hess_g)):
        if abs(ours - reference) > thermo_tol:
            fail("the chain %s differs from --hess by %.3e hartree "
                 "(tolerance %.1e)" % (name, ours - reference, thermo_tol),
                 output)
    print("chain: max dfreq %.4f cm^-1 (tol %.2f), dZPE %.3e, dH %.3e, "
          "dG %.3e hartree (tol %.1e)"
          % (max_dfreq, tol, vib_zpe - hess_zpe, vib_h - hess_h,
             vib_g - hess_g, thermo_tol))
    print("chain: frequencies %s cm^-1"
          % " ".join("%.2f" % v for v in vib_freq))


def check_chain_failure(exe, inp, rundir):
    """The chain must abort when a step fails and must not run the next one.

    Step 1 (the optimisation) succeeds and writes the geometry of the chain.
    Step 2 (the frequency analysis) is misconfigured - it is told to read a
    Hessian file that does not exist - and has to fail.  A runner that checks
    the exit status stops there; a runner that continues would run the
    reference step and leave `hessian`/`vibspectrum` behind, which is what the
    artifact assertion below detects.

    The rundir is removed first: the assertion is the *absence* of an
    artefact, and an artefact left behind by a previous run (a manual run, or
    the injected run of the validation report) would otherwise make this
    check fail although the product behaved correctly (M7 fix round finding
    F1).
    """
    rundir = os.path.abspath(rundir)
    shutil.rmtree(rundir, ignore_errors=True)
    os.makedirs(rundir, exist_ok=True)
    status, output = run(exe, inp, ["--gfn", "2", "--lopt", "--lopt-level",
                                    "tight", "--norestart"], rundir)
    if status != 0:
        fail("the first chain step exited %d" % status, output)
    geometry = os.path.join(rundir, "xtbopt.xyz")
    if not os.path.exists(geometry):
        fail("the first chain step did not write xtbopt.xyz", output)
    status, output = run(exe, geometry, ["--gfn", "2", "--lvib",
                                         "--lvib-hessian", "file",
                                         "--lvib-hessian-file",
                                         "chain_missing_hessian",
                                         "--norestart"], rundir)
    if status == 0:
        fail("a failed chain step exited 0; a caller cannot detect the "
             "failure", output)
    if "error" not in output.lower() and "cannot" not in output.lower():
        fail("the failed step did not report why it failed", output)
    for artifact in ("hessian", "vibspectrum", "g98.out"):
        if os.path.exists(os.path.join(rundir, artifact)):
            fail("the chain continued after a failed step (%s exists)"
                 % artifact, output)
    print("chain-failure: step 2 exited %d and the reference step did not run"
          % status)


def check_ljson(exe, inp, rundir, task):
    """Parse the JSON of one task, cross-assert its counters and pin its bytes."""
    if task not in TASK_ARGS:
        fail("unknown task '%s'" % task)
    doc_task = DOCUMENT_TASK.get(task, task)
    args = TASK_ARGS[task] + ["--norestart", "--ljson"]

    # 1. the document is valid JSON
    status, output = run(exe, inp, args, rundir)
    if status != 0:
        fail("--ljson %s exited %d" % (task, status), output)
    name = "landscape_%s.json" % doc_task
    path = os.path.join(rundir, name)
    if not os.path.exists(path):
        fail("%s was not written" % name, output)
    try:
        with open(path) as handle:
            document = json.load(handle)
    except ValueError as error:
        fail("%s is not valid JSON: %s" % (name, error), read(path))
    if "wrote " + name not in output:
        fail("the task did not announce %s" % name, output)
    if document.get("task") != doc_task:
        fail("%s declares task '%s'" % (name, document.get("task")), read(path))
    if document.get("status") != 0:
        fail("%s reports status %s" % (name, document.get("status")), read(path))
    if list(document) != ["format", "task", "status", "status_name",
                          "converged", "message", "n_atoms", "energy",
                          "geometry", "provenance", "counters", "values",
                          "flags", "history"]:
        fail("%s has an unexpected key order: %s" % (name, list(document)),
             read(path))

    # 2. every counter is the number the task printed on stdout
    counters = document["counters"]
    values = document["values"]
    expected = {
        "opt": {
            "n_steps": (r"^\s+steps:\s+(\d+)", int),
            "n_eval": (r"energy calls:\s+(\d+)", int),
            "n_cache_miss": (r"cache misses:\s+(\d+)", int),
        },
        "ts": {
            "n_steps": (r"^\s+steps:\s+(\d+)", int),
            "n_eval": (r"energy calls:\s+(\d+)", int),
            "n_cache_miss": (r"cache misses:\s+(\d+)", int),
            "n_imag": (r"saddle confirmation - imaginary modes:\s+(\d+)", int),
        },
        "scan": {
            # n_eval is the engine total of the run; the M7 validation round
            # found that this check did not cross-assert it (finding F4).
            "n_eval": (r"energy evaluations\s+(\d+)", int),
            "n_cache_miss": (r"energy evaluations\s+(\d+)", int),
        },
        "irc": {
            "n_fwd": (r"forward points:\s+(\d+)", int),
            "n_bwd": (r"backward points:\s+(\d+)", int),
        },
    }
    for key, (pattern, cast) in expected.get(task, {}).items():
        printed = parse_int(pattern, output, "the printed %s" % key, output)
        if counters.get(key) != printed:
            fail("%s: counter %s is %s but stdout says %d"
                 % (name, key, counters.get(key), printed), read(path))

    # per-task value cross-assertions
    if task == "opt":
        printed = parse_float(r"final energy:\s+(-?[\d.]+) Eh", output,
                              "the final energy")
        # the stdout table prints 12 decimals, the document 17 digits
        if abs(values["energy"] - printed) > 1.0e-9:
            fail("%s: energy is %r but stdout says %r"
                 % (name, values["energy"], printed), read(path))
    elif task in ("scan", "scan2d"):
        printed = parse_float(
            r"lowest grid energy:\s+(-?[\d.eE+-]+) Eh" if task == "scan2d"
            else r"lowest energy:\s+(-?[\d.eE+-]+) Eh",
            output, "the lowest energy")
        # the scan table prints 6 significant digits
        if abs(values["energy"] - printed) > 1.0e-5:
            fail("%s: energy is %r but stdout says %r"
                 % (name, values["energy"], printed), read(path))
        frames = read(os.path.join(rundir, "xtbscan.log")).count("energy:")
        if counters.get("n_points") != frames:
            fail("%s: n_points is %s but xtbscan.log holds %d frames"
                 % (name, counters.get("n_points"), frames), read(path))
        # The point bookkeeping of the 1D profile and of the 2D grid is the
        # same three counters (n_fail / n_unconverged / n_valid); the scan
        # task prints the written frame count and the omitted points, so the
        # document can be checked against the run (M7 fix round, findings
        # L4-L6: none of the three was cross-asserted before).  The frame
        # count is the number of points that carry an energy, which is
        # `n_valid` (converged *and* evaluated) plus the unconverged ones, so
        # the exact relation is a sandwich, not an equality.
        written, omitted = parse_pair(
            r"wrote xtbscan\.log with (\d+) frame\(s\), (\d+) failed "
            r"point\(s\) omitted", output, "the scan frame report", output)
        if counters.get("n_points") != written + omitted:
            fail("%s: the scan has %s points but %d frame(s) plus %d omitted "
                 "point(s)" % (name, counters.get("n_points"), written,
                                omitted), read(path))
        if counters.get("n_fail") != omitted:
            fail("%s: n_fail is %s but the run reports %d failed point(s)"
                 % (name, counters.get("n_fail"), omitted), read(path))
        if not counters.get("n_valid", -1) <= written <= \
                counters.get("n_valid", -1) + counters.get("n_unconverged", 0):
            fail("%s: %d frame(s) is outside [n_valid, n_valid + "
                 "n_unconverged] = [%s, %s]"
                 % (name, written, counters.get("n_valid"),
                    counters.get("n_valid", 0) +
                    counters.get("n_unconverged", 0)), read(path))
        if not 0 <= counters.get("n_unconverged", -1) <= \
                counters.get("n_points"):
            fail("%s: n_unconverged is %s out of %s points"
                 % (name, counters.get("n_unconverged"),
                    counters.get("n_points")), read(path))
        if task == "scan2d":
            n1, n2 = counters.get("n1"), counters.get("n2")
            if n1 * n2 != counters.get("n_points"):
                fail("%s: the grid is %s x %s but n_points is %s"
                     % (name, n1, n2, counters.get("n_points")), read(path))
            if not re.search(r"2D scan grid:\s+%d\s+x%d points" % (n1, n2),
                             output):
                fail("%s: stdout does not report a %s x %s grid"
                     % (name, n1, n2), output)
    elif task == "vib":
        zpe = parse_float(r"zero point energy:\s+(-?[\d.eE+-]+) Eh", output,
                          "the zero point energy")
        # the task prints the zero point energy with 6 significant digits
        if abs(values["zpe"] - zpe) > 1.0e-7:
            fail("%s: zpe is %r but stdout says %r" % (name, values["zpe"], zpe),
                 read(path))
        if counters.get("n_imag") != 0:
            fail("%s: the water minimum has %s imaginary modes"
                 % (name, counters.get("n_imag")), read(path))
    elif task == "irc":
        printed = parse_float(r"transition state energy:\s+(-?[\d.]+) Eh",
                              output, "the transition state energy")
        if abs(values["energy_ts"] - printed) > 1.0e-9:
            fail("%s: energy_ts is %r but stdout says %r"
                 % (name, values["energy_ts"], printed), read(path))
        if counters.get("n_imag") != 1:
            fail("%s: the transition state has %s imaginary modes"
                 % (name, counters.get("n_imag")), read(path))
    elif task == "ts":
        if counters.get("n_imag") != 1:
            fail("%s: the saddle has %s imaginary modes"
                 % (name, counters.get("n_imag")), read(path))
        # M7 fix round M-1: the printed counter and the engine report of the
        # *same* stdout have to be the same number.  Both are the complete
        # run, confirmation Hessian included; the dedicated mode
        # `ts-counters` additionally pins the size of that Hessian.
        calls = parse_int(r"landscape engine: energy evaluations\s+(\d+)",
                          output, "the engine evaluation count", output)
        if counters.get("n_cache_miss") != calls:
            fail("%s: the counter is %s but the engine of the same run "
                 "reports %d energy evaluations"
                 % (name, counters.get("n_cache_miss"), calls), output)

    # 3. determinism: a second run produces byte identical output
    again = rundir + "-again"
    status, output2 = run(exe, inp, args, again)
    if status != 0:
        fail("the second --ljson %s run exited %d" % (task, status), output2)
    first_bytes = read(path)
    second_bytes = read(os.path.join(again, name))
    if first_bytes != second_bytes:
        fail("%s is not reproducible between two runs" % name,
             "first:\n%s\nsecond:\n%s" % (first_bytes, second_bytes))

    # 4. the file name goes through the namespace
    namespaced = rundir + "-ns"
    status, output3 = run(exe, inp, args + ["--namespace", "e2ens"], namespaced)
    if status != 0:
        fail("--ljson %s --namespace exited %d" % (task, status), output3)
    ns_name = "e2ens." + name
    if not os.path.exists(os.path.join(namespaced, ns_name)):
        fail("--namespace did not prefix %s" % name, output3)
    print("ljson: %s valid, counters cross-checked, byte identical, "
          "namespace ok" % task)


# The `method` counter every TS / IRC document carries (the identifier of the
# settings type, see the doc comments of `write_json` in tasks/ts.f90 and
# tasks/irc.f90).
TS_METHODS = [("prfo", 1), ("bprfo", 2), ("dimer", 3), ("neb", 4),
              ("cineb", 5), ("autoneb", 6), ("string", 7), ("gsm", 8)]
IRC_METHODS = [("gs", 1), ("hpc", 2), ("lqa", 3), ("eulerpc", 4)]


def check_ljson_methods(exe, inp, rundir, product):
    """`--ljson` must produce a valid document for every task method.

    The `ts` task has eight methods and the `irc` task four; four of the eight
    `ts` methods (the band searches) take a different code path in the task and
    return early, so the document of those methods is exercised separately.
    The searches run with a small budget: the assertion is the document (task
    name, method identifier, parseability), not convergence, so a run that
    stops at the budget is fine as long as its document exists.
    """
    os.makedirs(rundir, exist_ok=True)
    product_file = os.path.join(rundir, "product.xyz")
    with open(product_file, "w") as handle:
        handle.write(product)

    def document(name):
        path = os.path.join(rundir, name)
        if not os.path.exists(path):
            fail("%s was not written" % name)
        try:
            with open(path) as handle:
                return json.load(handle)
        except ValueError as error:
            fail("%s is not valid JSON: %s" % (name, error), read(path))

    for method, identifier in TS_METHODS:
        args = ["--gfn", "2", "--lts", "--lts-method", method,
                "--lts-max-steps", "20", "--lts-level", "superloose",
                "--norestart", "--ljson", "--namespace", "mj"]
        if method in ("neb", "cineb", "autoneb", "string", "gsm"):
            args += ["--lts-product", product_file]
        run(exe, inp, args, rundir)
        doc = document("mj.landscape_ts.json")
        if doc.get("task") != "ts":
            fail("the %s document declares task %r" % (method, doc.get("task")),
                 read(os.path.join(rundir, "mj.landscape_ts.json")))
        if doc["counters"].get("method") != identifier:
            fail("the %s document reports method %r, expected %d"
                 % (method, doc["counters"].get("method"), identifier),
                 read(os.path.join(rundir, "mj.landscape_ts.json")))
        os.remove(os.path.join(rundir, "mj.landscape_ts.json"))

    for method, identifier in IRC_METHODS:
        args = ["--gfn", "2", "--lirc", "--lirc-method", method,
                "--lirc-max-steps", "4", "--norestart", "--ljson",
                "--namespace", "mj"]
        run(exe, inp, args, rundir)
        doc = document("mj.landscape_irc.json")
        if doc.get("task") != "irc":
            fail("the %s document declares task %r" % (method, doc.get("task")),
                 read(os.path.join(rundir, "mj.landscape_irc.json")))
        if doc["counters"].get("method") != identifier:
            fail("the %s document reports method %r, expected %d"
                 % (method, doc["counters"].get("method"), identifier),
                 read(os.path.join(rundir, "mj.landscape_irc.json")))
        os.remove(os.path.join(rundir, "mj.landscape_irc.json"))

    print("ljson-methods: 8 ts methods and 4 irc methods write their document")


def check_ts_counters(exe, inp, rundir):
    """The sampling point of the `ts` evaluation counter (PLAN section 12.11).

    The `ts` task samples `result%n_cache_miss` *after* the frequency
    confirmation (M7 fix round, finding M-1), so the counter of a run is the
    cost of the complete run and equals the `landscape engine: energy
    evaluations` line of the same stdout.  The check pins both properties:

      * for one run, the printed `cache misses`, the engine's own total and
        the JSON counter are one and the same number;
      * `--lts-no-verify` removes exactly the confirmation Hessian, i.e. whole
        6N finite difference evaluations (6N for N atoms, central
        differences of the 3N coordinates).  A counter sampled before the
        confirmation gives a difference of 0 here and a `cache misses` that
        disagrees with the engine report, so each of the three assertions
        below fails on its own.
    """
    document = "landscape_ts.json"
    counters = {}
    for label, extra in (("confirmed", []), ("unconfirmed", ["--lts-no-verify"])):
        status, output = run(exe, inp, ["--gfn", "2", "--lts", "--lts-method",
                                        "prfo", "--norestart", "--ljson"] + extra,
                             rundir + "-" + label)
        if status != 0:
            fail("the %s ts run exited %d" % (label, status), output)
        printed = parse_int(r"cache misses:\s+(\d+)", output,
                            "the printed cache misses", output)
        calls = parse_int(r"landscape engine: energy evaluations\s+(\d+)",
                          output, "the engine evaluation count", output)
        if printed != calls:
            fail("%s: the same stdout prints 'cache misses: %d' and the "
                 "engine reports %d energy evaluations"
                 % (label, printed, calls), output)
        try:
            with open(os.path.join(rundir + "-" + label, document)) as handle:
                doc = json.load(handle)
        except (OSError, ValueError) as error:
            fail("cannot read the %s document: %s" % (label, error))
        if doc["counters"].get("n_cache_miss") != printed:
            fail("%s: the JSON counter is %s but stdout prints %d"
                 % (label, doc["counters"].get("n_cache_miss"), printed))
        counters[label] = printed

    # 6N central difference evaluations per confirmation Hessian of N atoms
    atoms = int(open(os.path.abspath(inp)).readline().split()[0])
    hessian = 6 * atoms
    difference = counters["confirmed"] - counters["unconfirmed"]
    if difference <= 0 or difference % hessian != 0:
        fail("the confirmed ts run costs %d evaluations more than the run "
             "without --lts-no-verify; expected a positive multiple of one "
             "6N confirmation Hessian (%d evaluations for %d atoms)"
             % (difference, hessian, atoms))
    print("ts-counters: %d vs %d evaluations, difference %d = %d x 6N "
          "(N=%d); both agree with the engine report"
          % (counters["confirmed"], counters["unconfirmed"], difference,
             difference // hessian, atoms))


def check_scan_counters(exe, inp, rundir):
    """The point bookkeeping of a scan, converged and unconverged.

    The document of the `scan` task carries four point counters (`n_points`,
    `n_valid`, `n_fail`, `n_unconverged`).  They are cross-asserted against
    the two numbers the task prints about `xtbscan.log` (frames written,
    points omitted) in two runs:

      * a converged 1D scan - every point usable;
      * a scan whose inner relaxation never converges (`--lscan-max-steps 1`),
        which is the documented non-zero exit case and gives non-zero
        `n_unconverged` with `n_valid = 0`.
    """
    args = TASK_ARGS["scan"]
    document = "landscape_scan.json"

    status, output = run(exe, inp, args + ["--norestart", "--ljson"],
                         rundir + "-converged")
    if status != 0:
        fail("the converged scan exited %d" % status, output)
    with open(os.path.join(rundir + "-converged", document)) as handle:
        doc = json.load(handle)
    counters = doc["counters"]
    written, omitted = parse_pair(
        r"wrote xtbscan\.log with (\d+) frame\(s\), (\d+) failed point\(s\) "
        r"omitted", output, "the scan frame report", output)
    if doc["status"] != 0 or not doc["converged"]:
        fail("the converged scan reports status %s / converged %s"
             % (doc["status"], doc["converged"]))
    if counters["n_unconverged"] != 0 or counters["n_fail"] != 0:
        fail("the converged scan reports %d unconverged and %d failed points"
             % (counters["n_unconverged"], counters["n_fail"]), output)
    if counters["n_valid"] != written or counters["n_valid"] != counters["n_points"]:
        fail("the converged scan reports n_valid %s, n_points %s but wrote %d "
             "frame(s)" % (counters["n_valid"], counters["n_points"], written))
    summary_counts = dict(counters)

    status, output = run(exe, inp, args + ["--lscan-max-steps", "1",
                                           "--norestart", "--ljson"],
                         rundir + "-unconverged")
    if status == 0:
        fail("a scan without a converged point exited 0", output)
    with open(os.path.join(rundir + "-unconverged", document)) as handle:
        doc = json.load(handle)
    counters = doc["counters"]
    written, omitted = parse_pair(
        r"wrote xtbscan\.log with (\d+) frame\(s\), (\d+) failed point\(s\) "
        r"omitted", output, "the scan frame report", output)
    if doc["status"] != -1 or doc["status_name"] != "error":
        fail("the unconverged scan reports status %s / %s"
             % (doc["status"], doc["status_name"]))
    if counters["n_unconverged"] != counters["n_points"]:
        fail("the unconverged scan reports %d of %s points unconverged"
             % (counters["n_unconverged"], counters["n_points"]), output)
    if counters["n_valid"] != 0 or counters["n_fail"] != omitted:
        fail("the unconverged scan reports n_valid %s / n_fail %s "
             "(%d point(s) omitted)"
             % (counters["n_valid"], counters["n_fail"], omitted), output)
    if written != counters["n_points"]:
        fail("the unconverged scan wrote %d frame(s) for %s points"
             % (written, counters["n_points"]), output)
    print("scan-counters: converged run %d valid of %d points, unconverged "
          "run %d unconverged of %s points with status error and non-zero "
          "exit"
          % (summary_counts["n_valid"], summary_counts["n_points"],
             counters["n_unconverged"], counters["n_points"]))


#> The product of the ammonia inversion fixture: the same four atoms with the
#> nitrogen below the H3 plane.  Written here so the mode needs one input file.
NH3_DOWN = ("4\nammonia (N below the H3 plane)\n"
            "N  0.0000000  0.0000000 -0.3800000\n"
            "H  0.9400000  0.0000000  0.0000000\n"
            "H -0.4700000  0.8140000  0.0000000\n"
            "H -0.4700000 -0.8140000  0.0000000\n")

#> A linear HNC, i.e. the molecule's own minimum.  With the escape out of a
#> minimum switched off, the P-RFO search converges **onto that minimum**,
#> which is the H2 fixture of `xtb/landscape-ts-minimum`; written here so the
#> status mode needs no second input file.
HNC_LINEAR = ("3\nlinear HNC (a minimum)\n"
              "N  0.0000000  1.1560000  0.0000000\n"
              "C  0.0000000  0.0000000  0.0000000\n"
              "H  0.0000000 -1.0500000  0.0000000\n")


def check_ts_failure(exe, inp, rundir):
    """A failed `--lts` search keeps its own message and its own cost.

    The confirmation Hessian is only built for a **converged** search (report
    221, findings S2 and H1).  Two fixtures exercise both branches of the task:
    a single geometry P-RFO search stopped after one step, and a 9-image NEB
    band stopped after three iterations.  Each is run with and without
    `--lts-no-verify`; the two runs must cost exactly the same number of
    evaluations, because a failed search has nothing left to confirm.

    Every assertion below is red on the pre-221 task: the confirmation cost
    one 6N Hessian in both cases (50 -> 26 and 61 -> 31 evaluations for the
    four-atom fixture), the run printed the confirmation's verdict, and the
    document carried `n_imag` = 0 with `status` = error instead of the
    search's own message with `status` = incomplete.
    """
    os.makedirs(rundir, exist_ok=True)
    product = os.path.join(os.path.abspath(rundir), "product-nh3-down.xyz")
    with open(product, "w") as handle:
        handle.write(NH3_DOWN)

    atoms = int(open(os.path.abspath(inp)).readline().split()[0])
    images = 9
    steps = 3
    # The search-only cost of the two fixtures, derived rather than recorded:
    # the single geometry search spends one point on the initial gradient, the
    # 6N initial P-RFO Hessian and one step; the band spends one evaluation per
    # image on the initial band, `steps` sweeps over the **internal** images
    # (the two fixed endpoints are evaluated once, on the initial band, and
    # cached, report 223) and one on the candidate.  The 6N confirmation
    # Hessian is not in either number.
    single_cost = 1 + 6 * atoms + 1
    band_cost = 1 + images + steps * (images - 2)
    cases = (
        ("single", ["--lts", "--lts-method", "prfo", "--lts-max-steps", "1"],
         "maximum number of steps reached", single_cost),
        ("band", ["--lts", "--lts-method", "neb", "--lts-images", str(images),
                  "--lts-max-steps", str(steps), "--lts-product", product],
         "band reached the iteration limit without converging", band_cost),
    )
    # Every run below exits non-zero on purpose and is then judged by its
    # document, so a directory left by an earlier run would be read as this
    # run's result when the child died before writing - the stale-state
    # flakiness AGENTS.md section 6 records for `scan-failed-point`.
    for label, _args, _message, _cost in cases:
        for variant in ("confirmed", "unconfirmed"):
            shutil.rmtree(rundir + "-%s-%s" % (label, variant),
                          ignore_errors=True)
    for label, args, message, expected in cases:
        counters = {}
        for variant, extra in (("confirmed", []),
                               ("unconfirmed", ["--lts-no-verify"])):
            status, output = run(exe, inp, ["--gfn", "2"] + args +
                                 ["--norestart", "--ljson"] + extra,
                                 "%s-%s-%s" % (rundir, label, variant))
            if status == 0:
                fail("the %s %s search exited 0" % (label, variant), output)
            printed = parse_int(r"cache misses:\s+(\d+)", output,
                                "the printed cache misses", output)
            calls = parse_int(r"landscape engine: energy evaluations\s+(\d+)",
                              output, "the engine evaluation count", output)
            if printed != calls:
                fail("%s %s: stdout prints 'cache misses: %d' and the engine "
                     "reports %d energy evaluations"
                     % (label, variant, printed, calls), output)
            path = os.path.join(rundir + "-%s-%s" % (label, variant),
                                "landscape_ts.json")
            try:
                with open(path) as handle:
                    doc = json.load(handle)
            except (OSError, ValueError) as error:
                fail("cannot read the %s %s document: %s"
                     % (label, variant, error), output)
            counters[variant] = printed

            if doc["converged"]:
                fail("the %s %s search reports converged" % (label, variant),
                     output)
            if doc["status"] != 1:
                fail("%s %s: status is %s, expected 1 (incomplete); -1 is the "
                     "'converged on something that is not a saddle' verdict"
                     % (label, variant, doc["status"]), output)
            if doc["counters"].get("n_imag") != -1:
                fail("%s %s: n_imag is %s, expected -1 (the confirmation was "
                     "not performed)"
                     % (label, variant, doc["counters"].get("n_imag")), output)
            if message not in (doc.get("message") or ""):
                fail("%s %s: the search lost its own message, the document "
                     "says %r, expected %r"
                     % (label, variant, doc.get("message"), message), output)
            if "saddle confirmation" in output:
                fail("%s %s: the failed search still built the confirmation "
                     "Hessian" % (label, variant), output)
            if doc["counters"].get("n_cache_miss") != printed:
                fail("%s %s: the JSON counter is %s but stdout prints %d"
                     % (label, variant, doc["counters"].get("n_cache_miss"),
                        printed), output)

        if counters["confirmed"] != counters["unconfirmed"]:
            fail("%s: the failed search costs %d evaluations with the "
                 "confirmation and %d without it; a failed run has nothing "
                 "left to confirm"
                 % (label, counters["confirmed"], counters["unconfirmed"]))
        if counters["confirmed"] != expected:
            fail("%s: the failed search costs %d evaluations, the pinned "
                 "search-only cost is %d"
                 % (label, counters["confirmed"], expected))
    print("ts-failure: a failed single-geometry search and a failed band "
          "keep their own message and cost the search only (%d and %d "
          "evaluations)" % (cases[0][3], cases[1][3]))


def check_ts_status(exe, inp, rundir):
    """The two `--lts` failure statuses of the document, side by side.

    `landscape_json_status_incomplete` (1) means "the result is written but a
    criterion was not met"; `landscape_json_status_error` (-1) means "an error
    was reported through `env`" (`src/landscape/types/json.f90`).  The task has
    exactly two ways to fail a search and both must map onto that contract:

      * a search that ran out of steps keeps its own message and reports 1 -
        the search found nothing to confirm, so `n_imag` stays at its "not
        performed" value -1;
      * a search that **converged onto a minimum** reports -1 with the
        confirmation's message, because the task raises that verdict through
        `env` like every other error.  This branch is easy to lose: the
        confirmation clears `result%converged` when it rejects the point, so a
        status predicate written on `converged` never sees it and downgrades
        the document to 1, which is what report 224 (finding 1) found and this
        check pins.
    """
    os.makedirs(rundir, exist_ok=True)
    for name in ("exhausted", "rejected"):
        shutil.rmtree(os.path.join(rundir, name), ignore_errors=True)
    hnc = os.path.join(os.path.abspath(rundir), "hnc-linear.xyz")
    with open(hnc, "w") as handle:
        handle.write(HNC_LINEAR)

    def document(path, output):
        try:
            with open(path) as handle:
                return json.load(handle)
        except (OSError, ValueError) as error:
            fail("cannot read %s: %s" % (path, error), output)

    # (a) exhausted: the search's own message, status incomplete
    status, output = run(exe, inp, ["--gfn", "2", "--lts", "--lts-method",
                                    "prfo", "--lts-max-steps", "1",
                                    "--norestart", "--ljson"],
                         os.path.join(rundir, "exhausted"))
    if status == 0:
        fail("the exhausted search exited 0", output)
    doc = document(os.path.join(rundir, "exhausted", "landscape_ts.json"),
                   output)
    if doc.get("status") != 1 or doc.get("status_name") != "incomplete":
        fail("the exhausted search reports %s / %s, expected 1 / incomplete"
             % (doc.get("status"), doc.get("status_name")), output)
    if "maximum number of steps reached" not in (doc.get("message") or ""):
        fail("the exhausted search lost its own message, the document says %r"
             % doc.get("message"), output)
    if doc["counters"].get("n_imag") != -1:
        fail("the exhausted search reports n_imag %s, expected -1"
             % doc["counters"].get("n_imag"), output)

    # (b) converged onto a minimum: the confirmation's message, status error
    status, output = run(exe, hnc, ["--gfn", "2", "--lts", "--lts-method",
                                    "prfo", "--lts-escape-step", "0",
                                    "--norestart", "--ljson"],
                         os.path.join(rundir, "rejected"))
    if status == 0:
        fail("a minimum was accepted as a transition state", output)
    doc = document(os.path.join(rundir, "rejected", "landscape_ts.json"),
                   output)
    if doc.get("status") != -1 or doc.get("status_name") != "error":
        fail("the rejected point reports %s / %s, expected -1 / error; 1 is "
             "the verdict of a search that did not converge"
             % (doc.get("status"), doc.get("status_name")), output)
    if doc.get("converged"):
        fail("the rejected point still reports converged", output)
    if doc["counters"].get("n_imag") != 0:
        fail("the rejected point reports n_imag %s, expected 0"
             % doc["counters"].get("n_imag"), output)
    if "the located stationary point is a minimum" not in (doc.get("message")
                                                           or ""):
        fail("the rejected point lost the confirmation's message, the "
             "document says %r" % doc.get("message"), output)
    print("ts-status: an exhausted search is status 1 with its own message, a "
          "search that converged onto a minimum is status -1 with the "
          "confirmation's verdict")


def check_layers(xtb, tester, inp, rundir, ts_input=None):
    """Run L1, L2 and L3 and write a structured summary.

    The suite lists are the complete landscape unit suites of
    `test/unit/main.F90` (M7 fix round finding F7: the first version listed 6
    of the 14 suites and claimed L1+L2+L3 in one report).  The second L3 leg -
    the `ts -> irc` chain - needs the transition state input; it is run when
    the caller passes one.
    """
    rundir = os.path.abspath(rundir)
    os.makedirs(rundir, exist_ok=True)
    summary = {"layers": []}
    l1_suites = ["landscape-api", "landscape-settings", "landscape-coords",
                 "landscape-hessian-update", "landscape-optimize-math",
                 "landscape-cg", "landscape-hessian", "landscape-scan",
                 "landscape-saddle", "landscape-path", "landscape-irc"]
    l2_suites = ["landscape-optimize-xtb", "landscape-vib",
                 "landscape-failure", "landscape-e2e"]
    for suite in l1_suites:
        proc = subprocess.run([os.path.abspath(tester), suite], cwd=rundir,
                              capture_output=True, text=True)
        summary["layers"].append({"layer": "L1", "suite": suite,
                                  "status": proc.returncode,
                                  "passed": proc.returncode == 0})
    for suite in l2_suites:
        proc = subprocess.run([os.path.abspath(tester), suite], cwd=rundir,
                              capture_output=True, text=True)
        summary["layers"].append({"layer": "L2", "suite": suite,
                                  "status": proc.returncode,
                                  "passed": proc.returncode == 0})
    chain_dir = os.path.join(rundir, "chain")
    try:
        check_chain(xtb, inp, chain_dir)
        summary["layers"].append({"layer": "L3", "suite": "chain",
                                  "status": 0, "passed": True})
    except SystemExit:
        summary["layers"].append({"layer": "L3", "suite": "chain",
                                  "status": 1, "passed": False})
    if ts_input:
        try:
            summary["ts_irc"] = check_ts_irc(xtb, ts_input,
                                             os.path.join(rundir, "ts-irc"))
            summary["layers"].append({"layer": "L3", "suite": "ts-irc",
                                      "status": 0, "passed": True})
        except SystemExit:
            summary["layers"].append({"layer": "L3", "suite": "ts-irc",
                                      "status": 1, "passed": False})
    with open(os.path.join(rundir, "layers.json"), "w") as handle:
        json.dump(summary, handle, indent=1, sort_keys=True)
    if not all(entry["passed"] for entry in summary["layers"]):
        fail("at least one layer failed", json.dumps(summary, indent=1))
    print("layers: %d checks passed" % len(summary["layers"]))


def check_ts_irc(xtb, inp, rundir):
    """L3: the `ts -> irc` chain, returned as a topology record.

    Shared by the `layers` and `maple` modes (the second used to carry its own
    copy).  The assertion is topology: the search reports a first order saddle
    and the integration reaches a minimum in both directions.
    """
    rundir = os.path.abspath(rundir)
    status, output = run(xtb, inp, ["--gfn", "2", "--lts", "--lts-method",
                                    "prfo", "--norestart"], rundir)
    if status != 0:
        fail("the HCN transition state search failed", output)
    if not re.search(r"saddle confirmation - imaginary modes:\s+1", output):
        fail("the HCN candidate is not a first order saddle", output)
    status, output = run(xtb, inp, ["--gfn", "2", "--lirc", "--lirc-method",
                                    "hpc", "--norestart"], rundir)
    if status != 0:
        fail("the HCN IRC failed", output)
    endpoints = re.findall(r"endpoint energy:\s+(-?[\d.]+) Eh", output)
    if len(endpoints) != 2:
        fail("the HCN IRC reports %d endpoint(s)" % len(endpoints), output)
    return {
        "saddle_points": 1,
        "endpoint_minima": len(endpoints),
        "endpoint_energies_hartree": [float(v) for v in endpoints],
        "reached_minimum_forward": bool(re.search(r"forward reached minimum: T",
                                                  output)),
        "reached_minimum_backward": bool(re.search(r"backward reached minimum: T",
                                                   output)),
    }


# One committed fixture per task of the MAPLE example tree.  The suite cannot
# *run* any of them (their potentials are machine learning models whose
# checkpoints are not part of the checkout), so what is recorded is a property
# of the committed file: the number of frames, the atom count and the extrema
# of the committed energy series.  PLAN section 9.8 asks for this record "per
# task"; before the final fix round only the IRC fixture was read (finding L-8
# of `50-m7-independent-review.md`).
MAPLE_FIXTURES = [
    ("sp", "sp/sp.out", "single point output of the top level `sp` example"),
    ("opt", "opt/lbfgs/inp1_traj.xyz", "L-BFGS trajectory of `opt/lbfgs` example 1"),
    ("scan", "scan/C18_scan_final.xyz", "relaxed scan of the `scan` C18 example"),
    ("ts", "ts/dmf/inp1_prfo_ts.xyz", "P-RFO transition state of `ts/dmf` example 1"),
    ("freq", "freq/mw/inp1.xyz", "frequency example `freq/mw` 1"),
    ("irc", "irc/hpc/inp1_full.xyz", "HPC IRC path (`full` direction) example 1"),
]


def read_xyz_energies(path):
    """Atom count and the `Energy = ` series of a committed xyz trajectory.

    The MAPLE fixtures write one frame per image and put the energy of the
    frame into its comment line; a frame without an energy keeps `None` so the
    series length still equals the frame count.
    """
    with open(path) as handle:
        lines = handle.readlines()
    natoms = 0
    series = []
    i = 0
    while i < len(lines):
        fields = lines[i].split()
        if not fields:
            break
        per_frame = int(fields[0])
        if natoms == 0:
            natoms = per_frame
        match = re.search(r"Energy =\s*(-?[\d.]+)", lines[i + 1])
        series.append(float(match.group(1)) if match else None)
        i += per_frame + 2
    return natoms, series


def read_sp_output(path):
    """Atom count and total energy of a committed single point output."""
    natoms = 0
    energy = None
    with open(path) as handle:
        for line in handle:
            match = re.search(r"^Energy:\s*(-?[\d.]+)", line)
            if match:
                energy = float(match.group(1))
            if re.match(r"^\d+\s+[A-Z][a-z]?\s+-?[\d.]", line):
                natoms += 1
    return natoms, energy


def check_maple(xtb, inp, rundir):
    """Record the MAPLE comparison methodology and the comparable topology.

    MAPLE (MAchine-learning Potential for Landscape Exploration,
    `/home/owenchan/work_path/coding/python/MAPLE`) implements the same task
    matrix as this suite, but **every fixture of its example tree uses a
    machine learning potential** (ANI-1xnr, AIMNet2, MACE, ...), so its
    numbers are not a GFN2 baseline and PLAN N9/R-14 forbid a numerical
    comparison.  What can be compared is topology: the number of minima and
    saddle points a path connects.

    This mode therefore does three things:

      1. it **really probes** whether MAPLE can be used in this environment -
         import the dispatcher and instantiate the object a task run needs,
         not merely `import maple` (finding M-1 of the fix round, and F6 of
         the validation round: the shallow import returns 0 although
         `maple.function.dispatcher` cannot even be imported),
      2. it reads one *committed* fixture per task of the MAPLE example tree
         (`MAPLE_FIXTURES` below) and records the frames, the atoms and the
         extrema of the committed energy series, so the reference records of
         PLAN section 9.8 exist per task even though MAPLE itself cannot be
         run here, and
      3. it runs this suite's HCN/HNC `ts -> irc` chain and reports its
         topology.

    The mode never claims MAPLE was executed, and it asserts the probe
    instead of merely recording it: `maple_importable` is only true when the
    dispatcher was imported **and** an object was instantiated, a failing
    probe must carry a non-empty diagnosis naming the missing module, and a
    probe that succeeds here fails the check, because this checkout's whole
    MAPLE record (and the validation document) states that MAPLE was *not*
    run.  A false `true` next to an empty error string is the most dangerous
    kind of defect this project can ship: it reads as "validated against
    MAPLE" without any evidence.
    """
    maple_dir = os.environ.get("MAPLE_DIR",
                               "/home/owenchan/work_path/coding/python/MAPLE")
    record = {
        "maple_dir": maple_dir,
        "maple_executed": False,
        "maple_importable": False,
        "comparison": "topology only (MAPLE fixtures use ML potentials, "
                      "PLAN R-14 / N9)",
    }
    if os.path.isdir(maple_dir):
        record.update(probe_maple(maple_dir))
    else:
        record["maple_import_error"] = "checkout not found"
        record["maple_missing_dependencies"] = []

    # One committed fixture per task, read - never run.  The path is relative
    # to `<maple_dir>/examples`, so a reader of `maple.json` sees both the file
    # and what was extracted from it.
    fixtures = {}
    without_fixture = []
    for task, relative, fixture_note in MAPLE_FIXTURES:
        path = os.path.join(maple_dir, "examples", *relative.split("/"))
        entry = {"path": "examples/" + relative, "note": fixture_note}
        if not os.path.exists(path):
            entry["missing"] = True
            without_fixture.append(task)
        elif relative.endswith(".xyz"):
            natoms, series = read_xyz_energies(path)
            values = [value for value in series if value is not None]
            entry["n_atoms"] = natoms
            entry["n_frames"] = len(series)
            entry["energy_extrema"] = count_extrema(values) if values else None
        else:
            natoms, energy = read_sp_output(path)
            entry["n_atoms"] = natoms
            entry["energy_hartree"] = energy
        fixtures[task] = entry
    record["maple_fixtures"] = fixtures
    record["maple_fixtures_note"] = (
        "committed output files of the MAPLE example tree, read without "
        "running MAPLE; every one of them drives a machine learning "
        "potential, so only topology and order of magnitude are comparable "
        "(PLAN R-14 / N9).  `energy_extrema` counts the interior minima and "
        "maxima of the committed energy series.")
    if os.path.isdir(maple_dir) and without_fixture:
        fail("the MAPLE checkout at %s no longer holds the committed fixtures "
             "of %s; the per-task record of PLAN section 9.8 has to be "
             "revised" % (maple_dir, ", ".join(without_fixture)))

    # The IRC fixture keeps its own keys: the validation document cites it.
    irc_fixture = fixtures.get("irc", {})
    if irc_fixture.get("missing", False):
        record["maple_fixture"] = None
    else:
        record["maple_fixture"] = os.path.join(maple_dir, "examples", "irc",
                                               "hpc", "inp1_full.xyz")
        record["maple_fixture_npoints"] = irc_fixture["n_frames"]
        record["maple_fixture_extrema"] = irc_fixture["energy_extrema"]
        record["maple_fixture_note"] = (
            "the committed MAPLE IRC path of a 7-atom molecule (its own "
            "example, not HCN/HNC); the extrema count is a property of the "
            "committed file, not of a run of MAPLE on this machine")

    record["maple_not_executed_reason"] = (
        "no MAPLE job was run here.  The dispatcher itself cannot be "
        "imported in this environment (missing Python packages: %s; %s), so "
        "no task can be dispatched at all; and even with those packages the "
        "example tree drives machine learning potentials (ANI-1xnr, "
        "AIMNet2, MACE, ...) whose checkpoints and GPU stack are not part of "
        "the checkout.  Only the committed output files are used."
        % (", ".join(record.get("maple_missing_dependencies") or
                      ["unknown"]),
           record.get("maple_import_error") or "no import error recorded"))

    # The probe and the record must not disagree, and in this environment the
    # probe *must* fail (see the docstring).
    if record["maple_importable"] != (record.get("maple_import_error") == ""):
        fail("the MAPLE probe is not self-consistent: importable=%s with "
             "error %r" % (record["maple_importable"],
                           record.get("maple_import_error")))
    if record["maple_importable"]:
        fail("the deep MAPLE probe now reports MAPLE usable; this checkout's "
             "validation record says MAPLE was not run, so the record, the "
             "reason text and docs/xtb-landscape-M7-validation.md have to be "
             "rewritten by hand before this check may pass again")
    if not record.get("maple_import_error"):
        fail("the MAPLE probe failed but recorded no diagnosis (empty "
             "error); the previous round shipped exactly that as a false "
             "green")
    if record["maple_executed"]:
        fail("maple_executed is true, but this mode runs no MAPLE job")

    record["our_topology"] = check_ts_irc(xtb, inp, os.path.join(rundir,
                                                                 "ts-irc"))
    covered = [task for task in fixtures if not fixtures[task].get("missing")]
    fixtures_sentence = ("Committed MAPLE fixtures were read for the task(s) "
                         "%s." % (", ".join(covered) if covered else "none"))
    if without_fixture:
        fixtures_sentence += ("  The example tree holds no committed fixture "
                              "for %s." % ", ".join(without_fixture))
    record["conclusion"] = (
        "topology comparison only: this suite's HCN/HNC IRC starts from a "
        "first order saddle point (1 imaginary mode) and ends in two distinct "
        "minima.  " + fixtures_sentence + "  A numerical comparison is not "
        "possible (MAPLE uses ML potentials, this suite GFN2, PLAN R-14/N9); "
        "MAPLE itself was not executed here and could not be, its dispatcher "
        "does not import in this environment.")
    with open(os.path.join(rundir, "maple.json"), "w") as handle:
        json.dump(record, handle, indent=1, sort_keys=True)
    print("maple: executed=%s importable=%s missing=%s topology=%s"
          % (record["maple_executed"], record["maple_importable"],
             record.get("maple_missing_dependencies"), record["our_topology"]))


def probe_maple(maple_dir):
    """Really probe MAPLE: import the dispatcher and instantiate an object.

    `import maple` only executes the shallow `__init__` and returns 0 even
    when nothing in the package can be used (validation round finding F6).
    The probe below is therefore the import a job would need plus the
    instantiation of the dispatcher, and it records the missing modules
    separately so the reason text can name them.
    """
    deep = ("from maple.function.dispatcher import Dispatcher; "
            "from maple.function.dispatcher.jobABC import JobABC; "
            "Dispatcher")
    shallow = subprocess.run([sys.executable, "-c", "import maple"],
                             cwd=maple_dir, capture_output=True, text=True)
    proc = subprocess.run([sys.executable, "-c", deep], cwd=maple_dir,
                          capture_output=True, text=True)
    missing = []
    for module in ("ase", "torch", "fairchem", "numpy"):
        check = subprocess.run(
            [sys.executable, "-c", "import %s" % module], cwd=maple_dir,
            capture_output=True, text=True)
        if check.returncode != 0:
            missing.append(module)
    error = ""
    if proc.returncode != 0:
        lines = [line for line in proc.stderr.strip().splitlines() if line]
        error = lines[-1] if lines else "the dispatcher import failed"
    return {
        "maple_toplevel_importable": shallow.returncode == 0,
        "maple_importable": proc.returncode == 0,
        "maple_import_probe": deep,
        "maple_import_error": error,
        "maple_missing_dependencies": missing,
    }


def count_extrema(energies):
    """Count local minima and maxima of a committed energy series."""
    minima = 0
    maxima = 0
    for i in range(1, len(energies) - 1):
        if energies[i] < energies[i - 1] and energies[i] < energies[i + 1]:
            minima += 1
        if energies[i] > energies[i - 1] and energies[i] > energies[i + 1]:
            maxima += 1
    return {"minima": minima, "maxima": maxima}


# The uniform top level of a `--ljson` document (format 2).  Every task, suite
# or native, writes exactly these keys in exactly this order (report 117).
CORE_KEYS = ["format", "task", "status", "status_name", "converged",
             "message", "n_atoms", "energy", "geometry", "provenance",
             "counters", "values", "flags", "history"]
#> Keys of the geometry block; the unit is stated explicitly (bohr).
GEOMETRY_KEYS = ["unit", "symbols", "coordinates"]
#> Keys of the provenance block.
PROVENANCE_KEYS = ["oxtb_version", "engine_version", "commit", "command_line",
                   "method", "charge", "spin", "input_file"]


def parse_document(path, expected_task):
    """The single parser: read any task document and return the core block.

    This one function is the whole point of the unified schema.  It asserts
    the same top level keys, the same types and the same core block for every
    task - suite or native - and returns the parsed document so a caller can
    make further assertions.  A task that omits or renames a core key,
    changes a core type, drops the geometry, states the wrong unit or writes a
    different `format` fails here, whichever task it is.
    """
    try:
        with open(path) as handle:
            doc = json.load(handle)
    except (OSError, ValueError) as error:
        fail("%s is not valid JSON: %s" % (path, error),
             read(path) if os.path.exists(path) else "")
    if list(doc) != CORE_KEYS:
        fail("%s has an unexpected key order: %s" % (path, list(doc)),
             read(path))
    if isinstance(doc["format"], bool) or not isinstance(doc["format"], int) \
            or doc["format"] != 2:
        fail("%s does not declare format 2: %r" % (path, doc["format"]))
    if doc["task"] != expected_task:
        fail("%s declares task %r, expected %r"
             % (path, doc["task"], expected_task))
    if not isinstance(doc["status"], int) or isinstance(doc["status"], bool):
        fail("%s status is not an integer" % path)
    if not isinstance(doc["status_name"], str):
        fail("%s status_name is not a string" % path)
    if not isinstance(doc["converged"], bool):
        fail("%s converged is not a boolean" % path)
    if not isinstance(doc["message"], str):
        fail("%s message is not a string" % path)
    if not isinstance(doc["n_atoms"], int) or isinstance(doc["n_atoms"], bool) \
            or doc["n_atoms"] <= 0:
        fail("%s n_atoms is not a positive integer" % path)
    if isinstance(doc["energy"], bool) \
            or not isinstance(doc["energy"], (int, float)):
        fail("%s energy is not a number" % path)
    geometry = doc["geometry"]
    if not isinstance(geometry, dict) or list(geometry) != GEOMETRY_KEYS:
        fail("%s geometry has an unexpected shape: %s" % (path, geometry))
    if geometry["unit"] != "bohr":
        fail("%s geometry unit is %r, expected 'bohr'"
             % (path, geometry["unit"]))
    symbols = geometry["symbols"]
    coordinates = geometry["coordinates"]
    if not isinstance(symbols, list) or len(symbols) != doc["n_atoms"]:
        fail("%s geometry has %s symbols for %d atoms"
             % (path, len(symbols), doc["n_atoms"]))
    if any(not isinstance(s, str) or not s for s in symbols):
        fail("%s geometry symbols are not non-empty strings" % path)
    if not isinstance(coordinates, list) or len(coordinates) != doc["n_atoms"]:
        fail("%s geometry has %s coordinate rows for %d atoms"
             % (path, len(coordinates), doc["n_atoms"]))
    for row in coordinates:
        if (not isinstance(row, list) or len(row) != 3 or
                any(isinstance(v, bool) or not isinstance(v, (int, float))
                    for v in row)):
            fail("%s geometry row is not an [x, y, z] triple: %s"
                 % (path, row))
    provenance = doc["provenance"]
    if not isinstance(provenance, dict) or list(provenance) != PROVENANCE_KEYS:
        fail("%s provenance has an unexpected shape: %s" % (path, provenance))
    for key in ("oxtb_version", "engine_version", "commit", "command_line",
                "method", "input_file"):
        if not isinstance(provenance[key], str):
            fail("%s provenance %s is not a string" % (path, key))
    if not provenance["oxtb_version"] or not provenance["engine_version"]:
        fail("%s provenance carries no version" % path)
    if not provenance["command_line"]:
        fail("%s provenance carries no command line" % path)
    for key in ("charge", "spin"):
        if isinstance(provenance[key], bool) \
                or not isinstance(provenance[key], int):
            fail("%s provenance %s is not an integer" % (path, key))
    for key in ("counters", "values", "flags", "history"):
        if not isinstance(doc[key], dict):
            fail("%s %s is not an object" % (path, key))
    return doc


def check_ljson_core(exe, inp, rundir, product):
    """The single-parser contract: every task parses with one parser.

    Runs a representative matrix - the six suite tasks, all eight `ts`
    methods, all four `irc` methods and the native run types - and feeds
    every document to `parse_document`, the one parser.  This is the
    falsifiable form of "one parser reads every task": removing the core
    block from any single task makes the parser reject that task's document,
    so the check goes red.  The assertions are structural (keys, types,
    shape, unit), not numerical; the numerical cross-checks live in the
    `ljson` and `*-counters` modes.
    """
    rundir = os.path.abspath(rundir)
    os.makedirs(rundir, exist_ok=True)
    product_file = os.path.join(rundir, "product.xyz")
    with open(product_file, "w") as handle:
        handle.write(product)
    md_ctrl = os.path.join(rundir, "md.ctrl")
    with open(md_ctrl, "w") as handle:
        handle.write("$md\n  time=0.02\n  temp=300.0\n$end\n")

    # (label, args, task) - the native run types first, then the suite ones
    # (the `ts`/`irc` method loops below add the rest).
    matrix = [
        ("native-sp", ["--sp"], "sp"),
        ("native-grad", ["--grad"], "grad"),
        ("native-opt", ["--opt", "--method", "ancopt"], "opt"),
        ("native-hess", ["--hess"], "hess"),
        ("native-ohess", ["--ohess", "--method", "ancopt"], "ohess"),
        ("native-md", ["--md", "--input", md_ctrl], "md"),
        ("suite-opt", ["--lopt", "--lopt-level", "loose"], "opt"),
        ("suite-scan", TASK_ARGS["scan"], "scan"),
        ("suite-vib", ["--lvib"], "vib"),
    ]
    for method, _ in TS_METHODS:
        args = ["--lts", "--lts-method", method, "--lts-max-steps", "20",
                "--lts-level", "superloose"]
        if method in ("neb", "cineb", "autoneb", "string", "gsm"):
            args += ["--lts-product", product_file]
        matrix.append(("ts-%s" % method, args, "ts"))
    for method, _ in IRC_METHODS:
        matrix.append(("irc-%s" % method,
                       ["--lirc", "--lirc-method", method,
                        "--lirc-max-steps", "4"], "irc"))

    tasks = set()
    for label, args, task in matrix:
        run_dir = os.path.join(rundir, label)
        status, output = run(exe, inp,
                             ["--gfn", "2", "--norestart", "--ljson"] + args,
                             run_dir)
        path = os.path.join(run_dir, "landscape_%s.json" % task)
        if not os.path.exists(path):
            fail("the %s run (status %d) wrote no landscape_%s.json"
                 % (label, status, task), output)
        parse_document(path, task)
        tasks.add(task)
    print("ljson-core: %d documents of %d tasks (%s) parse with one parser"
          % (len(matrix), len(tasks), ", ".join(sorted(tasks))))


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    mode = sys.argv[1]
    if mode == "ljson":
        if len(sys.argv) < 6:
            print("usage: ljson <exe> <input> <rundir> <task>")
            sys.exit(2)
        check_ljson(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif mode == "chain":
        check_chain(sys.argv[2], sys.argv[3], sys.argv[4],
                    sys.argv[5] if len(sys.argv) > 5 else None)
    elif mode == "chain-failure":
        check_chain_failure(sys.argv[2], sys.argv[3], sys.argv[4])
    elif mode == "ljson-methods":
        check_ljson_methods(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif mode == "ljson-core":
        check_ljson_core(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif mode == "ts-counters":
        check_ts_counters(sys.argv[2], sys.argv[3], sys.argv[4])
    elif mode == "ts-failure":
        check_ts_failure(sys.argv[2], sys.argv[3], sys.argv[4])
    elif mode == "ts-status":
        check_ts_status(sys.argv[2], sys.argv[3], sys.argv[4])
    elif mode == "scan-counters":
        check_scan_counters(sys.argv[2], sys.argv[3], sys.argv[4])
    elif mode == "layers":
        check_layers(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5],
                     sys.argv[6] if len(sys.argv) > 6 else None)
    elif mode == "maple":
        check_maple(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        print("unknown mode '%s'" % mode)
        sys.exit(2)


if __name__ == "__main__":
    main()
