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
"""CLI level checks of the landscape suite.

The CMake build drives these checks with small `cmake -P` scripts
(`test/unit/CMakeLists.txt`); a `cmake -P` script cannot inspect the output of a
test, so the checks live here and are invoked as test programs.  Every
mode mirrors one CMake test one-to-one:

  energy <exe> <input> <rundir>       -> xtb/landscape-opt-energy
  acc    <exe> <input> <rundir>       -> xtb/landscape-opt-acc
  hessian <exe> <input> <rundir>      -> xtb/landscape-opt-hessian
  scf    <exe> <input> <rundir>       -> xtb/landscape-scf-failure
  scanfail <exe> <input> <rundir>     -> xtb/landscape-scan-failed-point
  vib    <exe> <input> <rundir>       -> xtb/landscape-vib-cli
  vibnamespace <exe> <input> <rundir> -> xtb/landscape-vib-namespace
  optnamespace <exe> <input> <rundir> -> xtb/landscape-opt-namespace
  ljsondefault <exe> <input> <rundir> -> xtb/landscape-ljson-default
  vibetemp <exe> <input> <rundir>     -> xtb/landscape-vib-etemp
  enumcase <exe> <input> <rundir> <second>
                                     -> xtb/landscape-enum-case
  strictunknown <exe> <input> <rundir> -> xtb/landscape-strict-unknown
  flagvalue <exe> <input> <rundir>    -> xtb/landscape-flag-value
  nativelevel <exe> <input> <rundir>  -> xtb/landscape-native-level-value
  parallel <exe> <input> <rundir>     -> xtb/landscape-parallel-value
  positional <exe> <input> <rundir>   -> xtb/landscape-positional-file
  color  <exe> <input> <rundir>       -> xtb/landscape-color
  gfn    <exe> <input> <rundir>       -> xtb/landscape-gfn-value
  dockgfn <exe> <inputA> <rundir> <inputB>
                                     -> xtb/landscape-dock-gfn-value
  echo <exe> <input> <rundir> <case>
                                     -> xtb/landscape-fc-echo,
                                         xtb/landscape-lopt-engine-echo,
                                         xtb/landscape-lirc-default-echo,
                                         xtb/landscape-lts-escape-echo
  fix    <exe> <input> <rundir>       -> xtb/landscape-fix-rfo
  scanblock <exe> <input> <rundir>    -> xtb/landscape-scan-xcontrol
  saddlefix <exe> <input> <rundir>    -> xtb/landscape-lts-fix
  rejectfix <exe> <input> <rundir>    -> xtb/landscape-lirc-fix
    thermo <exe> <input> <rundir>       -> xtb/landscape-vib-thermo
  thermoknobs <exe> <input> <rundir> <saddle>
                                     -> xtb/landscape-vib-thermo-knobs
  ts     <exe> <input> <rundir>       -> xtb/landscape-ts-prfo
  tsetemp <exe> <input> <rundir>      -> xtb/landscape-ts-etemp
  dimer  <exe> <input> <rundir>       -> xtb/landscape-ts-dimer
  tsminimum <exe> <input> <rundir>    -> xtb/landscape-ts-minimum
  tsnamespace <exe> <input> <rundir>  -> xtb/landscape-ts-namespace
  maxsteps <exe> <input> <rundir>     -> xtb/landscape-scan-max-steps
  scanhess <exe> <input> <rundir>     -> xtb/landscape-scan-hessian
  autonebdepth <exe> <input> <rundir> <product>
                                     -> xtb/landscape-path-autoneb-depth-limit
  irc    <exe> <input> <rundir>       -> xtb/landscape-irc-cli
  ircpc  <exe> <input> <rundir>       -> xtb/landscape-irc-pc
  ircknobs <exe> <input> <rundir>     -> xtb/landscape-irc-knobs
  irctermination <exe> <input> <rundir> -> xtb/landscape-irc-termination
  ircswap <exe> <input> <rundir>     -> xtb/landscape-irc-swap
  ircosc <exe> <input> <rundir>      -> xtb/landscape-irc-osc-window
  blasthreads <exe> <input> <rundir>  -> xtb/landscape-blas-threads
  blasthreads-timing <exe> <input> <rundir> -> xtb/landscape-blas-threads-timing
  banner <exe> <input> <rundir>       -> xtb/landscape-banner
  enginebanner <exe> <input> <rundir> -> xtb/landscape-engine-banner
  programname <exe> <input> <rundir>  -> xtb/landscape-program-name
  nonfinite <exe> <input> <rundir>    -> xtb/landscape-nonfinite-value
  irfatal <exe> <input> <rundir>      -> xtb/landscape-ir-fatal-status
  helptopic <exe> <input> <rundir> <main.F90>
                                     -> xtb/landscape-help-topic
  helpcolour <exe> <input> <rundir>
                                     -> xtb/landscape-help-colour
  error  <exe> <input> <rundir> <pattern> <arg>...  -> the `--lopt`/`--lscan`
                                          error branches of the CMake build
  method <exe> <input> <rundir> <hcn> <nh3> <nh3product>
                                     -> xtb/landscape-method
  methodnative <exe> <input> <rundir> <big>
                                     -> xtb/landscape-method-native,
                                         xtb/landscape-method-ohess,
                                         xtb/landscape-method-opt-reject
  nativeoptsym <exe> <input> <rundir>
                                     -> xtb/native-opt-lbfgs-symbols
  methodtable <exe> <input> <rundir> <hcn> <nh3> <nh3product>
                                     -> xtb/landscape-method-table (all 18
                                         suite names and the four `auto` runs)

The stage 2b checks of the CLI redesign (section 6.3) run more than one flag
combination and compare the output files of two runs byte by byte, so they are
implemented here **once** and the CMake build drives this script as well (the
precedent is `landscape_e2e_check.py`, which CMake runs through
`landscape_e2e_py.cmake`); their modes are

  optknobs <exe> <input> <rundir>    -> xtb/landscape-opt-knobs
  scanknobs <exe> <input> <rundir>   -> xtb/landscape-scan-knobs
  hessstepalias <exe> <input> <rundir> <hcn>
                                     -> xtb/landscape-hess-step-alias
  restartknobs <exe> <input> <rundir>
                                     -> xtb/landscape-restart-knobs
  hessused <exe> <input> <rundir> <hcn> <nh3> <nh3product>
                                     -> xtb/landscape-hess-step-used
  lvibrestart <exe> <input> <rundir> -> xtb/landscape-lvib-restart-warm

The stage 3 checks of the CLI redesign run `--opt`/`--ohess` as landscape
suite tasks and are implemented here once for both builds; their modes are

  optalias <exe> <input> <rundir>    -> xtb/landscape-opt-alias-once
  optnative <exe> <input> <rundir>   -> xtb/landscape-opt-native-exit
  ohesssuite <exe> <input> <rundir>  -> xtb/landscape-ohess-suite-once
  ohessnative <exe> <input> <rundir> -> xtb/landscape-ohess-native-hess
  ohessbridge <exe> <input> <rundir> -> xtb/landscape-ohess-fres-bridge
  ohesslevelaxes <exe> <input> <rundir>
                                     -> xtb/landscape-ohess-level-axes
  optguard <exe> <input> <rundir> <scanblock|pbc|xcontrol>
                                     -> xtb/landscape-opt-scanblock,
                                        xtb/landscape-opt-pbc,
                                        xtb/landscape-opt-xcontrol
  optxtbopt <exe> <input> <rundir>   -> xtb/landscape-opt-xtbopt-once
  optlevel / cycles / levelaxes <exe> <input> <rundir>
                                     -> xtb/landscape-opt-level-bridge,
                                        xtb/landscape-cycles-map,
                                        xtb/landscape-level-axes

The environment (XTBPATH) is inherited from the test environment, exactly as in
the rest of the suite.
"""

import json
import os
import re
import shutil
import subprocess
import sys
import time


def decode(raw):
    """Decode captured output without ever raising.

    xtb writes Fortran record markers into its formatted output in rare
    cases (measured with the 501 atom GFN-FF fixture of `check_method_native`:
    1/8/192 invalid bytes depending on the run - a pre-existing upstream
    defect, `docs/_simplify-recon/28-stage2-fix-report.md` T2).  Letting
    `UnicodeDecodeError` escape here would make every consumer of this module
    red for a reason that has nothing to do with the assertion under test,
    so the captured bytes are decoded with `errors='replace'`.
    """
    if isinstance(raw, str):
        return raw
    return raw.decode("utf-8", errors="replace")


def run(exe, inp, args, rundir):
    """Run xtb and return (status, combined output)."""
    os.makedirs(rundir, exist_ok=True)
    proc = subprocess.run(
        [exe, inp] + args,
        cwd=rundir,
        capture_output=True,
        text=False,
    )
    return proc.returncode, decode(proc.stdout) + decode(proc.stderr)


def run_split(exe, inp, args, rundir):
    """Run xtb and return (status, stdout, stderr), the streams kept apart.

    The termination lines of `src/mctc/error.f90` go to **stderr** while the
    banner and the results go to stdout, so a check of the program name in
    those lines has to read stderr alone; `run` above deliberately merges the
    two streams.
    """
    os.makedirs(rundir, exist_ok=True)
    proc = subprocess.run(
        [exe, inp] + args,
        cwd=rundir,
        capture_output=True,
        text=False,
    )
    return proc.returncode, decode(proc.stdout), decode(proc.stderr)


def fail(message, output=""):
    print("FAIL: " + message)
    if output:
        print(output)
    sys.exit(1)


def check_energy(exe, inp, rundir):
    """The TOTAL ENERGY line must be the energy written to xtbopt.xyz."""
    status, output = run(exe, inp, ["--gfn", "2", "--lopt", "--norestart"], rundir)
    if status != 0:
        fail("xtb --lopt terminated with status %d" % status, output)
    match = re.search(r"TOTAL ENERGY +(-?[0-9]+\.[0-9]+)", output)
    if not match:
        fail("no TOTAL ENERGY in the output", output)
    stdout_energy = match.group(1)
    with open(os.path.join(rundir, "xtbopt.xyz")) as handle:
        xyz = handle.read()
    match = re.search(r"energy: +(-?[0-9]+\.[0-9]+)", xyz)
    if not match:
        fail("no energy in xtbopt.xyz", xyz)
    if stdout_energy != match.group(1):
        fail(
            "TOTAL ENERGY (%s) differs from xtbopt.xyz (%s)"
            % (stdout_energy, match.group(1))
        )


def opt_result(exe, inp, args, rundir):
    """Run one minimisation and return (reported accuracy, xtbopt energy line)."""
    xyzfile = os.path.join(rundir, "xtbopt.xyz")
    if os.path.exists(xyzfile):
        os.remove(xyzfile)
    status, output = run(
        exe, inp, ["--gfn", "2", "--lopt", "--norestart"] + args, rundir
    )
    if status != 0:
        fail("xtb --lopt %s failed with status %d" % (" ".join(args), status), output)
    with open(xyzfile) as handle:
        match = re.search(r"energy: [^\n]*", handle.read())
    result = match.group(0) if match else ""
    match = re.search(r"landscape engine: SCC accuracy +([0-9.]+)", output)
    return (match.group(1) if match else "", result, output)


def check_acc(exe, inp, rundir):
    """`--acc` must reach the engine control and change the result."""
    default_acc, default_result, output = opt_result(exe, inp, [], rundir)
    acc_acc, acc_result, _ = opt_result(exe, inp, ["--acc", "0.05"], rundir)
    if default_acc != "1.000000":
        fail("the default run reported accuracy %r" % default_acc, output)
    if acc_acc != "0.050000":
        fail("--acc 0.05 did not reach the engine control (%r)" % acc_acc)
    if default_result == acc_result:
        fail("--acc 0.05 had no effect on the --lopt result (%s)" % default_result)


def check_hessian(exe, inp, rundir):
    """`--lopt-hessian` must reach the engine and stay comparable."""
    cases = [
        ("model_lbfgs", ["--lopt-engine", "lbfgs", "--lopt-hessian", "model"], "model"),
        ("num_lbfgs", ["--lopt-engine", "lbfgs", "--lopt-hessian", "numerical"], "numerical"),
        ("model_rfo", ["--lopt-engine", "rfo", "--lopt-hessian", "model"], "model"),
        ("num_rfo", ["--lopt-engine", "rfo", "--lopt-hessian", "numerical"], "numerical"),
    ]
    energies = {}
    for name, args, expected_source in cases:
        status, output = run(
            exe, inp, ["--gfn", "2", "--lopt", "--norestart"] + args, rundir
        )
        if status != 0:
            fail("xtb --lopt %s failed with status %d" % (" ".join(args), status), output)
        if "landscape: converged" not in output:
            fail("xtb --lopt %s did not converge" % " ".join(args), output)
        match = re.search(r"Hessian source: +([a-z]+)", output)
        found = match.group(1) if match else ""
        if found != expected_source:
            fail(
                "%s reported Hessian source %r, expected %r"
                % (name, found, expected_source)
            )
        with open(os.path.join(rundir, "xtbopt.xyz")) as handle:
            xyz = handle.read()
        match = re.search(r"energy: +(-?[0-9]+\.[0-9]{3})", xyz)
        if not match:
            fail("the %s run reported no energy" % name, xyz)
        energies[name] = match.group(1)
    for first, second in (("model_lbfgs", "num_lbfgs"), ("model_rfo", "num_rfo")):
        if energies[first] != energies[second]:
            fail(
                "%s and %s converge to different energies (%s vs %s)"
                % (first, second, energies[first], energies[second])
            )


def check_scf(exe, inp, rundir):
    """An unconverged SCC must abort the run with a non-zero status."""
    os.makedirs(rundir, exist_ok=True)
    xcontrol = os.path.join(rundir, "scf_fail.inp")
    with open(xcontrol, "w") as handle:
        handle.write("$scc\n maxiterations=1\n$end\n")
    status, output = run(
        exe,
        inp,
        ["--gfn", "2", "--lopt", "--input", xcontrol, "--norestart"],
        rundir,
    )
    if status == 0:
        fail("the unconverged SCC did not fail the run", output)
    if "landscape: converged" in output:
        fail("a converged landscape run was reported", output)


def check_scanfail(exe, inp, rundir):
    """A scan point without an energy must be reported and not written.

    A rigid distance scan with a 5 angstrom step puts the second point at
    5.96 angstrom, where the SCC of the single point fails.  The run has to
    end with a non-zero status and a readable diagnostic, and the scan log
    must not contain a frame with a fabricated zero energy (nor a frame per
    requested point, since a point without an energy has none).
    """
    args = [
        "--gfn", "2", "--lscan", "--lscan-mode", "rigid",
        "--lscan-coord", "distance", "--lscan-atoms", "1,2",
        "--lscan-step", "5.0", "--lscan-nstep", "2",
        "--namespace", "landscape3",
    ]
    status, output = run(exe, inp, args, rundir)
    if status == 0:
        fail("the failed scan points did not fail the run", output)
    if "scan points failed" not in output:
        fail("the failed scan points were not reported", output)
    with open(os.path.join(rundir, "landscape3.xtbscan.log")) as handle:
        scanlog = handle.read()
    if "energy: 0.000000000000" in scanlog:
        fail("a failed scan point was written with a zero energy", scanlog)
    frames = len(re.findall(r"energy:", scanlog))
    if frames >= 5:
        fail("the scan log holds %d frames, expected fewer than 5" % frames)



def check_vib(exe, inp, rundir):
    """`--lvib` must write the `--hess` files and reproduce them from disk.

    Mirrors the CMake script of `xtb/landscape-vib-cli`: the four output files
    have to exist, the task energy has to be the reported TOTAL ENERGY (the
    `res` write-back), and a re-run with `--lvib-hessian file` must reproduce
    the frequencies of the direct numerical Hessian.
    """
    status, output = run(exe, inp, ["--gfn", "2", "--lvib", "--norestart"], rundir)
    if status != 0:
        fail("xtb --lvib terminated with status %d" % status, output)
    if "landscape: harmonic frequencies" not in output:
        fail("the frequency header was not printed", output)
    if not re.search(r"Hessian source: +numerical", output):
        fail("the Hessian source was not reported", output)
    if "landscape task terminated: frequency analysis completed" not in output:
        fail("the frequency task did not complete", output)
    if "modes: 3 internal, 6 projected zero modes" not in output:
        fail("the mode count of water is wrong", output)
    if not re.search(r"linear molecule: +F", output):
        fail("water was classified as linear", output)
    for name in ("hessian", "vibspectrum", "vib_normal_modes", "g98.out"):
        path = os.path.join(rundir, name)
        if not os.path.exists(path):
            fail("%s was not written" % name)
        if os.path.getsize(path) == 0:
            fail("%s is empty" % name)
    match = re.search(r"energy: +(-?[0-9]+[.][0-9]+)", output)
    task_energy = match.group(1) if match else ""
    match = re.search(r"TOTAL ENERGY +(-?[0-9]+[.][0-9]+)", output)
    total_energy = match.group(1) if match else ""
    if not task_energy:
        fail("the frequency task reported no energy", output)
    if task_energy != total_energy:
        fail(
            "the task energy (%s) is not the reported TOTAL ENERGY (%s)"
            % (task_energy, total_energy)
        )
    match = re.search(r"eigval : +([-0-9. ]+)", output)
    vib_freqs = " ".join(match.group(1).split()) if match else ""
    # the direct run builds the numerical Hessian from 6N displaced single
    # points plus the gradient of the input geometry: 1 + 6*3 = 19 for water
    match = re.search(r"landscape engine: energy evaluations +(\d+)", output)
    if not match:
        fail("the direct run reports no evaluation counter", output)
    direct_eval = int(match.group(1))
    if direct_eval != 19:
        fail("the direct --lvib run used %d evaluations, expected 1+6N = 19"
             % direct_eval, output)
    status, file_output = run(
        exe,
        inp,
        ["--gfn", "2", "--lvib", "--lvib-hessian", "file",
         "--lvib-hessian-file", "hessian", "--norestart"],
        rundir,
    )
    if status != 0:
        fail("reading the Hessian file failed with status %d" % status, file_output)
    if not re.search(r"Hessian file: +hessian", file_output):
        fail("the file provider did not report its file", file_output)
    # the file provider must *save* the displaced single points: reading the
    # file needs the gradient of the input geometry only (19 -> 1 for water),
    # otherwise a provider that silently recomputed the Hessian would be green
    match = re.search(r"landscape engine: energy evaluations +(\d+)", file_output)
    if not match:
        fail("the file run reports no evaluation counter", file_output)
    file_eval = int(match.group(1))
    if file_eval != 1:
        fail("--lvib-hessian file used %d evaluations, expected 1; the "
             "provider did not reuse the file" % file_eval, file_output)
    if direct_eval <= file_eval:
        fail("the Hessian file did not save evaluations (%d vs %d)"
             % (direct_eval, file_eval))
    match = re.search(r"eigval : +([-0-9. ]+)", file_output)
    file_freqs = " ".join(match.group(1).split()) if match else ""
    if file_freqs != vib_freqs:
        fail(
            "the Hessian file reproduces different frequencies (%s vs %s)"
            % (file_freqs, vib_freqs)
        )


def check_vib_namespace(exe, inp, rundir):
    """The `hessian` file of `--lvib` carries the namespace exactly once.

    Mirrors the CMake script of `xtb/landscape-vib-namespace`: with
    `--namespace ns` the task writes `ns.hessian` (the name the `--hess` path
    produces as well), not `ns.ns.hessian`, every other output carries the
    prefix too (`ns.vibspectrum`, `ns.vib_normal_modes`, `ns.g98.out`), and the
    re-run with `--lvib-hessian file --lvib-hessian-file hessian --namespace
    ns` finds and reproduces the file.  The rundir is removed at entry, because
    the check asserts the *absence* of the unprefixed output names.
    """
    rundir = os.path.abspath(rundir)
    shutil.rmtree(rundir, ignore_errors=True)
    status, output = run(
        exe, inp, ["--gfn", "2", "--lvib", "--norestart", "--namespace", "ns"], rundir
    )
    if status != 0:
        fail("xtb --lvib --namespace ns terminated with status %d" % status, output)
    if not os.path.exists(os.path.join(rundir, "ns.hessian")):
        fail("the namespaced hessian file was not written", output)
    if os.path.exists(os.path.join(rundir, "ns.ns.hessian")):
        fail("the hessian file name carries the namespace twice")
    # every other output of the task carries the namespace as well
    for name in ("hessian", "vibspectrum", "vib_normal_modes", "g98.out"):
        if not os.path.exists(os.path.join(rundir, "ns." + name)):
            fail("--namespace ns wrote no ns.%s" % name, output)
        if os.path.exists(os.path.join(rundir, name)):
            fail("--namespace ns also wrote the unprefixed %s" % name)
    match = re.search(r"eigval : +([-0-9. ]+)", output)
    vib_freqs = " ".join(match.group(1).split()) if match else ""
    status, file_output = run(
        exe,
        inp,
        ["--gfn", "2", "--lvib", "--lvib-hessian", "file",
         "--lvib-hessian-file", "hessian", "--norestart", "--namespace", "ns"],
        rundir,
    )
    if status != 0:
        fail("reading the namespaced Hessian file failed with status %d" % status, file_output)
    match = re.search(r"eigval : +([-0-9. ]+)", file_output)
    file_freqs = " ".join(match.group(1).split()) if match else ""
    if file_freqs != vib_freqs:
        fail(
            "the namespaced Hessian file reproduces different frequencies (%s vs %s)"
            % (file_freqs, vib_freqs)
        )


def check_vib_etemp(exe, inp, rundir):
    """The CLI `--etemp` must reach the landscape engine.

    Mirrors the CMake script of `xtb/landscape-vib-etemp`: the engine control
    is the single source of truth for the electronic temperature, so a missing
    `main.F90` sync would silently run every landscape task at 300 K.  The
    stretched H2 input shifts its only mode by >30 % between 300 K and 5000 K
    (report 21 section 3.1), so the two runs have to report different effective
    temperatures and different frequencies.
    """
    status, cold = run(
        exe, inp, ["--gfn", "2", "--lvib", "--norestart", "--etemp", "300"], rundir
    )
    if status != 0:
        fail("xtb --lvib --etemp 300 terminated with status %d" % status, cold)
    if not re.search(r"landscape engine: etemp +300[.]000 K", cold):
        fail("the effective electronic temperature 300 K was not reported", cold)
    status, hot = run(
        exe, inp, ["--gfn", "2", "--lvib", "--norestart", "--etemp", "5000"], rundir
    )
    if status != 0:
        fail("xtb --lvib --etemp 5000 terminated with status %d" % status, hot)
    if not re.search(r"landscape engine: etemp +5000[.]000 K", hot):
        fail("--etemp 5000 did not reach the engine control (missing CLI sync)", hot)
    match = re.search(r"eigval : +([-0-9. ]+)", cold)
    cold_freqs = " ".join(match.group(1).split()) if match else ""
    match = re.search(r"eigval : +([-0-9. ]+)", hot)
    hot_freqs = " ".join(match.group(1).split()) if match else ""
    if not cold_freqs:
        fail("the 300 K run reported no frequencies", cold)
    if not hot_freqs:
        fail("the 5000 K run reported no frequencies", hot)
    if cold_freqs == hot_freqs:
        fail(
            "the electronic temperature did not change the frequencies (%s)" % cold_freqs
        )


def check_fix(exe, inp, rundir):
    """`$fix` must hold the frozen atoms for every minimiser.

    Mirrors `xtb/landscape-fix-rfo`.  A first order direction is built from
    the gradient, and the calculator zeroes the gradient of a fixed atom
    (`src/xtb/calculator.f90:304-306`), so SD/CG/L-BFGS/BFGS/GDIIS keep it in
    place by construction.  The RFO step is expanded in the eigenvectors of the
    Hessian, which mix the frozen and the free components, so it moved the
    frozen oxygen of this fixture by 0.0098/0.0109 angstrom in x/z before the
    direction was projected (docs/_simplify-recon/04-*.md).  `--lopt-hessian
    model` seeds BFGS with the Lindh-D2 matrix as well, so that combination is
    covered here too.  The tolerance is
    the xyz write/read round trip (the writer prints 14 decimals), not a
    measured spread: a projected direction leaves the coordinate unchanged.
    """
    xcontrol = os.path.join(rundir, "landscape_fix.inp")
    os.makedirs(rundir, exist_ok=True)
    rundir = os.path.abspath(rundir)
    xcontrol = os.path.abspath(xcontrol)
    with open(xcontrol, "w") as handle:
        handle.write("$fix\n   atoms: 1\n$end\n")
    geometry = os.path.join(rundir, "xtbopt.xyz")
    # The first order engines are only *provably* decoupled when their seed is
    # decoupled as well: `--lopt-hessian model` seeds BFGS with the Lindh-D2
    # curvature matrix (`optimize/driver.f90` -> `landscape_bfgs%set_hessian`),
    # whose frozen and free blocks are coupled, and the BFGS direction is
    # expanded in the eigenvectors of that matrix exactly like the RFO one
    # (`optimize/bfgs.f90`, `d = -sum_i (V_i.g)/max(w_i, hlow) V_i`).  The
    # projection is therefore not a no-op for `bfgs model`/`rfo model` either,
    # so the two combinations are covered here (finding P2-3 of
    # `05-w0-independent-review.md`).
    cases = [("sd", "sd", []), ("cg", "cg", []), ("lbfgs", "lbfgs", []),
             ("bfgs", "bfgs", []), ("diis", "diis", []), ("rfo", "rfo", []),
             ("bfgs model", "bfgs", ["--lopt-hessian", "model"]),
             ("rfo model", "rfo", ["--lopt-hessian", "model"])]
    for label, engine, extra in cases:
        if os.path.exists(geometry):
            os.remove(geometry)
        status, output = run(
            exe,
            inp,
            ["--gfn", "2", "--lopt", "--lopt-engine", engine,
             "--input", xcontrol, "--norestart"] + extra,
            rundir,
        )
        if status != 0:
            fail("xtb --lopt-engine %s terminated with status %d" % (label, status),
                 output)
        if "landscape: converged" not in output:
            fail("the --lopt-engine %s run did not converge" % label, output)
        if not os.path.exists(geometry):
            fail("the --lopt-engine %s run wrote no xtbopt.xyz" % label, output)
        with open(geometry) as handle:
            lines = handle.read().splitlines()
        fields = lines[2].split()
        if fields[0] != "O":
            fail("the first output atom of %s is not the frozen oxygen" % label,
                 lines[2])
        shift = max(abs(float(value)) for value in fields[1:4])
        if shift > 1.0e-10:
            fail(
                "--lopt-engine %s moved the frozen atom by %.5e angstrom (%s)"
                % (label, shift, " ".join(fields[1:4])),
                output,
            )
        print("fix projection: %-11s frozen atom displacement %.5e angstrom"
              % (label, shift))


def check_scanblock(exe, inp, rundir):
    """A landscape run type must reject the xcontrol `$scan` block.

    Mirrors `xtb/landscape-scan-xcontrol`.  `$scan` belongs to the ANCopt
    family - `relaxed_scan` runs it only on the `anyopt` path, and since stage
    3 of the CLI redesign the redeployed `--opt` is a landscape run type and no
    longer walks that path - and no landscape run type executes it.  Before
    the guard, `--lopt --input <file with $scan>` returned status 0 after a
    plain minimisation and wrote no `xtbscan.log`: the caller asked for a scan
    and silently got a different task.  The check pins the non-zero status, the
    actionable message and the absence of the scan log for all five run types,
    and that `--opt` now refuses the block as loudly as they do while the
    native exit (`--opt --method ancopt`) still executes it.
    """
    os.makedirs(rundir, exist_ok=True)
    rundir = os.path.abspath(rundir)
    xcontrol = os.path.join(rundir, "landscape_scan.inp")
    with open(xcontrol, "w") as handle:
        handle.write("$constrain\n  distance: 1, 2, auto\n"
                     "$scan\n  1: 0.9, 1.1, 3\n$end\n")
    scanlog = os.path.join(rundir, "xtbscan.log")
    for task in ("lopt", "lscan", "lts", "lirc", "lvib"):
        if task == "lscan":
            # a *valid* `--lscan` command line: the `$scan` block is rejected
            # even though the run type itself is fine
            extra = ["--lscan-coord", "distance", "--lscan-atoms", "1,2",
                     "--lscan-step", "0.05", "--lscan-nstep", "2"]
        else:
            extra = []
        if os.path.exists(scanlog):
            os.remove(scanlog)
        status, output = run(
            exe,
            inp,
            ["--gfn", "2", "--" + task] + extra
            + ["--input", xcontrol, "--norestart"],
            rundir,
        )
        if status == 0:
            fail("xtb --%s accepted the xcontrol $scan block" % task, output)
        if "does not execute the xcontrol $scan block" not in output:
            fail("xtb --%s did not explain the rejected $scan block" % task, output)
        if os.path.exists(scanlog):
            fail("xtb --%s wrote xtbscan.log although it refused the $scan block"
                 % task, output)
        print("scan block: --%s rejected the $scan block with status %d"
              % (task, status))
    # the redeployed `--opt` refuses the block as well (CLI redesign section
    # 3.8(a)); it is no longer the flag that executes it
    optdir = os.path.join(rundir, "opt")
    if os.path.exists(os.path.join(optdir, "xtbscan.log")):
        os.remove(os.path.join(optdir, "xtbscan.log"))
    status, output = run(
        exe, inp, ["--gfn", "2", "--opt", "--input", xcontrol, "--norestart"],
        optdir,
    )
    if status == 0:
        fail("xtb --opt accepted the xcontrol $scan block", output)
    # the guard message itself: the `[NOTICE]` line of the same run mentions
    # the native exit as well, so the bare substring proves nothing
    if "does not execute the xcontrol $scan block: use --opt --method ancopt" \
            not in output:
        fail("--opt did not refuse the $scan block through the guard", output)
    if os.path.exists(os.path.join(optdir, "xtbscan.log")):
        fail("--opt wrote xtbscan.log although it refused the $scan block", output)
    # ... and the native exit still executes it
    nativedir = os.path.join(rundir, "opt-native")
    status, output = run(
        exe, inp,
        ["--gfn", "2", "--opt", "--method", "ancopt", "--input", xcontrol,
         "--norestart"],
        nativedir,
    )
    if status != 0:
        fail("xtb --opt --method ancopt with a $scan block terminated with "
             "status %d" % status, output)
    if "RELAXED SCAN" not in output:
        fail("--opt --method ancopt did not run the $scan block", output)
    if not os.path.exists(os.path.join(nativedir, "xtbscan.log")):
        fail("--opt --method ancopt did not write xtbscan.log", output)


def check_strict_unknown(exe, inp, rundir):
    """`--strict` covers the parser, and an unknown flag stays a warning.

    Mirrors `xtb/landscape-strict-unknown` (CLI redesign section 17.5 5a).  The
    unknown-flag message is an `env%warning`; it is only fatal while
    `env%strict` is set, and that used to happen *after* the parsing loop, so
    `--strict --zzzz` exited 0.  The second half pins BC-1: without `--strict`
    the same flag only warns and the single point still runs.  The third run
    puts the unknown flag **before** `--strict`: the fatal decision is taken
    after the loop, so the order must not matter (review finding M5 of
    `20-stage1-review.md`; the position independent scheme is the reason the
    check moved to the parser tail).
    """
    os.makedirs(rundir, exist_ok=True)
    status, output = run(
        exe, inp, ["--strict", "--zzzz", "--gfn", "2", "--sp", "--norestart"],
        rundir,
    )
    if status == 0:
        fail("--strict --zzzz was accepted", output)
    if "Unknown option '--zzzz' provided" not in output:
        fail("the unknown flag was not reported", output)
    status, output = run(
        exe, inp, ["--zzzz", "--strict", "--gfn", "2", "--sp", "--norestart"],
        rundir,
    )
    if status == 0:
        fail("--zzzz --strict was accepted", output)
    if "Unknown option '--zzzz' provided" not in output:
        fail("the unknown flag before --strict was not reported", output)
    status, output = run(
        exe, inp, ["--zzzz", "--gfn", "2", "--sp", "--norestart"], rundir
    )
    if status != 0:
        fail("an unknown flag without --strict failed the run", output)
    if "Unknown option '--zzzz' provided" not in output:
        fail("the unknown flag warning disappeared", output)


def check_native_level(exe, inp, rundir):
    """An unknown native level keyword is reported, laxly and strictly.

    Mirrors `xtb/landscape-native-level-value`.  `--bhess`/`--metaopt` and
    the native `--opt`/`--ohess --method <engine>` hand their level to the
    upstream `optlevel2int`, whose `case default` maps any keyword it does not
    know to `normal`: `--bhess superloose` (a suite-only level) and
    `--bhess zzqq` used to run the default level without a message.  The check
    pins both directions of the lax/strict split of the parser tail: a
    documented keyword stays silent, an unknown keyword is named in a warning,
    and `--strict` turns the same keyword into a fatal error.  The xcontrol
    `$opt` route is pinned too, because it reaches the same parser through a
    second dispatcher.
    """
    os.makedirs(rundir, exist_ok=True)

    def level(args, label):
        directory = os.path.join(rundir, label)
        shutil.rmtree(directory, ignore_errors=True)
        return run(exe, inp, args + ["--norestart"], directory)

    # (1) the documented native keywords run without a report
    for flag in ("--bhess", "--metaopt"):
        for value in ("normal", "tight", "verytight", "lax"):
            status, output = level(["--gfn", "2", flag, value],
                                   "%s-%s" % (flag[2:], value))
            if status != 0:
                fail("%s %s failed with status %d" % (flag, value, status),
                     output)
            if "Unknown %s level" % flag in output:
                fail("%s %s was reported as unknown" % (flag, value), output)

    # (2) a keyword the native table does not know is reported and the run
    #     continues: the suite-only `superloose`, a misspelling, and
    #     `TIGHT`, because the upstream table is case sensitive
    for flag in ("--bhess", "--metaopt"):
        for value in ("superloose", "zzqq", "TIGHT"):
            status, output = level(["--gfn", "2", flag, value],
                                   "%s-bad-%s" % (flag[2:], value))
            if status != 0:
                fail("lax %s %s failed with status %d" % (flag, value, status),
                     output)
            if ("Unknown %s level '" % flag) not in output:
                fail("%s %s was not reported" % (flag, value), output)

    # (3) the same keyword is fatal with --strict, in either flag position
    for args, label in ((
            ["--gfn", "2", "--bhess", "superloose", "--strict"],
            "bhess-strict"), (
            ["--metaopt", "zzqq", "--gfn", "2", "--strict"],
            "metaopt-strict-before")):
        status, output = level(args, label)
        if status == 0:
            fail("--strict accepted %r" % " ".join(args), output)
        if "level '" not in output:
            fail("--strict did not report %r" % " ".join(args), output)

    # (4) the native `--opt`/`--ohess` keep the same oracle; the level has to
    #     follow the flag directly, the parser reads the next token
    for flag in ("--opt", "--ohess"):
        status, output = level(["--gfn", "2", flag, "superloose",
                                "--method", "ancopt"], "%s-native" % flag[2:])
        if status != 0:
            fail("%s superloose --method ancopt failed with status %d"
                 % (flag, status), output)
        if "Unknown --opt level 'superloose'" not in output:
            fail("%s superloose --method ancopt was not reported" % flag,
                 output)
        status, output = level(["--gfn", "2", flag, "normal", "--method",
                                "ancopt"], "%s-native-normal" % flag[2:])
        if status != 0:
            fail("%s normal --method ancopt failed with status %d"
                 % (flag, status), output)
        if "Unknown --opt level" in output:
            fail("%s normal --method ancopt was reported" % flag, output)

    # (5) the xcontrol `$opt` block reaches the same parser
    control = os.path.join(rundir, "native-level.inp")
    with open(control, "w") as handle:
        handle.write("$opt\n   optlevel=superloose\n")
    status, output = level(["--gfn", "2", "--opt", "--method", "ancopt",
                            "--input", control], "xcontrol")
    if status != 0:
        fail("the xcontrol $opt optlevel=superloose run failed with status %d"
             % status, output)
    if "Unknown optimization level 'superloose'" not in output:
        fail("the xcontrol $opt optlevel=superloose was not reported", output)
    print("native level: unknown keywords reported on every native route")


def parallel_team(output):
    """The OpenMP team size the run reported, as an integer."""
    match = re.search(r"omp threads +: +([0-9]+)", output)
    if not match:
        fail("no 'omp threads' line in the output", output)
    return int(match.group(1))


def total_energy(output):
    """The printed total energy in hartree, as the string of the output."""
    match = re.search(r"TOTAL ENERGY +(-?[0-9]+\.[0-9]+)", output)
    if not match:
        fail("no TOTAL ENERGY in the output", output)
    return match.group(1)


def check_parallel(exe, inp, rundir):
    """`-P`/`--parallel` validates its value instead of passing it through.

    Mirrors `xtb/landscape-parallel-value`.  The flag is one
    `omp_set_num_threads` call with no range check: `-P 0` and `-P -2` became
    one thread and `-P abc` the default team, both silently, `-P` as the last
    argument was consumed and dropped, and `-P 1000` really created a thousand
    threads on an eight processor machine (measured, report 189).  The check
    pins the project's lax/strict shape: a value that is not a positive integer
    is reported and the run keeps the old behaviour, or is fatal with
    `--strict`, decided in the parser tail so that the position of `--strict`
    on the command line does not matter.  A value above the processor count is
    a third case: it is reported too, but the value is applied as given and
    never clamped - a silent clamp would be the defect class being fixed.

    The accepted path has to stay unchanged: `-P 2` must still set the team to
    2, print nothing about the flag and give the same energy as a run without
    it.  The `-P 1000` half uses GFN-FF, whose shorter run keeps the 1000
    thread team cheap while still creating it.
    """
    os.makedirs(rundir, exist_ok=True)

    def single_point(args, label):
        directory = os.path.join(rundir, label)
        shutil.rmtree(directory, ignore_errors=True)
        return run(exe, inp, args, directory)

    # (1) the accepted path: `-P 2` and `--parallel 2` are one flag, they set
    #     the team, they are silent and they do not move the energy
    energies = {}
    teams = {}
    for label, args in (
            ("no-flag", ["--gfn", "2", "--sp", "--norestart"]),
            ("short", ["-P", "2", "--gfn", "2", "--sp", "--norestart"]),
            ("long", ["--parallel", "2", "--gfn", "2", "--sp",
                       "--norestart"])):
        status, output = single_point(args, label)
        if status != 0:
            fail("the %s run failed with status %d" % (label, status), output)
        if "-P/--parallel" in output:
            fail("the valid value of the %s run was reported" % label, output)
        energies[label] = total_energy(output)
        teams[label] = parallel_team(output)
    for label in ("short", "long"):
        if teams[label] != 2:
            fail("-P 2 did not set the OpenMP team: the %s run reports %d"
                 % (label, teams[label]))
        if energies[label] != energies["no-flag"]:
            fail("-P 2 moved the total energy: %s vs %s without the flag"
                 % (energies[label], energies["no-flag"]))
    print("parallel: -P 2 and --parallel 2 set a team of 2, silently")

    # (2) a value that is not a positive integer is reported and the run keeps
    #     the behaviour the flag had before the validation: a non-positive
    #     value reaches the runtime, which maps it to a single thread
    for value, label, expected_team in (("0", "zero", 1),
                                        ("-2", "negative", None),
                                        ("abc", "text", None)):
        status, output = single_point(
            ["-P", value, "--gfn", "2", "--sp", "--norestart"], label)
        if status != 0:
            fail("the lax run -P %s failed with status %d" % (value, status),
                 output)
        if ("-P/--parallel needs a positive integer, not '%s'" % value) \
                not in output:
            fail("lax -P %s was not reported" % value, output)
        if expected_team is not None and parallel_team(output) != expected_team:
            fail("lax -P %s did not keep the old single thread team" % value,
                 output)
    print("parallel: 0, -2 and abc are reported, the run continues")

    # (3) `-P` as the last argument used to be consumed and ignored
    status, output = single_point(
        ["--gfn", "2", "--sp", "--norestart", "-P"], "missing")
    if status != 0:
        fail("'-P' as the last argument failed with status %d" % status, output)
    if "-P/--parallel needs a positive integer, but no value was given" \
            not in output:
        fail("'-P' as the last argument was not reported", output)
    print("parallel: a missing value is reported instead of guessed")

    # (4) `--strict` turns the report into a fatal error, in either flag
    #     position, and covers the non-positive value as well
    for label, args in (
            ("strict-before", ["--strict", "-P", "abc", "--gfn", "2",
                                "--sp", "--norestart"]),
            ("strict-after", ["-P", "abc", "--strict", "--gfn", "2",
                               "--sp", "--norestart"]),
            ("strict-zero", ["--strict", "-P", "0", "--gfn", "2", "--sp",
                              "--norestart"])):
        status, output = single_point(args, label)
        if status == 0:
            fail("--strict accepted the %s command line" % label, output)
        if "-P/--parallel needs a positive integer" not in output:
            fail("--strict did not report the %s value" % label, output)
    print("parallel: --strict makes the report fatal, in either position")

    # (5) more threads than the machine has processors: reported, applied as
    #     given, and only a warning with `--strict` as well, because nothing
    #     was ignored
    for label, flags in (("oversubscribed", []),
                         ("oversubscribed-strict", ["--strict"])):
        status, output = single_point(
            flags + ["-P", "1000", "--gfn", "ff", "--sp", "--norestart"],
            label)
        if status != 0:
            fail("the %s run failed with status %d" % (label, status), output)
        if "-P/--parallel 1000 exceeds the" not in output:
            fail("the %s run did not report the thread count" % label, output)
        if parallel_team(output) != 1000:
            fail("the %s run clamped the value to %d"
                 % (label, parallel_team(output)), output)
    print("parallel: a value above the processor count warns, unclamped")


def check_flag_value(exe, inp, rundir):
    """An option-looking token outside the flag character set is reported.

    Mirrors `xtb/landscape-flag-value`.  `TArgument%isFlag` used to require
    every character to be a flag character (`[A-Za-z0-9_-]`), so `--gfn=0`
    carried a character the parser did not recognise as flag input, was
    classified as *not* a flag and was never offered to the unknown-option
    reporter: it vanished and the run silently computed the GFN2 default.
    `=` is deliberately not supported, so the token is unknown-option input
    and the check pins both directions of the lax/strict split.  The last
    section pins what must **not** change: a numeric option value, the
    discouraged `-gfn` spelling and a dash-leading file name.
    """
    os.makedirs(rundir, exist_ok=True)

    def energy(args, label):
        directory = os.path.join(rundir, label)
        shutil.rmtree(directory, ignore_errors=True)
        status, output = run(
            exe, inp, args + ["--sp", "--norestart"], directory
        )
        if status != 0:
            fail("'%s' failed with status %d" % (" ".join(args), status), output)
        match = re.search(r"TOTAL ENERGY +(-?[0-9]+\.[0-9]+)", output)
        if not match:
            fail("no TOTAL ENERGY for '%s'" % " ".join(args), output)
        return match.group(1), output

    # (1) lax: every option-looking token is warned about and the run continues
    for token in ("--gfn=0", "--gfn=zzz", "--gfn=", "--gfn:0"):
        label = "lax" + re.sub(r"[^A-Za-z0-9]", "_", token)
        _, output = energy([token], label)
        if ("Unknown option '%s' provided" % token) not in output:
            fail("the option '%s' was silently dropped" % token, output)

    # (2) the dropped token was not applied: `--gfn=0` gives the GFN2 default
    equals, output = energy(["--gfn=0"], "equals-zero")
    default, _ = energy([], "default")
    if "Unknown option '--gfn=0' provided" not in output:
        fail("--gfn=0 was not reported", output)
    if equals != default:
        fail("--gfn=0 was applied (%s) instead of reported (default %s)"
             % (equals, default), output)

    # (3) strict: the same token becomes fatal
    directory = os.path.join(rundir, "strict")
    shutil.rmtree(directory, ignore_errors=True)
    status, output = run(
        exe, inp, ["--strict", "--gfn=zzz", "--sp", "--norestart"], directory
    )
    if status == 0:
        fail("--strict --gfn=zzz was accepted", output)
    if "Unknown option '--gfn=zzz' provided" not in output:
        fail("--strict --gfn=zzz was not reported as an unknown option", output)

    # (4) the inputs that work today must stay untouched
    zero, output = energy(["--gfn", "0"], "gfn-zero")
    if "Unknown option" in output:
        fail("--gfn 0 was misreported as an unknown option", output)
    short, output = energy(["-gfn", "0"], "gfn-short")
    if short != zero:
        fail("-gfn 0 (%s) differs from --gfn 0 (%s)" % (short, zero))
    if "discouraged" not in output:
        fail("-gfn was not reported as the discouraged spelling", output)
    directory = os.path.join(rundir, "dash-file")
    shutil.rmtree(directory, ignore_errors=True)
    os.makedirs(directory, exist_ok=True)
    with open(os.path.join(directory, "-dash.xyz"), "w") as handle:
        handle.write("$sp\n")
    status, output = run(
        exe, inp,
        ["--gfn", "2", "--input", "-dash.xyz", "--sp", "--norestart"],
        directory,
    )
    if status != 0:
        fail("a dash-leading input file name was treated as an option", output)
    if "Unknown option" in output:
        fail("a dash-leading input file name was reported as an option", output)
    print("option-looking tokens reported, working inputs untouched")


def check_positional(exe, inp, rundir):
    """A positional argument the program cannot use aborts the run.

    `TArgument%isFile` is `.true.` only for a path that exists, and both
    `countFiles` and `nextFile` look at files alone, so `oxtb a.xyz b.xyz`
    counted one file and silently read `b.xyz`, and `oxtb water.xyz --lvib
    normal` dropped `normal`.  Both were a warning (`TArgParser%nextMissingFile`
    reports the second shape) until the owner ruling recorded in
    `docs/_simplify-recon/197`: the warning let a run continue with a command
    line the program had not honoured, and `--lvib normal` exited 0 having
    computed the frequencies of the **unrelaxed** geometry.

    Every case below therefore asserts the non-zero exit, the message that
    names the offending token, and that no result was produced at all; the
    single-file command line keeps working and a token that is a convergence
    level of the run type gets the flag that takes it as a hint.
    """
    def do(label, args, extra=None):
        directory = os.path.join(rundir, label)
        shutil.rmtree(directory, ignore_errors=True)
        os.makedirs(directory, exist_ok=True)
        for name in extra or []:
            shutil.copy(inp, os.path.join(directory, name))
        proc = subprocess.run(
            [exe] + args, cwd=directory, capture_output=True, text=False
        )
        return proc.returncode, decode(proc.stdout) + decode(proc.stderr)

    def energy(label, args, extra=None):
        status, output = do(label, args, extra)
        if status != 0:
            fail("'%s' failed with status %d" % (" ".join(args), status), output)
        match = re.search(r"TOTAL ENERGY +(-?[0-9]+\.[0-9]+)", output)
        if not match:
            fail("no TOTAL ENERGY for '%s'" % " ".join(args), output)
        return match.group(1), output

    def rejected(label, args, token, extra=None, hint=None):
        status, output = do(label, args, extra)
        if status == 0:
            fail("'%s' succeeded although the run could not honour it"
                 % " ".join(args), output)
        if "Input file '%s'" % token not in output:
            fail("'%s' does not name the unusable token %r"
                 % (" ".join(args), token), output)
        if "TOTAL ENERGY" in output:
            fail("'%s' reported a total energy although it aborted"
                 % " ".join(args), output)
        if hint is not None and hint not in output:
            fail("'%s' does not offer the actionable hint %r"
                 % (" ".join(args), hint), output)
        if hint is None and "did you mean" in output:
            fail("'%s' offers a hint for a token that is not a level"
                 % " ".join(args), output)
        return output

    common = ["--gfn", "2", "--sp", "--norestart"]

    # (1) a single existing file: the baseline, nothing is reported
    single, output = energy("single", [inp] + common)
    if "Input file '" in output:
        fail("a single input file was reported as unusable", output)

    # (2) two existing files: the run reads one, so the first is an error and
    #     the second is not read either - nothing is computed
    rejected("duplicate", [inp, "second.xyz"] + common, inp,
             extra=["second.xyz"])

    # (3) an existing file and a missing second one
    rejected("missing-second", [inp, "nope.xyz"] + common, "nope.xyz")

    # (4) a missing first file and an existing second one
    rejected("missing-first", ["nope.xyz", inp] + common, "nope.xyz")

    # (5) only a missing file: the input-file error stays, nothing runs
    status, output = do("only-missing", ["nope.xyz"] + common)
    if status == 0:
        fail("a run without an existing input file succeeded", output)
    if "No input file given" not in output:
        fail("only a missing file did not report the missing input", output)

    # (6) the dangerous case of the review: a level keyword passed as a
    #     positional argument used to exit 0 and compute the frequencies of
    #     the unrelaxed input geometry (`rejected` also pins that no result
    #     was reported at all)
    rejected("vib-level", [inp, "--gfn", "2", "--lvib", "normal",
                            "--norestart"], "normal")

    # (7) the hint: `normal` is a legal value of the level flag of `--lopt`
    rejected("lopt-level", [inp, "--gfn", "2", "--lopt", "normal",
                             "--norestart"], "normal",
             hint="did you mean --lopt-level normal?")

    print("positional arguments: an unusable one aborts the run (exit non "
          "zero, token named, no result), the level keyword gets its flag")


def check_color(exe, inp, rundir):
    """`--color` consumes its own argument (CLI redesign stage 4).

    Mirrors `xtb/landscape-color`.  The `case ('--color')` branch of the
    argument parser read `sec` without ever calling `args%nextArg(sec)`, so
    `sec` held whatever the *previous* flag had left there: `--sp --color
    never` aborted with "No color scheme provided for --color option" (the
    variable was never allocated), while `--charge 1 --color never` warned
    "Unknown color option '1'" and left the scheme unchanged.  The flag was
    therefore unusable.  The third run pins the genuine error path: a `--color`
    without a following argument still has to abort.
    """
    status, output = run(
        exe, inp, ["--gfn", "2", "--sp", "--norestart", "--color", "never"],
        rundir,
    )
    if status != 0:
        fail("--color never alone failed the run", output)
    if "No color scheme provided" in output:
        fail("the bare --color branch did not consume its argument", output)
    status, output = run(
        exe,
        inp,
        ["--gfn", "2", "--sp", "--charge", "1", "--norestart",
         "--color", "never"],
        rundir,
    )
    if status != 0:
        fail("--charge 1 --color never failed the run", output)
    if "Unknown color option" in output:
        fail("--color still reads the argument of the previous flag", output)
    status, output = run(
        exe, inp, ["--gfn", "2", "--sp", "--norestart", "--color"], rundir
    )
    if status == 0:
        fail("a --color without an argument was accepted", output)
    if "No color scheme provided" not in output:
        fail("the missing color scheme was not reported", output)


def check_gfn(exe, inp, rundir):
    """`--gfn` accepts exactly 0, 1, 2 and ff; every other value is an error.

    Mirrors `xtb/landscape-gfn-value`.  The argument of an unknown `--gfn`
    value used to reach only the integer parser of the xcontrol reader, whose
    "could not parse" warning left `set%gfn_method` at its default 2: `--gfn
    zzz`, `--gfn ffn`, `--gfn 3`, `--gfn 4` and `--gfn ''` all exited 0 and
    computed **GFN2**.  The same fall-through sat on the xcontrol route
    (`$gfn method=zzz`), so a run could answer a different question than the
    one that was asked without a single actionable message.  Both routes are
    pinned here against the two directions:

    * every illegal spelling has to abort with the list of valid values - the
      status alone is not enough, an unrelated crash must not pass the check;
    * every legal spelling has to run and produce the energy of that method,
      `ff` bitwise equal to `--gfnff` and the no-`--gfn` default bitwise equal
      to `--gfn 2`, so neither equivalence can be redefined silently later.

    Both routes of `set_gfn` share one literal decoder, so the same spellings
    are refused on the command line and in the xcontrol file.  The `02` and
    `+2` controls are the point of the xcontrol half: the old integer parser
    accepted them, i.e. the xcontrol route ran GFN2 while the command line
    already refused the same text.
    """
    valid_hint = ("valid values are 0 (GFN0-xTB), 1 (GFN1-xTB), "
                  "2 (GFN2-xTB) and ff (GFN-FF)")
    os.makedirs(rundir, exist_ok=True)

    def energy(args, label):
        directory = os.path.join(rundir, label)
        shutil.rmtree(directory, ignore_errors=True)
        status, output = run(
            exe, inp, args + ["--sp", "--norestart"], directory
        )
        if status != 0:
            fail("'%s' failed with status %d" % (" ".join(args), status), output)
        match = re.search(r"TOTAL ENERGY +(-?[0-9]+\.[0-9]+)", output)
        if not match:
            fail("no TOTAL ENERGY for '%s'" % " ".join(args), output)
        return match.group(1)

    # (1) every illegal spelling aborts and names the valid values
    for value in ("zzz", "ffn", "3", "4", ""):
        directory = os.path.join(rundir, "illegal-" + (value or "empty"))
        shutil.rmtree(directory, ignore_errors=True)
        status, output = run(
            exe, inp, ["--gfn", value, "--sp", "--norestart"], directory
        )
        if status == 0:
            fail("an illegal --gfn %r was accepted" % value, output)
        if "unknown method" not in output:
            fail("an illegal --gfn %r did not name the rejected value" % value,
                 output)
        if valid_hint not in output:
            fail("an illegal --gfn %r did not list the valid values" % value,
                 output)
        print("--gfn %r refused with status %d" % (value, status))

    # (2) the four legal spellings run and give the energy of their method
    energies = {}
    for value in ("0", "1", "2", "ff"):
        energies[value] = energy(["--gfn", value], "legal-" + value)

    # `ff` is the GFN-FF force field, i.e. `--gfnff`, case insensitively
    gfnff = energy(["--gfnff"], "gfnff")
    if energies["ff"] != gfnff:
        fail("--gfn ff (%s) differs from --gfnff (%s)"
             % (energies["ff"], gfnff))
    upper = energy(["--gfn", "FF"], "legal-FF")
    if upper != gfnff:
        fail("--gfn FF (%s) differs from --gfnff (%s)" % (upper, gfnff))

    # no --gfn at all is the documented default GFN2
    default = energy([], "default")
    if default != energies["2"]:
        fail("the default method (%s) differs from --gfn 2 (%s)"
             % (default, energies["2"]))
    print("gfn values: gfn0 %s, gfn1 %s, gfn2 %s, gfnff %s"
          % (energies["0"], energies["1"], energies["2"], gfnff))

    # (3) the xcontrol route validates the same way and, for every legal
    # literal, selects the method whose energy `--gfn` produced above.  The
    # illegal spellings include values that are integers to a lenient parser -
    # `02` and `+2` used to slip through on this route alone - so the two
    # spellings of `2` are the alignment this check exists for.
    def xcontrol(value, label):
        xcontrol_file = os.path.join(rundir, "gfn-%s.inp" % label)
        directory = os.path.join(rundir, "xcontrol-" + label)
        with open(xcontrol_file, "w") as handle:
            handle.write("$gfn\n   method=%s\n$end\n" % value)
        shutil.rmtree(directory, ignore_errors=True)
        return run(
            exe, inp, ["--sp", "--input", xcontrol_file, "--norestart"],
            directory,
        )

    for value in ("zzz", "02", "+2", "2.0", "' 2'"):
        label = "illegal-" + re.sub(r"[^A-Za-z0-9]+", "_", value)
        status, output = xcontrol(value, label)
        if status == 0:
            fail("an xcontrol $gfn method=%s was accepted" % value, output)
        if "unknown method" not in output or valid_hint not in output:
            fail("the xcontrol $gfn method=%s did not list the valid values"
                 % value, output)
        print("xcontrol $gfn method=%s refused with status %d"
              % (value, status))

    for value in ("0", "1", "2", "ff"):
        status, output = xcontrol(value, "legal-" + value)
        if status != 0:
            fail("an xcontrol $gfn method=%s failed with status %d"
                 % (value, status), output)
        match = re.search(r"TOTAL ENERGY +(-?[0-9]+\.[0-9]+)", output)
        if not match or match.group(1) != energies[value]:
            fail("an xcontrol $gfn method=%s did not match --gfn %s"
                 % (value, value), output)
        print("xcontrol $gfn method=%s matches --gfn %s" % (value, value))


def check_dock_gfn(exe, mol_a, rundir, mol_b):
    """`oxtb dock` validates its `--gfn` with the shared decoder.

    Mirrors `xtb/landscape-dock-gfn-value`.  `dock` is a submodule of its own:
    `src/prog/primary.f90` dispatches on `argv[1]` and never enters `xtbMain`,
    so the strict `--gfn` of the main program cannot see the docking command
    line.  `src/prog/dock.f90` had four bare `if` statements instead, without
    an `else`: every value they did not match left `optlvl` at its GFN2 default
    and the run silently computed GFN2 - including the legal `--gfn FF` of the
    main program, i.e. one flag with two meanings depending on the submodule.

    The two entries of the docking module (`--gfn <value>` and its documented
    spelling `--optlvl gfn <value>`) now call the same
    `gfn_method_from_literal`/`gfn_method_error_message` pair as the main
    program, and this check pins the aligned behaviour on both:

    * an unknown spelling aborts with the shared list of valid values;
    * `0` aborts with the docking message `GFN0 is not supported for docking`
      and *not* with the shared text: GFN0 is a docking restriction, not an
      unknown spelling;
    * `1`, `2` and `ff` select their method, case insensitively, so `FF` is
      GFN-FF exactly as in the main program, and the no-`--gfn` default stays
      GFN2.

    `--fast` shrinks the search (`maxgen 7`, `n_opt 4`, measured 0.5 s instead
    of 16 s) and is orthogonal to the flag under test; the method line is
    printed before the search starts.  Every run is bounded by a timeout, so a
    hang cannot block the suite.
    """
    valid_hint = ("valid values are 0 (GFN0-xTB), 1 (GFN1-xTB), "
                  "2 (GFN2-xTB) and ff (GFN-FF)")
    method_line = re.compile(r"Method for final opts\.\s*:\s*(\S+)")

    def dock(args, label):
        directory = os.path.join(rundir, label)
        shutil.rmtree(directory, ignore_errors=True)
        os.makedirs(directory, exist_ok=True)
        try:
            proc = subprocess.run(
                [exe, "dock", mol_a, mol_b] + args,
                cwd=directory,
                capture_output=True,
                text=False,
                timeout=300,
            )
        except subprocess.TimeoutExpired:
            fail("dock %s did not finish" % " ".join(args))
        return proc.returncode, decode(proc.stdout) + decode(proc.stderr)

    def method(output):
        match = method_line.search(output)
        return match.group(1) if match else None

    # (1) every illegal spelling is refused by both entries with the shared text
    for value in ("zzz", "ffn", "02", "2.0"):
        status, output = dock(["--gfn", value, "--fast"],
                              "illegal-gfn-" + re.sub(r"[^A-Za-z0-9]+", "_", value))
        if status == 0:
            fail("dock --gfn %r was accepted" % value, output)
        if "unknown method" not in output or valid_hint not in output:
            fail("dock --gfn %r did not print the shared refusal" % value, output)
        print("dock --gfn %r refused with status %d" % (value, status))

    status, output = dock(["--optlvl", "gfn", "zzz", "--fast"], "illegal-optlvl")
    if status == 0:
        fail("dock --optlvl gfn zzz was accepted", output)
    if "unknown method" not in output or valid_hint not in output:
        fail("dock --optlvl gfn zzz did not print the shared refusal", output)
    print("dock --optlvl gfn zzz refused with status %d" % status)

    # (2) the missing argument is an error on both entries, like in the main
    # program, and not a silent default either
    for args, label in ((["--fast", "--gfn"], "missing-gfn"),
                        (["--fast", "--optlvl", "gfn"], "missing-optlvl")):
        status, output = dock(args, label)
        if status == 0:
            fail("dock %s was accepted" % " ".join(args), output)
        if "No method provided for --gfn option" not in output:
            fail("dock %s did not name the missing method" % " ".join(args),
                 output)
        print("dock %s refused with status %d" % (" ".join(args), status))

    # (3) GFN0 keeps the docking message: the literal is a legal one, the
    # restriction is the docking module's own
    for args, label in ((["--gfn", "0", "--fast"], "gfn0"),
                        (["--optlvl", "gfn", "0", "--fast"], "optlvl-gfn0")):
        status, output = dock(args, label)
        if status == 0:
            fail("dock %s was accepted" % " ".join(args), output)
        if "GFN0 is not supported for docking" not in output:
            fail("dock %s lost the docking GFN0 message" % " ".join(args),
                 output)
        if "unknown method" in output:
            fail("dock %s answered GFN0 with the shared refusal"
                 % " ".join(args), output)
        print("dock %s refused with the docking GFN0 message" % " ".join(args))

    # (4) the legal spellings select their method; `Ff` pins the case
    # insensitivity of the shared decoder on the docking route (`ff` is the
    # only literal that has letters, `GFN2` is not a legal spelling on either
    # route).  The energy of the lowest interaction is read as well: on this
    # fixture GFN1, GFN2 and GFN-FF give -4.96, -5.17 and -6.47 kcal/mol, so
    # `FF` equal to `ff` and different from `2` proves that the four `if`
    # statements of the old code did not merely mislabel the result.
    energy_line = re.compile(r"Lowest Interaction Energy: +(-?[0-9]+\.[0-9]+)")
    energies = {}
    for args, expected, label in (([], "gfn2", "default"),
                                  (["--gfn", "1"], "gfn1", "gfn1"),
                                  (["--gfn", "2"], "gfn2", "gfn2"),
                                  (["--gfn", "ff"], "gfnff", "gfnff"),
                                  (["--gfn", "FF"], "gfnff", "gfnff-upper"),
                                  (["--gfn", "Ff"], "gfnff", "gfnff-mixed"),
                                  (["--optlvl", "gfn", "1"], "gfn1", "optlvl-gfn1")):
        status, output = dock(args + ["--fast"], "legal-" + label)
        if status != 0:
            fail("dock %s failed with status %d" % (" ".join(args), status),
                 output)
        selected = method(output)
        if selected != expected:
            fail("dock %s selected %r instead of %r"
                 % (" ".join(args), selected, expected), output)
        match = energy_line.search(output)
        if not match:
            fail("dock %s printed no lowest interaction energy"
                 % " ".join(args), output)
        energies[label] = match.group(1)
        print("dock %s -> %s (%s kcal/mol)"
              % (" ".join(args) or "(no --gfn)", selected, energies[label]))

    if energies["gfnff"] != energies["gfnff-upper"]:
        fail("dock --gfn ff (%s) and --gfn FF (%s) differ"
             % (energies["gfnff"], energies["gfnff-upper"]))
    if energies["gfn1"] == energies["gfn2"]:
        fail("GFN1 and GFN2 give the same energy on the fixture, so the method "
             "check cannot falsify a wrong dispatch")
    if energies["gfnff"] in (energies["gfn1"], energies["gfn2"]):
        fail("GFN-FF gives the energy of GFN1 or GFN2 on the fixture, so the "
             "method check cannot falsify a wrong dispatch")


def compare_enum_case(exe, inp, rundir, label, result, pattern, upper, lower):
    """Run one task with upper and lower case enumeration values.

    Used by `check_enum_case`; mirrors the `compare_enum_case` function of the
    CMake script of `xtb/landscape-enum-case`.
    """
    upperdir = os.path.join(rundir, "upper-" + label)
    lowerdir = os.path.join(rundir, "lower-" + label)
    shutil.rmtree(upperdir, ignore_errors=True)
    shutil.rmtree(lowerdir, ignore_errors=True)
    upper_status, upper_output = run(exe, inp, upper, upperdir)
    lower_status, lower_output = run(exe, inp, lower, lowerdir)
    if upper_status != lower_status:
        fail("%s: the upper case run exits %d, the lower case run %d"
             % (label, upper_status, lower_status), upper_output)
    if "Unknown landscape" in upper_output:
        fail("%s: the upper case value was rejected" % label, upper_output)
    if not re.search(pattern, upper_output):
        fail("%s: the value did not reach the task" % label, upper_output)
    paths = [os.path.join(directory, result)
             for directory in (upperdir, lowerdir)]
    for name in paths:
        if not os.path.exists(name):
            fail("%s: %s was not written" % (label, name), upper_output)
    contents = []
    for name in paths:
        with open(name, "rb") as handle:
            contents.append(handle.read())
    if contents[0] != contents[1]:
        fail("%s: '%s' of the upper case run differs from the lower case run"
             % (label, result), upper_output)


def check_enum_case(exe, inp, rundir, second):
    """The suite enumeration names are matched case insensitively.

    Mirrors `xtb/landscape-enum-case` (CLI redesign sections 7.5 and 17.5 5b).
    Every case runs the same task with upper and lower case values and
    requires the same exit status and a byte identical result file; the byte
    comparison is what makes a missing normalisation of a *stored* value
    visible (`--lopt-hessian MODEL` would pass the name check and then run the
    numerical Hessian, `--lscan-mode RIGID` the relaxed scan and
    `--lts-init-mode FORCE` the `h0` start).  `second` is the HCN fixture: the
    reaction path needs a transition state as its input, so its cases cannot
    use the water minimum.
    """
    os.makedirs(rundir, exist_ok=True)
    common = ["--gfn", "2", "--norestart"]
    compare_enum_case(
        exe, inp, rundir, "lopt", "xtbopt.xyz", r"Hessian source: +model",
        common + ["--lopt", "--lopt-engine", "BFGS", "--lopt-level", "TIGHT",
                  "--lopt-hessian", "MODEL"],
        common + ["--lopt", "--lopt-engine", "bfgs", "--lopt-level", "tight",
                  "--lopt-hessian", "model"],
    )
    compare_enum_case(
        exe, inp, rundir, "lscan", "xtbscan.log", "landscape: rigid scan",
        common + ["--lscan-coord", "DISTANCE", "--lscan-atoms", "1,2",
                  "--lscan-step", "0.05", "--lscan-nstep", "1",
                  "--lscan-oneway", "--lscan-mode", "RIGID",
                  "--lscan-method", "CG", "--lscan-restart", "WARM"],
        common + ["--lscan-coord", "distance", "--lscan-atoms", "1,2",
                  "--lscan-step", "0.05", "--lscan-nstep", "1",
                  "--lscan-oneway", "--lscan-mode", "rigid",
                  "--lscan-method", "cg", "--lscan-restart", "warm"],
    )
    compare_enum_case(
        exe, inp, rundir, "lts", "xtbts.xyz", "dimer initial axis: +force",
        common + ["--lts", "--lts-method", "dimer", "--lts-ninit", "FORCE",
                  "--lts-max-steps", "1"],
        common + ["--lts", "--lts-method", "dimer", "--lts-ninit", "force",
                  "--lts-max-steps", "1"],
    )
    compare_enum_case(
        exe, inp, rundir, "prfo", "xtbts.xyz", "prfo",
        common + ["--lts", "--lts-method", "PRFO", "--lts-init-mode", "FORCE",
                  "--lts-hessian", "MODEL", "--lts-escape-step", "0",
                  "--lts-max-steps", "1"],
        common + ["--lts", "--lts-method", "prfo", "--lts-init-mode", "force",
                  "--lts-hessian", "model", "--lts-escape-step", "0",
                  "--lts-max-steps", "1"],
    )
    # The remaining `landscape_lowercase` stores of `main.F90` (review
    # finding M3 of `20-stage1-review.md`, zero coverage entries Z4-Z8).  The
    # five Hessian sources are observed through the `Hessian source:` line the
    # task prints from the settings object (the provider factory lower cases
    # internally, so the stored spelling is only visible there); the two
    # restart policies fail loudly when the `select case` is not normalised
    # (`case default` reports `requires off, pin or warm`), so the status
    # comparison already catches the mutation.
    compare_enum_case(
        exe, inp, rundir, "lscanhess", "xtbscan.log",
        r"Hessian source: +model",
        common + ["--lscan-coord", "distance", "--lscan-atoms", "1,2",
                  "--lscan-step", "0.05", "--lscan-nstep", "1",
                  "--lscan-oneway", "--lscan-mode", "rigid",
                  "--lscan-hessian", "MODEL"],
        common + ["--lscan-coord", "distance", "--lscan-atoms", "1,2",
                  "--lscan-step", "0.05", "--lscan-nstep", "1",
                  "--lscan-oneway", "--lscan-mode", "rigid",
                  "--lscan-hessian", "model"],
    )
    compare_enum_case(
        exe, inp, rundir, "lvibhess", "vibspectrum",
        r"Hessian source: +model",
        common + ["--lvib", "--lvib-hessian", "MODEL"],
        common + ["--lvib", "--lvib-hessian", "model"],
    )
    compare_enum_case(
        exe, inp, rundir, "ltshess", "xtbts.xyz",
        r"Hessian source: +model",
        common + ["--lts", "--lts-method", "dimer", "--lts-ninit", "force",
                  "--lts-hessian", "MODEL", "--lts-max-steps", "1"],
        common + ["--lts", "--lts-method", "dimer", "--lts-ninit", "force",
                  "--lts-hessian", "model", "--lts-max-steps", "1"],
    )
    compare_enum_case(
        exe, inp, rundir, "ltsrestart", "xtbts.xyz",
        r"restart policy: +warm",
        common + ["--lts", "--lts-method", "prfo", "--lts-init-mode", "force",
                  "--lts-hessian", "model", "--lts-escape-step", "0",
                  "--lts-max-steps", "1", "--lts-restart", "WARM"],
        common + ["--lts", "--lts-method", "prfo", "--lts-init-mode", "force",
                  "--lts-hessian", "model", "--lts-escape-step", "0",
                  "--lts-max-steps", "1", "--lts-restart", "warm"],
    )
    compare_enum_case(
        exe, second, rundir, "lircrestart", "irc.log",
        r"restart policy: +warm",
        common + ["--lirc", "--lirc-method", "gs", "--lirc-max-steps", "1",
                  "--lirc-restart", "WARM"],
        common + ["--lirc", "--lirc-method", "gs", "--lirc-max-steps", "1",
                  "--lirc-restart", "warm"],
    )
    # The stored `lircset%hessian_source` used to have no observation point
    # (zero coverage entry R1 of `21-stage1-fix-report.md`): the provider
    # factory lower cases internally and the task did not print the source.
    # The startup banner of `--lirc` prints the effective source now, which
    # makes the stored spelling observable here.  The value is `numerical`
    # (the default) because the model Hessian of this transition state guess
    # has no imaginary mode and the run aborts before writing `irc.log`
    # (measured; the `WARM` restart is needed for the same reason as in
    # `lircrestart` above).
    compare_enum_case(
        exe, second, rundir, "lirchessian", "irc.log",
        r"Hessian source: +numerical",
        common + ["--lirc", "--lirc-method", "gs", "--lirc-max-steps", "1",
                  "--lirc-hessian", "NUMERICAL", "--lirc-restart", "WARM"],
        common + ["--lirc", "--lirc-method", "gs", "--lirc-max-steps", "1",
                  "--lirc-hessian", "numerical", "--lirc-restart", "warm"],
    )


def read_thermo(text):
    """Enthalpy and free energy (hartree) of a thermochemistry printout."""
    match = re.search(r"total free energy\s+(-?[\d.]+) Eh", text)
    if not match:
        fail("no total free energy in the output", text)
    gibbs = float(match.group(1))
    match = re.search(r"TOTAL ENTHALPY\s+(-?[\d.]+) Eh", text)
    if match:
        enthalpy = float(match.group(1))
    else:
        match = re.search(r"total enthalpy H\(T\):\s+(-?[\d.]+) Eh", text)
        if not match:
            fail("no total enthalpy in the output", text)
        enthalpy = float(match.group(1))
    return enthalpy, gibbs


def read_json(rundir, name):
    """Load a `--ljson` document of `rundir`."""
    path = os.path.join(rundir, name)
    if not os.path.exists(path):
        fail("no %s in %s" % (name, rundir))
    with open(path) as handle:
        return json.load(handle)


def check_thermo(exe, inp, rundir):
    """`--lvib` thermochemistry follows the `$thermo` priority rule.

    Mirrors `xtb/landscape-vib-thermo`.  Priority: explicit command line >
    xcontrol `$thermo` > builtin default.  Before the fix `--lvib` overwrote
    all four keys of the block unconditionally, so `$thermo temp=298.15,300,400`
    collapsed to the single 298.15 K the settings carried and `scale=0.9` was
    ignored, while `--hess` honoured both: the two routes disagreed by 7.5e-3
    hartree on one geometry whose validated agreement is 1e-6 (the projector
    spread, docs/_landscape-recon/49-*.md).  The comparisons below are made at
    a *tight* minimum, which is where that tolerance is calibrated.
    """
    os.makedirs(rundir, exist_ok=True)
    rundir = os.path.abspath(rundir)
    inputs = {
        "list": "$thermo\n   temp=298.15,300.0,400.0\n$end\n",
        "scale": "$thermo\n   scale=0.9000\n   temp=298.15\n$end\n",
        "temp298": "$thermo\n   temp=298.15\n$end\n",
        "temp400": "$thermo\n   temp=400.0\n$end\n",
    }
    for name, block in inputs.items():
        with open(os.path.join(rundir, name + ".inp"), "w") as handle:
            handle.write(block)

    # the thermochemistry table of `--hess` and `--lvib` has to agree, so the
    # geometry must be a (tight) minimum
    status, output = run(exe, inp, ["--gfn", "2", "--lopt", "--lopt-level", "tight",
                                    "--norestart"], rundir)
    if status != 0:
        fail("the thermochemistry optimisation terminated with status %d" % status,
             output)
    geometry = os.path.join(rundir, "xtbopt.xyz")
    if not os.path.exists(geometry):
        fail("the thermochemistry optimisation wrote no xtbopt.xyz", output)

    def thermo(task, args):
        status, output = run(exe, geometry, ["--gfn", "2", "--" + task,
                                             "--norestart"] + args, rundir)
        if status != 0:
            fail("xtb --%s %s terminated with status %d"
                 % (task, " ".join(args), status), output)
        return read_thermo(output), output

    for name in ("list", "scale"):
        (vh, vg), vib_out = thermo("lvib", ["-I", name + ".inp"])
        (hh, hg), hess_out = thermo("hess", ["-I", name + ".inp"])
        if abs(vg - hg) > 1.0e-6 or abs(vh - hh) > 1.0e-6:
            fail("the --lvib and --hess thermochemistry disagree on %s by "
                 "%.3e/%.3e hartree (H/G, tolerance 1e-6)\n--lvib:\n%s\n--hess:\n%s"
                 % (name, vh - hh, vg - hg, vib_out, hess_out))
        print("thermo: %-5s --lvib vs --hess  dH %.3e, dG %.3e hartree"
              % (name, vh - hh, vg - hg))

    # `temp=298.15,300.0,400.0` must produce one table row per temperature; the
    # summary reports the *last* one, here 400 K (the same value `--hess`
    # prints)
    (_, list_g), list_out = thermo("lvib", ["-I", "list.inp"])
    for temperature in ("298.15", "300.00", "400.00"):
        if not re.search(r"^ *%s +VIB" % re.escape(temperature), list_out, re.M):
            fail("the `$thermo temp` list lost the %s K row" % temperature, list_out)
    (_, temp400_g), _ = thermo("lvib", ["-I", "temp400.inp"])
    if abs(list_g - temp400_g) > 1.0e-9:
        fail("`temp=298.15,300.0,400.0` does not report the 400 K value "
             "(%.12f vs %.12f)" % (list_g, temp400_g))

    (_, default_g), _ = thermo("lvib", [])
    if abs(default_g - temp400_g) < 1.0e-9:
        fail("the xcontrol `temp=400` did not override the builtin 298.15 K")
    # explicit command line beats the xcontrol block
    (_, cli_g), _ = thermo("lvib", ["-I", "temp298.inp", "--lvib-temp", "400"])
    if abs(cli_g - temp400_g) > 1.0e-9:
        fail("--lvib-temp 400 did not override `$thermo temp=298.15` "
             "(%.12f vs %.12f)" % (cli_g, temp400_g))
    # ... and the xcontrol block beats the builtin default (same value here)
    (_, xc_g), _ = thermo("lvib", ["-I", "temp400.inp"])
    if abs(xc_g - cli_g) > 1.0e-9:
        fail("`$thermo temp=400` differs from `--lvib-temp 400` "
             "(%.12f vs %.12f)" % (xc_g, cli_g))
    print("thermo: CLI `--lvib-temp 400` == xcontrol `temp=400` == %.12f Eh" % cli_g)

    # the same for `scale`: the CLI value and the xcontrol value are the same
    # 0.9 on the same (builtin) temperature
    (scaled_h, scaled_g), scaled_out = thermo("lvib", ["-I", "scale.inp"])
    (cli_h, cli_g2), cli_out = thermo("lvib", ["--lvib-freq-scal", "0.9"])
    if abs(scaled_g - cli_g2) > 1.0e-9 or abs(scaled_h - cli_h) > 1.0e-9:
        fail("`$thermo scale=0.9` and `--lvib-freq-scal 0.9` differ "
             "(%.12f vs %.12f)\n%s\n%s" % (scaled_g, cli_g2, scaled_out, cli_out))
    if abs(scaled_g - default_g) < 1.0e-9:
        fail("the xcontrol `scale=0.9` did not change the thermochemistry")
    print("thermo: xcontrol `scale=0.9` == CLI `--lvib-freq-scal 0.9` == %.12f Eh"
          % scaled_g)

    # --- `--ljson` records the *effective* temperature, not the settings -----
    # default (finding MEDIUM-1 of `05-w0-independent-review.md`): before this
    # round `--lvib --ljson -I temp400.inp` wrote `temperature: 298.15` next to
    # the 400 K `gibbs` of the same document, because the JSON read the
    # `settings%` fields while the thermochemistry read the merged `$thermo`
    # values.  Two runs with the same `gibbs` were computed at the same
    # temperature, hence they must report the same `temperature` - and that
    # number has to be the 400 K the block asked for.
    def thermo_json(args):
        status, output = run(exe, geometry, ["--gfn", "2", "--lvib", "--ljson",
                                            "--norestart"] + args, rundir)
        if status != 0:
            fail("xtb --lvib --ljson %s terminated with status %d"
                 % (" ".join(args), status), output)
        return read_json(rundir, "landscape_vib.json"), read_thermo(output)[1], output

    xc_doc, xc_json_g, xc_json_out = thermo_json(["-I", "temp400.inp"])
    cli_doc, cli_json_g, cli_json_out = thermo_json(["--lvib-temp", "400"])
    if abs(xc_json_g - cli_json_g) > 1.0e-9:
        fail("the JSON runs of `$thermo temp=400` and `--lvib-temp 400` differ "
             "(%.12f vs %.12f)" % (xc_json_g, cli_json_g), xc_json_out)
    xc_temp = xc_doc["values"]["temperature"]
    if abs(xc_temp - 400.0) > 1.0e-9:
        fail("landscape_vib.json reports temperature %.6f K for a run whose "
             "gibbs is the 400 K value" % xc_temp, xc_json_out)
    if abs(xc_temp - cli_doc["values"]["temperature"]) > 1.0e-9:
        fail("the JSON temperature of the `$thermo temp=400` run (%.6f) and of "
             "the `--lvib-temp 400` run (%.6f) disagree"
             % (xc_temp, cli_doc["values"]["temperature"]), cli_json_out)
    print("thermo json: `$thermo temp=400` and `--lvib-temp 400` both report "
          "temperature %.3f K for gibbs %.12f Eh" % (xc_temp, xc_json_g))

    sc_doc, _, sc_json_out = thermo_json(["-I", "scale.inp"])
    cli_sc_doc, _, cli_sc_json_out = thermo_json(["--lvib-freq-scal", "0.9"])
    if abs(sc_doc["values"]["freq_scal"] - 0.9) > 1.0e-12:
        fail("landscape_vib.json reports freq_scal %.9f for a `$thermo "
             "scale=0.9` run" % sc_doc["values"]["freq_scal"], sc_json_out)
    if abs(sc_doc["values"]["freq_scal"]
           - cli_sc_doc["values"]["freq_scal"]) > 1.0e-12:
        fail("the JSON scaling factors of the `$thermo scale=0.9` and the "
             "`--lvib-freq-scal 0.9` run disagree (%r vs %r)"
             % (sc_doc["values"]["freq_scal"], cli_sc_doc["values"]["freq_scal"]),
             cli_sc_json_out)
    print("thermo json: `$thermo scale=0.9` == `--lvib-freq-scal 0.9` -> "
          "freq_scal %.3f" % sc_doc["values"]["freq_scal"])

    # --- unusable values of the `$thermo` block are rejected -----------------
    # The validation runs on the *merged* state, so a bad block value is caught
    # exactly like a bad `--lvib-*` flag and before the first single point
    # (finding P2-2 of `05-w0-independent-review.md`).  `set_thermo` checks
    # nothing, so `$thermo scale=-1` used to be multiplied into every
    # frequency (measured: it changed G by 1.9e-2 hartree instead of aborting).
    invalid = {
        "badscale": ("scale=-1.0", "the frequency scaling factor must be positive"),
        "badtemp": ("temp=-5.0", "the thermochemical temperature must be positive"),
        "badimagthr": ("imagthr=20.0",
                       "the imaginary frequency cutoff must not be positive"),
        "badsthr": ("sthr=-1.0",
                    "the rotor cutoff of the quasi-RRHO treatment must not be negative"),
    }
    for name, (key, pattern) in invalid.items():
        with open(os.path.join(rundir, name + ".inp"), "w") as handle:
            handle.write("$thermo\n   %s\n$end\n" % key)
        status, output = run(exe, geometry,
                             ["--gfn", "2", "--lvib", "--norestart",
                              "-I", name + ".inp"], rundir)
        if status == 0:
            fail("the `$thermo %s` block was accepted" % key, output)
        if pattern not in output:
            fail("the `$thermo %s` block was rejected without the documented "
                 "message" % key, output)
        print("thermo: `$thermo %-12s` rejected with status %d" % (key, status))


def check_thermo_knobs(exe, inp, rundir, saddle):
    """The rotor and imaginary cutoffs of `$thermo` take effect (`--lvib`).

    Mirrors `xtb/landscape-vib-thermo-knobs`.  The `$thermo` block has exactly
    four keys and the frequency settings carry one `*_was_set` guard per key
    (`landscape_vib_settings`); `xtb/landscape-vib-thermo` covers `temp` and
    `scale`, so this check closes `sthr` (`--lvib-rotor-cutoff`) and `imagthr`
    (`--lvib-imag-cutoff`).  Both are asserted twice: an *effect* (the number
    changes) and the *merge* (the xcontrol key and the CLI flag give the same
    result).  Deleting the two guards of `sthr`/`imagthr`
    (`landscape_vib_set_thermo_settings`) used to leave the whole suite green
    - the CLI flag had no test at all (`sthr`) or only its error branch
    (`imagthr`) - which is finding P2-1 of `05-w0-independent-review.md`.
    """
    os.makedirs(rundir, exist_ok=True)
    rundir = os.path.abspath(rundir)
    saddle_dir = os.path.join(rundir, "saddle")
    os.makedirs(saddle_dir, exist_ok=True)

    # --- sthr / --lvib-rotor-cutoff, on a tight water minimum ---------------
    # The quasi-RRHO treatment replaces the harmonic entropy of the modes below
    # the cutoff by the Grimme free-rotor interpolation, so the cutoff changes
    # G (measured 1.5e-5 hartree here, H is untouched).
    status, output = run(exe, inp, ["--gfn", "2", "--lopt", "--lopt-level",
                                    "tight", "--norestart"], rundir)
    if status != 0:
        fail("the rotor cutoff optimisation terminated with status %d" % status,
             output)
    geometry = os.path.join(rundir, "xtbopt.xyz")
    if not os.path.exists(geometry):
        fail("the rotor cutoff optimisation wrote no xtbopt.xyz", output)
    with open(os.path.join(rundir, "sthr.inp"), "w") as handle:
        handle.write("$thermo\n   sthr=2000.0\n$end\n")

    def water_thermo(args):
        status, output = run(exe, geometry,
                             ["--gfn", "2", "--lvib", "--norestart"] + args,
                             rundir)
        if status != 0:
            fail("xtb --lvib %s terminated with status %d"
                 % (" ".join(args), status), output)
        return read_thermo(output)[1]

    sthr_default = water_thermo([])
    sthr_cli = water_thermo(["--lvib-rotor-cutoff", "2000"])
    sthr_xc = water_thermo(["-I", "sthr.inp"])
    if abs(sthr_cli - sthr_default) < 1.0e-9:
        fail("--lvib-rotor-cutoff 2000 did not change the free energy "
             "(%.12f vs %.12f)" % (sthr_cli, sthr_default))
    if abs(sthr_cli - sthr_xc) > 1.0e-12:
        fail("`--lvib-rotor-cutoff 2000` and `$thermo sthr=2000` differ "
             "(%.12f vs %.12f)" % (sthr_cli, sthr_xc))
    print("thermo knobs: rotor cutoff  --lvib-rotor-cutoff 2000 == "
          "$thermo sthr=2000 = %.12f Eh (default %.12f, dG %.3e)"
          % (sthr_cli, sthr_default, sthr_cli - sthr_default))

    # --- imagthr / --lvib-imag-cutoff, on the bent HCN saddle --------------
    # A mode whose magnitude is below the cutoff is treated as an artifact of
    # the projected Hessian and `print_thermo` inverts its sign
    # (`if (vibs(i) < 0 .and. vibs(i) > ithr)`, `src/main/property.F90`).  The
    # fixture has one imaginary mode of -1077 cm^-1, so a cutoff of -10000 cm^-1
    # has to turn it into a real mode: the analysis counts `n_imag` 0 instead
    # of 1 (status `incomplete` -> `ok`) and the mode enters the partition
    # function (measured dG 2.4e-3 hartree).
    #
    # The two routes are compared on the *thermochemistry*, not on `n_imag`:
    # the suite's own mode analysis (`landscape_vib_analysis`) classifies with
    # `settings%imag_cutoff`, so an xcontrol `imagthr` reaches the reused model
    # (`set%thermo_ithr`) but not the analysis counter - the divergence is
    # visible in the warning of `landscape_vib_thermo_output` ("the
    # thermochemistry counted a different number of imaginary modes ...") and
    # is recorded as an open finding in `docs/_simplify-recon/06-*.md`.  What
    # this check pins is the *effect* of the CLI flag and that the block
    # reaches the same effective cutoff (same G, same printed cutoff).
    with open(os.path.join(saddle_dir, "imagthr.inp"), "w") as handle:
        handle.write("$thermo\n   imagthr=-10000.0\n$end\n")

    def saddle_json(args):
        status, output = run(exe, saddle,
                             ["--gfn", "2", "--lvib", "--ljson",
                              "--norestart"] + args, saddle_dir)
        if status != 0:
            fail("xtb --lvib --ljson %s terminated with status %d"
                 % (" ".join(args), status), output)
        return read_json(saddle_dir, "landscape_vib.json"), output

    default_doc, default_out = saddle_json([])
    cli_doc, cli_out = saddle_json(["--lvib-imag-cutoff", "-10000"])
    xc_doc, xc_out = saddle_json(["-I", "imagthr.inp"])
    if default_doc["counters"]["n_imag"] != 1:
        fail("the HCN fixture does not have one imaginary mode (%d)"
             % default_doc["counters"]["n_imag"], default_out)
    if cli_doc["counters"]["n_imag"] != 0:
        fail("--lvib-imag-cutoff -10000 did not invert the imaginary mode "
             "(n_imag %d)" % cli_doc["counters"]["n_imag"], cli_out)
    if cli_doc["status"] != 0 or cli_doc["status_name"] != "ok":
        fail("the run without an imaginary mode is not `ok` (%d/%s)"
             % (cli_doc["status"], cli_doc["status_name"]), cli_out)
    if abs(cli_doc["values"]["gibbs"]
           - default_doc["values"]["gibbs"]) < 1.0e-6:
        fail("the imaginary cutoff did not change the free energy "
             "(%.12f vs %.12f)"
             % (cli_doc["values"]["gibbs"], default_doc["values"]["gibbs"]),
             default_out)
    if abs(cli_doc["values"]["gibbs"] - xc_doc["values"]["gibbs"]) > 1.0e-12:
        fail("`--lvib-imag-cutoff -10000` and `$thermo imagthr=-10000` differ "
             "in G (%.12f vs %.12f)"
             % (cli_doc["values"]["gibbs"], xc_doc["values"]["gibbs"]), xc_out)
    match = re.search(r"imag[.] cutoff\s+(-?[\d.]+) cm", xc_out)
    if not match or abs(float(match.group(1)) + 10000.0) > 1.0e-9:
        fail("`$thermo imagthr=-10000` did not reach the thermochemistry",
             xc_out)
    if "inverting freq" not in xc_out:
        fail("`$thermo imagthr=-10000` did not invert the imaginary mode",
             xc_out)
    print("thermo knobs: imag cutoff   --lvib-imag-cutoff -10000 == "
          "$thermo imagthr=-10000 -> analytic n_imag 0, dG %.3e Eh (default "
          "n_imag %d, status %s)"
          % (cli_doc["values"]["gibbs"] - default_doc["values"]["gibbs"],
             default_doc["counters"]["n_imag"], default_doc["status_name"]))


def check_saddlefix(exe, inp, rundir):
    """`$fix` holds the frozen atom for the single ended saddle searches.

    Mirrors `xtb/landscape-lts-fix`; the input is the bent HCN transition state
    guess whose atom 1 is the terminal H at `(0, 1.05, 0)` angstrom, and
    `$fix atoms: 1` freezes it.  Before the projection the accepted step of
    every method moved it (measured: prfo 0.1171, bprfo 0.1475, dimer 0.1482
    angstrom - docs/_simplify-recon/04-*.md section 5.2).  After the fix the
    coordinate is bitwise the input value, so the tolerance is the xyz
    write/read round trip, not a measured spread.

    `prfo` does **not** converge on this constrained fixture (it exhausts the
    step budget, measured with 120/400/1200 steps - the same happens before the
    projection, so the fixture is the reason, not the projection).  The check
    therefore asserts the honest contract of either outcome: a converged run
    exits 0, a non-converged run exits non-zero and says so.  No convergence is
    faked and the step budget is not tuned to make the fixture pass.
    """
    # rebuilt: the check below asserts the *absence* of `xtbts.xyz` for the
    # method that does not converge, and a leftover file of an earlier run (a
    # manual run, or the injected run of the falsifiability experiment) would
    # satisfy that existence check (the same trap as finding F1 of
    # `50-m7-independent-review.md`)
    rundir = os.path.abspath(rundir)
    shutil.rmtree(rundir, ignore_errors=True)
    os.makedirs(rundir, exist_ok=True)
    xcontrol = os.path.join(rundir, "landscape_fix.inp")
    with open(xcontrol, "w") as handle:
        handle.write("$fix\n   atoms: 1\n$end\n")
    reference = None
    with open(inp) as handle:
        for line in handle.read().splitlines()[2:4]:
            fields = line.split()
            if fields and fields[0] == "H":
                reference = [float(value) for value in fields[1:4]]
    if reference is None:
        fail("the saddle fixture has no hydrogen at atom 1", inp)
    for method in ("prfo", "bprfo", "dimer"):
        geometry = os.path.join(rundir, "xtbts.xyz")
        if os.path.exists(geometry):
            os.remove(geometry)
        status, output = run(
            exe,
            inp,
            ["--gfn", "2", "--lts", "--lts-method", method,
             "--lts-max-steps", "120", "--input", xcontrol, "--norestart"],
            rundir,
        )
        converged = "landscape: converged" in output
        if not os.path.exists(geometry):
            fail("the --lts-method %s run wrote no xtbts.xyz" % method, output)
        with open(geometry) as handle:
            lines = handle.read().splitlines()
        fields = lines[2].split()
        if fields[0] != "H":
            fail("the first output atom of %s is not the frozen hydrogen" % method,
                 lines[2])
        shift = max(abs(float(value) - ref)
                    for value, ref in zip(fields[1:4], reference))
        if shift > 1.0e-10:
            fail(
                "--lts-method %s moved the frozen atom by %.5e angstrom (%s)"
                % (method, shift, " ".join(fields[1:4])),
                output,
            )
        if method == "prfo":
            # `prfo` does not converge on this constrained fixture (it exhausts
            # the step budget, measured with 120/400/1200 steps - the same
            # happens before the projection, so the fixture is the reason, not
            # the projection).  What the check pins is therefore **not**
            # "prfo must not converge": it is that the outcome is reported
            # honestly - a converged run exits 0, a non-converged one exits
            # non-zero with the message of `landscape_ts` ("transition state
            # search did not converge: <reason>").  Asserting the surfaced
            # reason instead ("maximum number of steps reached") would turn a
            # future improvement of the trust radius into a red test for the
            # wrong reason (finding P2-4 of `05-w0-independent-review.md`).
            if converged:
                if status != 0:
                    fail("the constrained prfo run reported convergence with "
                         "status %d" % status, output)
            else:
                if status == 0:
                    fail("the constrained prfo run exited 0 although it did "
                         "not converge", output)
                if "transition state search did not converge" not in output:
                    fail("the constrained prfo run failed without reporting "
                         "the missing convergence", output)
        else:
            if status != 0 or not converged:
                fail("--lts-method %s did not converge (status %d)"
                     % (method, status), output)
        print("saddle fix: %-5s status %d converged %s  frozen atom displacement "
              "%.5e angstrom" % (method, status, converged, shift))


def check_rejectfix(exe, inp, rundir):
    """The methods without `$fix` support must refuse the run, not move atoms.

    Mirrors `xtb/landscape-lirc-fix`.  The band searches of `--lts` and the
    reaction path `--lirc` have no defined meaning for a frozen coordinate yet
    (interpolated images / the mass weighted arc length contract), so they
    abort with an actionable message instead of running in a geometry that
    violates the constraint (docs/_simplify-recon/04-*.md sections 5.2/5.3).
    """
    # rebuilt: the check asserts the *absence* of the artefacts the refused run
    # would have written when the guard is injected (finding F1 of
    # `50-m7-independent-review.md`)
    rundir = os.path.abspath(rundir)
    shutil.rmtree(rundir, ignore_errors=True)
    os.makedirs(rundir, exist_ok=True)
    xcontrol = os.path.join(rundir, "landscape_fix.inp")
    with open(xcontrol, "w") as handle:
        handle.write("$fix\n   atoms: 1\n$end\n")
    cases = [
        ("--lirc", ["--lirc", "--lirc-method", "hpc", "--lirc-max-steps", "4"],
         "the reaction path (--lirc) does not support $fix"),
        ("--lts neb", ["--lts", "--lts-method", "neb", "--lts-max-steps", "4"],
         "the band searches (--lts-method neb|cineb|autoneb|string|gsm) do not "
         "support $fix"),
    ]
    for label, args, pattern in cases:
        status, output = run(
            exe, inp, ["--gfn", "2"] + args + ["--input", xcontrol, "--norestart"],
            rundir,
        )
        if status == 0:
            fail("%s accepted $fix although the method does not support it" % label,
                 output)
        if pattern not in output:
            fail("%s did not explain the rejected $fix block" % label, output)
        for artifact in ("xtbts.xyz", "irc_forward.xyz", "irc_path.xyz"):
            if os.path.exists(os.path.join(rundir, artifact)):
                fail("%s wrote %s although it refused the $fix block"
                     % (label, artifact), output)
        print("reject fix: %-9s rejected $fix with status %d" % (label, status))


def check_ts(exe, inp, rundir):
    """`--lts` locates a first order saddle point and writes its two files.

    Mirrors the CMake script of `xtb/landscape-ts-prfo`.  The input is a bent
    HCN, a guess for the HCN/HNC transition state; P-RFO converges in a few
    steps to a stationary point with exactly one imaginary mode, which is the
    definition of the requested saddle point.  The frozen electronic
    temperature is part of the output (M4 hard constraint 1).
    """
    for name in ("xtbts.xyz", "xtbts.log"):
        path = os.path.join(rundir, name)
        if os.path.exists(path):
            os.remove(path)
    status, output = run(
        exe, inp, ["--gfn", "2", "--lts", "--lts-method", "prfo", "--norestart"], rundir
    )
    if status != 0:
        fail("xtb --lts terminated with status %d" % status, output)
    if "landscape: converged" not in output:
        fail("the transition state search did not converge", output)
    if "saddle confirmation - imaginary modes: 1" not in output:
        fail("the located point is not a first order saddle point", output)
    if not re.search(r"landscape engine: etemp +300[.]000 K", output):
        fail("the frozen electronic temperature was not reported", output)
    for name in ("xtbts.xyz", "xtbts.log"):
        path = os.path.join(rundir, name)
        if not os.path.exists(path):
            fail("%s was not written" % name)
        if os.path.getsize(path) == 0:
            fail("%s is empty" % name)
    with open(os.path.join(rundir, "xtbts.log")) as handle:
        log = handle.read()
    if "confirmation imaginary modes: 1" not in log:
        fail("the log does not report the imaginary mode count", log)
    # the restricted step branch has to reach the same kind of stationary point
    status, output = run(
        exe, inp, ["--gfn", "2", "--lts", "--lts-method", "bprfo", "--norestart"], rundir
    )
    if status != 0:
        fail("xtb --lts --lts-method bprfo terminated with status %d" % status, output)
    if "landscape: converged" not in output:
        fail("the B-PRFO search did not converge", output)
    if "saddle confirmation - imaginary modes: 1" not in output:
        fail("the B-PRFO point is not a first order saddle point", output)


def check_dimer(exe, inp, rundir):
    """`--lts --lts-method dimer` locates a saddle point without a Hessian.

    Mirrors the CMake script of `xtb/landscape-ts-dimer`.  The search starts
    from the deterministic `force` initial axis (the largest Cartesian
    gradient component direction, PLAN section 6.3), converges to a first
    order saddle point of the bent HCN, writes `xtbts.xyz`/`xtbts.log` and
    prints its rotation/torque and curvature diagnostics together with the
    frozen electronic temperature.
    """
    for name in ("xtbts.xyz", "xtbts.log"):
        path = os.path.join(rundir, name)
        if os.path.exists(path):
            os.remove(path)
    status, output = run(
        exe, inp,
        ["--gfn", "2", "--lts", "--lts-method", "dimer", "--lts-ninit", "force",
         "--norestart"],
        rundir,
    )
    if status != 0:
        fail("xtb --lts --lts-method dimer terminated with status %d" % status, output)
    if "dimer initial axis: force" not in output:
        fail("the deterministic initial axis was not used", output)
    if "landscape: converged" not in output:
        fail("the dimer search did not converge", output)
    if "saddle confirmation - imaginary modes: 1" not in output:
        fail("the dimer point is not a first order saddle point", output)
    if not re.search(r"landscape engine: etemp +300[.]000 K", output):
        fail("the frozen electronic temperature was not reported", output)
    if not re.search(r"dimer rotation: iterations [0-9]+, torque rms", output):
        fail("the rotation/torque diagnostic is missing", output)
    if "dimer curvature:" not in output:
        fail("the curvature diagnostic is missing", output)
    for name in ("xtbts.xyz", "xtbts.log"):
        path = os.path.join(rundir, name)
        if not os.path.exists(path):
            fail("%s was not written" % name)
        if os.path.getsize(path) == 0:
            fail("%s is empty" % name)
    with open(os.path.join(rundir, "xtbts.log")) as handle:
        log = handle.read()
    if "confirmation imaginary modes: 1" not in log:
        fail("the log does not report the imaginary mode count", log)
    if "dimer" not in log:
        fail("the log does not identify the dimer method", log)


def check_ts_minimum(exe, inp, rundir):
    """A search that converges to a minimum is not a transition state.

    Mirrors the CMake script of `xtb/landscape-ts-minimum` (H2 of
    `29-m4-independent-review.md`): from a linear HNC the P-RFO search lands on
    the molecule's own minimum (0 imaginary modes) when the escape out of a
    minimum is disabled (`--lts-escape-step 0`); the run must abort with a
    non-zero status, print the diagnosis, and still write `xtbts.xyz`/
    `xtbts.log` for debugging.  With the default escape the same input finds
    the saddle point (`check_ts_escape`).

    The `--ljson` document of the same run must report `status` -1 (`error`),
    because it reports the confirmation's rejection of the point with a message
    that the task also raises through `env`; `1` (`incomplete`) is the verdict
    of a search that ran out of steps, which this fixture is not (report 224,
    finding 1).
    """
    for name in ("xtbts.xyz", "xtbts.log", "landscape_ts.json"):
        path = os.path.join(rundir, name)
        if os.path.exists(path):
            os.remove(path)
    status, output = run(
        exe, inp,
        ["--gfn", "2", "--lts", "--lts-method", "prfo", "--norestart",
         "--lts-escape-step", "0", "--ljson"],
        rundir,
    )
    if status == 0:
        fail("a minimum was accepted as a transition state", output)
    if "saddle confirmation - imaginary modes: 0" not in output:
        fail("the imaginary mode count is not reported", output)
    if "confirmation: minimum" not in output:
        fail("the minimum was not named as such", output)
    if "not a first order saddle point" not in output:
        fail("the failure was not reported", output)
    for name in ("xtbts.xyz", "xtbts.log"):
        if not os.path.exists(os.path.join(rundir, name)):
            fail("%s was not written despite the failure" % name)
    path = os.path.join(rundir, "landscape_ts.json")
    try:
        with open(path) as handle:
            doc = json.load(handle)
    except (OSError, ValueError) as error:
        fail("cannot read the document of the rejected point: %s" % error,
             output)
    if doc.get("status") != -1 or doc.get("status_name") != "error":
        fail("the rejected point reports status %s / %s, expected -1 / error"
             % (doc.get("status"), doc.get("status_name")), output)
    if doc.get("converged"):
        fail("the rejected point still reports converged", output)
    if doc["counters"].get("n_imag") != 0:
        fail("the document reports n_imag %s, expected 0"
             % doc["counters"].get("n_imag"), output)
    if "the located stationary point is a minimum" not in (
            doc.get("message") or ""):
        fail("the document lost the confirmation's message, it says %r"
             % doc.get("message"), output)


def check_ts_escape(exe, inp, rundir):
    """The escape walk turns a minimum start into a successful search.

    Mirrors the CMake script of `xtb/landscape-ts-escape` (MEDIUM-1 of
    `31-m4-review-round2.md`): from the same linear minimum as
    `check_ts_minimum`, the driver walks the geometry along the lowest
    internal mode until its curvature turns negative and the search then
    reaches the HCN/HNC saddle point with exactly one imaginary mode.
    """
    for name in ("xtbts.xyz", "xtbts.log"):
        path = os.path.join(rundir, name)
        if os.path.exists(path):
            os.remove(path)
    status, output = run(
        exe, inp, ["--gfn", "2", "--lts", "--lts-method", "prfo", "--norestart"], rundir
    )
    if status != 0:
        fail("the escaped search terminated with status %d" % status, output)
    if not re.search(r"escape from a minimum: +[0-9]+ push", output):
        fail("the escape walk out of the minimum was not reported", output)
    if not re.search(r"escape from a minimum: +0[.]200000 bohr", output):
        fail("the configured escape displacement was not reported", output)
    if not re.search(r"selected mode, rigid body overlap: +0[.]000000", output):
        fail("a rigid body direction was followed", output)
    if "saddle confirmation - imaginary modes: 1" not in output:
        fail("the escaped search did not reach a first order saddle", output)
    if "confirmation: first order saddle point" not in output:
        fail("the saddle confirmation is missing", output)
    for name in ("xtbts.xyz", "xtbts.log"):
        if not os.path.exists(os.path.join(rundir, name)):
            fail("%s was not written" % name)
    with open(os.path.join(rundir, "xtbts.log")) as handle:
        log = handle.read()
    if "confirmation imaginary modes: 1" not in log:
        fail("the log does not report the imaginary mode count", log)
    status, output = run(
        exe, inp, ["--gfn", "2", "--lts", "--lts-method", "bprfo", "--norestart"], rundir
    )
    if status != 0:
        fail("the escaped B-PRFO search terminated with status %d" % status, output)
    if not re.search(r"escape from a minimum: +[0-9]+ push", output):
        fail("the B-PRFO escape walk was not reported", output)
    if "saddle confirmation - imaginary modes: 1" not in output:
        fail("the escaped B-PRFO search did not reach a saddle", output)


def check_ts_namespace(exe, inp, rundir):
    """`--namespace ns --lts` writes `ns.xtbts.xyz`/`ns.xtbts.log`.

    Mirrors the CMake script of `xtb/landscape-ts-namespace` (M3 of
    `29-m4-independent-review.md`): the file names pass through
    `get_namespace` exactly like `--lscan` and `--lvib`, so a namespaced run
    does not collide with a plain one.
    """
    for name in ("landscape6.xtbts.xyz", "landscape6.xtbts.log"):
        path = os.path.join(rundir, name)
        if os.path.exists(path):
            os.remove(path)
    status, output = run(
        exe, inp,
        ["--gfn", "2", "--lts", "--lts-method", "prfo", "--norestart",
         "--namespace", "landscape6"],
        rundir,
    )
    if status != 0:
        fail("xtb --lts --namespace terminated with status %d" % status, output)
    for name in ("landscape6.xtbts.xyz", "landscape6.xtbts.log"):
        if not os.path.exists(os.path.join(rundir, name)):
            fail("%s was not written (namespace not applied)" % name)


def check_ts_etemp(exe, inp, rundir):
    """The CLI `--etemp` must reach the transition state engine.

    Mirrors the CMake script of `xtb/landscape-ts-etemp`: the settings carry the
    frozen temperature and `main.F90` has to synchronise them from `set%etemp`,
    so a missing sync would leave the task at 300 K.  Both the frozen value and
    the value read back from the calculator are checked against 5000 K.
    """
    status, output = run(
        exe, inp,
        ["--gfn", "2", "--lts", "--lts-method", "prfo", "--etemp", "5000",
         "--norestart"],
        rundir,
    )
    if status != 0:
        fail("xtb --lts --etemp 5000 terminated with status %d" % status, output)
    if not re.search(r"electronic temperature: +5000[.]000 K", output):
        fail("the frozen electronic temperature is not 5000 K", output)
    if not re.search(r"landscape engine: etemp +5000[.]000 K", output):
        fail("--etemp 5000 did not reach the engine control (missing CLI sync)", output)
    # MEDIUM-3: the result field itself (`etemp:` of `landscape_result%write_log`
    # and the frozen value of `xtbts.log`) must carry the effective temperature
    if not re.search(r"etemp: +5000[.]000 K", output):
        fail("the result object does not carry the effective etemp", output)
    with open(os.path.join(rundir, "xtbts.log")) as handle:
        log = handle.read()
    if not re.search(r"electronic temperature \(frozen\): +5000[.]000 K", log):
        fail("xtbts.log does not carry the effective etemp", log)


def check_maxsteps(exe, inp, rundir):
    """`--lscan-max-steps` reaches the inner relaxation (M2 leftovers M-1/M-2).

    Mirrors the CMake script of `xtb/landscape-scan-max-steps`: with a zero
    step budget no point converges, the run reports a missing minimum (NaN)
    and a readable diagnostic with a non-zero exit status, and every evaluated
    point still reaches the scan log.  With a real budget the same scan
    converges and exits 0.
    """
    args = [
        "--gfn", "2", "--lscan", "--lscan-coord", "angle",
        "--lscan-atoms", "2,1,3", "--lscan-step", "3", "--lscan-nstep", "2",
        "--lscan-max-steps", "0", "--norestart", "--namespace", "landscape4",
    ]
    status, output = run(exe, inp, args, rundir)
    # M7 (PLAN §12.11 INFO-3): a scan without a single converged point has no
    # usable result, so the run exits non-zero next to the NaN total energy
    # (before this round the branch only warned and the run exited 0).
    if status == 0:
        fail("a scan without a converged point exited 0", output)
    if "scan completed without a converged point" not in output:
        fail("the missing convergence was not reported", output)
    if not re.search(r"lowest energy: +NaN", output):
        fail("the unconverged scan reported a minimum", output)
    with open(os.path.join(rundir, "landscape4.xtbscan.log")) as handle:
        frames = len(re.findall(r"energy:", handle.read()))
    if frames != 5:
        fail("xtbscan.log holds %d frames, expected 5" % frames)
    args = [
        "--gfn", "2", "--lscan", "--lscan-coord", "angle",
        "--lscan-atoms", "2,1,3", "--lscan-step", "3", "--lscan-nstep", "2",
        "--lscan-max-steps", "256", "--lscan-level", "loose", "--norestart",
        "--namespace", "landscape4",
    ]
    status, output = run(exe, inp, args, rundir)
    if status != 0:
        fail("the reference scan failed with status %d" % status, output)
    if "landscape task terminated: scan completed" not in output:
        fail("the reference scan did not complete", output)
    if "without a converged point" in output:
        fail("the reference scan reported a missing minimum", output)


def check_autoneb_depth_limit(exe, inp, rundir, product):
    """An AutoNEB run that stops at the insertion cap must name the cap.

    Mirrors the CMake script of `xtb/landscape-path-autoneb-depth-limit`: the
    band of the last stage converges, but the requested path was not built, so
    the result message - and therefore `xtbts.log` - has to report the cap
    instead of a bare "converged" (report 38 L-2).
    """
    args = [
        "--gfn", "2", "--lts", "--lts-method", "autoneb",
        "--lts-product", product, "--lts-autoneb-max-depth", "1",
        "--lts-max-steps", "400", "--norestart", "--namespace", "mesonP",
    ]
    status, output = run(exe, inp, args, rundir)
    if status != 0:
        fail("the capped AutoNEB run terminated with status %d" % status, output)
    if "autoneb: budget limit reached" not in output:
        fail("the AutoNEB insertion cap was not reached", output)
    cap = "landscape task terminated: stopped at the insertion depth limit max_depth = 1"
    if cap not in output:
        fail("the cap did not reach the result message", output)
    with open(os.path.join(rundir, "mesonP.xtbts.log")) as handle:
        log = handle.read()
    if "termination: stopped at the insertion depth limit max_depth = 1" not in log:
        fail("the cap did not reach xtbts.log", log)
    if "termination: converged" in log:
        fail("xtbts.log claims convergence for a capped run", log)


def check_scanhess(exe, inp, rundir):
    """`--lscan-hessian` reaches the engine of the scan points."""
    args = [
        "--gfn", "2", "--lscan", "--lscan-method", "rfo",
        "--lscan-hessian", "model", "--lscan-coord", "angle",
        "--lscan-atoms", "2,1,3", "--lscan-step", "3", "--lscan-nstep", "1",
        "--lscan-level", "loose", "--norestart",
    ]
    status, output = run(exe, inp, args, rundir)
    if status != 0:
        fail("--lscan-hessian model failed with status %d" % status, output)
    if not re.search(r"Hessian source: +model", output):
        fail("the scan Hessian source was not reported", output)


def check_basin(label, energy, lower, upper):
    """Require an endpoint energy to lie inside the measured basin window."""
    value = float(energy)
    if value > upper or value < lower:
        fail("the %s endpoint %s is outside the basin window (%s .. %s)"
             % (label, energy, lower, upper))


def read_log(rundir, name):
    """Read one IRC log file of the run directory."""
    with open(os.path.join(rundir, name)) as handle:
        return handle.read()


def endpoint_energy(log, direction):
    """Endpoint energy of one direction of an IRC log."""
    match = re.search(r"%s endpoint energy: +(-?[0-9.]+)" % direction, log)
    if not match:
        fail("the IRC log reports no %s endpoint energy" % direction, log)
    return match.group(1)


def check_irc(exe, inp, rundir):
    """`--lirc` from a converged transition state (xtb/landscape-irc-cli).

    The input is the bent HCN transition state **guess**: the check first runs
    `--lts` to converge it, because an IRC started from a non-stationary
    geometry funnels both directions into the same basin (report 42).
    """
    status, output = run(
        exe, inp,
        ["--gfn", "2", "--lts", "--lts-method", "prfo", "--norestart",
         "--namespace", "irc0"],
        rundir,
    )
    if status != 0:
        fail("the preparation --lts run failed with status %d" % status, output)
    ts = os.path.join(rundir, "irc0.xtbts.xyz")
    if not os.path.exists(ts):
        fail("the preparation --lts run wrote no xtbts.xyz", output)
    irc_args = ["--gfn", "2", "--lirc", "--lirc-method", "gs", "--norestart",
                "--namespace", "irc1"]
    status, output = run(exe, ts, irc_args, rundir)
    if status != 0:
        fail("xtb --lirc terminated with status %d" % status, output)
    for pattern in [
        r"landscape: intrinsic reaction coordinate [(]gs[)]",
        r"transition state imaginary modes: +1",
        r"phase component: +[0-9.]+",
        r"landscape engine: etemp +300[.]000 K",
        r"forward reached minimum: +T",
        r"backward reached minimum: +T",
    ]:
        if not re.search(pattern, output):
            fail("missing %r in the --lirc output" % pattern, output)
    for name in ["irc1.irc_forward.xyz", "irc1.irc_backward.xyz",
                 "irc1.irc_path.xyz", "irc1.irc.log"]:
        path = os.path.join(rundir, name)
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            fail("%s was not written" % name, output)
    with open(os.path.join(rundir, "irc1.irc_path.xyz")) as handle:
        pathfile = handle.read()
    if " TS" not in pathfile:
        fail("the transition state is not marked in irc_path.xyz", pathfile)
    with open(os.path.join(rundir, "irc1.irc.log")) as handle:
        log = handle.read()
    match = re.search(r"forward endpoint energy: +(-?[0-9.]+)", log)
    if not match:
        fail("the IRC log reports no forward endpoint energy", log)
    fwd = match.group(1)
    match = re.search(r"backward endpoint energy: +(-?[0-9.]+)", log)
    if not match:
        fail("the IRC log reports no backward endpoint energy", log)
    if fwd == match.group(1):
        fail("both IRC directions reached the same energy (%s)" % fwd)
    bwd = match.group(1)
    # golden values of the Gonzalez-Schlegel baseline (report 44 section 3):
    # the measured minima +/- 1e-4 Eh, 2x the measured same-basin spread of
    # the four methods (5.0e-5 Eh) and 8x narrower than the distance of the
    # one method outside the basin (eulerpc, 8.0e-4 Eh)
    check_basin("forward", fwd, -5.47225988, -5.47205988)
    check_basin("backward", bwd, -5.50416622, -5.50396622)
    for pattern in [r"forward arc length: +3[.]949160",
                    r"backward arc length: +3[.]538733",
                    r"mean energy/gradient evaluations per point: +1[.]88",
                    r"mean energy/gradient evaluations per point: +2[.]14",
                    r"endpoint rmsd criterion: rmsd <= rmsd_end = 0[.]18897000"]:
        if not re.search(pattern, log):
            fail("missing %r in irc.log" % pattern, log)
    # the equality constrained Lagrange step of the sphere relaxation (report
    # 64) has to show up as **zero** pinned points on both paths.  Before the
    # fix 33 of 42 forward and 38 of 40 backward points were pinned and the
    # check asserted `>= 1` on purpose; the golden arc lengths and E/g of the
    # loop above were updated with the same fix (4.051234 -> 3.949160,
    # 3.695670 -> 3.538733, 16.9/20.2 -> 1.88/2.14 E/g).
    match = re.search(r"pinned at the micro cycle budget: forward ([0-9]+) of", log)
    if not match or int(match.group(1)) != 0:
        fail("the GS forward sphere relaxation is pinned again", log)
    match = re.search(r"backward ([0-9]+) of [0-9]+ points", log)
    if not match or int(match.group(1)) != 0:
        fail("the GS backward sphere relaxation is pinned again", log)
    with open(os.path.join(rundir, "irc1.irc_forward.xyz")) as handle:
        first = handle.read()
    status, output = run(exe, ts, irc_args, rundir)
    if status != 0:
        fail("the repeated --lirc run failed with status %d" % status, output)
    with open(os.path.join(rundir, "irc1.irc_forward.xyz")) as handle:
        second = handle.read()
    if first != second:
        fail("irc_forward.xyz is not bitwise reproducible")


def check_irc_knobs(exe, inp, rundir):
    """`--lirc-level`, `--lirc-rmsd-end` and the default method (M-1).

    Each of the three has to make the reported numbers move: the two flags were
    parsed and then never read before the M6 fix round, and the default method
    was `gs` although `hpc` is both more accurate and 8.5x cheaper.
    """
    status, output = run(
        exe, inp,
        ["--gfn", "2", "--lts", "--lts-method", "prfo", "--norestart",
         "--namespace", "ircknob0"],
        rundir,
    )
    if status != 0:
        fail("the preparation --lts run failed with status %d" % status, output)
    ts = os.path.join(rundir, "ircknob0.xtbts.xyz")
    if not os.path.exists(ts):
        fail("the preparation --lts run wrote no xtbts.xyz", output)
    # the default method is hpc (M6 behaviour change, PLAN section 12.15)
    status, output = run(
        exe, ts, ["--gfn", "2", "--lirc", "--norestart", "--namespace", "irc4"],
        rundir,
    )
    if status != 0:
        fail("xtb --lirc without --lirc-method failed with status %d" % status, output)
    log = read_log(rundir, "irc4.irc.log")
    if not re.search(r"intrinsic reaction coordinate: hpc", log):
        fail("the default IRC method is no longer hpc", log)
    check_basin("default forward", endpoint_energy(log, "forward"),
                -5.47225988, -5.47205988)
    check_basin("default backward", endpoint_energy(log, "backward"),
                -5.50416622, -5.50396622)
    if not re.search(r"pinned at the micro cycle budget: forward 0 of", log):
        fail("the default (hpc) run reports pinned micro cycles", log)
    # `--lirc-rmsd-end` attributes the endpoint (the forward endpoint is
    # 0.0099 bohr from its relaxed minimum, the backward one 1.6e-5)
    status, output = run(
        exe, ts,
        ["--gfn", "2", "--lirc", "--lirc-method", "gs", "--lirc-rmsd-end", "1e-3",
         "--norestart", "--namespace", "irc5"],
        rundir,
    )
    if status != 0:
        fail("xtb --lirc --lirc-rmsd-end 1e-3 failed with status %d" % status, output)
    log = read_log(rundir, "irc5.irc.log")
    if not re.search(r"endpoint rmsd criterion: rmsd <= rmsd_end = 0[.]00100000", log):
        fail("--lirc-rmsd-end 1e-3 did not reach the settings", log)
    if not re.search(r"forward reached minimum: +F", log):
        fail("--lirc-rmsd-end 1e-3 did not reject the forward endpoint", log)
    if not re.search(r"backward reached minimum: +T", log):
        fail("--lirc-rmsd-end 1e-3 rejected the backward endpoint", log)
    # `--lirc-level` drives the endpoint relaxation (measured: 1 step at
    # superloose to -5.47213843 Eh, 7 steps at vtight to -5.47215989 Eh)
    minima = {}
    for level in ["superloose", "vtight"]:
        status, output = run(
            exe, ts,
            ["--gfn", "2", "--lirc", "--lirc-method", "gs", "--lirc-level", level,
             "--norestart", "--namespace", "irc-" + level],
            rundir,
        )
        if status != 0:
            fail("xtb --lirc-level %s failed with status %d" % (level, status), output)
        log = read_log(rundir, "irc-%s.irc.log" % level)
        match = re.search(
            r"forward relaxed minimum: +(-?[0-9.]+) Eh [(]([0-9]+) relaxation steps[)]",
            log,
        )
        if not match:
            fail("the %s log reports no relaxed minimum" % level, log)
        minima[level] = (float(match.group(1)), int(match.group(2)))
    if minima["vtight"][1] <= minima["superloose"][1]:
        fail("--lirc-level vtight did not need more relaxation steps than superloose "
             "(%d vs %d)" % (minima["vtight"][1], minima["superloose"][1]))
    if minima["vtight"][0] >= minima["superloose"][0]:
        fail("--lirc-level vtight did not relax further than superloose (%g vs %g)"
             % (minima["vtight"][0], minima["superloose"][0]))
    # without a relaxation there is no relaxed minimum: the attribution falls
    # back to the termination signal and the log has to say so
    status, output = run(
        exe, ts,
        ["--gfn", "2", "--lirc", "--lirc-method", "gs", "--lirc-no-endpoint-relax",
         "--norestart", "--namespace", "irc6"],
        rundir,
    )
    if status != 0:
        fail("xtb --lirc --lirc-no-endpoint-relax failed with status %d" % status, output)
    log = read_log(rundir, "irc6.irc.log")
    if not re.search(r"forward reached minimum: +T", log):
        fail("--lirc-no-endpoint-relax lost the termination based attribution", log)
    if re.search(r"forward relaxed minimum:", log):
        fail("--lirc-no-endpoint-relax still reports a relaxed minimum", log)
    if not re.search(r"the rmsd_end criterion was skipped", log):
        fail("--lirc-no-endpoint-relax does not report the skipped rmsd criterion", log)


def check_irc_pc(exe, inp, rundir):
    """`--lirc-method hpc` and `eulerpc` (xtb/landscape-irc-pc).

    Both predictor-corrector integrators have to reach a minimum in both
    directions, write the four output files and keep the mean E/g per point at
    the predictor-corrector level (see the module header of `irc/hpc.f90`).
    """
    status, output = run(
        exe, inp,
        ["--gfn", "2", "--lts", "--lts-method", "prfo", "--norestart",
         "--namespace", "ircpc0"],
        rundir,
    )
    if status != 0:
        fail("the preparation --lts run failed with status %d" % status, output)
    ts = os.path.join(rundir, "ircpc0.xtbts.xyz")
    if not os.path.exists(ts):
        fail("the preparation --lts run wrote no xtbts.xyz", output)
    for method in ["hpc", "eulerpc"]:
        status, output = run(
            exe, ts,
            ["--gfn", "2", "--lirc", "--lirc-method", method, "--norestart",
             "--namespace", method],
            rundir,
        )
        if status != 0:
            fail("xtb --lirc --lirc-method %s failed with status %d" % (method, status), output)
        for pattern in [
            r"landscape: intrinsic reaction coordinate [(]%s[)]" % method,
            r"transition state imaginary modes: +1",
            r"termination forward: +(gradient|energy_rise)",
            r"termination backward: +(gradient|energy_rise)",
            r"forward reached minimum: +T",
            r"backward reached minimum: +T",
        ]:
            if not re.search(pattern, output):
                fail("missing %r in the %s output" % (pattern, method), output)
        for name in ["%s.irc_forward.xyz" % method, "%s.irc_backward.xyz" % method,
                     "%s.irc_path.xyz" % method, "%s.irc.log" % method]:
            path = os.path.join(rundir, name)
            if not os.path.exists(path) or os.path.getsize(path) == 0:
                fail("%s was not written" % name, output)
        with open(os.path.join(rundir, "%s.irc.log" % method)) as handle:
            log = handle.read()
        match = re.search(r"mean energy/gradient evaluations per point: +([0-9.]+)", log)
        if not match:
            fail("%s: the log reports no mean E/g" % method, log)
        mean_ev = float(match.group(1))
        if not 1.0 <= mean_ev <= 6.0:
            fail("%s: the mean E/g per point is %s, far above the predictor-corrector level" % (method, mean_ev), log)
        if not re.search(r"forward endpoint energy: +-5[.]47", log):
            fail("%s: the forward endpoint did not reach the HNC basin" % method, log)
        if not re.search(r"backward endpoint energy: +-5[.]50", log):
            fail("%s: the backward endpoint did not reach the HCN basin" % method, log)
    # the new flags reach the settings
    status, output = run(
        exe, ts,
        ["--gfn", "2", "--lirc", "--lirc-method", "eulerpc",
         "--lirc-loose-cycles", "0", "--lirc-erise-steps", "5",
         "--lirc-endpoint-relax", "--norestart", "--namespace", "ircpc1"],
        rundir,
    )
    if status != 0:
        fail("the new --lirc flags failed with status %d" % status, output)
    # a path cut off by the step budget is not assigned to a minimum
    status, output = run(
        exe, ts,
        ["--gfn", "2", "--lirc", "--lirc-method", "hpc", "--lirc-max-steps", "3",
         "--norestart", "--namespace", "ircpc2"],
        rundir,
    )
    if status != 0:
        fail("the --lirc-max-steps 3 run failed with status %d" % status, output)
    if not re.search(r"termination forward: +max_steps", output):
        fail("the step budget did not report signal (a)", output)
    if not re.search(r"forward reached minimum: +F", output):
        fail("a path cut off by the step budget was assigned to a minimum", output)


def check_irc_termination(exe, inp, rundir):
    """Termination signal (c) at the command line (xtb/landscape-irc-termination)."""
    status, output = run(
        exe, inp,
        ["--gfn", "2", "--lts", "--lts-method", "prfo", "--norestart",
         "--namespace", "ircterm0"],
        rundir,
    )
    if status != 0:
        fail("the preparation --lts run failed with status %d" % status, output)
    ts = os.path.join(rundir, "ircterm0.xtbts.xyz")
    if not os.path.exists(ts):
        fail("the preparation --lts run wrote no xtbts.xyz", output)
    status, output = run(
        exe, ts,
        ["--gfn", "2", "--lirc", "--lirc-method", "gs", "--lirc-grad-tol", "0",
         "--norestart", "--namespace", "ircterm1"],
        rundir,
    )
    if status != 0:
        fail("xtb --lirc --lirc-grad-tol 0 failed with status %d" % status, output)
    for pattern in [
        r"termination forward: +energy_rise",
        r"termination backward: +energy_rise",
        r"forward reached minimum: +T",
        r"backward reached minimum: +T",
        r"landscape engine: etemp +300[.]000 K",
    ]:
        if not re.search(pattern, output):
            fail("missing %r in the --lirc output" % pattern, output)


def check_irc_swap(exe, inp, rundir):
    """`--lirc-swap` exchanges the two directions (xtb/landscape-irc-swap).

    Mirrors the inline CMake script of the CMake build.  The input is the bent
    HCN transition state guess, converged with `--lts` first (see
    `check_irc`).  The swap check would be vacuous if both directions reached
    the same energy, so that is asserted as well.
    """
    status, output = run(
        exe, inp,
        ["--gfn", "2", "--lts", "--lts-method", "prfo", "--norestart",
         "--namespace", "ircswap0"],
        rundir,
    )
    if status != 0:
        fail("the preparation --lts run failed with status %d" % status, output)
    ts = os.path.join(rundir, "ircswap0.xtbts.xyz")
    if not os.path.exists(ts):
        fail("the preparation --lts run wrote no xtbts.xyz", output)

    def case(suffix, swap):
        args = ["--gfn", "2", "--lirc", "--lirc-method", "gs", "--norestart",
                "--namespace", suffix]
        if swap:
            args.append("--lirc-swap")
        status, output = run(exe, ts, args, rundir)
        if status != 0:
            fail("xtb --lirc %s failed with status %d" % (suffix, status), output)
        log = read_log(rundir, "%s.irc.log" % suffix)
        return endpoint_energy(log, "forward"), endpoint_energy(log, "backward")

    plain_fwd, plain_bwd = case("ircswap1", False)
    swapped_fwd, swapped_bwd = case("ircswap2", True)
    if plain_fwd == plain_bwd:
        fail("both directions reached the same energy (%s); the swap check "
             "would be vacuous" % plain_fwd)
    if plain_fwd != swapped_bwd:
        fail("--lirc-swap did not exchange the forward endpoint: %s vs %s"
             % (plain_fwd, swapped_bwd))
    if plain_bwd != swapped_fwd:
        fail("--lirc-swap did not exchange the backward endpoint: %s vs %s"
             % (plain_bwd, swapped_fwd))


def check_irc_osc(exe, inp, rundir):
    """GS macro step limit cycle guard (xtb/landscape-irc-osc-window).

    Mirrors the inline CMake script of the CMake build.  The Gonzalez-Schlegel
    macro step stalls in a two cycle when the pivot sits on a basin floor, so
    `--lirc-osc-window` (default 8) leaves the path after that many consecutive
    macro steps that cancel each other (report 67, report 69).  The unit case
    `gs-limit-cycle-guard` pins the truncation on the analytic Mueller-Brown
    surface; this check shows the flag is live end to end from the other side:
    on the HCN/HNC path the default window is inert, window 1 fires at the
    first cancelling pair, window 0 reproduces the default, and a negative
    window is rejected by the task.
    """
    status, output = run(
        exe, inp,
        ["--gfn", "2", "--lts", "--lts-method", "prfo", "--norestart",
         "--namespace", "ircosc0"],
        rundir,
    )
    if status != 0:
        fail("the preparation --lts run failed with status %d" % status, output)
    ts = os.path.join(rundir, "ircosc0.xtbts.xyz")
    if not os.path.exists(ts):
        fail("the preparation --lts run wrote no xtbts.xyz", output)

    def case(suffix, window, extra=None):
        args = ["--gfn", "2", "--lirc", "--lirc-method", "gs",
                "--lirc-grad-tol", "0", "--norestart", "--namespace", suffix]
        if window is not None:
            args += ["--lirc-osc-window", str(window)]
        if extra:
            args.append(extra)
        status, output = run(exe, ts, args, rundir)
        if status != 0:
            fail("xtb --lirc (window %s) failed with status %d"
                 % (window, status), output)
        points = {}
        for direction in ["forward", "backward"]:
            match = re.search(r"%s points: +([0-9]+)" % direction, output)
            if not match:
                fail("the run reports no %s point count" % direction, output)
            points[direction] = int(match.group(1))
        match = re.search(r"forward arc length: +([0-9.]+)", output)
        if not match:
            fail("the run reports no forward arc length", output)
        return points, match.group(1), output

    default, default_arc, default_out = case("ircosc1", None)
    off, off_arc, off_out = case("ircosc2", 0)
    one, one_arc, one_out = case("ircosc3", 1)

    if re.search(r"limit cycle", default_out):
        fail("the default GS limit cycle guard fired on the HCN path", default_out)
    if default["forward"] != 45 or default["backward"] != 39:
        fail("the default GS path has %d/%d points, expected 45/39"
             % (default["forward"], default["backward"]), default_out)
    if re.search(r"limit cycle", off_out):
        fail("--lirc-osc-window 0 still reports a limit cycle", off_out)
    if off != default or off_arc != default_arc:
        fail("--lirc-osc-window 0 changed the path: %s (%s) vs %s (%s)"
             % (off, off_arc, default, default_arc))
    if not re.search(
            r"limit cycle detected [(]the last 1 macro steps cancel[)] - terminating",
            one_out):
        fail("--lirc-osc-window 1 did not fire", one_out)
    if one["forward"] >= default["forward"] or one["backward"] >= default["backward"]:
        fail("--lirc-osc-window 1 did not truncate the path: %s vs %s"
             % (one, default), one_out)
    # the truncation is its own termination signal (d), not signal (c): `irc.log`
    # names it, the warning says the endpoint is not a converged path point and
    # `--ljson` carries it as `term_fwd = 4` with an incomplete, non-converged
    # run whose endpoint is not assigned to a minimum (round-2 review HIGH-1)
    log = read_log(rundir, "ircosc3.irc.log")
    if not re.search(r"termination forward: +limit_cycle", log):
        fail("irc.log does not name the limit cycle termination", log)
    if not re.search(
            r"one phase of the oscillation, not a converged path point", one_out):
        fail("the limit cycle truncation warning is missing", one_out)
    status, output = run(
        exe, ts,
        ["--gfn", "2", "--lirc", "--lirc-method", "gs", "--lirc-grad-tol", "0",
         "--lirc-osc-window", "1", "--lirc-no-endpoint-relax", "--ljson",
         "--norestart", "--namespace", "ircosc6"],
        rundir,
    )
    if status != 0:
        fail("the --ljson limit cycle run failed with status %d" % status, output)
    doc = read_json(rundir, "ircosc6.landscape_irc.json")
    if doc["counters"]["term_fwd"] != 4:
        fail("--ljson does not report term_fwd = 4 (limit_cycle): %s"
             % doc["counters"]["term_fwd"], output)
    if doc["status_name"] != "incomplete":
        fail("a truncated IRC is not reported as incomplete (%s)"
             % doc["status_name"], output)
    if doc["converged"]:
        fail("a truncated IRC still reports converged = true", output)
    if doc["flags"]["reached_min_fwd"]:
        fail("the truncated forward endpoint was assigned to a minimum", output)
    if "limit cycle guard" not in doc["message"]:
        fail("the truncated run message does not name the guard (%r)"
             % doc["message"], output)
    # the opt-in adaptive step (`--lirc-adaptive`, default off) reduces the step
    # at the first cancelling pair instead of ending the walk: it has to report
    # the reduction, change the path and leave no limit cycle message
    adaptive, adaptive_arc, adaptive_out = case("ircosc5", None, "--lirc-adaptive")
    if not re.search(r"reducing the effective step length", adaptive_out):
        fail("--lirc-adaptive did not reduce the step length", adaptive_out)
    if re.search(r"limit cycle", adaptive_out):
        fail("--lirc-adaptive still reports a limit cycle", adaptive_out)
    if adaptive_arc == default_arc:
        fail("--lirc-adaptive did not change the forward arc length (%s)"
             % adaptive_arc)
    # a negative window is rejected by the task, so the value really reaches
    # `check_irc_settings` and is not parsed and forgotten
    status, output = run(
        exe, ts,
        ["--gfn", "2", "--lirc", "--lirc-method", "gs", "--lirc-osc-window", "-1",
         "--norestart", "--namespace", "ircosc4"],
        rundir,
    )
    if status == 0:
        fail("a negative --lirc-osc-window was accepted", output)
    if not re.search(r"the IRC limit cycle window must not be negative", output):
        fail("the negative --lirc-osc-window error does not name the window", output)
    # a direction that stopped on the step budget is as unconverged as a
    # truncated one: `status` already reported it as `incomplete`, but
    # `converged` still claimed the run succeeded (the sibling of the limit
    # cycle signal, found while fixing round 2 HIGH-1).  `--lirc-max-steps 2`
    # forces signal (a) on both directions.
    status, output = run(
        exe, ts,
        ["--gfn", "2", "--lirc", "--lirc-method", "gs", "--lirc-grad-tol", "0",
         "--lirc-max-steps", "2", "--lirc-no-endpoint-relax", "--ljson",
         "--norestart", "--namespace", "ircosc7"],
        rundir,
    )
    if status != 0:
        fail("the --ljson step budget run failed with status %d" % status, output)
    doc = read_json(rundir, "ircosc7.landscape_irc.json")
    if doc["counters"]["term_fwd"] != 1:
        fail("--lirc-max-steps 2 did not terminate on the step budget (term_fwd %s)"
             % doc["counters"]["term_fwd"], output)
    if doc["status_name"] != "incomplete":
        fail("a step budget stop is not reported as incomplete (%s)"
             % doc["status_name"], output)
    if doc["converged"]:
        fail("a step budget stop still reports converged = true", output)


def check_opt_namespace(exe, inp, rundir):
    """`--lopt` geometry and trajectory honour `--namespace`.

    Mirrors the CMake script of `xtb/landscape-opt-namespace` (defect N-1 of
    the final fix round): `xtb_landscape_log::landscape_write_files` used the
    bare names `xtbopt.xyz`/`xtb.trj` although every other landscape writer and
    the native geometry writer (the `open_file` path of `--opt --method
    ancopt`) pass their names through `get_namespace`, so two minimisations
    with different namespaces in one directory overwrote each other.  Without
    `--namespace` the names must stay `xtbopt.xyz`/`xtb.trj`.

    `inp` is ignored: the check writes its own *stretched* water, because two
    runs of a geometry that is already at the minimum could produce identical
    files, and an identical pair cannot tell an overwrite from a correct pair.
    The rundir is removed at entry, because the check asserts the *absence* of
    the unprefixed `xtbopt.xyz`.
    """
    rundir = os.path.abspath(rundir)
    shutil.rmtree(rundir, ignore_errors=True)
    os.makedirs(os.path.join(rundir, "plain"))
    os.makedirs(os.path.join(rundir, "two"))
    stretched = os.path.join(rundir, "stretched_h2o.xyz")
    with open(stretched, "w") as handle:
        handle.write("3\nstretched water\n"
                     "O  0.0000000  0.0000000  0.0000000\n"
                     "H  0.0000000  0.0000000  1.3000000\n"
                     "H  1.3000000  0.0000000 -0.4000000\n")

    # 1. the default names (tool-chain compatibility)
    plain = os.path.join(rundir, "plain")
    status, output = run(
        exe, stretched,
        ["--gfn", "2", "--lopt", "--lopt-level", "loose", "--norestart"],
        plain,
    )
    if status != 0:
        fail("the default --lopt run failed with status %d" % status, output)
    for name in ("xtbopt.xyz", "xtb.trj"):
        path = os.path.join(plain, name)
        if not os.path.exists(path):
            fail("without --namespace the output is still %s, but it was not "
                 "written" % name, output)
        if os.path.getsize(path) == 0:
            fail("%s is empty" % name)

    # 2. two namespaces in one directory
    two = os.path.join(rundir, "two")
    status, output = run(
        exe, stretched,
        ["--gfn", "2", "--lopt", "--lopt-level", "loose", "--norestart",
         "--namespace", "c1"],
        two,
    )
    if status != 0:
        fail("the --namespace c1 run failed with status %d" % status, output)
    first = os.path.join(two, "c1.xtbopt.xyz")
    if not os.path.exists(first):
        fail("--namespace c1 wrote no c1.xtbopt.xyz", output)
    with open(first) as handle:
        first_text = handle.read()

    status, output = run(
        exe, stretched,
        ["--gfn", "2", "--lopt", "--lopt-level", "vtight", "--norestart",
         "--namespace", "c2"],
        two,
    )
    if status != 0:
        fail("the --namespace c2 run failed with status %d" % status, output)
    for name in ("c1.xtbopt.xyz", "c2.xtbopt.xyz", "c1.xtb.trj", "c2.xtb.trj"):
        if not os.path.exists(os.path.join(two, name)):
            fail("%s was not written" % name, output)
    if os.path.exists(os.path.join(two, "xtbopt.xyz")):
        fail("a namespaced run also wrote the unprefixed xtbopt.xyz")
    with open(first) as handle:
        if handle.read() != first_text:
            fail("the --namespace c2 run overwrote c1.xtbopt.xyz")
    with open(os.path.join(two, "c2.xtbopt.xyz")) as handle:
        second_text = handle.read()
    if first_text == second_text:
        fail("the two namespaces hold the same geometry; the check cannot "
             "tell an overwrite from a correct pair")


def check_ljson_default(exe, inp, rundir):
    """The bare default run type (`--sp`) writes its document.

    Report 117 retires the old contract "`--ljson` without a landscape run type
    writes nothing and warns": the bare default `set%runtyp` is 2 = `--sp`, and
    every documented run type now writes `landscape_<task>.json`.  This mode
    mirrors the CMake script of `xtb/landscape-ljson-default`: the run exits 0,
    writes `landscape_sp.json` with the uniform core block, and no warning
    about a missing run type is printed.
    """
    rundir = os.path.abspath(rundir)
    shutil.rmtree(rundir, ignore_errors=True)
    os.makedirs(rundir)
    status, output = run(
        exe, os.path.abspath(inp), ["--gfn", "2", "--ljson"], rundir
    )
    if status != 0:
        fail("xtb --ljson terminated with status %d" % status, output)
    if "writes nothing" in output:
        fail("the retired --ljson no-op warning is still printed", output)
    document = read_json(rundir, "landscape_sp.json")
    if document.get("task") != "sp":
        fail("the default --ljson document is not an sp result", output)
    for key in ("format", "n_atoms", "energy", "geometry", "provenance"):
        if key not in document:
            fail("the default --ljson document lacks the core key %s" % key,
                 output)
    print("ljson default: the bare --sp default writes landscape_sp.json")


def check_method(exe, inp, rundir, hcn, nh3, nh3down):
    """`--method` is the canonical name of the four method flags.

    Mirrors `xtb/landscape-method` (CLI redesign section 4).  Three
    properties are asserted, each with its own falsifiable mutation:

    * **alias**: every pair `--method <name>` / `--<task>-method-flag <name>`
      has to produce a byte identical result file, and the run has to reach
      the task with that name (the milestone line), so the alias cannot be a
      silent no-op;
    * **filter**: a name that belongs to another task has to abort with a
      message that names the flag of that task (`--lopt --method dimer` ->
      `--lts`), and an unknown name or a missing run type has to abort as
      well;
    * **conflict/auto/repeat**: two different values through two different
      flags are an error, `auto` keeps the default of the task (so it is
      byte identical to the run without `--method`) and a repeated flag keeps
      the last value.

    The extra fixtures are the HCN transition state guess (`--lts`/`--lirc`),
    the ammonia pair (band search) and the 501 atom GFN-FF cluster (native
    engine of `xtb/landscape-method-native`); `inp` is the water input.
    """

    def invoke(label, input_file, args, cwd=None):
        """Run one command in its own directory, return (status, text, dir)."""
        workdir = os.path.join(rundir, label) if cwd is None else cwd
        if cwd is None:
            shutil.rmtree(workdir, ignore_errors=True)
            os.makedirs(workdir)
        proc = subprocess.run(
            [exe, input_file] + args, cwd=workdir, capture_output=True, text=False
        )
        return proc.returncode, decode(proc.stdout) + decode(proc.stderr), workdir

    def identical(label, file_a, dir_a, file_b, dir_b):
        with open(os.path.join(dir_a, file_a), "rb") as handle:
            left = handle.read()
        with open(os.path.join(dir_b, file_b), "rb") as handle:
            right = handle.read()
        if left != right:
            fail("%s: %s differs between the two routes" % (label, file_a))

    def expect_milestone(label, args, input_file, milestone, skipped=()):
        status, output, workdir = invoke(label, input_file, args)
        if status != 0:
            fail("%s: exit status %d" % (label, status), output)
        if milestone not in output:
            fail("%s: the task did not report %r" % (label, milestone), output)
        for other in skipped:
            if other in output:
                fail("%s: the run reported %r as well" % (label, other), output)
        return workdir, output

    def expect_error(label, args, input_file, pattern, second=None):
        status, output, _ = invoke(label, input_file, args)
        if status == 0:
            fail("%s: expected a non-zero exit" % label, output)
        if not re.search(pattern, output):
            fail("%s: missing %r in the output" % (label, pattern), output)
        if second is not None and not re.search(second, output):
            fail("%s: missing %r in the output" % (label, second), output)

    # --- alias: opt, scan, ts (single ended and band) and irc ---------------
    opt_method, _ = expect_milestone(
        "alias-opt-method",
        ["--gfn", "2", "--lopt", "--method", "bfgs", "--norestart"],
        inp, "landscape: bfgs minimisation",
        skipped=("landscape: lbfgs minimisation",),
    )
    opt_engine, _ = expect_milestone(
        "alias-opt-engine",
        ["--gfn", "2", "--lopt", "--lopt-engine", "bfgs", "--norestart"],
        inp, "landscape: bfgs minimisation",
    )
    identical("alias-opt", "xtbopt.xyz", opt_method, "xtbopt.xyz", opt_engine)
    opt_upper, _ = expect_milestone(
        "alias-opt-upper",
        ["--gfn", "2", "--lopt", "--method", "BFGS", "--norestart"],
        inp, "landscape: bfgs minimisation",
    )
    identical("alias-opt-upper", "xtbopt.xyz", opt_method, "xtbopt.xyz", opt_upper)

    scan_args = ["--lscan-coord", "distance", "--lscan-atoms", "1,2",
                 "--lscan-step", "0.05", "--lscan-nstep", "2", "--lscan-oneway"]
    scan_method, _ = expect_milestone(
        "alias-scan-method",
        ["--gfn", "2", "--lscan", "--method", "cg"] + scan_args,
        inp, "landscape: relaxed scan",
    )
    scan_flag, _ = expect_milestone(
        "alias-scan-flag",
        ["--gfn", "2", "--lscan", "--lscan-method", "cg"] + scan_args,
        inp, "landscape: relaxed scan",
    )
    identical("alias-scan", "xtbscan.log", scan_method, "xtbscan.log", scan_flag)

    ts_common = ["--gfn", "2", "--lts-ninit", "force", "--norestart"]
    ts_method, _ = expect_milestone(
        "alias-ts-method", ["--lts", "--method", "dimer"] + ts_common,
        hcn, "landscape: converged",
    )
    ts_flag, _ = expect_milestone(
        "alias-ts-flag", ["--lts", "--lts-method", "dimer"] + ts_common,
        hcn, "landscape: converged",
    )
    identical("alias-ts", "xtbts.xyz", ts_method, "xtbts.xyz", ts_flag)

    # The climbing image side effect has to live in the shared setter: both
    # routes ignite the same image and write the same band result.
    band_common = ["--gfn", "2", "--lts-images", "5",
                   "--lts-product", nh3down, "--norestart", "--namespace", "p1"]
    cineb_method, cineb_output = expect_milestone(
        "alias-cineb-method", ["--lts", "--method", "cineb"] + band_common,
        nh3, "landscape: converged",
    )
    cineb_flag, cineb_flag_output = expect_milestone(
        "alias-cineb-flag", ["--lts", "--lts-method", "cineb"] + band_common,
        nh3, "landscape: converged",
    )
    for label, output in (("alias-cineb-method", cineb_output),
                          ("alias-cineb-flag", cineb_flag_output)):
        if not re.search(r"climbing image [0-9]+ ignited at iteration", output):
            fail("%s: the climbing image was never ignited" % label, output)
        if not re.search(r"climbing image [(]TS candidate[)]: [0-9]+", output):
            fail("%s: the climbing image was not reported" % label, output)
    identical("alias-cineb", "p1.path.xyz", cineb_method, "p1.path.xyz", cineb_flag)
    identical("alias-cineb", "p1.xtbts.xyz", cineb_method, "p1.xtbts.xyz", cineb_flag)

    # The IRC needs a converged transition state as input.
    status, output, workdir = invoke(
        "irc-prepare", hcn,
        ["--gfn", "2", "--lts", "--lts-method", "prfo", "--norestart",
         "--namespace", "irc0"],
    )
    if status != 0:
        fail("the --lirc preparation run failed with status %d" % status, output)
    ts_file = os.path.join(workdir, "irc0.xtbts.xyz")
    if not os.path.exists(ts_file):
        fail("the --lirc preparation run wrote no xtbts.xyz", output)
    irc_common = ["--gfn", "2", "--lirc-max-steps", "2", "--norestart",
                  "--namespace", "irc1"]
    irc_status_m, irc_output_m, irc_dir_m = invoke(
        "alias-irc-method", ts_file, ["--lirc", "--method", "gs"] + irc_common
    )
    irc_status_f, irc_output_f, irc_dir_f = invoke(
        "alias-irc-flag", ts_file, ["--lirc", "--lirc-method", "gs"] + irc_common
    )
    if irc_status_m != 0 or irc_status_f != 0:
        fail("the two --lirc routes exited %d and %d"
             % (irc_status_m, irc_status_f), irc_output_m + irc_output_f)
    identical("alias-irc", "irc1.irc.log", irc_dir_m, "irc1.irc.log", irc_dir_f)

    # --- filter: a name of another task names the flag of that task ---------
    expect_error("filter-lopt-dimer",
                 ["--gfn", "2", "--lopt", "--method", "dimer", "--norestart"],
                 inp, r"dimer is a transition-state method", r"--lts")
    expect_error("filter-lopt-hpc",
                 ["--gfn", "2", "--lopt", "--method", "hpc", "--norestart"],
                 inp, r"intrinsic reaction coordinate method", r"--lirc")
    expect_error("filter-lscan-dimer",
                 ["--gfn", "2", "--lscan", "--method", "dimer"],
                 inp, r"dimer is a transition-state method", r"--lts")
    expect_error("filter-lts-bfgs",
                 ["--gfn", "2", "--lts", "--method", "bfgs", "--norestart"],
                 inp, r"bfgs is a minimisation method", r"--lopt")
    expect_error("filter-lirc-bfgs",
                 ["--gfn", "2", "--lirc", "--method", "bfgs"],
                 inp, r"bfgs is a minimisation method", r"--lopt")
    expect_error("filter-lopt-native",
                 ["--gfn", "2", "--lopt", "--method", "ancopt", "--norestart"],
                 inp, r"only available with --opt or --ohess", r"--opt")
    expect_error("filter-lvib",
                 ["--gfn", "2", "--lvib", "--method", "bfgs"],
                 inp, r"cannot be used with --lvib", r"--lvib-hessian")
    expect_error("filter-sp",
                 ["--gfn", "2", "--sp", "--method", "bfgs"],
                 inp, r"needs a landscape run type")
    expect_error("filter-none",
                 ["--gfn", "2", "--method", "bfgs"],
                 inp, r"needs a landscape run type")
    expect_error("filter-unknown",
                 ["--gfn", "2", "--lopt", "--method", "bogus", "--norestart"],
                 inp, r"unknown landscape method 'bogus'")
    expect_error("filter-method-missing",
                 ["--gfn", "2", "--lopt", "--method"],
                 inp, r"requires a method name")
    expect_error("filter-native-none",
                 ["--gfn", "2", "--method", "ancopt", "--norestart"],
                 inp, r"only available with --opt or --ohess")

    # --- conflict: two different values through two different flags --------
    for label, args in (
        ("conflict-opt", ["--gfn", "2", "--lopt", "--method", "bfgs",
                          "--lopt-engine", "cg", "--norestart"]),
        ("conflict-opt-reverse", ["--gfn", "2", "--lopt", "--lopt-engine", "cg",
                                  "--method", "bfgs", "--norestart"]),
        ("conflict-scan", ["--gfn", "2", "--lscan", "--method", "cg",
                           "--lscan-method", "bfgs"]),
        ("conflict-ts", ["--gfn", "2", "--lts", "--method", "dimer",
                         "--lts-method", "prfo", "--norestart"]),
        ("conflict-irc", ["--gfn", "2", "--lirc", "--method", "gs",
                          "--lirc-method", "lqa"]),
    ):
        expect_error(label, args, inp, r"conflicting landscape methods", r"--method")

    # same value through both flags is idempotent
    status, output, workdir = invoke(
        "conflict-same",
        inp,
        ["--gfn", "2", "--lopt", "--method", "bfgs", "--lopt-engine", "bfgs",
         "--norestart"],
    )
    if status != 0:
        fail("the same value through both flags was rejected", output)
    if "landscape: bfgs minimisation" not in output:
        fail("the idempotent pair did not reach the task", output)

    # --- auto: the default of the selected task, a no-op -------------------
    plain, _ = expect_milestone(
        "auto-plain", ["--gfn", "2", "--lopt", "--norestart"], inp,
        "landscape: lbfgs minimisation",
    )
    auto, auto_output = expect_milestone(
        "auto-method", ["--gfn", "2", "--lopt", "--method", "auto", "--norestart"],
        inp, "landscape: lbfgs minimisation",
    )
    identical("auto-lopt", "xtbopt.xyz", plain, "xtbopt.xyz", auto)
    if "Engine: lbfgs (default)" not in auto_output:
        fail("`auto` did not report the builtin default engine", auto_output)
    explicit, _ = expect_milestone(
        "auto-explicit",
        ["--gfn", "2", "--lopt", "--lopt-engine", "bfgs", "--norestart"],
        inp, "landscape: bfgs minimisation",
    )
    auto_explicit, auto_explicit_output = expect_milestone(
        "auto-explicit-auto",
        ["--gfn", "2", "--lopt", "--lopt-engine", "bfgs", "--method", "auto",
         "--norestart"],
        inp, "landscape: bfgs minimisation",
    )
    identical("auto-explicit", "xtbopt.xyz", explicit, "xtbopt.xyz", auto_explicit)
    if "Engine: bfgs (explicit)" not in auto_explicit_output:
        fail("`auto` overrode the explicit engine", auto_explicit_output)
    ts_auto, _ = expect_milestone(
        "auto-ts", ["--gfn", "2", "--lts", "--method", "auto",
                    "--lts-ninit", "force", "--norestart"],
        hcn, "landscape: converged",
    )
    ts_plain, _ = expect_milestone(
        "auto-ts-plain", ["--gfn", "2", "--lts", "--lts-ninit", "force",
                          "--norestart"],
        hcn, "landscape: converged",
    )
    identical("auto-ts", "xtbts.xyz", ts_auto, "xtbts.xyz", ts_plain)

    # --- a repeated flag keeps the last value ------------------------------
    cg, cg_output = expect_milestone(
        "repeat-cg", ["--gfn", "2", "--lopt", "--method", "cg", "--norestart"],
        inp, "landscape: cg minimisation",
    )
    last, last_output = expect_milestone(
        "repeat-last",
        ["--gfn", "2", "--lopt", "--method", "bfgs", "--method", "cg",
         "--norestart"],
        inp, "landscape: cg minimisation",
        skipped=("landscape: bfgs minimisation",),
    )
    identical("repeat-method", "xtbopt.xyz", cg, "xtbopt.xyz", last)
    engine_last, _ = expect_milestone(
        "repeat-engine-last",
        ["--gfn", "2", "--lopt", "--lopt-engine", "bfgs", "--lopt-engine", "cg",
         "--norestart"],
        inp, "landscape: cg minimisation",
        skipped=("landscape: bfgs minimisation",),
    )
    identical("repeat-engine", "xtbopt.xyz", cg, "xtbopt.xyz", engine_last)

    # --- a repeated band flag is equivalent to a repeated `--method` --------
    # (stage 2 review L3)  The climbing image stage that `cineb` implies is a
    # side effect of the *final* method, not of every flag: otherwise
    # `--lts-method cineb --lts-method neb` kept the climbing stage while the
    # equivalent `--method cineb --method neb` did not.
    repeat_band = ["--gfn", "2", "--lts-images", "5",
                   "--lts-product", nh3down, "--norestart", "--namespace", "p2"]
    band_flag, band_flag_output = expect_milestone(
        "repeat-band-flag",
        ["--gfn", "2", "--lts", "--lts-method", "cineb",
         "--lts-method", "neb"] + repeat_band[2:],
        nh3, "landscape: transition state search (neb)",
    )
    band_method, band_method_output = expect_milestone(
        "repeat-band-method",
        ["--gfn", "2", "--lts", "--method", "cineb",
         "--method", "neb"] + repeat_band[2:],
        nh3, "landscape: transition state search (neb)",
    )
    for label, output in (("repeat-band-flag", band_flag_output),
                          ("repeat-band-method", band_method_output)):
        if re.search(r"climbing image [0-9]+ ignited", output):
            fail("%s: the discarded `cineb` still ignited the climbing image"
                 % label, output)
    identical("repeat-band", "p2.path.xyz", band_flag, "p2.path.xyz", band_method)
    identical("repeat-band", "p2.xtbts.xyz", band_flag, "p2.xtbts.xyz", band_method)

    # an explicit `--lts-climb` still switches the stage on
    explicit_climb = expect_milestone(
        "repeat-band-explicit-climb",
        ["--gfn", "2", "--lts", "--method", "neb", "--lts-climb"]
        + repeat_band[2:],
        nh3, "landscape: transition state search (neb)",
    )
    if not re.search(r"climbing image [0-9]+ ignited", explicit_climb[1]):
        fail("repeat-band-explicit-climb: `--lts-climb` was lost")


def check_method_native(exe, inp, rundir, big):
    """The four native engines behind `--method` (CLI redesign section 4.8).

    Mirrors `xtb/landscape-method-native`, `-ohess` and `-opt-reject`.  The
    discriminators are the measured stdout tokens of the native engines
    (`ANC micro-cycles` / `L-ANC optimizer` / `initial time step` /
    `current lattice:`); `landscape: ` appearing in the same output is the
    failure - the name has to select the native path, not the suite.

    `big` is a 501 atom GFN-FF cluster: `src/main/defaults.f90:42-46` gives
    such a system `p_engine_lbfgs` unless the engine is pinned, so it is the
    only fixture that can see the explicit assignment of `ancopt` (on a small
    system the default is `rf` as well).

    Scope note for `landscape_native_engine` (stage 2 review INFO-1, fix
    report T7): the value **is** read at runtime - the parser tail writes it
    into `set%opt_engine`, and the 501 atom fixture below is falsifiable
    against that assignment (without it the run reports `L-ANC optimizer`
    instead of `ANC micro-cycles`).  Since stage 3 it is also the dispatch
    flag of `--opt`/`--ohess`: the redeployed paths are asserted in
    `xtb/landscape-opt-alias-once`, `xtb/landscape-opt-native-exit` and the
    four `--ohess` tests of the same round, so this check keeps its stage 2
    scope (the four native exits) and no longer claims that the redeployment
    is untested.  The declaration default itself (the
    `landscape_native_engine = 0` line before the flag loop) stays a
    **static** check - `intent(out)` makes a missing initialisation read
    undefined memory that `-O2` usually happens to see as 0, which is what the
    design section 4.7.1 records.
    """

    def invoke(label, input_file, args):
        workdir = os.path.join(rundir, label)
        shutil.rmtree(workdir, ignore_errors=True)
        os.makedirs(workdir)
        proc = subprocess.run(
            [exe, input_file] + args, cwd=workdir, capture_output=True, text=False
        )
        return proc.returncode, decode(proc.stdout) + decode(proc.stderr), workdir

    for name, token in (("ancopt", "ANC micro-cycles"),
                        ("rf", "ANC micro-cycles"),
                        ("lancopt", "L-ANC optimizer"),
                        ("l-ancopt", "L-ANC optimizer"),
                        ("fire", "initial time step"),
                        ("inertial", "initial time step"),
                        ("pbc_lbfgs", "current lattice:")):
        status, output, _ = invoke(
            "native-" + name, inp,
            ["--gfn", "2", "--opt", "--method", name, "--norestart"],
        )
        if status != 0:
            fail("--opt --method %s exited %d" % (name, status), output)
        if token not in output:
            fail("--opt --method %s did not reach the native engine" % name, output)
        if "landscape: " in output:
            fail("--opt --method %s ran the landscape suite" % name, output)

    # the order of `--method` and the run type flag is irrelevant
    status, forward, dir_forward = invoke(
        "native-order-forward", inp,
        ["--gfn", "2", "--opt", "--method", "ancopt", "--norestart"],
    )
    status2, backward, dir_backward = invoke(
        "native-order-backward", inp,
        ["--gfn", "2", "--method", "ancopt", "--opt", "--norestart"],
    )
    if status != 0 or status2 != 0:
        fail("the two flag orders exited %d and %d" % (status, status2),
             forward + backward)
    for label, output in (("forward", forward), ("backward", backward)):
        if "ANC micro-cycles" not in output:
            fail("the %s order did not reach the native engine" % label, output)
    with open(os.path.join(dir_forward, "xtbopt.xyz"), "rb") as handle:
        left = handle.read()
    with open(os.path.join(dir_backward, "xtbopt.xyz"), "rb") as handle:
        right = handle.read()
    if left != right:
        fail("the two flag orders converge to different geometries")

    # an explicit `ancopt` must beat the >500 atom default `lbfgs`
    status, output, _ = invoke(
        "native-big", big,
        ["--gfnff", "--opt", "--method", "ancopt", "--norestart"],
    )
    if status != 0:
        fail("the 501 atom run exited %d" % status, output)
    if "ANC micro-cycles" not in output:
        fail("the explicit `ancopt` did not reach the native engine", output)
    if "L-ANC optimizer" in output:
        fail("the >500 atom default overrode the explicit engine", output)

    # `--ohess` dispatches the same native engines and keeps the native
    # Hessian/frequency path: `ancopt` is the default engine of a small
    # system, so the files have to be byte identical to the plain run.
    status, output, ohess_native = invoke(
        "native-ohess", inp,
        ["--gfn", "2", "--ohess", "--method", "ancopt", "--norestart"],
    )
    if status != 0:
        fail("--ohess --method ancopt exited %d" % status, output)
    if "ANC micro-cycles" not in output:
        fail("--ohess --method ancopt did not reach the native engine", output)
    if "landscape: " in output:
        fail("--ohess --method ancopt ran the landscape suite", output)
    if "Numerical Hessian" not in output:
        fail("--ohess --method ancopt lost the native Hessian block", output)
    # the plain `--ohess` is the redeployed one: it runs the suite chain and
    # its files are asserted in `xtb/landscape-ohess-suite-once`; here only
    # the dispatch (native side versus suite side) is pinned
    status, output, ohess_plain = invoke(
        "native-ohess-plain", inp,
        ["--gfn", "2", "--ohess", "--norestart"],
    )
    if status != 0:
        fail("the plain --ohess run exited %d" % status, output)
    if "landscape: harmonic frequencies" not in output:
        fail("the plain --ohess run did not reach the landscape suite", output)
    if "ANC micro-cycles" in output:
        fail("the plain --ohess run also ran the native ANCopt", output)
    for name in ("hessian", "vibspectrum"):
        if not os.path.exists(os.path.join(ohess_plain, name)):
            fail("the plain --ohess run did not write %s" % name)

    # a landscape method name under `--opt`/`--ohess` reaches the suite
    # settings object the run reads; the stage 2 refusal (R39 of the redesign)
    # is gone with the redeployment
    for label, args, milestone in (
        ("stage3-opt-suite-method", ["--gfn", "2", "--opt", "--method", "bfgs"],
         "landscape: bfgs minimisation"),
        ("stage3-ohess-suite-method", ["--gfn", "2", "--ohess", "--method", "bfgs"],
         "landscape: bfgs minimisation"),
    ):
        status, output, _ = invoke(label, inp, args)
        if status != 0:
            fail("%s: exit status %d" % (label, status), output)
        if milestone not in output:
            fail("%s: the landscape method name did not reach the suite" % label,
                 output)

    # --- U15: the command line wins over the xcontrol, but not silently ----
    # (CLI redesign section 4.7.1, fix report T10)  `$opt engine=` is ignored
    # as soon as `--method <native>` pinned the engine; the conflicting
    # request has to be reported instead of disappearing, and a request that
    # agrees with the command line must stay quiet.  Since stage 3 the block
    # has to arrive with a native exit, because a bare `--opt` plus an `$opt`
    # block is exactly what the C7 guard refuses; the two cases that used to
    # run the block alone are therefore spelled with `--method ancopt` (the
    # engine they exercise is `rf`), and the case that used to test "the
    # xcontrol selected the engine" is covered by `u15-xcontrol-twice`.
    rf_control = os.path.join(rundir, "opt_engine_rf.ctl")
    with open(rf_control, "w") as handle:
        handle.write("$opt\n   engine=rf\n")
    lbfgs_control = os.path.join(rundir, "opt_engine_lbfgs.ctl")
    with open(lbfgs_control, "w") as handle:
        handle.write("$opt\n   engine=lbfgs\n")
    same_control = os.path.join(rundir, "opt_engine_twice.ctl")
    with open(same_control, "w") as handle:
        handle.write("$opt\n   engine=rf\n$end\n$opt\n   engine=lbfgs\n")

    for label, args, warned, token in (
        ("u15-conflict", ["--gfn", "2", "--opt", "--method", "lancopt",
                          "--norestart", "--input", rf_control], True,
         "L-ANC optimizer"),
        ("u15-agree", ["--gfn", "2", "--opt", "--method", "lancopt",
                       "--norestart", "--input", lbfgs_control], False,
         "L-ANC optimizer"),
        ("u15-cli-only", ["--gfn", "2", "--opt", "--method", "lancopt",
                          "--norestart"], False, "L-ANC optimizer"),
        ("u15-xcontrol-only", ["--gfn", "2", "--opt", "--method", "ancopt",
                               "--norestart", "--input", rf_control], False,
         "ANC micro-cycles"),
        # a repeated `$opt engine=` with two different values was silently
        # ignored before as well; it is the same conflict and is reported too
        ("u15-xcontrol-twice", ["--gfn", "2", "--opt", "--method", "ancopt",
                                "--norestart", "--input", same_control], True,
         "ANC micro-cycles"),
    ):
        status, output, _ = invoke(label, inp, args)
        if status != 0:
            fail("%s: exit status %d" % (label, status), output)
        ignored = "is ignored: an optimizer engine was already selected" in output
        if ignored != warned:
            fail("%s: expected %s warning for the ignored $opt engine"
                 % (label, "a" if warned else "no"), output)
        if warned and "[WARNING]" not in output:
            fail("%s: the conflict was reported without a [WARNING] line"
                 % label, output)
        # the command line selection survives in every case
        if token not in output:
            fail("%s: the command line engine was overridden by the $opt block"
                 % label, output)


def check_native_opt_symbols(exe, inp, rundir):
    """The native `--opt` L-BFGS path must write the element symbols intact.

    Mirrors `xtb/native-opt-lbfgs-symbols`.  The L-ANCopt engine used to copy
    the geometry with an intrinsic `TMolecule` assignment
    (`src/relaxation_engine.f90:615` and `:779`); gfortran 12.2.0 corrupts
    `sym(2:)` of a `character(len=:), allocatable :: sym(:)` component under
    intrinsic assignment while every length stays correct (a stride bug in
    deferred-length character array components - the suite's own
    `src/landscape/types/log.f90` already works around it).  Measured
    consequence: with the L-ANCopt engine
    the symbol column of `xtbopt.xyz` **and** of the `final structure:` block
    of stdout came out as NUL bytes, while the numbers were correct.  The `rf`
    engine was never affected, which is why the native `--opt` (default engine
    `rf` for a small system) looked healthy in the BC-4 fixtures.

    The engine is selected with `--method lancopt` and not with an
    `$opt engine=lbfgs` xcontrol anymore: since stage 3 of the CLI redesign
    `--opt` without `--method` is a landscape run type, and an `$opt` block
    next to it is refused (C7).  The selection is the same one - `--method
    lancopt` pins `p_engine_lbfgs` in the parser, `$opt engine=lbfgs` set it in
    the control file - so the engine under test is unchanged.

    The assertion is on the **bytes**, not on the text: a NUL byte is a valid
    UTF-8 code point (U+0000), so a strict `decode("utf-8")` alone stays green
    on the corrupted output - and so does the tolerant `errors='replace'`
    decoding that stage 2 introduced for the unrelated Fortran record markers
    of the 501 atom fixture (fix report T2).  Therefore the check scans both
    streams for control bytes and compares the symbol column of the xyz file
    with the input geometry as well.
    """
    workdir = os.path.join(rundir, "native-opt-lbfgs")
    shutil.rmtree(workdir, ignore_errors=True)
    os.makedirs(workdir)

    proc = subprocess.run(
        [exe, inp, "--gfn", "2", "--opt", "normal", "--method", "lancopt",
         "--norestart"],
        cwd=workdir, capture_output=True, text=False,
    )
    if proc.returncode != 0:
        fail("the native L-BFGS `--opt` exited %d" % proc.returncode,
             decode(proc.stdout) + decode(proc.stderr))
    output = decode(proc.stdout) + decode(proc.stderr)
    if "L-ANC optimizer" not in output:
        fail("the L-BFGS run did not reach the L-ANCopt engine", output)

    # the input fixture is the reference for the atom count and the symbols
    with open(inp, "rb") as handle:
        reference = handle.read().decode("utf-8").splitlines()
    natoms = int(reference[0].split()[0])
    expected = [line.split()[0] for line in reference[2:2 + natoms]]

    with open(os.path.join(workdir, "xtbopt.xyz"), "rb") as handle:
        raw_xyz = handle.read()
    raw_stdout = proc.stdout + proc.stderr

    for label, raw in (("stdout", raw_stdout), ("xtbopt.xyz", raw_xyz)):
        try:
            raw.decode("utf-8")
        except UnicodeDecodeError as error:
            fail("%s is not valid UTF-8: %s" % (label, error))
        for offset, byte in enumerate(raw):
            if byte < 0x20 and byte not in (0x09, 0x0a, 0x0d):
                fail("%s carries the control byte 0x%02x at offset %d"
                     % (label, byte, offset))

    # the symbol column of the written geometry
    xyz = raw_xyz.decode("utf-8").splitlines()
    if int(xyz[0].split()[0]) != natoms:
        fail("xtbopt.xyz reports %s atoms instead of %d"
             % (xyz[0].split()[0], natoms))
    written = [line.split()[0] for line in xyz[2:2 + natoms]]
    if written != expected:
        fail("xtbopt.xyz carries the symbols %r instead of %r"
             % (written, expected))

    # the same symbols in the `final structure:` block of stdout
    block = re.search(r"final structure:\s*\n=+\n\d+\n[^\n]*\n"
                      r"((?:\S+\s+[-0-9.]+\s+[-0-9.]+\s+[-0-9.]+\n)+)",
                      decode(raw_stdout))
    if not block:
        fail("the `final structure:` block is missing from stdout", output)
    printed = [line.split()[0] for line in block.group(1).splitlines()]
    if printed != expected:
        fail("the `final structure:` block carries the symbols %r instead "
             "of %r" % (printed, expected))


def check_method_table(exe, inp, rundir, hcn, nh3, nh3product):
    """Every suite method name of `--method` has to reach its task (M1).

    Mirrors `xtb/landscape-method-table`.  `check_method` (mode `method`) pins
    the alias/conflict/`auto` grammar and runs 6 of the 18 suite names; this
    table runs **all 18** plus the four `auto` runs, so a name that is accepted
    by the parser and then never reaches the task cannot stay green.

    The assertion per name is deliberately two-sided where the task writes the
    resolved name into a file:

      * `opt`  (6): stdout milestone `landscape: <name> minimisation` + a
        non-empty `xtbopt.xyz`.  No product of this task carries the name, so
        the file check only proves the task ran; the milestone is the
        dispatch assertion (measured: `sd` and `cg` coincide on GFN2 water by
        design of the CG restart, so a trajectory comparison cannot separate
        all six names).
      * `ts`   (8): milestone `landscape: transition state search (<name>)`
        **and** `<ns>.xtbts.log`, whose first line is written from the resolved
        method (`src/landscape/tasks/ts.f90:685`).
      * `irc`  (4): milestone **and** `<ns>.irc.log`
        (`src/landscape/tasks/irc.f90:996`), which carries the same information.

    Failures are collected and reported together: one mutation of the name
    table has to name every affected name, not only the first one.

    `hcn` is a transition state guess (`--lts`/`--lirc` parent), `nh3` with
    `nh3product` is the band fixture, `inp` is the water input.
    """
    misses = []

    def invoke(label, input_file, args):
        workdir = os.path.join(rundir, label)
        shutil.rmtree(workdir, ignore_errors=True)
        os.makedirs(workdir)
        proc = subprocess.run(
            [exe, input_file] + args, cwd=workdir, capture_output=True, text=False
        )
        return proc.returncode, decode(proc.stdout) + decode(proc.stderr), workdir

    def inspect(label, input_file, args, milestone, product=None, token=None):
        status, output, workdir = invoke(label, input_file, args)
        if status != 0:
            misses.append("%s: exit status %d" % (label, status))
        if milestone is not None and milestone not in output:
            misses.append("%s: the run did not report %r" % (label, milestone))
        if product is not None:
            path = os.path.join(workdir, product)
            if not os.path.exists(path):
                misses.append("%s: no %s was written" % (label, product))
            else:
                with open(path, "rb") as handle:
                    body = decode(handle.read())
                if token not in body:
                    misses.append("%s: %s does not report %r" % (label, product, token))
        return workdir, output

    # --- the six minimisers, on the water fixture --------------------------
    for name in ("sd", "cg", "lbfgs", "bfgs", "diis", "rfo"):
        _, output = inspect(
            "table-opt-" + name, inp,
            ["--gfn", "2", "--lopt", "--method", name, "--norestart"],
            "landscape: %s minimisation" % name,
        )
        # the engine banner says `(explicit)` only when a name was given
        if "Engine: %s (explicit)" % name not in output:
            misses.append("table-opt-%s: the name did not count as an explicit engine" % name)

    # --- the eight saddle searches -----------------------------------------
    # single ended: P-RFO / B-PRFO / dimer
    for name in ("prfo", "bprfo", "dimer"):
        inspect(
            "table-ts-" + name, hcn,
            ["--gfn", "2", "--lts", "--method", name, "--lts-ninit", "force",
             "--norestart", "--namespace", "nm"],
            "landscape: transition state search (%s)" % name,
            product="nm.xtbts.log",
            token="transition state search: %s" % name,
        )
    # chain of states: NEB / CI-NEB / AutoNEB (band) and string / GSM (growing
    # string).  The growing string ignores `--lts-images` (it is configured by
    # `--lts-string-nnodes`, `settings.f90:271-275`), so it gets its own
    # argument list.
    band = ["--lts-images", "5", "--lts-product", nh3product,
            "--norestart", "--namespace", "nm"]
    string = ["--lts-string-nnodes", "5", "--lts-product", nh3product,
              "--norestart", "--namespace", "nm"]
    # the climbing stage of the family: the band ignites a climbing image, the
    # string refines a climbing node; `neb` must not climb at all
    climbing = {
        "neb": None,
        "cineb": r"climbing image [0-9]+ ignited at iteration",
        "autoneb": r"climbing image [0-9]+ ignited at iteration",
        "string": r"climbing image [(]TS candidate[)]: [0-9]+",
        "gsm": r"climbing image [(]TS candidate[)]: [0-9]+",
    }
    for name in ("neb", "cineb", "autoneb"):
        _, output = inspect(
            "table-ts-" + name, nh3,
            ["--gfn", "2", "--lts", "--method", name] + band,
            "landscape: transition state search (%s)" % name,
            product="nm.xtbts.log",
            token="transition state search: %s" % name,
        )
        found = re.search(climbing[name], output) if climbing[name] else None
        if climbing[name] and found is None:
            misses.append("table-ts-%s: the climbing stage was never reported" % name)
        if climbing[name] is None and found is not None:
            misses.append("table-ts-%s: `neb` ignited a climbing image" % name)
    for name in ("string", "gsm"):
        _, output = inspect(
            "table-ts-" + name, nh3,
            ["--gfn", "2", "--lts", "--method", name] + string,
            "landscape: transition state search (%s)" % name,
            product="nm.xtbts.log",
            token="transition state search: %s" % name,
        )
        if re.search(climbing[name], output) is None:
            misses.append("table-ts-%s: the climbing stage was never reported" % name)

    # --- the four reaction path integrators --------------------------------
    status, output, workdir = invoke(
        "table-irc-prepare", hcn,
        ["--gfn", "2", "--lts", "--lts-method", "prfo", "--norestart",
         "--namespace", "irc0"],
    )
    ts_file = os.path.join(workdir, "irc0.xtbts.xyz")
    if status != 0 or not os.path.exists(ts_file):
        fail("the IRC preparation run failed with status %d" % status, output)
    for name in ("gs", "lqa", "hpc", "eulerpc"):
        inspect(
            "table-irc-" + name, ts_file,
            ["--gfn", "2", "--lirc", "--method", name,
             "--lirc-max-steps", "2", "--norestart", "--namespace", "nm"],
            "landscape: intrinsic reaction coordinate (%s)" % name,
            product="nm.irc.log",
            token="intrinsic reaction coordinate: %s" % name,
        )

    # --- `auto` is the default of every task -------------------------------
    # (stage 2 review section 9 listed `--lscan --method auto` and
    # `--lirc --method auto` as uncovered; `lopt`/`lts` are in mode `method`)
    scan_args = ["--lscan-coord", "distance", "--lscan-atoms", "1,2",
                 "--lscan-step", "0.05", "--lscan-nstep", "2",
                 "--lscan-oneway"]
    scan_plain, _ = inspect(
        "table-scan-auto-plain", inp, ["--gfn", "2", "--lscan"] + scan_args,
        "landscape: relaxed scan",
    )
    scan_auto, _ = inspect(
        "table-scan-auto", inp,
        ["--gfn", "2", "--lscan", "--method", "auto"] + scan_args,
        "landscape: relaxed scan",
    )
    with open(os.path.join(scan_plain, "xtbscan.log"), "rb") as handle:
        plain_log = handle.read()
    with open(os.path.join(scan_auto, "xtbscan.log"), "rb") as handle:
        auto_log = handle.read()
    if plain_log != auto_log:
        misses.append("table-scan-auto: `--method auto` changed xtbscan.log")
    irc_plain, _ = inspect(
        "table-irc-auto-plain", ts_file,
        ["--gfn", "2", "--lirc", "--lirc-max-steps", "2", "--norestart",
         "--namespace", "nm"],
        None,
    )
    irc_auto, _ = inspect(
        "table-irc-auto", ts_file,
        ["--gfn", "2", "--lirc", "--method", "auto", "--lirc-max-steps", "2",
         "--norestart", "--namespace", "nm"],
        None,
    )
    with open(os.path.join(irc_plain, "nm.irc.log"), "rb") as handle:
        plain_irc = handle.read()
    with open(os.path.join(irc_auto, "nm.irc.log"), "rb") as handle:
        auto_irc = handle.read()
    if plain_irc != auto_irc:
        misses.append("table-irc-auto: `--method auto` changed irc.log")

    # --- the inner minimiser of a relaxed scan is observable ---------------
    # (stage 2 review L2 / fix report T8)  `landscape_scan_settings%method` is
    # handed to the constraint relaxation, but no output file carries it and
    # the milestone above is mode only, so a name that was stored and then
    # ignored was invisible.  The task now echoes the effective inner
    # minimiser, and two different inner minimisers have to change
    # `xtbscan.log`, which is the part the echo alone cannot prove.  The
    # fixture is written here so that the check does not depend on a fixture
    # of the build system; it is the stretched water of
    # `LANDSCAPE-H2O-STRETCH` in `test/unit/CMakeLists.txt`.
    stretch = os.path.join(rundir, "stretch.xyz")
    with open(stretch, "w") as handle:
        handle.write("3\nstretched water\n"
                     "O  0.0000000  0.0000000  0.0000000\n"
                     "H  0.0000000  0.0000000  1.3000000\n"
                     "H  1.3000000  0.0000000 -0.4000000\n")
    inner_log = {}
    for name in ("lbfgs", "cg", "sd"):
        workdir, output = inspect(
            "table-scan-inner-" + name, stretch,
            ["--gfn", "2", "--lscan", "--method", name] + scan_args,
            "landscape: relaxed scan",
        )
        banner = re.search(r"inner minimiser: +(\S+)", output)
        if banner is None or banner.group(1) != name:
            misses.append("table-scan-inner-%s: the task did not echo the "
                          "effective inner minimiser" % name)
        with open(os.path.join(workdir, "xtbscan.log"), "rb") as handle:
            inner_log[name] = handle.read()
    for name in ("cg", "sd"):
        if inner_log[name] == inner_log["lbfgs"]:
            misses.append("table-scan-inner-%s: the inner minimiser did not "
                          "change xtbscan.log" % name)
    # a rigid scan builds the geometry and has no inner minimiser
    _, rigid_output = inspect(
        "table-scan-rigid", stretch,
        ["--gfn", "2", "--lscan", "--lscan-mode", "rigid"] + scan_args,
        "landscape: rigid scan",
    )
    if "inner minimiser" in rigid_output:
        misses.append("table-scan-rigid: a rigid scan echoed an inner minimiser")

    if misses:
        fail("%d of the method names did not reach their task:\n  %s"
             % (len(misses), "\n  ".join(misses)))


def check_error(exe, inp, rundir, pattern, args):
    """An invalid request must abort with a non-zero status and a message."""
    status, output = run(exe, inp, args, rundir)
    if status == 0:
        fail("expected a non-zero exit for '%s'" % " ".join(args), output)
    if not re.search(pattern, output):
        fail("missing %r in the output of '%s'" % (pattern, " ".join(args)), output)


#: The dangerous builtin defaults that a suite task echoes as a `[NOTICE]`
#: line (CLI redesign sections 8.3 and 17.3).  Every run names the milestone
#: line that proves the task reached its print block, whether the banner has
#: to be there, and optionally the `--ljson` document that still has to parse.
#:
#: `banner` carries the **value** the task prints next to the flag name, not
#: only the flag prefix: the notice is documented as "the value read from the
#: settings object", so a hard coded literal (or a wrong enum name) has to
#: make the default route fail.  `absent` is the value free prefix and is used
#: for the suppression direction - a value carrying pattern would miss a
#: notice that prints the *explicit* value (the bug that check has to catch).
#: The milestone line is asserted in every route, so a task that aborted
#: before its print block cannot look like a correctly suppressed notice.
ECHO_CASES = {
    "fc": {
        "banner": "[NOTICE] --lscan-fc default 0.50 hartree/bohr^2",
        "absent": "[NOTICE] --lscan-fc default",
        "runs": [
            {
                "label": "default",
                "milestone": "landscape: relaxed scan",
                "banner": True,
                "json": "landscape_scan.json",
                "args": ["--gfn", "2", "--lscan", "--lscan-coord",
                         "distance", "--lscan-atoms", "1,2",
                         "--lscan-step", "0.1", "--lscan-nstep", "2",
                         "--norestart", "--ljson"],
            },
            {
                "label": "explicit",
                "milestone": "landscape: relaxed scan",
                "banner": False,
                "json": "landscape_scan.json",
                "args": ["--gfn", "2", "--lscan", "--lscan-coord",
                         "distance", "--lscan-atoms", "1,2",
                         "--lscan-step", "0.1", "--lscan-nstep", "2",
                         "--lscan-fc", "0.05", "--norestart", "--ljson"],
            },
            {
                "label": "rigid",
                "milestone": "landscape: rigid scan",
                "banner": False,
                "args": ["--gfn", "2", "--lscan", "--lscan-mode", "rigid",
                         "--lscan-coord", "distance", "--lscan-atoms",
                         "1,2", "--lscan-step", "0.1", "--lscan-nstep",
                         "2", "--norestart"],
            },
        ],
    },
    "lopt": {
        "banner": "Engine: lbfgs (default)",
        "absent": "Engine: lbfgs (",
        "runs": [
            {
                "label": "default",
                "milestone": "landscape: lbfgs minimisation",
                "banner": True,
                "json": "landscape_opt.json",
                "json_task": "opt",
                "args": ["--gfn", "2", "--lopt", "--norestart",
                         "--ljson"],
            },
            {
                "label": "explicit",
                "milestone": "landscape: cg minimisation",
                "banner": False,
                "args": ["--gfn", "2", "--lopt", "--lopt-engine", "cg",
                         "--norestart"],
            },
        ],
    },
    "lirc": {
        "banner": "IRC integrator: hpc (default)",
        "absent": "IRC integrator: hpc (",
        "runs": [
            {
                "label": "default",
                "milestone": "landscape: intrinsic reaction coordinate (",
                "banner": True,
                "require": ["Hessian source:           numerical"],
                "args": ["--gfn", "2", "--lirc", "--lirc-max-steps", "1",
                         "--norestart"],
            },
            {
                "label": "explicit",
                "milestone": "landscape: intrinsic reaction coordinate (",
                "banner": False,
                "require": ["Hessian source:           numerical"],
                "args": ["--gfn", "2", "--lirc", "--lirc-method", "gs",
                         "--lirc-max-steps", "1", "--norestart"],
            },
            {
                "label": "hessian",
                "milestone": "landscape: intrinsic reaction coordinate (",
                "banner": True,
                "require": ["Hessian source:           model"],
                "args": ["--gfn", "2", "--lirc", "--lirc-hessian", "model",
                         "--lirc-max-steps", "1", "--norestart"],
            },
        ],
    },
    "lts": {
        "banner": "[NOTICE] --lts-escape-step default 0.20 bohr",
        "absent": "[NOTICE] --lts-escape-step default",
        "runs": [
            {
                "label": "default",
                "milestone": "landscape: transition state search (",
                "banner": True,
                "args": ["--gfn", "2", "--lts", "--lts-method", "prfo",
                         "--lts-max-steps", "1", "--norestart"],
            },
            {
                "label": "explicit",
                "milestone": "landscape: transition state search (",
                "banner": False,
                "args": ["--gfn", "2", "--lts", "--lts-method", "prfo",
                         "--lts-escape-step", "0.3", "--lts-max-steps", "1",
                         "--norestart"],
            },
        ],
    },
}


def check_echo(exe, inp, rundir, case):
    """The dangerous builtin defaults are echoed, the explicit values not.

    Mirrors `xtb/landscape-fc-echo`, `xtb/landscape-lirc-default-echo`,
    `xtb/landscape-lts-escape-echo` and `xtb/landscape-lopt-engine-echo` (CLI
    redesign sections 8.3 and 17.3).  Every route of the case is run once and
    the milestone line is asserted in **all** of them, so a task that aborted
    before its print block cannot look like a correctly suppressed notice.
    The `--ljson` routes parse the document as well: the notice goes to
    stdout, the structured result is a file, and the check pins that the two
    do not interfere.

    The two engine cases (`lopt`, `lirc`) no longer read a `[NOTICE]` line:
    their default/explicit announcement is the startup engine banner
    (`xtb_prog_main::landscape_engine_banner`), so the route asserts the box
    text instead.

    A route may add `require`: patterns that have to appear in the output of
    that route.  The `lirc` case uses it for the effective Hessian source
    (`--lirc-hessian`), an effective value rather than a dangerous default:
    it is printed for the builtin default (`numerical`) and for the explicit
    spelling alike, and a hard coded literal breaks one of the two.
    """
    spec = ECHO_CASES.get(case)
    if spec is None:
        fail("unknown echo case %r" % case)
    for route in spec["runs"]:
        label = route["label"]
        route_dir = os.path.join(rundir, case + "-" + label)
        shutil.rmtree(route_dir, ignore_errors=True)
        status, output = run(exe, inp, route["args"], route_dir)
        if route["milestone"] not in output:
            fail("%s (%s): the task did not reach its print block"
                 % (case, label), output)
        present = spec["banner"] in output
        if route["banner"] and not present:
            fail("%s (%s): the builtin default was not echoed"
                 % (case, label), output)
        if not route["banner"] and spec["absent"] in output:
            fail("%s (%s): the notice was printed although the flag was "
                 "given" % (case, label), output)
        for pattern in route.get("require", []):
            if pattern not in output:
                fail("%s (%s): %r is missing from the output"
                     % (case, label, pattern), output)
        if "json" in route:
            document = read_json(route_dir, route["json"])
            if document.get("task") != route.get("json_task", "scan"):
                fail("%s (%s): the document is not a %s result"
                     % (case, label, route.get("json_task", "scan")), output)


# ---------------------------------------------------------------------------
# CLI redesign stage 2b: the flags of section 6.3 that had no command line
# ---------------------------------------------------------------------------
#
# The design lists a block of "unified name" flags (`--lopt-max-steps`,
# `--l*-hess-step`, `--l*-sccacc`, `--lopt/--lvib-restart`).  Stage 2 found
# that none of them had landed (finding R40) and moved the block to stage 2b.
# Every one of them needs a test that proves it takes effect, not one that
# proves it parses: the assertions below all compare **observable results** of
# two runs (the reported energy or budget, or the bytes of an output file).


def knobs_run(exe, inp, args, rundir, label):
    """Run one flag combination in its own directory.

    Every combination gets a fresh directory, so the output files of two runs
    can be compared byte by byte without a `--namespace`.  Returns the exit
    status, the combined output and the directory.
    """
    workdir = os.path.join(rundir, label)
    shutil.rmtree(workdir, ignore_errors=True)
    os.makedirs(workdir)
    status, output = run(exe, inp, args, workdir)
    return status, output, workdir


def knobs_file(workdir, name):
    """Read one output file of a `knobs_run` directory as bytes."""
    path = os.path.join(workdir, name)
    if not os.path.exists(path):
        fail("%s was not written" % path)
    with open(path, "rb") as handle:
        return handle.read()


def knobs_final_energy(output, label):
    """The `final energy:` value of an `--lopt` run."""
    match = re.search(r"final energy: +(-?[0-9.]+) Eh", output)
    if match is None:
        fail("%s: the run does not report a final energy" % label, output)
    return match.group(1)


def check_opt_knobs(exe, inp, rundir):
    """The four new `--lopt` knobs reach the minimisation.

    Mirrors `xtb/landscape-opt-knobs` (CLI redesign sections 6.3 and 9.2):

      * `--lopt-max-steps 1` exhausts the budget in the first macro iteration,
        so the run echoes the budget, stays unconverged, exits non-zero and
        stops at a different energy than the default run;
      * `--lopt-hess-step` and `--lopt-sccacc` act on the numerical Hessian,
        which the default L-BFGS never builds (measured: the default run's
        `landscape engine: Hessians` counter is 0), so they are checked with
        `--lopt-engine rfo --lopt-hessian numerical` and have to move the
        final energy;
      * `--lopt-restart` reaches `landscape_engine%apply_restart`: the engine
        echo line is printed from the control object the task handed to the
        engine, and `pin`/`warm` change the trajectory (byte comparison of
        `xtb.trj`) against the default `off`.
    """
    common = ["--gfn", "2", "--norestart"]
    opt = common + ["--lopt"]
    status, output, default = knobs_run(exe, inp, opt, rundir, "optknobs-default")
    if status != 0:
        fail("the default minimisation failed with status %d" % status, output)
    # Since stage 3c the default budget is the one the default `normal` level
    # derives for the molecule: `clamp(n_atoms * 3, 200, 10000)` - 200 for the
    # three atom water fixture, where the lower clamp binds (the fixed default
    # of 256 this test pinned before is gone, CLI redesign section 3.5/R5).
    if not re.search(r"maximum number of steps: +200\s*$", output, re.M):
        fail("the default step budget is not clamp(3n, 200, 10000) = 200",
             output)
    reference = knobs_final_energy(output, "optknobs-default")

    status, output, budget = knobs_run(
        exe, inp, opt + ["--lopt-max-steps", "1"], rundir, "optknobs-max-steps")
    if status == 0:
        fail("a minimisation that exhausts its budget exited 0", output)
    if not re.search(r"maximum number of steps: +1\s*$", output, re.M):
        fail("--lopt-max-steps 1 did not reach the minimisation", output)
    if not re.search(r"^ steps: +1\s*$", output, re.M):
        fail("the run did not stop after one step", output)
    if not re.search(r"converged: +F", output):
        fail("the run claims convergence after one step", output)
    if "maximum number of steps reached" not in output:
        fail("the exhausted budget was not reported", output)
    if knobs_final_energy(output, "optknobs-max-steps") == reference:
        fail("the step budget did not change the result", output)

    hessian_alg = opt + ["--lopt-engine", "rfo", "--lopt-hessian",
                         "numerical", "--lopt-max-steps", "2"]
    status, output, rfo = knobs_run(exe, inp, hessian_alg, rundir, "optknobs-rfo")
    if status != 0:
        fail("the RFO minimisation failed with status %d" % status, output)
    rfo_energy = knobs_final_energy(output, "optknobs-rfo")
    status, output, workdir = knobs_run(
        exe, inp, hessian_alg + ["--lopt-hess-step", "0.2"], rundir,
        "optknobs-hess-step")
    if status != 0:
        fail("the RFO minimisation with --lopt-hess-step failed", output)
    if knobs_final_energy(output, "optknobs-hess-step") == rfo_energy:
        fail("--lopt-hess-step did not change the minimisation", output)
    status, output, workdir = knobs_run(
        exe, inp, hessian_alg + ["--lopt-sccacc", "10.0"], rundir,
        "optknobs-sccacc")
    if status != 0:
        fail("the RFO minimisation with --lopt-sccacc failed", output)
    if knobs_final_energy(output, "optknobs-sccacc") == rfo_energy:
        fail("--lopt-sccacc did not change the minimisation", output)

    trajectories = {}
    for mode in ("off", "pin", "warm"):
        status, output, workdir = knobs_run(
            exe, inp, opt + ["--lopt-restart", mode], rundir,
            "optknobs-restart-" + mode)
        if status != 0:
            fail("--lopt-restart %s failed with status %d" % (mode, status),
                 output)
        if not re.search(r"landscape engine: restart policy +%s\s*$" % mode,
                         output, re.M):
            fail("the engine did not read back the %s policy" % mode, output)
        trajectories[mode] = knobs_file(workdir, "xtb.trj")
    for mode in ("pin", "warm"):
        if trajectories[mode] == trajectories["off"]:
            fail("--lopt-restart %s did not change the trajectory" % mode)


def check_scan_knobs(exe, inp, rundir):
    """The two new `--lscan` knobs reach the scan points.

    Mirrors `xtb/landscape-scan-knobs` (CLI redesign sections 6.3 and 9.2):
    `--lscan-hess-step` and `--lscan-sccacc` act on the numerical Hessian of
    the inner relaxation, so the scan runs with `--lscan-method rfo
    --lscan-hessian numerical`, and each knob has to change the byte image of
    `xtbscan.log` (the de-biased energies and the drift report of the three
    points of a one step bidirectional sweep).
    """
    scan = ["--gfn", "2", "--norestart", "--lscan", "--lscan-method", "rfo",
            "--lscan-hessian", "numerical", "--lscan-coord", "angle",
            "--lscan-atoms", "2,1,3", "--lscan-step", "3",
            "--lscan-nstep", "1", "--lscan-max-steps", "2"]
    logs = {}
    for label, extra in (("default", []),
                         ("hess-step", ["--lscan-hess-step", "0.2"]),
                         ("sccacc", ["--lscan-sccacc", "10.0"])):
        status, output, workdir = knobs_run(
            exe, inp, scan + extra, rundir, "scanknobs-" + label)
        if status != 0:
            fail("the scan with %s failed with status %d" % (label, status),
                 output)
        logs[label] = knobs_file(workdir, "xtbscan.log")
    for label in ("hess-step", "sccacc"):
        if logs[label] == logs["default"]:
            fail("--lscan-%s did not change the scan" % label)


#: The progress line the numerical Hessian provider prints before it starts
#: its displaced single points (`src/landscape/types/hessian.f90`).
HESS_PROGRESS_LINE = re.compile(
    r"^ (computing|rebuilding) the numerical Hessian: (\d+) displaced "
    r"single points \(6N, step +5[.]00E-03 bohr\)$")
#: One row of the live per-step table of a minimisation
#: (`landscape_print_step`, `src/landscape/types/results.f90`): iteration
#: number, the 12-decimal energy and the 4-decimal dE column.  The
#: `convergence history` table printed after the loop has a `es15.6` gnorm in
#: that column instead, so the two tables cannot be confused.
HESS_STEP_LINE = re.compile(r"^ +(\d{1,6}) +-?\d+[.]\d{12} +-?\d+[.]\d{4} ")


def hess_events(output):
    """The Hessian progress lines and the macro steps of a minimisation.

    Returns the events between the settings block and the `converged:`
    summary in output order, as `("build", "computing"|"rebuilding")` and
    `("step", n)`.  The cut at the summary drops the post-loop
    `convergence history` table, which repeats the step numbers.
    """
    stop = re.search(r"(?m)^ converged:", output)
    text = output[:stop.start()] if stop else output
    events = []
    for line in text.splitlines():
        build = HESS_PROGRESS_LINE.match(line)
        if build is not None:
            events.append(("build", build.group(1)))
            continue
        step = HESS_STEP_LINE.match(line)
        if step is not None:
            events.append(("step", int(step.group(1))))
    return events


def hess_builds(output):
    """The progress lines of a run, in order, as `(verb, n_eval)` pairs."""
    return [(match.group(1), int(match.group(2)))
            for match in (HESS_PROGRESS_LINE.match(line)
                          for line in output.splitlines())
            if match is not None]


def check_hess_progress(exe, inp, rundir, stretch):
    """The Hessian build explains its pause before the pause happens.

    Mirrors `xtb/landscape-hess-progress`.  The numerical provider is 6N
    displaced single points at `printlevel = -1`, so the run used to go quiet
    for the whole build (`docs/_simplify-recon/201`); the line names the
    source and the exact number of evaluations the provider will make, it is
    printed **before** the first step row, and a periodic rebuild announces
    itself the same way.  A `model` Hessian costs no evaluation and therefore
    prints nothing at all.
    """
    common = ["--gfn", "2", "--norestart", "--lopt", "--lopt-engine", "rfo"]
    status, output, _ = knobs_run(exe, stretch, common, rundir,
                                  "hessprogress-default")
    if status != 0:
        fail("the default RFO minimisation failed with status %d" % status,
             output)
    builds = hess_builds(output)
    if [verb for verb, _ in builds] != ["computing", "rebuilding"]:
        fail("the run did not announce its Hessian builds (%r)" % (builds,),
             output)
    if any(n_eval != 18 for _, n_eval in builds):
        fail("the progress line does not name the 6N = 18 evaluations the "
             "provider makes (%r)" % (builds,), output)
    events = hess_events(output)
    if events[:1] != [("build", "computing")]:
        fail("the progress line is not the first event after the settings "
             "block (%r)" % (events[:3],), output)
    try:
        rebuild = events.index(("build", "rebuilding"))
    except ValueError:
        fail("the rebuild announced no progress line (%r)" % (events,), output)
    before = [event for event in events[:rebuild] if event[0] == "step"]
    after = [event for event in events[rebuild + 1:] if event[0] == "step"]
    if not before or not after:
        fail("the run has no step rows around the rebuild (%r)" % (events,),
             output)
    if before[-1] != ("step", 5) or after[0] != ("step", 6):
        fail("the rebuild line does not sit between step 5 and step 6 (%r)"
             % (events,), output)

    status, output, _ = knobs_run(
        exe, stretch, common + ["--lopt-hessian", "model"], rundir,
        "hessprogress-model")
    if status != 0:
        fail("the model Hessian run failed with status %d" % status, output)
    if hess_builds(output):
        fail("a model Hessian (0 evaluations) printed a progress line",
             output)


def check_opt_recalc(exe, inp, rundir, stretch):
    """`--lopt-recalc` sets the Hessian rebuild interval of the RFO step.

    Mirrors `xtb/landscape-opt-recalc`.  The field was fixed at 5 with no
    flag reaching it (`src/landscape/types/settings.f90`); the check pins the
    semantics: 0 disables the rebuild (one Hessian, the step-6 pause is
    gone), 2 rebuilds four times, the explicit 5 is the default (byte
    identical geometry), and a value that is not a non-negative integer is
    reported (fatal with `--strict`).
    """
    common = ["--gfn", "2", "--norestart", "--lopt", "--lopt-engine", "rfo"]

    def interval(value):
        label = "optrecalc-%s" % value if value else "optrecalc-default"
        extra = ["--lopt-recalc", value] if value else []
        status, output, workdir = knobs_run(exe, stretch, common + extra,
                                           rundir, label)
        if status != 0:
            fail("--lopt-recalc %s failed with status %d"
                 % (value or "(default)", status), output)
        if not re.search(r"Hessian rebuild interval: +%s steps"
                         % (value if value else "5"), output):
            fail("--lopt-recalc %s did not reach the minimisation"
                 % (value or "(default)"), output)
        return output, workdir

    # Measured on the stretched water: interval 0 converges in 12 steps after
    # one Hessian, the default interval 5 in 10 steps after two (the rebuild
    # at step 6), interval 2 in 5 steps after three (rebuilds at steps 2 and
    # 4).  The shorter trajectory of the tighter interval is the reason the
    # counts are not simply `steps / interval + 1`.
    for value, expected_builds, expected_hessians in (("0", 1, 1),
                                                      ("5", 2, 2),
                                                      ("2", 3, 3)):
        output, workdir = interval(value)
        verbs = [verb for verb, _ in hess_builds(output)]
        if len(verbs) != expected_builds:
            fail("--lopt-recalc %s announced %d builds, expected %d (%r)"
                 % (value, len(verbs), expected_builds, verbs), output)
        reported = re.search(r"landscape engine: Hessians +(\d+)\s*$",
                             output, re.M)
        if reported is None or int(reported.group(1)) != expected_hessians:
            fail("--lopt-recalc %s reported %r Hessian builds, expected %d"
                 % (value, reported.group(1) if reported else None,
                    expected_hessians), output)
        if value == "0":
            if not re.search(r"(?m)^ steps: +12\s*$", output):
                fail("--lopt-recalc 0 did not change the trajectory", output)
            if "rebuilding the numerical Hessian" in output:
                fail("--lopt-recalc 0 still rebuilt the Hessian", output)
        elif value == "5":
            if not re.search(r"(?m)^ steps: +10\s*$", output):
                fail("the default interval (5) does not stop after 10 steps",
                     output)
            if "rebuilding the numerical Hessian" not in output:
                fail("--lopt-recalc 5 did not rebuild the Hessian", output)

    # An explicit 5 is the default: the same number of builds and the same
    # geometry.  The two runs are compared by file, because stdout carries
    # machine dependent timing lines.
    _, default_dir = interval("")
    _, five_dir = interval("5")
    if knobs_file(default_dir, "xtbopt.xyz") != knobs_file(five_dir,
                                                           "xtbopt.xyz"):
        fail("--lopt-recalc 5 is not the default interval")

    status, output, _ = knobs_run(
        exe, stretch, common + ["--lopt-recalc", "x"], rundir,
        "optrecalc-invalid")
    if status == 0:
        fail("--lopt-recalc x was accepted", output)
    if "The --lopt-recalc flag requires an integer" not in output:
        fail("the unparsable value was not reported", output)

    negative = common + ["--lopt-recalc", "-3"]
    status, output, _ = knobs_run(exe, stretch, negative, rundir,
                                  "optrecalc-negative")
    if status != 0:
        fail("--lopt-recalc -3 aborted without --strict (status %d)" % status,
             output)
    if "[WARNING]" not in output:
        fail("--lopt-recalc -3 was not reported in lax mode", output)
    if "--lopt-recalc needs 0 (no rebuild)" not in output:
        fail("the lax report does not name the flag", output)
    status, output, _ = knobs_run(exe, stretch, negative + ["--strict"],
                                  rundir, "optrecalc-negative-strict")
    if status == 0:
        fail("--lopt-recalc -3 did not abort with --strict", output)
    if "--lopt-recalc needs 0 (no rebuild)" not in output:
        fail("the --strict report does not name the flag", output)


def check_hess_step_alias(exe, inp, rundir, second):
    """The canonical `--l*-hess-step` names and their older spellings.

    Mirrors `xtb/landscape-hess-step-alias` (CLI redesign sections 6.2, 6.3,
    7.2 and 9.2).  Three field names were inconsistent (`--lts-step`,
    `--lirc-step-hess` and `--lvib-step`), so each pair gets a canonical
    spelling and the old name stays as an alias.  The check pins both halves
    of the alias contract:

      * the two spellings given the **same** value are equivalent - the output
        files are byte identical (and each of them differs from the default
        run, so the comparison cannot pass by both flags being ignored);
      * the two spellings given **different** values abort with a non-zero
        status and a message that names both flags;
      * a repeated flag keeps the last value.

    `second` is the bent HCN fixture: the `--lts` and `--lirc` runs need a
    transition state as their input and cannot use the water minimum.
    """
    common = ["--gfn", "2", "--norestart"]

    def compare(label, fixture, args, files, old, new):
        status, output, _ = knobs_run(exe, fixture, args, rundir, label + "-default")
        if status != 0:
            fail("%s: the default run failed with status %d" % (label, status), output)
        status, output, olddir = knobs_run(
            exe, fixture, args + [old, "0.02"], rundir, label + "-old")
        if status != 0:
            fail("%s: %s failed with status %d" % (label, old, status), output)
        status, output, newdir = knobs_run(
            exe, fixture, args + [new, "0.02"], rundir, label + "-new")
        if status != 0:
            fail("%s: %s failed with status %d" % (label, new, status), output)
        defaultdir = os.path.join(rundir, label + "-default")
        for name in files:
            if knobs_file(olddir, name) != knobs_file(newdir, name):
                fail("%s: %s and %s %s are not equivalent"
                     % (label, old, new, name))
            if knobs_file(newdir, name) == knobs_file(defaultdir, name):
                fail("%s: %s did not change %s" % (label, new, name))

    compare("alias-ts", second, common + ["--lts", "--lts-method", "prfo"],
            ("xtbts.log", "xtbts.xyz"), "--lts-step", "--lts-hess-step")
    compare("alias-lirc", second,
            common + ["--lirc", "--lirc-method", "gs", "--lirc-max-steps", "1"],
            ("irc.log", "irc_path.xyz"), "--lirc-step-hess", "--lirc-hess-step")
    compare("alias-lvib", inp, common + ["--lvib"], ("hessian",),
            "--lvib-step", "--lvib-hess-step")

    # the same value through both spellings is not a conflict, also when the
    # two spellings are written differently (`1.0e-2` vs `0.01`)
    status, output, _ = knobs_run(
        exe, inp, common + ["--lvib", "--lvib-step", "0.01",
                            "--lvib-hess-step", "1.0e-2"],
        rundir, "alias-lvib-same")
    if status != 0:
        fail("the same step through both spellings was rejected", output)
    if not re.search(r"Hessian step: +1[.]000000E-02 bohr", output):
        fail("the shared value did not reach the task", output)
    # a repeated flag keeps the last value
    status, output, _ = knobs_run(
        exe, inp, common + ["--lvib", "--lvib-step", "0.01",
                            "--lvib-step", "0.03"],
        rundir, "alias-lvib-repeat")
    if status != 0:
        fail("the repeated flag was rejected", output)
    if not re.search(r"Hessian step: +3[.]000000E-02 bohr", output):
        fail("the repeated flag did not keep the last value", output)

    # a different value through the two spellings is a loud error naming both
    for runflag, old, new in (("--lts", "--lts-step", "--lts-hess-step"),
                              ("--lirc", "--lirc-step-hess", "--lirc-hess-step"),
                              ("--lvib", "--lvib-step", "--lvib-hess-step")):
        status, output = run(
            exe, inp, common + [runflag, old, "0.02", new, "0.03"], rundir)
        if status == 0:
            fail("the conflicting spellings of the Hessian step exited 0", output)
        if "conflicting Hessian finite difference steps" not in output:
            fail("the conflict was not reported for %s" % runflag, output)
        if new not in output or old not in output:
            fail("the conflict message does not name both flags", output)


def check_restart_knobs(exe, inp, rundir):
    """The two new restart flags reach the engine's restart policy.

    Mirrors `xtb/landscape-restart-knobs` (CLI redesign sections 6.3 and 9.2).
    `--lopt-restart` and `--lvib-restart` close the two tasks that had no
    policy flag (`--lscan-restart`, `--lts-restart` and `--lirc-restart`
    existed).  Each mode is checked twice: the engine report has to read back
    the requested policy (it prints the control object the task handed over)
    and the result has to change - the trajectory for the minimisation, the
    `hessian` file for the frequency task.  Measured magnitude of the change:
    the `off`/`pin` geometries differ in the 12th decimal, which is a real,
    reproducible difference (a repeated `pin` run is byte identical).
    """
    common = ["--gfn", "2", "--norestart"]
    for task, files in (("--lopt", "xtb.trj"), ("--lvib", "hessian")):
        results = {}
        for mode in ("off", "pin", "warm"):
            label = "restart-%s-%s" % (task[2:], mode)
            status, output, workdir = knobs_run(
                exe, inp, common + [task, "%s-restart" % task, mode], rundir,
                label)
            if status != 0:
                fail("%s %s-restart %s failed with status %d"
                     % (task, task, mode, status), output)
            if not re.search(r"restart policy: +%s\s*$" % mode, output, re.M):
                fail("%s did not echo the %s policy" % (task, mode), output)
            if not re.search(r"landscape engine: restart policy +%s\s*$" % mode,
                             output, re.M):
                fail("the engine did not read back the %s policy" % mode, output)
            results[mode] = knobs_file(workdir, files)
        for mode in ("pin", "warm"):
            if results[mode] == results["off"]:
                fail("%s-restart %s did not change %s" % (task, mode, files))


def check_hessian_step_used(exe, inp, rundir, hcn, nh3, nh3product):
    """The effectivity echo of the finite difference parameters.

    Mirrors `xtb/landscape-hess-step-used` (CLI redesign sections 7.3 and 9.2,
    fix report `35-stage2b-fix.md` section T1).  `--l*-hess-step` and
    `--l*-sccacc` are knobs of the numerical Hessian, and whether the selected
    algorithm and Hessian source consume them differs per task.  The tasks now
    echo the effective value and mark the parameters as unused otherwise, and
    this check pins both halves for all four tasks that own the knobs:

      * a combination that consumes the step echoes the value **from the
        settings object** (the flag value, not the default) and carries no
        `(not used:` line;
      * a combination that does not consume it echoes the same value and the
        one line reason.

    The predicate of `--lopt`, `--lscan` and `--lirc` is the measured one of
    the fact table (only an RFO minimiser and, for `--lirc`, the initial
    Hessian of every method builds a numerical Hessian).  `--lts` and `--lvib`
    do not decide anything before the run: their notice is derived from the
    Hessian builds the engine actually performed, so it is the same fact as
    the `landscape engine: Hessians` counter and cannot contradict it
    (`37-stage2b-fix2.md` sections T1 and T2).  That is why the `--lts` cases
    below assert the observed run, not a method name: `--lts-method autoneb
    --lts-no-verify` runs to its natural end and **does** build a Hessian (the
    candidate modes of `path_ts_list.xyz`), so the old method-name guess
    reported it as unused next to a non-zero counter.

    Two cases exist to catch a predicate that is written in the wrong
    direction: `--lts --lts-method prfo --lts-no-verify` still consumes the
    step through the search Hessian and must stay silent, and a non-numerical
    Hessian source marks the step unused even for RFO.

    The step and the accuracy flag are **one family** (the numerical provider
    consumes both, `model`/`file` ignore both), so the `--lts`/`--lvib` notice
    appears when the user named **either** of them and this run read neither.
    The `*-sccacc-only` cases give the accuracy flag alone - no
    `--l*-hess-step` - which the gate that only looked at `hess_step_source`
    left silent (`38-stage2b-fix2-review.md` LOW-1).
    """
    common = ["--gfn", "2", "--norestart"]
    step_text = "2.000000E-02"
    acc_text = "1.000000E+01"

    def expect(label, fixture, args, hint, used=True, require=None,
               step=step_text):
        status, output, _ = knobs_run(exe, fixture, args, rundir, label)
        if require is not None and not re.search(require, output):
            fail("%s: the run did not reach the expected state (%s)"
                 % (label, require), output)
        if not re.search(r"Hessian step: +%s bohr" % re.escape(step),
                         output):
            fail("%s: the run does not echo the effective step %s"
                 % (label, step), output)
        if not re.search(r"Hessian SCC accuracy: +%s\s*$" % re.escape(acc_text),
                         output, re.M):
            fail("%s: the run does not echo the effective SCC accuracy %s"
                 % (label, acc_text), output)
        # the echo alone does not prove the claim: a used combination has to
        # build at least one Hessian through the engine (the counter counts
        # every provider, so it is only asserted on the used side)
        if used:
            counter = re.search(r"landscape engine: Hessians +([0-9]+)", output)
            if counter is None or int(counter.group(1)) < 1:
                fail("%s: the run claims to consume the step but built no "
                     "Hessian" % label, output)
        match = re.findall(r"^\s*\(not used: (.*)\)\s*$", output, re.M)
        if hint is None:
            # a run can carry `(not used: ...)` notices that have nothing to do
            # with the finite difference knobs (`--lvib` names the energy
            # tolerance it never reads, a rigid scan names the missing inner
            # relaxation), so only a notice about the Hessian counts here - the
            # assertion is about *this* step, not about the absence of every
            # unused notice
            hessian_notices = [reason for reason in match
                               if "Hessian" in reason or "finite" in reason]
            if hessian_notices:
                fail("%s: the step is consumed but the run reports it as "
                     "unused (%r)" % (label, hessian_notices), output)
        else:
            if not match:
                fail("%s: the run does not report the unused step" % label,
                     output)
            # a run can carry more than one `(not used: ...)` notice (a rigid
            # scan names the missing inner relaxation as well), so the hint has
            # to appear in one of them instead of in the first one
            if not any(hint in reason for reason in match):
                fail("%s: the unused reasons are %r, expected one with %r"
                     % (label, match, hint), output)
        return status, output

    # --- opt ---------------------------------------------------------------
    opt = common + ["--lopt", "--lopt-max-steps", "2"]
    expect("used-opt-rfo", inp,
           opt + ["--lopt-engine", "rfo", "--lopt-hessian", "numerical",
                  "--lopt-hess-step", "0.02", "--lopt-sccacc", "10.0"],
           None)
    expect("unused-opt-lbfgs", inp,
           opt + ["--lopt-hess-step", "0.02", "--lopt-sccacc", "10.0"],
           "lbfgs builds no numerical Hessian", used=False)
    expect("unused-opt-model", inp,
           opt + ["--lopt-engine", "rfo", "--lopt-hessian", "model",
                  "--lopt-hess-step", "0.02", "--lopt-sccacc", "10.0"],
           "the Hessian source 'model' has no finite differences", used=False)

    # --- scan --------------------------------------------------------------
    scan = common + ["--lscan", "--lscan-coord", "angle", "--lscan-atoms",
                     "2,1,3", "--lscan-step", "3", "--lscan-nstep", "1",
                     "--lscan-max-steps", "2"]
    expect("used-scan-rfo", inp,
           scan + ["--lscan-method", "rfo", "--lscan-hessian", "numerical",
                   "--lscan-hess-step", "0.02", "--lscan-sccacc", "10.0"],
           None)
    expect("unused-scan-lbfgs", inp,
           scan + ["--lscan-method", "lbfgs", "--lscan-hess-step", "0.02",
                   "--lscan-sccacc", "10.0"],
           "lbfgs builds no numerical Hessian", used=False)
    expect("unused-scan-rigid", inp,
           scan + ["--lscan-mode", "rigid", "--lscan-method", "rfo",
                   "--lscan-hess-step", "0.02", "--lscan-sccacc", "10.0"],
           "a rigid scan builds no Hessian", used=False)

    # --- ts ----------------------------------------------------------------
    # `--lts` does not decide before the run: whether a build happens depends
    # on the method, on the confirmation and on the part of the path the run
    # reaches, so the notice comes from the engine counters (T1).
    ts = common + ["--lts", "--lts-max-steps", "2"]
    expect("used-ts-prfo", hcn,
           ts + ["--lts-method", "prfo", "--lts-hess-step", "0.02",
                 "--lts-sccacc", "10.0"], None)
    # the search Hessian of P-RFO keeps consuming the step without the
    # confirmation, so this must not be marked unused
    expect("used-ts-prfo-no-verify", hcn,
           ts + ["--lts-method", "prfo", "--lts-no-verify",
                 "--lts-hess-step", "0.02", "--lts-sccacc", "10.0"], None)
    # AutoNEB builds no Hessian during the band relaxation, but it does build
    # one for the candidate modes of `path_ts_list.xyz`, so a run that reaches
    # its natural end reports a non-zero counter and must carry no notice.
    expect("used-ts-autoneb-no-verify", nh3,
           common + ["--lts", "--lts-method", "autoneb", "--lts-images", "5",
                     "--lts-product", nh3product, "--lts-no-verify",
                     "--lts-hess-step", "0.02", "--lts-sccacc", "10.0"],
           None, require=r"autoneb: termination autoneb completed")
    expect("unused-ts-dimer-no-verify", hcn,
           ts + ["--lts-method", "dimer", "--lts-no-verify",
                 "--lts-hess-step", "0.02", "--lts-sccacc", "10.0"],
           "no numerical Hessian was built by this run", used=False)
    expect("unused-ts-neb-no-verify", nh3,
           ts + ["--lts-method", "neb", "--lts-images", "5", "--lts-product",
                 nh3product, "--lts-no-verify", "--lts-hess-step", "0.02",
                 "--lts-sccacc", "10.0"],
           "no numerical Hessian was built by this run", used=False)
    expect("unused-ts-model", hcn,
           ts + ["--lts-method", "prfo", "--lts-hessian", "model",
                 "--lts-hess-step", "0.02", "--lts-sccacc", "10.0"],
           "the Hessian source 'model' has no finite differences", used=False)

    # --- the accuracy knob of the same family (LOW-1) ----------------------
    # `--l*-sccacc` alone (no `--l*-hess-step`), so the echo shows the builtin
    # step.  A run that reads neither knob has to carry the notice, and one
    # that reads the accuracy flag (it builds a numerical Hessian) must not.
    default_step = "5.000000E-03"
    expect("unused-ts-sccacc-only", hcn,
           ts + ["--lts-method", "dimer", "--lts-no-verify",
                 "--lts-sccacc", "10.0"],
           "no numerical Hessian was built by this run", used=False,
           step=default_step)
    expect("used-ts-sccacc-only", hcn,
           ts + ["--lts-method", "prfo", "--lts-sccacc", "10.0"], None,
           step=default_step)

    # --- irc ---------------------------------------------------------------
    irc = common + ["--lirc", "--lirc-method", "gs", "--lirc-max-steps", "1"]
    expect("used-irc-gs", hcn,
           irc + ["--lirc-hess-step", "0.02", "--lirc-sccacc", "10.0"], None)
    # measured: the model Hessian of this fixture has no imaginary mode, so the
    # run aborts at the transition state confirmation after the banner; the
    # step is not consumed either way
    status, output = expect("unused-irc-model", hcn,
                            irc + ["--lirc-hessian", "model",
                                   "--lirc-hess-step", "0.02",
                                   "--lirc-sccacc", "10.0"],
                            "the Hessian source 'model' has no finite differences",
                            used=False)
    if status == 0:
        fail("the model Hessian of the IRC fixture did not abort the run",
             output)

    # --- lvib --------------------------------------------------------------
    # The same-source notice for the frequency task (T2): the source decides,
    # not the method, and the flag has to be visible in the file it claims to
    # change.
    lvib = common + ["--lvib", "--lvib-sccacc", "10.0"]
    expect("unused-lvib-model", inp,
           lvib + ["--lvib-hessian", "model", "--lvib-hess-step", "0.02"],
           "the Hessian source 'model' has no finite differences", used=False)
    expect("used-lvib-numerical", inp,
           lvib + ["--lvib-hessian", "numerical", "--lvib-hess-step", "0.02"],
           None)
    # the same pair for the accuracy flag alone (LOW-1)
    expect("unused-lvib-sccacc-only", inp,
           common + ["--lvib", "--lvib-hessian", "model", "--lvib-sccacc",
                     "10.0"],
           "the Hessian source 'model' has no finite differences", used=False,
           step=default_step)
    expect("used-lvib-sccacc-only", inp,
           common + ["--lvib", "--lvib-hessian", "numerical",
                     "--lvib-sccacc", "10.0"], None, step=default_step)
    for source, changes in (("numerical", True), ("model", False)):
        files = {}
        for step in ("0.005", "0.02"):
            _, _, workdir = knobs_run(
                exe, inp, common + ["--lvib", "--lvib-sccacc", "10.0",
                                    "--lvib-hessian", source,
                                    "--lvib-hess-step", step],
                rundir, "lvib-hess-step-%s-%s" % (source, step))
            files[step] = knobs_file(workdir, "hessian")
        if changes and files["0.005"] == files["0.02"]:
            fail("--lvib-hess-step did not change the %s hessian" % source)
        if not changes and files["0.005"] != files["0.02"]:
            fail("--lvib-hess-step changed the %s hessian although the source "
                 "has no finite differences" % source)


def check_lvib_restart_warm(exe, inp, rundir):
    """`--lvib-restart warm` is not distinguishable from `pin`.

    Mirrors `xtb/landscape-lvib-restart-warm` (fix report
    `35-stage2b-fix.md` section T2).  The finite difference sequence of the
    frequency task resets the wavefunction to the state it entered the Hessian
    with before every displacement, so `warm` (keep the incoming state) and
    `pin` (reuse the recorded state) hand the provider the same start: the
    `hessian` files are byte identical, while `off` records the cold start
    wavefunction and differs.  The equivalence is a documented property, not a
    defect, and this check pins it so that a later change of the restart policy
    cannot make the documentation wrong silently.
    """
    results = {}
    for mode in ("off", "pin", "warm"):
        status, output, workdir = knobs_run(
            exe, inp, ["--gfn", "2", "--norestart", "--lvib",
                       "--lvib-restart", mode],
            rundir, "lvib-restart-" + mode)
        if status != 0:
            fail("--lvib-restart %s failed with status %d" % (mode, status),
                 output)
        if not re.search(r"landscape engine: restart policy +%s\s*$" % mode,
                         output, re.M):
            fail("the engine did not read back the %s policy" % mode, output)
        results[mode] = knobs_file(workdir, "hessian")
    if results["pin"] == results["off"]:
        fail("--lvib-restart pin did not change the hessian file")
    if results["warm"] != results["pin"]:
        fail("--lvib-restart warm differs from pin: the documented "
             "equivalence of `--help` no longer holds")


# ---------------------------------------------------------------------------
# CLI redesign stage 3: `--opt`/`--ohess` are executed by the landscape suite
# ---------------------------------------------------------------------------
#
# `--opt`/`--ohess` keep their run types and their identity, but the run they
# start is the suite task `opt` - plus the frequency analysis for `--ohess` -
# and the native ANCopt family is reachable through `--method
# ancopt|lancopt|fire|pbc_lbfgs` only.  Every check below is written against a
# token the *other* path prints, so "the wrong path ran" and "both paths ran"
# are both failures:
#
#   suite minimisation      `landscape: <algorithm> minimisation`
#   suite frequencies       `landscape: harmonic frequencies`
#   native ANCopt           `ANC micro-cycles`
#   native L-ANCopt         `L-ANC optimizer`
#   native FIRE             `initial time step`
#   native pbc L-BFGS       `current lattice:`
#   native Hessian          `Numerical Hessian`
#   native frequency table  `Frequency Printout`
#
# The design is `docs/_simplify-recon/11-cli-redesign.md` (sections 3.8, 3.11,
# 3.11.1); the implementation is `docs/_simplify-recon/42-stage3a-impl.md`.

TOTAL_RE = re.compile(r"TOTAL ENTHALPY\s+([-0-9.eE+]+) Eh")
GIBBS_RE = re.compile(r"TOTAL FREE ENERGY\s+([-0-9.eE+]+) Eh")
TOTAL_ENERGY_RE = re.compile(r"TOTAL ENERGY\s+([-0-9.eE+]+) Eh")
LVIB_H_RE = re.compile(r"total enthalpy H\(T\):\s+([-0-9.eE+]+) Eh")
LVIB_G_RE = re.compile(r"total free energy G\(T\):\s+([-0-9.eE+]+) Eh")


def one(pattern, output, what):
    """Return the single match of `pattern` in `output` (or fail)."""
    found = pattern.findall(output)
    if len(found) != 1:
        fail("expected exactly one %s, found %d" % (what, len(found)), output)
    return float(found[0])


def check_opt_alias_once(exe, inp, rundir):
    """`--opt` runs the landscape suite, exactly once.

    Mirrors `xtb/landscape-opt-alias-once` (CLI redesign sections 3.11 and
    4.7.1, R24/M-3).  Two independent signals prove "once": the suite prints
    its algorithm line while the native engine prints `ANC micro-cycles`
    (`src/optimizer.f90:392`), and the suite writes `landscape_opt.json`
    although `--opt` was a run type the `--ljson` switch warned about.  The
    document is compared with the one of `--lopt`, because "the redeployed
    `--opt` is the same task" is the point of the alias.
    """
    rundir = os.path.abspath(rundir)
    shutil.rmtree(rundir, ignore_errors=True)
    optdir = os.path.join(rundir, "opt")
    loptdir = os.path.join(rundir, "lopt")
    status, output = run(
        exe, inp, ["--gfn", "2", "--opt", "--ljson", "--norestart"], optdir
    )
    if status != 0:
        fail("--opt exited %d" % status, output)
    if "landscape: lbfgs minimisation" not in output:
        fail("--opt did not run the landscape suite", output)
    if "ANC micro-cycles" in output:
        fail("--opt ran the native ANCopt as well", output)
    if "[NOTICE]" not in output:
        fail("the redeployed --opt did not announce the semantic change", output)
    if "writes nothing without a landscape run type" in output:
        fail("--ljson reported the redeployed --opt as a run type that knows "
             "nothing about the switch", output)
    status2, output2 = run(
        exe, inp, ["--gfn", "2", "--lopt", "--ljson", "--norestart"], loptdir
    )
    if status2 != 0:
        fail("--lopt exited %d" % status2, output2)
    document = read_json(optdir, "landscape_opt.json")
    reference = read_json(loptdir, "landscape_opt.json")
    if document.get("task") != "opt":
        fail("the --opt document is not an opt result", output)
    if document.get("counters", {}).get("n_steps") != \
            reference.get("counters", {}).get("n_steps"):
        fail("--opt and --lopt report a different number of steps: %s vs %s"
             % (document.get("counters"), reference.get("counters")),
             output + output2)
    print("opt alias: --opt ran the suite (%s steps, same as --lopt)"
          % document.get("counters", {}).get("n_steps"))


def check_opt_native_exit(exe, inp, rundir):
    """The four native engines behind `--method` keep `--opt` native.

    Mirrors `xtb/landscape-opt-native-exit` (CLI redesign sections 3.8 and
    4.7.1).  The three guards, the suppression of the native optimisation block
    and the `--ljson` criterion all read the same flag, so a run that names a
    native engine has to keep the native path **and** must not write a
    landscape document.
    """
    rundir = os.path.abspath(rundir)
    for name, token in (("ancopt", "ANC micro-cycles"),
                        ("lancopt", "L-ANC optimizer"),
                        ("fire", "initial time step"),
                        ("pbc_lbfgs", "current lattice:")):
        workdir = os.path.join(rundir, "native-" + name)
        shutil.rmtree(workdir, ignore_errors=True)
        status, output = run(
            exe, inp,
            ["--gfn", "2", "--opt", "--method", name, "--ljson", "--norestart"],
            workdir,
        )
        if status != 0:
            fail("--opt --method %s exited %d" % (name, status), output)
        if token not in output:
            fail("--opt --method %s did not reach the native engine" % name,
                 output)
        if "landscape: " in output:
            fail("--opt --method %s ran the landscape suite" % name, output)
        if "[NOTICE]" in output:
            fail("--opt --method %s announced a redeployment" % name, output)
        # Report 117: the native run types write the unified document too, so
        # the old "writes nothing without a landscape run type" contract is
        # retired.  The native engine must still leave the suite untouched
        # (asserted above), but it now writes the same core document every
        # task writes.
        document = read_json(workdir, "landscape_opt.json")
        if document.get("task") != "opt":
            fail("--opt --method %s wrote a non-opt document" % name, output)
        if document.get("status_name") not in ("ok", "incomplete", "error"):
            fail("--opt --method %s wrote a status_name %r"
                 % (name, document.get("status_name")), output)
        for key in ("format", "n_atoms", "energy", "geometry", "provenance"):
            if key not in document:
                fail("--opt --method %s document lacks the core key %s"
                     % (name, key), output)
        print("opt native exit: --method %-10s -> %s" % (name, token))


def check_ohess_suite_once(exe, inp, rundir):
    """`--ohess` is the suite minimisation followed by the suite frequencies.

    Mirrors `xtb/landscape-ohess-suite-once` (CLI redesign section 3.11.1,
    T1/T2/T3).  The native Hessian block and the native frequency printout are
    suppressed through the same flag as the native optimisation, and their
    tokens (`Numerical Hessian`, `Frequency Printout`, `ANC micro-cycles`) are
    the falsifiable part; the container they fill is covered by
    `xtb/landscape-ohess-fres-bridge`.  Each stage of the chain must run
    exactly once, which is what the two milestone counts and the single
    `landscape engine: Hessians           1` line of the frequency report show.
    """
    rundir = os.path.abspath(rundir)
    shutil.rmtree(rundir, ignore_errors=True)
    status, output = run(exe, inp, ["--gfn", "2", "--ohess", "--norestart"], rundir)
    if status != 0:
        fail("--ohess exited %d" % status, output)
    for token, expected in (("Numerical Hessian", 0), ("Frequency Printout", 0),
                            ("ANC micro-cycles", 0),
                            ("analytical hessian", 0),
                            ("landscape: lbfgs minimisation", 1),
                            ("landscape: harmonic frequencies", 1),
                            ("landscape engine: Hessians           1", 1)):
        count = output.count(token)
        if count != expected:
            fail("--ohess printed %r %d times instead of %d"
                 % (token, count, expected), output)
    for name in ("hessian", "vibspectrum", "g98.out", "xtbopt.xyz"):
        if not os.path.exists(os.path.join(rundir, name)):
            fail("--ohess wrote no %s" % name, output)
    print("ohess suite: one minimisation + one frequency analysis, "
          "no native Hessian/frequency block")


def check_ohess_level_axes(exe, inp, rundir):
    """`--ohess <LEVEL>` drives both stages of the chain (stage 3c close-out).

    Mirrors `xtb/landscape-ohess-level-axes`.  `--ohess` is a single run made
    of two stages - the suite minimisation and the frequency analysis of the
    minimised geometry - and the level used to reach the minimisation only
    (`lset`): the frequency half kept its own defaults, so one run printed
    `landscape engine: SCC accuracy 0.200000` for the minimisation and
    `1.000000` (level `normal`) for the frequency analysis.  Both halves now
    read the same level, and the axes are observable in the banner of each
    half, which is the falsifier of the fix:

    * `--ohess tight` reports 0.2 and 1e-6 twice (once per stage),
    * `--ohess tight --acc 0.5` reports 0.5 twice: the explicit `--acc` still
      wins over the level in **both** halves,
    * the plain `--ohess` reports the level `normal` of both halves (1.0, 5e-6),
      so a fix that unconditionally overwrote the frequency defaults (instead
      of following the same level) fails here.
    """
    rundir = os.path.abspath(rundir)
    shutil.rmtree(rundir, ignore_errors=True)

    acc = r"landscape engine: SCC accuracy\s+(\S+)"
    etol = r"energy tolerance:\s+(\S+)"

    def axes(label, args):
        workdir = os.path.join(rundir, label)
        status, output = run(exe, inp, ["--gfn", "2", "--ohess"] + args +
                             ["--norestart"], workdir)
        if status != 0:
            fail("--ohess %s exited %d" % (" ".join(args), status), output)
        return re.findall(acc, output), re.findall(etol, output)

    tight, tight_etol = axes("tight", ["tight"])
    if tight != ["0.200000", "0.200000"]:
        fail("--ohess tight reported the SCC accuracies %r instead of two "
             "times 0.200000: the level does not reach both stages" % (tight,))
    if tight_etol != ["1.000000E-06", "1.000000E-06"]:
        fail("--ohess tight reported the energy tolerances %r instead of two "
             "times 1.000000E-06" % (tight_etol,))

    explicit, _ = axes("tight-acc", ["tight", "--acc", "0.5"])
    if explicit != ["0.500000", "0.500000"]:
        fail("--ohess tight --acc 0.5 reported the SCC accuracies %r instead "
             "of two times 0.500000: the explicit --acc does not win in both "
             "stages" % (explicit,))

    plain, plain_etol = axes("plain", [])
    if plain != ["1.000000", "1.000000"]:
        fail("--ohess without a level reported the SCC accuracies %r instead "
             "of two times 1.000000 (level normal)" % (plain,))
    if plain_etol != ["5.000000E-06", "5.000000E-06"]:
        fail("--ohess without a level reported the energy tolerances %r "
             "instead of two times 5.000000E-06" % (plain_etol,))

    print("ohess level: both stages of the chain report the same SCC accuracy "
          "and energy tolerance (level and explicit --acc)")


def check_ohess_native_hess(exe, inp, rundir):
    """`--ohess --method ancopt` keeps the native Hessian and frequency path.

    Mirrors `xtb/landscape-ohess-native-hess`, the reverse direction of the
    suppression: a criterion that is always false (for instance one that only
    looks at the suite) would silently drop the native Hessian here, and the
    missing `Numerical Hessian` line and the missing files are the failure.
    """
    rundir = os.path.abspath(rundir)
    shutil.rmtree(rundir, ignore_errors=True)
    status, output = run(
        exe, inp, ["--gfn", "2", "--ohess", "--method", "ancopt", "--norestart"],
        rundir,
    )
    if status != 0:
        fail("--ohess --method ancopt exited %d" % status, output)
    for token, expected in (("Numerical Hessian", 1), ("Frequency Printout", 1),
                            ("analytical hessian", 1),
                            ("landscape: harmonic frequencies", 0),
                            ("landscape: lbfgs minimisation", 0)):
        count = output.count(token)
        if count != expected:
            fail("--ohess --method ancopt printed %r %d times instead of %d"
                 % (token, count, expected), output)
    for name in ("hessian", "vibspectrum", "g98.out"):
        if not os.path.exists(os.path.join(rundir, name)):
            fail("--ohess --method ancopt wrote no %s" % name, output)
    print("ohess native: native Hessian and frequency path kept")


def check_ohess_fres_bridge(exe, inp, rundir):
    """The `--ohess` frequency bridge reproduces the native console contract.

    Mirrors `xtb/landscape-ohess-fres-bridge` (CLI redesign section 3.11.1,
    T3).  `anyhess` stays the criterion of the three readers of the frequency
    container, so the redeployed `--ohess` has to *fill* it: without the bridge
    the property block loses `TOTAL ENTHALPY`/`TOTAL FREE ENERGY`, and a half
    fix (container printed, not filled) prints `TOTAL ENERGY 0.000000000000`.
    The totals are compared with `--lvib` on the minimised geometry, so a
    bridge that fills the wrong field fails as well.
    """
    rundir = os.path.abspath(rundir)
    shutil.rmtree(rundir, ignore_errors=True)
    ohessdir = os.path.join(rundir, "ohess")
    vibdir = os.path.join(rundir, "vib")
    status, output = run(exe, inp, ["--gfn", "2", "--ohess", "--norestart"], ohessdir)
    if status != 0:
        fail("--ohess exited %d" % status, output)
    energy = one(TOTAL_ENERGY_RE, output, "TOTAL ENERGY line")
    enthalpy = one(TOTAL_RE, output, "TOTAL ENTHALPY line")
    gibbs = one(GIBBS_RE, output, "TOTAL FREE ENERGY line")
    if abs(energy) < 1.0e-8:
        fail("--ohess reports TOTAL ENERGY 0: the frequency container was "
             "printed but not filled", output)
    status2, output2 = run(
        exe, os.path.join(ohessdir, "xtbopt.xyz"),
        ["--gfn", "2", "--lvib", "--norestart"], vibdir,
    )
    if status2 != 0:
        fail("--lvib on the minimised geometry exited %d" % status2, output2)
    reference_h = one(LVIB_H_RE, output2, "total enthalpy H(T) line")
    reference_g = one(LVIB_G_RE, output2, "total free energy G(T) line")
    if abs(enthalpy - reference_h) > 1.0e-6:
        fail("--ohess TOTAL ENTHALPY %.12f differs from --lvib %.12f"
             % (enthalpy, reference_h), output + output2)
    if abs(gibbs - reference_g) > 1.0e-6:
        fail("--ohess TOTAL FREE ENERGY %.12f differs from --lvib %.12f"
             % (gibbs, reference_g), output + output2)
    print("ohess bridge: TOTAL ENERGY %.12f, H(T) %.12f (--lvib %.12f)"
          % (energy, enthalpy, reference_h))


def check_opt_guard(exe, inp, rundir, case):
    """The three guards of the redeployed `--opt` (CLI redesign section 3.8).

    Mirrors `xtb/landscape-opt-scanblock`, `xtb/landscape-opt-pbc` and
    `xtb/landscape-opt-xcontrol`.  Each case asserts both directions: the
    redeployed `--opt` refuses the input with a message that names the native
    exit, and the same input through `--opt --method ancopt` still works -
    a guard that is too wide would fail the second half, a guard that is
    missing would fail the first.  The scan case additionally pins the absence
    of `xtbscan.log` (the block must not run at all) and the xcontrol case the
    full scope of C7 (`--lopt` refuses the same block).
    """
    rundir = os.path.abspath(rundir)
    shutil.rmtree(rundir, ignore_errors=True)
    os.makedirs(rundir, exist_ok=True)

    def invoke(label, input_file, args):
        workdir = os.path.join(rundir, label)
        os.makedirs(workdir)
        return run(exe, input_file, args, workdir) + (workdir,)

    if case == "scanblock":
        control = os.path.join(rundir, "scan.inp")
        with open(control, "w") as handle:
            handle.write("$constrain\n  distance: 1, 2, auto\n"
                         "$scan\n  1: 0.9, 1.1, 3\n$end\n")
        status, output, workdir = invoke(
            "guard", inp, ["--gfn", "2", "--opt", "--input", control, "--norestart"]
        )
        if status == 0:
            fail("--opt accepted the xcontrol $scan block", output)
        # the message of the guard itself; the `[NOTICE]` line of the same run
        # also contains `--method ancopt`, so the bare substring is not enough
        if "does not execute the xcontrol $scan block: use --opt --method ancopt" \
                not in output:
            fail("the refusal is not the $scan guard of the redeployed --opt",
                 output)
        if os.path.exists(os.path.join(workdir, "xtbscan.log")):
            fail("--opt wrote xtbscan.log although it refused the block", output)
        status, output, workdir = invoke(
            "baseline", inp,
            ["--gfn", "2", "--opt", "--method", "ancopt", "--input", control,
             "--norestart"],
        )
        if status != 0:
            fail("the native exit refused the $scan block (exit %d)" % status,
                 output)
        if "RELAXED SCAN" not in output:
            fail("the native exit did not execute the $scan block", output)
        if not os.path.exists(os.path.join(workdir, "xtbscan.log")):
            fail("the native exit wrote no xtbscan.log", output)
        print("opt guard scanblock: refused with %r, native exit runs it"
              % "--method ancopt")
    elif case == "pbc":
        status, output, _ = invoke(
            "guard", inp, ["--gfnff", "--opt", "--norestart"]
        )
        if status == 0:
            fail("--opt accepted a periodic system", output)
        # the guard, not the generic refusal of the driver and not the
        # `[NOTICE]` line (both of which also mention the native exit)
        if "runs the landscape suite, which does not support periodic " \
                "boundary conditions: use --opt --method ancopt" not in output:
            fail("the refusal is not the periodic guard of the redeployed "
                 "--opt", output)
        status, output, _ = invoke(
            "baseline", inp,
            ["--gfnff", "--opt", "--method", "ancopt", "--norestart"],
        )
        if status != 0:
            fail("the native exit refused the periodic system (exit %d)" % status,
                 output)
        if "current lattice:" not in output:
            fail("the native periodic optimiser did not run", output)
        print("opt guard pbc: refused with %r, native exit runs pbc_lbfgs"
              % "--method ancopt")
    elif case == "xcontrol":
        control = os.path.join(rundir, "opt.inp")
        with open(control, "w") as handle:
            handle.write("$opt\n   optlevel=tight\n")
        status, output, _ = invoke(
            "guard", inp, ["--gfn", "2", "--opt", "--input", control, "--norestart"]
        )
        if status == 0:
            fail("--opt accepted the xcontrol $opt block", output)
        if "the xcontrol $opt block does nothing for --opt" not in output:
            fail("the refusal is not the $opt guard of the redeployed --opt",
                 output)
        status, output, _ = invoke(
            "baseline", inp,
            ["--gfn", "2", "--opt", "--method", "ancopt", "--input", control,
             "--norestart"],
        )
        if status != 0:
            fail("the native exit refused the $opt block (exit %d)" % status,
                 output)
        # C7 covers every suite run type, not only the redeployed `--opt`
        status, output, _ = invoke(
            "lopt", inp, ["--gfn", "2", "--lopt", "--input", control, "--norestart"]
        )
        if status == 0:
            fail("--lopt accepted the xcontrol $opt block", output)
        if "the xcontrol $opt block does nothing for --lopt" not in output:
            fail("the --lopt refusal is not the $opt guard", output)
        print("opt guard xcontrol: refused for --opt and --lopt, "
              "native exit honours it")
    else:
        fail("unknown opt guard case %r" % case)


def check_opt_xtbopt_once(exe, inp, rundir):
    """`--opt` writes `xtbopt.xyz` once and labels its timer correctly.

    Mirrors `xtb/landscape-opt-xtbopt-once` (CLI redesign read points 13/14).
    The native optimisers have a second geometry writer of their own next to
    the one of the suite, and the timing block names the timer that ran; for
    the redeployed `--opt` both have to follow the suite.
    """
    rundir = os.path.abspath(rundir)
    shutil.rmtree(rundir, ignore_errors=True)
    plain = os.path.join(rundir, "plain")
    status, output = run(exe, inp, ["--gfn", "2", "--opt", "--norestart"], plain)
    if status != 0:
        fail("--opt exited %d" % status, output)
    if output.count("optimized geometry written to:") != 0:
        fail("--opt printed the native geometry message as well", output)
    if output.count("wrote xtbopt.xyz") != 1:
        fail("--opt wrote xtbopt.xyz %d times instead of 1"
             % output.count("wrote xtbopt.xyz"), output)
    if "landscape optimizer" not in output:
        fail("the timing block does not name the landscape optimiser", output)
    if "ANC optimizer" in output:
        fail("the timing block still names the native optimiser", output)
    namespaced = os.path.join(rundir, "namespaced")
    status, output = run(
        exe, inp, ["--gfn", "2", "--opt", "--namespace", "c1", "--norestart"],
        namespaced,
    )
    if status != 0:
        fail("--opt --namespace c1 exited %d" % status, output)
    if not os.path.exists(os.path.join(namespaced, "c1.xtbopt.xyz")):
        fail("--opt --namespace c1 wrote no c1.xtbopt.xyz", output)
    if os.path.exists(os.path.join(namespaced, "xtbopt.xyz")):
        fail("--opt --namespace c1 also wrote the unprefixed xtbopt.xyz", output)
    print("opt xtbopt: one writer, landscape timer, namespace honoured")


def count_xyz_frames(path):
    """Number of complete XMOL frames in `path`."""
    with open(path) as fh:
        lines = fh.read().splitlines()
    n = 0
    i = 0
    while i < len(lines):
        try:
            nat = int(lines[i].split()[0])
        except (ValueError, IndexError):
            break
        if nat <= 0 or i + nat + 2 > len(lines):
            break
        n += 1
        i += nat + 2
    return n


def check_file_contract_opt(exe, inp, rundir):
    """A converged `--opt` writes `.xtboptok` and the `xtbopt.log` trajectory.

    Mirrors `xtb/landscape-file-contract-opt` (CLI redesign section 3.9).  The
    upstream optimiser creates `.xtboptok` only on success
    (`src/geoopt_driver.f90:186`) and opens `xtbopt.log` as a multi frame XMOL
    file (`:121`, frames at `src/optimizer.f90:801`); the redeployed `--opt`
    has to follow both.
    """
    rundir = os.path.abspath(rundir)
    shutil.rmtree(rundir, ignore_errors=True)
    run_dir = os.path.join(rundir, "run")
    status, output = run(exe, inp, ["--gfn", "2", "--opt", "--norestart"], run_dir)
    if status != 0:
        fail("--opt exited %d" % status, output)
    ok_path = os.path.join(run_dir, ".xtboptok")
    if not os.path.exists(ok_path):
        fail("a converged --opt wrote no .xtboptok", output)
    if os.path.getsize(ok_path) != 0:
        fail(".xtboptok is not empty", output)
    if os.path.exists(os.path.join(run_dir, "NOT_CONVERGED")):
        fail("a converged --opt also wrote NOT_CONVERGED", output)
    log_path = os.path.join(run_dir, "xtbopt.log")
    if not os.path.exists(log_path):
        fail("--opt wrote no xtbopt.log", output)
    frames = count_xyz_frames(log_path)
    if frames < 2:
        fail("xtbopt.log is not a multi frame trajectory (%d frames)" % frames,
             output)
    print("file contract opt: .xtboptok + xtbopt.log (%d frames)" % frames)


def check_file_contract_opt_fail(exe, inp, rundir):
    """A failed `--opt` writes `NOT_CONVERGED` and `xtblast.xyz`, not `.xtboptok`.

    Mirrors `xtb/landscape-file-contract-opt-fail` (section 3.9).  Upstream
    (`src/geoopt_driver.f90:183`) creates `NOT_CONVERGED` instead of
    `.xtboptok` on failure, and the native `--opt` writes the last geometry to
    `xtblast.<ext>` with the message `last geometry written to: <name>`
    (`src/prog/main.F90:1024-1034`; the block moved when stage 3b added the
    sentinel handling above it).
    """
    rundir = os.path.abspath(rundir)
    shutil.rmtree(rundir, ignore_errors=True)
    run_dir = os.path.join(rundir, "run")
    status, output = run(
        exe, inp, ["--gfn", "2", "--opt", "--cycles", "2", "--norestart"],
        run_dir,
    )
    if status == 0:
        fail("a two step --opt converged unexpectedly", output)
    if os.path.exists(os.path.join(run_dir, ".xtboptok")):
        fail("a failed --opt wrote .xtboptok", output)
    if not os.path.exists(os.path.join(run_dir, "NOT_CONVERGED")):
        fail("a failed --opt wrote no NOT_CONVERGED", output)
    last_path = os.path.join(run_dir, "xtblast.xyz")
    if not os.path.exists(last_path):
        fail("a failed --opt wrote no xtblast.xyz", output)
    if count_xyz_frames(last_path) != 1:
        fail("xtblast.xyz is not a single geometry", output)
    if "last geometry written to:" not in output:
        fail("the last geometry message is missing", output)
    print("file contract opt-fail: NOT_CONVERGED + xtblast.xyz, no .xtboptok")


def check_file_contract_vib(exe, inp, rundir):
    """`--lvib` writes the legacy `compliance.dat` (CLI redesign section 3.9).

    Mirrors `xtb/landscape-file-contract-vib`.  The native `--hess` runs
    `compliance_driver` after the numerical Hessian (`src/hessian.F90:331`)
    under the condition `mol%n > 1 .and. freezeset%n == 0 .and.
    .not. set%periodic`; the suite has to produce the same file.  It is
    written under the bare name (`src/compliance.f90:287`).  The frequency
    half of the redeployed `--ohess` is the second writer of the same file.
    """
    rundir = os.path.abspath(rundir)
    shutil.rmtree(rundir, ignore_errors=True)
    run_dir = os.path.join(rundir, "run")
    status, output = run(exe, inp, ["--gfn", "2", "--lvib"], run_dir)
    if status != 0:
        fail("--lvib exited %d" % status, output)
    path = os.path.join(run_dir, "compliance.dat")
    if not os.path.exists(path):
        fail("--lvib wrote no compliance.dat", output)
    with open(path) as fh:
        text = fh.read()
    if "compliance.dat" not in text:
        fail("compliance.dat does not name itself", text)
    # The `--ohess` chain runs the same frequency half, so it writes the same
    # file (`landscape_ohess_chain` -> `tasks/vib.f90`); the assertion is here
    # and not in the failure-path check, because it is a property of the
    # successful chain.
    ohess_dir = os.path.join(rundir, "ohess")
    status, output = run(
        exe, inp, ["--gfn", "2", "--ohess", "--norestart"], ohess_dir)
    if status != 0:
        fail("--ohess exited %d" % status, output)
    if not os.path.exists(os.path.join(ohess_dir, "hessian")):
        fail("--ohess wrote no hessian", output)
    if not os.path.exists(os.path.join(ohess_dir, "compliance.dat")):
        fail("--ohess wrote no compliance.dat in its frequency half", output)
    print("file contract vib: compliance.dat written by --lvib and by --ohess")


def check_vib_compliance_skip(exe, inp, rundir):
    """The converse of the `compliance.dat` contract.

    Mirrors `xtb/landscape-vib-compliance-skip`.  `compliance_driver` is called
    only for a nonperiodic system of more than one atom whose `$fix` mask is
    empty, so a single atom and a run with `$fix atoms:` must not write the
    file.  Without this case the condition is only covered from the side where
    it is true, and deleting either clause would stay green.
    """
    rundir = os.path.abspath(rundir)
    shutil.rmtree(rundir, ignore_errors=True)
    # 1. a single atom (`mol%n > 1` is false)
    atom = os.path.join(rundir, "landscape_he.xyz")
    os.makedirs(rundir, exist_ok=True)
    with open(atom, "w") as handle:
        handle.write("1\nhelium\nHe  0.0000000  0.0000000  0.0000000\n")
    single_dir = os.path.join(rundir, "single")
    status, output = run(exe, atom, ["--gfn", "2", "--lvib"], single_dir)
    if status != 0:
        fail("--lvib on a single atom exited %d" % status, output)
    if os.path.exists(os.path.join(single_dir, "compliance.dat")):
        fail("a single atom wrote compliance.dat", output)
    # 2. `$fix atoms:` (the `$fix` mask is not empty)
    xcontrol = os.path.join(rundir, "landscape_vib_fix.inp")
    with open(xcontrol, "w") as handle:
        handle.write("$fix\n   atoms: 1\n$end\n")
    fixed_dir = os.path.join(rundir, "fixed")
    status, output = run(
        exe, inp, ["--gfn", "2", "--lvib", "--input", xcontrol], fixed_dir)
    if status != 0:
        fail("--lvib with $fix atoms: exited %d" % status, output)
    if os.path.exists(os.path.join(fixed_dir, "compliance.dat")):
        fail("$fix atoms: still wrote compliance.dat", output)
    if not os.path.exists(os.path.join(fixed_dir, "hessian")):
        fail("the frozen run produced no hessian at all", output)
    print("vib compliance skip: no compliance.dat for a single atom / $fix atoms:")


def check_file_contract_namespace(exe, inp, rundir):
    """The legacy names carry the namespace, except `compliance.dat`.

    Mirrors `xtb/landscape-file-contract-namespace` (section 3.9).  Upstream
    routes `.xtboptok`/`xtbopt.log`/`xtblast.xyz` through `open_file`/
    `touch_file`, which adds the prefix (`.c1.xtboptok`, `c1.xtbopt.log`),
    while `compliance.dat` uses a bare Fortran `open` and stays unprefixed.
    """
    rundir = os.path.abspath(rundir)
    shutil.rmtree(rundir, ignore_errors=True)
    opt_dir = os.path.join(rundir, "opt")
    status, output = run(
        exe, inp, ["--gfn", "2", "--opt", "--namespace", "c1", "--norestart"],
        opt_dir,
    )
    if status != 0:
        fail("--opt --namespace c1 exited %d" % status, output)
    if not os.path.exists(os.path.join(opt_dir, ".c1.xtboptok")):
        fail("--opt --namespace c1 wrote no .c1.xtboptok", output)
    if not os.path.exists(os.path.join(opt_dir, "c1.xtbopt.log")):
        fail("--opt --namespace c1 wrote no c1.xtbopt.log", output)
    if os.path.exists(os.path.join(opt_dir, ".xtboptok")):
        fail("--opt --namespace c1 also wrote the bare .xtboptok", output)
    if os.path.exists(os.path.join(opt_dir, "xtbopt.log")):
        fail("--opt --namespace c1 also wrote the bare xtbopt.log", output)
    vib_dir = os.path.join(rundir, "vib")
    status, output = run(
        exe, inp, ["--gfn", "2", "--lvib", "--namespace", "c1"], vib_dir)
    if status != 0:
        fail("--lvib --namespace c1 exited %d" % status, output)
    if not os.path.exists(os.path.join(vib_dir, "compliance.dat")):
        fail("--lvib --namespace c1 wrote no bare compliance.dat", output)
    if os.path.exists(os.path.join(vib_dir, "c1.compliance.dat")):
        fail("compliance.dat must not carry the namespace prefix", output)
    print("file contract namespace: legacy names prefixed, compliance.dat bare")


def check_ohess_fail(exe, inp, rundir):
    """A failed `--ohess` aborts with the opt error, without the cascade.

    Mirrors `xtb/landscape-ohess-fail` (CLI redesign section 3.9, T2).  The
    frequency half runs only on a converged geometry, so a failed minimisation
    reports its single error, writes `NOT_CONVERGED`, and produces neither the
    frequency files nor the misleading `Could not generate lattice points` /
    `Setup of Coulomb evaluator failed` cascade of the old chain.
    """
    rundir = os.path.abspath(rundir)
    shutil.rmtree(rundir, ignore_errors=True)
    fail_dir = os.path.join(rundir, "fail")
    status, output = run(
        exe, inp, ["--gfn", "2", "--ohess", "--cycles", "2", "--norestart"],
        fail_dir,
    )
    if status == 0:
        fail("a two step --ohess converged unexpectedly", output)
    for needle in ("Could not generate lattice points",
                   "Setup of Coulomb evaluator failed"):
        if needle in output:
            fail("the misleading cascade is back: %s" % needle, output)
    if os.path.exists(os.path.join(fail_dir, ".xtboptok")):
        fail("a failed --ohess wrote .xtboptok", output)
    if not os.path.exists(os.path.join(fail_dir, "NOT_CONVERGED")):
        fail("a failed --ohess wrote no NOT_CONVERGED", output)
    for name in ("hessian", "vibspectrum"):
        if os.path.exists(os.path.join(fail_dir, name)):
            fail("a failed --ohess wrote %s" % name, output)
    # The contrast: a converged --ohess runs both halves.
    ok_dir = os.path.join(rundir, "ok")
    status, output = run(exe, inp, ["--gfn", "2", "--ohess", "--norestart"], ok_dir)
    if status != 0:
        fail("a converged --ohess exited %d" % status, output)
    if not os.path.exists(os.path.join(ok_dir, ".xtboptok")):
        fail("a converged --ohess wrote no .xtboptok", output)
    if not os.path.exists(os.path.join(ok_dir, "hessian")):
        fail("a converged --ohess wrote no hessian", output)
    print("ohess fail: single opt error, NOT_CONVERGED, no cascade; "
          "converged run writes the frequency files")


def check_native_runtypes_unaffected(exe, inp, rundir):
    """`--omd`/`--screen`/`--metaopt` still run their native ANCopt first.

    Mirrors `xtb/landscape-native-runtypes-unaffected` (review HIGH-1 of
    `52-stage3b-review.md`, a regression of stage 3a).  The suppression of the
    native optimisation block reads `landscape_is_suite_run`, not
    `landscape_native_engine /= 0`: these three run types are never redeployed
    and cannot take `--method` at all (`main.F90` refuses it with "only
    available with --opt or --ohess"), so an `/= 0` test makes the block
    unreachable for them - `--omd` then runs its dynamics on the **input**
    geometry.  The geometry golden values below are the ones the pre-stage-3a
    binary `/tmp/xtb-land/build/xtb` produced for this fixture
    (`--gfn 2 --omd|--screen|--metaopt`), measured in the fix round; the
    optimisation is deterministic (only the dynamics draws the run's random
    seed), so the comparison is exact.

    `--screen` needs a `screen` file in the working directory and aborts
    without one **after** the optimisation; that abort is upstream behaviour
    and identical before and after the fix (`norm`-ed stdout of the two
    binaries differs only in backtrace addresses), so the expected status is
    pinned to it instead of being skipped.
    """
    golden_energy = -5.070544363296
    cases = [("--omd", "--omd", 0), ("--screen", "--screen", 1),
             ("--metaopt", "--metaopt", 0)]
    rundir = os.path.abspath(rundir)
    shutil.rmtree(rundir, ignore_errors=True)
    for label, flag, want_status in cases:
        case_dir = os.path.join(rundir, label.strip("-"))
        status, output = run(
            exe, inp, ["--gfn", "2", flag, "--norestart"], case_dir)
        if status != want_status:
            fail("%s exited %d instead of the native %d"
                 % (label, status, want_status), output)
        if "ANC micro-cycles" not in output:
            fail("%s did not run the native ANCopt (no 'ANC micro-cycles')"
                 % label, output)
        if "optimized geometry written to:" not in output:
            fail("%s did not report an optimized geometry" % label, output)
        path = os.path.join(case_dir, "xtbopt.xyz")
        if not os.path.exists(path):
            fail("%s wrote no xtbopt.xyz" % label, output)
        with open(path) as handle:
            lines = handle.read().splitlines()
        # frame comment of `writeMolecule`: ' energy: <E> gnorm: <g> xtb: ...'
        energy = float(lines[1].split("energy:")[1].split()[0])
        if abs(energy - golden_energy) > 1.0e-9:
            fail("%s ended at %s instead of the native %s"
                 % (label, energy, golden_energy), "\n".join(lines))
        # the input O-H distance is 1.3000 angstrom; the native minimisation
        # takes it to 0.96, so a geometry that was never optimised is caught
        oxygen = [float(v) for v in lines[2].split()[1:4]]
        longest = 0.0
        for line in lines[3:5]:
            hydrogen = [float(v) for v in line.split()[1:4]]
            longest = max(longest, sum((a - b) ** 2
                                       for a, b in zip(oxygen, hydrogen)) ** 0.5)
        if not 0.9 < longest < 1.05:
            fail("%s left the O-H distance at %.4f angstrom" % (label, longest),
                 "\n".join(lines))
    print("native runtypes unaffected: --omd/--screen/--metaopt optimise first, "
          "energy %s, O-H %.4f" % (golden_energy, longest))


def check_opt_level_bridge(exe, inp, rundir):
    """`--opt <LEVEL>` drives the suite convergence level (redesign 3.5).

    Mirrors `xtb/landscape-opt-level-bridge` (stage 3a fix round, finding B and
    the destructive change F3).  Before the fix the upstream level keyword
    reached the native ANCopt only, which `--opt` left in stage 3, so `--opt`,
    `--opt tight` and `--opt crude` were the same run.  The bridge makes
    `--opt tight` the same run as `--lopt --lopt-level tight` (byte-identical
    geometry and trajectory), `--opt crude` a *different* one (it maps to
    `superloose`), an unknown keyword an error that names the flag, and the
    level flag of the suite unusable with the `--opt` run type.
    """
    rundir = os.path.abspath(rundir)
    shutil.rmtree(rundir, ignore_errors=True)
    opt = os.path.join(rundir, "opt-tight")
    lopt = os.path.join(rundir, "lopt-tight")
    status, output = run(exe, inp,
                         ["--gfn", "2", "--opt", "tight", "--norestart"], opt)
    if status != 0:
        fail("--opt tight exited %d" % status, output)
    if "landscape convergence level: tight" not in output:
        fail("--opt tight did not report the bridged level", output)
    status2, output2 = run(
        exe, inp,
        ["--gfn", "2", "--lopt", "--lopt-level", "tight", "--norestart"],
        lopt,
    )
    if status2 != 0:
        fail("--lopt --lopt-level tight exited %d" % status2, output2)
    for name in ("xtbopt.xyz", "xtb.trj"):
        with open(os.path.join(opt, name), "rb") as handle:
            bridged = handle.read()
        with open(os.path.join(lopt, name), "rb") as handle:
            reference = handle.read()
        if bridged != reference:
            fail("--opt tight and --lopt --lopt-level tight wrote different %s"
                 % name, output + output2)

    def steps(label, args):
        workdir = os.path.join(rundir, label)
        run(exe, inp, ["--gfn", "2"] + args + ["--ljson", "--norestart"], workdir)
        return read_json(workdir, "landscape_opt.json")["counters"]["n_steps"]

    crude = steps("opt-crude", ["--opt", "crude"])
    tight = steps("opt-tight-json", ["--opt", "tight"])
    if crude == tight:
        fail("--opt crude and --opt tight produced the same run (%d steps)"
             % crude)

    status, output = run(exe, inp,
                         ["--gfn", "2", "--opt", "banana", "--norestart"],
                         os.path.join(rundir, "banana"))
    if status == 0:
        fail("--opt banana was accepted silently", output)
    if "Unknown --opt level 'banana'" not in output:
        fail("the unknown --opt level was not named", output)

    status, output = run(
        exe, inp,
        ["--gfn", "2", "--opt", "tight", "--lopt-level", "loose",
         "--norestart"],
        os.path.join(rundir, "conflict"),
    )
    if status == 0:
        fail("--opt tight with --lopt-level loose was accepted", output)
    print("opt level bridge: --opt tight == --lopt --lopt-level tight, "
          "crude %d != tight %d steps" % (crude, tight))


def check_cycles_mapping(exe, inp, rundir):
    """`--cycles N` drives the suite step budget of the run type (3.10b).

    Mirrors `xtb/landscape-cycles-map`, `-lvib` and `-conflict` (stage 3a fix
    round, finding C).  `--cycles` used to set the native
    `set%optset%maxoptcycle` that no landscape run type reads, so it was a dead
    knob.  The water fixture needs five steps by default, so a budget of four
    binds and the document reports exactly four steps - the same number as
    `--lopt-max-steps 4`, and fewer than the default.  A budget the
    minimisation cannot meet exits non-zero, which is expected; the document
    is written in either case, so only `counters.n_steps` is asserted.
    """
    rundir = os.path.abspath(rundir)
    shutil.rmtree(rundir, ignore_errors=True)

    def steps(label, args):
        workdir = os.path.join(rundir, label)
        run(exe, inp, ["--gfn", "2"] + args + ["--ljson", "--norestart"],
            workdir)
        return read_json(workdir, "landscape_opt.json")["counters"]["n_steps"]

    default = steps("default", ["--lopt"])
    cycles = steps("cycles", ["--lopt", "--cycles", "4"])
    flag = steps("flag", ["--lopt", "--lopt-max-steps", "4"])
    opt = steps("opt", ["--opt", "--cycles", "4"])
    if cycles != flag:
        fail("--lopt --cycles 4 reports %d steps, --lopt --lopt-max-steps 4 %d"
             % (cycles, flag))
    if cycles >= default:
        fail("--lopt --cycles 4 does not bind (default %d, cycles %d)"
             % (default, cycles))
    if opt != cycles:
        fail("--opt --cycles 4 reports %d steps, --lopt --cycles 4 %d"
             % (opt, cycles))
    same = steps("same", ["--lopt", "--cycles", "4", "--lopt-max-steps", "4"])
    if same != cycles:
        fail("--cycles and --lopt-max-steps of the same value disagree")

    status, output = run(
        exe, inp,
        ["--gfn", "2", "--lopt", "--cycles", "4", "--lopt-max-steps", "5",
         "--norestart"],
        os.path.join(rundir, "conflict"),
    )
    if status == 0:
        fail("--cycles 4 with --lopt-max-steps 5 was accepted", output)
    if "--lopt-max-steps" not in output:
        fail("the step budget conflict does not name --lopt-max-steps", output)

    status, output = run(
        exe, inp, ["--gfn", "2", "--lvib", "--cycles", "4", "--norestart"],
        os.path.join(rundir, "lvib"),
    )
    if status == 0:
        fail("--lvib --cycles 4 was accepted", output)
    if "--cycles has no meaning for --lvib" not in output:
        fail("--lvib --cycles did not report the missing step budget", output)
    print("cycles map: --cycles 4 == --lopt-max-steps 4 (%d steps, default %d), "
          "lvib and conflict refused" % (cycles, default))


def blas_backend(output, requested):
    """Check one `--blas-threads` read-back line and return the backend.

    Three renderings are documented (`docs/_simplify-recon/78-blas-threads-
    flag.md` section 1.1) and any of them is correct for the build under test:

    * `blas threads: N (OpenBLAS, read back M)` - the OpenBLAS form,
    * `blas threads: N (MKL, read back M)` - the MKL form with its getter,
    * `blas threads: N (MKL, not readable)` - a setter without a getter,
    * the warning `no OpenBLAS or MKL symbol found` - a reference BLAS, or a
      statically linked BLAS whose symbol never enters the dynamic scope.

    Whenever a count is printed back it must equal the requested one: the test
    environment pins `OMP_NUM_THREADS=1`, so the pool already is one and only a
    setter that ran can report the request.  The function returns the backend
    label of the read-back form (`'OpenBLAS'`, `'MKL'`, or `None` for the
    warning) so the caller can require the same backend in a second run.
    """
    match = re.search(
        r"blas threads: (\d+) \((OpenBLAS|MKL), read back (\d+)\)", output
    )
    if match:
        if int(match.group(1)) != requested:
            fail("the read-back line reports a request of %s, not %d"
                 % (match.group(1), requested), output)
        if int(match.group(3)) != requested:
            fail("--blas-threads %d did not read back %d (backend reported "
                 "%s): with OMP_NUM_THREADS=1 the pool already is one, so "
                 "only a setter that ran can report the requested value"
                 % (requested, requested, match.group(3)), output)
        return match.group(2)
    if re.search(r"blas threads: %d \(MKL, not readable\)" % requested, output):
        return "MKL"
    if re.search(r"no OpenBLAS or MKL symbol found", output):
        return None
    fail("--blas-threads %d printed none of the documented lines (an OpenBLAS "
         "or MKL read-back, `(MKL, not readable)`, or the no-symbol warning)"
         % requested, output)


def check_blas_threads(exe, inp, rundir):
    """`--blas-threads N` sets the BLAS pool and the setting is observable.

    Mirrors `xtb/landscape-blas-threads`.  `--blas-threads` goes through
    `dlopen(NULL, RTLD_NOW)` + `dlsym`, so the printed line is the proof that
    the runtime setter of the linked backend really ran:

    * a plain `--gfn 2 --sp` prints no `blas threads:` line at all - the
      default leaves the backend untouched (no banner, no behaviour change),

    * `--sp --blas-threads 1` stays a single point: the flag is a machine
      setting and must not go through `select_runtyp`, where it would be
      reported as a run type conflict,

    * the requested value is read back from the backend, so
      `blas threads: 1 (OpenBLAS, read back 1)` has to appear,

    * the falsifier: the test environment pins `OMP_NUM_THREADS=1`, so the
      process pool starts at one and a request of three can only read back
      three if the setter was called.  An implementation that resolves the
      symbol but skips the call reads back one and fails here,

    * an invalid value names the flag instead of claiming that nothing was
      provided (the branch had no coverage before).

    `blas_backend` accepts the three documented fallbacks so that a build
    linked against MKL or the reference BLAS does not turn the registered
    check red; the build under test links OpenBLAS and therefore exercises the
    strict read-back form.
    """
    rundir = os.path.abspath(rundir)
    shutil.rmtree(rundir, ignore_errors=True)

    status, plain = run(exe, inp, ["--gfn", "2", "--sp", "--norestart"], rundir)
    if status != 0:
        fail("xtb --gfn 2 --sp terminated with status %d" % status, plain)
    if re.search(r"blas threads:", plain):
        fail("--gfn 2 --sp without --blas-threads printed a 'blas threads:' "
             "line; the default must not touch the backend", plain)

    status, one = run(
        exe, inp,
        ["--gfn", "2", "--sp", "--blas-threads", "1", "--norestart"],
        rundir,
    )
    if status != 0:
        fail("--sp --blas-threads 1 terminated with status %d (run type "
             "conflict?)" % status, one)
    backend = blas_backend(one, 1)

    status, three = run(
        exe, inp,
        ["--gfn", "2", "--sp", "--blas-threads", "3", "--norestart"],
        rundir,
    )
    if status != 0:
        fail("--blas-threads 3 terminated with status %d" % status, three)
    if blas_backend(three, 3) != backend:
        fail("the symbol lookup changed between two runs of the same binary "
             "(%r then %r)" % (backend, blas_backend(three, 3)), three)

    check_error(exe, inp, rundir,
                "--blas-threads needs at least one thread",
                ["--gfn", "2", "--sp", "--blas-threads", "0", "--norestart"])
    check_error(exe, inp, rundir,
                "--blas-threads needs at least one thread",
                ["--gfn", "2", "--sp", "--blas-threads", "-1", "--norestart"])
    for bad in ("abc", "1.5"):
        check_error(exe, inp, rundir,
                    "--blas-threads needs an integer >= 1",
                    ["--gfn", "2", "--sp", "--blas-threads", bad,
                     "--norestart"])

    print("blas threads: silent by default, no conflict with --sp, the backend "
          "reports the requested value back, and invalid values name the flag")


def check_blas_threads_timing(exe, inp, rundir):
    """`--blas-threads 1` removes the OpenMP x BLAS oversubscription.

    Mirrors `xtb/landscape-blas-threads-timing`.  The test runs with
    `OMP_NUM_THREADS=8` and - crucially - an explicitly pinned
    `OPENBLAS_NUM_THREADS=8`: that is the variable the `--help` block tells
    users to export, and an inherited value of one would make both sides of the
    ratio serial and fail the assertion for an environment reason.  The test is
    `RUN_SERIAL` in CMake because the ratio is a wall clock measurement that a
    parallel `ctest -j` neighbour would distort.  With the pool at eight the
    GFN2 single point of a system just over the gemm threading threshold
    (anthracene, 24 atoms, 66 basis functions, 66^3 > 262144) is oversubscribed
    8x8.  Measured: 14.15-14.95x faster with `--blas-threads 1`
    (`docs/_simplify-recon/80-blas-threads-validation.md`); the assertion asks
    for 3x, a margin wide enough for a loaded machine.  A machine with fewer
    than four CPUs has almost no pool to oversubscribe, so the check is skipped
    there.
    """
    if (os.cpu_count() or 1) < 4:
        print("blas threads timing: skipped, the machine has fewer than four "
              "CPUs")
        return
    rundir = os.path.abspath(rundir)
    shutil.rmtree(rundir, ignore_errors=True)

    def measure(args, repeats=2):
        best = None
        for _ in range(repeats):
            start = time.monotonic()
            status, output = run(
                exe, inp,
                ["--gfn", "2", "--sp", "--norestart"] + args,
                rundir,
            )
            elapsed = time.monotonic() - start
            if status != 0:
                fail("the timing run %s failed with status %d"
                     % (" ".join(args) or "(default)", status), output)
            best = elapsed if best is None else min(best, elapsed)
        return best

    measure(["--blas-threads", "1"])  # warm up the page cache
    default = measure([])
    limited = measure(["--blas-threads", "1"])
    if limited * 3.0 > default:
        fail("--blas-threads 1 took %.3f s against the default %.3f s: less "
             "than the required 3x speedup (OMP_NUM_THREADS=8, anthracene)"
             % (limited, default))
    print("blas threads timing: %.3f s default vs %.3f s with "
          "--blas-threads 1 (%.1fx)" % (default, limited, default / limited))


def check_level_axes(exe, inp, rundir):
    """A convergence level drives `e_tol`, `acc` and the step budget (3.5/R5).

    Mirrors `xtb/landscape-level-axes` (stage 3c).  Three axes are asserted on
    the stdout banner, which is the only place they are observable:

    * `--lopt --lopt-level loose` reports `SCC accuracy 2.000000`; before the
      guard of `acc_was_set` the global `--acc` default 1.0 overwrote the
      level-derived value and the axis was a dead knob, so this line is the
      falsifier of T2.
    * `--lopt --lopt-level loose --acc 1.0` reports 1.0: an explicit `--acc`
      still wins over the level.
    * `energy tolerance` follows the level (5e-5 loose, 1e-6 tight, 5e-4
      superloose); the line did not exist before stage 3c, so the axis had no
      observable at all.
    * `maximum number of steps` is `clamp(n_atoms * multiplier, 200, 10000)`
      for water (3 atoms, so the lower clamp 200 binds at every level), and an
      explicit `--lopt-max-steps 8` or `--cycles 4` keeps its value.

    The scaling with the atom count (the multiplier itself) is asserted in the
    unit test `landscape-settings / level-step-budget`, which can build a 2000
    atom settings object without paying for a 2000 atom SCF.
    """
    rundir = os.path.abspath(rundir)
    shutil.rmtree(rundir, ignore_errors=True)

    def banner(label, args):
        workdir = os.path.join(rundir, label)
        status, output = run(exe, inp,
                             ["--gfn", "2"] + args + ["--norestart"], workdir)
        return output

    def value(output, pattern, label):
        match = re.search(pattern, output)
        if match is None:
            fail("the banner does not report %s" % label, output)
        return match.group(1)

    acc = r"landscape engine: SCC accuracy\s+(\S+)"
    etol = r"energy tolerance:\s+(\S+)"
    budget = r"maximum number of steps:\s+(\S+)"

    default = banner("default", ["--lopt"])
    if value(default, etol, "energy tolerance:") != "5.000000E-06":
        fail("the default level does not report e_tol 5e-6", default)
    if value(default, acc, "the SCC accuracy of the engine") != "1.000000":
        fail("the default level does not report acc 1.0", default)
    if value(default, budget, "maximum number of steps:") != "200":
        fail("the default level does not report clamp(3n, 200, 10000) = 200",
             default)

    loose = banner("loose", ["--lopt", "--lopt-level", "loose"])
    if value(loose, acc, "the SCC accuracy of the engine") != "2.000000":
        fail("--lopt-level loose does not report acc 2.0 (dead knob, T2)",
             loose)
    if value(loose, etol, "energy tolerance:") != "5.000000E-05":
        fail("--lopt-level loose does not report e_tol 5e-5 (T3)", loose)

    explicit = banner("loose-acc",
                      ["--lopt", "--lopt-level", "loose", "--acc", "1.0"])
    if value(explicit, acc, "the SCC accuracy of the engine") != "1.000000":
        fail("--acc 1.0 did not win over --lopt-level loose", explicit)

    tight = banner("tight", ["--lopt", "--lopt-level", "tight"])
    if value(tight, etol, "energy tolerance:") != "1.000000E-06":
        fail("--lopt-level tight does not report e_tol 1e-6 (T3)", tight)
    if value(tight, budget, "maximum number of steps:") != "200":
        fail("--lopt-level tight does not report clamp(5n, 200, 10000) = 200",
             tight)

    superloose = banner("superloose", ["--lopt", "--lopt-level", "superloose"])
    if value(superloose, etol, "energy tolerance:") != "5.000000E-04":
        fail("--lopt-level superloose does not report e_tol 5e-4 (T3)",
             superloose)

    explicit_budget = banner("tight-budget",
                             ["--lopt", "--lopt-level", "tight",
                              "--lopt-max-steps", "8"])
    if value(explicit_budget, budget, "maximum number of steps:") != "8":
        fail("--lopt-max-steps 8 did not win over the level-derived budget",
             explicit_budget)

    cycles = banner("tight-cycles",
                    ["--lopt", "--lopt-level", "tight", "--cycles", "4"])
    if value(cycles, budget, "maximum number of steps:") != "4":
        fail("--cycles 4 did not win over the level-derived budget", cycles)
    print("level axes: acc loose 2.0 (--acc 1.0 wins), e_tol loose/tight/"
          "superloose echoed, budget 200 = clamp, explicit 8 and 4 win")


# The startup banner (`src/header.f90`) is the only place where the product
# version, the engine version and the legal notices are printed together, and
# nothing covered it before this check (`xtb/Argparser_print_version` only ran
# `--version` and looked at the exit status).  `xtb_header` centres every line
# in one block of `BANNER_WIDTH` characters and indents the block by
# `BANNER_INDENT`, so the check pins the exact indentation of the wordmark and
# of every text line, both version numbers, the build date line (commit hash
# included), the legal notice, the order of the block - and that the old boxed
# `x T B` banner is gone.  The banner is four semantic blocks (wordmark,
# product description, authors, version and legal notice) separated by one
# blank line each; the blank line is part of the exact block below, so moving a
# line across a boundary, removing a line and rewording a line all fail.  The
# wordmark is compared **with** its trailing blanks, because those are part of
# the drawing.
#
# Three properties are asserted on their own, not only through the exact
# block: every printed banner line fits the 80 columns of a standard terminal,
# the wordmark is exactly the drawing below, and the build date line carries
# the short commit hash of the build, taken from the `provenance` block of the
# same binary's `--ljson` document.
BANNER_WIDTH = 58
BANNER_INDENT = 11
BANNER_ART_WIDTH = 24
BANNER_LIMIT = 80
BANNER_WORDMARK = [
    r" ___  __    _____ ____  ",
    r"/ / \ \ \_/  | |  | |_) ",
    r"\_\_/ /_/ \  |_|  |_|_) ",
]
BANNER_DESCRIPTION = [
    "OxTB - An Omni-functional PES Exploration and Exploitation",
    "Toolkit Based on xTB",
]
BANNER_AFFILIATION = [
    "Original Author - S. Grimme",
    "Mulliken Center for Theoretical Chemistry",
    "University of Bonn",
]
BANNER_TEXT = BANNER_DESCRIPTION + BANNER_AFFILIATION
BANNER_VERSION = [
    "* OxTB version 0.3.2 *",
    "* xTB version 6.7.1 *",
]
BANNER_LEGAL = [
    "Copyright (C) 2017-2020 Stefan Grimme,",
    "2019-2020 Sebastian Ehlert and contributors",
    "LGPL-3.0-or-later - NO WARRANTY;",
    "see `oxtb --license` for the full text",
]
# The dashed box of the banner this check replaced.
BANNER_OLD = "-" * 59
# What `cmake/CMakeLists.txt` bakes in when `git` cannot be asked for a commit.
BANNER_NO_COMMIT = "unknown-commit"


def centred_banner_line(text):
    """The line `xtb_header`/`banner_line` prints for `text`: centred in
    `BANNER_WIDTH` characters and indented by `BANNER_INDENT`."""
    offset = BANNER_INDENT + max(0, (BANNER_WIDTH - len(text.rstrip())) // 2)
    return " " * offset + text


def centred_wordmark(row):
    """The line `xtb_header` prints for one wordmark `row`.

    The rows are padded to `BANNER_ART_WIDTH` and share one offset, so the
    three of them start in the same column; the row keeps its trailing
    blanks."""
    offset = BANNER_INDENT + max(0, (BANNER_WIDTH - BANNER_ART_WIDTH) // 2)
    return " " * offset + row


def check_banner(exe, inp, rundir):
    """The OxTB banner: wordmark, versions, build date and legal notice.

    Mirrors `xtb/landscape-banner`.  `--version` prints the banner alone, so
    the whole block is compared line by line (including its indentation); a
    normal run has to print the same banner *and* still print the full
    disclaimer (`disclaimer`, the renamed `disclamer`), and `--license` has to
    keep printing the licence text.
    """
    status, output = run(exe, inp, ["--version"], rundir)
    if status != 0:
        fail("--version terminated with status %d" % status, output)
    lines = output.replace("\r\n", "\n").split("\n")

    compiled = [line.strip() for line in lines
                if line.strip().startswith("Compiled on ")]
    if len(compiled) != 1:
        fail("the banner has %d 'Compiled on' lines, expected exactly one"
             % len(compiled), output)
    match = re.fullmatch(
        r"Compiled on ([0-9]{4}/[0-9]{2}/[0-9]{2})(?: from (\S+))?",
        compiled[0])
    if not match:
        fail("the build date line %r is not 'Compiled on <date>' or "
             "'Compiled on <date> from <hash>'" % compiled[0], output)
    date, banner_commit = match.group(1), match.group(2)

    block = ([centred_wordmark(row) for row in BANNER_WORDMARK] + [""]
             + [centred_banner_line(text) for text in BANNER_DESCRIPTION]
             + [""]
             + [centred_banner_line(text) for text in BANNER_AFFILIATION]
             + [""]
             + [centred_banner_line(text) for text in BANNER_VERSION]
             + [centred_banner_line(compiled[0]), ""]
             + [centred_banner_line(text) for text in BANNER_LEGAL])
    # The product description line is stable while the wordmark is not, so it
    # is the anchor of the block: the wordmark starts `len(BANNER_WORDMARK)`
    # lines above it, one more for the blank line between the two blocks.
    anchor = centred_banner_line(BANNER_DESCRIPTION[0])
    try:
        start = lines.index(anchor) - len(BANNER_WORDMARK) - 1
    except ValueError:
        fail("the banner product line is missing (expected %r)" % anchor,
             output)
    if start < 0:
        fail("the banner wordmark is missing above %r" % anchor, output)

    # (b) The wordmark is exactly the drawing above, row for row.  Comparing
    # the rows verbatim is what catches an accidental edit to the drawing, a
    # stray leading blank on one row for instance, because it changes that
    # row's text as well as its alignment.
    #
    # An earlier version of this check instead required all three rows to
    # start in the same column.  That is wrong as a rule: the glyph tops are
    # deliberately inset or widened - the bar of the "O" sits inside its
    # sides, while the bar of the "T" is wider than its stem - so the rows do
    # not share a first column in a correct drawing.
    wordmark = [centred_wordmark(row) for row in BANNER_WORDMARK]
    if lines[start:start + len(wordmark)] != wordmark:
        fail("the wordmark is not the expected drawing\nexpected: %r\ngot:      %r"
             % (wordmark, lines[start:start + len(wordmark)]), output)

    # (a) Every printed banner line fits a standard terminal.
    for line in lines[start:start + len(block)]:
        if len(line) > BANNER_LIMIT:
            fail("the banner line %r is %d columns wide, more than %d"
                 % (line, len(line), BANNER_LIMIT), output)

    if lines[start:start + len(block)] != block:
        fail("the banner is not the expected centred block\nexpected: %r\ngot:      %r"
             % (block, lines[start:start + len(block)]), output)
    for line in lines:
        if line.strip() == BANNER_OLD:
            fail("the old boxed xTB banner is still printed", output)
    print("--version banner: %d centred lines, product 0.3.2, engine 6.7.1, "
          "build date %s, commit %s"
          % (len(block), date, banner_commit or "none"))

    status, output = run(exe, inp, ["--gfn", "2", "--sp", "--norestart",
                                    "--ljson"], rundir)
    if status != 0:
        fail("the GFN2 single point terminated with status %d" % status,
             output)
    normal = output.split("\n")
    for row in BANNER_WORDMARK:
        if centred_wordmark(row) not in normal:
            fail("the banner wordmark row %r is missing from a normal run"
                 % row, output)
    for text in BANNER_TEXT + BANNER_VERSION:
        if centred_banner_line(text) not in normal:
            fail("the banner line %r is missing from a normal run" % text,
                 output)
    if "Compiled on " not in output:
        fail("the banner of a normal run has no build date", output)
    if "GNU Lesser General Public License" not in output:
        fail("the disclaimer is no longer printed at startup", output)

    # The licence notice names the program, not the engine: it said `xtb`
    # while the line above it said `see 'oxtb --license'`.  Only the two
    # records that begin with the name are pinned, not the notice as a whole,
    # and the old spelling is matched as a whole line because `xtb is free
    # software:` is a substring of `oxtb is free software:`.
    licence = [line.strip() for line in normal]
    for text in ("oxtb is free software: you can redistribute it and/or "
                 "modify it under",
                 "oxtb is distributed in the hope that it will be useful,"):
        if text not in licence:
            fail("the licence notice does not name the program `oxtb`: %r is "
                 "missing" % text, output)
    for text in ("xtb is free software: you can redistribute it and/or "
                 "modify it under",
                 "xtb is distributed in the hope that it will be useful,"):
        if text in licence:
            fail("the licence notice still names the engine `xtb`: %r" % text,
                 output)

    # The build date line carries the commit hash of the build.  The same
    # binary reports the hash in the `provenance` block of its `--ljson`
    # document, so the check compares the two instead of trusting the format;
    # a build without a hash (no `git` at configure time) must fall back to
    # the date alone rather than print an empty token.
    build_commit = read_json(rundir, "landscape_sp.json")["provenance"]["commit"]
    if build_commit and build_commit != BANNER_NO_COMMIT:
        expected = build_commit[:7]
        if banner_commit != expected:
            fail("the build date line carries %r, expected the commit hash %r"
                 % (banner_commit, expected), output)
    elif banner_commit is not None:
        fail("the build carries no commit hash but the banner claims %r"
             % banner_commit, output)

    status, output = run(exe, inp, ["--license"], rundir)
    if status != 0:
        fail("--license terminated with status %d" % status, output)
    if "but WITHOUT ANY WARRANTY" not in output:
        fail("--license no longer prints the warranty disclaimer", output)
    if "GNU Lesser General Public License" not in output:
        fail("--license no longer prints the licence", output)
    if "oxtb is free software:" not in output:
        fail("--license still names the program `xtb`", output)
    print("banner also printed by a normal run, disclaimer and --license "
          "intact")


# The box lines of the startup engine banner.  Every run that has an engine
# prints exactly one of them and a run without an engine prints none.
ENGINE_BANNER_LABELS = ("Engine: ", "Saddle search: ", "IRC integrator: ")

# (label, arguments, the box line the run has to print): every engine name of
# every selector, the builtin default of each task, and the native engines.
# The runs are deliberately cut short (`--l*max-steps 1`) where the fixture is
# not a real minimum/saddle: the banner is printed before the first geometry is
# evaluated, which is the property under test, so the outcome of the run is not
# asserted here.
ENGINE_BANNER_CASES = (
    ("lopt-default", ["--lopt"], "Engine: lbfgs (default)"),
    ("lopt-sd", ["--lopt", "--lopt-engine", "sd"], "Engine: sd (explicit)"),
    ("lopt-cg", ["--lopt", "--lopt-engine", "cg"], "Engine: cg (explicit)"),
    ("lopt-lbfgs", ["--lopt", "--lopt-engine", "lbfgs"],
     "Engine: lbfgs (explicit)"),
    ("lopt-bfgs", ["--lopt", "--lopt-engine", "bfgs"],
     "Engine: bfgs (explicit)"),
    ("lopt-diis", ["--lopt", "--lopt-engine", "diis"],
     "Engine: diis (explicit)"),
    ("lopt-rfo", ["--lopt", "--lopt-engine", "rfo", "--lopt-max-steps", "1"],
     "Engine: rfo (explicit)"),
    ("scan-default", ["--lscan", "--lscan-coord", "distance", "--lscan-atoms",
                       "1,2", "--lscan-step", "0.1", "--lscan-nstep", "1"],
     "Engine: lbfgs (default)"),
    ("scan-cg", ["--lscan", "--lscan-method", "cg", "--lscan-coord",
                  "distance", "--lscan-atoms", "1,2", "--lscan-step", "0.1",
                  "--lscan-nstep", "1"], "Engine: cg (explicit)"),
    ("ts-default", ["--lts", "--lts-max-steps", "1"],
     "Saddle search: prfo (default)"),
    ("ts-prfo", ["--lts", "--lts-method", "prfo", "--lts-max-steps", "1"],
     "Saddle search: prfo (explicit)"),
    ("ts-bprfo", ["--lts", "--lts-method", "bprfo", "--lts-max-steps", "1"],
     "Saddle search: bprfo (explicit)"),
    ("ts-dimer", ["--lts", "--lts-method", "dimer", "--lts-max-steps", "1"],
     "Saddle search: dimer (explicit)"),
    ("ts-neb", ["--lts", "--lts-method", "neb", "--lts-max-steps", "1"],
     "Saddle search: neb (explicit)"),
    ("ts-cineb", ["--lts", "--lts-method", "cineb", "--lts-max-steps", "1"],
     "Saddle search: cineb (explicit)"),
    ("ts-autoneb", ["--lts", "--lts-method", "autoneb", "--lts-max-steps",
                     "1"], "Saddle search: autoneb (explicit)"),
    ("ts-string", ["--lts", "--lts-method", "string", "--lts-max-steps",
                    "1"], "Saddle search: string (explicit)"),
    ("ts-gsm", ["--lts", "--lts-method", "gsm", "--lts-max-steps", "1"],
     "Saddle search: gsm (explicit)"),
    ("irc-default", ["--lirc", "--lirc-max-steps", "1"],
     "IRC integrator: hpc (default)"),
    ("irc-gs", ["--lirc", "--lirc-method", "gs", "--lirc-max-steps", "1"],
     "IRC integrator: gs (explicit)"),
    ("irc-hpc", ["--lirc", "--lirc-method", "hpc", "--lirc-max-steps", "1"],
     "IRC integrator: hpc (explicit)"),
    ("irc-lqa", ["--lirc", "--lirc-method", "lqa", "--lirc-max-steps", "1"],
     "IRC integrator: lqa (explicit)"),
    ("irc-eulerpc", ["--lirc", "--lirc-method", "eulerpc",
                      "--lirc-max-steps", "1"],
     "IRC integrator: eulerpc (explicit)"),
    ("native-ancopt", ["--opt", "--method", "ancopt", "--norestart"],
     "Engine: ancopt (explicit)"),
    ("native-lancopt", ["--opt", "--method", "lancopt", "--norestart"],
     "Engine: lancopt (explicit)"),
    ("native-fire", ["--opt", "--method", "fire", "--norestart"],
     "Engine: fire (explicit)"),
    ("native-pbc-lbfgs", ["--opt", "--method", "pbc_lbfgs", "--norestart"],
     "Engine: pbc_lbfgs (explicit)"),
    ("native-ohess", ["--ohess", "--method", "fire", "--norestart"],
     "Engine: fire (explicit)"),
)

# Run types with no engine at all: they must print no box line.
ENGINE_FREE_CASES = (
    ("sp", ["--sp"]),
    ("hess", ["--hess"]),
    ("lvib", ["--lvib"]),
    ("rigid-scan", ["--lscan", "--lscan-mode", "rigid", "--lscan-coord",
                    "distance", "--lscan-atoms", "1,2", "--lscan-step",
                    "0.1", "--lscan-nstep", "1"]),
)


def check_engine_banner(exe, inp, rundir):
    """One engine announcement per run, for every engine, early and unique.

    The engine used to be echoed by each task's print block, after the
    always-done single point, and only as a `[NOTICE]` for the builtin
    default - an explicitly named engine was never announced and `--lts-method`
    had no echo at all.  The startup banner replaces that: this check pins

    (1) every name of the three selectors and every native `--method` engine,
        with `(explicit)` when the command line named it and `(default)` when
        the builtin default runs - so a name that stops reaching the banner,
        and a run that stops distinguishing the two, both turn it red;
    (2) exactly one box line per run, so the announcement cannot be doubled
        again (the retired `[NOTICE] --lopt-engine`/`--lirc-method` lines are
        asserted absent as well);
    (3) nothing at all for a run with no engine (`--sp`, `--hess`, `--lvib`,
        rigid scan);
    (4) the position: the box is printed before the first geometry is
        evaluated, i.e. before the `SETUP` block and the first `TOTAL ENERGY`
        of the run.
    """
    for label, args, expected in ENGINE_BANNER_CASES:
        workdir = os.path.join(rundir, label)
        shutil.rmtree(workdir, ignore_errors=True)
        # the cut-short searches are allowed to fail: the banner is printed
        # before any geometry is evaluated
        _, output = run(exe, inp, args + ["--gfn", "2"], workdir)
        if expected not in output:
            fail("%s: the engine banner does not say %r" % (label, expected),
                 output)
        lines = [line for line in output.split("\n")
                 if any(text in line for text in ENGINE_BANNER_LABELS)]
        if len(lines) != 1:
            fail("%s: %d engine banner lines, expected one\n%r"
                 % (label, len(lines), lines), output)
        for retired in ("[NOTICE] --lopt-engine default",
                        "[NOTICE] --lirc-method default"):
            if retired in output:
                fail("%s: the retired engine notice is back (%r)"
                     % (label, retired), output)

    for label, args in ENGINE_FREE_CASES:
        workdir = os.path.join(rundir, label)
        shutil.rmtree(workdir, ignore_errors=True)
        status, output = run(exe, inp, args + ["--gfn", "2", "--norestart"],
                             workdir)
        if status != 0:
            fail("%s terminated with status %d" % (label, status), output)
        for text in ENGINE_BANNER_LABELS:
            if text in output:
                fail("%s has no engine but announces one (%r)"
                     % (label, text), output)

    # The banner has to precede the first evaluated geometry.  A `--lopt` run
    # is the case that used to print its echo 16 lines below the first SCC
    # summary: the box must be above both the `SETUP` block and the first
    # `TOTAL ENERGY` of that run.
    workdir = os.path.join(rundir, "early")
    shutil.rmtree(workdir, ignore_errors=True)
    status, output = run(exe, inp, ["--gfn", "2", "--lopt", "--norestart"],
                         workdir)
    if status != 0:
        fail("the early-position probe terminated with status %d" % status,
             output)
    lines = output.split("\n")
    banner = [i for i, line in enumerate(lines) if "Engine: " in line]
    setup = [i for i, line in enumerate(lines) if "SETUP" in line]
    energy = [i for i, line in enumerate(lines) if "TOTAL ENERGY" in line]
    if len(banner) != 1 or not setup or not energy:
        fail("the early-position probe printed %d banners, %d SETUP and %d "
             "energy lines" % (len(banner), len(setup), len(energy)), output)
    if not banner[0] < min(setup[0], energy[0]):
        fail("the engine banner is printed after the first geometry "
             "evaluation (banner %d, SETUP %d, energy %d)"
             % (banner[0], setup[0], energy[0]), output)

    print("engine banner: %d engine cases, %d engine-free cases, announced "
          "before the first geometry"
          % (len(ENGINE_BANNER_CASES), len(ENGINE_FREE_CASES)))


def check_programname(exe, inp, rundir):
    """Both runtime termination messages name the program `oxtb`.

    `src/prog/primary.f90` hands the program name to `mctc_init`, and
    `src/mctc/error.f90` prints it as the `normal termination of ...` /
    `abnormal termination of ...` line on **stderr**.  Both messages said
    `xtb` until the rename, and nothing covered them (`--help` and the banner
    print other text), so reverting the one literal in `primary.f90` left every
    test green.  This check pins both directions: a successful run has to print
    `normal termination of oxtb`, a deliberately failing run (a missing input
    file) `abnormal termination of oxtb`, and neither may print the old name.
    """
    status, stdout, stderr = run_split(exe, inp, ["--gfn", "2", "--sp",
                                                  "--norestart"], rundir)
    if status != 0:
        fail("the GFN2 single point terminated with status %d" % status,
             stdout + stderr)
    if "normal termination of oxtb" not in stderr:
        fail("the successful run does not report 'normal termination of "
             "oxtb' on stderr", stderr)
    if "termination of xtb" in stderr:
        fail("the successful run still names the program 'xtb' on stderr",
             stderr)
    print("successful run: 'normal termination of oxtb' on stderr")

    # A missing input file makes `xtbMain` raise through the environment, which
    # `src/mctc/error.f90` turns into the abnormal termination message; the
    # file is removed first so a leftover from an earlier run cannot hide it.
    missing = os.path.join(os.path.abspath(rundir), "oxtb-no-such-input.xyz")
    if os.path.exists(missing):
        os.remove(missing)
    status, stdout, stderr = run_split(exe, missing, ["--gfn", "2", "--sp",
                                                      "--norestart"], rundir)
    if status == 0:
        fail("a run with a missing input file exited 0", stdout + stderr)
    if "abnormal termination of oxtb" not in stderr:
        fail("the failing run does not report 'abnormal termination of "
             "oxtb' on stderr", stderr)
    if "termination of xtb" in stderr:
        fail("the failing run still names the program 'xtb' on stderr",
             stderr)
    print("failing run: 'abnormal termination of oxtb' on stderr")


#: A value no finite-number parser rejects but no comparison catches either:
#: `read(sec,*,iostat=)` accepts `nan` and `inf`, `<`/`>` are false against NaN,
#: and an SCC accuracy or a force constant of NaN poisons every number the run
#: prints.  The landscape suite therefore refuses a non-finite value where the
#: setting is consumed (the xcontrol block merges into the same state after the
#: command line, so a parse-site check would miss `$thermo temp=nan`).
NONFINITE_CASES = [
    ("--lvib --lvib-temp nan", ["--lvib", "--lvib-temp", "nan"], "--lvib-temp"),
    ("--lvib --lvib-freq-scal inf",
     ["--lvib", "--lvib-freq-scal", "inf"], "--lvib-freq-scal"),
    ("--lvib --lvib-imag-cutoff nan",
     ["--lvib", "--lvib-imag-cutoff", "nan"], "--lvib-imag-cutoff"),
    ("--lvib --lvib-hess-step nan",
     ["--lvib", "--lvib-hess-step", "nan"], "--lvib-hess-step"),
    ("--lscan --lscan-fc nan",
     ["--lscan", "--lscan-coord", "distance", "--lscan-atoms", "1,2",
      "--lscan-step", "0.1", "--lscan-nstep", "2", "--lscan-fc", "nan"],
     "--lscan-fc"),
]


#: Labels of the settings echo: the task prints the value it was given, so a
#: NaN here is the *input* being named, not a number the run computed.
ECHO_LABELS = ("program call", "Hessian step:", "Hessian SCC accuracy:")


def nan_lines(output, extra=()):
    """Result lines that carry a NaN/Inf number.

    The diagnostic itself names the value (`NaN`), so the `-1-`/`-2-` lines of
    the environment log, the echoed `program call` line, the settings echo of
    `ECHO_LABELS` and the termination report are filtered out; what remains is
    a number the run computed.
    """
    hits = []
    for line in output.splitlines():
        if not re.search(r"(?i)\b(nan|inf)\b", line):
            continue
        stripped = line.lstrip()
        if stripped.startswith("-1-") or stripped.startswith("-2-"):
            continue
        if "ERROR" in line.upper():
            continue
        if any(label in line for label in ECHO_LABELS):
            continue
        if any(token in line for token in extra):
            continue
        hits.append(line)
    return hits


def check_nonfinite(exe, inp, rundir):
    """A non-finite value for a numeric flag is reported, never computed with.

    Mirrors `xtb/landscape-nonfinite-value`.  Every flag below accepts `nan` or
    `inf` at the parser, and in a release build the value used to reach the
    calculation: `--lvib-temp nan` printed a complete table of NaN
    (`H(0)-H(T)`, `T*S`, `G(T)`) and exited **0** - a silent wrong answer - and
    `--acc nan` stored a NaN SCC accuracy.  The five landscape flags are
    refused where the merged setting is consumed, so the run stops with a
    message that names flag and value; `--acc` keeps the default and is a
    warning, fatal under `--strict`.
    """
    for label, args, flag in NONFINITE_CASES:
        directory = os.path.join(rundir, "nonfinite" +
                                 re.sub(r"[^A-Za-z0-9]", "_", label))
        shutil.rmtree(directory, ignore_errors=True)
        status, output = run(exe, inp, ["--gfn", "2", "--norestart"] + args,
                             directory)
        if status == 0:
            fail("'%s' exited 0 instead of refusing the value" % label, output)
        if flag not in output:
            fail("'%s' does not name the flag %s" % (label, flag), output)
        hits = nan_lines(output, extra=(flag, "nan", "inf"))
        if hits:
            fail("'%s' computed with the value: %s"
                 % (label, hits[0].strip()), output)

    # The same flags with a usable value must still run, so the check cannot be
    # satisfied by rejecting the flag outright.
    directory = os.path.join(rundir, "finite")
    shutil.rmtree(directory, ignore_errors=True)
    status, output = run(exe, inp, ["--gfn", "2", "--norestart", "--lvib",
                                    "--lvib-temp", "300"], directory)
    if status != 0:
        fail("--lvib --lvib-temp 300 failed", output)

    # `--acc nan`: lax keeps the default accuracy (the energy is the one of a
    # run without the flag) and warns; `--strict` makes the same value fatal.
    energy_re = re.compile(r"TOTAL ENERGY +(-?[0-9]+\.[0-9]+)")
    plain_dir = os.path.join(rundir, "acc-plain")
    shutil.rmtree(plain_dir, ignore_errors=True)
    status, output = run(exe, inp, ["--gfn", "2", "--sp", "--norestart"],
                         plain_dir)
    if status != 0:
        fail("the reference single point failed", output)
    reference = energy_re.search(output)
    if not reference:
        fail("no TOTAL ENERGY in the reference single point", output)
    lax_dir = os.path.join(rundir, "acc-lax")
    shutil.rmtree(lax_dir, ignore_errors=True)
    status, output = run(exe, inp, ["--gfn", "2", "--sp", "--norestart",
                                    "--acc", "nan"], lax_dir)
    if status != 0:
        fail("--acc nan is fatal without --strict", output)
    if "--acc needs a finite accuracy" not in output:
        fail("--acc nan was not reported", output)
    energy = energy_re.search(output)
    if not energy:
        fail("no TOTAL ENERGY with --acc nan", output)
    if energy.group(1) != reference.group(1):
        fail("--acc nan changed the energy (%s, default %s)"
             % (energy.group(1), reference.group(1)), output)
    strict_dir = os.path.join(rundir, "acc-strict")
    shutil.rmtree(strict_dir, ignore_errors=True)
    status, output = run(exe, inp, ["--gfn", "2", "--sp", "--norestart",
                                    "--strict", "--acc", "nan"], strict_dir)
    if status == 0:
        fail("--strict --acc nan exited 0", output)
    print("non-finite values reported; the finite ones still run")


def check_ir_fatal(exe, inp, rundir):
    """A fatal path of the `ir` subcommand exits non-zero.

    Mirrors `xtb/landscape-ir-fatal-status`.  `oxtb ir` without `--born` used
    an **unallocated** `intent(out)` file name (`src/prog/irmod.f90`), and the
    corrupted-file branch printed `... aborting` and executed a bare `stop`,
    which returns exit status **0**: a fatal error that a script cannot see.
    The name now defaults to `born.out`, a missing file is reported before the
    reader fails on its own, and the unreadable-file branch goes through
    `env%error` + `env%checkpoint`, whose termination status is 1.
    """
    def run_ir(directory):
        """Run the `ir` subcommand: the subcommand word comes first."""
        os.makedirs(directory, exist_ok=True)
        proc = subprocess.run(
            [exe, "ir", inp],
            cwd=directory,
            capture_output=True,
            text=False,
        )
        return proc.returncode, decode(proc.stdout) + decode(proc.stderr)

    # (1) a truncated dipole derivative file: the read fails and the run aborts
    broken = os.path.join(rundir, "broken")
    shutil.rmtree(broken, ignore_errors=True)
    os.makedirs(broken, exist_ok=True)
    with open(os.path.join(broken, "born.out"), "w") as handle:
        handle.write("1.0 0.0 0.0\n")
    status, output = run_ir(broken)
    if status == 0:
        fail("a truncated born.out exited 0", output)
    if "born.out" not in output or "corrupted" not in output:
        fail("the truncated born.out was not reported", output)

    # (2) no born.out at all: reported against the flag, not by a raw abort
    missing = os.path.join(rundir, "missing")
    shutil.rmtree(missing, ignore_errors=True)
    os.makedirs(missing, exist_ok=True)
    status, output = run_ir(missing)
    if status == 0:
        fail("a missing born.out exited 0", output)
    if "--born" not in output:
        fail("a missing born.out does not name --born", output)
    print("ir: both fatal paths exit non-zero")


#: The category keys of the grouped `--help` index, in print order; the
#: unknown-topic error has to list every one of them.
HELP_CATEGORY_KEYS = (
    "runtypes", "suite", "levels", "minimisation", "scan", "saddle",
    "band", "irc", "frequencies", "hamiltonian", "electronic", "external",
    "solvation", "analysis", "output", "workflow", "performance",
)

#: One sampled topic per category, as (topic, first line, distinctive phrase
#: of the full text).  `--method` is the longest option block of the page.
HELP_TOPIC_SAMPLES = (
    ("--method", "--method NAME", "canonical name of --lopt-engine"),
    ("--hess", "--hess", "numerical hessian calculation on input geometry"),
    ("--lopt", "--lopt", "minimize the geometry with the landscape task suite"),
    ("--lopt-engine", "--lopt-engine NAME",
     "minimisation algorithm of the landscape suite"),
    ("--lscan-coord", "--lscan-coord TYPE",
     "type of the first scanned coordinate"),
    ("--lts-trust", "--lts-trust VALUE", "initial trust radius of the step"),
    ("--lts-product", "--lts-product NAME",
     "product endpoint of a chain-of-states search"),
    ("--lirc-step", "--lirc-step VALUE", "step length in bohr*amu^1/2"),
    ("--lvib-temp", "--lvib-temp VALUE",
     "temperature of the thermodynamic table in K"),
    ("--gfn", "--gfn INT", "specify parametrisation of GFN-xTB"),
    ("--chrg", "-c, --chrg, --charge INT", "specify molecular charge as INT"),
    ("--tm", "--tm", "external Turbomole program for the energy"),
    ("--alpb", "--alpb SOLVENT [STATE]", "analytical linearized Poisson-Boltzmann"),
    ("--pop", "--pop", "Mulliken population analysis"),
    ("--json", "--json", "write xtbout.json file"),
    ("--input", "-I, --input FILE", "input source for xcontrol(7) instructions"),
    ("--parallel", "-P, --parallel INT", "number of OpenMP threads of the run"),
)


def run_help(exe, rundir, args):
    """Run `oxtb` without a geometry and return (status, combined output)."""
    os.makedirs(rundir, exist_ok=True)
    proc = subprocess.run([exe] + args, cwd=rundir, capture_output=True,
                          text=False)
    return proc.returncode, decode(proc.stdout) + decode(proc.stderr)


def index_spellings(text):
    """The spelling set of a rendered grouped index page.

    An index entry is one line `"  " + spellings + "  " + summary`; the
    spellings never contain two consecutive spaces, so the first double space
    separates them from the gloss.  The token and `--[no]x` rules are the ones
    of `landscape_help_check.py`, so the rendered set is comparable with the
    parser's set.
    """
    from landscape_help_check import TOKEN, SPLIT

    flags = set()
    for line in text.splitlines():
        if not line.startswith("  -"):
            continue
        body = line[2:]
        cut = body.find("  ")
        spec = body[:cut] if cut >= 0 else body
        for token in SPLIT.split(spec.strip()):
            if not TOKEN.match(token):
                continue
            if "[no]" in token:
                flags.add(token.replace("[no]", ""))
                flags.add(token.replace("[no]", "no"))
            else:
                flags.add(token)
    return flags


def check_help_topic(exe, inp, rundir, main_src):
    """The grouped `--help` index and its `--help <topic>` dispatch.

    The expected spelling set is taken from the **parser** (`src/prog/main.F90`,
    through the extractor of `landscape_help_check.py`), not from the help
    source, so the check is falsifiable: dropping a spelling from the index
    turns it red even though the help gate and the index would still agree with
    each other.  The samples cover one option of every category, `--method`
    (the longest block) included, and the unknown topic has to name the
    categories and exit non-zero.
    """
    from landscape_help_check import HELP_ALLOWLIST, implemented_flags

    with open(main_src) as handle:
        implemented = implemented_flags(handle.read())
    if len(implemented) != 249:
        fail("the parser implements %d spellings, expected 249"
             % len(implemented))
    # `--siman` and `--gmd` are parsed but are not registered run types; they
    # are deliberately absent from the help table (see the allowlist in
    # `landscape_help_check.py`), so the index must not show them and must not
    # be expected to.
    expected = implemented - set(HELP_ALLOWLIST)

    status, page = run_help(exe, rundir, ["--help"])
    if status != 0:
        fail("oxtb --help terminated with status %d" % status, page)
    shown = index_spellings(page)
    missing = sorted(expected - shown)
    phantom = sorted(shown - expected)
    if missing:
        fail("the index misses %d of the %d implemented spellings: %s"
             % (len(missing), len(expected), " ".join(missing)), page)
    if phantom:
        fail("the index shows %d spellings the parser does not implement: %s"
             % (len(phantom), " ".join(phantom)), page)
    print("index: %d spellings, the %d implemented ones, none phantom"
          % (len(shown), len(expected)))

    for topic, first, phrase in HELP_TOPIC_SAMPLES:
        status, out = run_help(exe, rundir, ["--help", topic])
        if status != 0:
            fail("oxtb --help %s terminated with status %d" % (topic, status),
                 out)
        lines = out.splitlines()
        if not lines or lines[0].strip() != first:
            fail("oxtb --help %s did not start with %r" % (topic, first), out)
        if phrase not in out:
            fail("oxtb --help %s misses the text %r" % (topic, phrase), out)
        if len(lines) > 40:
            fail("oxtb --help %s printed %d lines, the whole index?"
                 % (topic, len(lines)), out)
        if "Native run types" in out or "Landscape suite tasks" in out:
            fail("oxtb --help %s printed the whole index" % topic, out)
    print("topics: %d sampled options print their own full text"
          % len(HELP_TOPIC_SAMPLES))

    # A bare name and an upper case spelling are the same topic.
    status, bare = run_help(exe, rundir, ["--help", "gfn"])
    status2, upper = run_help(exe, rundir, ["--help", "--GFN"])
    if status != 0 or status2 != 0 or bare != upper:
        fail("--help gfn and --help --GFN differ", bare + upper)
    print("topics: a bare and an upper case spelling resolve alike")

    # One category: its heading and its options, nothing of the other groups.
    status, cat = run_help(exe, rundir, ["--help", "saddle"])
    if status != 0:
        fail("oxtb --help saddle terminated with status %d" % status, cat)
    lines = cat.splitlines()
    if not lines or lines[0].strip() != "Saddle search options (--lts)":
        fail("oxtb --help saddle did not start with the heading", cat)
    for flag in ("--lts-method", "--lts-trust", "--lts-escape-step"):
        if flag not in cat:
            fail("oxtb --help saddle misses %s" % flag, cat)
    if "Native run types" in cat or "--gfn" in cat:
        fail("oxtb --help saddle printed another category", cat)
    print("categories: --help saddle prints one group (%d lines)" % len(lines))

    # An unknown topic is an error that lists the categories.
    status, bad = run_help(exe, rundir, ["--help", "nonsense"])
    if status == 0:
        fail("oxtb --help nonsense exited 0", bad)
    if "Unknown help topic" not in bad or "nonsense" not in bad:
        fail("the unknown topic error does not name the topic", bad)
    for key in HELP_CATEGORY_KEYS:
        if key not in bad:
            fail("the unknown topic error does not list %r" % key, bad)
    print("unknown topic: status %d, the 17 categories listed" % status)


def run_help_raw(exe, rundir, args, overrides):
    """Run `oxtb` with a controlled colour environment; return stderr-free
    stdout as raw bytes.

    The page goes to stdout alone (`normal termination` and the error trace
    of an unknown topic go to stderr), so the byte comparison reads stdout.
    `CLICOLOR_FORCE` and `NO_COLOR` are removed from the inherited environment
    first, so an ambient setting cannot turn the check red or green for a
    foreign reason; `overrides` then sets the ones the case needs.
    """
    os.makedirs(rundir, exist_ok=True)
    env = dict(os.environ)
    env.pop("CLICOLOR_FORCE", None)
    env.pop("NO_COLOR", None)
    env.update(overrides)
    proc = subprocess.run([exe] + args, cwd=rundir, capture_output=True,
                          text=False, env=env)
    return proc.returncode, proc.stdout


#: A basic ANSI SGR run: the help page uses no 256-colour codes (see
#: `docs/_simplify-recon/124-help-colour.md`).
HELP_SGR = re.compile(rb"\x1b\[[0-9;]*m")


def check_help_colour(exe, inp, rundir):
    """The `--help` page is coloured by `auto`, and the colour is additive.

    A pipe is not a terminal, so the default page must carry no escape byte -
    that is what keeps every byte-for-byte capture of `--help` in the repo
    valid.  `CLICOLOR_FORCE` forces colour even into a pipe and must change
    **only** escapes: the stripped coloured page is byte-identical to the
    plain one, for the index, for a topic page and for the unknown-topic
    error.  The precedence is `CLICOLOR_FORCE` > `NO_COLOR` > terminal.
    """
    cases = (
        ("the index", ["--help"]),
        ("a topic page", ["--help", "--method"]),
        ("a category", ["--help", "saddle"]),
        ("the unknown-topic error", ["--help", "nonsense"]),
    )
    for label, args in cases:
        status, plain = run_help_raw(exe, rundir, args, {})
        if b"\x1b" in plain:
            fail("%s: the plain page carries an escape byte" % label,
                 decode(plain))
        forced_status, coloured = run_help_raw(
            exe, rundir, args, {"CLICOLOR_FORCE": "1"})
        if forced_status != status:
            fail("%s: forcing colour changed the exit status (%d -> %d)"
                 % (label, status, forced_status))
        escapes = len(HELP_SGR.findall(coloured))
        if escapes == 0:
            fail("%s: CLICOLOR_FORCE=1 produced no escape sequence" % label,
                 decode(coloured))
        stripped = HELP_SGR.sub(b"", coloured)
        if stripped != plain:
            fail("%s: stripping the colour is not the plain page" % label,
                 decode(stripped))
        print("%s: %d SGR runs, the stripped page is byte-identical"
              % (label, escapes))

    # The precedence: CLICOLOR_FORCE wins over NO_COLOR (the deliberate
    # override beats the plain opt-out, as `uv` documents it), NO_COLOR
    # disables colour without a pipe, and a CLICOLOR_FORCE of `0` or an
    # empty value is not a force.
    _, both = run_help_raw(exe, rundir, ["--help"],
                           {"NO_COLOR": "1", "CLICOLOR_FORCE": "1"})
    if b"\x1b" not in both:
        fail("CLICOLOR_FORCE=1 must win over NO_COLOR=1")
    _, nocolor = run_help_raw(exe, rundir, ["--help"], {"NO_COLOR": "1"})
    if b"\x1b" in nocolor:
        fail("NO_COLOR=1 must disable colour")
    _, zero = run_help_raw(exe, rundir, ["--help"], {"CLICOLOR_FORCE": "0"})
    if b"\x1b" in zero:
        fail("CLICOLOR_FORCE=0 is not a force")
    _, empty = run_help_raw(exe, rundir, ["--help"], {"CLICOLOR_FORCE": ""})
    if b"\x1b" in empty:
        fail("an empty CLICOLOR_FORCE is not a force")
    print("precedence: CLICOLOR_FORCE > NO_COLOR > terminal, asserted")


def main():
    if len(sys.argv) < 5:
        fail("usage: landscape_cli_check.py MODE EXE INPUT RUNDIR [ARGS]")
    mode, exe, inp, rundir = sys.argv[1:5]
    rest = sys.argv[5:]
    if mode == "energy":
        check_energy(exe, inp, rundir)
    elif mode == "acc":
        check_acc(exe, inp, rundir)
    elif mode == "hessian":
        check_hessian(exe, inp, rundir)
    elif mode == "scf":
        check_scf(exe, inp, rundir)
    elif mode == "scanfail":
        check_scanfail(exe, inp, rundir)
    elif mode == "vib":
        check_vib(exe, inp, rundir)
    elif mode == "vibnamespace":
        check_vib_namespace(exe, inp, rundir)
    elif mode == "vibetemp":
        check_vib_etemp(exe, inp, rundir)
    elif mode == "fix":
        check_fix(exe, inp, rundir)
    elif mode == "scanblock":
        check_scanblock(exe, inp, rundir)
    elif mode == "saddlefix":
        check_saddlefix(exe, inp, rundir)
    elif mode == "rejectfix":
        check_rejectfix(exe, inp, rundir)
    elif mode == "thermo":
        check_thermo(exe, inp, rundir)
    elif mode == "thermoknobs":
        if len(rest) < 1:
            fail("the thermoknobs mode needs the saddle geometry")
        check_thermo_knobs(exe, inp, rundir, rest[0])
    elif mode == "optnamespace":
        check_opt_namespace(exe, inp, rundir)
    elif mode == "ljsondefault":
        check_ljson_default(exe, inp, rundir)
    elif mode == "ts":
        check_ts(exe, inp, rundir)
    elif mode == "tsetemp":
        check_ts_etemp(exe, inp, rundir)
    elif mode == "tsminimum":
        check_ts_minimum(exe, inp, rundir)
    elif mode == "tsescape":
        check_ts_escape(exe, inp, rundir)
    elif mode == "tsnamespace":
        check_ts_namespace(exe, inp, rundir)
    elif mode == "dimer":
        check_dimer(exe, inp, rundir)
    elif mode == "maxsteps":
        check_maxsteps(exe, inp, rundir)
    elif mode == "scanhess":
        check_scanhess(exe, inp, rundir)
    elif mode == "irc":
        check_irc(exe, inp, rundir)
    elif mode == "ircknobs":
        check_irc_knobs(exe, inp, rundir)
    elif mode == "ircpc":
        check_irc_pc(exe, inp, rundir)
    elif mode == "ircosc":
        check_irc_osc(exe, inp, rundir)
    elif mode == "ircswap":
        check_irc_swap(exe, inp, rundir)
    elif mode == "irctermination":
        check_irc_termination(exe, inp, rundir)
    elif mode == "autonebdepth":
        if len(rest) < 1:
            fail("the autonebdepth mode needs the product geometry")
        check_autoneb_depth_limit(exe, inp, rundir, rest[0])
    elif mode == "enumcase":
        if len(rest) < 1:
            fail("the enumcase mode needs the second (transition state) input")
        check_enum_case(exe, inp, rundir, rest[0])
    elif mode == "strictunknown":
        check_strict_unknown(exe, inp, rundir)
    elif mode == "flagvalue":
        check_flag_value(exe, inp, rundir)
    elif mode == "nativelevel":
        check_native_level(exe, inp, rundir)
    elif mode == "parallel":
        check_parallel(exe, inp, rundir)
    elif mode == "positional":
        check_positional(exe, inp, rundir)
    elif mode == "enginebanner":
        check_engine_banner(exe, inp, rundir)
    elif mode == "color":
        check_color(exe, inp, rundir)
    elif mode == "gfn":
        check_gfn(exe, inp, rundir)
    elif mode == "dockgfn":
        if len(rest) < 1:
            fail("the dockgfn mode needs the second geometry")
        check_dock_gfn(exe, inp, rundir, rest[0])
    elif mode == "echo":
        if len(rest) < 1:
            fail("the echo mode needs a case name")
        check_echo(exe, inp, rundir, rest[0])
    elif mode == "method":
        if len(rest) < 3:
            fail("the method mode needs the HCN, NH3 and NH3 product fixtures")
        check_method(exe, inp, rundir, rest[0], rest[1], rest[2])
    elif mode == "methodnative":
        if len(rest) < 1:
            fail("the methodnative mode needs the 501 atom fixture")
        check_method_native(exe, inp, rundir, rest[0])
    elif mode == "methodtable":
        if len(rest) < 3:
            fail("the methodtable mode needs the HCN, NH3 and NH3 product fixtures")
        check_method_table(exe, inp, rundir, rest[0], rest[1], rest[2])
    elif mode == "nativeoptsym":
        check_native_opt_symbols(exe, inp, rundir)
    elif mode == "optknobs":
        check_opt_knobs(exe, inp, rundir)
    elif mode == "hessprogress":
        if len(rest) < 1:
            fail("the hessprogress mode needs the stretched water fixture")
        check_hess_progress(exe, inp, rundir, rest[0])
    elif mode == "optrecalc":
        if len(rest) < 1:
            fail("the optrecalc mode needs the stretched water fixture")
        check_opt_recalc(exe, inp, rundir, rest[0])
    elif mode == "scanknobs":
        check_scan_knobs(exe, inp, rundir)
    elif mode == "hessstepalias":
        if len(rest) < 1:
            fail("the hessstepalias mode needs the transition state input")
        check_hess_step_alias(exe, inp, rundir, rest[0])
    elif mode == "restartknobs":
        check_restart_knobs(exe, inp, rundir)
    elif mode == "hessused":
        if len(rest) < 3:
            fail("the hessused mode needs the HCN, NH3 and NH3 product "
                 "fixtures")
        check_hessian_step_used(exe, inp, rundir, rest[0], rest[1], rest[2])
    elif mode == "lvibrestart":
        check_lvib_restart_warm(exe, inp, rundir)
    elif mode == "optalias":
        check_opt_alias_once(exe, inp, rundir)
    elif mode == "optnative":
        check_opt_native_exit(exe, inp, rundir)
    elif mode == "ohesssuite":
        check_ohess_suite_once(exe, inp, rundir)
    elif mode == "ohesslevelaxes":
        check_ohess_level_axes(exe, inp, rundir)
    elif mode == "ohessnative":
        check_ohess_native_hess(exe, inp, rundir)
    elif mode == "ohessbridge":
        check_ohess_fres_bridge(exe, inp, rundir)
    elif mode == "optguard":
        if len(rest) < 1:
            fail("the optguard mode needs the guard case "
                 "(scanblock, pbc or xcontrol)")
        check_opt_guard(exe, inp, rundir, rest[0])
    elif mode == "optxtbopt":
        check_opt_xtbopt_once(exe, inp, rundir)
    elif mode == "filecontractopt":
        check_file_contract_opt(exe, inp, rundir)
    elif mode == "filecontractoptfail":
        check_file_contract_opt_fail(exe, inp, rundir)
    elif mode == "filecontractvib":
        check_file_contract_vib(exe, inp, rundir)
    elif mode == "vibcomplianceskip":
        check_vib_compliance_skip(exe, inp, rundir)
    elif mode == "filecontractns":
        check_file_contract_namespace(exe, inp, rundir)
    elif mode == "ohessfail":
        check_ohess_fail(exe, inp, rundir)
    elif mode == "nativetypes":
        check_native_runtypes_unaffected(exe, inp, rundir)
    elif mode == "optlevel":
        check_opt_level_bridge(exe, inp, rundir)
    elif mode == "cycles":
        check_cycles_mapping(exe, inp, rundir)
    elif mode == "levelaxes":
        check_level_axes(exe, inp, rundir)
    elif mode == "blasthreads":
        check_blas_threads(exe, inp, rundir)
    elif mode == "blastiming":
        check_blas_threads_timing(exe, inp, rundir)
    elif mode == "banner":
        check_banner(exe, inp, rundir)
    elif mode == "programname":
        check_programname(exe, inp, rundir)
    elif mode == "nonfinite":
        check_nonfinite(exe, inp, rundir)
    elif mode == "irfatal":
        check_ir_fatal(exe, inp, rundir)
    elif mode == "helptopic":
        if len(rest) < 1:
            fail("the helptopic mode needs the path of src/prog/main.F90")
        check_help_topic(exe, inp, rundir, rest[0])
    elif mode == "helpcolour":
        check_help_colour(exe, inp, rundir)
    elif mode == "error":
        if len(rest) < 2:
            fail("the error mode needs a pattern and the xtb arguments")
        check_error(exe, inp, rundir, rest[0], rest[1:])
    else:
        fail("unknown mode %r" % mode)


if __name__ == "__main__":
    main()
