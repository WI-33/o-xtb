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
"""Precision oracle: the delivered numbers, asserted at the delivered precision.

The suite's existing gates are regression tripwires, not precision oracles: the
``--lvib``/``--hess`` frequency comparison allows 0.2 cm^-1 against a measured
0.0100, the installer's water check allows 1e-6 Eh, and three of the four
pinned water energies of ``AGENTS.md`` section 6 (GFN1 and GFNff among them) were
recorded datums with no test at all.  A patch that degraded one of those numbers
by a relative 1e-5 would have stayed green everywhere.  This check closes that
gap for two quantities that are cheap to produce and expensive to fake:

**The four water single points.**  The reference geometry is the one the GFN0
test uses (``test/unit/CMakeLists.txt``, ``LANDSCAPE-GFN0-WATER``), a staggered
water at the experimental geometry.  Each of ``--gfn 0/1/2/ff`` is run with
``--ljson`` and the document's canonical ``energy`` field (18 decimals) is
compared against the delivered value, which is that field for the release build
of commit ``2f91a0d`` - not the 12-decimal ``total energy`` line, whose own
rounding would spend part of the tolerance.  The tolerance is **1e-12 Eh**:

  * at |E| <= 5.8 Eh one unit in the last place is 8.9e-16 Eh, so 1e-12 Eh is
    ~1120 ulp - a budget, not a bit pattern;
  * with ``OMP_NUM_THREADS=1`` (set by the test registration) no OpenMP
    reduction can reassociate, so the run-to-run spread of a fixed binary is
    **zero**;
  * the measured cross-build spread is also zero - the four energies are
    bit-identical before and after the 2026-10-16 performance PRs
    (``AGENTS.md`` section 2) - so the 1120 ulp are a margin for a *different*
    BLAS or compiler, not for this build.  Measured cross-compiler: the
    gfortran 14.4.0 build of the gating CI reproduces GFN0, GFN1 and GFNff
    exactly and GFN2 to +2.7e-15 Eh, so the whole budget is available to a
    compiler difference.  A change that moves any of the four by more than
    1e-12 Eh is a real numerical change and has to be re-baselined in a report,
    which is exactly what this check exists to force.

**The finite difference Hessian of the same water.**  ``--hess`` writes the raw
matrix to ``hessian`` (a ``$hessian`` header and 9x9 values with 10 decimals).
The check asserts three things: the shape and header of the file, that two runs
of the same binary are **byte identical**, and that every entry reproduces the
matrix delivered at the commit this oracle was written against to within
**1e-9 Eh/bohr^2** (about 10x the file's own print resolution).  The stored
matrix is 81 numbers, so it needs no golden file next to the source.

Both assertions are falsifiable in the obvious way: perturb one stored energy or
one stored matrix entry and the check fails on that quantity alone.

Usage::

    landscape_precision_check.py <oxtb> <rundir> <water.xyz>

Exit status is non-zero when an energy leaves its tolerance, a matrix entry
leaves its tolerance, the file shape changed, or the two ``--hess`` runs of the
same binary disagree by one byte.
"""

import json
import subprocess
import sys
from pathlib import Path

#: The delivered water single points (hartree), in the order of ``METHODS``,
#: taken from the full-precision ``energy`` field of the ``--ljson`` document
#: (the field carries 18 decimals, not the 12 the ``total energy`` line
#: prints).  The 12-decimal values of ``AGENTS.md`` section 6 are
#: -4.361986999510 / -5.764209595203 / -5.065772968305 / -0.319465230453 Eh;
#: the extra digits below are what those lines round.  Measured on the release
#: build of commit 2f91a0d with gfortran 12.2.0, they reproduce on the gfortran
#: 14.4.0 build of the gating CI to 0 (GFN0/1/ff) and +2.7e-15 Eh (GFN2).
#: GFN0 and GFN2 were gated before (a 10-digit prefix and a 1e-6 installer
#: check); GFN1 and GFNff were datums only.
METHODS = [("0", -4.36198699950954), ("1", -5.76420959520329),
           ("2", -5.065772968305037), ("ff", -0.31946523045327296)]
#: Energy tolerance (hartree): 1e-12 Eh, ~1120 ulp at 5.8 Eh.  The reference
#: is the unrounded document value, so the whole budget is cross-build margin
#: (measured cross-compiler spread above: 2.7e-15 Eh) instead of partly the
#: rounding of a 12-decimal datum, which alone spent 5.4e-13 of it (report
#: 224, finding 2).
ENERGY_TOL = 1.0e-12

#: The raw ``--hess`` matrix of the same water at GFN2 (hartree/bohr^2), rows of
#: the 9x9 matrix, 10 decimals each - the precision the file itself carries.
HESSIAN = """
     0.9450203445  0.0000000000 -0.0000000000 -0.4725101723  0.0000000000  0.3141032052 -0.4725101723  0.0000000000 -0.3141032052
     0.0000000000  0.0000000000  0.0000000000 -0.0000000000  0.0000000000 -0.0000000000 -0.0000000000  0.0000000000  0.0000000000
    -0.0000000000  0.0000000000  0.4561871336  0.2865951422  0.0000000000 -0.2280935668 -0.2865951422  0.0000000000 -0.2280935668
    -0.4725101723 -0.0000000000  0.2865951422  0.5041375861  0.0000000000 -0.3003491737 -0.0316274138  0.0000000000  0.0137540315
     0.0000000000  0.0000000000  0.0000000000  0.0000000000  0.0000000000 -0.0000000000 -0.0000000000  0.0000000000 -0.0000000000
     0.3141032052 -0.0000000000 -0.2280935668 -0.3003491737 -0.0000000000  0.2184475237 -0.0137540315 -0.0000000000  0.0096460431
    -0.4725101723 -0.0000000000 -0.2865951422 -0.0316274138 -0.0000000000 -0.0137540315  0.5041375861 -0.0000000000  0.3003491737
     0.0000000000  0.0000000000  0.0000000000  0.0000000000  0.0000000000 -0.0000000000 -0.0000000000  0.0000000000 -0.0000000000
    -0.3141032052  0.0000000000 -0.2280935668  0.0137540315 -0.0000000000  0.0096460431  0.3003491737 -0.0000000000  0.2184475237
"""
#: Matrix tolerance (hartree/bohr^2): 10x the file's print resolution (1e-10).
HESSIAN_TOL = 1.0e-9
#: The 3x3-matrix count of the water fixture.
HESSIAN_VALUES = 81


def fail(message, output=""):
    sys.stderr.write("landscape_precision_check.py: %s\n" % message)
    if output:
        sys.stderr.write(output if output.endswith("\n") else output + "\n")
    raise SystemExit(1)


def run(binary, args, rundir):
    """Run the binary in ``rundir``; return (exit status, stdout+stderr)."""
    rundir.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run([str(binary), *args], cwd=rundir,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          text=True, timeout=1800)
    return proc.returncode, proc.stdout


def read_hessian(path):
    """The values of a ``$hessian`` file, in file order."""
    values = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("$"):
                continue
            values += [float(token) for token in line.split()]
    return values


def check_energies(binary, water, rundir):
    """The four water single points of the reference geometry."""
    for method, reference in METHODS:
        run_dir = rundir / ("sp-gfn%s" % method)
        status, output = run(binary, [str(water), "--gfn", method, "--sp",
                                      "--ljson", "--norestart"], run_dir)
        if status != 0:
            fail("the --gfn %s single point exited %d" % (method, status), output)
        document = run_dir / "landscape_sp.json"
        if not document.is_file():
            fail("the --gfn %s run wrote no %s" % (method, document.name), output)
        with open(document) as handle:
            energy = json.load(handle)["energy"]
        delta = energy - reference
        printed = "?"
        for line in output.splitlines():
            if "total energy" in line and " Eh " in line:
                tokens = line.split()
                printed = tokens[tokens.index("Eh") - 1]
                break
        print("  --gfn %-2s total energy %s Eh  delta %+.3e" % (method, printed,
                                                                delta))
        if abs(delta) > ENERGY_TOL:
            fail("the --gfn %s single point is %.17g Eh, the delivered "
                 "value is %.17g Eh: the delta %+.3e exceeds the %.0e Eh "
                 "tolerance"
                 % (method, energy, reference, delta, ENERGY_TOL), output)


def check_hessian(binary, water, rundir):
    """The raw finite difference Hessian: shape, determinism, stored matrix."""
    names = ("first", "second")
    for name in names:
        status, output = run(binary, [str(water), "--gfn", "2", "--hess",
                                      "--norestart"], rundir / name)
        if status != 0:
            fail("the %s --hess run exited %d" % (name, status), output)
    first = (rundir / "first" / "hessian").read_bytes()
    second = (rundir / "second" / "hessian").read_bytes()
    if first != second:
        fail("two --hess runs of the same binary wrote different hessian files "
             "(%d vs %d bytes)" % (len(first), len(second)))
    print("  --hess hessian file %d bytes, byte identical over two runs"
          % len(first))

    values = read_hessian(rundir / "first" / "hessian")
    if len(values) != HESSIAN_VALUES:
        fail("the hessian file carries %d values, expected %d"
             % (len(values), HESSIAN_VALUES))
    reference = [float(token) for token in HESSIAN.split()]
    if len(reference) != HESSIAN_VALUES:
        fail("the stored reference matrix carries %d values, expected %d"
             % (len(reference), HESSIAN_VALUES))
    worst_index = max(range(len(values)),
                      key=lambda i: abs(values[i] - reference[i]))
    worst = abs(values[worst_index] - reference[worst_index])
    print("  --hess matrix: largest deviation %.3e Eh/bohr^2 at element %d"
          % (worst, worst_index + 1))
    if worst > HESSIAN_TOL:
        fail("hessian element %d is %.10f, the delivered value is %.10f: the "
             "deviation %.3e exceeds the %.0e tolerance"
             % (worst_index + 1, values[worst_index], reference[worst_index],
                worst, HESSIAN_TOL))


def main(argv):
    if len(argv) != 4:
        sys.stderr.write("usage: landscape_precision_check.py <oxtb> <rundir> "
                         "<water.xyz>\n")
        return 2
    binary, rundir, water = Path(argv[1]).resolve(), Path(argv[2]), Path(argv[3]).resolve()
    if not binary.is_file():
        fail("no executable at %s" % binary)
    if not water.is_file():
        fail("no water fixture at %s" % water)
    rundir.mkdir(parents=True, exist_ok=True)

    check_energies(binary, water, rundir)
    check_hessian(binary, water, rundir)
    print("precision: four water energies within %.0e Eh, the water Hessian "
          "within %.0e Eh/bohr^2" % (ENERGY_TOL, HESSIAN_TOL))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
