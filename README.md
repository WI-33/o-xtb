# o-xtb

`o-xtb` is a fork of the semiempirical extended tight-binding program package
[`xtb`](https://github.com/grimme-lab/xtb) (version 6.7.1), extended with a
potential-energy-surface task suite.  The electronic-structure methods, the
parameter files and the original implementation come from that upstream
project; this repository adds the `sp`, `opt`, `scan`, `ts`, `freq` and `irc`
tasks on top of the same engine.  The command the build installs is `oxtb`.

## Table of Contents
- [What it is](#what-it-is)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Command line](#command-line)
- [Machine-readable output](#machine-readable-output)
- [Documentation](#documentation)
- [Contributing](#contributing)
- [Upstream](#upstream)
- [Citations](#citations)
- [License](#license)

## What it is

`o-xtb` adds a complete potential-energy-surface (PES) task suite to the xTB
engine.  The suite covers six tasks:

| Task | Run type | Covers |
|------|----------|--------|
| Single point | `--sp` | energy and gradient at the input geometry |
| Optimisation | `--opt`, `--lopt` | minimisation (sd, cg, lbfgs, bfgs, diis, rfo) |
| Scan | `--lscan` | relaxed or rigid scan of one or two internal coordinates |
| Transition state | `--lts` | PRFO, B-PRFO, dimer, NEB, CI-NEB, AutoNEB, string, GSM |
| Frequencies | `--lvib` | numerical Hessian, normal modes and thermochemistry |
| IRC | `--lirc` | GS, HPC, LQA, EulerPC |

The suite exists because the optimiser that xTB ships is task-incomplete: it
only minimises.  Upstream xTB, and therefore the engine in this tree, has no
transition-state optimiser, no NEB or dimer search and no intrinsic reaction
coordinate, so nothing above a minimum is reachable with the native engines
alone (`--method ancopt` is the legacy minimiser, still available).  Every
suite task is built on the same xTB calculator (`class(TCalculator)`) and
reuses its energy, gradient and numerical Hessian, so the underlying GFN
methods and parameters are unchanged.

Product version `0.3.2`, engine xTB `6.7.1`.

Names: the command and the project are `oxtb` / `o-xtb`, but output file names
keep their historical xTB spellings (`xtbopt.xyz`, `xtbscan.log`,
`xtbout.json`, `xtb.trj`).  That is deliberate.

## Installation

### One command: the static release

`install.sh` configures, builds, tests, installs, registers and self-checks a
fully static executable in one step, and is idempotent:

```bash
./install.sh
```

It stops at the first failing step and leaves that step's exit status behind.
The defaults can be overridden with environment variables: `XTB_PREFIX` (the
install prefix, default `$HOME/programs/o-xtb`), `XTB_BIN_DIR` (the
registration directory, default `$HOME/.local/bin`), `XTB_CMAKE` and
`XTB_JOBS`.

The static release builds with `WITH_OBJECT=OFF`, so it deliberately ships no
C API: a fully static executable cannot be linked against the shared library
`liboxtb.so`.  Build with `WITH_OBJECT=ON` (a dynamic build) if you need
`liboxtb.so` and `xtb.h`.

### Manual CMake build

```bash
cmake -S . -B build \
      -DWITH_TBLITE=OFF -DWITH_CPCMX=OFF -DWITH_TESTS=ON \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_Fortran_FLAGS_RELEASE="-O2 -DNDEBUG"
cmake --build build -j
ctest --test-dir build -R landscape
```

`WITH_TBLITE=OFF` and `WITH_CPCMX=OFF` are required: with CMake 4.3.2 the
`toml-f` dependency fails to build otherwise.  The executable is
`build/oxtb`.

Compilers and dependencies:

- GNU Fortran (`gfortran`) is the supported compiler.  `-O2` is the supported
  optimisation level and this checkout is not supported at `-O3` (CMake's
  `Release` default is `-O3` for GNU), so `CMAKE_Fortran_FLAGS_RELEASE` pins
  `-O2` explicitly.
- CMake 3.17 or newer (`cmake_minimum_required` in `CMakeLists.txt`) and a
  BLAS/LAPACK implementation are required.

### Windows (MSYS2)

The Windows build is native and lives in the **MSYS2 UCRT64** environment;
MSYS2 documents MINGW64 as deprecated, so UCRT64 is the environment to use.
Install MSYS2, open the *MSYS2 UCRT64* shell and install the toolchain:

```bash
pacman -S mingw-w64-ucrt-x86_64-gcc-fortran mingw-w64-ucrt-x86_64-openblas \
          mingw-w64-ucrt-x86_64-cmake mingw-w64-ucrt-x86_64-ninja \
          mingw-w64-ucrt-x86_64-git mingw-w64-ucrt-x86_64-python
```

OpenBLAS is LAPACK complete, so the separate `lapack` package is not needed -
a second implementation would only make the math library lookup depend on
search order.  `python` is in the line because `WITH_TESTS` defaults to `ON`
and the test target needs a Python 3 interpreter; `git` is what fetches
`mctc-lib` and `test-drive`, neither of which MSYS2 packages.

Then configure and build from the repository root:

```bash
cmake -S . -B build -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_Fortran_FLAGS_RELEASE="-O2 -DNDEBUG" \
      -DWITH_TBLITE=OFF -DWITH_CPCMX=OFF -DWITH_OBJECT=OFF
cmake --build build -j
```

The executable is `build/oxtb.exe`.  Do not pass the Windows specifics by
hand: the `if(WIN32)` blocks of `CMakeLists.txt` already supply `_WIN32`
(gfortran's preprocessor defines no target macro at all), the static OpenMP
archives, `-static` and the larger stack.  This is the flow the Windows CI job
runs on `windows-latest`; that job also passes `-DBLA_VENDOR=OpenBLAS` and
`-DWITH_TESTS=ON -DTEST_DRIVE_BUILD_TESTING=OFF` to build and run a bounded
test subset.

What a Windows build does differently:

- `WITH_OBJECT=OFF` is required: a `-static` executable cannot host the shared
  library, so there is no C API and no `liboxtb.dll` (upstream's Windows build
  has the same limitation).
- No parameter files are needed: GFN1-xTB, GFN2-xTB, GFN-FF and, since 0.3.0,
  GFN0-xTB are compiled into the binary.
- The console: the program sets its own console output code page to UTF-8 at
  startup, so a Chinese-locale console shows `·`, `²`, `Δ` and `⁻¹` correctly
  without `chcp`; redirected output is not translated and keeps the raw UTF-8
  bytes, and the console window stays on UTF-8 after the run.
- `--blas-threads` is inert in a static build and says so; set
  `OPENBLAS_NUM_THREADS` instead.
- The man pages and the shell completions are Unix install steps and do not
  apply.
- Status: this is the path CI builds on `windows-latest` with MSYS2/UCRT64, and
  a cross-built `oxtb.exe` was smoke-tested on a real Windows machine (the four
  `--gfn` methods, `--hess` and `--lopt`).  Nothing beyond that is claimed.

## Quick start

All examples use water; save this as `water.xyz`:

```
3
water
O   0.00000000000000   0.00000000000000   0.00000000000000
H   0.75860200000000   0.00000000000000  -0.50428400000000
H  -0.75860200000000   0.00000000000000  -0.50428400000000
```

A single point (`--sp`, GFN2-xTB):

```bash
oxtb water.xyz --sp --gfn 2
```

```
         :: total energy              -5.065772968305 Eh    ::
```

A geometry optimisation (`--opt`, which runs the suite L-BFGS minimiser):

```bash
oxtb water.xyz --opt --gfn 2
```

```
 landscape: lbfgs minimisation
 ...
 landscape: converged
 steps:          5
 final energy:        -5.070544444843 Eh
 final gradient:    6.984305E-05 Eh/bohr
 ...
 wrote xtbopt.xyz
```

A suite task, a relaxed scan of the O-H1 distance (`--lscan`):

```bash
oxtb water.xyz --lscan --gfn 2 \
   --lscan-coord distance --lscan-atoms 1,2 --lscan-step 0.1 --lscan-nstep 4
```

```
    point      value/Ang reached/Ang drift/Ang         energy/Eh      residual      steps  conv
        1       0.91092234       0.92781628       0.01689394      -5.069638201154   5.9247E-05      8 T
        2       1.01092234       0.99578699      -0.01513535      -5.069515287587   1.9459E-04      9 T
        3       1.11092234       1.07514461      -0.03577772      -5.061998815351   4.2842E-05      8 T
        4       1.21092234       1.16343009      -0.04749224      -5.048656668183   1.7950E-04      8 T
        5       1.31092234       1.25793646      -0.05298587      -5.031593674422   1.4239E-04      8 T
 ...
 wrote xtbscan.log with 9 frame(s), 0 failed point(s) omitted
```

The scan sweeps forward and back and reports the hysteresis between the two
branches; the five rows above are the forward half of a sweep that finishes
with four points back down the coordinate and writes nine frames to
`xtbscan.log`.  Energies are in hartree, the coordinate in ångström.

## Command line

`oxtb --help` prints a grouped index of every option in 17 categories
(`runtypes`, `suite`, `levels`, `minimisation`, `scan`, `saddle`, `band`,
`irc`, `frequencies`, `hamiltonian`, `electronic`, `external`, `solvation`,
`analysis`, `output`, `workflow`, `performance`).  The full text of one option
and one whole category are reachable too:

```bash
oxtb --help              # the grouped index
oxtb --help --method     # full text of one option
oxtb --help suite        # one category
```

An unknown topic exits 1 with an error that names the category keys.  Colour
is automatic: a terminal gets the coloured page, a pipe or a redirect gets the
plain page with no escape bytes.  `CLICOLOR_FORCE=1` forces colour even into a
pipe and wins over everything else; `NO_COLOR=1` disables it.

`--method` is the canonical selector for the algorithm of the current task.
It is the single name of `--lopt-engine`, `--lscan-method`, `--lts-method`
and `--lirc-method`, is not case-sensitive, and a value that does not fit the
run type aborts the run.

## Machine-readable output

Pass `--ljson` to write one `landscape_<TASK>.json` document per run, for
every task, whether it is a native run type or a suite task: `--sp` writes
`landscape_sp.json`, `--lvib` writes `landscape_vib.json`.  Every document
carries `format: 2` and the same top-level core keys in the same order:

```
format, task, status, status_name, converged, message, n_atoms, energy,
geometry (unit = bohr, symbols, coordinates), provenance,
counters, values, flags, history
```

Task-specific data stays inside `counters`, `values`, `flags` and `history`,
so a single parser reads every task.  The older `--json` option is unchanged
and still writes the legacy `xtbout.json` (which the task suite does not
touch).

## Documentation

- `man 1 oxtb` and `man 7 xcontrol` — the manual pages.  The AsciiDoc sources
  `man/oxtb.1.adoc` / `man/xcontrol.7.adoc` are installed under `share/man`, and
  `man/tools/adoc2man.py` renders them into the roff pages `man/oxtb.1` /
  `man/xcontrol.7`, which are committed and installed alongside the sources
  (the drift gate `xtb/landscape-man-render` keeps the two in step).

The documentation site is published at <https://owenchan-0216.github.io/o-xtb/>
and is rebuilt on every push to `main`.  It holds the user guide, the generated
command-line reference and a method-notes chapter; its sources live in
`docs/docs-site/`, and both `mkdocs build --strict` and the dependency-free site
gate (`docs/docs-site/tools/check_site.py`) run in CI.

## Contributing

Please read our [contributing guidelines](CONTRIBUTING.md)
before contributing to this project.

## Upstream

`o-xtb` is a fork of the [`xtb`](https://github.com/grimme-lab/xtb) program
package by the Grimme group at the University of Bonn.  The scientific methods,
the parameter files and the original implementation come from that project; the
history of this repository includes its commits.

The contributors of the upstream project are listed on the upstream repository
(<https://github.com/grimme-lab/xtb/graphs/contributors>) and in its own README.
This repository does not restate them, because it is not the upstream project.

## Citations

General Reference to `xtb` and the implemented GFN methods:
- C. Bannwarth, E. Caldeweyher, S. Ehlert, A. Hansen, P. Pracht, J. Seibert, S. Spicher, S. Grimme
  *WIREs Comput. Mol. Sci.*, **2020**, 11, e01493.
  DOI: [10.1002/wcms.1493](https://doi.org/10.1002/wcms.1493)

for GFN-xTB:
- S. Grimme, C. Bannwarth, P. Shushkov, *J. Chem. Theory Comput.*, **2017**, 13, 1989-2009.
  DOI: [10.1021/acs.jctc.7b00118](https://dx.doi.org/10.1021/acs.jctc.7b00118)
- C. Bannwarth, S. Ehlert and S. Grimme., *J. Chem. Theory Comput.*, **2019**, 15, 1652-1671.
  DOI: [10.1021/acs.jctc.8b01176](https://dx.doi.org/10.1021/acs.jctc.8b01176)
- P. Pracht, E. Caldeweyher, S. Ehlert, S. Grimme, *ChemRxiv*, **2019**, preprint.
  DOI: [10.26434/chemrxiv.8326202.v1](https://dx.doi.org/10.26434/chemrxiv.8326202.v1)

for GFN-FF:
- S. Spicher and S. Grimme, *Angew. Chem. Int. Ed.*, **2020**, 59, 15665–15673
  DOI: [10.1002/anie.202004239](https://doi.org/10.1002/anie.202004239)

for PTB:
- S. Grimme, M. Müller, A. Hansen, *J. Chem. Phys.* **2023**, 158, 124111.
  DOI: [10.1063/5.0137838](https://doi.org/10.1063/5.0137838)

for GBSA and ALPB implicit solvation:
- S. Ehlert, M. Stahn, S. Spicher, S. Grimme,
  *J. Chem. Theory Comput.*, **2021**, 17, 4250-4261
  DOI: [10.1021/acs.jctc.1c00471](https://doi.org/10.1021/acs.jctc.1c00471)

for ddCOSMO and CPCM-X implicit solvation:
- M.Stahn, S. Ehlert, S. Grimme,
  *J. Phys. Chem. A*, **2023**, XX, XXX-XXX
  DOI: [10.1021/acs.jpca.3c04382](https://doi.org/10.1021/acs.jpca.3c04382)

for DFT-D4:
- E. Caldeweyher, C. Bannwarth and S. Grimme, *J. Chem. Phys.*, **2017**, 147, 034112.
  DOI: [10.1063/1.4993215](https://dx.doi.org/10.1063/1.4993215)
- E. Caldeweyher, S. Ehlert, A. Hansen, H. Neugebauer, S. Spicher, C. Bannwarth and S. Grimme, *J. Chem. Phys.*,
  **2019**, 150, 154122. DOI: [10.1063/1.5090222](https://dx.doi.org/10.1063/1.5090222)
- E. Caldeweyher, J.-M. Mewes, S. Ehlert and S. Grimme, *Phys. Chem. Chem. Phys.*, **2020**, 22, 8499-8512. 
  DOI: [10.1039/D0CP00502A](https://dx.doi.org/10.1039/D0CP00502A) 

for sTDA-xTB:
- S. Grimme and C. Bannwarth, *J. Chem. Phys.*, **2016**, 145, 054103.
  DOI: [10.1063/1.4959605](https://dx.doi.org/10.1063/1.4959605)

in the mass-spec context:
- V. Asgeirsson, C. Bauer and S. Grimme, *Chem. Sci.*, **2017**, 8, 4879.
  DOI: [10.1039/c7sc00601b](https://dx.doi.org/10.1039/c7sc00601b)
- J. Koopman and S. Grimme, *ACS Omega*, **2019**, 4, 12, 15120-15133.
  DOI: [10.1021/acsomega.9b02011](https://dx.doi.org/10.1021/acsomega.9b02011)
- J. Koopman and S. Grimme, *J. Am. Soc. Mass Spectrom.*, **2021**, 32, 7, 1735-1751.
  DOI: [10.1021/jasms.1c00098](https://dx.doi.org/10.1021/jasms.1c00098)

for metadynamics refer to:
- S. Grimme, *J. Chem. Theory Comput.*, **2019**, 155, 2847-2862.
  DOI: [10.1021/acs.jctc.9b00143](https://dx.doi.org/10.1021/acs.jctc.9b00143)
  
for SPH calculations refer to:
- S. Spicher and S. Grimme, *J. Chem. Theory Comput.*, **2021**, 17, 1701–1714.
  DOI: [10.1021/acs.jctc.0c01306](https://doi.org/10.1021/acs.jctc.0c01306) 

for ONIOM refer to:
- C. Plett, A. Katbashev, S. Ehlert, S. Grimme, M. Bursch, *Phys. Chem. Chem. Phys.*, **2023**, 25, 17860-17868.
  DOI: [10.1039/D3CP02178E](https://doi.org/10.1039/D3CP02178E)

for O1NumHess refer to:
- B. Wang, S. Luo, Z. Wang, W. Liu, *J. Chem. Theory Comput.* **2025**, 21, 10893-10909.
  DOI: 10.1021/acs.jctc.5c01354

All references are available in [bibtex format](./assets/references.bib).


## License

`xtb` is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

`xtb` is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the
GNU Lesser General Public License for more details.

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in `xtb` by you, as defined in the
GNU Lesser General Public license, shall be licensed as above, without any
additional terms or conditions.
