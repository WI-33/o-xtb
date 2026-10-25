# o-xtb

**o-xtb** is a fork of [xTB 6.7.1](https://github.com/grimme-lab/xtb) that adds a
complete **potential energy surface (PES) task suite** on top of the GFN2-xTB
semiempirical engine.  The xTB engine itself is unchanged; the suite is an
additional layer of tasks written in Fortran under `src/landscape/`.

The command is called `oxtb`; the shared library is `liboxtb.*`.  The xTB
module and symbol names, the `param_*.txt` data files, the environment variable
`XTBPATH` and the output file names (`xtbopt.xyz`, `xtbscan.log`, ...) are
deliberately **not** renamed.

## Two version numbers

| Number | Meaning |
|---|---|
| **o-xtb 0.3.2** | the product version of this fork |
| **xTB 6.7.1** | the engine version this fork tracks |

Both are printed by the startup banner:

```text
 * OxTB version 0.3.2 *
 * xTB version 6.7.1 *
```

The banner also prints the commit the binary was built from.  A run started
with `--ljson` records the same two numbers in the `provenance` block of its
JSON document (see the [JSON output](reference/ljson.md) chapter).

## Install

The supported route is the one-command static release, which configures,
builds, tests, installs and self-checks:

```bash
./install.sh
```

It installs into `~/programs/o-xtb` (executable `~/programs/o-xtb/bin/oxtb`)
and registers `~/.local/bin/oxtb` as a symlink.  See
[Getting started](getting-started.md) for the manual build and the first run.

## The task suite

The suite covers the whole PES task matrix with one task driver per family:

| Task | Command | Writes |
|---|---|---|
| Single point / gradient | `--sp`, `--grad` | scalar results; `xtbout.json` with `--json` |
| Geometry optimisation | `--opt`, `--lopt` | `xtbopt.xyz`, `xtbopt.log`, `xtb.trj` |
| Surface scan | `--lscan` | `xtbscan.log` |
| Transition state | `--lts` | `xtbts.xyz`, `xtbts.log`, `path.xyz` |
| Frequencies / thermochemistry | `--hess`, `--lvib` | `hessian`, `vibspectrum`, `g98.out` |
| Intrinsic reaction coordinate | `--lirc` | `irc_forward.xyz`, `irc_backward.xyz`, `irc_path.xyz`, `irc.log` |

Each family has a page in the [user guide](user-guide/index.md); the exact
options are in the generated [command-line reference](reference/cli.md).

## Map of this site

* **[Getting started](getting-started.md)** - build, install and a first run.
* **[User guide](user-guide/index.md)** - one page per task family.
* **[Reference](reference/cli.md)** - the command-line reference (generated
  from the program itself), the `--ljson` document format and the output files.
* **[Method notes](method-notes/index.md)** - how selected mechanisms work.
* **[API reference](api/index.md)** - the Fortran API, not published yet.

!!! note "This site is built from the repository"
    Every command, file name and number on these pages is taken from a build of
    this repository.  Pages that would need content that does not exist yet say
    so explicitly instead of guessing.
