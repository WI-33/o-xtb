# JSON output (`--ljson`)

`--ljson` is the **single machine-readable entry point for every task**.  The
six suite tasks and the native run types all write one
`landscape_<TASK>.json` document, so a single parser reads every task.  The
historical `landscape_` prefix is kept for the native tasks too.

```bash
oxtb water.xyz --gfn 2 --sp --ljson          # writes landscape_sp.json
oxtb water.xyz --gfn 2 --lopt --ljson        # writes landscape_opt.json
```

The schema version is **`format: 2`**.  A parser branches on `format`.
`xtbout.json` (written by the separate `--json` flag) is a different, older
format and is unchanged by this document format.

## Core keys

Every document carries the same top-level keys, always in this order, whatever
the task:

<!-- gate:core-keys -->
```text
format
task
status
status_name
converged
message
n_atoms
energy
geometry
provenance
counters
values
flags
history
```

| Key | Type | Meaning |
|---|---|---|
| `format` | integer | schema version, `2` |
| `task` | string | the task name, e.g. `sp`, `opt`, `scan`, `ts`, `vib`, `irc` |
| `status` | integer | numeric status code |
| `status_name` | string | the status code as a name |
| `converged` | bool | whether the task converged |
| `message` | string | a human-readable message (empty when there is none) |
| `n_atoms` | integer | number of atoms |
| `energy` | number | the canonical total energy in hartree (duplicates the task-specific energy in `values`) |
| `geometry` | object | the geometry block, see below |
| `provenance` | object | the run metadata, see below |
| `counters` | object | task-specific counters (energies, gradient norms, counts) |
| `values` | object | task-specific values |
| `flags` | object | task-specific boolean flags |
| `history` | object | per-step convergence series |

The task-specific data stays inside `counters`, `values`, `flags` and
`history`, so a parser can ignore what it does not know without failing.

## The `geometry` block

<!-- gate:geometry-keys -->
```text
unit
symbols
coordinates
```

`unit` is always the string `"bohr"`; `symbols` is the list of element symbols,
one per atom; `coordinates` is the list of `[x, y, z]` positions, in bohr, in
the same order as `symbols`.  The geometry block has the same shape for every
task.

## The `provenance` block

<!-- gate:provenance-keys -->
```text
oxtb_version
engine_version
commit
command_line
method
charge
spin
input_file
```

| Key | Meaning |
|---|---|
| `oxtb_version` | the product version of this fork, e.g. `0.3.2` |
| `engine_version` | the xTB engine version, e.g. `6.7.1` |
| `commit` | the git commit the binary was built from |
| `command_line` | the caller's own command line (not a property of the process) |
| `method` | the electronic-structure method, e.g. `GFN2-xTB` |
| `charge` | the total charge |
| `spin` | the number of unpaired electrons |
| `input_file` | the input file name |

## Frozen properties

* the top-level keys are always the same and always in the same order;
* counters and values are written in the order the task added them, so a diff
  of two files of the same task is stable;
* real numbers are written with 17 significant digits, enough for a bitwise
  double round trip, and the conversion is deterministic, so two runs of the
  same input at the **same `OMP_NUM_THREADS`** produce byte-identical files;
  at different thread counts they differ in the last digits, because the
  OpenMP reduction order is not fixed;
* no wall time, no path and no random quantity is written;
* a missing value is JSON `null`, never a fabricated number.

!!! note "The key lists above are checked"
    The three key lists in this page are the contract the site gate verifies:
    `docs/docs-site/tools/check_site.py` compares them with the enforced
    constants of the end-to-end checker (`test/unit/landscape_e2e_check.py`), so
    a task that renames or reorders a core key fails the build until this page
    is updated with it.
