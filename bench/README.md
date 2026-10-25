# `bench/` — wall-clock benchmark of the PES tasks

`test/` certifies *numbers*; nothing in this tree measures *time*.  Every
performance statement in `docs/` was taken with a hand-written command line on
one machine, which is why none of them can be reproduced from the repository
alone.  `bench/tools/bench_landscape.py` is that command line.

It times one fixed set of tasks on fixtures that are **in the repository**
(nothing is generated, nothing is downloaded) and prints a table of wall-clock
times.  It prints numbers, not verdicts: timings depend on the machine, the
load, the compiler and the BLAS, so **no gate runs this tool** and it has no
`--check` mode.  Compare two runs of the same machine and the same binary; never
compare a run against a stored number.

## Running it

```bash
# after a normal configure (the target exists when a Python 3 interpreter is found)
cmake --build <build-dir> --target benchmark

# or directly, which also works for an installed binary
bench/tools/bench_landscape.py <build-dir>/oxtb
bench/tools/bench_landscape.py <build-dir>/oxtb --runs 5          # more repetitions
bench/tools/bench_landscape.py <build-dir>/oxtb --only lvib,lirc-hpc
bench/tools/bench_landscape.py --list                            # the case table
```

The machine model is pinned: `OMP_NUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1`
are forced into the child environment, one process per timed run, a **fresh
rundir per repetition** and `--norestart`.  `XTBPATH` is inherited, not set — the
build tree finds its own parameter directory (see `AGENTS.md` §5, exception #1).

## Cases

| case | task | fixture | atoms |
|---|---|---|---|
| `sp-gfn2` | `--sp` | `examples/sp/sp-01-gfn2-single-point-json/input.xyz` | 30 |
| `lopt-lbfgs` | `--lopt --lopt-engine lbfgs` | `examples/lopt/lopt-01-lbfgs-from-distorted/input.xyz` | 33 |
| `lopt-rfo` | `--lopt --lopt-engine rfo` | same input | 33 |
| `lvib` | `--lvib` | `examples/hess/hess-01-numerical/input.xyz` | 30 |
| `lts-band-neb` | `--lts --lts-method neb` (5 images, 5 steps) | `examples/structures/benzoic_acid_dimer/reactant.xyz` -> `ts_guess.xyz` | 30 |
| `lirc-hpc` | `--lirc --lirc-method hpc --lirc-max-steps 40` | `examples/irc/irc-01-four-methods/input.xyz` | 30 |
| `hess-taxol` | `--hess` | `assets/inputs/xyz/taxol.xyz` | 113 |

The default repetition count is 3, except for the two cases that cost of the
order of 10 s per run (`lopt-rfo`, `hess-taxol`), which default to 1; `--runs`
overrides both.  The whole default set takes about 75 s on the reference machine.

**The `evals` column** is the landscape engine's energy-evaluation counter
(`landscape engine: energy evaluations`), which is the quantity most of the
speedup candidates in `docs/_simplify-recon/220-performance-campaign.md` move.
The native `--sp` and `--hess` routes do not print it, so their rows show `-`.

**A non-zero exit status** in the table is not a defect of the tool: the band
case stops at its 5-step budget on a real reaction path (reactant -> TS guess,
5 images), so it reports "band reached the iteration limit without converging"
and exits 1 like any other bounded `--lts` run.  It is deliberately the
propagating case: a band that converges after zero iterations (a degenerate
reactant/product pair) would only time the initial interpolation and could not
see a change of the band kernel.

## Baseline at `4f1a0ef`

Release `-O2`, `gfortran 12.2.0`, one thread, idle 8-core machine, median of the
`n` runs, recorded by the campaign's step 0 (report
`docs/_simplify-recon/221-step0-benchmark-and-s2.md`).  The band row is the one
number report 221 changed: at `4f1a0ef` a failed search still paid for a
confirmation Hessian, which the five-step budget case shows as 211 evaluations
and 1.96 s; after the S2 fix (only a converged search is confirmed) the same
command line costs 31 evaluations and 0.26 s.  The endpoint cache of the next
commit (`landscape_band_endpoint_cache`, report 223) removes the two fixed
endpoints from every sweep after the first, which the same case shows as **20
evaluations and 0.19 s** (measured on the fixed tree of report 225, gfortran
12.2.0).  Every other row is unchanged.

| case | n | median / s | min / s | max / s | evals |
|---|---|---|---|---|---|
| `sp-gfn2` | 3 | 0.021 | 0.021 | 0.022 | – |
| `lopt-lbfgs` | 3 | 0.626 | 0.618 | 0.636 | 43 |
| `lopt-rfo` | 1 | 14.804 | 14.804 | 14.804 | 1015 |
| `lvib` | 3 | 1.098 | 1.094 | 1.131 | 181 |
| `lts-band-neb` | 3 | 1.963 | 1.943 | 1.995 | 211 |
| `lirc-hpc` | 3 | 2.275 | 2.271 | 2.341 | 221 |
| `hess-taxol` | 1 | 55.112 | 55.112 | 55.112 | – |
| **total (medians)** | | **75.9** | | | |
