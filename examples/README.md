# Worked examples for the o-xtb landscape task suite

Sixteen runnable examples of the potential-energy-surface task suite (`--sp`, `--lopt`, `--lscan`,
`--hess`/`--lvib`, `--lts`, `--lirc`), two per task (seven tasks, fourteen cases), plus two
multi-task chains.  Every example is a complete little case: the input geometry, the exact command
line, the captured output, the natural output files, and a README that says what it shows and how to
read the result.

Every system is a **30-33 atom GFN2-xTB** job that finishes in seconds to a couple of minutes at
`OMP_NUM_THREADS=1`.  Every number quoted in a case README is a line in that case's captured log -
the logs are the record, and they ship with the case.

The examples were produced with the installed release:

```
OxTB version 0.1.0   (xTB engine 6.7.1)   Compiled on 2026/10/16 from 8a0ad67
```

The current release is **0.3.0**.  Since the examples were produced the version string has moved
three times (`0.1.0` -> `0.2.0` -> `0.2.1` -> `0.3.0`), which changes only the `OxTB version` line of
the banner and the `oxtb_version` field of the `--ljson` provenance block.  The changes under `src/`
in the range `8a0ad67..HEAD` are one guard in `src/disp/dftd4.F90` (`zeta`/`dzeta`: `qmod.lt.0` ->
`qmod.le.0`), which fires only when a Mulliken charge is **exactly** zero, and the help/CLI plumbing
of the machine-readable `--help-values` value sets (`src/xhelp.f90`, `src/prog/main.F90`), which does
not touch the energy, gradient or Hessian code.  The guard was measured not to move the record: the
relaxed distance scan `lscan-01` - the case that reaches that boundary geometry - reproduces its
shipped `xtbscan.log` (1225 lines, 35 per-point energies, `diff`-clean after normalising the build
hash in the header), the `sp-01` total energy is `-52.373771244067 Eh` as shipped, and the four
pinned water single points are unchanged.  So the recorded numbers still stand.  The features added
since the examples were produced - the man-page render pipeline, the extended coverage gate, the
`--help-values` value sets and the shell completions - are command-line surface; the examples were
**not** re-run for them.

## How to run a case

```bash
cd examples/lopt/lopt-01-lbfgs-from-distorted
export XTBPATH=/path/to/o-xtb/source          # or use an installed tree; --gfn 2 needs no file
export OMP_NUM_THREADS=1                       # the examples were measured single-threaded
bash run.sh                                    # the exact command the README documents
```

`run.sh` prints the command it runs and writes the same output file names you find in the folder, so
you can diff your run against the shipped one.  The cases are deterministic: every case was re-run
twice from clean directories and every data file (`xtbopt.xyz`, `xtbscan.log`, `hessian`,
`vibspectrum`, `path.xyz`, `irc_*.xyz`, `landscape_<task>.json`, ...) came out byte for byte
identical.  Only the captured `*.log` transcripts (run date, random seed, wall-time table) and the
`xtbtopo.mol` header line (its date stamp) differ between two runs.

## Layout

```
examples/
  README.md                 this index
  structures/               the shared 30-36 atom systems and their reaction coordinates
  sp/    sp-01-... sp-02-...       single points and the machine-readable document
  lopt/  lopt-01-... lopt-02-...   minimisation, two engines
  lscan/ lscan-01-... lscan-02-... relaxed scans, 1D and 2D
  hess/  hess-01-... hess-02-...   numerical Hessian, and the model Hessian against it
  ts/    ts-01-... ts-02-...       saddle searches: PRFO/B-PRFO/dimer, and the band family
  vib/   vib-01-... vib-02-...     thermochemistry of a minimum, imaginary mode of a saddle
  irc/   irc-01-... irc-02-...     four IRC integrations, and the direction convention
  chains/ chain-01-... chain-02-...  multi-task workflows on one structure
```

Each case folder holds `input.xyz` (plus any extra inputs), `run.sh`, one or more captured logs,
the output files with their natural names (`xtbopt.xyz`, `xtbscan.log`, `hessian`, `vibspectrum`,
`path.xyz`, `irc_forward.xyz`, `landscape_<task>.json`, ...), and `README.md`.

The category `.gitignore` files re-include those output names, which the repository-wide
`.gitignore` excludes elsewhere - that is what makes the examples self-contained.

## The cases

| category | case | system (atoms) | what it shows |
|---|---|---|---|
| `sp` | `sp-01-gfn2-single-point-json` | benzoic acid dimer (30) | one GFN2 single point with `--ljson`: energy, gap, and the unified document |
| `sp` | `sp-02-method-comparison` | benzoic acid dimer (30) | the same geometry at `--gfn 0/1/2/ff` - and why the four numbers are not comparable |
| `lopt` | `lopt-01-lbfgs-from-distorted` | anthracene + maleic anhydride (33) | L-BFGS from a deliberately distorted geometry: 41 steps, -0.956420 Eh |
| `lopt` | `lopt-02-rfo-from-distorted` | anthracene + maleic anhydride (33) | RFO with the model Hessian on the same start: 24 steps, 5 Hessian builds |
| `lscan` | `lscan-01-distance-relaxed` | anthracene + maleic anhydride (33) | bidirectional relaxed scan of the forming C9-C28 bond: 35 frames, hysteresis 9.7e-05 Eh |
| `lscan` | `lscan-02-dihedral-2d` | benzoic acid dimer (30) | relaxed 2D scan of the two proton transfers: concerted (6.05 kcal/mol) versus stepwise (23.9) |
| `hess` | `hess-01-numerical` | benzoic acid dimer (30) | numerical Hessian of a minimum: 84 modes, 6 projected zero modes, ZPE and G(T) |
| `hess` | `hess-02-model-vs-numerical` | benzoic acid dimer (30) | Lindh-D2 model Hessian against the numerical one: 1 versus 181 energy evaluations |
| `ts` | `ts-01-prfo-bprfo-dimer` | benzoic acid dimer (30) | three saddle searches agreeing to 8e-6 Eh, and their shipped-default failures |
| `ts` | `ts-02-band-family` | anthracene + maleic anhydride (33) | NEB and CI-NEB on a working band, and a documented failure of the same family |
| `vib` | `vib-01-minimum-thermochemistry` | benzoic acid dimer (30) | frequencies, ZPE, H(T) and G(T) of a minimum, and the temperature dependence |
| `vib` | `vib-02-imaginary-mode-at-the-ts` | benzoic acid dimer (30) | the single imaginary mode of an optimised saddle and the non-minimum flags |
| `irc` | `irc-01-four-methods` | benzoic acid dimer (30) | one saddle, four integrations: `hpc`, `lqa`, `eulerpc`, `gs` |
| `irc` | `irc-02-direction-and-swap` | benzoic acid dimer (30) | `--lirc-oneway` and `--lirc-swap`: what forward and backward actually mean |

## The chains

A chain carries one structure through several tasks, so you can see which file each step consumes
from the previous one.

| chain | steps |
|---|---|
| `chain-01-lopt-scan-dimer-ts-irc` | `00-structure` -> `01-lopt` -> `02-lscan` -> `03-dimer` -> `04-ts` -> `05-irc` |
| `chain-02-band-to-irc` | `00-endpoints` -> `01-band` -> `02-refine` -> `03-irc` |

Every step is itself a complete case folder with its own `run.sh`, log and README; the chain README
ties them together and gives the energy at each step.  `chain-01` shows the dimer and PRFO routes to
the same saddle; `chain-02` shows how to take the saddle out of a band (`path_ts_list.xyz`) and feed
it to the IRC.

## The command each case runs

This index is generated from every `run.sh` under `examples/`; the case READMEs still hold the
captured logs and the numbers.  Regenerate it with `python3 examples/tools/gen_command_index.py`
(the gate `xtb/landscape-examples-index` fails when it is stale).

<!-- BEGIN GENERATED: command index -->
<!-- Generated by `examples/tools/gen_command_index.py` from every `run.sh`; do not edit by hand. -->

Every block below reproduces the `export` preamble the case READMEs show and the exact
`oxtb` invocation(s) its `run.sh` runs, continuation lines joined and shell redirects kept.
The index lists **22 cases** (14 standalone cases plus 8 chain steps) in **2 chains**.

### `sp-01-gfn2-single-point-json`

category `sp` · system `structures/benzoic_acid_dimer`

```bash
export OMP_NUM_THREADS=1
export XTBPATH=<repository root>
oxtb input.xyz --gfn 2 --sp --ljson --norestart > oxtb.log 2>&1
```

### `sp-02-method-comparison`

category `sp` · system `structures/benzoic_acid_dimer`

```bash
export OMP_NUM_THREADS=1
export XTBPATH=<repository root>
oxtb input.xyz --gfn 0  --sp --norestart > oxtb-gfn0.log  2>&1
oxtb input.xyz --gfn 1  --sp --norestart > oxtb-gfn1.log  2>&1
oxtb input.xyz --gfn 2  --sp --norestart > oxtb-gfn2.log  2>&1
oxtb input.xyz --gfn ff --sp --norestart > oxtb-gfn-ff.log 2>&1
```

### `lopt-01-lbfgs-from-distorted`

category `lopt` · system `structures/anthracene_maleic_anhydride`

```bash
export OMP_NUM_THREADS=1
export XTBPATH=<repository root>
oxtb input.xyz --gfn 2 --lopt --lopt-engine lbfgs --lopt-level normal  --lopt-max-steps 300 --norestart > oxtb.log 2>&1
```

### `lopt-02-rfo-from-distorted`

category `lopt` · system `structures/anthracene_maleic_anhydride`

```bash
export OMP_NUM_THREADS=1
export XTBPATH=<repository root>
oxtb input.xyz --gfn 2 --lopt --lopt-engine rfo --lopt-hessian model  --lopt-level normal --lopt-max-steps 300 --norestart > oxtb.log 2>&1
```

### `lscan-01-distance-relaxed`

category `lscan` · system `structures/anthracene_maleic_anhydride`

```bash
export OMP_NUM_THREADS=1
export XTBPATH=<repository root>
oxtb input.xyz --gfn 2 --lscan --lscan-coord distance --lscan-atoms 7,28  --lscan-step 0.1 --lscan-nstep 17 --lscan-level normal  --lscan-max-steps 150 --norestart > oxtb.log 2>&1
```

### `lscan-02-dihedral-2d`

category `lscan` · system `structures/benzoic_acid_dimer`

```bash
export OMP_NUM_THREADS=1
export XTBPATH=<repository root>
oxtb input.xyz --gfn 2 --lscan --lscan-ndim 2  --lscan-coord  distance --lscan-atoms  9,10  --lscan-step  0.125 --lscan-nstep  8  --lscan-coord2 distance --lscan-atoms2 24,25 --lscan-step2 0.125 --lscan-nstep2 8  --lscan-level normal --lscan-max-steps 150 --norestart > oxtb.log 2>&1
```

### `hess-01-numerical`

category `hess` · system `structures/benzoic_acid_dimer`

```bash
export OMP_NUM_THREADS=1
export XTBPATH=<repository root>
oxtb input.xyz --gfn 2 --lvib --norestart > oxtb.log 2>&1
```

### `hess-02-model-vs-numerical`

category `hess` · system `structures/benzoic_acid_dimer`

```bash
export OMP_NUM_THREADS=1
export XTBPATH=<repository root>
oxtb input.xyz --gfn 2 --lvib --lvib-hessian model --norestart > oxtb-model.log 2>&1
oxtb input.xyz --gfn 2 --lvib --lvib-hessian numerical --norestart > oxtb-numerical.log 2>&1
```

### `ts-01-prfo-bprfo-dimer`

category `ts` · system `structures/benzoic_acid_dimer`

```bash
export OMP_NUM_THREADS=1
export XTBPATH=<repository root>
oxtb input.xyz --gfn 2 --lts --lts-method prfo --lts-level normal --lts-max-steps 60  --namespace prfo --norestart > prfo.log 2>&1
oxtb input.xyz --gfn 2 --lts --lts-method bprfo --lts-level normal --lts-max-steps 30  --namespace bprfo-default --norestart > bprfo-default.log 2>&1
oxtb input-displaced.xyz --gfn 2 --lts --lts-method bprfo --lts-level normal --lts-max-steps 150  --namespace bprfo --norestart > bprfo.log 2>&1
oxtb input.xyz --gfn 2 --lts --lts-method dimer --lts-level normal --lts-max-steps 120  --namespace dimer-default --norestart > dimer-default.log 2>&1
oxtb input-displaced.xyz --gfn 2 --lts --lts-method dimer --lts-level normal --lts-max-steps 150  --lts-delta 0.1 --lts-ninit force --lts-rot-max-iter 30  --namespace dimer --norestart > dimer.log 2>&1
```

### `ts-02-band-family`

category `ts` · system `structures/anthracene_maleic_anhydride`

```bash
export OMP_NUM_THREADS=1
export XTBPATH=<repository root>
oxtb product.xyz --gfn 2 --lts --lts-method neb --lts-product ts_guess.xyz --lts-images 9  --lts-level normal --lts-max-steps 200 --namespace local-neb --norestart > neb-local.log 2>&1
oxtb product.xyz --gfn 2 --lts --lts-method cineb --lts-product ts_guess.xyz --lts-images 9  --lts-level normal --lts-max-steps 200 --namespace local-cineb --norestart > cineb-local.log 2>&1
oxtb product.xyz --gfn 2 --lts --lts-method string --lts-product ts_guess.xyz  --lts-level normal --lts-max-steps 300 --namespace local-string --norestart > string-local.log 2>&1
oxtb reactant.xyz --gfn 2 --lts --lts-method neb --lts-product product.xyz --lts-images 9  --lts-level normal --lts-max-steps 100 --namespace full-neb --norestart > neb-full.log 2>&1
oxtb precomplex.xyz --gfn 2 --lts --lts-method neb --lts-product product.xyz --lts-images 11  --lts-level normal --lts-max-steps 300 --namespace precomplex-neb --norestart > neb-precomplex.log 2>&1
oxtb precomplex.xyz --gfn 2 --lopt --lopt-engine lbfgs --lopt-level loose --lopt-max-steps 50  --namespace precomplex-opt --norestart > precomplex-opt.log 2>&1
oxtb reactant.xyz --gfn 2 --lts --lts-method neb --lts-product product.xyz --lts-images 15  --lts-level normal --lts-max-steps 300 --namespace full-neb15 --norestart > neb-full15.log 2>&1
oxtb reactant.xyz --gfn 2 --lts --lts-method cineb --lts-product product.xyz --lts-images 15  --lts-level normal --lts-max-steps 300 --namespace full-cineb15 --norestart > cineb-full15.log 2>&1
oxtb reactant.xyz --gfn 2 --lts --lts-method neb --lts-product product.xyz --lts-images 15  --lts-no-align --lts-level normal --lts-max-steps 300 --namespace full-neb15-na --norestart > neb-full15-noalign.log 2>&1
oxtb reactant.xyz --gfn 2 --lts --lts-method gsm --lts-product product.xyz  --lts-level normal --lts-max-steps 300 --namespace full-gsm --norestart > gsm-full.log 2>&1
oxtb reactant.xyz --gfn 2 --lts --lts-method autoneb --lts-product product.xyz  --lts-level normal --lts-max-steps 300 --namespace full-autoneb --norestart > autoneb-full.log 2>&1
```

### `vib-01-minimum-thermochemistry`

category `vib` · system `structures/benzoic_acid_dimer`

```bash
export OMP_NUM_THREADS=1
export XTBPATH=<repository root>
oxtb input.xyz --gfn 2 --lvib --lvib-temp 298.15 --norestart > oxtb.log 2>&1
oxtb input.xyz --gfn 2 --lvib --lvib-temp 1000 --namespace t1000 --norestart > vib1000.log 2>&1
```

### `vib-02-imaginary-mode-at-the-ts`

category `vib` · system `structures/benzoic_acid_dimer`

```bash
export OMP_NUM_THREADS=1
export XTBPATH=<repository root>
oxtb input.xyz --gfn 2 --lvib --lvib-temp 298.15 --norestart > oxtb.log 2>&1
```

### `irc-01-four-methods`

category `irc` · system `structures/benzoic_acid_dimer`

```bash
export OMP_NUM_THREADS=1
export XTBPATH=<repository root>
oxtb input.xyz --gfn 2 --lirc --lirc-method hpc --lirc-max-steps 40 --norestart > oxtb.log 2>&1
oxtb input.xyz --gfn 2 --lirc --lirc-method "$m" --lirc-max-steps 40  --namespace "$m" --norestart > "$m.log" 2>&1
```

### `irc-02-direction-and-swap`

category `irc` · system `structures/benzoic_acid_dimer`

```bash
export OMP_NUM_THREADS=1
export XTBPATH=<repository root>
oxtb input.xyz --gfn 2 --lirc --lirc-method hpc --lirc-max-steps 40 --norestart > oxtb.log 2>&1
oxtb input.xyz --gfn 2 --lirc --lirc-method hpc --lirc-oneway --lirc-max-steps 40  --namespace oneway --norestart > oneway.log 2>&1
oxtb input.xyz --gfn 2 --lirc --lirc-method hpc --lirc-swap --lirc-max-steps 40  --namespace swap --norestart > swap.log 2>&1
```

### Chains

The two `chains/` directories carry one structure through several tasks; every numbered
step is a separate case folder with its own `run.sh`, listed here in the order the chain
runs them.

#### `chain-01-lopt-scan-dimer-ts-irc`

system `structures/benzoic_acid_dimer` · step order `01-lopt` → `02-lscan` → `03-dimer` → `04-ts` → `05-irc`

**`01-lopt`** (category `chains`)

```bash
export OMP_NUM_THREADS=1
export XTBPATH=<repository root>
oxtb input.xyz --gfn 2 --lopt --lopt-engine lbfgs --lopt-level normal  --lopt-max-steps 300 --norestart > oxtb.log 2>&1
```

**`02-lscan`** (category `chains`)

```bash
export OMP_NUM_THREADS=1
export XTBPATH=<repository root>
oxtb input.xyz --gfn 2 --lscan --lscan-coord distance --lscan-atoms 9,10  --lscan-step 0.05 --lscan-nstep 25 --lscan-level normal --lscan-max-steps 150  --norestart > oxtb.log 2>&1
```

**`03-dimer`** (category `chains`)

```bash
export OMP_NUM_THREADS=1
export XTBPATH=<repository root>
oxtb input.xyz --gfn 2 --lts --lts-method dimer --lts-level normal --lts-max-steps 150  --lts-delta 0.1 --lts-ninit force --lts-rot-max-iter 30 --norestart > oxtb.log 2>&1
```

**`04-ts`** (category `chains`)

```bash
export OMP_NUM_THREADS=1
export XTBPATH=<repository root>
oxtb input.xyz --gfn 2 --lts --lts-method prfo --lts-level normal --lts-max-steps 60  --norestart > oxtb.log 2>&1
```

**`05-irc`** (category `chains`)

```bash
export OMP_NUM_THREADS=1
export XTBPATH=<repository root>
oxtb input.xyz --gfn 2 --lirc --lirc-method hpc --lirc-max-steps 60 --lirc-grad-tol 1e-5  --norestart > oxtb.log 2>&1
```

#### `chain-02-band-to-irc`

system `structures/anthracene_maleic_anhydride` · step order `01-band` → `02-refine` → `03-irc`

**`01-band`** (category `chains`)

```bash
export OMP_NUM_THREADS=1
export XTBPATH=<repository root>
oxtb input.xyz --gfn 2 --lts --lts-method cineb --lts-product ../00-endpoints/ts_guess.xyz  --lts-images 9 --lts-level normal --lts-max-steps 200 --norestart > oxtb.log 2>&1
```

**`02-refine`** (category `chains`)

```bash
export OMP_NUM_THREADS=1
export XTBPATH=<repository root>
oxtb input.xyz --gfn 2 --lts --lts-method prfo --lts-level normal --lts-max-steps 80  --norestart > oxtb.log 2>&1
```

**`03-irc`** (category `chains`)

```bash
export OMP_NUM_THREADS=1
export XTBPATH=<repository root>
oxtb input.xyz --gfn 2 --lirc --lirc-method hpc --lirc-max-steps 80  --lirc-step 0.4 --lirc-grad-tol 1e-4 --norestart > oxtb.log 2>&1
```
<!-- END GENERATED: command index -->

## The shared structures

`structures/README.md` documents each system: atom count, chemistry, the reaction coordinate with
its atom indices, the validated command lines, and which tasks it suits.  Read it before writing a
new case - the atom numbering there is the numbering every case README uses.

| system | atoms | reaction |
|---|---|---|
| `anthracene_maleic_anhydride` | 33 | Diels-Alder across the meso carbons; **distinct** endpoints, the only band-capable system |
| `benzoic_acid_dimer` | 30 | concerted double proton transfer; **degenerate** (product is the same species) |
| `ptoluic_acid_dimer` | 36 | the same reaction on a larger monomer; also degenerate |

## Caveats worth knowing before you copy a recipe

* **A band needs endpoints that differ by a smooth motion.**  `ts-02` documents what happens when
  they do not: the separated fragments sit 12.6 A apart and coplanar while the product needs an
  out-of-plane rotation, so the Cartesian/IDPP interpolation walks the images through the ring and
  the SCC blows up (the pathological "highest image" barriers are 818 / 501 / 918 kcal/mol in the
  shipped logs - a non-converged, chaotic number that moves by hundreds of kcal/mol between builds).
  The same case shows the fix: start from a pre-reactive complex.
* **Some saddle searches need tuning on these systems.**  `ts-01` reports, with logs, that B-PRFO
  and the dimer method reach a higher-order saddle from the shipped `ts_guess.xyz` with the shipped
  defaults, and what it took to reach the first-order saddle instead.
* **The model Hessian is a cheap approximation, not a replacement** - `hess-02` quantifies both the
  cost (1 versus 181 energy evaluations) and the damage (an O-H frequency off by ~1000 cm-1).
* **The acid dimers are degenerate reactions.**  That is ideal for a saddle search and for showing
  that the IRC directions are what the flags say, and useless for a band: a band between identical
  endpoints has nothing to interpolate.
* **The logs ship with the cases on purpose.**  If a README claims a number, that number is a line
  in the log next to it; if you edit a case, re-run it and update both.
