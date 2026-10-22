<!-- The development review reports cited below (the numbered reports under docs/_landscape-recon/ and docs/_simplify-recon/) and AGENTS.md are not part of this repository. -->
# The xtb CLI redesign direction — overall delivery report

> **Reader**: a developer/reviewer who needs to decide "what was delivered, why it can be trusted, and what risk remains" without reading ~60 staged reports.
> **Evidence density over length**: every conclusion is followed by a verifiable location (commit id / `file:line` / report name / command + expected output).
> **Annotation convention**:
> - `[measured · this report]` = a command I **ran myself** while writing this report (command and output in §8 and appendix A).
> - `[measured · report N]` = a measurement quoted from a staged report; I give the report name and `file:line`.
> - `[unverified]` = I did **not** independently recompute it, I only faithfully relayed it (collected in appendix B).
> - This report **changes no code/test/AGENTS.md**, it only adds this file; `docs/` is already ignored by `.gitignore` and has been untracked from git, so there is **no `git add`** at all.

---

## 1. Summary

### 1.1 What this direction solved

The repository is a `grimme-lab/xtb` 6.7.1 base plus a complete PES task suite already landed in `src/landscape/` (M0..M7: `sp/opt/scan/ts/hessian/freq/irc`, see `AGENTS.md` §2). The suite itself works, but it **exposed two independent problems**:

1. **The native tasks are incomplete**: xtb's default optimiser ANCopt only **minimises**; there is no TS search, no NEB/dimer, no IRC. The suite implements them, but **only hangs them off the `--l*` private flags**, while the user's default path (`--opt`) still runs ANCopt.
2. **The CLI surface is chaotic**: the level semantics of `--opt`/`--ohess`, the naming split between `--lopt-engine` vs `--lscan-method`, the word `-step` with three meanings and two units, conflicting combinations silently winning, required-item errors not naming the flag — see the survey report `docs/_simplify-recon/02-cli-audit.md` (886 lines).

The output of this direction (W0 + W1 design + stages 1/2/2b/3a/3b/3c) is: **re-dispatch `--opt`/`--ohess` to the suite** (option C), add a **unified method selector `--method`**, and ground the "never silent" principle in **18 breaking changes**; at the same time **fix a batch of real defects** (including 3 W0 defects that produced silently wrong results and 1 upstream compiler defect).

### 1.2 Net effect at a glance

| Item | Value | Source |
|---|---|---|
| direction commits | **9** (plus 1 `cdd9ffb` that only moves `docs/` out of git tracking, a meta-commit) | `git log --oneline` |
| code change size (`src/`+`test/`+`cmake/`+`man/`, `24c069f..HEAD`) | 38 files, **+9018 / −483** lines | `git diff --shortstat 24c069f..HEAD -- src/ test/ cmake/ man/` `[measured · this report]` |
| of which `src/` | 25 files, +3107 / −396 | same as above |
| of which `test/` | 13 files, +5911 / −87 | same as above |
| landscape ctest target count | 95 (M7 delivery) → **147** | `ctest -N -R landscape` = 147 `[measured · this report]`; the M7 95 is in `docs/xtb-landscape-M7-validation.md` §2 |
| breaking changes | **18** (A3 + B0 + C8 + D1 + E1 + F4 + G1), all tested | design `11-cli-redesign.md` §10 summary row |
| core-layer changes (`src/type\|xtb\|mctc\|param\|io\|disp\|coulomb`) | **0** | `git diff --name-only 24c069f..HEAD -- <7 dirs>` empty `[measured · this report]` |
| `use xtb_landscape_*` outside `main.F90`/`src/landscape/` | **0** | `grep -rln "use xtb_landscape" src/ \| grep -v main.F90 \| grep -v src/landscape/` empty `[measured · this report]` |
| landscape tests at the close of this direction | **147/147 passed** (44.4 s, `OMP_NUM_THREADS=1`) | `ctest -R landscape` `[measured · this report]`; the tree carries **171** `add_test(` entries at `d7181ea` (2026-10-22), up from the 147 registered here |
| current regression subset | **12/12 pass** | `ctest -R 'Singlepoint\|GFN\|model-hessian\|optimizer\|Argparser\|CAPI'` `[measured · this report]` |
| work tree | clean, `git diff --cached` empty | `git status --short` empty `[measured · this report]` |

### 1.3 One-sentence conclusion

**Deliverable**: every registered breaking change has a falsifiable test, the final conclusion of all three independent reviews is **accept** (one round each for 3a/3b/3c, see `48/54/57` and `58`), the core layer has 0 changes, and the re-dispatch of `--opt`/`--ohess` has 4 native exits (`--method ancopt`/`lancopt`/`fire`/`pbc_lbfgs`). **But**: stage 4 (documentation close-out) is not done, the rc channel (`.xtbrc`) is still silently ignored by the suite, and the meson entries cannot be verified on this machine — see §7.

---

## 2. Delivery list (9 commits)

In chronological order (`git log --oneline` reversed goes from `fa969ce` to `bdd8eb8`). Size column = **code** (`src/`+`test/`+`cmake/`+`man/`) / **total** (including reports, while `docs/` was still tracked).

| # | Commit | Subject | Size (code / total) | What landed |
|---|---|---|---|---|
| W0 | `bdd8eb8` | Switch the build to `-O2` and stop three silent wrong results | 1458+ / 73− ; total 5471+ / 75− | Three silently wrong results (`$fix`, `$thermo`, `$scan`) + build switched to `-O2`; added `types/constraints.f90`, `test/unit/landscape_cli_check.py`, `landscape_testpot.f90` |
| W1 | `a845def` | Add the survey and the CLI redesign | **0 code** / total 4602+ | 4-way capability-boundary survey (`08`) + **design document `11-cli-redesign.md`** (2851 lines of body) + 4 rounds of review records |
| 1 | `f5b070f` | conflict check, flag names in errors, echo risky defaults | 2256+ / 229− | conflict detection (C1/C2/C8), naming required items, echoing risky defaults, two-way `--help` consistency, `--strict` moved earlier, enum names case-insensitive |
| 2 | `0406baf` | unified `--method` selector | 1364+ / 44− | `--method <name>` + 4 old flags become aliases + native exits (error during the transition period); `src/set_module.f90` gains `xcontrol_supplied` |
| — | `3f155cc` | Fix an intrinsic `TMolecule` assignment in the relaxation engine | 125+ / 3− | **Independent upstream bug fix**: `src/relaxation_engine.f90` 2 lines `molopt = mol` → `call molopt%copy(mol)` (see §5.2) |
| 2b | `2b49ef8` | Expose the missing task knobs and name Hessian flags consistently | 1279+ / 18− | `--lopt-max-steps`, `--l*-hess-step` (5 spellings sharing a setter), `--l*-sccacc`, `--l*-restart`, `--l*-hessian` naming unified |
| 3a | `b822c1f` | Run the landscape suite for `--opt` and `--ohess` | 1461+ / 147− | `--opt`/`--ohess` re-dispatch + `[NOTICE]` + 3 guard rails (`$scan`/periodic/`$opt`) + 4 native exits + `fres` bridge + native block suppression |
| 3b | `8b85c6d` | Write the legacy files again and keep a failed `--ohess` from cascading | 659+ / 11− | `.xtboptok`/`NOT_CONVERGED`/`xtbopt.log`/`xtblast.xyz`/`compliance.dat` contract; stop on failure (removes the misleading cascade); `--namespace` reaches `xtbopt.xyz` |
| 3c | `fa969ce` | Let the convergence level drive all four axes | 538+ / 36− | The level drives 4 axes at once (4 thresholds + `e_tol` + `engine%acc` + `max_steps = clamp(3n,200,10000)`) + `e_tol` echo + the two halves of `--ohess <LEVEL>` agree |
| meta | `cdd9ffb` | Stop tracking the development and review material under `docs/` | 0 | only changes `.gitignore` + `git rm -r --cached docs/` (`docs/` was always session state, not a deliverable) |

---

## 3. Stage by stage

> The "breaking changes" column of each stage lists only **which §10 items that stage landed**; the complete per-item status of all 18 is in **§4**.

### 3.1 W0 — stop three silently wrong results + switch to `-O2`

**What changed** (`04-o2-and-silent-defects-fix.md`, 922 lines)
- `$fix` was violated by every second-order step → added `src/landscape/types/constraints.f90`, projecting the accepted step onto the free subspace.
- `--lvib` unconditionally overwrote the xcontrol `$thermo` → `vib/thermo.f90` changed to the four-key priority "explicit CLI > `$thermo` > builtin default".
- `--lopt` silently ignored `$scan` → `main.F90` gains a loud refusal.
- The build switched from `-O3` to `-O2` (`cmake/CMakeLists.txt` gains an `-O3` detection warning): under `-O3` gfortran 12.2.0 once produced a **wrong Hessian** for edcfbbe (water frequencies `-931.59/3329.78/3345.94` vs the correct `1557.41/3677.84/3684.01`).

**Breaking changes**: none (W0 is a bug fix, it does not enter the §10 list).

**Acceptance evidence**
- `ctest -R landscape` **100/100** (`04` §5.4) → close-out round **101/101** (`06-w0-closeout-report.md` §0).
- Behaviour after the three fixes: `$fix` all 8/8 displacements `0.00000e+00 Å`; `--lvib` vs `--hess` dG `1.4e-8 Eh`; `--lopt`+`$scan` exit 1 (`06` §3.3).
- Injection: delete the two guards of `landscape_vib_set_thermo_settings` → `xtb/landscape-vib-thermo-knobs` **red**; roll back the two lines of `write_json` → `xtb/landscape-vib-thermo` **red** (`06` §1.4/§2.1).

### 3.2 W1 — 4-way capability-boundary survey + full CLI redesign

**What changed** (`08-four-way-capability-boundary.md`, `11-cli-redesign.md`)
- Read-only survey: `02-cli-audit.md` audits all 122 suite flags item by item (read point / test coverage / dead-knob risk), `03-build-and-toolchain.md` audits the build chain.
- Design document `11-cli-redesign.md` (2851 lines of body + 41 revisions `R1..R45`): `--opt`/`--ohess` semantic decision (option C), the `--method` mechanism, the **18-item breaking-change list §10**, the 122-spelling mapping table §6, the backward-compatibility matrix §5, the test strategy §9.2, the stage split §9.1.
- The design went through **4 rounds of review** (`12/14/15/16`), revision record `13-cli-redesign-revision.md`.

**Breaking changes**: W1 defined and froze all 18 (the count = A3 + B0 + C8 + D1 + E1 + F4 + G1; C8 was registered later in stage 1, F4 landed in stage 1, the count went 16→17→18).

**Acceptance evidence**: the 4-way boundary report of `08` + the design self-check (the "command / expected / actual" tables of `11` §13/§16.4).

### 3.3 Stage 1 — non-breaking improvements

**What changed** (`17-stage1-impl-a.md`, `18-stage1-impl-b.md`, followed by three fix rounds `20/21/22/24`)
- Conflicting-combination detection: `select_runtyp` (`main.F90`) raises `env%error` on the second request for a **different** runtyp and **names the second flag**.
- Required-parameter errors name the flag (`--lscan-step`/`--lscan-atoms`/`--lscan-nstep`).
- Risky-default echo (`--lscan-fc`/IRC method/`--lopt-engine`/`--lts-escape-step`; D4 = the parent session ruled that the `--lscan` default bidirectional sweep is **not changed in code**).
- Two-way `--help` consistency script `test/unit/landscape_help_check.py` (set S difference = 0 as the gate).
- `--strict` moved earlier + suite enum values **case-insensitive** (`BFGS` == `bfgs`).

**Breaking changes (§10)**: **C1**, **C2**, **C3**, **C4**, **C5**, **C8**, **F2**, **F4**.

**Acceptance evidence**
- `ctest -R landscape` 101 → **109** (group A) → **115** (group B) → **118/118** (final); regression subset 12/12 (`21` §0, `22` §0).
- 8 test assertions for the conflict family; `--sp --lopt-engine bfgs` goes from exit 0 → **1** (`11` §10 C8's `[measured]`).
- Injection: `ts-can-tighten-case` reverse mutation red ("the upper case numerical source did not tighten"); the `.false.` mutation of the D3 banner **was still green at 115/115** → the fix round added falsifiable coverage for `algorithm_was_set` (`21` §3).
- `--help` set S: 124/124 (at stage 1) difference 0 (`19` §4).

### 3.4 Stage 2 — unified `--method` selector

**What changed** (`25-stage2-impl.md` + `26/27/28/29`)
- `--method <name>` (6 minimisers + 8 ts + 4 irc + 4 native engines + `auto`); `--lopt-engine`/`--lscan-method`/`--lts-method`/`--lirc-method` become **permanent aliases** (not deprecated, no warning printed).
- Two-segment normalisation at the end of parsing (suite set + native-exit set); conflicts (different values) error; **un-re-dispatched runtypes** error (`--metaopt --method bfgs` etc.).
- `src/set_module.f90` gains the `set%xcontrol_supplied` flag.

**Breaking changes (§10)**: **C6**.

**Acceptance evidence**
- `ctest -R landscape` 118 → **120** → **121/121** (fix round, `28` §0).
- `xtb/landscape-method-table`: 18 suite names + 4 `auto`, asserted effective name by name (`11` §9.2); shifting the name table ⇒ 38 lines turn red one by one.
- BC-4: the 4 rewritten old flags (`--lopt-engine`/`--lscan-method`/`--lts-method`/`--lirc-method`) produce output files that are **bit-identical** under `cmp` (`11` §9.1) - measured against the then-current builds; cross-build byte identity became round-off-level in the 2026-10-16 performance-PR round (§8.4).
- **A real experiment that genuinely turned red**: first-round `ctest` = 119/120, the only red item `xtb/landscape-help-coverage` (a new flag not in `--help`); after adding `xhelp.f90` 120/120 (`25` §5) — the set S gate really does stop it.

### 3.5 Stage 2b — fill in the missing task knobs + unify Hessian flag naming

**What changed** (`32-stage2b-impl.md` + `33/35/37/39`)
- Added `--lopt-max-steps N` (previously **missing**), `--lopt/lscan-hess-step`, `--lopt/lscan-sccacc`, `--lopt/lvib-restart`.
- The three names of `hess_step` (`--lts-step`/`--lvib-step`/`--lirc-step-hess`) unified to `--l*-hess-step`, 5 spellings going through **the same setter** (with different-value conflict detection + a `hess_step_source` provenance field).
- The hint threshold widened from "only watch the step-length flag" to "any flag of the same family" (`hess_acc_source` companion field).

**Breaking changes**: none (pure additions + aliases, BC-3 "the set is not smaller").

**Acceptance evidence**
- `ctest -R landscape` 121 → **125** → **127/127** (fix round, `35` §0; `39` §0).
- `xtb/landscape-hess-step-alias`/`-opt-knobs`/`-scan-knobs`/`-restart-knobs`: each knob is **double-asserted** (alias equivalence + a changed value must change the result).
- Injection: change the consumer back to the old constant ⇒ the corresponding assertion turns red (`--lopt-hess-step did not change the minimisation` etc., `32` §4 E1/E3/E4/E5).
- BC: 11 invocation groups, old/new binary output files **bit-identical** (then-current builds; round-off-level from the 2026-10-16 performance-PR round, §8.4); stdout gains only 2 lines `(not used: …)` (`39` §T1.4).
- `--help` set S: 135/135, difference 0 (`33` §5).

### 3.6 Stage 3a — `--opt` / `--ohess` re-dispatch (the breaking core)

**What changed** (`42-stage3a-impl.md` + `43/44/45/47/48`)
- `--opt [LEVEL]` / `-o` → the suite default optimiser (`p_run_lopt`), printing `[NOTICE]`; `--ohess [LEVEL]` → the suite opt+vib chain.
- **4 native exits**: `--method ancopt`/`rf`, `lancopt`/`l-ancopt`, `fire`/`inertial`, `pbc_lbfgs`.
- **3 guard rails**: `--opt` + xcontrol `$scan` / periodic system / `$opt` block ⇒ a loud error + a pointer to `--method ancopt`.
- Native block suppression (the `anyopt` branch of `main.F90`) + `fres` bridge (preserving the two native lines `TOTAL ENTHALPY`/`TOTAL FREE ENERGY`) + native Hessian/frequency path suppression (T1 `numhess`/T2 `main_freq`) + `xtbopt.xyz` second-write suppression + timing-label suppression.
- `--opt <LEVEL>` bridge (mapped to the suite level) + `--cycles N` bridge (mapped to `max_steps`) + `--opt banana` ⇒ `env%error`.
- `src/setparam.f90` gains `set%xcontrol_supplied` (correcting the field path of design §18.1).

**Breaking changes (§10)**: **A1**, **A2** (partial), **A3**, **C7**, **D1**, **F1**, **F3**, **G1** (the `check_scanblock` 492-503 section is rewritten).

**Acceptance evidence**
- `ctest -R landscape` 127 → **136** → **138/138** (fix round, `45` §3.3; `47` §0); regression subset 12/12.
- **20 native-branch CLI assertions** (`11` §4.7.1), including `--ohess` running only once (`Numerical Hessian` = 0, `landscape: harmonic frequencies` = 1, `landscape engine: Hessians` = 1).
- Fix of the false timing-label claim: before, `--ohess` printed `analytical hessian: … 0.000 sec / NaN speedup`; after, `grep -c 'analytical hessian'` = **0**, and the reverse `--ohess --method ancopt` = **1** (`45` §2 T1).
- Injection: hard-wire `--opt crude` and `--opt tight` to the same level ⇒ `FAIL: --opt crude and --opt tight produced the same run` (`45` §2 T2).
- **Design premises refuted** (honestly registered): an explicit `--lopt-level` "higher than" `--opt <LEVEL>` is **unreachable** on the CLI (`--lopt-level` itself triggers a runtype conflict and exits 1 before the bridge); `--cycles 8` does **not bind** for the water fixture (default 5 steps) ⇒ the test switched to `--cycles 4` (`45` §4.3/§4.4).

### 3.7 Stage 3b — file contract + stop on failure (and fix the regression introduced by 3a)

**What changed** (`49-stage3b-recon.md`, `50-stage3b-impl.md` + `51/52/53/54`)
- The suite `--lopt`/`--ohess` writes the full upstream semantics of `.xtboptok` (only on convergence, empty file)/`NOT_CONVERGED` (on failure)/`xtbopt.log` (multi-frame)/`xtblast.<ext>` (on failure).
- The suite `--lvib` writes `compliance.dat` (bare name, no namespace, same as upstream).
- **Stop on failure**: `main.F90` guards the whole vib half with `if (lres%converged)` ⇒ removes the misleading cascade (see §5.4).
- `--namespace` reaches `xtbopt.xyz`/`xtb.trj`/`xtbopt.log` (`types/log.f90`).
- **Fix of the HIGH-2 introduced by 3a**: the help/man texts of `--bhess`/`--omd`/`--metaopt` reverted to the facts (see §5.3).
- **Fix of the HIGH-1 inherited from HEAD**: the `anyopt` predicate + the `murks` initial value of `main.F90` (see §5.3).

**Breaking changes**: no additions (class B stays 0; the change of `--lvib`'s stdout/file set is registered as the third intentional difference of BC-4, `11` §10 class B `[registered T4]`).

**Acceptance evidence**
- `ctest -R landscape` 138 → **143** → **145/145** (fix round, `53` §0; `54` review conclusion = **accept, can be committed**).
- The measured three states of the `--omd` regression: pre-3a exit 0 / `ANC micro-cycles` 1 / `xtbopt.xyz` = −5.070544363296; 3a+3b before the fix exit 0 / **0** / −4.962807380732 (= the input geometry); after the fix exit 0 / **1** / −5.070544363296 (`53` §1.2).
- Injection: restore the `main.F90` predicate to `anyopt .and. landscape_native_engine /= 0` ⇒ `xtb/landscape-native-runtypes-unaffected` **red** (`FAIL: --omd did not run the native ANCopt`) (`53` §1.3).
- BC-4: `--omd`/`--screen`/`--metaopt`'s `xtbopt.xyz`/`xtbopt.log` **bit-identical** to pre-3a (`53` §1.4) - then-current builds; round-off-level from the 2026-10-16 performance-PR round (§8.4).
- 5 file-contract tests (positive + negative: single atom/`$fix atoms:` does not write `compliance.dat`).

### 3.8 Stage 3c — the level drives four axes (the only stage that moved numbers)

**What changed** (`55-stage3c-impl.md` + `56/57/58`)
- `--l*-level` drives **4 axes** at once: 4 thresholds + `e_tol` + `engine%acc` + `max_steps = clamp(n_atoms × multiplier, 200, 10000)` (the `max_steps` of the default `normal` level goes from 256 → 200@H2O).
- Explicit `--acc`/`--l*-max-steps`/`--cycles` take **priority** (`max_steps_source` + `acc_was_set` provenance fields).
- `e_tol` gains a stdout echo (otherwise that axis is not falsifiable).
- `acc_was_set` changed from a module-level mutable variable to an `intent(out)` dummy of `parseArguments` (review M1).
- Four `(default 256)` strings in `xhelp.f90` changed to level-derived wording (review M2).
- **Fix of the two-halves inconsistency of `--ohess <LEVEL>`**: the bridge now also sets `lvset%set_conv_level` (see §5.7).

**Breaking changes (§10)**: **A2** (complete), **E1**.

**Acceptance evidence**
- `ctest -R landscape` 145 → **146** → **147/147** (close-out round, `58` §5); regression subset 12/12; `--help` consistency PASS; no duplicate `add_test`.
- The four-axis numbers checked value by value: `settings.f90:1165-1199` vs upstream `optimizer.f90:59-105` (all 5 levels × 3 axes agree, `57` §1.1-1.2).
- The new test `xtb/landscape-level-axes` pins three axes + explicit priority; the existing `check_acc` reverse-pins "an explicit `--acc` must take effect".
- **Re-baseline conclusion: no accepted number moved** (M3 dfreq 0.0100 cm⁻¹, M4 prfo-vs-bprfo 6.3022e-06 bohr, M5 LEPS 8.9985 kcal/mol, M6 four-method MEP deviation 0.03995/0.00657/0.00215/0.02574 bohr); the only intentionally updated golden = the default budget pinned by `xtb/landscape-opt-knobs`, 256 → 200.
- Injection: (a) change `if (acc_was_set)` back to unconditional ⇒ the acc of `--lopt-level loose` returns to 1.000000 (the dead knob reappears) ⇒ `xtb/landscape-level-axes` **red**; (b) bridge `--ohess` with the pre-fix binary ⇒ `['0.200000','1.000000']` exit 1 / after the fix `['0.200000','0.200000']` exit 0 (`56` §1; `58` §3.3).
- Final review conclusion: **conditionally accepted** (M1/M2 fixed) → after the close-out round `58` records it as deliverable (`57` §0).

---

## 4. Per-item status of the breaking changes (design §10)

**Authoritative count** (design `11-cli-redesign.md` §10 summary row): **A 3 + B 0 + C 8 + D 1 + E 1 + F 4 + G 1 = 18 items**.

| Item | Change | Status | Landed in stage | Verifiable evidence |
|---|---|---|---|---|
| **A1** | the `--opt` engine and the numerical result change (ANCopt → the suite) | landed | 3a | `xtb/landscape-opt-alias` (stdout contains `[NOTICE]` + `landscape: lbfgs`; `xtbopt.xyz` bit-identical to `--lopt`) |
| **A2** | the criterion of `--opt normal` changes (four criteria + level-dependent `e_tol`/`acc`/`max_steps`) | landed | 3a+3c | `xtb/landscape-opt-level-map`, `xtb/landscape-level-axes` |
| **A3** | `-o` re-dispatched together with all levels | landed | 3a | `xtb/landscape-opt-level-bridge` (`--opt tight` == `--lopt --lopt-level tight`, files bit-identical) |
| **class B** | file contract changes | **0 items** (R4 cleared it: the suite writes all 4 legacy files) | 3b | `xtb/landscape-file-contract-opt`/`-vib`/`-namespace` |
| **C1** | `--lopt --lscan` silently first-option-wins → error | landed | 1 | `xtb/landscape-runtype-conflict` |
| **C2** | `--opt --lopt` / `--opt --hess` → error | landed | 1 | `xtb/landscape-opt-lopt-conflict` |
| **C3** | `--lts --lopt-engine bfgs` → error (naming `--lts-method`) | landed | 1 | `xtb/landscape-subflag-mismatch` |
| **C4** | `--hess --lvib-temp 400` → error | landed | 1 | `xtb/landscape-vib-temp-mismatch` |
| **C5** | any native runtype + suite sub-option (the whole class) → error | landed | 1 | the general rule of `xtb/landscape-subflag-mismatch` |
| **C6** | an un-re-dispatched runtype + `--method` → error | landed | 2 | `xtb/landscape-method` (`--hess --method bfgs` etc.) |
| **C7** | a suite runtype + xcontrol `$opt` block → error | landed | 3a | `xtb/landscape-opt-xcontrol-scope` (`--lopt --input <with $opt>` non-zero + naming `--method ancopt`) |
| **C8** | two or more mutually exclusive runtype flags → error (naming the second) | landed | 1 | `xtb/landscape-runtype-conflict`; `--sp --lopt-engine bfgs` 0 → 1 |
| **D1** | `--cycles INT` mapped to the suite `max_steps` (`--lvib` errors) | landed | 3a | `xtb/landscape-cycles-map`, `-cycles-lvib`, `-cycles-conflict` |
| **E1** | `--l*-level` drives `e_tol`/`acc`/`max_steps` at the same time | landed | 3c | `xtb/landscape-level-axes`, unit case `level-step-budget` |
| **F1** | `-o <level>` → the suite | landed | 3a | same as A3 |
| **F2** | `--strict --zzzz` exit 0 → non-zero (position-independent) | landed | 1 | `xtb/landscape-strict-unknown` |
| **F3** | `--opt banana` silently normal → error | landed | 3a | `xtb/landscape-opt-level-bridge` (`Unknown --opt level 'banana'`) |
| **F4** | suite enum values case-insensitive (`BFGS` == `bfgs`) | landed | 1 | the `enumcase` unit-case series + `xtb/landscape-method-table` |
| **G1** | `xtb/landscape-scan-xcontrol` turns red (a test-layer consequence, must be rewritten in step) | handled | 3a | `check_scanblock` lines 492-503 changed to "assert a loud error + name it in the text" |

**Open item (not blocking, honestly registered)**: **U17** — under the suite `--ohess`, `--o1nh` (ODLR) is a dead knob (`set%o1numhess` vs `settings%engine%odlr` have different sources). If adopted the count goes 18 → 19 (candidate **C9**, C8 is already taken by the mutually exclusive runtype flags). See design §12 U17.

---

## 5. Real defects fixed (the most valuable section of this report)

> Each entry gives: **defect / symptom / root cause / fix / regression test**. This section lists only "code changed and proven by a test", not pure documentation corrections.

### 5.1 The three **silently wrong results** of W0

| # | Defect | Symptom (measured before the fix) | Root cause | Fix | Regression test |
|---|---|---|---|---|---|
| 1 | **`$fix` violated by every second-order step** | For water with `$fix atoms: 1` (O at the origin), the RFO step moved the fixed atom **0.0098 Å (x) / 0.0109 Å (z)**; SD/CG/L-BFGS/BFGS/GDIIS were spared because their gradient is zeroed | The calculator layer only zeroes the gradient of `$fix` atoms (`src/xtb/calculator.f90:304-306`), whereas the RFO step expands in the Hessian eigenbasis and mixes in the frozen components | Added `types/constraints.f90`: every accepted step is projected onto the free subspace (`fixset` as the single source of truth) | `xtb/landscape-fix-rfo` etc.; after the fix all 6 engines drift `0.00000e+00 Å` (`AGENTS.md` §6 measured section) |
| 2 | **`--lvib` unconditionally overwrites the xcontrol `$thermo`** | `$thermo temp=298.15,300.0,400.0` collapsed to the single 298.15 K, `scale=0.9` ignored, `--lvib` vs `--hess` differed by **7.5e-3 Eh (4.7 kcal/mol)** | The old `vib/thermo.f90` copied all four keys | Changed to "explicit CLI > `$thermo` > builtin default" (four `*_was_set` guards); validation moved to the **merged** effective state | `xtb/landscape-vib-thermo`/`-thermo-knobs`; deleting the guards by injection ⇒ red |
| 3 | **`--lopt` silently ignores the xcontrol `$scan`** | `--lopt --input <with $scan>` exit **0**, writes no `xtbscan.log`, prints nothing | `relaxed_scan` runs only on the `anyopt` path (`src/prog/main.F90:914-916`) | `--lopt/--lscan/--lts/--lirc/--lvib` + `$scan` ⇒ `env%error`, pointing at `--opt` (or `--lscan --lscan-coord`) | `xtb/landscape-scan-xcontrol` (lines 466-491 keep the suite runtype refusal assertion; lines 492-503 are rewritten to a loud error) |

**The "prevent the next one" artefact of W0**: `test/unit/landscape_cli_check.py` (602 lines) + `landscape_testpot.f90` — from then on every suite option has a CLI-level "prove it takes effect" assertion.

### 5.2 The intrinsic-assignment defect of upstream `relaxation_engine.f90` (an independent upstream bug, not among the 18)

- **Defect**: the L-ANCopt engine copied the geometry with an **intrinsic derived-type assignment** (`src/relaxation_engine.f90:615 molopt = mol` / `:779 mol = molopt`).
- **Symptom**: gfortran 12.2.0 corrupts `sym(2:)` of a `TMolecule` containing `character(len=:), allocatable :: sym(:)` (a stride defect, `AGENTS.md` §6 item 5). When built with this machine's compiler, the **L-BFGS path of native `--opt`** (`$opt engine=lbfgs`, or the default engine for >500 atoms) wrote **garbage bytes** (NUL/random bytes) into the symbol column of `xtbopt.xyz` / `xtbopt.log` / `xtbtopo.mol` and into the `final structure:` / bond table of stdout. **The numbers were all correct**, so it went unnoticed.
- **Root cause**: a gfortran 12.2.0 compiler defect (reportable to GCC Bugzilla); the `rf`/FIRE engines already wrote `call molopt%copy(mol)`, so they were clean.
- **Fix**: **2 lines** changed to an explicit `%copy` (same spelling and same signature `copyMolecule(self, mol0)` as the `rf` path).
- **Regression test**: `xtb/native-opt-lbfgs-symbols` (byte level: strict UTF-8 + control-byte scan + the symbol column against the input geometry; restoring the two lines by injection ⇒ red, restoring ⇒ green).
- **Numerical-neutrality proof**: `xtbopt.xyz` 12 / `xtbopt.log` 60 / `xtbtopo.mol` 10 / normalised stdout 389 numeric-token sequences are **identical**, only the symbol column goes from garbage → correct (`11` §10 R42 section, `30-upstream-relaxation-fix.md`).
- **Counting convention**: it **does not belong** to the 18 CLI breaking changes, it is an independent upstream fix. **Note**: the same defect still exists in `src/screening.f90:107` (`mol = mol0`), which this direction **did not fix** (out of scope, `AGENTS.md` §6 item 5, last paragraph).

### 5.3 The `--omd`/`--screen`/`--metaopt` regression introduced by 3a and fixed by 3b

- **Defect**: 3a wrote the suppression predicate of the native optimiser block as `if (anyopt .and. landscape_native_engine /= 0)` (`main.F90`).
- **Symptom**: `anyopt` covers `--omd`/`--screen`/`--metaopt`, but these three runtypes were **never re-dispatched** and `--method` **errors** for them ⇒ `landscape_native_engine` is always 0 for them ⇒ **the native optimiser block is unreachable**. `--omd` ran MD directly from the **input geometry** (no preceding optimisation); `--metaopt` read the **uninitialised** `murks` at `main.F90:1524` (UB, which could randomly trigger `Optimization did not converge, aborting`).
- **Root cause**: the meaning of `anyopt` ("any runtype of the optimiser family") is wider than "re-dispatched to the suite".
- **Fix**: the predicate changed to `if (anyopt .and. .not. landscape_is_suite_run(set%runtyp, landscape_native_engine))` (reusing the existing predicate, no second criterion) + `murks = .false.` initial value.
- **Regression test**: `xtb/landscape-native-runtypes-unaffected` (`--omd` exit 0 + `ANC micro-cycles` appearing 1 time + `xtbopt.xyz` golden energy 1e-9 + O–H ∈ (0.9,1.05) Å, **not just the exit code**).
- **Measured three states** (`53` §1.2): pre-3a `-5.070544363296` / 3a+3b before the fix `-4.962807380732` (= the input geometry) / after the fix `-5.070544363296`.
- **HIGH-2 fixed at the same time**: 3b had changed the help/man text of `--bhess`/`--omd`/`--metaopt` to "the geometry is optimised by the suite (pass `--method ancopt`)" — **neither statement holds** (`--method ancopt` errors directly for these three runtypes, `main.F90:4368-4379`). Review `52`'s P1-1 caught and reverted it; `src/xhelp.f90` is **0 diff** relative to HEAD.

### 5.4 The **misleading cascade** of the `--ohess` failure path ("Could not generate lattice points")

- **Defect**: after a failed minimisation the suite `--ohess` **still ran the vib half** and reused the same contaminated `calc`/`env`.
- **Symptom** (measured before the fix, `49` §D12): `--ohess --cycles 2` printed a 7-layer error stack, with the root cause `landscape minimisation did not converge` buried at the bottom, while the top showed the **misleading** `type_latticepoint_update: Could not generate lattice points` and `Setup of Coulomb evaluator failed`; exit 1; **no** `hessian`/`vibspectrum` written.
- **Root cause**: two coupled defects — (i) the chain ran the vib half without checking `lres%converged`; (ii) `latticepoint` raised its own error while `env` already had a pending error (`src/type/latticepoint.f90:198`), turning the second stage's failure into "lattice points".
- **Fix**: `main.F90` guards **the whole vib half** with `if (lres%converged)` (structurally sealing the cascade: the vib half is the only consumer that reuses the same `calc`). **It does not change** `src/type/latticepoint.f90` (core layer, forbidden by AGENTS §5).
- **Regression test**: `xtb/landscape-ohess-fail` (non-zero exit, no cascade wording, no `hessian`/`vibspectrum`, `NOT_CONVERGED` present and no `.xtboptok`).
- **Intentional deviation**: native `--ohess --method ancopt --cycles 1` **exits 0** and still writes `hessian`/`vibspectrum`/`.xtboptok` (measured `49` §D12(b)); the suite exits 1 + `NOT_CONVERGED` + writes no frequency file. The reason (single-cause message / the meaning of frequencies at a non-stationary point is dubious / `NOT_CONVERGED` is an existing upstream signal) is registered in design §3.9 — **see §7.3**.

### 5.5 **4 output-naming defects** (the `--namespace` series)

This class recurred, and is the latent trap of "two conventions left in one code base" (`AGENTS.md` working-agreement section):

| # | Milestone | Defect | Regression test |
|---|---|---|---|
| 1 | M3 | `--namespace` **double-prefixed** `hessian` (`ns.ns.hessian`) | `xtb/landscape-vib-namespace` (write→read round trip) |
| 2 | M4 | `--namespace` did not reach the output of `--lts` | `xtb/landscape-ts-namespace` |
| 3 | M6 | the output of `--lirc-*` was not prefixed | the namespace assertion of `xtb/landscape-irc-*` |
| 4 | M7 | `--namespace` **did not reach the geometry/trajectory of `--lopt`**: two minimisations with different namespaces in the same directory **overwrote each other's** `xtbopt.xyz` | `xtb/landscape-opt-namespace` + injection xiii (`docs/xtb-landscape-M7-validation.md` §5) |

The 4th was found by **the parent session's own end-to-end run**; at the time the suite's chain test used an **unprefixed** name and therefore **bypassed** it (`53-m7-final-fix-report.md` §1). The fix is at the write point `src/landscape/types/log.f90`; the default names (`xtbopt.xyz`/`xtb.trj`) are unchanged and asserted.

### 5.6 Dead knobs ("parsed but read by nobody" or "claimed effective but not effective")

> **Honest note on the counting convention**: design §3.5 records "this project has caught **4** of them" (M3/M4/M6×2); `48-m7-e2e-closeout-report.md` records "three earlier dead knob defects". The 3a review (`44` LOW-1) caught **2** more (`--opt <LEVEL>`, `--cycles N`), and 3c caught **1** more (the acc axis of `--l*-level`). **I cannot find a source for the "5" given in the task brief**; the table below lists **7** by **verifiable source**.

| # | Knob | Symptom | Fix/evidence |
|---|---|---|---|
| 1 | M3 `engine%etemp` | the field was dead: `xtb_single::singlepoint` dropped the `et` argument | the engine pushes `etemp` into `TxTBCalculator%etemp`; `xtb/landscape-vib-etemp` (stretched H2 300 K vs 5000 K, the only mode 579.9 → 397.6 cm⁻¹) |
| 2 | M4 `results%etemp` | the `etemp` in the result had a different source than the effective value | **read back** from the calculator; `xtb/landscape-ts-etemp` |
| 3 | M6 `--lirc-level` | parsed, but `relax_endpoint` only copied `write_files/printlevel/engine/etemp`, the four criteria were not passed down | `xtb/landscape-irc-knobs` (1 step vs 7 steps); injection red |
| 4 | M6 `--lirc-rmsd-end` | set the field but `endpoint_reached` did not read it | `xtb/landscape-irc-knobs` (1e-3 ⇒ F/T) + unit case `endpoint-rmsd-criterion` |
| 5 | 3a `--opt <LEVEL>` (`set%optset%optlev`) | after re-dispatch **no reader at all** ⇒ `--opt tight` equivalent to the suite default level (**silently ineffective**) | 3a fix round T2: bridged to `lset%set_conv_level`; `xtb/landscape-opt-level-bridge` |
| 6 | 3a `--opt --cycles N` (`set%optset%maxoptcycle`) | same as above ⇒ silently discarded | 3a fix round T3: bridged to `max_steps`; `xtb/landscape-cycles-map` |
| 7 | 3c the **acc axis** of `--l*-level` | the level-derived `acc` was overwritten by the unconditional `engine%acc = set%acc` (= 1.0) ⇒ always 1.0 | 3c T2: 6 write points changed to `if (acc_was_set)`; `xtb/landscape-level-axes`; injection red (`55` §T2, `56` §1) |

**The defensive artefact of this direction**: the `hess_step_source`/`hess_acc_source`/`max_steps_source` provenance fields of `types/settings.f90` + the `--help` echo — making every switch that "claims to take effect" **require** a test that can falsify it (the lesson of the `AGENTS.md` working-agreement section).

### 5.7 Other real defects (item by item)

| Defect | Symptom | Root cause | Fix | Regression test |
|---|---|---|---|---|
| **The false claim `analytical hessian: 0.000 sec`** | at the end of the suite `--ohess`, `analytical hessian: … 0.000 sec / ratio c/w: NaN speedup` was printed, but the suite **never** does an analytic Hessian | the timing label of `finalize_xtb` was decided by a three-way test on the **bare runtyp** (`p_run_hess/ohess/bhess`) | slot 5 changed to read `landscape_native_hessian_path` (the **same criterion** as the suppression) | `xtb/landscape-ohess-suite-once` (`analytical hessian` appears **0** times) + `-ohess-native-hess` (reverse: **1** time) |
| **`--opt tight` silently ineffective** | after re-dispatch `--opt tight` and giving no level **produce the same result** | see §5.6 #5 | 3a fix round T2 | `xtb/landscape-opt-level-bridge` |
| **`--ohess <LEVEL>`'s two halves disagree** | in one `--ohess tight` run `landscape engine: SCC accuracy` is printed as **0.2** (opt half) and **1.0** (vib half) — a user-visible self-contradiction | the bridge only did `call lset%set_conv_level(...)`, `lvset` was untouched | the bridge now also sets `lvset%set_conv_level` for `p_run_ohess` | `xtb/landscape-ohess-level-axes` (3 × `--ohess`, asserting that the two halves' acc/e_tol each appear 2 times and are equal) |
| **`--raman` bypasses `select_runtyp`** | `--raman 0.5 1000 --lopt` reports no conflict and `--lopt` is silently dropped (`--raman` goes through `set_elprop('raman')` → `set_runtyp('hess')`, bypassing conflict detection — an instance of "validated on one side, broken on the consuming side") | the `--raman` bypass in `main.F90` | `select_runtyp` gains an optional `apply` argument (`apply=.false.` only records, does not call `set_runtyp`); `--raman` uses it | `xtb/landscape-raman-conflict` (non-zero exit + naming `--raman` and `--lopt` + `run type conflict`); a reverse mutation that deletes the call ⇒ the error becomes `Raman activities are not implemented` (`21` §1.2) |
| **The two-way gap between `--help` and the implementation** | stage 4 not done: many native flags are still "implemented but undocumented" | `xhelp.f90` not synchronised | stage 1 establishes the set S gate (suite-flag difference = 0), the native set N is only printed | `xtb/landscape-help-coverage` `[measured · this report]`: S(main)=135=S(help), difference **0** |

---

## 6. Methodology (the most distinctive part of this project)

### 6.1 The review loop

Every stage goes through **implement → independent validation → independent review (fresh-context, read-only) → fix → re-review**. This direction had ~20 rounds of review in total, and the three final stage conclusions: 3a = conditionally accepted (`44`, LOW×4/INFO×7 → `48` fix-round review), 3b = **accept** (`54`, CRITICAL 0/HIGH 0/MEDIUM 0/LOW 3/INFO 4), 3c = **conditionally accepted** (`57`, MEDIUM 2 → `58` close-out fix).

The review caught problems **that both the implementer and the validator missed**: the HIGH-2 false text of 3b (`52` P1-1), the module-level mutable variable of 3c (`57` M1), and the refutation of C7's blast radius (design §18.2: the real channel is `XTBPATH` → `.xtbrc` at the repo root, not the test files). **"Text-class regressions have no automatic gate"** (`landscape_help_check.py` only compares flag sets) is explicitly registered by review `54` INFO-2 as a residual risk — that is exactly how HIGH-2 escaped validation.

### 6.2 Injection experiments (every assertion must be able to turn red)

Rule: **every new test must prove that it is red under the old code**. This direction's injection recipe: back up the source → single-point text mutation → incremental build → run the corresponding test (expect red) → restore byte by byte → re-green (`sha256sum -c` asserts the restore). Every stage report pastes the real **red-before / green-after** output.

**Lesson (the 4th static conclusion refuted by measurement, `34-m4-coverage-fix-report.md` §): an injection must first prove the mutation is **semantically equivalent** to the old logic. A mutation that is "different but happens to be a correct alternative implementation" **stays green**, and in that case "green" proves nothing.

### 6.3 The same-source principle (hint and count from the same source ⇒ contradiction impossible by construction)

**The same pure predicate/provenance field** drives all read points; a second criterion never appears:
- `effective_landscape_task` / `landscape_is_suite_run` / `landscape_ohess_chain` / `landscape_native_hessian_path` serve at once the guards, the native block suppression, `--ljson`, the timing labels, `[NOTICE]`, and the `xtbopt.xyz` second write.
- All six write points of 3c's `acc` go through `acc_was_set`; `max_steps_source` is written by the 4 flags + the `--cycles` bridge and read by `set_conv_level`.
- The `--ljson` counters and stdout have the **same source** (contradiction impossible by construction), with 17-digit reals, `null` for non-finite values, and byte determinism (same `OMP_NUM_THREADS`).

### 6.4 Instances of "the tests are all green but still wrong" (≥10, each verifiable)

| # | Instance | Why green | Source |
|---|---|---|---|
| 1 | M0/M1: L-BFGS returned an uphill direction, the line search rejected it and fell back to SD, **converging as usual** | the old assertion only looked at `converged` | `09-h3-test-falsifiability-report.md`, `08-independent-review-m0-m1.md:89` |
| 2 | M2: the rigid-dihedral bug was **encoded as expected behaviour** | the test treated the bug as the expectation | `AGENTS.md` §2 M2 section |
| 3 | M2: the task-level `result%converged`/exit code contradicts its own docstring | the CLI assertion was hit by a substring | `17-m2-independent-review.md:150` |
| 4 | M3: extra near-zero modes → uninitialised index (HIGH-1) | the green suite **never took** that branch | `22-m3-independent-review.md:60` |
| 5 | M3: the `water-vs-trproj` unit case was a **no-op** (reprojecting an already projected matrix) | the assertion is always true | `AGENTS.md` §2 M3 section |
| 6 | M4: non-stationary geometry exit 0 | the CLI test only ran the positive HCN case | `29-m4-independent-review.md:100` |
| 7 | M4: the H3 injection was **not a semantically equivalent** mutation | the mutation was a "correct alternative implementation" | `34-m4-coverage-fix-report.md` |
| 8 | M5: `grow_frontier` argument misordered (`other` passed as `budget`), the first node toward the coordinate origin | the harmonic-oscillator fixture could not see it, and end-to-end it was absorbed by the spherical relaxation | `38-m5-independent-review.md:66` |
| 9 | M6: splitting `gs.f90` silently changed the walk | only the endpoint test caught it | `AGENTS.md` §2 M6 section |
| 10 | M7: `--namespace` did not reach the `--lopt` geometry | the chain test used an **unprefixed** name, bypassing that path | `53-m7-final-fix-report.md` §1 |
| 11 | stage 1: changing the D3 banner condition to `.false.` ⇒ **still fully green at 115/115** | no test referenced `algorithm_was_set` | `19-stage1-validation.md:375`, `20-stage1-review.md:218` |
| 12 | stage 1: `--lscan-nstep-missing`'s pattern was hit by the **command-line echo** (false green) | the assertion only asserted the flag name | `18-stage1-impl-b.md` §8.1 |
| 13 | stage 3b: the help text said the opposite | `landscape-help-coverage` only compares the flag **set**, not the semantics | `52-stage3b-review.md` P1-1, `54` INFO-2 |

### 6.5 The list of "facts corrected by measurement"

**The discipline of this direction**: static conclusions must be measured. The stage reports give the count themselves: the stage-1 fix round records **"the 13th of this direction (6 of them the parent session)"** (`21-stage1-fix-report.md:81`); the 3c close-out records **"the 14th"** (`58-stage3c-closeout.md:205`). **I cannot find a direct source for the "7 parent-session times" given in the task brief**; the table below lists by verifiable source (including the parent session's):

| # | Refuted conclusion | Who raised it | Measured truth | Source |
|---|---|---|---|---|
| 1 | "CI will go red: 4 new tests do not pin `XTBPATH`" | review | under `env -u XTBPATH` 4/4 pass (GFN2 parameters are compiled into the binary) | `21` §1.1 |
| 2 | "`--lscan-nstep-missing`'s flag-name assertion is falsifiable" | the implementer themself | it is unconditionally hit by the command-line echo (false green) | `18` §8.1 |
| 3 | "an explicit `--lopt-level` is higher than `--opt <LEVEL>` (different values error)" | design | unreachable (the runtype conflict fires before the bridge) | `45` §4.3 |
| 4 | "`--cycles 8` is below the default" | task brief | the water fixture default is 5 steps, 8 does **not bind** | `45` §4.4 |
| 5 | "BC-4 = stdout is also line-identical" | design | not robust (timing lines without a percent sign are missed) | `45` §4.5 |
| 6 | "the parent session's three `--ohess` energies are not reproducible" | the parent session | the measurement method was wrong (`head -3` caught the initial geometry) | `45` §4.1 |
| 7 | "`--opt --method ancopt` in two orders gives bit-identical stdout" | design | even two runs of the same binary are not bit-identical (24 timing/seed lines) | `25` §6 F3 |
| 8 | "`case (p_run_lopt, p_run_opt, p_run_ohess)` all enter the suite setter" | design | it would create a dead knob; changed to an error during the transition period | `25` §6 F4 |
| 9 | "C7's blast radius = 0 (0 hits in test files)" | design | the real channel is `XTBPATH` → `.xtbrc` (143/216 tests carry `XTBPATH`) | design §18.2 |
| 10 | "`set%optset%xcontrol_supplied`" | design | it is actually `set%xcontrol_supplied` (the line number the design cites points into `TSet`) | design §18.1 |
| 11 | "review INFO-1: the two halves of `--ohess <LEVEL>` having different accuracy is **not** a defect and should stay as is" | review | **Refuted**: printing 0.2/1.0 in one run is a user-visible self-contradiction | `58` §6.1 |
| 12 | "`--bhess`/`--omd`/`--metaopt` have their geometry optimised by the suite (pass `--method ancopt`)" | 3b implementation | `--method ancopt` errors directly for these three runtypes | `52` P1-1 |
| 13 | "the `compliance.dat` condition is isomorphic to upstream" | implementation/design | upstream reads `freezeset` (`$fix freeze:`), the suite reads `fixset` (`$fix atoms:`) — two keys | `52` P2-1 |
| 14 | "`--lopt-level loose`'s acc = 2.0 takes effect" | before 3c | it is always 1.0 (a dead knob) | `55` §T2 |

**The parent session's own error types** (summarised from the items above): **measurement-method errors** such as relative paths/`head -3` (#6), **second-hand information** (#12), **an old binary/old build directory** (`29-stage2-fix-review.md:105` states explicitly that the directory "can no longer be used as evidence"), and **treating a static inference as a conclusion** (#9/#13). **The operating discipline was therefore upgraded**: always absolute paths, confirm the binary/path before verifying, and compare against a baseline after running.

---

## 7. Residual risk and known limitations (**honesty is the whole value of this section**)

### 7.1 The rc channel (`.xtbrc`) is still **silently ignored** by the suite, and it is a **whole class**

- **Fact**: a `.xtbrc` in the directory pointed to by `XTBPATH` is read in by `rdcontrol`. The suite **consumes** `$scc` (`$scc temp=5000` ⇒ `landscape engine: etemp 5000.000 K` for `--lopt`) and `$thermo` (`$thermo temp=500.0` ⇒ `H(T)` of `--lvib` changes), but **does not consume** `$opt`/`$md`/`$siman`/`$hess`/`$modef`/`$cube`/`$symmetry` (`grep -rn 'set%optset|set%mdset|…' src/landscape/` = 1 hit and it is a comment).
- **Ruling**: C7's guard decides by **source** — it takes effect only for `$opt` in the **input file** and is **silent** for `$opt` in `.xtbrc` (user/parent-session ruling (A): "least hostile to the user"; the upstream release ships its own `.xtbrc` whose `$opt` block, if fatal, would hard-fail all such users).
- **Registered**: the whole **class** "rc-channel blocks not consumed by the suite" is registered as **to be re-evaluated** (re-evaluate whether it should give a loud error once there are real users). See design §18.3.
- **Note**: "the rc is always ignored" is **wrong** (`$scc`/`$thermo` are consumed through the global parameter state).

### 7.2 The meson entries are gone

`meson`/`ninja` are not installed (`AGENTS.md` §4). All new tests are registered in `test/unit/CMakeLists.txt`, the one build system of the tree; when this was written a hand-synchronised `test/unit/meson.build` existed as well and was never executed, and it was deleted together with the whole Meson build in `e714161`. Substitute checks: test names `uniq -d` = 0, bracket balance = 0, the python checker itself executed directly and passing (`25` §7, `39` §5, `54` §5). The Meson build no longer exists, so there is nothing left to verify there.

### 7.3 The `--ohess` failure path and an **intentional deviation** from native

| | native `--ohess --method ancopt` | suite `--ohess` |
|---|---|---|
| exit code when the minimisation does not converge | **0** | **1** |
| does it still write the frequency files | **yes** (`hessian`+`vibspectrum`+`.xtboptok`) | **no** (writes `NOT_CONVERGED`) |
| measured | `49` §D12(b) | `53` §1.3, `xtb/landscape-ohess-fail` |

The reason is registered in design §3.9: a single-cause message / the dubious meaning of a frequency analysis on a non-stationary geometry / `NOT_CONVERGED` being an existing upstream signal for "not converged". **This is an intentional difference, not a bug**, but external scripts that rely on "a native failure still produces frequencies" must be aware of it.

### 7.4 The 3c re-baseline deviations (accepted by a parent-session ruling, item by item)

| Milestone | Item | Ruling and reason | Evidence level |
|---|---|---|---|
| M2 | the 77 hysteresis files (`docs/_landscape-recon/m2-hysteresis-data/`) are **not re-run** | accepted: that batch was measured under the **old level semantics**, re-running it would give a different experiment | historical record, **not used** as a current golden reference |
| M4 | dimer-vs-prfo `7.7397e-05` | ~~accepted as a deviation~~ **correction: this item is not a deviation, it has been measured and closed** — `test/unit/test_landscape_saddle.f90:2004-2007` **does print** `rmsd(prfo)` (`es12.4`); `56-stage3c-validation.md:122/:177` measured `7.7397E-05` (**unchanged**) and noted "this item can be struck". The first draft of this report inherited "the print point no longer exists" from `58-stage3c-closeout.md` §4 — that is a **false statement** (design `11-cli-redesign.md:925` has been corrected in step) | **measured, unchanged** (not a static judgement) |
| M3 | CH3OH `dfreq 0.0600 cm⁻¹` is **not re-run** | accepted: a manual experiment, not in ctest; the H2O `dfreq 0.0100` of the same chain is already covered by the registered test `xtb/landscape-e2e-chain` and measured unchanged | the H2O of the same chain has registered-test coverage |

Source: `58-stage3c-closeout.md` §4 + the registration table of design §3.10. **The only intentionally updated golden** = the default budget pinned by `xtb/landscape-opt-knobs`, 256 → 200.

### 7.5 Stage 4 (documentation close-out) is **not done**

`xhelp.f90` still has documentation gaps for native flags. `[measured · this report]` `python3 test/unit/landscape_help_check.py src/prog/main.F90 src/xhelp.f90`:

```
dialect A: S(main) = 135, S(help) = 135, difference = 0        <- suite flags: gate PASS
dialect A: N(implemented, undocumented) = 33                   <- implemented, undocumented
dialect A: N(documented, unimplemented) = 4                    <- --ithr --silent --xparam -s
dialect B: N(implemented, undocumented) = 37 ; N(documented, unimplemented) = 6
PASS: set S is documented by `--help` in both dialects
```

- The task brief's "30 implemented-but-undocumented native flags + 3 documented-but-unimplemented" differs from the measured **33 / 4 (dialect A)**: the audit (`02` §2.2) counted 30 and did not include `-1`/`-2`/`-3` (numeric level aliases), which the measurement script does count; the "3" = `--ithr`/`--xparam`/`--silent` (4 spellings including `-s`). **The set S gate PASSes**, so this is a **registered debt**, not a regression.
- Stage 4 still has to do: `man/` and `docs/` synchronisation, and the deprecation policy (design §9.3).

### 7.6 The core layer's 0 changes — as a **boundary constraint**

`[measured · this report]` `git diff --name-only 24c069f..HEAD -- src/type src/xtb src/mctc src/param src/io src/disp src/coulomb` = **empty**; `grep -rln "use xtb_landscape" src/ | grep -v main.F90 | grep -v src/landscape/` = **empty**.

This means: **the suite's only integration points are `src/prog/main.F90` + `src/xhelp.f90` + `src/set_module.f90` + `src/setparam.f90` + `src/relaxation_engine.f90` (the 2-line fix of §5.2)**. The benefit is zero regression risk in the core layer; the cost is that **some defects can only be worked around in the suite layer, not fixed at the root** (e.g. the misleading `latticepoint` message of §5.4 can only be avoided by "not letting the error arise"; the same-class defect at `src/screening.f90:107` of `AGENTS.md` §6 item 5 is not fixed).

### 7.7 Other residuals (registered by the various reviews)

- **Text-class regressions have no automatic gate**: `landscape-help-coverage` only compares flag sets, not semantics (`54` INFO-2). That is exactly how HIGH-2 escaped validation.
- **The spherical (Lagrange) relaxation fix of IRC `gs`**: **intentionally deferred** (open item of PLAN §6.5); the current behaviour is pinned by the micro-cycle budget, measured cost 16.9/20.3 E/g vs 2 for HPC (a 3-4× gain).
- **MAPLE was not executed as a numerical reference**: the dispatcher is not importable (`ase`/`torch`/`fairchem` missing); only committed fixtures were read, comparing topology/frame count. **MAPLE numerical cross-validation must not be claimed** (`docs/xtb-landscape-M7-validation.md` §7.2).
- **IR/Raman intensities are 0** (the Hessian provider carries no dipole derivatives); the **index convention of `vib_normal_modes`** (rows = Cartesian components, columns = modes) is the opposite of the common reading.
- **`--lscan/--lts/--lirc/--lvib` do not reset the sentinels** (only the opt halves of `--lopt`/`--ohess` delete/create `.xtboptok`/`NOT_CONVERGED` at the start); a stale `.xtboptok` left by a previous `--lopt` may mislead an external pipeline that only looks at the sentinel (`51` §L1, still registered after handover to 3c).

---

## 8. Reproduction guide

### 8.1 Build (**`-O2` required**)

```bash
cd /home/owenchan/work_path/coding/fortran/xtb
cmake -S . -B /tmp/xtb-o2/build \
      -DWITH_TBLITE=OFF -DWITH_CPCMX=OFF -DWITH_TESTS=ON -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_Fortran_FLAGS_RELEASE="-O2 -DNDEBUG"
cmake --build /tmp/xtb-o2/build -j
```

- `WITH_CPCMX=OFF` is **required** (the toml-f target clash of cmake 4.3.2).
- `-O2` is **required**: CMake's `Release` is `-O3` for GNU, and under `-O3` gfortran 12.2.0 once produced a **wrong Hessian** for edcfbbe (water frequencies `-931.59/3329.78/3345.94` vs the correct `1557.41/3677.84/3684.01`); `-O3 -fno-tree-vectorize` is wrong too, the mechanism is unknown. `cmake/CMakeLists.txt` detects `-O3` and warns.

### 8.2 Tests (**`-R` filtering required**; a full `ctest` will time out)

```bash
export OMP_NUM_THREADS=1
# the whole landscape set (currently 147; ~45 s)
timeout 600 ctest --test-dir /tmp/xtb-o2/build -R landscape
# regression subset (12; <1 s)
timeout 600 ctest --test-dir /tmp/xtb-o2/build \
  -R 'Singlepoint|GFN|model-hessian|optimizer|Argparser|CAPI'
# duplicate-name check (the silent-drop trap of AGENTS §4)
grep -oE 'add_test\("[^"]+"' test/unit/CMakeLists.txt | sort | uniq -d    # expected: empty
# two-way `--help` consistency (set S difference = 0)
python3 test/unit/landscape_help_check.py src/prog/main.F90 src/xhelp.f90
```

⚠️ **Never run the full `ctest`** (84+ targets; `mctc-lib/ncoord` alone takes ~156 s, which exceeds a 30-minute agent timeout).
⚠️ **Always `OMP_NUM_THREADS=1` for scan/Hessian experiments** (measured 340× speedup: the same 2-point scan 0.13 s vs 44.12 s).

### 8.3 Re-running the key experiments

```bash
XTBPATH=$PWD /tmp/xtb-o2/build/xtb water.xyz --gfn2 --opt --lopt-engine lbfgs --lopt-level normal
# expected: stdout contains "[NOTICE] --opt now runs the landscape suite optimiser"
#       + "landscape: lbfgs minimisation"; writes .xtboptok / xtbopt.log / xtb.trj / xtbopt.xyz

# the native exit (the key of BC-4)
XTBPATH=$PWD /tmp/xtb-o2/build/xtb water.xyz --gfn2 --opt --method ancopt --norestart
# expected: stdout contains "ANC micro-cycles" and not "landscape:"; exit code 0

# the three guard rails
XTBPATH=$PWD /tmp/xtb-o2/build/xtb water.xyz --gfn2 --lopt --input <xcontrol with $opt>
# expected: exit 1, naming "--method ancopt"

# --ohess runs only once (no native Hessian/frequency block)
XTBPATH=$PWD /tmp/xtb-o2/build/xtb water.xyz --gfn2 --ohess --norestart
# expected: grep -c 'Numerical Hessian' = 0; grep -c 'analytical hessian' = 0;
#       grep -c 'landscape: harmonic frequencies' = 1

# the four axes of the level
XTBPATH=$PWD /tmp/xtb-o2/build/xtb water.xyz --gfn2 --lopt --lopt-level loose --norestart
# expected: maximum number of steps 200; energy tolerance 5.000000E-05; SCC accuracy 2.000000

# the upstream symbol-column defect (§5.2)
XTBPATH=$PWD /tmp/xtb-o2/build/xtb water.xyz --gfn2 --opt --method ancopt -I <$opt engine=lbfgs> --norestart
# expected: the symbol column of xtbopt.xyz = O/H/H (strict UTF-8, no control bytes)
```

### 8.4 The correct criterion for BC-4 (backward compatibility)

**Output files identical to round-off** + the stdout differences are **only** the registered non-determinism (version string `(<hash>|unknown-commit)`, date, random seed, `program call`, traceback addresses, **timing lines without a percent sign** `SCC setup/iterations … 0 min, X sec`). **"stdout is line-identical" is no longer claimed** (design §19.5).

The byte-level wording this section carried until 2026-10-16 is no longer true across builds, because
the fork adopted three upstream performance PRs (#1465 integrals, #1466 SCC linear algebra, #1464
dispersion gradient).  They remove O(N³) work from the SCC loop and are **numerically round-off only**:
measured on this tree, the four pinned H2O single-point energies stay bit-identical, full-precision
energies differ by 3e-15…5.6e-14 Eh, gradients by ≤1.4e-15 Eh/bohr, the `--hess` raw matrix by **0**,
and an `--opt` geometry by ≤1.5e-12 Å - but the SCC trace `ΔE` column, signed zeros and one Wiberg
ordering tie now differ, and `--molden` rotates degenerate orbitals.  The criterion above is therefore
the one to use; the BC-4 records below that were measured as `cmp`-identical were measured against the
then-current builds and remain valid as history.  See `AGENTS.md` §5 (exceptions #7-#9) and §6 for the
measurements.

```bash
git worktree add --detach /tmp/xtb-bc4/base 24c069f   # or any un-re-dispatched baseline
# build both sides, then cmp file by file for un-re-dispatched runtypes (--sp/--hess/--bhess/--md/--lscan/--lts/--lirc/--lvib)
```

---

## 9. Warnings for whoever comes next

1. **An injection experiment must first prove the mutation is semantically equivalent to the old logic**. A mutation that is "different but happens to be a correct alternative implementation" stays green, and in that case "green" proves nothing (`34-m4-coverage-fix-report.md`). This is the 4th lesson of the same kind in this direction.
2. **Confirm the binary and the path before verifying**. This direction tripped over: using an old binary/old build directory as evidence (`29` §105), mixing relative paths with `cd` (5 times), `test/unit/landscape_h2o.xyz` **not existing** in the source tree (it is generated by the `file(WRITE ...)` of `CMakeLists.txt:162-163`; the correct path is `<build>/test/unit/landscape_h2o.xyz`, `58` §6.3), and `head -3` catching the initial geometry rather than the final energy (`45` §4.1). **Always absolute paths.**
3. **An authoritative document can also contradict itself**: the design document `11-cli-redesign.md` was corrected at least **3 times** by implementation measurements (the field path §18.1, C7's blast radius §18.2, the PBC baseline sentence §18.5), plus §18.6 (`--ohess --ljson` writing two documents). **Corrections leave a trace, the history is not silently rewritten**: keep the original sentence + insert a correction pointer.
4. **`-O3` is not trustworthy, `-O2` is the supported level**. gfortran 12.2.0's `-O3` produced a wrong Hessian on a historical source state; the mechanism is unknown and the risk cannot be bounded. **Always use `-O2`.**
5. **`ctest` test names must be unique**: two `add_test` entries with the same name make CMake **silently drop** one (this once hid 19 scan unit cases). After adding tests, run `uniq -d`.
6. **A switch that "claims to take effect" must have a test that can falsify it**. The root cause of this direction's 7 dead knobs (§5.6) is always "parsed but read by nobody" or "unconditionally overwritten by a higher priority". The defence: provenance fields (`*_source`) + echo + an assertion that "changing the value must change the result".
7. **The core layer is a read-only boundary**. 0 changes to `src/type|xtb|mctc|param|io|disp|coulomb` is a hard constraint; when a core-layer defect is found, **work around it in the suite layer and register it** (§5.2's `screening.f90:107` is a known, unfixed instance of the same class).
8. **Never use intrinsic assignment between two `TMolecule` objects**, always `call a%copy(b)` (the stride defect of gfortran 12.2.0, `AGENTS.md` §6 item 5).
9. **The rc channel and the input file are two channels**: `$opt`/`$md`/… in `.xtbrc` are silent for the suite (§7.1), but `$scc`/`$thermo` **do** take effect. Do not treat "the rc is always ignored" as a fact.
10. **Do not run the full `ctest`**; do not run xtb without a `timeout`; do not `find /`.

---

## Appendix A — numbers I measured and verified myself in this report

| Command | Result |
|---|---|
| `git log --oneline` / `git status --short` | 9 direction commits + 1 meta-commit; the work tree is **clean**; `git diff --cached` empty |
| `git diff --shortstat 24c069f..HEAD -- src/ test/ cmake/ man/` | 38 files, +9018 / −483 |
| `ctest --test-dir /tmp/xtb-3cf/build -N -R landscape` | **Total Tests: 147** |
| `ctest --test-dir /tmp/xtb-3cf/build -R landscape` (`OMP_NUM_THREADS=1`) | **147/147 passed** (44.36 s) |
| `ctest -R 'Singlepoint\|GFN\|model-hessian\|optimizer\|Argparser\|CAPI'` | **12/12 passed** (0.29 s) |
| `grep -oE 'add_test\("[^"]+"' test/unit/CMakeLists.txt \| sort \| uniq -d` | **empty** (no duplicates) |
| `python3 test/unit/landscape_help_check.py src/prog/main.F90 src/xhelp.f90` | S(main)=S(help)=135, difference 0; N(impl,undoc)=33 (A)/37 (B); N(doc,unimpl)=4 (A)/6 (B); **PASS** |
| `git diff --name-only 24c069f..HEAD -- src/type src/xtb src/mctc src/param src/io src/disp src/coulomb` | **empty** (core layer 0 changes) |
| `grep -rln "use xtb_landscape" src/ \| grep -v main.F90 \| grep -v src/landscape/` | **empty** |
| `grep -c 'add_test("xtb/landscape' test/unit/CMakeLists.txt` | 127; plus the 6 `--ljson` tests expanded by `foreach(task IN ITEMS opt scan scan2d vib ts irc)`, plus the **14** `landscape-*` test-drive suites in the `tests` list = **147** (consistent with `ctest -N -R landscape`; all names start with `xtb/landscape-`) |

## Appendix B — items I **cannot** verify and only faithfully relay

| Item | Why not verified | Relay source |
|---|---|---|
| the historical test counts of each stage (101/109/115/118/120/121/125/127/136/138/143/145/146) | the intermediate work trees corresponding to those build directories **no longer exist** (only the final state `/tmp/xtb-3cf` can be run) | the stage reports (see the report numbers in the "acceptance evidence" column of §3 for each stage) |
| the "red-before / green-after" output of all injection experiments | would require re-mutating + rebuilding + restoring, beyond the scope of this report | the stage reports (`06`/`21`/`32`/`45`/`53`/`55`/`56`/`58`) |
| all M2..M7 acceptance numbers (dfreq, RMSD, MEP deviation, LEPS barrier, etc.) | would require re-running the corresponding manual experiments | `docs/xtb-landscape-M7-validation.md`, `AGENTS.md` §2, `58` §4 |
| the three re-baseline deviations of 3c that were not re-run (§7.4) | the report itself marks them "not measured / static judgement" | `58-stage3c-closeout.md` §4 |
| the **exact count** of "the 13th/14th time a fact was corrected by measurement in this direction (6/7 times the parent session)" | the reports give only "the 13th, 6 times the parent session" for stage 1 and "the 14th" for 3c; there is **no** unified "7 times" count | `21-stage1-fix-report.md:81`, `58-stage3c-closeout.md:205` |
| the counting convention of "5 dead knobs" | the reports give "4" (design §3.5) / "three earlier" (`48`); I list **7** by verifiable source | the table and sources of §5.6 |
| whether the meson-side tests are executable | no meson/ninja on this machine | `25` §7, `39` §5, `54` §5 |
| MAPLE numerical cross-validation | the dispatcher is not importable | `docs/xtb-landscape-M7-validation.md` §7.2 |

---

*This report only adds one file, `docs/xtb-cli-redesign-delivery.md`; it does not `git add`, does not commit, and does not change `src/`, `test/`, or `AGENTS.md`.*
