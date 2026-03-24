# Session State & TODO

## TODO

### Manuscript Datasets (3 models × 3 datasets)

For each model, 3 datasets: pseudo-RAD (stdpopsim, 10K×100bp), pseudo-WGS
(stdpopsim, 10K×100kb), and real WGS (published data, real fragment sizes).

**Vaquita2Epoch** (1-pop bottleneck, 3 params):
- Real data: Robinson et al. 2022 (SRA PRJNA751981, 15 unrelated vaquitas)
- Pipeline: `tests/Vaquita2Epoch/data/wgs_10Kx100kb/robinson_pipeline.sh`
- Status: pipeline step 3 complete, 5/15 ENA samples + 9 SRA BAMs, 10 ENA failed

**PonAbe TwoSpecies_2L11** (2-pop IM):
- Real data: Han et al. 2025 great ape genomes (PHAIDRA)
- 10 P. pygmaeus (pop1) + 10 P. abelii (pop2)
- VCF ready: `/mnt/sda/ponabe_realdata/pipemaster/ponabe_neutral_loci.vcf`
- 4.65M variants, 206K neutral fragments (median 1.4kb, 0.44 Gb)

**OutOfAfrica_3G09** (3-pop OoA):
- Real data: 1000 Genomes NYGC 30x (Byrska-Bishop et al. 2022)
- 20 YRI (pop1) + 20 CEU (pop2) + 20 CHB (pop3)
- VCF ready: `/mnt/sda/ooa_realdata/pipemaster/ooa_neutral_loci.vcf`
- 2.09M variants, 206K neutral fragments (same hg38 masks as PonAbe)

Neutral regions (PonAbe & OoA share hg38 masks: CDS±10kb, repeats, CpG):
206,079 fragments, 0.44 Gb, 99.6% are 1-5 kb, 35 fragments >100kb.

Priors are standardized per parameter class within each model (same range
for all Ne0, all Ne1, all mig, all time). Different models have different
prior ranges appropriate to their demographic history.

Existing data inventory (`tests/` directory):

| | Pseudo-RAD (10K×100bp) | Pseudo-WGS (10K×100kb) | Real WGS |
|---|---|---|---|
| **Vaquita** | | | |
| model | `model.RData` ✅ | `model.RData` ✅ | — |
| observed | `obs_allstats_100x_mean.txt` ✅ | `obs_allstats.txt` ✅ | waiting on ENA download |
| reftable | `SIMS_reftable_100K.txt` 100K ✅ | `SIMS_reftable_10K.txt` 10K ✅ | — |
| **PonAbe** | | | |
| model | `model.RData` ✅ | `model.RData` ✅ | `model.RData` ✅ |
| observed | `obs_allstats_100x_mean.txt` ✅ | `obs_allstats_wgs.txt` ✅ | `obs_allstats_real_wgs.txt` ✅ (stringent) |
| reftable | `SIMS_reftable_100K.txt` 100K ✅ | `SIMS_reftable_10K.txt` 10K ✅ | — |
| **OoA** | | | |
| model | `model.RData` ✅ (3-pop, standardized) | `model.RData` ✅ (3-pop, standardized) | `model.RData` ✅ |
| observed | pseudoobs PHYLIP exist | `obs_allstats_wgs.txt` ✅ | `obs_allstats_real_wgs.txt` ✅ |
| reftable | `SIMS_reftable_100K.txt` 102K ✅ (3-pop) | `SIMS_reftable_10K.txt` 10K ✅ | — |

All reftables (RAD + WGS) include Ne.anc params for multi-pop models.

Remaining TODO:
- [ ] Vaquita real WGS — pipeline step 3 done, 5/15 ENA + 9 SRA BAMs, 10 ENA failed
- [ ] Real WGS reftables for all 3 models (scrm_stats_multi_call ready)
- [ ] ABC/SML runs for PonAbe WGS, OoA RAD, OoA WGS, all real WGS

### CNN Posterior Estimation

Completed:
- Quantile regression — pinball loss, well-calibrated, matches ABC posteriors
- MC Dropout — too narrow, only captures model uncertainty, not recommended
- Bootstrap retraining — too narrow, only captures training variability
- Conformal prediction — guaranteed coverage, narrower CIs than ABC, no retraining
- WGAN-GP — conditional generator produces full posterior samples, @tf.function compiled

Remaining:
- [ ] WGAN-GP production test — rerun vaquita with 20 configs, 2000 max epochs
- [ ] Mixture Density Network (MDN) — Gaussian mixture params, skewed/multimodal
- [ ] Deep Ensembles — N models with different seeds, prediction spread = uncertainty
- [ ] Neural Posterior Estimation (NPE) — normalizing flows, state-of-the-art
- [ ] Evidential Deep Learning — Normal-Inverse-Gamma output, single forward pass

### Simulation-Based Inference

Two approaches for gradient-based parameter optimization:

**BOLFI** (Bayesian Optimization for Likelihood-Free Inference):
- GP surrogate for discrepancy, LCB acquisition, MCMC on GP pseudo-likelihood
- Very sample-efficient (500-2000 sims) but scales poorly above ~10-15 params
- New deps: DiceKriging/hetGP, DEoptim. Status: not started.

**Neural Emulator + backprop** (partially implemented):
- NN forward model θ→S, backprop loss=||f(θ)-S_obs||² to update θ
- Scales to high dimensions, reuses existing infrastructure
- Status: core pipeline implemented (train/optimize/ABC/MCMC)

### Other TODO Items

- [ ] Verify monomorphic class (bin 0) in simulated and observed SFS
- [ ] Remove scrm dependency (replace with built-in ms in sim.stacked.segsites)
- [ ] Check PonAbe reftable — verify 10,240 rows, inspect stats/SFS
- [ ] Convergence validation: compare point estimates across 10 search models
- [ ] `emulator.MCMC()`: add `start = "random"` option — LHS over prior
  bounds in log-space for better prior coverage when posterior is multimodal.
  Needs: `.lhs.sample()` helper, input validation, new branch in starting
  points section. See `add_random_start_to_emulator_MCMC.md` for details.
  Applies to both keras (`R/emulator_posterior.R`) and torch
  (`R/torch_posterior.R`) versions.

### hABC Workflow (Ubelmania)

- [ ] ABC reference tables (102K sims per species, 32 cores) — running on segovia
- [ ] Run ABC per species (posterior estimation)
- [ ] Build hModel for hierarchical analysis
- [ ] Hierarchical ABC analysis

### Port TF/Keras to R Torch

- [ ] `R/tune.nn.R` — model architectures, Hyperband, training loop
- [ ] `R/emulator.R` — train.emulator, emulator.optimize, emulator.ABC, Jacobian
- [ ] `R/emulator_posterior.R` — emulator.MCMC Python loop → pure R+torch
- [ ] `R/OOD.diagnose.R` — model inference calls

(Torch rewrite files already exist — see Emulator Pipeline (Torch Rewrite)
section in CLAUDE.md.)

### Papers & Shipping Plan

3 papers planned. Each has its own scope for what gets shipped vs gitignored.

**Paper 1 — PipeMaster v2 (shipped package, CRAN):**
- Simulation: `sim.all.stats()`, `sim.scrm.stats()`, `sim.sumstat()`, etc.
- Observation: `observed.sumstats()`, `observed.pw.distances()`
- SML: `tune.nn()`, `nn.predict()`, `save/load.tune.result()` (keras + torch backends)
- ABC: `abc.rejection()` (standalone, extracted from emulator.R)
- OOD: `OOD.diagnose()` (adapted for both keras and torch models)
- Model building: GUI, `PlotModel()`, all other existing exports
- DESCRIPTION: `torch`, `keras`, `tensorflow`, `reticulate` in Suggests

**Paper 2 — Neural emulator (not shipped, gitignored):**
- `R/emulator.R` — `train.emulator()`, `emulator.optimize()`, `emulator.ABC()`,
  `emulator.sensitivity()`, `emulator.diagnose.stats()`
- `R/emulator_posterior.R` — `emulator.MCMC()`, S3 methods
- `R/torch_emulator.R`, `R/torch_posterior.R`, `R/torch_OOD.R`

**Paper 3 — hABC (not shipped, gitignored):**
- `R/grid.search.R` — `grid.search()`, `update.priors.from.grid()`
- `sim.sumstat.lhs()`

### Torch Integration into tune.nn (in progress, 2026-03-13)

Changes needed to ship torch backend in PipeMaster v2:

1. `tune.nn()` — add `backend = c("auto", "torch", "keras")` param
2. `nn.predict()` — detect model type (keras vs torch), dispatch
3. `OOD.diagnose()` — Check 4 model disagreement: torch dispatch
4. `save/load.tune.result()` — handle torch model serialization
5. Extract `abc.rejection()` from `emulator.R` → `R/abc.rejection.R`
6. Track: `R/torch_modules.R`, `R/torch_training.R`, `R/abc.rejection.R`
7. Gitignore: `R/grid.search.R` (add to existing gitignore list)
8. NAMESPACE: add `abc.rejection`, `OOD.diagnose` exports

---

# Session State (2026-03-23)

## Commits pushed (through 2026-03-16)
- `4f0d1dd` — Remove old tutorials, binaries, images, deprecated files
- `9291b0f` — Add vendored scrm 1.7.5 C++ backend
- `f6b07a5` — Update R source: Ne.anc, stat names, tune.nn torch
- `0470d13` — GUI overhaul: model plot, migration UI, visualization, logo
- `3a5013d` — Condition list system, sample.w.cond rewrite, GUI and model fixes
- `05501bc` — GUI plot overhaul, condition cleanup, dead code removal

## Uncommitted Code Changes (2026-03-17 through 2026-03-19)

**`R/tune.nn.R`** — Major refactor:
- `nn.predict()` rewritten: 3 methods (`"point"`, `"bootstrap"`, `"ABC_NN_regression"`)
  - `"bootstrap"` — locally-weighted residual bootstrap (Epanechnikov kernel)
  - `"ABC_NN_regression"` — Beaumont et al. 2002 regression adjustment
  - Added `pls = FALSE, n.pls = 20L` params for PLS distance projection
  - Both methods: no retraining, single forward pass
- Removed conformal, MC dropout, retrain-bootstrap (~1,150 lines removed)
- Removed `mc_dropout` from `.build.nn()`, `.build.resnet()`, `.build.cnn1d()`, `.build.cnn2d()`
- S3 methods updated for new result structure (`$bootstrap`, `$abc`)

**`R/OOD.diagnose.R`**:
- Per-stat percentile check: 3-way classification (`"ok"`, `"outlier"`, `"uninformative"`)
- Uninformative = zero-variance bins (prevents false flags from constant SFS entries)
- Plot: grey bars for uninformative, legend with 3 categories

**`R/reduce.sfs.R`** — Unified function:
- Single `reduce.sfs()` accepts R objects OR file paths
- File path: batched I/O via C reader (`src/read_tsv.c`), 19x faster than `read.table()`
- R object: in-memory with 2GB size check, progress bar
- Methods: `"marginal2d"` (pairwise 2D marginals) and `"binned"` (coarsened nD)
- Removed `reduce.sfs.file()` — merged into `reduce.sfs()`

**`R/pls.transform.R`**:
- `pls.fit()`: added `max.rows = 10000L` — subsamples for NIPALS fitting,
  centering/scaling from full data. Prevents OOM on large reftables.
- Exported: `pls.fit`, `pls.project`

**`R/parameter_samplers.R`** — Orphan worker prevention:
- `.pm.parent.alive(pid_file)` — cross-platform parent check via `tools::pskill(pid, 0)`
- `.pm.register.parent(pid_file, worker_pids_env)` — writes PID file, registers
  `on.exit()` cleanup that kills all workers

**`R/sim.sumstat.msABC.R`** + **`R/sim.scrm.R`**:
- Added `.parent.pid.file` parameter
- Parent writes `.PM_parent.pid`, tracks worker PIDs, kills on exit
- Workers check parent alive before each block

**`R/shiny_menu.R`** — GUI plot fixes:
- Migration arrow direction: swapped `from`/`to` in ms parser to match
  biological direction (source → destination). ms `-m i j` = into i from j.
- Migration rate display: converted from ms-internal M to 4Nm scale using
  receiving pop's Ne (present-day for `-m`, time-appropriate for `-em`)
- Migration label lookup: adjusted for swapped direction
- Time labels: join names get priority over tied Ne-change names
- ms command window: shows conditions alongside sampled parameters

**`R/torch_modules.R`**: Removed `mc_dropout` from `.build.nn.torch()`,
removed `.torch.predict.mc()`

**`R/torch_training.R`**: Removed conformal/bootstrap worker scripts

**`src/read_tsv.c`** (NEW): Fast C file reader for tab-separated numeric data.
Parses doubles with `strtod`, skips unrequested columns, 19x faster than
`read.table()`. Registered as `read_tsv_call` in `src/init.c`.

**NAMESPACE**: `reduce.sfs`, `pls.fit`, `pls.project` exported.
Removed `reduce.sfs.file`.

**Manuscript** (`pipemaster_v2_draft.Rmd`):
- Updated `nn.predict()` methods section (residual bootstrap + ABC-NN regression)
- Added prior specification rules with formulas
- Updated OOD diagnostics (3-way classification)

## Test Folder Reorganization (2026-03-18)
```
tests/<Model>/data/<seqtype>/
  observed/      # obs stats, pseudoobs, VCF (shared)
  informative/   # model.RData, reftable (true/3 to true*3)
  diffuse/       # model.RData, reftable (true/10 to true*10)
tests/<Model>/{ABC,SML}/<seqtype>/
  informative/   # results
  diffuse/       # results (existing data moved here)
```

## Prior Specification Rules
- Ne/time: informative [true/3, true×3], diffuse [true/10, true×10]
- Migration (4Nm): informative [0, max(true×2.5, 1)], diffuse [0, max(true×5, 2)]
- Minimum migration upper bound: 1 Nm (informative), 2 Nm (diffuse)
- Build scripts: parametric `{informative|diffuse}` argument
- All 6 model.RData files built (3 models × 2 priors)
- Models match GUI-compatible structure: `labels`, `conds` (list format),
  `sum_anc_ne`, `tree` with `;`

## RAD Reftables (10K × 100bp) — COMPLETE
All 6 reftables generated on segovia (102K sims each, 32 cores):
- Vaquita: informative 38MB, diffuse 39MB
- PonAbe: informative 219MB, diffuse 227MB
- OoA: informative 14GB (gzipped 419MB), diffuse 14GB (gzipped 649MB)

OoA SFS reductions (on segovia):
- `SIMS_reftable_reduced.txt`: marginal2d, 5043 SFS cols, ~1.2GB each
- `SIMS_reftable_binned.txt`: binned bin_size=4, 1331 SFS cols, ~457MB each

## ABC Results (RAD, informative priors)
- **Vaquita**: done, results in `ABC/rad_10Kx100bp/informative/`
- **PonAbe**: done, results in `ABC/rad_10Kx100bp/informative/`
- **OoA marginal2d**: done, results in `ABC/rad_10Kx100bp/informative/`
- **OoA binned**: running comparison test
- ABC not very good on RADseq — expected (fewer seg sites, sparser SFS)

## SML Status (RAD)
- **Vaquita informative**: done
- **All others**: ready to run via `tests/run_sml_rad_all.sh`
- PonAbe diffuse had stale `tune_result` (7 params, missing Ne.anc_1_2 from
  pre-Ne.anc reftable) — must delete before retraining
- All scripts use full 100K reftable (no subsampling)
- Hyperband params (all models): `max_epochs=100, eta=3, n_searches=20,
  top_k=3, cores=20, gpus=0, backend="keras"`
- PLS in nn.predict: Vaquita `pls=FALSE`, PonAbe `n.pls=15`, OoA `n.pls=30`

## OoA SFS Reduction Comparison
- **marginal2d** (5043 bins): better results, needs PLS for distances
- **binned bin_size=4** (1331 bins): worse results, loses low-frequency resolution
- Both reduction scripts in `tests/OutOfAfrica_3G09/data/rad_10Kx100bp/`:
  `reduce_sfs.R` (marginal2d) and `reduce_sfs_binned.R` (binned)
- OoA ABC script uses marginal2d + PLS inside abc.rejection()

## Vaquita Real WGS Pipeline (segovia) — RESTARTED
- **Full restart launched** (2026-03-23): `robinson_pipeline.sh --clean`
  - Unified script: aria2c download (8 connections) + `gzip -t` verification
  - Re-downloads corrupt FASTQs automatically (up to 3 retries)
  - `samtools quickcheck` on sorted BAMs before MarkDuplicates
  - Root cause of prior failures: ENA FTP downloads without gzip verification
    left corrupt FASTQs that looked OK (right size, wrong content)
  - Script: `tests/Vaquita2Epoch/data/real_wgs/robinson_pipeline.sh`
  - Log: `/mnt/sda/vaquita_robinson/pipeline.log`
  - BAMs cleaned (freed 273GB → 485GB free), now verifying+redownloading FASTQs
  - Pipeline uses 8 threads (bwa/samtools/bbduk), SML uses 20 cores = 28 total

## WGS Reftable Generation (segovia) — COMPLETE
- All 6 reftables done (2026-03-21), 10,240 sims each, 20 cores
- **Bug found+fixed**: `source()` in WGS `build_model.R` clobbered `out_dir`,
  writing reftables to RAD dirs. Fix: removed nested `source()`, RAD models
  must be pre-built. Reftables moved to correct WGS dirs.
- Launcher: `tests/generate_wgs_reftables.sh`
- Timings: Vaquita 13min+1h, PonAbe 2.7h+12.3h, OoA 5.9h+22.5h (inf+diff)

| Model | Informative | Diffuse |
|---|---|---|
| Vaquita (71 cols) | 9.4 MB | 9.5 MB |
| PonAbe (558 cols) | 35 MB | 36 MB |
| OoA (69,112 cols) | 1.7 GB | 1.8 GB |

## Orphan Worker Prevention
- Added to `sim.sumstats()`, `sim.scrm.sumstats()`, and `.launch.rscript.pool()`
- Parent writes `.PM_parent.pid`, workers check with `tools::pskill(pid, 0)`
- `on.exit()` in parent kills all workers + removes PID file
- Hyperband workers check parent before search and before retraining
- Fixed empty `.done` file crash (length-0 exit_code → NA_integer_)

## Session Work (2026-03-23)

**Done this session:**
- Killed tracker-miner-fs (100% CPU for 93h), permanently disabled
- Fixed WGS reftable misplacement: `source()` in WGS `build_model.R` clobbered
  `out_dir`, writing reftables to RAD dirs. Fix: removed nested `source()`.
  Moved 6 reftables to correct WGS dirs, rebuilt all RAD+WGS model.RData files.
- Verified Ne.anc params present in all multi-pop reftables (RAD + WGS)
- Wrote unified `robinson_pipeline.sh` — aria2c + gzip verification + full
  pipeline in one script. Launched with `--clean` on segovia (Step 1 in progress,
  finding and fixing corrupt FASTQs)
- Wrote `tests/run_sml_rad_all.sh` — runs all 6 RAD SML jobs sequentially
- Updated all SML Hyperband params to `max_epochs=100, n_searches=20, cores=20`
- PonAbe diffuse has stale `tune_result` (7 params, missing Ne.anc_1_2) — must
  delete before running: `rm -rf tests/PonAbe_TwoSpecies/SML/rad_10Kx100bp/diffuse/tune_result`

**Running on segovia:**
- Robinson pipeline: Step 1 (FASTQ verify+download), 8/24 SRRs done, 6 corrupt
  files re-downloaded successfully. Monitor: `tail -30 /mnt/sda/vaquita_robinson/pipeline.log`
- SML RAD: **not yet launched** — user will launch manually

**PipeMaster on segovia:** installed 2026-03-19 (has nn.predict rewrite, Ne.anc,
all uncommitted changes through that date). `Date:` field in DESCRIPTION is
stale (2026-02-15) but `Built:` timestamp is 2026-03-19.

## Disk Space (local)
- 143GB free (Dropbox unsync of projects_finished complete)
- OoA reftables gzipped on segovia (14GB → ~500MB each)
- GPU: RTX A2000 4GB — does not fit OoA SML (5K inputs), all SML runs use CPU

## Known Issues
- OoA migration: 8 independent rates, original Gutenkunst model has 4
  symmetric pairs. Kept independent for now (stronger test, noted in paper)
- `pls.fit()` captures only ~25% variance on OoA 5K-stat reftable — PLS
  maximizes param covariance not total variance, so this is expected
- ABC not very good on RADseq data — expected, motivates WGS comparison

## Next Steps
- [ ] Launch SML RAD: `tests/run_sml_rad_all.sh` (delete PonAbe diffuse tune_result first)
- [ ] Wait for Vaquita pipeline (Step 1 downloading, then ~6 days alignment)
- [ ] Run ABC/SML on WGS reftables (all 6 available)
- [ ] OoA WGS: reduce SFS (69K cols) before ABC/SML
- [ ] Run diffuse prior ABC after informative results reviewed
- [ ] Commit all uncommitted changes
