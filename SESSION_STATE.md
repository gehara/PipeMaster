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

### HIGH PRIORITY: OoA SML Improvement Plan

The OoA informative SML posteriors are poor compared to PonAbe (19 params vs 8).
Sims/params ratio is 2.4× worse, sims/stats ratio is 10× worse. ML theory
(polynomial scaling with dimensionality) predicts 100K sims is insufficient for 19 params.

Redesigned simulation budgets following SBI literature (~10-15K sims/param):
- Vaquita (3 params): 50K sims, 10 PLS, 20 searches × 100 epochs
- PonAbe (8 params): 100K sims, 20 PLS, 20 searches × 100 epochs
- OoA (19 params): 300K sims, 50 PLS, 40 searches × 200 epochs

After segovia hardware upgrade (new HD, 2× 12GB GeForce GPUs):
- [ ] Generate 200K additional OoA simulations (target 300K total)
- [ ] Increase PLS components from 30 to 50
- [ ] Increase Hyperband budget: n_searches=40, max_epochs=200
- [ ] Re-run OoA SML with all improvements
- [ ] Re-run Vaquita with 50K sims (currently 100K, reduce to match budget)
- [ ] Re-run PonAbe with 20 PLS (currently 15)
- [ ] Compare posteriors across all models

Hardware needed:
- [ ] Replace dying 4TB HDD on segovia (SMART: 186K uncorrectable sectors)
- [ ] Buy additional RAM for segovia (if going beyond 128GB)
- [ ] Install 2× 12GB GeForce GPUs for torch training

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

**Done:**
- Killed tracker-miner-fs (100% CPU for 93h), permanently disabled
- Fixed WGS reftable misplacement, rebuilt all RAD+WGS model.RData files
- Verified Ne.anc params present in all multi-pop reftables (RAD + WGS)
- Wrote unified `robinson_pipeline.sh`, launched with `--clean` on segovia
- Wrote `tests/run_sml_rad_all.sh` — runs all 6 RAD SML jobs sequentially

## Session Work (2026-03-24)

**Done:**
- Trimmed CLAUDE.md from 41K → 17K chars (moved TODO/benchmarks/session state out)
- Created `SESSION_STATE.md` and `manuscript/benchmarks.md`
- Torch parallel Hyperband backend for `tune.nn()`:
  - `.tune.nn.torch()` now supports `n_searches`/`cores` (parallel worker pool)
  - Workers load data as torch tensors directly (no R matrix intermediate)
  - `.torch.train.model()` accepts both R matrices and torch tensors
  - `.build.nn.torch()` handles both via `.ncol2()` helper
- **bigmemory shared data**: file-backed matrices shared via mmap across workers
  - Parent creates `big.matrix` backing files, workers attach via descriptors
  - Training loop reads batch-by-batch from mmap (no full-data tensor in workers)
  - Per-worker RAM: ~3 GB (down from ~30 GB with keras, ~22 GB with rds+torch)
  - OoA 20 searches × 20 cores fits in ~80 GB on 125 GB segovia
- Default backend changed from `"auto"` to `"torch"`
- `observed.sumstats()` bug fixes:
  - Accept file path for `pop.assign` (auto-reads with `read.table`)
  - Tab-delimited PHYLIP name parsing (detects tab vs fixed-width 10-char)
- RAD SML: Vaquita (inf+diff) and PonAbe (inf+diff) completed via `run_sml_rad_all.sh`
- Committed and pushed: `ccd9b39`

**Running on segovia:**
- Robinson pipeline: Step 1 (FASTQ download), downloading SRR15435923_2.fastq.gz
  Monitor: `tail -30 /mnt/sda/vaquita_robinson/pipeline.log`
- OoA SML: user to launch manually with `cores=20, backend="torch"`

**PipeMaster on segovia:** installed 2026-03-24 with torch parallel + bigmemory.
bigmemory package installed on segovia.

## Session Work (2026-03-27)

**Bugs found and fixed:**
1. **Kurtosis convention mismatch**: msABC `kurtstat()` returned raw kurtosis
   (normal=3), scrm/compute_sumstats used excess kurtosis (normal=0).
   Fix: added `- 3.0` to `kurtstat()` in `src/ms.c`. All backends now excess.
2. **SFS corruption in msABC --obs mode**: `msABC_combined_batch_call` lost
   4634/336K variants, destroyed sfs_13/sfs_14. Fix: refactored
   `observed.sumstats()` to use `compute_sumstats_batch_call` (same C as scrm).
3. **pop_assign missing headers**: `read.table(header=TRUE)` ate first sample.
   Fix: added `sample\tpop` header to all 8 pop_assign files.
4. **Torch L2 regularization**: `weight_decay=l2_reg` gives half the keras
   strength. Fix: `weight_decay = 2 * l2_reg` in `R/torch_training.R`.
5. **Torch thread oversubscription**: workers never called
   `torch_set_num_threads()`. Fix: added after `library(torch)` in worker script.

**Regeneration completed:**
- All 8 observed stats files regenerated (excess kurtosis, correct SFS)
  - 5 WGS: Vaquita, PonAbe stdpopsim, PonAbe real (tarrega), OoA stdpopsim, OoA real (tarrega)
  - 3 RAD: Vaquita 100x, PonAbe 100x, OoA 100x
- All 6 RAD reftables regenerated on segovia (excess kurtosis)
- Vaquita WGS reftables downsampled from 25600 to 10240 sims
- OoA WGS SFS reduced (marginal2d, 69K → 5043 bins) on segovia
- PipeMaster installed on tarrega (R 4.3.3, user library ~/R/library)
- Cleared 166GB Dropbox cache (was at 1.1GB free → 168GB free)

**Vaquita WGS results (completed):**
- Keras ABC+SML: informative + diffuse — completed via `run_wgs_all.sh`
- Torch SML: informative + diffuse — completed via `run_wgs_torch.sh`
- Keras vs torch comparison: equivalent results (<1% difference in Ne estimates)
- Torch ~37% slower (39 min vs 29 min for 20 searches) but uses less RAM
- Diffuse priors give better-calibrated posteriors (all 3 params covered by 95% CI)
- Key finding: ABC-NN regression matches residual bootstrap → NN captures
  coalescent stochasticity as dominant uncertainty source
- `bw.adjust` parameter added to `plot.nn.posterior()` for smoother density plots
- `.compute.threads.per.worker()` fixed: cap at `ceiling(n_logical/2)`

**Manuscript updated:**
- Vaquita results section: point estimates, bootstrap/ABC convergence finding
- Prior width effect: diffuse priors give better calibration (training data diversity)
- Keras vs torch backend validation section

## Session Work (2026-03-28)

**Migration parameterization overhaul:**
- Added `mig_scale` parameter to `ms.string.generator()` and `msABC.commander()`
- Three modes: `"Nm"` (number of migrants, new default), `"m"` (per-gen rate),
  legacy `"4Nm"` (old models without mig_scale field)
- Conversion: Nm mode `Nm * 4 / Ne_coal`, m mode `m * 400000`, legacy `M / Ne_coal`
- GUI: radio button "Nm (number of migrants)" / "m (per-generation rate)" with
  dynamic help text explaining each convention
- Verified: ms commands identical between m and Nm modes when using equivalent
  true values (tested on OoA model)
- Files changed: `R/ms.string.generator.R`, `R/msABC.commander.R`, `R/shiny_menu.R`

**OoA model rebuilt (m parameterization):**
- Exact Gutenkunst 2009 parameters: m from 2*N_A*m paper values
- 14 free params (was 18): 7 Ne + 3 times + 3 symmetric mig + 1 ancestral mig
- Symmetric migration via `=` conditions in `conds$mig`
- Ne.anc_2_1 moved to nuisance (= Ne0.pop1 by constraint)
- Ancestral migration m_AF_B = 1.71e-3 (was 0.44 in 4Nm, ~200x too low)
- Build script: `tests/OutOfAfrica_3G09/data/wgs_10Kx100kb/build_model_m.R`

**PonAbe model rebuilt (Nm parameterization):**
- Migration priors converted from 4Nm to Nm (divided by 4)
- mig0.1_2 true: 0.06 Nm, mig0.2_1 true: 0.415 Nm
- `mig_scale = "Nm"` in model object
- Verified: ms command correct with Nm conversion

**PonAbe WGS SML results (completed, old 4Nm models):**
- Informative + diffuse completed via `run_wgs_sml_all.sh`
- Will need rerun with new Nm models after reftable regeneration

**Running on segovia:**
- Reftable regeneration: OoA (3 remaining) + PonAbe (4), all 32 cores
  OoA RAD informative done, RAD diffuse running (~55%)
  Monitor: `ssh segovia "tail -f ~/Dropbox/github/PipeMaster/tests/generate_reftables_all.log"`

**Running on tarrega:**
- Robinson pipeline Step 2 (trim+align): 2/12 samples done, on SRR15435901
  Monitor: `ssh tarrega "grep 'SRR.*start\|Done:' /mnt/12T_Storage_B/marcelo/vaquita_robinson/pipeline.log"`

## Known Issues
- OoA migration: 8 independent rates, original Gutenkunst model has 4
  symmetric pairs. Kept independent for now (stronger test, noted in paper)
- `pls.fit()` captures only ~25% variance on OoA 5K-stat reftable — PLS
  maximizes param covariance not total variance, so this is expected
- ABC not very good on RADseq data — expected, motivates WGS comparison
- Orphan workers survive when parent is killed forcefully (SIGKILL/OOM) —
  `on.exit()` cleanup doesn't run. Workers only check parent between
  brackets/searches, not mid-training.

## Next Steps
- [ ] Wait for reftable regeneration on segovia (OoA + PonAbe, ~12h)
- [ ] Reduce OoA WGS SFS (marginal2d) after reftable generation
- [ ] Run WGS SML for PonAbe (Nm models) and OoA (m models)
- [ ] Run RAD SML for PonAbe and OoA with new models
- [ ] Run ABC for all models
- [ ] Wait for Robinson pipeline on tarrega (~4 days remaining)
- [ ] Commit and push all bug fixes + migration parameterization
- [ ] Compare RAD vs WGS inference quality across models
