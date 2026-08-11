# Edit list: harmonise `ms.c` accumulation with `compute_sumstats.c` / `scrm_stats.cpp`

**Status: PROPOSED — nothing applied.** For review before touching the C source.

## Background

All three summary-statistic paths call the **same per-locus kernels** from
`data_sumstat.c` (`theta_pi`, `theta_w`, `tajD`, `ZnS`, `dvstat`,
`calculations`, `Fst_*`). The statistics themselves are computed by identical
code. The paths diverge **only in the accumulation layer**:

| Path | Uninformative locus | Denominator |
|---|---|---|
| `ms.c` (msABC) | dropped from that statistic's denominator | per-statistic / per-pop / per-pair |
| `scrm_stats.cpp:98` | `if (nsegsites == 0) return;` -> contributes 0 | uniform `nloci` |
| `compute_sumstats.c:461` | `if (isNull(hapmat)...) continue;` -> contributes 0 | uniform `nloci` |

msABC computes the mean over *informative* loci; the other two compute the
mean over *all* loci scoring uninformative ones as zero. These are different
estimators and they diverge by exactly the informative fraction.

Invisible at WGS scale (nearly every 100 kb locus is informative for every
statistic). Severe at RAD scale: with mean `segs` ~ 0.37 per 100 bp locus,
~69% of loci are monomorphic and only a small minority carry >= 2 segregating
sites.

Measured against simulations at the true parameters (PonAbe RAD, informative):

| Stat | obs/sim | z |
|---|---|---|
| segs, pi, thetaW | 0.995-1.008 | < 0.06 |
| Hd | 1.001-1.011 | < 0.06 |
| Fst | 0.552 | -2.63 |
| nhap | 0.669-0.752 | -4.6 to -6.9 |
| ZnS | 0.058-0.19 | -5.3 to -6.1 |

Statistics well-defined as zero on a monomorphic locus (segs, pi, thetaW, Hd)
agree. Every statistic that is *undefined* there diverges.

## Design decision

Harmonise **onto the `compute_sumstats.c` convention** (uniform denominator,
uninformative locus contributes 0), i.e. change `ms.c` only.

The alternative -- adopting msABC's per-statistic denominators in
`compute_sumstats.c`/`scrm_stats.cpp` -- would invalidate every observed file
and every WGS reference table, and it is the worse convention: a conditional
mean discards the informative-fraction signal, and makes a statistic's meaning
depend on how many loci happened to be informative, which is itself a function
of the parameters being estimated.

For simulation-based inference a summary statistic need not estimate anything;
it must be the *same deterministic function* of the data on both sides.

## Group A -- remove per-locus denominator decrements

Emitted statistics only. Each site currently subtracts the locus from that
statistic's denominator; after the edit the locus stays in the denominator and
contributes 0 (the accumulation is already skipped, so no `+= 0` is needed).

| Line | Case | Statistic | Trigger | Action |
|---|---|---|---|---|
| 1790 | 3 | `tajd` per-pop | `else if(pritruesegsites[ipop] == 0)` | delete decrement |
| 1815 | 3 | `tajd` global | `else if(truesegsites == 0)` | delete decrement |
| 1836 | 4 | `ZnS` per-pop | `else if(pritruesegsites[ipop] <= 1)` | delete decrement |
| 1847 | 4 | `ZnS` global | `else` on `truesegsites > 1 && ZnS(...)` | delete decrement |
| 1885 | 5 | `Fst` global | `else` on Fst calculation | delete decrement |
| 1926 | 6 | `shared` global | `else` | delete decrement |
| 1969 | 7 | `private` global | `else` | delete decrement |
| 2012 | 8 | `fixed` global | `else` | delete decrement |
| 2065 | 9 | pairwise `Fst` global | `else` | delete decrement |

Where the decrement is the sole statement of an `else`/`else if`, remove the
whole clause. Where it is inside a braced block, remove just the statement.

## Group B -- per-pair decrements: split the condition, do NOT delete

Sites 1921, 1964, 2007 (shared / private / fixed per pair) and 2054, 2058
(pairwise Fst) use a **mixed** trigger:

```c
if( (pars.cp.config)[npopi] >= MINSEQ && (pars.cp.config)[npopj] >= MINSEQ &&
    !isnan(shared[npopi][npopj]) ){
      ...accumulate...
} else
      sstats_popdenpair_weights[ii][npopi][npopj] -= sstats_weights[ii];
```

- `config[...] >= MINSEQ` is **population-level and constant across loci** --
  a population genuinely too small to support the statistic. This decrement
  must be KEPT.
- `!isnan(...)` is **per-locus** -- the pair carries no informative data at
  this locus. This must NOT decrement; it should contribute 0, matching the
  2026-05-26 `compute_sumstats.c` fix ("loci with piT_ij=0 contribute 0 to
  match the zero-init behavior of other empty-locus stats").

Required restructuring at each of the five sites:

```c
if( (pars.cp.config)[npopi] >= MINSEQ && (pars.cp.config)[npopj] >= MINSEQ ){
    if( !isnan(shared[npopi][npopj]) ){
        ...accumulate...
    }
    /* NaN at this locus: contributes 0, denominator unchanged */
} else {
    sstats_popdenpair_weights[ii][npopi][npopj] -= sstats_weights[ii];
}
```

## Group C -- `nhap` on monomorphic loci

Independent of the denominators. `dvstat` returns `nhap = 1, Hd = 0` for a
monomorphic locus. `ms.c` accepts it and adds 1; `compute_sumstats.c` skips
the locus and adds 0. This is why observed `s_mean_nhap_1` = 0.898 -- below 1,
which is impossible for a per-locus haplotype count -- while `Hd` matches
exactly (0 either way, from the same `dvstat` call).

Guard both accumulation sites in case 11 with `segsites > 0`:

- ~line 2185, per-pop: `dvstat(list, ..., samples_b[ipop], samples_e[ipop], &pridvk[ipop], &pridvh[ipop])`
- ~line 2203, global: `if( dvstat(list, pars.cp.nsam, segsites, 0, pars.cp.nsam, &dvk, &dvh))`

Do **not** decrement the denominator here -- monomorphic loci must remain in
it and contribute 0. `Hd` is unaffected (it contributed 0 before and after).

## Group D -- variance definition (optional; ~0.01% at n = 10,000)

```c
double varstat( double sumsq, double ave, double den){
  return (sumsq/den - ave*ave)/(den - 1)* den;   /* Bessel-corrected */
}
```

`compute_sumstats.c` uses the uncorrected population variance
`sum2/nl - mean*mean`. Note `skewstat`/`kurtstat` in `ms.c` already use the
uncorrected `m2`, so `ms.c` is internally inconsistent too.

Recommended: drop the correction in `varstat` so all three paths agree exactly.
Impact at n = 10,000 loci is a factor 1.0001 -- immaterial numerically, but it
removes the last known divergence.

## Sites to LEAVE UNCHANGED

| Line | Reason |
|---|---|
| 1723 | per-pop `pi`, trigger `config[ipop] >= MINSEQ` -- population-level, constant |
| 2198 | per-pop `nhap`/`Hd`, trigger `config[ipop] >= MINSEQ` -- population-level, constant |
| 2132, 2175 | Fay-Wu H -- computed but never emitted |
| 2232, 2255 | Thomson estimator/variance -- computed but never emitted |

Emitted statistics confirmed from a live RAD reftable header: `segs`, `pi`,
`thetaW`, `tajd`, `ZnS`, `nhap`, `Hd`, `Fst`, plus `shared_i_j`,
`private_i_j`, `fixed_i_j`, `Fst_i_j` and per-population variants.

## Non-issues (checked, no action)

- **Locus length weighting.** `sstats_weights[ii] = (double)flength[i]` at line
  1338 is gated by `if(weightmode == 1)`, commented "this is not used now";
  the default at line 271 is `1.`. Weights are uniformly 1, so the weighted
  mean reduces to the plain mean. (Would become a real divergence if
  `weightmode` were ever enabled with variable-length loci, since `var_*`
  accumulates *unweighted* squares while the denominator would be weighted.)
- **skew / kurtosis formulas.** `skewstat` and `kurtstat` are algebraically
  identical to `compute_sumstats.c`, including the `- 3.0` excess-kurtosis
  convention harmonised on 2026-03-27.
- **Per-locus kernels.** Shared via `data_sumstat.c` in all three paths.

## Consequences

- **Invalidates:** the four RAD reference tables (PonAbe + OoA, informative +
  diffuse). Regeneration ~10 min each; retraining ~1 h each.
- **Does not touch:** any observed file (RAD observed already uses the target
  convention), any WGS reference table, or the scrm path.
- **Performance:** none. Group A/B remove arithmetic; Group C adds one integer
  comparison per locus; Group D removes two operations.

## Verification plan

1. Rebuild, simulate one small RAD dataset through `sim.sumstats()`.
2. Write the same simulated data out and recompute with `observed.sumstats()`.
3. The two vectors must agree to floating-point tolerance on all
   ~104 statistics. Today they diverge on `ZnS` (5-17x), `nhap` (1.3-1.5x),
   `Fst` (1.8x) and `tajd`, and agree on `segs`/`pi`/`thetaW`/`Hd`.

This round-trip is the regression test that should have existed and did not.
Worth adding to the package tests regardless of the outcome here.
