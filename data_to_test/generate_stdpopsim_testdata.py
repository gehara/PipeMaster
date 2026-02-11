#!/usr/bin/env python3
"""
Generate test data for PipeMaster using stdpopsim.

Produces 4 datasets under known demographic models:
1. Africa_1T12 (1 pop)       -> 1D SFS
2. OutOfAfrica_2T12 (2 pop)  -> 2D joint SFS
3. OutOfAfrica_3G09 (2 pop subset) -> summary statistics
4. Vaquita2Epoch_1R22 (1 pop bottleneck) -> summary statistics

Output format matches PipeMaster conventions:
  - Tab-separated files
  - SFS: sfs_0, sfs_1, ... (1-pop) or sfs_i_j (2-pop)
  - Summary stats: s_average_*, s_variance_* (msABC naming)
"""

import stdpopsim
import numpy as np
import os
import sys

SEED = 42
OUTDIR = os.path.dirname(os.path.abspath(__file__))


# ============================================================================
# Helper functions
# ============================================================================

def compute_1d_sfs(ts, pop_id=None):
    """Compute 1D SFS from a tree sequence for a given population."""
    if pop_id is not None:
        sample_ids = ts.samples(population=pop_id)
    else:
        sample_ids = ts.samples()
    n = len(sample_ids)
    sfs = np.zeros(n + 1, dtype=int)
    for var in ts.variants(samples=sample_ids):
        count = int(np.sum(var.genotypes))
        sfs[count] += 1
    return sfs


def compute_joint_sfs_2d(ts, pop1_id, pop2_id):
    """Compute 2D joint SFS from a tree sequence."""
    samples_p1 = ts.samples(population=pop1_id)
    samples_p2 = ts.samples(population=pop2_id)
    n1, n2 = len(samples_p1), len(samples_p2)
    joint_sfs = np.zeros((n1 + 1, n2 + 1), dtype=int)

    all_samples = np.concatenate([samples_p1, samples_p2])
    for var in ts.variants(samples=all_samples):
        g = var.genotypes
        c1 = int(np.sum(g[:n1]))
        c2 = int(np.sum(g[n1:]))
        joint_sfs[c1, c2] += 1
    return joint_sfs


def _pop_basic_stats(ts, samp, n, pop_id):
    """Compute basic per-locus stats for one population sample."""
    from collections import Counter
    seg_sites = 0
    for var in ts.variants(samples=samp):
        freq = int(np.sum(var.genotypes))
        if 0 < freq < n:
            seg_sites += 1

    pi = float(np.squeeze(ts.diversity(sample_sets=[samp])))

    if seg_sites > 0:
        a1 = sum(1.0 / k for k in range(1, n))
        theta_w = seg_sites / a1
    else:
        theta_w = 0.0

    tajd = float(np.squeeze(ts.Tajimas_D(sample_sets=[samp])))
    if not np.isfinite(tajd):
        tajd = 0.0

    # Fay & Wu's H (pi - theta_H)
    sfs_pop = compute_1d_sfs(ts, pop_id=pop_id)
    theta_H = 0.0
    if n > 1:
        for freq_bin in range(1, n):
            theta_H += (freq_bin ** 2) * sfs_pop[freq_bin]
        theta_H = (2.0 / (n * (n - 1))) * theta_H
    fwh = pi - theta_H

    # DVK / DVH
    geno_matrix = ts.genotype_matrix(samples=samp)
    hap_list = []
    if geno_matrix.shape[0] > 0:
        for col_idx in range(n):
            hap_list.append(tuple(geno_matrix[:, col_idx]))
    dvk = len(set(hap_list))
    if dvk > 1 and n > 1:
        hc = Counter(hap_list)
        dvh = 1.0 - sum((c / n) ** 2 for c in hc.values())
        dvh = dvh * n / (n - 1)
    else:
        dvh = 0.0

    return seg_sites, pi, theta_w, tajd, fwh, dvk, dvh


def compute_perlocus_sumstats(ts, pop_ids):
    """
    Compute per-locus summary statistics matching PipeMaster sim.sumstat column naming.

    Single-pop: s_segs, s_pi, s_w, s_tajd, s_FayWuH, s_dvk, s_dvh
    Multi-pop:  per-pop (_1, _2), overall, pairwise (Fst, shared, private, fixed_dif, pairwise_fst)

    Column order matches sim.sumstat output:
      segs(per-pop, overall), pi(...), w(...), tajd(...), Fst, shared, private, fixed_dif,
      pairwise_fst, fwh(per-pop, overall), dvk(per-pop, overall), dvh(per-pop, overall)
    """
    from collections import Counter
    npop = len(pop_ids)
    stats = {}

    all_pop_samples = []
    all_pop_n = []
    pop_results = []

    for i, pid in enumerate(pop_ids):
        samp = ts.samples(population=pid)
        n = len(samp)
        all_pop_samples.append(samp)
        all_pop_n.append(n)
        pop_results.append(_pop_basic_stats(ts, samp, n, pid))

    if npop == 1:
        # Single-pop: no suffix, use FayWuH (not fwh)
        segs, pi, w, tajd, fwh, dvk, dvh = pop_results[0]
        stats["s_segs"] = segs
        stats["s_pi"] = pi
        stats["s_w"] = w
        stats["s_tajd"] = tajd
        stats["s_FayWuH"] = fwh
        stats["s_dvk"] = dvk
        stats["s_dvh"] = dvh
    else:
        # Multi-pop: PipeMaster groups by stat type: stat_1, stat_2, stat_overall
        # Compute overall stats
        all_samp = np.concatenate(all_pop_samples)
        n_total = len(all_samp)

        seg_total = 0
        for var in ts.variants(samples=all_samp):
            freq = int(np.sum(var.genotypes))
            if 0 < freq < n_total:
                seg_total += 1

        pi_total = float(np.squeeze(ts.diversity(sample_sets=[all_samp])))

        if seg_total > 0:
            a1_total = sum(1.0 / k for k in range(1, n_total))
            w_total = seg_total / a1_total
        else:
            w_total = 0.0

        tajd_total = float(np.squeeze(ts.Tajimas_D(sample_sets=[all_samp])))
        if not np.isfinite(tajd_total):
            tajd_total = 0.0

        # Group 1: segs per-pop + overall, pi per-pop + overall, w, tajd
        overall_basic = [seg_total, pi_total, w_total, tajd_total]
        basic_names = ["segs", "pi", "w", "tajd"]
        for si, sname in enumerate(basic_names):
            for i in range(npop):
                stats[f"s_{sname}_{i+1}"] = pop_results[i][si]
            stats[f"s_{sname}"] = overall_basic[si]

        # Group 2: Fst, shared, private, fixed_dif, pairwise_fst
        for i in range(npop):
            for j in range(i + 1, npop):
                samp_i = all_pop_samples[i]
                samp_j = all_pop_samples[j]
                ni, nj = all_pop_n[i], all_pop_n[j]
                pi_label = f"{i+1}_{j+1}"

                fst = float(np.squeeze(ts.Fst(sample_sets=[samp_i, samp_j])))
                stats["s_Fst"] = fst if np.isfinite(fst) else 0.0

                shared = 0
                private = 0
                fixed_dif = 0
                all_pair = np.concatenate([samp_i, samp_j])
                for var in ts.variants(samples=all_pair):
                    g = var.genotypes
                    c_i = int(np.sum(g[:ni]))
                    c_j = int(np.sum(g[ni:]))
                    poly_i = (0 < c_i < ni)
                    poly_j = (0 < c_j < nj)
                    if poly_i and poly_j:
                        shared += 1
                    elif poly_i and not poly_j and c_j == 0:
                        private += 1
                    elif poly_j and not poly_i and c_i == 0:
                        private += 1
                    if (c_i == ni and c_j == 0) or (c_j == nj and c_i == 0):
                        fixed_dif += 1

                stats[f"s_shared_{pi_label}"] = shared
                stats[f"s_private_{pi_label}"] = private
                stats[f"s_fixed_dif_{pi_label}"] = fixed_dif
                stats[f"s_pairwise_fst_{pi_label}"] = stats["s_Fst"]

        # Group 3: fwh per-pop + overall FayWuH
        for i in range(npop):
            stats[f"s_fwh_{i+1}"] = pop_results[i][4]

        sfs_total = compute_1d_sfs(ts)
        theta_H_total = 0.0
        if n_total > 1:
            for freq_bin in range(1, n_total):
                theta_H_total += (freq_bin ** 2) * sfs_total[freq_bin]
            theta_H_total = (2.0 / (n_total * (n_total - 1))) * theta_H_total
        stats["s_FayWuH"] = pi_total - theta_H_total

        # Group 4: dvk per-pop + overall, dvh per-pop + overall
        for i in range(npop):
            stats[f"s_dvk_{i+1}"] = pop_results[i][5]
            stats[f"s_dvh_{i+1}"] = pop_results[i][6]

        geno_all = ts.genotype_matrix(samples=all_samp)
        hap_list_all = []
        if geno_all.shape[0] > 0:
            for col_idx in range(n_total):
                hap_list_all.append(tuple(geno_all[:, col_idx]))
        stats["s_dvk"] = len(set(hap_list_all))
        if stats["s_dvk"] > 1 and n_total > 1:
            hc = Counter(hap_list_all)
            dvh_all = 1.0 - sum((c / n_total) ** 2 for c in hc.values())
            stats["s_dvh"] = dvh_all * n_total / (n_total - 1)
        else:
            stats["s_dvh"] = 0.0

    return stats


def simulate_loci(species, model, contig_length, n_loci, samples_dict, seed=42):
    """Simulate n_loci independent contigs and return list of tree sequences."""
    engine = stdpopsim.get_engine("msprime")
    rng = np.random.RandomState(seed)
    tree_seqs = []

    mu = model.mutation_rate if model.mutation_rate is not None else None
    for i in range(n_loci):
        contig = species.get_contig(length=contig_length, mutation_rate=mu)
        locus_seed = int(rng.randint(1, 2**31))
        ts = engine.simulate(model, contig, samples_dict, seed=locus_seed)
        tree_seqs.append(ts)
        if (i + 1) % 50 == 0:
            print(f"  Simulated {i+1}/{n_loci} loci", file=sys.stderr)

    return tree_seqs


def write_tsv(data_rows, fieldnames, filepath):
    """Write list of dicts to tab-separated file (PipeMaster convention)."""
    with open(filepath, "w") as f:
        f.write("\t".join(fieldnames) + "\n")
        for row in data_rows:
            vals = [str(row.get(k, "NA")) for k in fieldnames]
            f.write("\t".join(vals) + "\n")
    print(f"  Wrote {filepath} ({len(data_rows)} rows)")


def get_pop_id(ts, name_hint, fallback_id):
    """Find population ID by name hint, falling back to numeric ID."""
    for p in ts.populations():
        pname = p.metadata.get("name", "")
        if name_hint in pname:
            return p.id
    return fallback_id


def compute_moments(all_locus_stats, stat_keys):
    """
    Compute 4 moments (mean, var, skew, kurt) across loci for each statistic,
    matching PipeMaster sim.sumstat naming: s_mean_X, s_var_X, s_skew_X, s_kurt_X.
    """
    from scipy import stats as sp_stats
    observed = {}
    for key in stat_keys:
        vals = np.array([s[key] for s in all_locus_stats], dtype=float)
        clean = vals[np.isfinite(vals)]
        suffix = key.replace('s_', '', 1)
        if len(clean) > 0:
            observed[f"s_mean_{suffix}"] = float(np.mean(clean))
            observed[f"s_var_{suffix}"] = float(np.var(clean))
            observed[f"s_skew_{suffix}"] = float(sp_stats.skew(clean))
            observed[f"s_kurt_{suffix}"] = float(sp_stats.kurtosis(clean))
        else:
            observed[f"s_mean_{suffix}"] = 0.0
            observed[f"s_var_{suffix}"] = 0.0
            observed[f"s_skew_{suffix}"] = 0.0
            observed[f"s_kurt_{suffix}"] = 0.0
    return observed


def export_fasta_files(tree_seqs, outdir, pop_sample_map, contig_length):
    """
    Export FASTA files for each locus tree sequence.

    Args:
        tree_seqs: list of tree sequences (one per locus)
        outdir: output directory for FASTA files
        pop_sample_map: list of (pop_name_hint, pop_fallback_id, pop_number, n_diploid)
                        tuples defining pop ordering. pop_number is 1-based for pop.assign.
        contig_length: length of each contig (bp)

    Writes:
        - One .fas file per locus in outdir (locus_000.fas, locus_001.fas, ...)
        - A pop_assign_<basename>.txt file alongside outdir
    """
    os.makedirs(outdir, exist_ok=True)

    # Build ordered sample indices and names from first tree sequence
    ts0 = tree_seqs[0]
    ordered_indices = []
    sample_names = []
    pop_numbers = []
    counter = 1
    for name_hint, fallback_id, pop_num, n_dip in pop_sample_map:
        pid = get_pop_id(ts0, name_hint, fallback_id)
        samp = ts0.samples(population=pid)
        for s in samp:
            ordered_indices.append(s)
            sample_names.append(f"sample_{counter}")
            pop_numbers.append(pop_num)
            counter += 1

    n_samples = len(ordered_indices)

    for locus_idx, ts in enumerate(tree_seqs):
        filepath = os.path.join(outdir, f"locus_{locus_idx:03d}.fas")

        # Initialize all-A sequences
        seqs = [['A'] * contig_length for _ in range(n_samples)]

        # Rebuild sample index mapping for this locus tree sequence
        locus_indices = []
        for name_hint, fallback_id, pop_num, n_dip in pop_sample_map:
            pid = get_pop_id(ts, name_hint, fallback_id)
            samp = ts.samples(population=pid)
            locus_indices.extend(samp)

        # Build index lookup: global position -> position in locus_indices
        idx_to_pos = {int(s): i for i, s in enumerate(locus_indices)}

        # Apply variants
        all_samp = np.array(locus_indices)
        for var in ts.variants(samples=all_samp):
            pos = int(var.site.position)
            if pos >= contig_length:
                continue
            for i in range(n_samples):
                if var.genotypes[i] != 0:
                    seqs[i][pos] = 'T'

        # Write FASTA
        with open(filepath, 'w') as f:
            for i in range(n_samples):
                f.write(f'>{sample_names[i]}\n')
                f.write(''.join(seqs[i]) + '\n')

    print(f"  Wrote {len(tree_seqs)} FASTA files to {outdir}")

    # Write pop_assign file
    basename = os.path.basename(outdir).replace('fasta_', '')
    pop_assign_path = os.path.join(os.path.dirname(outdir), f"pop_assign_{basename}.txt")
    with open(pop_assign_path, 'w') as f:
        for i in range(n_samples):
            f.write(f'{sample_names[i]}\t{pop_numbers[i]}\n')
    print(f"  Wrote {pop_assign_path} ({n_samples} samples)")


# ============================================================================
# 1. Africa_1T12 - Single population 1D SFS
# ============================================================================
print("=== 1. Africa_1T12: 1D SFS (single population) ===")

species_hs = stdpopsim.get_species("HomSap")
model_1 = species_hs.get_demographic_model("Africa_1T12")
n_samples_1 = 20  # diploid individuals -> 40 haplotypes
n_loci_1 = 200
contig_len_1 = 1000  # 1kb per locus

tree_seqs_1 = simulate_loci(
    species_hs, model_1, contig_len_1, n_loci_1,
    {"AFR": n_samples_1}, seed=SEED
)

# Sum SFS across all loci
n_hap_1 = n_samples_1 * 2  # 40 haplotypes
total_sfs_1 = np.zeros(n_hap_1 + 1, dtype=int)
for ts in tree_seqs_1:
    total_sfs_1 += compute_1d_sfs(ts)

# Write observed SFS (tab-separated, matching PipeMaster sfs_N naming)
sfs_names = [f"sfs_{i}" for i in range(len(total_sfs_1))]
row = {name: int(val) for name, val in zip(sfs_names, total_sfs_1)}
write_tsv([row], sfs_names, os.path.join(OUTDIR, "observed_sfs_Africa_1T12.txt"))

meta = {
    "model": "Africa_1T12",
    "species": "HomSap",
    "n_diploid": n_samples_1,
    "n_haploid": n_hap_1,
    "n_loci": n_loci_1,
    "contig_length": contig_len_1,
    "mutation_rate": model_1.mutation_rate,
    "description": "Single African population, 3-epoch (Tennessen et al. 2012)"
}
write_tsv([meta], list(meta.keys()), os.path.join(OUTDIR, "meta_Africa_1T12.txt"))

# Export FASTA files for Africa_1T12
print("  Exporting FASTA files for Africa_1T12...")
export_fasta_files(
    tree_seqs_1,
    os.path.join(OUTDIR, "fasta_Africa_1T12"),
    pop_sample_map=[
        ("AFR", 0, 1, n_samples_1),  # pop 1: AFR
    ],
    contig_length=contig_len_1,
)


# ============================================================================
# 2. OutOfAfrica_2T12 - Two population joint 2D SFS
# ============================================================================
print("\n=== 2. OutOfAfrica_2T12: 2D joint SFS (AFR + EUR) ===")

model_2 = species_hs.get_demographic_model("OutOfAfrica_2T12")
n_afr_2, n_eur_2 = 20, 20  # diploid
n_loci_2 = 200
contig_len_2 = 1000

tree_seqs_2 = simulate_loci(
    species_hs, model_2, contig_len_2, n_loci_2,
    {"AFR": n_afr_2, "EUR": n_eur_2}, seed=SEED + 1
)

# Sum joint SFS across loci
nh_afr, nh_eur = n_afr_2 * 2, n_eur_2 * 2  # 40, 40
total_joint_sfs = np.zeros((nh_afr + 1, nh_eur + 1), dtype=int)
for ts in tree_seqs_2:
    afr_id = get_pop_id(ts, "AFR", 0)
    eur_id = get_pop_id(ts, "EUR", 1)
    total_joint_sfs += compute_joint_sfs_2d(ts, afr_id, eur_id)

# Flatten to vector using expand.grid order (PipeMaster convention: pop1 varies fastest)
# PipeMaster uses: expand.grid(0:n1, 0:n2) -> pop1 cycles first
# This means row-major with pop1 as inner loop
idx_grid = [(i, j) for j in range(nh_eur + 1) for i in range(nh_afr + 1)]
sfs_names_2d = [f"sfs_{i}_{j}" for i, j in idx_grid]
flat_vals = [int(total_joint_sfs[i, j]) for i, j in idx_grid]
row_2d = dict(zip(sfs_names_2d, flat_vals))
write_tsv([row_2d], sfs_names_2d, os.path.join(OUTDIR, "observed_sfs_OutOfAfrica_2T12.txt"))

# Also save as matrix for plotting with plot.2D.sfs
np.savetxt(os.path.join(OUTDIR, "observed_joint_sfs_matrix_OutOfAfrica_2T12.txt"),
           total_joint_sfs, delimiter="\t", fmt="%d")
print(f"  Wrote joint SFS matrix ({total_joint_sfs.shape})")

meta_2 = {
    "model": "OutOfAfrica_2T12",
    "species": "HomSap",
    "n_diploid_AFR": n_afr_2,
    "n_diploid_EUR": n_eur_2,
    "n_haploid_AFR": nh_afr,
    "n_haploid_EUR": nh_eur,
    "n_loci": n_loci_2,
    "contig_length": contig_len_2,
    "mutation_rate": model_2.mutation_rate,
    "description": "African + European Out-of-Africa (Tennessen et al. 2012)"
}
write_tsv([meta_2], list(meta_2.keys()), os.path.join(OUTDIR, "meta_OutOfAfrica_2T12.txt"))

# Export FASTA files for OutOfAfrica_2T12
print("  Exporting FASTA files for OutOfAfrica_2T12...")
export_fasta_files(
    tree_seqs_2,
    os.path.join(OUTDIR, "fasta_OutOfAfrica_2T12"),
    pop_sample_map=[
        ("AFR", 0, 1, n_afr_2),   # pop 1: AFR
        ("EUR", 1, 2, n_eur_2),   # pop 2: EUR
    ],
    contig_length=contig_len_2,
)


# ============================================================================
# 3. OutOfAfrica_3G09 - Two-pop summary stats (YRI + CEU)
# ============================================================================
print("\n=== 3. OutOfAfrica_3G09: Summary statistics (YRI + CEU) ===")

model_3 = species_hs.get_demographic_model("OutOfAfrica_3G09")
n_yri, n_ceu = 20, 20
n_loci_3 = 200
contig_len_3 = 1000

tree_seqs_3 = simulate_loci(
    species_hs, model_3, contig_len_3, n_loci_3,
    {"YRI": n_yri, "CEU": n_ceu, "CHB": 0}, seed=SEED + 2
)

# Compute per-locus summary stats
all_stats_3 = []
for k, ts in enumerate(tree_seqs_3):
    yri_id = get_pop_id(ts, "YRI", 0)
    ceu_id = get_pop_id(ts, "CEU", 1)
    stats = compute_perlocus_sumstats(ts, pop_ids=[yri_id, ceu_id])
    stats["locus"] = k + 1
    all_stats_3.append(stats)

# Per-locus file
stat_keys = [k for k in all_stats_3[0].keys() if k != "locus"]
fields_locus = ["locus"] + stat_keys
write_tsv(all_stats_3, fields_locus, os.path.join(OUTDIR, "perlocus_sumstats_OutOfAfrica_3G09.txt"))

# Observed: mean and variance across loci (msABC naming)
observed_3 = compute_moments(all_stats_3, stat_keys)
write_tsv([observed_3], list(observed_3.keys()), os.path.join(OUTDIR, "observed_sumstats_OutOfAfrica_3G09.txt"))

meta_3 = {
    "model": "OutOfAfrica_3G09",
    "species": "HomSap",
    "pops_sampled": "YRI,CEU",
    "n_diploid_YRI": n_yri,
    "n_diploid_CEU": n_ceu,
    "n_loci": n_loci_3,
    "contig_length": contig_len_3,
    "mutation_rate": model_3.mutation_rate,
    "description": "Out-of-Africa YRI+CEU (Gutenkunst et al. 2009)"
}
write_tsv([meta_3], list(meta_3.keys()), os.path.join(OUTDIR, "meta_OutOfAfrica_3G09.txt"))

# Export FASTA files for OutOfAfrica_3G09
print("  Exporting FASTA files for OutOfAfrica_3G09...")
export_fasta_files(
    tree_seqs_3,
    os.path.join(OUTDIR, "fasta_OutOfAfrica_3G09"),
    pop_sample_map=[
        ("YRI", 0, 1, n_yri),   # pop 1: YRI
        ("CEU", 1, 2, n_ceu),   # pop 2: CEU
    ],
    contig_length=contig_len_3,
)


# ============================================================================
# 4. Vaquita2Epoch_1R22 - Single pop bottleneck, summary stats
# ============================================================================
print("\n=== 4. Vaquita2Epoch_1R22: Summary statistics (bottleneck) ===")

species_ps = stdpopsim.get_species("PhoSin")
model_4 = species_ps.get_demographic_model("Vaquita2Epoch_1R22")
n_vaq = 20
n_loci_4 = 200
contig_len_4 = 1000

tree_seqs_4 = simulate_loci(
    species_ps, model_4, contig_len_4, n_loci_4,
    {"Vaquita": n_vaq}, seed=SEED + 3
)

# Per-locus summary stats
all_stats_4 = []
for k, ts in enumerate(tree_seqs_4):
    stats = compute_perlocus_sumstats(ts, pop_ids=[0])
    stats["locus"] = k + 1
    all_stats_4.append(stats)

stat_keys_4 = [k for k in all_stats_4[0].keys() if k != "locus"]
fields_locus_4 = ["locus"] + stat_keys_4
write_tsv(all_stats_4, fields_locus_4, os.path.join(OUTDIR, "perlocus_sumstats_Vaquita2Epoch.txt"))

# Observed moments (msABC naming)
observed_4 = compute_moments(all_stats_4, stat_keys_4)
write_tsv([observed_4], list(observed_4.keys()), os.path.join(OUTDIR, "observed_sumstats_Vaquita2Epoch.txt"))

# Also compute 1D SFS
nh_vaq = n_vaq * 2
total_sfs_4 = np.zeros(nh_vaq + 1, dtype=int)
for ts in tree_seqs_4:
    total_sfs_4 += compute_1d_sfs(ts)

sfs_names_4 = [f"sfs_{i}" for i in range(len(total_sfs_4))]
row_4 = {name: int(val) for name, val in zip(sfs_names_4, total_sfs_4)}
write_tsv([row_4], sfs_names_4, os.path.join(OUTDIR, "observed_sfs_Vaquita2Epoch.txt"))

meta_4 = {
    "model": "Vaquita2Epoch_1R22",
    "species": "PhoSin",
    "n_diploid": n_vaq,
    "n_haploid": nh_vaq,
    "n_loci": n_loci_4,
    "contig_length": contig_len_4,
    "mutation_rate": model_4.mutation_rate,
    "description": "Vaquita 2-epoch bottleneck (Robinson et al. 2022)"
}
write_tsv([meta_4], list(meta_4.keys()), os.path.join(OUTDIR, "meta_Vaquita2Epoch.txt"))

# Export FASTA files for Vaquita2Epoch
print("  Exporting FASTA files for Vaquita2Epoch...")
export_fasta_files(
    tree_seqs_4,
    os.path.join(OUTDIR, "fasta_Vaquita2Epoch"),
    pop_sample_map=[
        ("Vaquita", 0, 1, n_vaq),  # pop 1: Vaquita
    ],
    contig_length=contig_len_4,
)


# ============================================================================
# Summary
# ============================================================================
print("\n=== Done! Generated files: ===")
for f in sorted(os.listdir(OUTDIR)):
    if f.endswith(".txt") and ("observed" in f or "meta" in f or "perlocus" in f or "pop_assign" in f):
        fpath = os.path.join(OUTDIR, f)
        print(f"  {f}  ({os.path.getsize(fpath)} bytes)")
for d in sorted(os.listdir(OUTDIR)):
    dpath = os.path.join(OUTDIR, d)
    if os.path.isdir(dpath) and d.startswith("fasta_"):
        n_fas = len([x for x in os.listdir(dpath) if x.endswith(".fas")])
        print(f"  {d}/  ({n_fas} FASTA files)")
