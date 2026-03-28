#!/usr/bin/env python3
"""
Generate RAD pseudo-observed VCFs for PipeMaster example data.

Three models:
  1. Vaquita2Epoch_1R22 (1 pop, 15 diploid)
  2. PonAbe TwoSpecies_2L11 (2 pops, 10+10 diploid)
  3. OutOfAfrica_3G09 (3 pops, 20+20+20 diploid)

Each: 10,000 loci x 100bp, single replicate, VCF + pop_assign.txt output.

Requires: stdpopsim, msprime, numpy
Usage:
    python3 inst/extdata/generate_rad_vcfs.py
"""

import os
import sys
import numpy as np
import stdpopsim
import msprime

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)))
N_LOCI = 10000
LOCUS_LEN = 100
SEED = 42


def write_vcf(vcf_path, pop_assign_path, species_id, model_id,
              pop_samples, pop_names=None):
    """Simulate N_LOCI independent loci and write a multi-locus VCF.

    pop_samples: dict mapping population name -> n_diploid
    pop_names: optional dict mapping population name -> pop number (1-based)
    """
    species = stdpopsim.get_species(species_id)
    model = species.get_demographic_model(model_id)
    engine = stdpopsim.get_engine("msprime")
    mu = model.mutation_rate

    # Ordered population list
    ordered_pops = list(pop_samples.keys())
    if pop_names is None:
        pop_names = {p: i + 1 for i, p in enumerate(ordered_pops)}

    # Build sample name list: ind001..indN per pop, ordered by pop
    sample_names = []
    pop_assign = []
    idx = 1
    for pop in ordered_pops:
        n_dip = pop_samples[pop]
        for _ in range(n_dip):
            name = f"ind{idx:03d}"
            sample_names.append(name)
            pop_assign.append((name, pop_names[pop]))
            idx += 1

    n_diploid_total = sum(pop_samples.values())
    n_haploid = 2 * n_diploid_total

    rng = np.random.RandomState(SEED)
    seeds = rng.randint(1, 2**31, size=N_LOCI)

    total_snps = 0
    print(f"  Model: {model_id}")
    print(f"  mu = {mu:.2e}")
    print(f"  {N_LOCI} loci x {LOCUS_LEN} bp")
    print(f"  Samples: {pop_samples}")

    with open(vcf_path, "w") as vcf:
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write(f"##source=stdpopsim_{model_id}_msprime\n")
        vcf.write(f"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        for i in range(N_LOCI):
            vcf.write(f"##contig=<ID=locus{i+1},length={LOCUS_LEN}>\n")
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                  + "\t".join(sample_names) + "\n")

        for i in range(N_LOCI):
            if (i + 1) % 2000 == 0:
                print(f"    {i+1}/{N_LOCI} loci...", flush=True)

            contig = species.get_contig(length=LOCUS_LEN, mutation_rate=mu)
            locus_seed = int(seeds[i])
            ts = engine.simulate(model, contig, pop_samples, seed=locus_seed)

            # Get samples ordered by population
            all_samp = []
            for pop in ordered_pops:
                pop_id = None
                for p in ts.populations():
                    if p.metadata.get("name", "") == pop:
                        pop_id = p.id
                        break
                if pop_id is None:
                    pop_id = ordered_pops.index(pop)
                all_samp.extend(ts.samples(population=pop_id))

            chrom = f"locus{i+1}"
            for var in ts.variants(samples=all_samp):
                if len(var.alleles) != 2:
                    continue
                pos = int(var.site.position) + 1
                gts = []
                for d in range(n_diploid_total):
                    a1 = var.genotypes[2 * d]
                    a2 = var.genotypes[2 * d + 1]
                    gts.append(f"{a1}|{a2}")
                vcf.write(f"{chrom}\t{pos}\t.\tA\tT\t.\tPASS\t.\tGT\t"
                          + "\t".join(gts) + "\n")
                total_snps += 1

    # Pop assign with header
    with open(pop_assign_path, "w") as f:
        f.write("sample\tpop\n")
        for name, pop_num in pop_assign:
            f.write(f"{name}\t{pop_num}\n")

    size_kb = os.path.getsize(vcf_path) / 1024
    print(f"  {total_snps} SNPs, {size_kb:.0f} KB")
    print(f"  VCF: {vcf_path}")
    print(f"  Pop assign: {pop_assign_path}")


if __name__ == "__main__":
    print("=" * 50)
    print("  Generating RAD pseudo-observed VCFs")
    print("=" * 50)

    # --- 1. Vaquita ---
    print("\n=== Vaquita2Epoch_1R22 ===")
    write_vcf(
        os.path.join(OUTDIR, "Vaquita2Epoch_RAD.vcf"),
        os.path.join(OUTDIR, "Vaquita2Epoch_pop_assign.txt"),
        "PhoSin", "Vaquita2Epoch_1R22",
        {"Vaquita": 15}
    )

    # --- 2. PonAbe ---
    print("\n=== PonAbe TwoSpecies_2L11 ===")
    write_vcf(
        os.path.join(OUTDIR, "PonAbe_TwoSpecies_RAD.vcf"),
        os.path.join(OUTDIR, "PonAbe_TwoSpecies_pop_assign.txt"),
        "PonAbe", "TwoSpecies_2L11",
        {"Bornean": 10, "Sumatran": 10}
    )

    # --- 3. OoA ---
    print("\n=== OutOfAfrica_3G09 ===")
    write_vcf(
        os.path.join(OUTDIR, "OutOfAfrica_3G09_RAD.vcf"),
        os.path.join(OUTDIR, "OutOfAfrica_3G09_pop_assign.txt"),
        "HomSap", "OutOfAfrica_3G09",
        {"YRI": 20, "CEU": 20, "CHB": 20}
    )

    print("\n" + "=" * 50)
    print("  Done!")
    print("=" * 50)
