#!/usr/bin/env Rscript
# Test script for sim.sfs() function
# Validates SFS simulation for 1-pop, 2-pop, folded, and custom mu.rates

suppressMessages(library(PipeMaster))

# Determine script directory robustly
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if(length(file_arg) > 0) {
  test_dir <- normalizePath(dirname(sub("^--file=", "", file_arg)))
} else {
  test_dir <- normalizePath("data_to_test")
}
setwd(test_dir)

cat("=== Loading test data ===\n")
load("test_models.RData")

pass <- 0
fail <- 0

check <- function(desc, condition) {
  if(isTRUE(condition)) {
    cat(paste0("  PASS: ", desc, "\n"))
    pass <<- pass + 1
  } else {
    cat(paste0("  FAIL: ", desc, "\n"))
    fail <<- fail + 1
  }
}

# Create temporary directory for sim output
tmpdir <- tempdir()

# ============================================================================
# Test 1: Vaquita2Epoch (1-pop bottleneck)
# ============================================================================
cat("\n=== Test 1: Vaquita2Epoch (1-pop, sim.sfs) ===\n")

sim.sfs(model = Vaquita2Epoch,
        nsim.blocks = 1,
        block.size = 20,
        path = tmpdir,
        use.alpha = FALSE,
        output.name = "test_vaq",
        ncores = 1)

outfile_vaq <- file.path(tmpdir, "SIM_SFS_test_vaq.txt")
check("Vaquita: output file exists", file.exists(outfile_vaq))

sim_vaq <- read.table(outfile_vaq, header = TRUE, sep = "\t")
check("Vaquita: 20 simulation rows", nrow(sim_vaq) == 20)

# Check parameter columns exist
check("Vaquita: has Ne0.pop1 column", "Ne0.pop1" %in% colnames(sim_vaq))
check("Vaquita: has Ne1.pop1 column", "Ne1.pop1" %in% colnames(sim_vaq))
check("Vaquita: has t.Ne1.pop1 column", "t.Ne1.pop1" %in% colnames(sim_vaq))
check("Vaquita: has mean.rate column", "mean.rate" %in% colnames(sim_vaq))
check("Vaquita: has sd.rate column", "sd.rate" %in% colnames(sim_vaq))

# Check SFS columns
# scrm -oSFS produces nsam-1 bins (freq classes 1..nsam-1), labeled sfs_0..sfs_{nsam-2}
sfs_cols_vaq <- grep("^sfs_", colnames(sim_vaq), value = TRUE)
check("Vaquita: 39 SFS columns (sfs_0..sfs_38)", length(sfs_cols_vaq) == 39)
check("Vaquita: first SFS col is sfs_0", sfs_cols_vaq[1] == "sfs_0")
check("Vaquita: last SFS col is sfs_38", sfs_cols_vaq[39] == "sfs_38")

# Check SFS values
sfs_mat_vaq <- as.matrix(sim_vaq[, sfs_cols_vaq])
check("Vaquita: all SFS values non-negative", all(sfs_mat_vaq >= 0))
row_totals_vaq <- rowSums(sfs_mat_vaq)
check("Vaquita: most rows have SNPs (>75%)", mean(row_totals_vaq > 0) > 0.75)

# Compare mean simulated SFS with observed reference
# sim.sfs produces 39 bins (freq 1..39); obs.sfs produces 41 bins (freq 0..40)
# Align by taking obs bins 2:40 (freq 1..39) to match sim bins
mean_sfs_vaq <- colMeans(sfs_mat_vaq)
ref_vaq <- as.numeric(observed_sfs_Vaquita2Epoch[1, 2:40])
nonzero <- (mean_sfs_vaq + ref_vaq) > 0
if(sum(nonzero) > 2) {
  r <- cor(mean_sfs_vaq[nonzero], ref_vaq[nonzero])
  cat(paste0("  Vaquita mean sim SFS vs observed correlation: r = ", round(r, 4), "\n"))
  check("Vaquita: SFS correlation > 0.5", r > 0.5)
} else {
  cat("  WARNING: too few non-zero bins for correlation\n")
}

# Cleanup
file.remove(outfile_vaq)


# ============================================================================
# Test 2: Africa_1T12 (1-pop, 3-epoch)
# ============================================================================
cat("\n=== Test 2: Africa_1T12 (1-pop, 3-epoch, sim.sfs) ===\n")

sim.sfs(model = Africa_1T12,
        nsim.blocks = 1,
        block.size = 20,
        path = tmpdir,
        use.alpha = FALSE,
        output.name = "test_afr",
        ncores = 1)

outfile_afr <- file.path(tmpdir, "SIM_SFS_test_afr.txt")
check("Africa: output file exists", file.exists(outfile_afr))

sim_afr <- read.table(outfile_afr, header = TRUE, sep = "\t")
check("Africa: 20 simulation rows", nrow(sim_afr) == 20)

# Check parameter columns (3-epoch: Ne0, Ne1, Ne2, t.Ne1, t.Ne2)
check("Africa: has Ne0.pop1 column", "Ne0.pop1" %in% colnames(sim_afr))
check("Africa: has Ne1.pop1 column", "Ne1.pop1" %in% colnames(sim_afr))
check("Africa: has Ne2.pop1 column", "Ne2.pop1" %in% colnames(sim_afr))
check("Africa: has t.Ne1.pop1 column", "t.Ne1.pop1" %in% colnames(sim_afr))
check("Africa: has t.Ne2.pop1 column", "t.Ne2.pop1" %in% colnames(sim_afr))

# Check time ordering constraint (t.Ne1 < t.Ne2 always)
check("Africa: t.Ne1.pop1 < t.Ne2.pop1 (constraint enforced)",
      all(sim_afr$t.Ne1.pop1 < sim_afr$t.Ne2.pop1))

# Check SFS columns (39 bins: scrm -oSFS for nsam=40)
sfs_cols_afr <- grep("^sfs_", colnames(sim_afr), value = TRUE)
check("Africa: 39 SFS columns", length(sfs_cols_afr) == 39)

sfs_mat_afr <- as.matrix(sim_afr[, sfs_cols_afr])
check("Africa: all SFS values non-negative", all(sfs_mat_afr >= 0))

# Compare with observed (align bins: obs cols 2:40 = freq 1..39)
mean_sfs_afr <- colMeans(sfs_mat_afr)
ref_afr <- as.numeric(observed_sfs_Africa_1T12[1, 2:40])
nonzero <- (mean_sfs_afr + ref_afr) > 0
if(sum(nonzero) > 2) {
  r <- cor(mean_sfs_afr[nonzero], ref_afr[nonzero])
  cat(paste0("  Africa mean sim SFS vs observed correlation: r = ", round(r, 4), "\n"))
  check("Africa: SFS correlation > 0.5", r > 0.5)
}

file.remove(outfile_afr)


# ============================================================================
# Test 3: OutOfAfrica_2T12 (2-pop joint SFS)
# ============================================================================
cat("\n=== Test 3: OutOfAfrica_2T12 (2-pop joint SFS, sim.sfs) ===\n")

sim.sfs(model = OutOfAfrica_2T12,
        nsim.blocks = 1,
        block.size = 10,
        path = tmpdir,
        use.alpha = FALSE,
        output.name = "test_ooa2",
        ncores = 1)

outfile_ooa2 <- file.path(tmpdir, "SIM_SFS_test_ooa2.txt")
check("OoA_2T12: output file exists", file.exists(outfile_ooa2))

sim_ooa2 <- read.table(outfile_ooa2, header = TRUE, sep = "\t")
check("OoA_2T12: 10 simulation rows", nrow(sim_ooa2) == 10)

# Check parameter columns
check("OoA_2T12: has Ne0.pop1 column", "Ne0.pop1" %in% colnames(sim_ooa2))
check("OoA_2T12: has Ne0.pop2 column", "Ne0.pop2" %in% colnames(sim_ooa2))
check("OoA_2T12: has join1 column", "join1" %in% colnames(sim_ooa2))
check("OoA_2T12: has mig0.1_2 column", "mig0.1_2" %in% colnames(sim_ooa2))

# Check SFS columns (41*41 = 1681)
sfs_cols_ooa2 <- grep("^sfs_", colnames(sim_ooa2), value = TRUE)
check("OoA_2T12: 1681 SFS columns (41x41)", length(sfs_cols_ooa2) == 1681)
check("OoA_2T12: first SFS col is sfs_0_0", sfs_cols_ooa2[1] == "sfs_0_0")

# Check expand.grid ordering: pop1 varies fastest
# Second column should be sfs_1_0 (pop1=1, pop2=0)
check("OoA_2T12: second SFS col is sfs_1_0 (pop1 varies fastest)",
      sfs_cols_ooa2[2] == "sfs_1_0")

sfs_mat_ooa2 <- as.matrix(sim_ooa2[, sfs_cols_ooa2])
check("OoA_2T12: all SFS values non-negative", all(sfs_mat_ooa2 >= 0))
check("OoA_2T12: all rows have SNPs", all(rowSums(sfs_mat_ooa2) > 0))

file.remove(outfile_ooa2)


# ============================================================================
# Test 4: Folded SFS (1-pop)
# ============================================================================
cat("\n=== Test 4: Folded SFS (sim.sfs) ===\n")

sim.sfs(model = Vaquita2Epoch,
        nsim.blocks = 1,
        block.size = 10,
        path = tmpdir,
        use.alpha = FALSE,
        output.name = "test_folded",
        folded = TRUE,
        ncores = 1)

outfile_fold <- file.path(tmpdir, "SIM_SFS_test_folded.txt")
check("Folded: output file exists", file.exists(outfile_fold))

sim_fold <- read.table(outfile_fold, header = TRUE, sep = "\t")
check("Folded: 10 simulation rows", nrow(sim_fold) == 10)

sfs_cols_fold <- grep("^sfs_fold_", colnames(sim_fold), value = TRUE)
# scrm produces 39 unfolded bins -> fold_sfs: floor(39/2)=19 + 1 middle (odd) = 20
check("Folded: 20 SFS columns", length(sfs_cols_fold) == 20)
check("Folded: first col is sfs_fold_0", sfs_cols_fold[1] == "sfs_fold_0")
check("Folded: last col is sfs_fold_19", sfs_cols_fold[20] == "sfs_fold_19")

sfs_mat_fold <- as.matrix(sim_fold[, sfs_cols_fold])
check("Folded: all values non-negative", all(sfs_mat_fold >= 0))

file.remove(outfile_fold)


# ============================================================================
# Test 5: Custom mutation rates (mu.rates)
# ============================================================================
cat("\n=== Test 5: Custom mu.rates (sim.sfs) ===\n")

sim.sfs(model = Vaquita2Epoch,
        nsim.blocks = 1,
        block.size = 10,
        path = tmpdir,
        use.alpha = FALSE,
        output.name = "test_murates",
        mu.rates = list("runif", 200, 1e-9, 1e-8),
        ncores = 1)

outfile_mu <- file.path(tmpdir, "SIM_SFS_test_murates.txt")
check("mu.rates: output file exists", file.exists(outfile_mu))

sim_mu <- read.table(outfile_mu, header = TRUE, sep = "\t")
check("mu.rates: 10 simulation rows", nrow(sim_mu) == 10)

sfs_cols_mu <- grep("^sfs_", colnames(sim_mu), value = TRUE)
check("mu.rates: 39 SFS columns", length(sfs_cols_mu) == 39)
check("mu.rates: all SFS values non-negative",
      all(as.matrix(sim_mu[, sfs_cols_mu]) >= 0))

file.remove(outfile_mu)


# ============================================================================
# Summary
# ============================================================================
cat(paste0("\n=== Results: ", pass, " passed, ", fail, " failed ===\n"))
if(fail > 0) {
  cat("SOME TESTS FAILED!\n")
  quit(status = 1)
} else {
  cat("ALL TESTS PASSED!\n")
}
