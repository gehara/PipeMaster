#!/usr/bin/env Rscript
# Test script for obs.sfs() function
# Validates observed SFS computation from PHYLIP data against Python reference SFS

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

# ============================================================================
# Test 1: Vaquita2Epoch (1-pop, 40 haplotypes, PHYLIP input)
# ============================================================================
cat("\n=== Test 1: Vaquita2Epoch (1-pop, PHYLIP) ===\n")

pop_assign_vaq <- read.table("pop_assign_Vaquita2Epoch.txt", header = FALSE)
obs_vaq <- obs.sfs(model = Vaquita2Epoch,
                   path.to.phylip = file.path(test_dir, "phylip_Vaquita2Epoch.phy"),
                   pop.assign = pop_assign_vaq)

check("Vaquita: 1 row", nrow(obs_vaq) == 1)
check("Vaquita: 41 columns (sfs_0..sfs_40)", ncol(obs_vaq) == 41)
check("Vaquita: column names start with sfs_", all(grepl("^sfs_", colnames(obs_vaq))))
check("Vaquita: first col is sfs_0", colnames(obs_vaq)[1] == "sfs_0")
check("Vaquita: last col is sfs_40", colnames(obs_vaq)[41] == "sfs_40")
check("Vaquita: all values non-negative", all(obs_vaq >= 0))
check("Vaquita: total SNPs > 0", sum(obs_vaq) > 0)

# Compare with Python reference SFS
ref_vaq <- read.table("observed_sfs_Vaquita2Epoch.txt", header = TRUE, sep = "\t")
check("Vaquita: same dimensions as reference", ncol(obs_vaq) == ncol(ref_vaq))
check("Vaquita: same column names as reference", all(colnames(obs_vaq) == colnames(ref_vaq)))

# Correlation test (major-allele polarization may differ from true ancestral)
obs_vec <- as.numeric(obs_vaq[1, ])
ref_vec <- as.numeric(ref_vaq[1, ])
# Use only non-zero bins for correlation (avoid inflating with shared zeros)
nonzero <- (obs_vec + ref_vec) > 0
if(sum(nonzero) > 2) {
  r <- cor(obs_vec[nonzero], ref_vec[nonzero])
  cat(paste0("  Vaquita SFS correlation (non-zero bins): r = ", round(r, 4), "\n"))
  check("Vaquita: SFS correlation > 0.9", r > 0.9)
} else {
  cat("  WARNING: too few non-zero bins to compute correlation\n")
}

cat(paste0("  Vaquita total SNPs (obs): ", sum(obs_vaq),
           ", (ref): ", sum(ref_vaq), "\n"))
check("Vaquita: total SNP count matches reference", sum(obs_vaq) == sum(ref_vaq))


# ============================================================================
# Test 2: Africa_1T12 (1-pop, 40 haplotypes, PHYLIP input)
# ============================================================================
cat("\n=== Test 2: Africa_1T12 (1-pop, PHYLIP) ===\n")

pop_assign_afr <- read.table("pop_assign_Africa_1T12.txt", header = FALSE)
obs_afr <- obs.sfs(model = Africa_1T12,
                   path.to.phylip = file.path(test_dir, "phylip_Africa_1T12.phy"),
                   pop.assign = pop_assign_afr)

check("Africa: 1 row", nrow(obs_afr) == 1)
check("Africa: 41 columns", ncol(obs_afr) == 41)
check("Africa: all values non-negative", all(obs_afr >= 0))
check("Africa: total SNPs > 0", sum(obs_afr) > 0)

ref_afr <- read.table("observed_sfs_Africa_1T12.txt", header = TRUE, sep = "\t")
check("Africa: same dimensions as reference", ncol(obs_afr) == ncol(ref_afr))

obs_vec <- as.numeric(obs_afr[1, ])
ref_vec <- as.numeric(ref_afr[1, ])
nonzero <- (obs_vec + ref_vec) > 0
if(sum(nonzero) > 2) {
  r <- cor(obs_vec[nonzero], ref_vec[nonzero])
  cat(paste0("  Africa SFS correlation (non-zero bins): r = ", round(r, 4), "\n"))
  check("Africa: SFS correlation > 0.9", r > 0.9)
}

check("Africa: total SNP count matches reference", sum(obs_afr) == sum(ref_afr))


# ============================================================================
# Test 3: OutOfAfrica_2T12 (2-pop joint SFS, 40+40 haplotypes, PHYLIP)
# ============================================================================
cat("\n=== Test 3: OutOfAfrica_2T12 (2-pop joint SFS, PHYLIP) ===\n")

pop_assign_ooa2 <- read.table("pop_assign_OutOfAfrica_2T12.txt", header = FALSE)
obs_ooa2 <- obs.sfs(model = OutOfAfrica_2T12,
                    path.to.phylip = file.path(test_dir, "phylip_OutOfAfrica_2T12.phy"),
                    pop.assign = pop_assign_ooa2)

check("OoA_2T12: 1 row", nrow(obs_ooa2) == 1)
check("OoA_2T12: 1681 columns (41*41)", ncol(obs_ooa2) == 1681)
check("OoA_2T12: first col is sfs_0_0", colnames(obs_ooa2)[1] == "sfs_0_0")
check("OoA_2T12: all values non-negative", all(obs_ooa2 >= 0))
check("OoA_2T12: total SNPs > 0", sum(obs_ooa2) > 0)

# Check sfs_matrix attribute
sfs_mat <- attr(obs_ooa2, "sfs_matrix")
check("OoA_2T12: has sfs_matrix attribute", !is.null(sfs_mat))
check("OoA_2T12: matrix is 41x41", all(dim(sfs_mat) == c(41, 41)))
check("OoA_2T12: matrix sum equals data frame sum", sum(sfs_mat) == sum(obs_ooa2))

# Compare with Python reference
ref_ooa2 <- read.table("observed_sfs_OutOfAfrica_2T12.txt", header = TRUE, sep = "\t")
check("OoA_2T12: same dimensions as reference", ncol(obs_ooa2) == ncol(ref_ooa2))
check("OoA_2T12: same column names as reference",
      all(colnames(obs_ooa2) == colnames(ref_ooa2)))

obs_vec <- as.numeric(obs_ooa2[1, ])
ref_vec <- as.numeric(ref_ooa2[1, ])
nonzero <- (obs_vec + ref_vec) > 0
if(sum(nonzero) > 2) {
  r <- cor(obs_vec[nonzero], ref_vec[nonzero])
  cat(paste0("  OoA_2T12 joint SFS correlation (non-zero bins): r = ", round(r, 4), "\n"))
  check("OoA_2T12: joint SFS correlation > 0.9", r > 0.9)
}

check("OoA_2T12: total SNP count matches reference", sum(obs_ooa2) == sum(ref_ooa2))

# Compare matrix with Python reference matrix
ref_mat <- as.matrix(read.table("observed_joint_sfs_matrix_OutOfAfrica_2T12.txt"))
check("OoA_2T12: matrix dimensions match reference", all(dim(sfs_mat) == dim(ref_mat)))


# ============================================================================
# Test 4: OutOfAfrica_3G09 (2-pop, 40+40 haplotypes, PHYLIP)
# ============================================================================
cat("\n=== Test 4: OutOfAfrica_3G09 (2-pop, PHYLIP) ===\n")

pop_assign_ooa3 <- read.table("pop_assign_OutOfAfrica_3G09.txt", header = FALSE)
obs_ooa3 <- obs.sfs(model = OutOfAfrica_3G09,
                    path.to.phylip = file.path(test_dir, "phylip_OutOfAfrica_3G09.phy"),
                    pop.assign = pop_assign_ooa3)

check("OoA_3G09: 1 row", nrow(obs_ooa3) == 1)
check("OoA_3G09: 1681 columns (41*41)", ncol(obs_ooa3) == 1681)
check("OoA_3G09: all values non-negative", all(obs_ooa3 >= 0))
check("OoA_3G09: total SNPs > 0", sum(obs_ooa3) > 0)
check("OoA_3G09: has sfs_matrix attribute", !is.null(attr(obs_ooa3, "sfs_matrix")))


# ============================================================================
# Test 5: Folded SFS (1-pop, PHYLIP)
# ============================================================================
cat("\n=== Test 5: Folded SFS (PHYLIP) ===\n")

obs_vaq_folded <- obs.sfs(model = Vaquita2Epoch,
                          path.to.phylip = file.path(test_dir, "phylip_Vaquita2Epoch.phy"),
                          pop.assign = pop_assign_vaq,
                          folded = TRUE)

check("Folded: 1 row", nrow(obs_vaq_folded) == 1)
# Folded SFS: nsam=40 -> 41 bins unfolded, fold_sfs: n=41, half=20,
# plus middle entry since 41 is odd -> 21 total
check("Folded: 21 columns", ncol(obs_vaq_folded) == 21)
check("Folded: column names start with sfs_fold_",
      all(grepl("^sfs_fold_", colnames(obs_vaq_folded))))
check("Folded: all values non-negative", all(obs_vaq_folded >= 0))


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
