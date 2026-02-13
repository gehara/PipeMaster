#!/usr/bin/env Rscript
# Create PipeMaster model objects matching the stdpopsim test data.
#
# Models:
#   1. Africa_1T12        - 1 pop, 3-epoch (expansion)             -> for sim.sfs
#   2. OutOfAfrica_2T12   - 2 pop, AFR+EUR with migration          -> for sim.sfs (joint SFS)
#   3. OutOfAfrica_3G09   - 3 pop, YRI+CEU+CHB with migration      -> for sim.sumstat + sim.sfs
#   4. Vaquita2Epoch      - 1 pop, bottleneck                      -> for sim.sumstat + sim.sfs
#   5. PonAbe_TwoSpecies  - 2 pop, isolation-with-migration+growth -> for sim.sfs + sim.sumstat
#
# All models use 10000 loci x 100bp, 20 diploid (40 haploid) per pop.
# Priors are set as uniform ranges bracketing the true stdpopsim values.

n_loci   <- 10000
locus_bp <- 100
n_dip    <- 20     # diploid individuals per pop
n_hap    <- 40     # haploid samples per pop

# ============================================================================
# Helper: build loci matrix (shared rate across loci)
# ============================================================================
make_loci <- function(n_loci, bp, mu_min, mu_max, inheritance=1) {
  loci <- matrix(NA, nrow=n_loci, ncol=6)
  for (i in 1:n_loci) {
    loci[i,] <- c(paste0("rate"), as.character(bp), as.character(inheritance),
                  as.character(mu_min), as.character(mu_max), "runif")
  }
  loci
}

# ============================================================================
# Helper: build I matrix (sample configuration)
# ============================================================================
make_I <- function(n_loci, npop, samples_per_pop) {
  ncols <- 3 + npop
  I <- matrix(NA, nrow=n_loci, ncol=ncols)
  for (i in 1:n_loci) {
    I[i,] <- c(paste0("locus",i), "-I", as.character(npop),
               as.character(samples_per_pop))
  }
  I
}

# ============================================================================
# Helper: build condition matrices
# ============================================================================
make_size_conds <- function(ne_names) {
  n <- length(ne_names)
  m <- matrix(NA, nrow=n, ncol=n)
  diag(m) <- "0"
  rownames(m) <- colnames(m) <- ne_names
  m
}

make_mig_conds <- function(mig_names) {
  n <- length(mig_names)
  if (n == 0) return(matrix(NA, nrow=0, ncol=0))
  m <- matrix(NA, nrow=n, ncol=n)
  rownames(m) <- colnames(m) <- mig_names
  m
}

make_time_conds <- function(time_names) {
  n <- length(time_names)
  if (n == 0) return(matrix(NA, nrow=0, ncol=0))
  m <- matrix(NA, nrow=n, ncol=n)
  diag(m) <- "0"
  rownames(m) <- colnames(m) <- time_names
  m
}

# ============================================================================
# 1. Africa_1T12: Single pop, 3-epoch
#    True params: Ne_present=432125 (growth), Ne_middle=14474, Ne_ancient=7310
#    Epoch boundaries: 205 gen, 5920 gen
# ============================================================================
cat("Building Africa_1T12 model...\n")

# For single pop with sim.sfs, we approximate the 3-epoch model:
#   Ne0.pop1 ~ current Ne (14474, end of growth phase)
#   Ne1.pop1 ~ ancestral Ne (7310)
#   The exponential growth in the most recent epoch can be captured with use.alpha=TRUE
Africa_1T12 <- list(
  loci = make_loci(n_loci, locus_bp, mu_min="1e-09", mu_max="5e-08"),
  I    = make_I(n_loci, npop=1, samples_per_pop=n_hap),
  flags = list(
    n = matrix(c(
      "Ne0.pop1", "-n", "1", "5000", "600000", "runif"
    ), nrow=1, ncol=6, byrow=TRUE),
    en = list(
      size = matrix(c(
        "Ne1.pop1", "-en", "1", "5000", "30000", "runif",
        "Ne2.pop1", "-en", "1", "1000", "20000", "runif"
      ), nrow=2, ncol=6, byrow=TRUE),
      time = matrix(c(
        "t.Ne1.pop1", "-en", "1", "50",  "500",   "runif",
        "t.Ne2.pop1", "-en", "1", "2000", "10000", "runif"
      ), nrow=2, ncol=6, byrow=TRUE)
    ),
    ej = NULL
  ),
  conds = list(
    size.matrix = make_size_conds(c("Ne0.pop1","Ne1.pop1","Ne2.pop1")),
    mig.matrix  = matrix(NA, nrow=0, ncol=0),
    time.matrix = {
      tm <- make_time_conds(c("t.Ne1.pop1","t.Ne2.pop1"))
      tm["t.Ne1.pop1","t.Ne2.pop1"] <- "<"
      tm["t.Ne2.pop1","t.Ne1.pop1"] <- ">"
      tm
    }
  ),
  tree = "(1);"
)
class(Africa_1T12) <- "Model"

# Store true parameter values as attribute for reference
attr(Africa_1T12, "true_params") <- list(
  Ne0.pop1 = 432125, Ne1.pop1 = 14474, Ne2.pop1 = 7310,
  t.Ne1.pop1 = 205, t.Ne2.pop1 = 5920,
  mutation_rate = 2.36e-08, model = "Africa_1T12"
)

# ============================================================================
# 2. OutOfAfrica_2T12: 2 pop (AFR + EUR) with migration
#    True params:
#      Ne_AFR_present=432125, Ne_AFR_middle=14474, Ne_AFR_ancient=7310
#      Ne_EUR_present=501436, Ne_EUR_920=9279, Ne_EUR_bot=1861
#      Split: 2040 gen, size changes at 205, 920, 5920
#      Migration: 2.5e-5 per gen (recent), 1.5e-4 (old)
# ============================================================================
cat("Building OutOfAfrica_2T12 model...\n")

OutOfAfrica_2T12 <- list(
  loci = make_loci(n_loci, locus_bp, mu_min="1e-09", mu_max="5e-08"),
  I    = make_I(n_loci, npop=2, samples_per_pop=c(n_hap, n_hap)),
  flags = list(
    n = matrix(c(
      "Ne0.pop1", "-n", "1", "5000",  "600000", "runif",
      "Ne0.pop2", "-n", "2", "5000",  "600000", "runif"
    ), nrow=2, ncol=6, byrow=TRUE),
    m = matrix(c(
      "mig0.1_2", "-m", "1 2", "0", "100", "runif",
      "mig0.2_1", "-m", "2 1", "0", "100", "runif"
    ), nrow=2, ncol=6, byrow=TRUE),
    en = list(
      size = matrix(c(
        "Ne1.pop1", "-en", "1", "5000",  "30000",  "runif",
        "Ne1.pop2", "-en", "2", "500",   "5000",   "runif",
        "Ne2.pop1", "-en", "1", "1000",  "20000",  "runif"
      ), nrow=3, ncol=6, byrow=TRUE),
      time = matrix(c(
        "t.Ne1.pop1", "-en", "1", "50",   "500",   "runif",
        "t.Ne1.pop2", "-en", "2", "400",  "2000",  "runif",
        "t.Ne2.pop1", "-en", "1", "2000", "10000", "runif"
      ), nrow=3, ncol=6, byrow=TRUE)
    ),
    ej = matrix(c(
      "join1", "-ej", "1 2", "1000", "5000", "runif"
    ), nrow=1, ncol=6, byrow=TRUE)
  ),
  conds = list(
    size.matrix = make_size_conds(c("Ne0.pop1","Ne0.pop2","Ne1.pop1","Ne1.pop2","Ne2.pop1")),
    mig.matrix  = make_mig_conds(c("mig0.1_2","mig0.2_1")),
    time.matrix = {
      tm <- make_time_conds(c("join1","t.Ne1.pop1","t.Ne1.pop2","t.Ne2.pop1"))
      # Ne1.pop1 (205 gen) < Ne1.pop2 (920 gen) < join1 (2040 gen) < Ne2.pop1 (5920 gen)
      tm["t.Ne1.pop1","t.Ne1.pop2"] <- "<"
      tm["t.Ne1.pop2","t.Ne1.pop1"] <- ">"
      tm["t.Ne1.pop2","join1"]      <- "<"
      tm["join1","t.Ne1.pop2"]      <- ">"
      tm["join1","t.Ne2.pop1"]      <- "<"
      tm["t.Ne2.pop1","join1"]      <- ">"
      tm
    }
  ),
  tree = "(1,2);"
)
class(OutOfAfrica_2T12) <- "Model"

attr(OutOfAfrica_2T12, "true_params") <- list(
  Ne0.pop1 = 432125, Ne0.pop2 = 501436,
  Ne1.pop1 = 14474, Ne1.pop2 = 1861, Ne2.pop1 = 7310,
  t.Ne1.pop1 = 205, t.Ne1.pop2 = 920, t.Ne2.pop1 = 5920,
  join1 = 2040,
  mig0.1_2 = 4 * 432125 * 2.5e-5,   # 4Nm = 43.21 (stdpopsim m = 2.5e-5)
  mig0.2_1 = 4 * 501436 * 2.5e-5,   # 4Nm = 50.14 (stdpopsim m = 2.5e-5)
  mutation_rate = 2.36e-08, model = "OutOfAfrica_2T12"
)

# ============================================================================
# 3. OutOfAfrica_3G09: Full 3 pop (YRI + CEU + CHB)
#    True params (Gutenkunst et al. 2009):
#      Ne_YRI=12300, Ne_CEU_present=29725 (growth), Ne_CHB_present=54090 (growth)
#      Ne_CEU_bot=1000, Ne_CHB_bot=510, Ne_ancestral=7300
#      Split YRI-OoA: 5600 gen, Split CEU-CHB: 848 gen
#      T_AF (ancestral Ne change): 8800 gen
#      Migration: YRI<->CEU, YRI<->CHB, CEU<->CHB
# ============================================================================
cat("Building OutOfAfrica_3G09 model (YRI+CEU+CHB)...\n")

OutOfAfrica_3G09 <- list(
  loci = make_loci(n_loci, locus_bp, mu_min="1e-09", mu_max="5e-08"),
  I    = make_I(n_loci, npop=3, samples_per_pop=c(n_hap, n_hap, n_hap)),
  flags = list(
    n = matrix(c(
      "Ne0.pop1", "-n", "1", "5000",  "25000",  "runif",
      "Ne0.pop2", "-n", "2", "5000",  "60000",  "runif",
      "Ne0.pop3", "-n", "3", "5000",  "100000", "runif"
    ), nrow=3, ncol=6, byrow=TRUE),
    m = matrix(c(
      "mig0.1_2", "-m", "1 2", "0", "8",  "runif",
      "mig0.2_1", "-m", "2 1", "0", "8",  "runif",
      "mig0.1_3", "-m", "1 3", "0", "5",  "runif",
      "mig0.3_1", "-m", "3 1", "0", "5",  "runif",
      "mig0.2_3", "-m", "2 3", "0", "20", "runif",
      "mig0.3_2", "-m", "3 2", "0", "20", "runif"
    ), nrow=6, ncol=6, byrow=TRUE),
    en = list(
      size = matrix(c(
        "Ne1.pop2", "-en", "2", "200",  "3000",  "runif",
        "Ne1.pop3", "-en", "3", "100",  "2000",  "runif",
        "Ne1.pop1", "-en", "1", "2000", "15000", "runif"
      ), nrow=3, ncol=6, byrow=TRUE),
      time = matrix(c(
        "t.Ne1.pop2", "-en", "2", "200",  "2000",  "runif",
        "t.Ne1.pop3", "-en", "3", "200",  "2000",  "runif",
        "t.Ne1.pop1", "-en", "1", "2000", "15000", "runif"
      ), nrow=3, ncol=6, byrow=TRUE)
    ),
    ej = matrix(c(
      "join2_3", "-ej", "2 3", "200",  "2000",  "runif",
      "join1",   "-ej", "1 3", "2000", "10000", "runif"
    ), nrow=2, ncol=6, byrow=TRUE)
  ),
  conds = list(
    size.matrix = make_size_conds(c("Ne0.pop1","Ne0.pop2","Ne0.pop3","Ne1.pop2","Ne1.pop3","Ne1.pop1")),
    mig.matrix  = make_mig_conds(c("mig0.1_2","mig0.2_1","mig0.1_3","mig0.3_1","mig0.2_3","mig0.3_2")),
    time.matrix = {
      tm <- make_time_conds(c("join2_3","join1","t.Ne1.pop2","t.Ne1.pop3","t.Ne1.pop1"))
      # CEU-CHB bottleneck times ~ split time (848 gen)
      # join2_3 (CEU-CHB split) ~ 848 gen
      # join1 (YRI-OoA split) ~ 5600 gen
      # t.Ne1.pop1 (ancestral) ~ 8800 gen
      tm["join2_3","join1"]         <- "<"
      tm["join1","join2_3"]         <- ">"
      tm["join1","t.Ne1.pop1"]      <- "<"
      tm["t.Ne1.pop1","join1"]      <- ">"
      # Bottleneck times must be <= join2_3 (they happen at split)
      tm["t.Ne1.pop2","join2_3"]    <- "<"
      tm["join2_3","t.Ne1.pop2"]    <- ">"
      tm["t.Ne1.pop3","join2_3"]    <- "<"
      tm["join2_3","t.Ne1.pop3"]    <- ">"
      tm
    }
  ),
  tree = "((1,2),3);"
)
class(OutOfAfrica_3G09) <- "Model"

# True values from Gutenkunst et al. 2009
# Ne_EU_present = 1000 * exp(0.004 * 848) = 29725
# Ne_AS_present = 510  * exp(0.0055 * 848) = 54090
attr(OutOfAfrica_3G09, "true_params") <- list(
  Ne0.pop1 = 12300, Ne0.pop2 = 29725, Ne0.pop3 = 54090,
  Ne1.pop2 = 1000, Ne1.pop3 = 510, Ne1.pop1 = 7300,
  t.Ne1.pop2 = 848, t.Ne1.pop3 = 848, t.Ne1.pop1 = 8800,
  join2_3 = 848, join1 = 5600,
  mig0.1_2 = 4 * 12300 * 3e-5,   # 4Nm = 1.476
  mig0.2_1 = 4 * 29725 * 3e-5,   # 4Nm = 3.567
  mig0.1_3 = 4 * 12300 * 1.9e-5, # 4Nm = 0.935
  mig0.3_1 = 4 * 54090 * 1.9e-5, # 4Nm = 4.111
  mig0.2_3 = 4 * 29725 * 9.6e-5, # 4Nm = 11.41
  mig0.3_2 = 4 * 54090 * 9.6e-5, # 4Nm = 20.77
  mutation_rate = 2.35e-08, model = "OutOfAfrica_3G09"
)

# ============================================================================
# 4. Vaquita2Epoch: 1 pop, simple bottleneck
#    True params: Ne_present=2807, Ne_ancient=4485, time=2162 gen
# ============================================================================
cat("Building Vaquita2Epoch model...\n")

Vaquita2Epoch <- list(
  loci = make_loci(n_loci, locus_bp, mu_min="1e-10", mu_max="1e-08"),
  I    = make_I(n_loci, npop=1, samples_per_pop=n_hap),
  flags = list(
    n = matrix(c(
      "Ne0.pop1", "-n", "1", "500", "10000", "runif"
    ), nrow=1, ncol=6, byrow=TRUE),
    en = list(
      size = matrix(c(
        "Ne1.pop1", "-en", "1", "1000", "15000", "runif"
      ), nrow=1, ncol=6, byrow=TRUE),
      time = matrix(c(
        "t.Ne1.pop1", "-en", "1", "500", "5000", "runif"
      ), nrow=1, ncol=6, byrow=TRUE)
    ),
    ej = NULL
  ),
  conds = list(
    size.matrix = {
      sm <- make_size_conds(c("Ne0.pop1","Ne1.pop1"))
      # Present Ne < ancestral Ne (bottleneck: decline going forward)
      sm["Ne0.pop1","Ne1.pop1"] <- "<"
      sm["Ne1.pop1","Ne0.pop1"] <- ">"
      sm
    },
    mig.matrix  = matrix(NA, nrow=0, ncol=0),
    time.matrix = make_time_conds(c("t.Ne1.pop1"))
  ),
  tree = "(1);"
)
class(Vaquita2Epoch) <- "Model"

attr(Vaquita2Epoch, "true_params") <- list(
  Ne0.pop1 = 2807, Ne1.pop1 = 4485,
  t.Ne1.pop1 = 2162,
  mutation_rate = 5.83e-09, model = "Vaquita2Epoch_1R22"
)

# ============================================================================
# 5. PonAbe TwoSpecies_2L11: 2 pop isolation-with-migration + growth
#    True params (Locke et al. 2011):
#      Na=17934, N_S(present)=37661, N_B(present)=8805
#      N_S(at split)=7299, N_B(at split)=10635
#      Split: 20157 gen
#      Growth: r_S=0.0000453/gen, r_B=-0.0000828/gen
#      Migration: m_S_B=1.099e-5, m_B_S=6.646e-6
# ============================================================================
cat("Building PonAbe TwoSpecies_2L11 model...\n")

PonAbe_TwoSpecies <- list(
  loci = make_loci(n_loci, locus_bp, mu_min="1e-09", mu_max="5e-08"),
  I    = make_I(n_loci, npop=2, samples_per_pop=c(n_hap, n_hap)),
  flags = list(
    n = matrix(c(
      "Ne0.pop1", "-n", "1", "10000", "80000", "runif",
      "Ne0.pop2", "-n", "2", "2000",  "25000", "runif"
    ), nrow=2, ncol=6, byrow=TRUE),
    m = matrix(c(
      "mig0.1_2", "-m", "1 2", "0", "5", "runif",
      "mig0.2_1", "-m", "2 1", "0", "5", "runif"
    ), nrow=2, ncol=6, byrow=TRUE),
    en = list(
      size = matrix(c(
        "Ne1.pop1", "-en", "1", "5000", "40000", "runif"
      ), nrow=1, ncol=6, byrow=TRUE),
      time = matrix(c(
        "t.Ne1.pop1", "-en", "1", "5000", "50000", "runif"
      ), nrow=1, ncol=6, byrow=TRUE)
    ),
    ej = matrix(c(
      "join1", "-ej", "1 2", "5000", "50000", "runif"
    ), nrow=1, ncol=6, byrow=TRUE)
  ),
  conds = list(
    size.matrix = make_size_conds(c("Ne0.pop1","Ne0.pop2","Ne1.pop1")),
    mig.matrix  = make_mig_conds(c("mig0.1_2","mig0.2_1")),
    time.matrix = {
      tm <- make_time_conds(c("join1","t.Ne1.pop1"))
      # Ne change at split time (they're the same event in Locke model)
      tm["t.Ne1.pop1","join1"] <- "<"
      tm["join1","t.Ne1.pop1"] <- ">"
      tm
    }
  ),
  tree = "(1,2);"
)
class(PonAbe_TwoSpecies) <- "Model"

# True values from Locke et al. 2011
# Pop 1 = Sumatran (grows from 7299 to 37661 over 20157 gen)
# Pop 2 = Bornean (declines from 10635 to 8805 over 20157 gen)
# Ne1.pop1 = ancestral Na = 17934 (after join)
attr(PonAbe_TwoSpecies, "true_params") <- list(
  Ne0.pop1 = 37661, Ne0.pop2 = 8805,
  Ne1.pop1 = 17934,
  join1 = 20157, t.Ne1.pop1 = 20157,
  mig0.1_2 = 4 * 37661 * 1.099e-5,  # 4Nm = 1.656 (Sumatran -> Bornean)
  mig0.2_1 = 4 * 8805 * 6.646e-6,   # 4Nm = 0.234 (Bornean -> Sumatran)
  mutation_rate = 2e-08, model = "TwoSpecies_2L11"
)

# ============================================================================
# Save all models and observed data
# ============================================================================
cat("Loading observed data...\n")

observed_sfs_Africa_1T12 <- tryCatch(
  read.table("observed_sfs_Africa_1T12.txt", header=TRUE, sep="\t"),
  error = function(e) { cat("  (Africa_1T12 SFS not yet generated)\n"); NULL })

observed_sfs_OutOfAfrica_2T12 <- tryCatch(
  read.table("observed_sfs_OutOfAfrica_2T12.txt", header=TRUE, sep="\t"),
  error = function(e) { cat("  (OutOfAfrica_2T12 SFS not yet generated)\n"); NULL })

observed_joint_sfs_matrix_OutOfAfrica_2T12 <- tryCatch(
  as.matrix(read.table("observed_joint_sfs_matrix_OutOfAfrica_2T12.txt", header=FALSE, sep="\t")),
  error = function(e) { cat("  (OutOfAfrica_2T12 joint SFS matrix not yet generated)\n"); NULL })

observed_sumstats_OutOfAfrica_3G09 <- tryCatch(
  read.table("observed_sumstats_OutOfAfrica_3G09.txt", header=TRUE, sep="\t"),
  error = function(e) { cat("  (OutOfAfrica_3G09 sumstats not yet generated)\n"); NULL })

observed_sumstats_Vaquita2Epoch <- tryCatch(
  read.table("observed_sumstats_Vaquita2Epoch.txt", header=TRUE, sep="\t"),
  error = function(e) { cat("  (Vaquita2Epoch sumstats not yet generated)\n"); NULL })

observed_sfs_Vaquita2Epoch <- tryCatch(
  read.table("observed_sfs_Vaquita2Epoch.txt", header=TRUE, sep="\t"),
  error = function(e) { cat("  (Vaquita2Epoch SFS not yet generated)\n"); NULL })

# PonAbe observed data
observed_sfs_PonAbe <- tryCatch(
  read.table("observed_sfs_PonAbe.txt", header=TRUE, sep="\t"),
  error = function(e) { cat("  (PonAbe SFS not yet generated)\n"); NULL })

observed_sumstats_PonAbe <- tryCatch(
  read.table("observed_sumstats_PonAbe.txt", header=TRUE, sep="\t"),
  error = function(e) { cat("  (PonAbe sumstats not yet generated)\n"); NULL })

# 3-pop OoA SFS
observed_sfs_OutOfAfrica_3G09 <- tryCatch(
  read.table("observed_sfs_OutOfAfrica_3G09.txt", header=TRUE, sep="\t"),
  error = function(e) { cat("  (3-pop OoA SFS not yet generated)\n"); NULL })

cat("Saving to test_models.RData...\n")
save_list <- c("Africa_1T12", "OutOfAfrica_2T12", "OutOfAfrica_3G09",
               "Vaquita2Epoch", "PonAbe_TwoSpecies")
obs_objects <- c("observed_sfs_Africa_1T12", "observed_sfs_OutOfAfrica_2T12",
                 "observed_joint_sfs_matrix_OutOfAfrica_2T12",
                 "observed_sumstats_OutOfAfrica_3G09",
                 "observed_sumstats_Vaquita2Epoch", "observed_sfs_Vaquita2Epoch",
                 "observed_sfs_PonAbe", "observed_sumstats_PonAbe",
                 "observed_sfs_OutOfAfrica_3G09")
for (obj in obs_objects) {
  if (!is.null(get(obj))) save_list <- c(save_list, obj)
}

save(list = save_list, file = "test_models.RData")

cat("\nDone! Saved test_models.RData with:\n")
cat("  Models: Africa_1T12, OutOfAfrica_2T12, OutOfAfrica_3G09, Vaquita2Epoch, PonAbe_TwoSpecies\n")
cat("  Observed: observed_sfs_Africa_1T12, observed_sfs_OutOfAfrica_2T12,\n")
cat("            observed_joint_sfs_matrix_OutOfAfrica_2T12,\n")
cat("            observed_sumstats_OutOfAfrica_3G09, observed_sfs_OutOfAfrica_3G09,\n")
cat("            observed_sumstats_Vaquita2Epoch, observed_sfs_Vaquita2Epoch,\n")
cat("            observed_sfs_PonAbe, observed_sumstats_PonAbe\n")

# Print summaries
cat("\n=== Model summaries ===\n")
for (mname in c("Africa_1T12","OutOfAfrica_2T12","OutOfAfrica_3G09","Vaquita2Epoch","PonAbe_TwoSpecies")) {
  m <- get(mname)
  cat(sprintf("\n%s:\n", mname))
  cat(sprintf("  Tree: %s\n", m$tree))
  cat(sprintf("  Loci: %d x %s bp\n", nrow(m$loci), m$loci[1,2]))
  cat(sprintf("  Pops: %s, samples: %s\n", m$I[1,3],
              paste(m$I[1,4:ncol(m$I)], collapse=",")))
  cat(sprintf("  Ne params: %s\n", paste(m$flags$n[,1], collapse=", ")))
  if (!is.null(m$flags$m)) cat(sprintf("  Mig params: %s\n", paste(m$flags$m[,1], collapse=", ")))
  if (!is.null(m$flags$en$size)) cat(sprintf("  En size: %s\n", paste(m$flags$en$size[,1], collapse=", ")))
  if (!is.null(m$flags$en$time)) cat(sprintf("  En time: %s\n", paste(m$flags$en$time[,1], collapse=", ")))
  if (!is.null(m$flags$ej)) cat(sprintf("  Join: %s\n", paste(m$flags$ej[,1], collapse=", ")))
  tp <- attr(m, "true_params")
  cat(sprintf("  True values: %s\n",
              paste(names(tp), "=", tp, collapse=", ")))
}
