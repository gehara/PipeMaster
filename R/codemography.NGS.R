#' Simulation of codemographic models with ngs data
#' @description Simulation of codemographic models. Accepts either an hModel object
#'   (created by \code{hModel()} or \code{h.menu.gui()}) or individual parameters.
#' @param nsims_or_hmodel Either: (1) an hModel object, in which case \code{nsims} is required;
#'   or (2) a numeric value specifying the total number of simulations (legacy mode).
#' @param nsims Total number of simulations (required when first argument is an hModel).
#' @param var.zeta Variation on zeta parameter. Can be "FREE" to vary or be set to a specific value (between 0-1). Overrides hModel value if provided.
#' @param coexp.prior Uniform prior for the coexpansion time. Vector of two numbers with the lower and upper boundary of the prior. Overrides hModel value if provided.
#' @param th Threshold. Minimum time difference between Ts, time of simultaneous change and population specific times. Overrides hModel value if provided.
#' @param Ne.prior Data frame with the prior values for the Ne of each population. Overrides hModel value if provided.
#' @param NeA.prior Data frame with the prior values for the ancestral Ne of each population. Overrides hModel value if provided.
#' @param time.prior Data frame with parameter values for the priors of the time of demographic change of each population. Overrides hModel value if provided.
#' @param gene.prior Data frame with parameter values for the priors of the mutation rate of each species. Overrides hModel value if provided.
#' @param alpha logical. If TRUE all demographic changes are exponential. If FALSE sudden changes. Overrides hModel value if provided.
#' @param append.sims logical. If TRUE simulations are appended to the simulations file. Default is FALSE.
#' @param mu.rates list. A list of parameters of mutation rate prior distribution. First parameter: the type of distribution. Second element should be NA. Third and remaining elements are parameters of the distribution which should be given in order. For instance if using the rnorm distribution the third parameter is the mean and the fourth the SD. Overrides hModel value if provided.
#' @param block.size Number of simulations to process per batch. Default is 100.
#' @param path Path to the directory to write the simulations. Default is the working directory.
#' @details To simulate the model of Chan et al. (2014) the th parameter should be set to zero and the time.prior should have the same value of the coexp.prior.
#' @details To simulate the Threshold model the th argument need to be higher than zero. To simulate the Narrow Coexpansion Time model the th argument need to be higher than zero and the boundaries of coexp.time should be narrower than the time.prior values.
#' See references for more details. Use the sim.coexp2 function to simulate the partitioned time model.
#' @references Gehara M., Garda A.A., Werneck F.P. et al. (2017) Estimating synchronous demographic changes across populations using hABC and its application for a herpetological community from northeastern Brazil. Molecular Ecology, 26, 4756--4771.
#' @references Chan Y.L., Schanzenbach D., & Hickerson M.J. (2014) Detecting concerted demographic response across community assemblages using hierarchical approximate Bayesian computation. Molecular Biology and Evolution, 31, 2501--2515.
#' @noRd
sim.coexp.ngs<-function(nsims_or_hmodel,
                    nsims = NULL,
                    var.zeta = NULL,
                    coexp.prior = NULL,
                    th = NULL,
                    Ne.prior = NULL,
                    NeA.prior = NULL,
                    time.prior = NULL,
                    gene.prior = NULL,
                    mu.rates = NULL,
                    alpha = NULL,
                    model = NULL,
                    append.sims=FALSE,
                    block.size=100,
                    path=getwd())
{
  # --- Dispatch: hModel vs legacy numeric ---
  if (inherits(nsims_or_hmodel, "hModel")) {
    hm <- nsims_or_hmodel
    if (is.null(nsims)) stop("'nsims' is required when using an hModel object")
    # Extract from hModel, allow explicit overrides
    if (is.null(var.zeta))   var.zeta   <- hm$var.zeta
    if (is.null(coexp.prior)) coexp.prior <- hm$coexp.prior
    if (is.null(th))         th         <- hm$th
    if (is.null(Ne.prior))   Ne.prior   <- hm$Ne.prior
    if (is.null(NeA.prior))  NeA.prior  <- hm$NeA.prior
    if (is.null(time.prior)) time.prior <- hm$time.prior
    if (is.null(gene.prior)) gene.prior <- hm$gene.prior
    if (is.null(mu.rates))   mu.rates   <- hm$mu.rates
    if (is.null(alpha))      alpha      <- hm$alpha
    if (is.null(model))      model      <- if (!is.null(hm$model)) hm$model else "single"
  } else if (is.numeric(nsims_or_hmodel)) {
    # Legacy mode: first arg is nsims
    nsims <- nsims_or_hmodel
    if (is.null(var.zeta) || is.null(coexp.prior) || is.null(th) ||
        is.null(Ne.prior) || is.null(NeA.prior) || is.null(time.prior) ||
        is.null(gene.prior) || is.null(mu.rates))
      stop("All prior arguments are required in legacy mode (when first argument is numeric nsims)")
    if (is.null(alpha)) alpha <- FALSE
    if (is.null(model)) model <- "single"
  } else {
    stop("First argument must be an hModel object or a numeric nsims value")
  }
  WD <- getwd()
  setwd(path)
  on.exit(setwd(WD))

  ## --- Column header for the simulation output ---------------------------
  ## Panel: zeta, ts, E(t), DI                                  (4 hyperparams)
  ## + per-species 4 within-species moments x 7 stats           (28 features/species)
  ## + 8 across-species summaries                               x 28 = 224 hyper-stats
  ## Naming:
  ##   within (msABC native): s_{mean,var,skew,kurt}_{segs,pi,thetaW,tajd,ZnS,nhap,Hd}
  ##   across (this aggregator): h_{mean,var,q25,q50,q75,mad,skew,kurt}_<within>
  WITHIN_STATS    <- c("segs", "pi", "thetaW", "tajd", "ZnS", "nhap", "Hd")
  WITHIN_MOMENTS  <- c("mean", "var", "skew", "kurt")
  ACROSS_FEATURES <- c("hmean", "hvar", "hq25", "hq50", "hq75", "hmad", "hskew", "hkurt")
  within_names <- as.vector(outer(WITHIN_MOMENTS, WITHIN_STATS,
                                  function(m, s) paste0("s_", m, "_", s)))
  hyper_names  <- as.vector(outer(ACROSS_FEATURES, within_names,
                                  function(a, w) paste0(a, "_", w)))
  N_HYPER <- length(hyper_names)   # 8 * 28 = 224

  if(!append.sims){
    nsp_hdr <- nrow(Ne.prior)
    sp_lbl  <- if (!is.null(Ne.prior$species)) as.character(Ne.prior$species)
               else as.character(Ne.prior[, 1])
    if (length(sp_lbl) != nsp_hdr) sp_lbl <- sprintf("pop_%d", seq_len(nsp_hdr))
    t_cols <- paste0("t_", sp_lbl)
    header_row <- c("K", "zeta", "ts", "E(t)", "DI", t_cols, hyper_names)
    write.table(matrix(header_row, nrow = 1), file = "simulations.txt",
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  }

  nsim.blocks <- ceiling(nsims / block.size)

  ## LHS-stratify K for the bins model: pre-shuffle a vector of length nsims
  ## that contains balanced counts of each K in {1..nspecies}. Each block
  ## consumes the next bs entries. For non-bins models, K is unused at the
  ## sampler level (zeta is sampled instead).
  nspecies <- nrow(Ne.prior)
  K_lhs <- if (model == "bins") {
    sample(rep(seq_len(nspecies), length.out = nsims))
  } else NULL

  cat(sprintf("PipeMaster:: sim.coexp.ngs model=%s\n", model))
  for(j in 1:nsim.blocks){
    bs <- min(block.size, nsims - (j-1) * block.size)
    K_block <- if (!is.null(K_lhs)) K_lhs[((j-1) * block.size + 1):((j-1) * block.size + bs)] else NULL

    x <- coexp.sample.pars.msABC(nruns=bs, var.zeta=var.zeta, coexp.prior=coexp.prior, th=th,
                                  Ne.prior=Ne.prior, NeA.prior=NeA.prior, time.prior=time.prior,
                                  gene.prior=gene.prior, model=model,
                                  K_block=K_block)

    y <- coexp.msABC.batch(MS.par=x$MS.par, gene.prior=gene.prior, alpha=alpha,
                           pop.par=x$pop.par, mu.rates=mu.rates)

    simulations <- cbind(x$coexp.par, y)

    write.table(simulations, file="simulations.txt", quote=FALSE, row.names=FALSE,
                col.names=FALSE, append=TRUE, sep="\t")

    total.sims <- min(j * block.size, nsims)
    cat(sprintf("PipeMaster:: %d/%d sims\n", total.sims, nsims))
  }
}

# Internal function of the sim.coexp function
# @description sample the parameters of the models from prior distributions.
# @export
coexp.sample.pars.msABC<-function(nruns,
                            var.zeta,
                            coexp.prior,
                            th,
                            Ne.prior,
                            NeA.prior,
                            time.prior,
                            gene.prior,
                            model = c("single", "PT", "bins"),
                            K_block = NULL) {

  model <- match.arg(model)

  MS.par<-list(NULL)
  pop.par<-list(NULL)
  nspecies<-nrow(Ne.prior)
  ## coexp.par columns: K, zeta, Ts, E(t), DI, t_pop_1, ..., t_pop_nspecies
  ##   K, zeta, Ts, E(t), DI: hyperparams (5)
  ##   t_pop_i              : per-species realized time (nspecies)
  ## Per-species times let downstream ABC pool times across accepted sims
  ## for multimodal "pulse" detection in the time posterior.
  coexp.par<-matrix(nrow=nruns,ncol=5L + nspecies)

  for(i in 1:nspecies){
    mat<-matrix(nrow=nruns,ncol=4)
    MS.par[[i]]<-mat
    pop.par[[i]]<-mat
  }

  for(j in 1:nruns) {
    time.prior.B<-time.prior

    ## --- Time assignment per species --------------------------------------
    ## single: Ts uniform from coexp.prior; non-coexp species get times from
    ##         their own per-species time.prior.
    ## PT    : divide coexp.prior into nspecies equal sub-bins, draw one time
    ##         from each (the e.t pool), pick Ts from the pool, non-coexp
    ##         species share the REMAINING pool entries (sampled without
    ##         replacement). Allows multiple distinct shared times.
    ## bins  : K-uniform partition model (Variant 3b) with LHS-stratified K
    ##   and UNIFORM (linear) time bins (flat density over coexp.prior).
    ##   1. Sample K from the LHS strata (K_block if supplied, else uniform).
    ##   2. Pick K time bins (without replacement) from nspecies equal-width
    ##      sub-bins of coexp.prior; one time uniformly from each. Uniform
    ##      bins -> flat prior density on time, so any posterior clustering
    ##      is data-driven, not prior-induced.
    ##   3. Sample K group sizes uniformly from compositions of nspecies
    ##      into K positive parts (stars-and-bars on nspecies-1 gaps).
    ##   4. Randomly assign species to groups respecting sizes.
    ##   K is reported uniformly on {1..nspecies}. zeta is the fraction of
    ##   species in the largest group; Ts is that group's shared time.
    if (model == "bins") {
      K <- if (!is.null(K_block)) K_block[j] else sample.int(nspecies, 1L)

      range <- coexp.prior[2] - coexp.prior[1]
      step  <- range / nspecies
      cuts  <- coexp.prior[1] + step * (0:nspecies)
      chosen_bins <- sort(sample.int(nspecies, K))
      bin_times   <- vapply(chosen_bins,
                            function(k) runif(1, cuts[k], cuts[k + 1L]),
                            numeric(1))

      ## Uniform sample of K positive group sizes summing to nspecies
      if (K == 1L) {
        group_sizes <- nspecies
      } else {
        divider_pos <- sort(sample.int(nspecies - 1L, K - 1L))
        group_sizes <- diff(c(0L, divider_pos, nspecies))
      }

      ## Random species-to-group assignment respecting sizes
      group_indices  <- rep.int(seq_len(K), group_sizes)
      species_perm   <- sample.int(nspecies)
      species_to_grp <- integer(nspecies)
      species_to_grp[species_perm] <- group_indices

      for (i in seq_len(nspecies)) {
        t <- bin_times[species_to_grp[i]]
        time.prior.B[i, c(3, 4)] <- t
      }

      zeta.b <- max(group_sizes)
      zeta   <- zeta.b / nspecies
      Ts     <- bin_times[which.max(group_sizes)]

    } else {
      ## single or PT: sample zeta (number of co-changing species) first.
      if (var.zeta == "FREE") {
        zeta.space <- (1:nspecies) / nspecies
        zeta       <- sample(zeta.space, 1)
        zeta.b     <- nspecies * zeta
      } else {
        zeta   <- var.zeta
        zeta.b <- nspecies * zeta
      }

    if (model == "PT") {
      ## Build the e.t pool: one time per sub-bin of coexp.prior
      range  <- coexp.prior[2] - coexp.prior[1]
      step   <- range / nspecies
      cuts   <- coexp.prior[1] + step * (0:nspecies)
      e.t    <- vapply(seq_len(nspecies),
                       function(k) runif(1, cuts[k], cuts[k + 1]),
                       numeric(1))
      ## Pick Ts from the pool; remove it from the remaining pool
      Ts_idx <- sample.int(nspecies, 1)
      Ts     <- e.t[Ts_idx]
      pool   <- e.t[-Ts_idx]

      if (zeta.b == nspecies) {
        time.prior.B[, 3] <- Ts
      } else {
        coexp.sp <- sort(sample.int(nspecies, zeta.b))
        time.prior.B[coexp.sp, c(3, 4)] <- Ts
        ## non-coexp species draw without replacement from the pool
        non_coexp <- setdiff(seq_len(nspecies), coexp.sp)
        picks     <- sample(pool, length(non_coexp))
        time.prior.B[non_coexp, c(3, 4)] <- picks
      }

    } else {  ## "single"
      if (zeta.b == 1) {
        for (u in 1:nrow(time.prior.B)) {
          time.prior.B[u, 3] <- runif(1, time.prior.B[u, 3], time.prior.B[u, 4])
        }
        Ts <- runif(1, coexp.prior[1], coexp.prior[2])
      } else {
        Ts <- runif(1, coexp.prior[1], coexp.prior[2])
        coexp.sp <- sort(sample.int(nspecies, zeta.b))
        time.prior.B[coexp.sp, c(3, 4)] <- Ts
        for (u in 1:nrow(time.prior.B)) {
          if (!(u %in% coexp.sp)) {
            time <- runif(1, time.prior.B[u, 3], time.prior.B[u, 4])
            while (abs(time - Ts) < th) {
              time <- runif(1, time.prior.B[u, 3], time.prior.B[u, 4])
            }
            time.prior.B[u, 3:4] <- time
          }
        }
      }
    }

      ## K for single/PT: deterministic from zeta.b
      ##   1 main shared time + (nspecies - zeta.b) other distinct times
      K <- nspecies - zeta.b + 1L
    }  ## end of single/PT else branch

    Et<-mean(time.prior.B[,3])
    Disp.index<-var(time.prior.B[,3])/Et

    coexp.par[j,]<-c(K, zeta, Ts, Et, Disp.index,
                     as.numeric(time.prior.B[, 3]))

    for(i in 1:nspecies){
      ms.par<-NULL

      Ne <- runif(1, Ne.prior[i,3], Ne.prior[i,4])
      Ne.EXP.t <- time.prior.B[i,3]/time.prior[i,5] #corrects by generations
      theta.A.ratio <- runif(1, NeA.prior[i,3], NeA.prior[i,4])# thetaA (NeA) ratio
      NeA <- Ne*theta.A.ratio
      scalar = 4*Ne # generates the coalescent scalar
      EXP.time=Ne.EXP.t/scalar # expansion time / 4N

      g.rate=-log(NeA/Ne)/Ne.EXP.t

      ms.par<-cbind(0,EXP.time,theta.A.ratio,g.rate)
      po.par<-c(Ne,time.prior.B[i,3],NeA,0)
      MS.par[[i]][j,]<-ms.par
      pop.par[[i]][j,]<-po.par
    }
    #print(j)
  }
  #write.table(coexp.par,file="sim.par.txt", quote=F,row.names=F, col.names=F, append=T, sep="\t")
  #return(MS.par)
  pars<-list(NULL,NULL,NULL)
  names(pars)<-c("coexp.par","MS.par","pop.par")
  pars$coexp.par<-coexp.par
  pars$MS.par<-MS.par
  pars$pop.par<-pop.par
  return(pars)
}

# internal function
# @description Batch version of coexp.msABC using msABC_batch_call
# @return matrix (block_size x 24) of hypersummary statistics
coexp.msABC.batch<-function(MS.par, gene.prior, alpha, pop.par, mu.rates) {

  nspecies <- length(MS.par)
  block_size <- nrow(MS.par[[1]])

  # Per-species batch simulation
  species_stats <- vector("list", nspecies)

  for(sp in 1:nspecies) {
    locfile <- gene.prior[[sp]]
    nloci <- nrow(locfile)
    nsamp <- gene.prior[[sp]][1, 2]

    # Per-species mu.rates: support both list-of-lists (hModel) and single list (legacy)
    if (is.list(mu.rates[[1]])) {
      sp_mu <- mu.rates[[sp]]
    } else {
      sp_mu <- mu.rates
    }

    # Build commands and mu_mat for entire block
    commands <- character(block_size)
    mu_mat <- matrix(0, nrow=nloci, ncol=block_size)

    for(i in 1:block_size) {
      # Sample mu rates for this sim using this species' distribution
      sp_mu[[2]] <- nloci
      mu_mat[, i] <- do.call(sp_mu[[1]], args=sp_mu[2:length(sp_mu)])

      # Build command string
      if(alpha) {
        commands[i] <- paste(nsamp, 1,
                             "-G", MS.par[[sp]][i, 4],
                             "-eG", MS.par[[sp]][i, 2], 0,
                             "-eN", MS.par[[sp]][i, 2], MS.par[[sp]][i, 3],
                             "--frag-begin --finp .1locfile.txt --N", pop.par[[sp]][i, 1],
                             "--frag-end")
      } else {
        commands[i] <- paste(nsamp, 1,
                             "-eN", MS.par[[sp]][i, 2], MS.par[[sp]][i, 3],
                             "--frag-begin --finp .1locfile.txt --N", pop.par[[sp]][i, 1],
                             "--frag-end")
      }
    }

    # Write locfile for msABC (needs mu and rec columns even though mu is overridden)
    msabc_locfile <- data.frame(locfile, mu = 0, rec = 0, stringsAsFactors = FALSE)
    write.table(msabc_locfile, ".1locfile.txt", row.names=FALSE, col.names=TRUE,
                quote=FALSE, sep=" ")

    # Single batch C call for all sims of this species
    outputs <- .Call("msABC_batch_call", commands, mu_mat, NULL,
                     PACKAGE = "PipeMaster")

    # Parse outputs -> extract s_mean_* columns
    ## Keep all 4 within-species moments (s_mean_, s_var_, s_skew_, s_kurt_)
    ## for the 7 retained stats: segs, pi, thetaW, tajd, ZnS, nhap, Hd.
    ## Drops thomson_est, thomson_var, FayWuH.
    ## Returns a sims-block-row x 28 matrix per species.
    species_stats[[sp]] <- t(sapply(strsplit(outputs, "\n", fixed=TRUE), function(lines) {
      tab <- read.table(text=lines, header=TRUE, sep="\t")
      cols_rm <- c(grep("thomson", colnames(tab)),
                   grep("FayWuH",  colnames(tab)))
      if(length(cols_rm) > 0) tab <- tab[, -cols_rm]
      moment_cols <- grep("^s_(mean|var|skew|kurt)_", colnames(tab))
      tab <- tab[, moment_cols]
      as.numeric(tab)
    }))
  }

  # Clean up locfile
  if(file.exists(".1locfile.txt")) file.remove(".1locfile.txt")

  ## Compute hypersummary for each sim across species.
  ## sp_mat is (nspecies x 28). Apply 8 across-species summaries -> 224 hyper-stats.
  n_features <- ncol(species_stats[[1]])
  n_hyper    <- 8L * n_features   # 8 across-species summaries
  result <- matrix(0, nrow=block_size, ncol=n_hyper)

  ## Detect mad availability (base::mad exists; e1071::kurtosis/skewness used elsewhere here)
  for(i in 1:block_size) {
    sp_mat <- do.call(rbind, lapply(species_stats, function(m) m[i, ]))

    hmean <- colMeans(sp_mat, na.rm=TRUE)
    hvar  <- apply(sp_mat, 2, stats::var, na.rm=TRUE)
    hq25  <- apply(sp_mat, 2, function(x) stats::quantile(x, probs=0.25, na.rm=TRUE))
    hq50  <- apply(sp_mat, 2, function(x) stats::quantile(x, probs=0.50, na.rm=TRUE))
    hq75  <- apply(sp_mat, 2, function(x) stats::quantile(x, probs=0.75, na.rm=TRUE))
    hmad  <- apply(sp_mat, 2, function(x) stats::mad(x, na.rm=TRUE))
    hskew <- apply(sp_mat, 2, e1071::skewness, na.rm=TRUE)
    hkurt <- apply(sp_mat, 2, e1071::kurtosis, na.rm=TRUE)

    result[i, ] <- c(hmean, hvar, hq25, hq50, hq75, hmad, hskew, hkurt)
  }

  result
}

# internal function (legacy, kept for reference)
# @description control ms simulations - single-call-per-sim version
# @return ms simulations
coexp.msABC.legacy<-function(MS.par, gene.prior, alpha, pop.par, mu.rates) {

  nspecies <- length(MS.par)
  nruns<-nrow(MS.par[[1]])
  sim<-list(NULL)

  for(ii in 1:nruns){

    for(xx in 1:nspecies) {

      locfile <- gene.prior[[xx]]

      # Per-species mu.rates: support both list-of-lists (hModel) and single list (legacy)
      if (is.list(mu.rates[[1]])) {
        sp_mu <- mu.rates[[xx]]
      } else {
        sp_mu <- mu.rates
      }

      sp_mu[[2]] <- nrow(locfile)
      mu_vals <- do.call(sp_mu[[1]], args=sp_mu[2:length(sp_mu)])

      # Build locfile with mu and rec columns for msABC
      msabc_locfile <- data.frame(locfile, mu = mu_vals, rec = 0, stringsAsFactors = FALSE)

      write.table(msabc_locfile, paste(".locfile.txt", sep=""), row.names = F, col.names = T, quote = F, sep=" ")

      sims <- read.table(text=run.msABC(paste(gene.prior[[xx]][1,2],1,"-eN",MS.par[[xx]][ii,2],MS.par[[xx]][ii,3],
                                              "--frag-begin --finp .locfile.txt --N",pop.par[[xx]][ii,1],"--frag-end")), sep="\t", header=T)

      sims <- sims[,grep("^s_mean_",colnames(sims))]
      sims <- sims[,-grep("thomson",colnames(sims))]
      sims <- sims[,-grep("FayWuH",colnames(sims))]

      sim[[xx]] <- sims
    }

    sim <- matrix(unlist(sim), ncol = 6, byrow = TRUE)
    colnames(sim) <- colnames(sims)
  }
  average <- colMeans(sim)
  vari <- apply(sim, 2, var, na.rm = T)
  kur <- apply(sim, 2, kurtosis, na.rm = T)
  skew <- apply(sim, 2, skewness, na.rm = T)

  names(skew)<-gsub("s_mean_","s_skew_",names(skew))
  names(vari)<-gsub("s_mean_","s_var_",names(vari))
  names(kur)<-gsub("s_mean_","s_kurt_",names(kur))

  return(c(average,vari,kur,skew))
}






