#' Shiny GUI for Model Builder
#' @description Launches a Shiny web app for building demographic models.
#'   Uses a shinydashboard layout with sidebar navigation.
#' @param input An optional existing Model object to load as template. Default NULL.
#' @return A Model object when the user clicks "Build Model", or NULL if cancelled.
#' @examples
#' \dontrun{
#' my.model <- main.menu.gui()
#' sim.sumstats(my.model, nsim.blocks = 10, block.size = 10, output.name = "test")
#' }
#' @export
main.menu.gui <- function(input = NULL) {

  # Save the template before Shiny shadows 'input'
  template_model <- input

  # Increase upload limit to 500 MB (default is 5 MB)
  old_limit <- getOption("shiny.maxRequestSize")
  options(shiny.maxRequestSize = 500 * 1024^2)
  on.exit(options(shiny.maxRequestSize = old_limit), add = TRUE)

  # Check for required packages
  for (pkg in c("shinydashboard", "shinyjs", "DT")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(paste0("Package '", pkg, "' is required. Install with: install.packages('", pkg, "')"))
  }

  # ===========================================================================
  # HELPER FUNCTIONS
  # ===========================================================================

  # Validate a Newick tree string
  validate_tree <- function(tree_str) {
    if (tree_str == "1") return(list(valid = TRUE, msg = "Single population model"))
    if (tree_str == "") return(list(valid = FALSE, msg = "Empty tree string"))
    chars <- strsplit(tree_str, "")[[1]]
    n_open  <- sum(chars == "(")
    n_close <- sum(chars == ")")
    n_comma <- sum(chars == ",")
    if (n_open != n_close)
      return(list(valid = FALSE, msg = paste0("Mismatched parentheses: ", n_open, " '(' vs ", n_close, " ')'")))
    if (n_close != n_comma)
      return(list(valid = FALSE, msg = "Missing comma or non-bifurcating tree"))
    return(list(valid = TRUE, msg = "Valid tree"))
  }

  # Parse a Newick tree and extract junction pairs (generalized for multi-char labels)
  parse_tree <- function(tree_str) {
    if (tree_str == "1" || tree_str == "") {
      return(list(npops = 1, joints = NULL, ej = NULL, pop_labels = "1"))
    }

    # Clean: remove branch lengths and semicolons
    clean <- gsub(":[0-9.eE+-]+", "", tree_str)
    clean <- gsub(";", "", clean)
    clean <- trimws(clean)

    # Extract leaf labels
    leaf_str <- gsub("[(),]", " ", clean)
    leaves <- strsplit(trimws(leaf_str), "\\s+")[[1]]
    leaves <- leaves[nchar(leaves) > 0]
    npops <- length(leaves)
    if (npops < 1) return(list(npops = 1, joints = NULL, ej = NULL, pop_labels = "1"))

    # If labels are not simple integers 1:n, map them to numeric
    numeric_check <- suppressWarnings(as.integer(leaves))
    if (any(is.na(numeric_check)) || !setequal(numeric_check, 1:npops)) {
      mapping <- setNames(as.character(1:npops), leaves)
      # Replace longer labels first to avoid partial matches
      for (lab in names(sort(nchar(names(mapping)), decreasing = TRUE))) {
        clean <- gsub(paste0("(?<=[,(])", gsub("([.|()\\^{}+$*?])", "\\\\\\1", lab), "(?=[,)])"),
                       mapping[lab], clean, perl = TRUE)
      }
      pop_labels <- leaves
    } else {
      pop_labels <- leaves
    }

    # Extract junctions by progressively collapsing innermost (X,Y)
    t_str <- clean
    joints <- character()

    while (grepl("(", t_str, fixed = TRUE)) {
      # Match innermost parenthesized pair (no nested parens inside)
      m <- regmatches(t_str, regexpr("\\([^()]+\\)", t_str))
      if (length(m) == 0) break
      inner <- m[1]
      content <- substr(inner, 2, nchar(inner) - 1)
      parts <- trimws(strsplit(content, ",")[[1]])
      if (length(parts) != 2) break

      joints <- c(joints, paste(parts[1], parts[2]))
      t_str <- sub(inner, parts[2], t_str, fixed = TRUE)
    }

    if (length(joints) == 0) {
      return(list(npops = npops, joints = NULL, ej = NULL, pop_labels = pop_labels))
    }

    # Build ej matrix (6 columns, same structure as PipeMaster's .e$ej)
    tot_join_par <- sapply(joints, function(j) {
      paste0("join", paste(strsplit(j, " ")[[1]], collapse = "_"))
    }, USE.NAMES = FALSE)

    ej <- matrix(nrow = length(joints), ncol = 6)
    ej[, 1] <- tot_join_par
    ej[, 2] <- "-ej"
    ej[, 3] <- joints
    ej[, 4] <- "500000"
    ej[, 5] <- "1500000"
    ej[, 6] <- "uniform"

    return(list(npops = npops, joints = joints, ej = ej, pop_labels = pop_labels))
  }

  # Build current Ne parameter matrix
  make_cur_ne <- function(npops, dist = "uniform") {
    nms <- paste0("Ne0.pop", 1:npops)
    n <- matrix(nrow = npops, ncol = 6)
    n[, 1] <- nms
    n[, 2] <- "-n"
    n[, 3] <- as.character(1:npops)
    n[, 4] <- "100000"
    n[, 5] <- "500000"
    n[, 6] <- dist
    n
  }

  # Build ancestral Ne matrices (size + time) for per-pop Ne changes
  make_anc_ne <- function(ne_changes, dist = "uniform") {
    anc_ne_par <- NULL
    time_anc_ne_par <- NULL
    pop_vec <- NULL
    for (i in seq_along(ne_changes)) {
      nc <- ne_changes[i]
      if (!is.na(nc) && nc > 0) {
        for (j in 1:nc) {
          anc_ne_par <- c(anc_ne_par, paste0("Ne", j, ".pop", i))
          time_anc_ne_par <- c(time_anc_ne_par, paste0("t.Ne", j, ".pop", i))
          pop_vec <- c(pop_vec, as.character(i))
        }
      }
    }
    if (is.null(anc_ne_par)) return(NULL)

    en_size <- matrix(nrow = length(anc_ne_par), ncol = 6)
    en_size[, 1] <- anc_ne_par
    en_size[, 2] <- "-en"
    en_size[, 3] <- pop_vec
    en_size[, 4] <- "1000"
    en_size[, 5] <- "10000"
    en_size[, 6] <- dist

    en_time <- matrix(nrow = length(anc_ne_par), ncol = 6)
    en_time[, 1] <- time_anc_ne_par
    en_time[, 2] <- "-en"
    en_time[, 3] <- pop_vec
    en_time[, 4] <- "10000"
    en_time[, 5] <- "100000"
    en_time[, 6] <- dist

    list(size = en_size, time = en_time)
  }

  # Build ancestral Ne matrices for join nodes (Ne.anc_X_Y at each -ej event)
  make_anc_ne_at_joins <- function(ej, dist = "uniform") {
    if (is.null(ej) || !is.matrix(ej) || nrow(ej) == 0) return(NULL)

    n_joins <- nrow(ej)
    en_size <- matrix(nrow = n_joins, ncol = 6)
    en_time <- matrix(nrow = n_joins, ncol = 6)

    for (i in 1:n_joins) {
      pair <- strsplit(ej[i, 3], " ")[[1]]
      surviving_pop <- pair[2]
      par_name <- paste0("Ne.anc_", paste(pair, collapse = "_"))
      time_name <- paste0("t.Ne.anc_", paste(pair, collapse = "_"))

      en_size[i, ] <- c(par_name, "-en", surviving_pop, ej[i, 4], ej[i, 5], dist)
      en_time[i, ] <- c(time_name, "-en", surviving_pop, ej[i, 4], ej[i, 5], dist)
    }

    list(size = en_size, time = en_time)
  }

  # Build migration parameter matrix
  # Build ancestral migration change parameters from current migration
  # mig_changes is a named vector: number of changes per migration parameter
  # e.g., c(1, 2) for 1 change on mig0.1_2 and 2 changes on mig0.2_1
  make_anc_mig <- function(m, mig_changes = NULL, dist = "uniform") {
    if (is.null(m) || !is.matrix(m) || nrow(m) == 0) return(NULL)
    n_mig <- nrow(m)
    if (is.null(mig_changes)) mig_changes <- rep(1, n_mig)
    em_size_rows <- list()
    em_time_rows <- list()
    for (i in 1:n_mig) {
      nc <- mig_changes[i]
      if (is.na(nc) || nc < 1) next
      pair <- m[i, 3]  # e.g., "1 2"
      pair_key <- gsub(" ", "_", pair)
      pop_recv <- strsplit(pair, " ")[[1]][1]
      for (j in 1:nc) {
        em_size_rows[[length(em_size_rows) + 1]] <- c(
          paste0("mig", j, ".", pair_key), "-em", pair, "0.1", "1", dist)
        em_time_rows[[length(em_time_rows) + 1]] <- c(
          paste0("t.mig", j, ".", pair_key), "-em", pop_recv, "10000", "100000", dist)
      }
    }
    if (length(em_size_rows) == 0) return(NULL)
    list(size = do.call(rbind, em_size_rows),
         time = do.call(rbind, em_time_rows))
  }

  make_mig_par <- function(npops, dist = "uniform") {
    mig_par <- NULL
    pops <- NULL
    for (i in 1:npops) {
      for (j in 1:npops) {
        if (i != j) {
          mig_par <- c(mig_par, paste0("mig0.", i, "_", j))
          pops <- c(pops, paste(i, j))
        }
      }
    }
    m <- matrix(nrow = length(mig_par), ncol = 6)
    m[, 1] <- mig_par
    m[, 2] <- "-m"
    m[, 3] <- pops
    m[, 4] <- "0.1"
    m[, 5] <- "1"
    m[, 6] <- dist
    m
  }

  # Build condition matrices (size, mig, time) from current parameters
  build_condition_matrices <- function(n, ej, en, m, em, en_joins = NULL) {
    size_names <- n[, 1]
    time_names <- NULL
    mig_names  <- NULL

    if (!is.null(ej) && is.matrix(ej)) time_names <- ej[, 1]
    if (!is.null(m)  && is.matrix(m))  mig_names  <- m[, 1]
    if (!is.null(en)) {
      size_names <- c(size_names, en$size[, 1])
      time_names <- c(time_names, en$time[, 1])
    }
    if (!is.null(en_joins)) {
      size_names <- c(size_names, en_joins$size[, 1])
      time_names <- c(time_names, en_joins$time[, 1])
    }
    if (!is.null(em)) {
      mig_names  <- c(mig_names, em$size[, 1])
      time_names <- c(time_names, em$time[, 1])
    }

    make_cond_mat <- function(nms) {
      if (is.null(nms) || length(nms) == 0) return(NULL)
      mat <- matrix(NA_character_, nrow = length(nms), ncol = length(nms))
      colnames(mat) <- nms
      rownames(mat) <- nms
      diag(mat) <- "0"
      mat
    }

    list(
      size.matrix = make_cond_mat(size_names),
      mig.matrix  = make_cond_mat(mig_names),
      time.matrix = make_cond_mat(time_names)
    )
  }

  # Mirror helper for condition matrices
  inv_mirror_lower <- function(x) {
    x1 <- t(x)[lower.tri(x, diag = FALSE)]
    for (i in seq_along(x1)) {
      if (is.na(x1[i]) || !x1[i] %in% c("<", ">")) next
      x1[i] <- if (x1[i] == "<") ">" else "<"
    }
    x[lower.tri(x, diag = FALSE)] <- x1
    return(x)
  }

  # Convert a 6-column parameter matrix to a display data.frame
  mat_to_df <- function(mat, dist) {
    if (is.null(mat)) return(NULL)
    df <- data.frame(
      Parameter = mat[, 1],
      Value1 = as.numeric(mat[, 4]),
      Value2 = as.numeric(mat[, 5]),
      stringsAsFactors = FALSE
    )
    if (dist == "uniform") {
      colnames(df) <- c("Parameter", "Min", "Max")
    } else {
      colnames(df) <- c("Parameter", "Mean", "SD")
    }
    df
  }

  # Assemble the Model object (equivalent to get.model)
  assemble_model <- function(rv) {
    convert_dist <- function(x) {
      if (is.null(x)) return(NULL)
      x <- gsub("normal", "rtnorm", x)
      x <- gsub("uniform", "runif", x)
      x
    }

    ej <- convert_dist(rv$ej)
    n  <- convert_dist(rv$n)
    m  <- convert_dist(rv$m)
    en <- NULL
    em <- NULL

    # Merge per-pop Ne changes with ancestral Ne at joins
    en_size <- NULL
    en_time <- NULL
    if (!is.null(rv$en)) {
      en_size <- convert_dist(rv$en$size)
      en_time <- convert_dist(rv$en$time)
    }
    if (!is.null(rv$en_joins)) {
      en_size <- rbind(en_size, convert_dist(rv$en_joins$size))
      en_time <- rbind(en_time, convert_dist(rv$en_joins$time))
    }
    if (!is.null(en_size) && !is.null(nrow(en_size))) {
      en <- list(size = en_size, time = en_time)
    }
    if (!is.null(rv$em)) {
      em <- list(size = convert_dist(rv$em$size), time = convert_dist(rv$em$time))
      if (is.null(nrow(em$size))) em <- NULL
    }
    if (!is.null(m)  && is.null(nrow(m)))  m  <- NULL
    if (!is.null(ej) && is.null(nrow(ej))) ej <- NULL

    loci <- convert_dist(rv$loci)

    model <- list(NULL, NULL, NULL, NULL, NULL)
    names(model) <- c("loci", "I", "flags", "conds", "tree")
    model$loci <- loci
    model$I    <- rv$I

    flags <- list(NULL, NULL, NULL, NULL, NULL)
    names(flags) <- c("n", "m", "en", "em", "ej")
    flags$n  <- n
    flags$m  <- m
    flags$en <- en
    flags$em <- em
    flags$ej <- ej
    model$flags <- flags

    conds <- build_condition_matrices(rv$n, rv$ej, rv$en, rv$m, rv$em, rv$en_joins)
    if (!is.null(rv$size_matrix)) conds$size.matrix <- rv$size_matrix
    if (!is.null(rv$mig_matrix))  conds$mig.matrix  <- rv$mig_matrix
    if (!is.null(rv$time_matrix)) conds$time.matrix  <- rv$time_matrix
    model$conds <- conds
    model$tree  <- rv$tree_string
    model$sum_anc_ne <- rv$sum_anc_ne

    class(model) <- "Model"
    return(model)
  }

  # Revert distribution names from model format to menu format
  revert_dist <- function(x) {
    if (is.null(x)) return(NULL)
    x <- gsub("rtnorm", "normal", x)
    x <- gsub("runif", "uniform", x)
    x
  }

  # Apply condition matrix constraints to a named numeric vector
  # "<" : row param < col param  (swap if violated)
  # ">" : row param > col param  (swap if violated)
  # "=" : row param = col param  (copy col to row)
  apply_cond <- function(vals, cond_mat) {
    if (is.null(cond_mat) || is.null(vals) || length(vals) == 0) return(vals)
    nms <- intersect(names(vals), rownames(cond_mat))
    if (length(nms) == 0) return(vals)
    sub <- cond_mat[nms, nms, drop = FALSE]

    # "<" conditions: ensure row < col
    lt_idx <- which(sub == "<", arr.ind = TRUE)
    if (nrow(lt_idx) > 0) {
      for (k in 1:nrow(lt_idx)) {
        n1 <- nms[lt_idx[k, 1]]; n2 <- nms[lt_idx[k, 2]]
        if (vals[n1] > vals[n2]) {
          tmp <- vals[n1]; vals[n1] <- vals[n2]; vals[n2] <- tmp
        }
      }
    }
    # ">" conditions: ensure row > col
    gt_idx <- which(sub == ">", arr.ind = TRUE)
    if (nrow(gt_idx) > 0) {
      for (k in 1:nrow(gt_idx)) {
        n1 <- nms[gt_idx[k, 1]]; n2 <- nms[gt_idx[k, 2]]
        if (vals[n1] < vals[n2]) {
          tmp <- vals[n1]; vals[n1] <- vals[n2]; vals[n2] <- tmp
        }
      }
    }
    # "=" conditions: copy col value to row
    eq_idx <- which(sub == "=", arr.ind = TRUE)
    if (nrow(eq_idx) > 0) {
      for (k in 1:nrow(eq_idx)) {
        n1 <- nms[eq_idx[k, 1]]; n2 <- nms[eq_idx[k, 2]]
        vals[n1] <- vals[n2]
      }
    }
    vals
  }

  # Sample once from a parameter matrix, returning named numeric vector
  sample_par_vec <- function(mat) {
    if (is.null(mat) || !is.matrix(mat) || nrow(mat) == 0) return(NULL)
    vals <- numeric(nrow(mat))
    for (i in 1:nrow(mat)) {
      v1 <- as.numeric(mat[i, 4]); v2 <- as.numeric(mat[i, 5])
      if (grepl("uniform", mat[i, 6], ignore.case = TRUE)) {
        vals[i] <- stats::runif(1, v1, v2)
      } else {
        vals[i] <- msm::rtnorm(1, v1, v2, lower = 0)
      }
    }
    setNames(vals, mat[, 1])
  }

  # Sample all model parameters once (from priors), grouped and conditioned
  sample_all_conditioned <- function(rv) {
    size_v <- c(sample_par_vec(rv$n),
                if (!is.null(rv$en)) sample_par_vec(rv$en$size),
                if (!is.null(rv$en_joins)) sample_par_vec(rv$en_joins$size))
    time_v <- c(if (!is.null(rv$ej) && is.matrix(rv$ej)) sample_par_vec(rv$ej),
                if (!is.null(rv$en)) sample_par_vec(rv$en$time),
                if (!is.null(rv$en_joins)) sample_par_vec(rv$en_joins$time),
                if (!is.null(rv$em)) sample_par_vec(rv$em$time))
    mig_v  <- c(if (!is.null(rv$m) && is.matrix(rv$m)) sample_par_vec(rv$m),
                if (!is.null(rv$em)) sample_par_vec(rv$em$size))
    size_v <- apply_cond(size_v, rv$size_matrix)
    time_v <- apply_cond(time_v, rv$time_matrix)
    mig_v  <- apply_cond(mig_v,  rv$mig_matrix)
    c(size_v, time_v, mig_v)
  }

  # Get conditioned parameter means by sampling many times and averaging
  # This produces correct conditional means (e.g. E[A | A < B] != E[A])
  conditioned_means <- function(rv, n_samples = 500) {
    all_samp <- replicate(n_samples, sample_all_conditioned(rv))
    if (is.null(all_samp)) return(NULL)
    if (is.matrix(all_samp)) return(rowMeans(all_samp))
    # single parameter case
    setNames(mean(all_samp), names(sample_all_conditioned(rv)))
  }

  # Custom demographic model plot (tree-like with rectangles for populations)
  plot_model_tree <- function(rv, use_avg = TRUE, font_scale = 1.0, n_samples = 200,
                              palette = "Tableau", bg_theme = "Dark", spacing = 1.5,
                              font_family = "Palatino", plot_title = NULL) {
    npops <- rv$npops
    if (is.null(npops) || npops < 1) {
      plot.new(); text(0.5, 0.5, "No populations defined", col = "white", cex = 1.5)
      return(invisible(NULL))
    }

    # Use msABC.commander() for parameter sampling — same code path as simulations
    model <- tryCatch(assemble_model(rv), error = function(e) NULL)
    if (is.null(model)) {
      plot.new(); text(0.5, 0.5, "Model not ready", col = "white", cex = 1.5)
      return(invisible(NULL))
    }
    # Ensure loci exist (dummy if not configured yet)
    if (is.null(model$loci) || !is.matrix(model$loci) || nrow(model$loci) == 0) {
      model$loci <- matrix(c("loc1", "100", "1", "1e-8", "1e-8", "runif"), nrow = 1)
    }
    # Ensure I exists
    if (is.null(model$I) || !is.matrix(model$I)) {
      samp_per_pop <- rep("10", npops)
      model$I <- matrix(c("I", paste(npops), as.character(npops), samp_per_pop),
                         nrow = 1)
    }

    # Sample parameters via msABC.commander (handles conditions, Ne.anc, sum)
    sample_once <- function() {
      cmd <- tryCatch(
        PipeMaster:::msABC.commander(model, use.alpha = FALSE, arg = "/tmp/plot/"),
        error = function(e) NULL)
      if (is.null(cmd)) return(NULL)
      p <- cmd[[2]]
      vals <- as.numeric(p[2, ])
      names(vals) <- p[1, ]
      vals
    }

    if (use_avg) {
      all_samp <- replicate(n_samples, sample_once())
      if (is.null(all_samp) || !is.matrix(all_samp)) {
        plot.new(); text(0.5, 0.5, "Could not sample parameters", col = "white", cex = 1.5)
        return(invisible(NULL))
      }
      cond_vals <- rowMeans(all_samp, na.rm = TRUE)
    } else {
      cond_vals <- sample_once()
      if (is.null(cond_vals)) {
        plot.new(); text(0.5, 0.5, "Could not sample parameters", col = "white", cex = 1.5)
        return(invisible(NULL))
      }
    }
    lookup <- function(name) unname(cond_vals[name])

    # Current Ne per population
    ne_current <- rep(100000, npops)
    if (!is.null(rv$n) && is.matrix(rv$n)) {
      for (i in 1:nrow(rv$n)) {
        pop_id <- as.integer(rv$n[i, 3])
        if (!is.na(pop_id) && pop_id >= 1 && pop_id <= npops)
          ne_current[pop_id] <- lookup(rv$n[i, 1])
      }
    }

    # Collect time-ordered events (joins + Ne changes)
    # Ne.anc events are NOT added here — they're handled directly at merge time
    events <- list()
    if (!is.null(rv$ej) && is.matrix(rv$ej) && nrow(rv$ej) > 0) {
      for (i in 1:nrow(rv$ej)) {
        pair <- strsplit(rv$ej[i, 3], " ")[[1]]
        pair_key <- paste(pair, collapse = "_")
        ne_anc_name <- paste0("Ne.anc_", pair_key)
        ne_anc_val <- lookup(ne_anc_name)
        events[[length(events) + 1]] <- list(
          type = "join", time = lookup(rv$ej[i, 1]),
          src = as.integer(pair[1]), tgt = as.integer(pair[2]),
          ne_anc = if (!is.na(ne_anc_val)) ne_anc_val else NA)
      }
    }
    if (!is.null(rv$en) && !is.null(rv$en$size) && is.matrix(rv$en$size)) {
      for (i in 1:nrow(rv$en$size)) {
        # Skip Ne.anc entries (handled at merge time)
        if (grepl("^Ne\\.anc_", rv$en$size[i, 1])) next
        events[[length(events) + 1]] <- list(
          type = "ne_change", time = lookup(rv$en$time[i, 1]),
          pop = as.integer(rv$en$size[i, 3]),
          ne = lookup(rv$en$size[i, 1]))
      }
    }
    # Sort: ne_changes before joins at equal times
    if (length(events) > 0)
      events <- events[order(
        sapply(events, `[[`, "time"),
        sapply(events, function(e) if (e$type == "join") 1 else 0)
      )]

    # Track state and build drawing primitives
    # Position populations by tree leaf order (keeps merging pops adjacent)
    leaf_order <- 1:npops
    if (!is.null(rv$tree_string) && rv$tree_string != "1" &&
        rv$tree_string != "non tree-like model") {
      parsed_tree <- parse_tree(rv$tree_string)
      if (!is.null(parsed_tree$pop_labels)) {
        lo <- suppressWarnings(as.integer(parsed_tree$pop_labels))
        if (!any(is.na(lo)) && length(lo) == npops) leaf_order <- lo
      }
    }
    # Spacing must prevent overlap between daughter bars and merged ancestral bars
    # at the midpoint: spacing/2 > 2*hw_max + gap, so spacing > 4*hw_max + 2*gap
    # hw_max = 0.38 is defined below but is a fixed constant
    leaf_spacing <- spacing
    x_pos  <- setNames(1 + (seq_along(leaf_order) - 1) * leaf_spacing, as.character(leaf_order))
    x_pos_orig <- x_pos  # Save original positions for bottom labels
    alive  <- setNames(rep(TRUE, npops), as.character(1:npops))
    cur_ne <- setNames(ne_current, as.character(1:npops))
    last_t <- setNames(rep(0, npops), as.character(1:npops))
    segs   <- list()
    merges <- list()

    for (ev in events) {
      if (ev$type == "ne_change") {
        p <- as.character(ev$pop)
        if (alive[p]) {
          segs[[length(segs) + 1]] <- list(
            pop = p, x = x_pos[p], t0 = last_t[p], t1 = ev$time, ne = cur_ne[p])
          cur_ne[p] <- ev$ne
          last_t[p] <- ev$time
        }
      } else if (ev$type == "join") {
        s <- as.character(ev$src)
        g <- as.character(ev$tgt)
        if (alive[s] && alive[g]) {
          # Close daughter segments at join time
          segs[[length(segs) + 1]] <- list(
            pop = s, x = x_pos[s], t0 = last_t[s], t1 = ev$time, ne = cur_ne[s])
          segs[[length(segs) + 1]] <- list(
            pop = g, x = x_pos[g], t0 = last_t[g], t1 = ev$time, ne = cur_ne[g])
          new_x <- mean(c(x_pos[s], x_pos[g]))
          # Ancestor Ne: use Ne.anc from msABC.commander if available, else sum daughters
          ne_anc <- if (!is.na(ev$ne_anc)) ev$ne_anc else cur_ne[s] + cur_ne[g]
          merges[[length(merges) + 1]] <- list(
            time = ev$time,
            x_src = x_pos[s], ne_src = cur_ne[s],
            x_tgt = x_pos[g], ne_tgt = cur_ne[g],
            x_new = new_x, ne_new = ne_anc)
          alive[s] <- FALSE
          x_pos[g] <- new_x
          cur_ne[g] <- ne_anc
          last_t[g] <- ev$time
        }
      }
    }

    # Close remaining active lineages (root extension)
    all_t1  <- if (length(segs) > 0) sapply(segs, `[[`, "t1") else 0
    max_t   <- max(c(all_t1, 1))
    root_ext <- max_t * 0.25
    # Island populations (never involved in any merge) stop at max_t;
    # only tree lineages (participated in merges) get root extension
    merged_pops <- character()
    for (ev in events) {
      if (ev$type == "join") merged_pops <- c(merged_pops, as.character(ev$src), as.character(ev$tgt))
    }
    has_joins <- length(merged_pops) > 0
    for (p in names(alive)) {
      if (alive[p]) {
        # Island pops stop at 80% of max_t, but only when there are tree lineages;
        # if no joins at all (single pop, pure island), everyone gets root extension
        is_tree_lineage <- p %in% merged_pops
        if (!has_joins || is_tree_lineage) {
          t_end <- max_t + root_ext
        } else {
          t_end <- max_t * 0.8
        }
        segs[[length(segs) + 1]] <- list(
          pop = p, x = x_pos[p], t0 = last_t[p],
          t1 = t_end, ne = cur_ne[p])
      }
    }

    # Scaling
    all_ne <- sapply(segs, `[[`, "ne")
    y_max  <- max(sapply(segs, `[[`, "t1"))
    ne_max <- max(all_ne); ne_min <- min(all_ne)
    hw_max <- 0.55; hw_min <- 0.18
    ne2hw <- function(ne) {
      if (ne_max == ne_min) return((hw_max + hw_min) / 2)
      hw_min + (hw_max - hw_min) * (ne - ne_min) / (ne_max - ne_min)
    }

    palettes <- list(
      "Tableau"    = c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2",
                       "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7",
                       "#9C755F", "#BAB0AC"),
      "Colorblind" = c("#0072B2", "#E69F00", "#009E73", "#CC79A7",
                       "#56B4E9", "#D55E00", "#F0E442", "#999999",
                       "#000000", "#FFFFFF"),
      "Pastel"     = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
                       "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
                       "#D9D9D9", "#BC80BD"),
      "Vivid"      = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                       "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
                       "#66C2A5", "#999999")
    )
    pal <- palettes[[if (!is.null(palette)) palette else "Tableau"]]
    pcol <- function(p) pal[((as.integer(p) - 1) %% length(pal)) + 1]

    # Background themes
    bg_themes <- list(
      "Dark"       = list(bg = "#222222", fg = "white", axis = "#cccccc",
                          grid = "#666666", annot = "#aaaaaa", merge = "#888888",
                          ne_text = "white", empty = "gray50"),
      "Light"      = list(bg = "#FFFFFF", fg = "black", axis = "#333333",
                          grid = "#999999", annot = "#555555", merge = "#aaaaaa",
                          ne_text = "black", empty = "gray60"),
      "Slate"      = list(bg = "#2F4F4F", fg = "white", axis = "#cccccc",
                          grid = "#5F7F7F", annot = "#bbbbbb", merge = "#7F9F9F",
                          ne_text = "white", empty = "gray60")
    )
    th <- bg_themes[[if (!is.null(bg_theme)) bg_theme else "Dark"]]

    delta <- y_max * 0.015

    # Plot area
    all_x <- unlist(lapply(segs, function(s) c(s$x - ne2hw(s$ne), s$x + ne2hw(s$ne))))
    x_lo <- min(all_x) - 0.3; x_hi <- max(all_x) + 0.3

    par(mar = c(3.5, 7.5, 2.5, 1.5), bg = th$bg, fg = th$fg,
        col.axis = th$axis, col.lab = th$axis, col.main = th$fg,
        family = font_family)
    plot(NULL, xlim = c(x_lo, x_hi), ylim = c(-y_max * 0.07, y_max),
         xlab = "", ylab = "",
         main = if (!is.null(plot_title) && nchar(plot_title) > 0) plot_title else "Demographic Model",
         axes = FALSE, cex.lab = 1.3 * font_scale, cex.main = 1.5 * font_scale)
    ax_ticks <- axTicks(2)
    axis(2, at = ax_ticks, labels = ax_ticks / 1000, las = 1,
         col = th$grid, col.ticks = th$grid, col.axis = th$axis, cex.axis = 1.1 * font_scale)
    mtext("Time (x1000 generations ago)", side = 2, line = 5.5,
          cex = 1.3 * font_scale, col = th$axis)

    # Draw merge connectors first (behind rectangles)
    for (mg in merges) {
      hw_s <- ne2hw(mg$ne_src); hw_t <- ne2hw(mg$ne_tgt); hw_n <- ne2hw(mg$ne_new)
      if (mg$x_src < mg$x_tgt) {
        lx <- mg$x_src; lhw <- hw_s; rx <- mg$x_tgt; rhw <- hw_t
      } else {
        lx <- mg$x_tgt; lhw <- hw_t; rx <- mg$x_src; rhw <- hw_s
      }
      nx <- mg$x_new; nhw <- hw_n
      # Left trapezoid
      polygon(x = c(lx - lhw, lx + lhw, nx, nx - nhw),
              y = c(mg$time, mg$time, mg$time + delta, mg$time + delta),
              col = adjustcolor(th$merge, 0.5), border = NA)
      # Right trapezoid
      polygon(x = c(rx - rhw, rx + rhw, nx + nhw, nx),
              y = c(mg$time, mg$time, mg$time + delta, mg$time + delta),
              col = adjustcolor(th$merge, 0.5), border = NA)
    }

    # Draw rectangles (widest first so daughter bars are drawn on top of ancestors)
    seg_order <- order(sapply(segs, function(s) ne2hw(s$ne)), decreasing = TRUE)
    for (idx in seg_order) {
      seg <- segs[[idx]]
      hw <- ne2hw(seg$ne)
      col <- pcol(seg$pop)
      rect(seg$x - hw, seg$t0, seg$x + hw, seg$t1,
           col = adjustcolor(col, alpha.f = 0.75),
           border = adjustcolor(col, alpha.f = 0.9), lwd = 1.5)
    }

    # Population labels (positioned by tree leaf order)
    labels <- if (!is.null(rv$pop_labels)) rv$pop_labels else as.character(1:npops)
    for (i in 1:npops) {
      text(x_pos_orig[as.character(i)], -y_max * 0.04, paste0("Pop ", labels[i]),
           col = pcol(as.character(i)), cex = 1.4 * font_scale, font = 2)
    }

    # Ne annotations inside rectangles
    for (seg in segs) {
      seg_h <- seg$t1 - seg$t0
      if (seg_h > y_max * 0.05) {
        text(seg$x, (seg$t0 + seg$t1) / 2,
             paste0(round(seg$ne / 1000), "k"),
             col = th$ne_text, cex = 1.0 * font_scale, srt = 90)
      }
    }

    # Time annotations at merge points (positioned next to the node)
    for (mg in merges) {
      hw_n <- ne2hw(mg$ne_new)
      text(mg$x_new + hw_n + 0.15, mg$time,
           formatC(mg$time, format = "e", digits = 2),
           col = th$annot, cex = 0.9 * font_scale, adj = c(0, 0.5))
    }

    # Migration arrows between populations
    if (!is.null(rv$m) && is.matrix(rv$m) && nrow(rv$m) > 0) {
      # Find the bottom-most segments for each leaf pop (the ones starting at t=0)
      leaf_segs <- list()
      for (seg in segs) {
        if (seg$t0 == 0) leaf_segs[[seg$pop]] <- seg
      }

      for (i in 1:nrow(rv$m)) {
        pair <- as.integer(strsplit(rv$m[i, 3], " ")[[1]])
        from_p <- as.character(pair[1])
        to_p   <- as.character(pair[2])
        if (is.null(leaf_segs[[from_p]]) || is.null(leaf_segs[[to_p]])) next

        seg_from <- leaf_segs[[from_p]]
        seg_to   <- leaf_segs[[to_p]]
        hw_from  <- ne2hw(seg_from$ne)
        hw_to    <- ne2hw(seg_to$ne)

        # Place arrow at 1/3 height of the shorter segment
        arr_y <- min(seg_from$t1, seg_to$t1) * 0.33

        # Arrow goes from the edge of one rectangle to the edge of the other
        if (seg_from$x < seg_to$x) {
          x0 <- seg_from$x + hw_from
          x1 <- seg_to$x   - hw_to
        } else {
          x0 <- seg_from$x - hw_from
          x1 <- seg_to$x   + hw_to
        }

        # Offset vertically so bidirectional arrows don't overlap
        y_off <- if (pair[1] < pair[2]) y_max * 0.012 else -y_max * 0.012

        arrows(x0, arr_y + y_off, x1, arr_y + y_off,
               col = "#FFD700", lwd = 1.8, length = 0.08, code = 2)

        # Migration rate label at midpoint
        mig_val <- lookup(rv$m[i, 1])
        text((x0 + x1) / 2, arr_y + y_off + y_max * 0.018,
             sprintf("%.2g", mig_val),
             col = "#FFD700", cex = 0.8 * font_scale)
      }
    }
  }

  # ===========================================================================
  # UI
  # ===========================================================================
  ui <- shinydashboard::dashboardPage(
    skin = "blue",

    # --- Header ---
    shinydashboard::dashboardHeader(
      title = "PipeMaster - Model Builder",
      titleWidth = 300
    ),

    # --- Sidebar ---
    shinydashboard::dashboardSidebar(
      width = 280,
      shiny::uiOutput("sidebar_logo"),
      shinydashboard::sidebarMenu(
        id = "sidebar",
        shinydashboard::menuItem("Getting Started",      tabName = "start",      icon = shiny::icon("play-circle")),
        shinydashboard::menuItem("Population Structure",  tabName = "structure",  icon = shiny::icon("sitemap")),
        shinydashboard::menuItem("Demography (Ne)",       tabName = "demography", icon = shiny::icon("chart-line")),
        shinydashboard::menuItem("Migration",             tabName = "migration",  icon = shiny::icon("exchange-alt")),
        shinydashboard::menuItem("Time Priors",           tabName = "time",       icon = shiny::icon("clock")),
        shinydashboard::menuItem("Conditions",            tabName = "conditions", icon = shiny::icon("filter")),
        shinydashboard::menuItem("Gene Setup",            tabName = "genes",      icon = shiny::icon("dna")),
        shinydashboard::menuItem("Visualize & Export",        tabName = "export",     icon = shiny::icon("download"))
      ),
      shiny::hr(),
      shinydashboard::box(
        title = "Model Status", width = 12, solidHeader = FALSE,
        shinydashboard::valueBoxOutput("status_pops",  width = 12),
        shinydashboard::valueBoxOutput("status_nodes", width = 12),
        shinydashboard::valueBoxOutput("status_loci",  width = 12)
      )
    ),

    # --- Body ---
    shinydashboard::dashboardBody(
      shiny::tags$head(shiny::tags$style(shiny::HTML("
        .content-wrapper, .right-side { background-color: #ecf0f5; }
        .box { border-radius: 5px; }
        .dataTable { font-size: 13px; }
        .error-message { color: #dd4b39; font-weight: bold; padding: 10px;
                         background-color: #f2dede; border-radius: 5px; }
        .success-message { color: #00a65a; font-weight: bold; padding: 10px;
                           background-color: #d4edda; border-radius: 5px; }
      "))),
      shinyjs::useShinyjs(),

      shinydashboard::tabItems(

        # =====================================================================
        # TAB: Getting Started
        # =====================================================================
        shinydashboard::tabItem(tabName = "start",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Welcome to PipeMaster GUI", width = 12,
              status = "primary", solidHeader = TRUE,
              shiny::h3("Build your coalescent model visually"),
              shiny::p("This graphical interface replaces the text-based main.menu()."),
              shiny::hr(),
              shiny::h4("Workflow:"),
              shiny::tags$ol(
                shiny::tags$li(shiny::strong("Population Structure:"), " Define your population tree topology"),
                shiny::tags$li(shiny::strong("Demography:"), " Set effective population size (Ne) priors"),
                shiny::tags$li(shiny::strong("Migration:"), " Configure gene flow between populations (optional)"),
                shiny::tags$li(shiny::strong("Time Priors:"), " Set temporal parameters for divergences and Ne changes"),
                shiny::tags$li(shiny::strong("Conditions:"), " Add parameter constraints (optional)"),
                shiny::tags$li(shiny::strong("Gene Setup:"), " Configure loci manually or load from data (PHYLIP/VCF)"),
                shiny::tags$li(shiny::strong("Visualize & Export:"), " Generate model object for sim.sumstats()")
              ),
              shiny::hr(),
              shiny::h4("Quick Start"),
              shiny::tags$pre(
                "model <- main.menu.gui()\n",
                "sim.sumstats(model, nsim.blocks = 100, block.size = 10,\n",
                "             output.name = \"my_sims\")"
              ),
              shiny::hr(),
              shiny::h4("Tree examples:"),
              shiny::tags$ul(
                shiny::tags$li(shiny::strong("Single population:"), shiny::code("1")),
                shiny::tags$li(shiny::strong("Two populations:"), shiny::code("(1,2)")),
                shiny::tags$li(shiny::strong("Three populations:"), shiny::code("(1,(2,3))")),
                shiny::tags$li(shiny::strong("Four populations:"), shiny::code("((1,2),(3,4))")),
                shiny::tags$li(shiny::strong("Named tips:"), shiny::code("((A,B),(C,D))"))
              )
            )
          )
        ),

        # =====================================================================
        # TAB: Population Structure
        # =====================================================================
        shinydashboard::tabItem(tabName = "structure",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Population Structure", width = 12,
              status = "primary", solidHeader = TRUE,
              shiny::p("Define the topology of your population tree in Newick format, or enter '1' for a single population."),
              shiny::fluidRow(
                shiny::column(8,
                  shiny::textInput("txt_tree",
                    "Tree Topology (Newick format):",
                    value = "1",
                    placeholder = "e.g., (1,(2,3)) or ((A,B),(C,D))")
                ),
                shiny::column(4,
                  shiny::br(),
                  shiny::actionButton("btn_validate_tree", "Validate & Apply Tree",
                    icon = shiny::icon("check-circle"), class = "btn-success")
                )
              ),
              shiny::uiOutput("tree_validation_message"),
              shiny::hr(),
              shiny::fluidRow(
                shiny::column(6, shiny::h4("Model Visualization:")),
                shiny::column(3, shiny::helpText(
                  shiny::icon("info-circle"),
                  "More options in the Visualize & Export tab."
                )),
                shiny::column(3, shiny::downloadButton("btn_save_model_pdf", "Save as PDF",
                  class = "btn-sm btn-default", style = "margin-top: 5px;"))
              ),
              shiny::plotOutput("plot_model_structure", height = "450px")
            )
          ),
          shiny::fluidRow(
            shinydashboard::box(
              title = "Detected Junctions (Divergence Events)", width = 6,
              status = "info", solidHeader = TRUE,
              DT::DTOutput("table_junctions"),
              shiny::hr(),
              shiny::h4("Island Model"),
              shiny::p("Remove divergence nodes to simulate an island model (populations without shared ancestry, connected only by migration)."),
              shiny::uiOutput("island_node_checkboxes"),
              shiny::actionButton("btn_remove_nodes", "Remove Selected Nodes",
                icon = shiny::icon("scissors"), class = "btn-danger",
                style = "margin-top: 5px;")
            ),
            shinydashboard::box(
              title = "Quick Examples", width = 6,
              status = "warning", solidHeader = TRUE,
              shiny::actionButton("btn_ex_2pop",  "(1,2)",          class = "btn-default", style = "margin: 2px;"),
              shiny::actionButton("btn_ex_3pop",  "(1,(2,3))",      class = "btn-default", style = "margin: 2px;"),
              shiny::actionButton("btn_ex_4pop",  "((1,2),(3,4))",  class = "btn-default", style = "margin: 2px;"),
              shiny::actionButton("btn_ex_1pop",  "1 (single pop)", class = "btn-default", style = "margin: 2px;")
            )
          )
        ),

        # =====================================================================
        # TAB: Demography (Ne)
        # =====================================================================
        shinydashboard::tabItem(tabName = "demography",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Ne Prior Distribution", width = 12,
              status = "primary", solidHeader = TRUE,
              shiny::fluidRow(
                shiny::column(6,
                  shiny::selectInput("select_ne_dist", "Prior Distribution:",
                    choices = c("Uniform" = "uniform", "Normal" = "normal"),
                    selected = "uniform")
                ),
                shiny::column(6,
                  shiny::checkboxInput("check_ancestral_ne",
                    "Include Ne changes through time", value = FALSE),
                  shiny::checkboxInput("check_sum_anc_ne",
                    "Ancestral Ne = sum of daughter Ne at join events", value = TRUE)
                )
              )
            )
          ),
          shiny::fluidRow(
            shinydashboard::box(
              title = "Current Ne Priors", width = 12,
              status = "info", solidHeader = TRUE,
              DT::DTOutput("table_current_ne"),
              shiny::br(),
              shiny::helpText("Click a cell to edit. The Parameter column is read-only.")
            )
          ),
          shiny::conditionalPanel(
            condition = "input.check_sum_anc_ne == false",
            shiny::fluidRow(
              shinydashboard::box(
                title = "Ancestral Ne at Join Events", width = 12,
                status = "success", solidHeader = TRUE,
                shiny::p("Each divergence node has an independent ancestral Ne parameter."),
                DT::DTOutput("table_anc_ne_joins"),
                shiny::br(),
                shiny::helpText("Click a cell to edit. When 'Sum daughter Ne' is checked, these values are computed automatically.")
              )
            )
          ),
          shiny::conditionalPanel(
            condition = "input.check_ancestral_ne == true",
            shiny::fluidRow(
              shinydashboard::box(
                title = "Ne Changes Per Population", width = 12,
                status = "warning", solidHeader = TRUE,
                shiny::uiOutput("ne_changes_per_pop_ui"),
                shiny::actionButton("btn_apply_ne_changes", "Apply Ne Changes",
                  icon = shiny::icon("sync"), class = "btn-primary")
              )
            ),
            shiny::fluidRow(
              shinydashboard::box(
                title = "Ancestral Ne Priors", width = 12,
                status = "info", solidHeader = TRUE,
                DT::DTOutput("table_ancestral_ne"),
                shiny::br(),
                shiny::helpText("Click a cell to edit.")
              )
            )
          )
        ),

        # =====================================================================
        # TAB: Migration
        # =====================================================================
        shinydashboard::tabItem(tabName = "migration",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Gene Flow Configuration", width = 12,
              status = "primary", solidHeader = TRUE,
              shiny::fluidRow(
                shiny::column(4,
                  shiny::checkboxInput("check_migration",
                    "Enable migration between populations", value = FALSE)
                ),
                shiny::column(4,
                  shiny::conditionalPanel(
                    condition = "input.check_migration == true",
                    shiny::checkboxInput("check_anc_mig",
                      "Enable ancestral migration changes", value = FALSE)
                  )
                ),
                shiny::column(4,
                  shiny::conditionalPanel(
                    condition = "input.check_migration == true",
                    shiny::selectInput("select_mig_dist", "Migration Prior Distribution:",
                      choices = c("Uniform" = "uniform", "Normal" = "normal"),
                      selected = "uniform")
                  )
                )
              )
            )
          ),
          shiny::conditionalPanel(
            condition = "input.check_migration == true",
            shiny::fluidRow(
              shinydashboard::box(
                title = "Current Migration Rates (4Nm)", width = 12,
                status = "info", solidHeader = TRUE,
                shiny::helpText(
                  "Migration is specified as 4Nm (population-scaled rate).",
                  shiny::tags$b("mig0.i_j"), "= pop i receives migrants from pop j.",
                  "For example,", shiny::tags$b("mig0.1_2"), "= gene flow into pop 1 from pop 2."
                ),
                shiny::tags$div(
                  class = "alert alert-warning", style = "margin-top: 10px;",
                  shiny::icon("exclamation-triangle"),
                  shiny::tags$b(" Ne changes affect migration rates."),
                  " Migration rates (M = 4Nm) are scaled by the current Ne (Ne0) of the",
                  " receiving population. Internally, the coalescent simulator stores the",
                  " per-generation rate m = M/Ne0, which remains fixed even if Ne changes",
                  " over time via ancestral Ne. This means the effective number of migrants",
                  " per generation (Nm) changes proportionally to the Ne change of the",
                  " receiving population.",
                  shiny::tags$br(), shiny::tags$br(),
                  "For example, if Ne doubles at some time in the past, the effective Nm",
                  " also doubles from that point backward, even though the user-specified M",
                  " was set relative to the current Ne0.",
                  shiny::tags$br(), shiny::tags$br(),
                  shiny::tags$b("To maintain constant Nm through time:"),
                  " add an ancestral migration change with the same M value at the same",
                  " time as the Ne change. PipeMaster will automatically rescale the",
                  " migration rate by the Ne at that time point, keeping Nm = M/4 constant.",
                  " In the Conditions tab, set:",
                  shiny::tags$ul(
                    shiny::tags$li(shiny::tags$code("t.mig1.i_j = t.Ne1.popi"),
                      " (sync migration change time to Ne change time)"),
                    shiny::tags$li(shiny::tags$code("mig1.i_j = mig0.i_j"),
                      " (same M value, PipeMaster rescales by ancestral Ne)")
                  ),
                  "This works for both msABC and scrm backends."
                ),
                DT::DTOutput("table_migration"),
                shiny::br(),
                shiny::helpText("Click a cell to edit."),
                shiny::hr(),
                shiny::h5("Remove Migration Parameters:"),
                shiny::uiOutput("mig_remove_checkboxes"),
                shiny::actionButton("btn_remove_mig", "Remove Selected",
                  icon = shiny::icon("trash"), class = "btn-danger btn-sm")
              )
            ),
            shiny::conditionalPanel(
              condition = "input.check_migration == true && input.check_anc_mig == true",
              shiny::fluidRow(
                shinydashboard::box(
                  title = "Migration Changes Per Parameter", width = 12,
                  status = "warning", solidHeader = TRUE,
                  shiny::uiOutput("mig_changes_per_pair_ui"),
                  shiny::actionButton("btn_apply_mig_changes", "Apply Migration Changes",
                    icon = shiny::icon("sync"), class = "btn-primary")
                )
              ),
              shiny::fluidRow(
                shinydashboard::box(
                  title = "Ancestral Migration Rates (4Nm)", width = 12,
                  status = "info", solidHeader = TRUE,
                  shiny::helpText(
                    "Ancestral migration changes (-em) allow migration rates to change",
                    " at a specific time in the past. To maintain constant Nm when Ne",
                    " changes, set", shiny::tags$code("t.mig1.i_j = t.Ne1.popi"),
                    "and", shiny::tags$code("mig1.i_j = mig0.i_j"),
                    "in the Conditions tab."
                  ),
                  DT::DTOutput("table_anc_mig_rates"),
                  shiny::br(),
                  shiny::helpText("Click a cell to edit.")
                )
              )
            )
          )
        ),

        # =====================================================================
        # TAB: Time Priors
        # =====================================================================
        shinydashboard::tabItem(tabName = "time",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Temporal Parameters", width = 12,
              status = "primary", solidHeader = TRUE,
              shiny::helpText(
                shiny::icon("info-circle"),
                " All times must be specified in", shiny::tags$b("generations"),
                "(not years). To convert from years, divide by generation time."
              ),
              shiny::selectInput("select_time_dist", "Time Prior Distribution:",
                choices = c("Uniform" = "uniform", "Normal" = "normal"),
                selected = "uniform"),
              shiny::tabsetPanel(
                id = "time_tabs",
                shiny::tabPanel("Divergence Times",
                  shiny::br(),
                  shiny::uiOutput("div_times_ui"),
                  DT::DTOutput("table_div_times"),
                  shiny::br(),
                  shiny::helpText("Click a cell to edit. Available when tree has junctions.")
                ),
                shiny::tabPanel("Ne Change Times",
                  shiny::br(),
                  shiny::conditionalPanel(
                    condition = "input.check_ancestral_ne == true",
                    DT::DTOutput("table_ne_change_times"),
                    shiny::br(),
                    shiny::helpText("Click a cell to edit.")
                  ),
                  shiny::conditionalPanel(
                    condition = "input.check_ancestral_ne == false",
                    shiny::p("Enable Ne changes through time in the Demography tab first.")
                  )
                ),
                shiny::tabPanel("Migration Change Times",
                  shiny::br(),
                  shiny::conditionalPanel(
                    condition = "input.check_migration == true && input.check_anc_mig == true",
                    DT::DTOutput("table_anc_mig_times"),
                    shiny::br(),
                    shiny::helpText("Click a cell to edit.")
                  ),
                  shiny::conditionalPanel(
                    condition = "input.check_migration == false",
                    shiny::p("Enable migration in the Migration tab first.")
                  ),
                  shiny::conditionalPanel(
                    condition = "input.check_migration == true && input.check_anc_mig == false",
                    shiny::p("Enable ancestral migration changes in the Migration tab first.")
                  )
                )
              )
            )
          )
        ),

        # =====================================================================
        # TAB: Conditions
        # =====================================================================
        shinydashboard::tabItem(tabName = "conditions",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Parameter Constraints", width = 12,
              status = "primary", solidHeader = TRUE,
              shiny::p("Add constraints between parameters (e.g., Ne0.pop1 > Ne0.pop2)."),
              shiny::tabsetPanel(
                id = "cond_tabs",
                shiny::tabPanel("Size Conditions",
                  shiny::br(),
                  shiny::fluidRow(
                    shiny::column(4, shiny::uiOutput("cond_size_par1_ui")),
                    shiny::column(2, shiny::selectInput("cond_size_op", "Operator:",
                      choices = c("<", ">", "="), selected = "<")),
                    shiny::column(4, shiny::uiOutput("cond_size_par2_ui")),
                    shiny::column(2, shiny::br(),
                      shiny::actionButton("btn_add_size_cond", "Add", class = "btn-success"))
                  ),
                  shiny::hr(),
                  shiny::h5("Current Size Condition Matrix:"),
                  shiny::verbatimTextOutput("txt_size_matrix"),
                  shiny::actionButton("btn_clear_size_cond", "Clear All", class = "btn-warning btn-sm")
                ),
                shiny::tabPanel("Time Conditions",
                  shiny::br(),
                  shiny::fluidRow(
                    shiny::column(4, shiny::uiOutput("cond_time_par1_ui")),
                    shiny::column(2, shiny::selectInput("cond_time_op", "Operator:",
                      choices = c("<", ">", "="), selected = "<")),
                    shiny::column(4, shiny::uiOutput("cond_time_par2_ui")),
                    shiny::column(2, shiny::br(),
                      shiny::actionButton("btn_add_time_cond", "Add", class = "btn-success"))
                  ),
                  shiny::hr(),
                  shiny::h5("Current Time Condition Matrix:"),
                  shiny::verbatimTextOutput("txt_time_matrix"),
                  shiny::actionButton("btn_clear_time_cond", "Clear All", class = "btn-warning btn-sm")
                )
              )
            )
          )
        ),

        # =====================================================================
        # TAB: Gene Setup
        # =====================================================================
        shinydashboard::tabItem(tabName = "genes",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Genetic Data Configuration", width = 12,
              status = "primary", solidHeader = TRUE,
              shiny::fluidRow(
                shiny::column(4,
                  shiny::selectInput("select_data_type", "Data Source:",
                    choices = c("Manual" = "manual", "From data (PHYLIP/VCF)" = "fromdata"),
                    selected = "manual")
                ),
                shiny::column(4,
                  shiny::conditionalPanel(
                    condition = "input.select_data_type == 'manual'",
                    shiny::numericInput("num_loci", "Number of Loci:", value = 1, min = 1, max = 10000)
                  ),
                  shiny::conditionalPanel(
                    condition = "input.select_data_type == 'fromdata'",
                    shiny::selectInput("select_file_type", "File Format:",
                      choices = c("PHYLIP" = "phylip", "VCF" = "vcf"), selected = "phylip")
                  )
                ),
                shiny::column(4,
                  shiny::conditionalPanel(
                    condition = "input.select_data_type == 'manual'",
                    shiny::selectInput("select_mut_dist", "Mutation Rate Distribution:",
                      choices = c("Uniform" = "uniform", "Normal" = "normal"), selected = "uniform")
                  )
                )
              )
            )
          ),
          # From data: file uploads
          shiny::conditionalPanel(
            condition = "input.select_data_type == 'fromdata'",
            shiny::fluidRow(
              shinydashboard::box(
                title = "Data Files", width = 12,
                status = "info", solidHeader = TRUE,
                shiny::fluidRow(
                  shiny::column(6,
                    shiny::conditionalPanel(
                      condition = "input.select_file_type == 'phylip'",
                      shiny::fileInput("file_phylip", "PHYLIP file:",
                        accept = c(".phy", ".phylip", ".txt"))
                    ),
                    shiny::conditionalPanel(
                      condition = "input.select_file_type == 'vcf'",
                      shiny::fileInput("file_vcf", "VCF file:",
                        accept = c(".vcf", ".vcf.gz"))
                    )
                  ),
                  shiny::column(6,
                    shiny::fileInput("file_pop_assign", "Population assignment (CSV/TSV):",
                      accept = c(".csv", ".tsv", ".txt"))
                  )
                ),
                shiny::conditionalPanel(
                  condition = "input.select_file_type == 'vcf'",
                  shiny::fluidRow(
                    shiny::column(6,
                      shiny::fileInput("file_chrom_sizes", "Chromosome sizes (CSV/TSV):",
                        accept = c(".csv", ".tsv", ".txt"))
                    )
                  )
                ),
                shiny::fluidRow(
                  shiny::column(12,
                    shiny::actionButton("btn_load_data", "Load Data Structure",
                      icon = shiny::icon("upload"), class = "btn-primary"),
                    shiny::br(), shiny::br(),
                    shiny::uiOutput("data_load_status")
                  )
                ),
                shiny::hr(),
                shiny::h5("File Format Examples:"),
                shiny::tags$p(shiny::strong("Population assignment"), " (CSV/TSV, no header):"),
                shiny::tags$pre(
                  "sample1\t1\n",
                  "sample2\t1\n",
                  "sample3\t2\n",
                  "sample4\t2"
                ),
                shiny::tags$p("First column: sample names (must match PHYLIP/VCF). ",
                  "Second column: population number (1, 2, ...)."),
                shiny::conditionalPanel(
                  condition = "input.select_file_type == 'vcf'",
                  shiny::tags$p(shiny::strong("Chromosome sizes"), " (CSV/TSV, no header):"),
                  shiny::tags$pre(
                    "chr1\t100000\n",
                    "chr2\t85000\n",
                    "chr3\t120000"
                  ),
                  shiny::tags$p("First column: chromosome/scaffold name (must match VCF CHROM). ",
                    "Second column: total callable base pairs.")
                )
              )
            )
          ),
          # Manual: loci table
          shiny::conditionalPanel(
            condition = "input.select_data_type == 'manual'",
            shiny::fluidRow(
              shinydashboard::box(
                title = "Locus Parameters", width = 12,
                status = "info", solidHeader = TRUE,
                DT::DTOutput("table_loci"),
                shiny::br(),
                shiny::helpText("Click a cell to edit locus length and mutation rate priors."),
                shiny::helpText(
                  shiny::icon("info-circle"),
                  " Mutation and recombination rates are not set in the model builder.",
                  " They are specified as arguments in the simulation functions",
                  " (e.g.,", shiny::tags$code("mu.rates"), "and",
                  shiny::tags$code("rec.rates"), "in",
                  shiny::tags$code("sim.scrm.stats()"), "or",
                  shiny::tags$code("sim.all.stats()"), ").",
                  " Rates must be", shiny::tags$b("per site per generation"),
                  "(not per year). To convert from per-year rates, multiply by generation time."
                )
              )
            ),
            shiny::fluidRow(
              shinydashboard::box(
                title = "Sample Sizes per Population", width = 12,
                status = "info", solidHeader = TRUE,
                DT::DTOutput("table_samples"),
                shiny::br(),
                shiny::helpText("Click a cell to edit sample sizes (haploid) per population per locus.")
              )
            )
          ),
          # From data: preview tables (shown after loading)
          shiny::conditionalPanel(
            condition = "input.select_data_type == 'fromdata'",
            shiny::fluidRow(
              shinydashboard::box(
                title = "Loaded Data Structure", width = 12,
                status = "success", solidHeader = TRUE,
                DT::DTOutput("table_loaded_loci"),
                shiny::br(),
                DT::DTOutput("table_loaded_samples")
              )
            )
          )
        ),

        # =====================================================================
        # TAB: Visualize & Export
        # =====================================================================
        shinydashboard::tabItem(tabName = "export",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Model Summary", width = 12,
              status = "primary", solidHeader = TRUE,
              shiny::verbatimTextOutput("txt_model_summary")
            )
          ),
          shiny::fluidRow(
            shinydashboard::box(
              title = "Model Visualization", width = 12,
              status = "info", solidHeader = TRUE,
              shiny::fluidRow(
                shiny::column(2,
                  shiny::textInput("plot_title", "Title:",
                    value = "Demographic Model")
                ),
                shiny::column(2,
                  shiny::checkboxInput("check_avg_priors",
                    "Use average of priors (uncheck to sample)", value = TRUE)
                ),
                shiny::column(1,
                  shiny::sliderInput("plot_font_size", "Font size:", min = 0.5, max = 3.0,
                    value = 1.0, step = 0.1)
                ),
                shiny::column(1,
                  shiny::sliderInput("plot_spacing", "Spacing:", min = 0.5, max = 3.0,
                    value = 1.5, step = 0.1)
                ),
                shiny::column(2,
                  shiny::selectInput("select_plot_palette", "Color scheme:",
                    choices = c("Tableau", "Colorblind", "Pastel", "Vivid"),
                    selected = "Colorblind")
                ),
                shiny::column(1,
                  shiny::selectInput("select_plot_font", "Font:",
                    choices = c("Palatino", "Helvetica", "AvantGarde", "Bookman",
                                "Times", "NewCenturySchoolbook"),
                    selected = "Palatino")
                ),
                shiny::column(1,
                  shiny::selectInput("select_plot_bg", "Background:",
                    choices = c("Dark", "Light", "Slate"),
                    selected = "Light")
                ),
                shiny::column(2,
                  shiny::actionButton("btn_plot_model", "Generate Model Plot",
                    icon = shiny::icon("chart-area"), class = "btn-info")
                ),
                shiny::column(2,
                  shiny::downloadButton("btn_save_diagram_pdf", "Save as PDF",
                    class = "btn-sm btn-default", style = "margin-top: 5px;")
                )
              ),
              shiny::hr(),
              shiny::plotOutput("plot_model_diagram", height = "500px")
            )
          ),
          shiny::fluidRow(
            shinydashboard::box(
              title = "Prior Distributions", width = 12,
              status = "warning", solidHeader = TRUE,
              shiny::fluidRow(
                shiny::column(4,
                  shiny::numericInput("num_prior_samples", "Number of samples:",
                    value = 1000, min = 100, max = 50000, step = 100)
                ),
                shiny::column(4,
                  shiny::actionButton("btn_plot_priors", "Plot Priors",
                    icon = shiny::icon("chart-line"), class = "btn-warning")
                )
              ),
              shiny::hr(),
              shiny::uiOutput("plot_priors_ui")
            )
          ),
          shiny::fluidRow(
            shinydashboard::box(
              title = "Export", width = 12,
              status = "success", solidHeader = TRUE,
              shiny::fluidRow(
                shiny::column(4,
                  shiny::actionButton("btn_build_model",
                    "Build Model & Return to R",
                    icon = shiny::icon("check"), class = "btn-success btn-lg")
                ),
                shiny::column(4,
                  shiny::downloadButton("btn_download_model",
                    "Download Model (.rds)", class = "btn-primary btn-lg")
                ),
                shiny::column(4,
                  shiny::actionButton("btn_cancel", "Cancel",
                    icon = shiny::icon("times"), class = "btn-danger btn-lg")
                )
              ),
              shiny::hr(),
              shiny::h5("After building, use in R:"),
              shiny::verbatimTextOutput("txt_usage_code")
            )
          )
        )
      ) # end tabItems
    ) # end dashboardBody
  ) # end dashboardPage

  # ===========================================================================
  # SERVER
  # ===========================================================================
  server <- function(input, output, session) {

    # Serve logo from inst/www
    shiny::addResourcePath("pmlogo", system.file("www", package = "PipeMaster"))

    output$sidebar_logo <- shiny::renderUI({
      shiny::tags$div(
        style = "text-align: center; padding: 10px 0;",
        shiny::tags$img(src = "pmlogo/logo.png", height = "200px")
      )
    })

    # --- Reactive values (replaces .e environment) ---
    rv <- shiny::reactiveValues(
      tree_string   = "1",
      tree_obj      = NULL,
      npops         = 1,
      pop_labels    = "1",
      joints        = NULL,
      ej            = NULL,
      n             = make_cur_ne(1),
      m             = NULL,
      en            = NULL,
      em            = NULL,
      mig_changes   = NULL,
      loci          = NULL,
      I             = NULL,
      ne_changes    = NULL,
      sum_anc_ne    = TRUE,
      en_joins      = NULL,
      size_matrix   = NULL,
      mig_matrix    = NULL,
      time_matrix   = NULL
    )

    # --- Load template if provided ---
    shiny::observe({
      if (!is.null(template_model) && is.list(template_model) && !is.null(template_model$flags)) {
        tmpl <- template_model
        rv$n  <- revert_dist(tmpl$flags$n)
        rv$ej <- revert_dist(tmpl$flags$ej)
        rv$m  <- revert_dist(tmpl$flags$m)
        if (!is.null(tmpl$flags$en)) {
          rv$en <- list(size = revert_dist(tmpl$flags$en$size),
                        time = revert_dist(tmpl$flags$en$time))
        }
        if (!is.null(tmpl$flags$em)) {
          rv$em <- list(size = revert_dist(tmpl$flags$em$size),
                        time = revert_dist(tmpl$flags$em$time))
        }
        rv$loci <- revert_dist(tmpl$loci)
        rv$I    <- tmpl$I
        rv$tree_string <- tmpl$tree
        rv$npops <- nrow(rv$n)
        rv$sum_anc_ne <- if (!is.null(tmpl$sum_anc_ne)) tmpl$sum_anc_ne else TRUE

        # Separate Ne.anc join entries from per-pop en entries
        if (!is.null(tmpl$flags$en)) {
          all_size <- revert_dist(tmpl$flags$en$size)
          all_time <- revert_dist(tmpl$flags$en$time)
          is_join <- grepl("^Ne\\.anc_", all_size[, 1])
          if (any(is_join)) {
            rv$en_joins <- list(size = all_size[is_join, , drop = FALSE],
                                time = all_time[is_join, , drop = FALSE])
          }
          if (any(!is_join)) {
            rv$en <- list(size = all_size[!is_join, , drop = FALSE],
                          time = all_time[!is_join, , drop = FALSE])
          }
        }

        if (!is.null(tmpl$conds$size.matrix)) rv$size_matrix <- tmpl$conds$size.matrix
        if (!is.null(tmpl$conds$mig.matrix))  rv$mig_matrix  <- tmpl$conds$mig.matrix
        if (!is.null(tmpl$conds$time.matrix)) rv$time_matrix <- tmpl$conds$time.matrix

        shiny::updateTextInput(session, "txt_tree",
          value = if (!is.null(tmpl$tree)) tmpl$tree else "1")
        if (!is.null(rv$m))  shiny::updateCheckboxInput(session, "check_migration", value = TRUE)
        if (!is.null(rv$en)) shiny::updateCheckboxInput(session, "check_ancestral_ne", value = TRUE)
        shiny::updateCheckboxInput(session, "check_sum_anc_ne", value = rv$sum_anc_ne)
        if (!is.null(rv$n) && rv$n[1, 6] == "normal")
          shiny::updateSelectInput(session, "select_ne_dist", selected = "normal")
        if (!is.null(rv$ej) && rv$ej[1, 6] == "normal")
          shiny::updateSelectInput(session, "select_time_dist", selected = "normal")
        if (!is.null(rv$m) && rv$m[1, 6] == "normal")
          shiny::updateSelectInput(session, "select_mig_dist", selected = "normal")
      }
    }) |> shiny::bindEvent(session$clientData$url_protocol, once = TRUE)

    # =======================================================================
    # POPULATION STRUCTURE
    # =======================================================================

    # Quick example buttons
    shiny::observeEvent(input$btn_ex_1pop, { shiny::updateTextInput(session, "txt_tree", value = "1") })
    shiny::observeEvent(input$btn_ex_2pop, { shiny::updateTextInput(session, "txt_tree", value = "(1,2)") })
    shiny::observeEvent(input$btn_ex_3pop, { shiny::updateTextInput(session, "txt_tree", value = "(1,(2,3))") })
    shiny::observeEvent(input$btn_ex_4pop, { shiny::updateTextInput(session, "txt_tree", value = "((1,2),(3,4))") })

    # Validate & apply tree
    shiny::observeEvent(input$btn_validate_tree, {
      tree_str <- trimws(input$txt_tree)
      if (tree_str == "") {
        shiny::showNotification("Please enter a tree string.", type = "error")
        return()
      }

      result <- validate_tree(tree_str)
      if (!result$valid) {
        shiny::showNotification(paste("Invalid tree:", result$msg), type = "error", duration = 8)
        return()
      }

      parsed <- parse_tree(tree_str)
      rv$tree_string <- tree_str
      rv$npops       <- parsed$npops
      rv$joints      <- parsed$joints
      rv$ej          <- parsed$ej
      rv$pop_labels  <- parsed$pop_labels

      # Rebuild Ne for new npops
      cur_dist <- input$select_ne_dist
      rv$n <- make_cur_ne(parsed$npops, cur_dist)

      # Reset ancestral Ne (per-pop changes)
      if (isTRUE(input$check_ancestral_ne)) {
        rv$ne_changes <- rep(1, parsed$npops)
        rv$en <- make_anc_ne(rv$ne_changes, cur_dist)
      } else {
        rv$en <- NULL
      }

      # Create ancestral Ne at join nodes
      rv$en_joins <- make_anc_ne_at_joins(parsed$ej, cur_dist)

      # Reset migration
      if (isTRUE(input$check_migration) && parsed$npops > 1) {
        rv$m <- make_mig_par(parsed$npops, input$select_mig_dist)
      } else {
        rv$m <- NULL
      }

      # Reset condition matrices
      rv$size_matrix <- NULL
      rv$mig_matrix  <- NULL
      rv$time_matrix <- NULL

      # Auto-generate time conditions from tree topology:
      # Enforce time ordering only between dependent joins (ancestor-descendant)
      # Join Y depends on join X if Y's pair contains the surviving pop from X
      if (!is.null(parsed$ej) && nrow(parsed$ej) > 1) {
        conds <- build_condition_matrices(rv$n, parsed$ej, rv$en, rv$m, rv$em, rv$en_joins)
        tm <- conds$time.matrix
        ej <- parsed$ej
        for (i in 1:(nrow(ej) - 1)) {
          surv_i <- strsplit(ej[i, 3], " ")[[1]][2]  # surviving pop of join i
          for (j in (i + 1):nrow(ej)) {
            pops_j <- strsplit(ej[j, 3], " ")[[1]]
            if (surv_i %in% pops_j) {
              tm[ej[i, 1], ej[j, 1]] <- "<"  # join i must happen before join j
            }
          }
        }
        tm <- inv_mirror_lower(tm)
        rv$time_matrix <- tm
      }

      # Try to build an ape tree object for visualization
      rv$tree_obj <- tryCatch({
        tree_for_ape <- tree_str
        if (!grepl(";", tree_for_ape)) tree_for_ape <- paste0(tree_for_ape, ";")
        if (tree_str == "1") NULL else ape::read.tree(text = tree_for_ape)
      }, error = function(e) NULL)

      shiny::showNotification(
        paste("Tree applied:", parsed$npops, "populations,",
              if (is.null(parsed$ej)) 0 else nrow(parsed$ej), "junctions"),
        type = "message")
    })

    # Tree validation message
    output$tree_validation_message <- shiny::renderUI({
      shiny::req(input$txt_tree)
      result <- validate_tree(trimws(input$txt_tree))
      if (result$valid) {
        shiny::tags$div(class = "success-message",
          shiny::icon("check-circle"), " ", result$msg,
          sprintf(" (%d populations detected)", rv$npops))
      } else {
        shiny::tags$div(class = "error-message",
          shiny::icon("exclamation-triangle"), " ", result$msg)
      }
    })

    # Model structure plot (reactive to structural changes)
    output$plot_model_structure <- shiny::renderPlot({
      # Reactive dependencies: tree, joins, Ne, en, en_joins, migrations, sum_anc_ne
      rv$tree_string; rv$ej; rv$n; rv$en; rv$en_joins; rv$m; rv$sum_anc_ne; rv$npops
      pal <- if (!is.null(input$select_plot_palette)) input$select_plot_palette else "Colorblind"
      bg_th <- if (!is.null(input$select_plot_bg)) input$select_plot_bg else "Light"
      bg_cols <- c("Dark" = "#222222", "Light" = "#FFFFFF", "Slate" = "#2F4F4F")
      tryCatch({
        plot_model_tree(rv, use_avg = TRUE, font_scale = 1.0, palette = pal, bg_theme = bg_th,
                            spacing = if (!is.null(input$plot_spacing)) input$plot_spacing else 1.5,
                            font_family = if (!is.null(input$select_plot_font)) input$select_plot_font else "Palatino",
                            plot_title = if (!is.null(input$plot_title)) input$plot_title else "Demographic Model")
      }, error = function(e) {
        par(bg = bg_cols[bg_th])
        plot.new()
        fg <- if (bg_th == "Light") "gray40" else "gray50"
        if (rv$npops == 1) {
          text(0.5, 0.5, "Single population", cex = 1.5, col = fg)
        } else {
          text(0.5, 0.5, "Set up population tree to see model",
               cex = 1.2, col = fg)
        }
      })
    }, bg = "#FFFFFF")

    # Save model plot as PDF
    output$btn_save_model_pdf <- shiny::downloadHandler(
      filename = function() {
        paste0("model_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
      },
      content = function(file) {
        pal <- if (!is.null(input$select_plot_palette)) input$select_plot_palette else "Colorblind"
        bg_th <- if (!is.null(input$select_plot_bg)) input$select_plot_bg else "Light"
        bg_cols <- c("Dark" = "#222222", "Light" = "#FFFFFF", "Slate" = "#2F4F4F")
        pdf(file, width = 10, height = 7, bg = bg_cols[bg_th])
        tryCatch(
          plot_model_tree(rv, use_avg = TRUE, font_scale = 1.0, palette = pal, bg_theme = bg_th,
                            spacing = if (!is.null(input$plot_spacing)) input$plot_spacing else 1.5,
                            font_family = if (!is.null(input$select_plot_font)) input$select_plot_font else "Palatino",
                            plot_title = if (!is.null(input$plot_title)) input$plot_title else "Demographic Model"),
          error = function(e) {
            plot.new()
            text(0.5, 0.5, paste("Error:", e$message), col = "red")
          })
        dev.off()
      }
    )

    # Junction table
    output$table_junctions <- DT::renderDT({
      ej <- rv$ej
      if (is.null(ej)) return(DT::datatable(data.frame(Message = "No junctions (single population)")))
      df <- data.frame(
        Junction = ej[, 1],
        Populations = ej[, 3],
        `Default Min` = ej[, 4],
        `Default Max` = ej[, 5],
        Distribution = ej[, 6],
        stringsAsFactors = FALSE, check.names = FALSE
      )
      DT::datatable(df, rownames = FALSE, selection = "none",
        options = list(dom = "t", paging = FALSE))
    })

    # Island model: dynamic checkboxes for node removal
    output$island_node_checkboxes <- shiny::renderUI({
      ej <- rv$ej
      if (is.null(ej) || nrow(ej) == 0) {
        return(shiny::p(shiny::em("No divergence nodes to remove."),
                        style = "color: gray;"))
      }
      choices <- setNames(seq_len(nrow(ej)),
                          paste0(ej[, 1], " (pops ", ej[, 3], ")"))
      shiny::checkboxGroupInput("chk_remove_nodes", "Select nodes to remove:",
                                choices = choices)
    })

    # Island model: remove selected nodes
    shiny::observeEvent(input$btn_remove_nodes, {
      ej <- rv$ej
      if (is.null(ej)) {
        shiny::showNotification("No nodes to remove.", type = "warning")
        return()
      }
      sel <- as.integer(input$chk_remove_nodes)
      if (length(sel) == 0) {
        shiny::showNotification("No nodes selected.", type = "warning")
        return()
      }
      if (length(sel) >= nrow(ej)) {
        # Remove all nodes
        rv$ej <- NULL
        rv$joints <- NULL
      } else {
        rv$ej <- ej[-sel, , drop = FALSE]
        rv$joints <- rv$joints[-sel]
      }
      rv$tree_string <- "non tree-like model"
      rv$tree_obj <- NULL

      # Update en_joins to match remaining ej
      if (!is.null(rv$en_joins) && !is.null(rv$ej)) {
        remaining_pairs <- rv$ej[, 3]
        keep <- sapply(1:nrow(rv$en_joins$size), function(i) {
          pair_key <- sub("^Ne\\.anc_", "", rv$en_joins$size[i, 1])
          any(gsub(" ", "_", remaining_pairs) == pair_key)
        })
        rv$en_joins$size <- rv$en_joins$size[keep, , drop = FALSE]
        rv$en_joins$time <- rv$en_joins$time[keep, , drop = FALSE]
        if (nrow(rv$en_joins$size) == 0) rv$en_joins <- NULL
      } else if (is.null(rv$ej)) {
        rv$en_joins <- NULL
      }

      # Rebuild time condition matrix for remaining parameters (only dependent joins)
      if (!is.null(rv$ej) && nrow(rv$ej) > 1) {
        conds <- build_condition_matrices(rv$n, rv$ej, rv$en, rv$m, rv$em, rv$en_joins)
        tm <- conds$time.matrix
        ej_rem <- rv$ej
        for (i in 1:(nrow(ej_rem) - 1)) {
          surv_i <- strsplit(ej_rem[i, 3], " ")[[1]][2]
          for (j in (i + 1):nrow(ej_rem)) {
            pops_j <- strsplit(ej_rem[j, 3], " ")[[1]]
            if (surv_i %in% pops_j) {
              tm[ej_rem[i, 1], ej_rem[j, 1]] <- "<"
            }
          }
        }
        tm <- inv_mirror_lower(tm)
        rv$time_matrix <- tm
      } else {
        rv$time_matrix <- NULL
      }

      # Add migration between populations that were connected by removed nodes
      # (without joins, disconnected pops would have infinite coalescent time)
      if (rv$npops > 1) {
        dist <- if (!is.null(input$select_mig_dist)) input$select_mig_dist else "uniform"
        removed_ej <- ej[sel, , drop = FALSE]
        new_mig_rows <- list()
        for (k in 1:nrow(removed_ej)) {
          pair <- as.integer(strsplit(removed_ej[k, 3], " ")[[1]])
          p1 <- pair[1]; p2 <- pair[2]
          # Add both directions (i→j and j→i)
          for (dir in list(c(p1, p2), c(p2, p1))) {
            par_name <- paste0("mig0.", dir[1], "_", dir[2])
            # Skip if already exists in rv$m
            if (!is.null(rv$m) && par_name %in% rv$m[, 1]) next
            row <- matrix(c(par_name, "-m", paste(dir[1], dir[2]), "0.1", "1", dist),
                          nrow = 1)
            new_mig_rows[[length(new_mig_rows) + 1]] <- row
          }
        }
        if (length(new_mig_rows) > 0) {
          new_m <- do.call(rbind, new_mig_rows)
          rv$m <- rbind(rv$m, new_m)
          shiny::updateCheckboxInput(session, "check_migration", value = TRUE)
          added_names <- paste(sapply(new_mig_rows, `[`, 1), collapse = ", ")
          shiny::showNotification(
            paste("Migration added for removed node(s):", added_names),
            type = "message", duration = 6)
        }
      }

      shiny::showNotification(
        paste(length(sel), "node(s) removed.",
              if (is.null(rv$ej)) "All nodes removed — pure island model."
              else paste(nrow(rv$ej), "node(s) remaining.")),
        type = "message")
    })

    # =======================================================================
    # DEMOGRAPHY (Ne)
    # =======================================================================

    # Ne distribution change
    shiny::observeEvent(input$select_ne_dist, {
      dist <- input$select_ne_dist
      if (!is.null(rv$n)) rv$n[, 6] <- dist
      if (!is.null(rv$en)) rv$en$size[, 6] <- dist
      if (!is.null(rv$en_joins)) rv$en_joins$size[, 6] <- dist
    })

    # Current Ne table
    output$table_current_ne <- DT::renderDT({
      mat <- rv$n
      if (is.null(mat)) return(NULL)
      df <- mat_to_df(mat, input$select_ne_dist)
      DT::datatable(df, editable = list(target = "cell", disable = list(columns = 0)),
        selection = "none", rownames = FALSE,
        options = list(dom = "t", paging = FALSE, ordering = FALSE))
    })

    shiny::observeEvent(input$table_current_ne_cell_edit, {
      info <- input$table_current_ne_cell_edit
      if (info$col == 1) rv$n[info$row, 4] <- as.character(info$value)
      if (info$col == 2) rv$n[info$row, 5] <- as.character(info$value)
    })

    # Ne changes per population UI
    output$ne_changes_per_pop_ui <- shiny::renderUI({
      np <- rv$npops
      if (np < 1) return(NULL)
      cur <- rv$ne_changes
      inputs <- lapply(1:np, function(i) {
        val <- if (!is.null(cur) && length(cur) >= i) cur[i] else 1
        shiny::numericInput(paste0("ne_ch_pop_", i),
          paste("Pop", i, "- number of Ne changes:"),
          value = val, min = 0, max = 10, step = 1)
      })
      do.call(shiny::tagList, inputs)
    })

    # Apply Ne changes
    shiny::observeEvent(input$btn_apply_ne_changes, {
      np <- rv$npops
      ne_ch <- sapply(1:np, function(i) {
        val <- input[[paste0("ne_ch_pop_", i)]]
        if (is.null(val) || is.na(val)) 1 else val
      })
      rv$ne_changes <- ne_ch
      rv$en <- make_anc_ne(ne_ch, input$select_ne_dist)
      shiny::showNotification("Ancestral Ne parameters updated.", type = "message")
    })

    # Ancestral Ne table
    output$table_ancestral_ne <- DT::renderDT({
      if (is.null(rv$en)) return(NULL)
      df <- mat_to_df(rv$en$size, input$select_ne_dist)
      DT::datatable(df, editable = list(target = "cell", disable = list(columns = 0)),
        selection = "none", rownames = FALSE,
        options = list(dom = "t", paging = FALSE, ordering = FALSE))
    })

    shiny::observeEvent(input$table_ancestral_ne_cell_edit, {
      info <- input$table_ancestral_ne_cell_edit
      if (info$col == 1) rv$en$size[info$row, 4] <- as.character(info$value)
      if (info$col == 2) rv$en$size[info$row, 5] <- as.character(info$value)
    })

    # Handle ancestral Ne toggle off
    shiny::observeEvent(input$check_ancestral_ne, {
      if (!isTRUE(input$check_ancestral_ne)) {
        rv$en <- NULL
        rv$ne_changes <- NULL
      }
    })

    # Sum ancestral Ne checkbox
    shiny::observeEvent(input$check_sum_anc_ne, {
      rv$sum_anc_ne <- isTRUE(input$check_sum_anc_ne)
    })

    # Ancestral Ne at joins table (only shown when sum is OFF)
    output$table_anc_ne_joins <- DT::renderDT({
      ej_ne <- rv$en_joins
      if (is.null(ej_ne)) return(DT::datatable(data.frame(
        Message = "No join events. Set a multi-population tree first.")))
      df <- mat_to_df(ej_ne$size, input$select_ne_dist)
      DT::datatable(df, editable = list(target = "cell", disable = list(columns = 0)),
        selection = "none", rownames = FALSE,
        options = list(dom = "t", paging = FALSE, ordering = FALSE))
    })

    shiny::observeEvent(input$table_anc_ne_joins_cell_edit, {
      info <- input$table_anc_ne_joins_cell_edit
      if (info$col == 1) rv$en_joins$size[info$row, 4] <- as.character(info$value)
      if (info$col == 2) rv$en_joins$size[info$row, 5] <- as.character(info$value)
    })

    # =======================================================================
    # MIGRATION
    # =======================================================================

    shiny::observeEvent(input$check_migration, {
      if (isTRUE(input$check_migration) && rv$npops > 1) {
        # Only create migration if none exists yet (don't overwrite partial migrations)
        if (is.null(rv$m)) {
          rv$m <- make_mig_par(rv$npops, input$select_mig_dist)
        }
      } else {
        rv$m <- NULL
      }
    })

    shiny::observeEvent(input$select_mig_dist, {
      if (!is.null(rv$m)) rv$m[, 6] <- input$select_mig_dist
    })

    output$table_migration <- DT::renderDT({
      if (is.null(rv$m)) return(NULL)
      dist <- if (!is.null(input$select_mig_dist)) input$select_mig_dist else "uniform"
      df <- mat_to_df(rv$m, dist)
      DT::datatable(df, editable = list(target = "cell", disable = list(columns = 0)),
        selection = "none", rownames = FALSE,
        options = list(dom = "t", paging = FALSE, ordering = FALSE))
    })

    shiny::observeEvent(input$table_migration_cell_edit, {
      info <- input$table_migration_cell_edit
      if (info$col == 1) rv$m[info$row, 4] <- as.character(info$value)
      if (info$col == 2) rv$m[info$row, 5] <- as.character(info$value)
    })

    # Migration removal checkboxes
    output$mig_remove_checkboxes <- shiny::renderUI({
      if (is.null(rv$m) || nrow(rv$m) == 0) return(NULL)
      choices <- setNames(1:nrow(rv$m), rv$m[, 1])
      shiny::checkboxGroupInput("chk_remove_mig", NULL, choices = choices, inline = TRUE)
    })

    # Remove selected migration parameters
    shiny::observeEvent(input$btn_remove_mig, {
      if (is.null(rv$m)) return()
      sel <- as.integer(input$chk_remove_mig)
      if (length(sel) == 0) {
        shiny::showNotification("No migration parameters selected.", type = "warning")
        return()
      }
      if (length(sel) >= nrow(rv$m)) {
        shiny::showNotification(
          "Cannot remove all migration parameters — at least one is needed.",
          type = "error", duration = 5)
        return()
      }
      rv$m <- rv$m[-sel, , drop = FALSE]
      shiny::showNotification(
        paste(length(sel), "migration parameter(s) removed.",
              nrow(rv$m), "remaining."),
        type = "message")
    })

    # =======================================================================
    # TIME PRIORS
    # =======================================================================

    shiny::observeEvent(input$select_time_dist, {
      dist <- input$select_time_dist
      if (!is.null(rv$ej)) rv$ej[, 6] <- dist
      if (!is.null(rv$en)) rv$en$time[, 6] <- dist
      if (!is.null(rv$en_joins)) rv$en_joins$time[, 6] <- dist
    })

    # Divergence times info
    output$div_times_ui <- shiny::renderUI({
      if (is.null(rv$ej)) {
        shiny::p("No junctions. Set a multi-population tree in Population Structure first.")
      } else {
        NULL
      }
    })

    # Divergence times table
    output$table_div_times <- DT::renderDT({
      if (is.null(rv$ej)) return(NULL)
      df <- mat_to_df(rv$ej, input$select_time_dist)
      DT::datatable(df, editable = list(target = "cell", disable = list(columns = 0)),
        selection = "none", rownames = FALSE,
        options = list(dom = "t", paging = FALSE, ordering = FALSE))
    })

    shiny::observeEvent(input$table_div_times_cell_edit, {
      info <- input$table_div_times_cell_edit
      if (info$col == 1) rv$ej[info$row, 4] <- as.character(info$value)
      if (info$col == 2) rv$ej[info$row, 5] <- as.character(info$value)
    })

    # Ne change times table
    output$table_ne_change_times <- DT::renderDT({
      if (is.null(rv$en)) return(NULL)
      df <- mat_to_df(rv$en$time, input$select_time_dist)
      DT::datatable(df, editable = list(target = "cell", disable = list(columns = 0)),
        selection = "none", rownames = FALSE,
        options = list(dom = "t", paging = FALSE, ordering = FALSE))
    })

    shiny::observeEvent(input$table_ne_change_times_cell_edit, {
      info <- input$table_ne_change_times_cell_edit
      if (info$col == 1) rv$en$time[info$row, 4] <- as.character(info$value)
      if (info$col == 2) rv$en$time[info$row, 5] <- as.character(info$value)
    })

    # Ancestral migration changes
    shiny::observeEvent(input$check_anc_mig, {
      if (isTRUE(input$check_anc_mig) && !is.null(rv$m)) {
        if (is.null(rv$em)) {
          dist <- if (!is.null(input$select_mig_dist)) input$select_mig_dist else "uniform"
          rv$mig_changes <- rep(1, nrow(rv$m))
          rv$em <- make_anc_mig(rv$m, rv$mig_changes, dist)
        }
      } else {
        rv$em <- NULL
        rv$mig_changes <- NULL
      }
    })

    # Migration changes per pair UI
    output$mig_changes_per_pair_ui <- shiny::renderUI({
      if (is.null(rv$m) || nrow(rv$m) == 0) return(NULL)
      cur <- rv$mig_changes
      inputs <- lapply(1:nrow(rv$m), function(i) {
        val <- if (!is.null(cur) && length(cur) >= i) cur[i] else 1
        shiny::numericInput(paste0("mig_ch_", i),
          paste0(rv$m[i, 1], " - number of changes:"),
          value = val, min = 0, max = 10, step = 1)
      })
      do.call(shiny::tagList, inputs)
    })

    # Apply migration changes
    shiny::observeEvent(input$btn_apply_mig_changes, {
      if (is.null(rv$m)) return()
      n_mig <- nrow(rv$m)
      mig_ch <- sapply(1:n_mig, function(i) {
        val <- input[[paste0("mig_ch_", i)]]
        if (is.null(val) || is.na(val)) 1 else val
      })
      rv$mig_changes <- mig_ch
      dist <- if (!is.null(input$select_mig_dist)) input$select_mig_dist else "uniform"
      rv$em <- make_anc_mig(rv$m, mig_ch, dist)
      shiny::showNotification("Ancestral migration parameters updated.", type = "message")
    })

    output$table_anc_mig_times <- DT::renderDT({
      if (is.null(rv$em)) return(NULL)
      df <- mat_to_df(rv$em$time, input$select_time_dist)
      DT::datatable(df, editable = list(target = "cell", disable = list(columns = 0)),
        selection = "none", rownames = FALSE,
        options = list(dom = "t", paging = FALSE, ordering = FALSE),
        caption = "Migration Change Times")
    })

    shiny::observeEvent(input$table_anc_mig_times_cell_edit, {
      info <- input$table_anc_mig_times_cell_edit
      if (info$col == 1) rv$em$time[info$row, 4] <- as.character(info$value)
      if (info$col == 2) rv$em$time[info$row, 5] <- as.character(info$value)
    })

    output$table_anc_mig_rates <- DT::renderDT({
      if (is.null(rv$em)) return(NULL)
      dist <- if (!is.null(input$select_mig_dist)) input$select_mig_dist else "uniform"
      df <- mat_to_df(rv$em$size, dist)
      DT::datatable(df, editable = list(target = "cell", disable = list(columns = 0)),
        selection = "none", rownames = FALSE,
        options = list(dom = "t", paging = FALSE, ordering = FALSE),
        caption = "Ancestral Migration Rates (4Nm)")
    })

    shiny::observeEvent(input$table_anc_mig_rates_cell_edit, {
      info <- input$table_anc_mig_rates_cell_edit
      if (info$col == 1) rv$em$size[info$row, 4] <- as.character(info$value)
      if (info$col == 2) rv$em$size[info$row, 5] <- as.character(info$value)
    })

    # =======================================================================
    # CONDITIONS
    # =======================================================================

    # Reactively build condition parameter lists
    size_params <- shiny::reactive({
      nms <- rv$n[, 1]
      if (!is.null(rv$en)) nms <- c(nms, rv$en$size[, 1])
      if (!is.null(rv$en_joins)) nms <- c(nms, rv$en_joins$size[, 1])
      nms
    })

    time_params <- shiny::reactive({
      nms <- character()
      if (!is.null(rv$ej)) nms <- c(nms, rv$ej[, 1])
      if (!is.null(rv$en)) nms <- c(nms, rv$en$time[, 1])
      if (!is.null(rv$en_joins)) nms <- c(nms, rv$en_joins$time[, 1])
      nms
    })

    # Condition UI dropdowns
    output$cond_size_par1_ui <- shiny::renderUI({
      shiny::selectInput("cond_size_par1", "Parameter 1:", choices = size_params())
    })
    output$cond_size_par2_ui <- shiny::renderUI({
      shiny::selectInput("cond_size_par2", "Parameter 2:", choices = size_params())
    })
    output$cond_time_par1_ui <- shiny::renderUI({
      shiny::selectInput("cond_time_par1", "Parameter 1:", choices = time_params())
    })
    output$cond_time_par2_ui <- shiny::renderUI({
      shiny::selectInput("cond_time_par2", "Parameter 2:", choices = time_params())
    })

    # Ensure condition matrices exist
    ensure_size_matrix <- function() {
      nms <- size_params()
      if (is.null(rv$size_matrix) || !identical(sort(rownames(rv$size_matrix)), sort(nms))) {
        mat <- matrix(NA_character_, nrow = length(nms), ncol = length(nms))
        rownames(mat) <- nms; colnames(mat) <- nms; diag(mat) <- "0"
        rv$size_matrix <- mat
      }
    }

    ensure_time_matrix <- function() {
      nms <- time_params()
      if (length(nms) == 0) { rv$time_matrix <- NULL; return() }
      if (is.null(rv$time_matrix) || !identical(sort(rownames(rv$time_matrix)), sort(nms))) {
        mat <- matrix(NA_character_, nrow = length(nms), ncol = length(nms))
        rownames(mat) <- nms; colnames(mat) <- nms; diag(mat) <- "0"
        # Preserve existing conditions
        if (!is.null(rv$time_matrix)) {
          shared <- intersect(rownames(rv$time_matrix), nms)
          if (length(shared) > 0) mat[shared, shared] <- rv$time_matrix[shared, shared]
        }
        # Re-apply tree-derived join ordering (only dependent joins)
        if (!is.null(rv$ej) && nrow(rv$ej) > 1) {
          ej <- rv$ej
          for (i in 1:(nrow(ej) - 1)) {
            if (!ej[i, 1] %in% nms) next
            surv_i <- strsplit(ej[i, 3], " ")[[1]][2]
            for (j in (i + 1):nrow(ej)) {
              if (!ej[j, 1] %in% nms) next
              pops_j <- strsplit(ej[j, 3], " ")[[1]]
              if (surv_i %in% pops_j) {
                if (is.na(mat[ej[i, 1], ej[j, 1]]) ||
                    mat[ej[i, 1], ej[j, 1]] == "0") {
                  mat[ej[i, 1], ej[j, 1]] <- "<"
                }
              }
            }
            mat <- inv_mirror_lower(mat)
          }
        }
        rv$time_matrix <- mat
      }
    }

    # Add size condition
    shiny::observeEvent(input$btn_add_size_cond, {
      ensure_size_matrix()
      p1 <- input$cond_size_par1; p2 <- input$cond_size_par2; op <- input$cond_size_op
      if (is.null(p1) || is.null(p2) || p1 == p2) {
        shiny::showNotification("Select two different parameters.", type = "warning"); return()
      }
      rv$size_matrix[p1, p2] <- op
      rv$size_matrix <- inv_mirror_lower(rv$size_matrix)
      shiny::showNotification(paste(p1, op, p2, "added"), type = "message")
    })

    # Add time condition
    shiny::observeEvent(input$btn_add_time_cond, {
      ensure_time_matrix()
      p1 <- input$cond_time_par1; p2 <- input$cond_time_par2; op <- input$cond_time_op
      if (is.null(p1) || is.null(p2) || p1 == p2) {
        shiny::showNotification("Select two different parameters.", type = "warning"); return()
      }
      rv$time_matrix[p1, p2] <- op
      rv$time_matrix <- inv_mirror_lower(rv$time_matrix)
      shiny::showNotification(paste(p1, op, p2, "added"), type = "message")
    })

    # Clear conditions
    shiny::observeEvent(input$btn_clear_size_cond, { rv$size_matrix <- NULL; ensure_size_matrix() })
    shiny::observeEvent(input$btn_clear_time_cond, { rv$time_matrix <- NULL; ensure_time_matrix() })

    # Display condition matrices
    output$txt_size_matrix <- shiny::renderPrint({
      ensure_size_matrix()
      if (!is.null(rv$size_matrix)) print(rv$size_matrix) else cat("No size parameters.")
    })
    output$txt_time_matrix <- shiny::renderPrint({
      ensure_time_matrix()
      if (!is.null(rv$time_matrix)) print(rv$time_matrix) else cat("No time parameters.")
    })

    # =======================================================================
    # GENE SETUP
    # =======================================================================

    # Manual mode: generate loci/I matrices from user inputs
    shiny::observe({
      dtype <- input$select_data_type
      if (is.null(dtype) || dtype != "manual") return()

      np <- rv$npops
      ng <- input$num_loci
      if (is.null(ng) || is.na(ng) || ng < 1) ng <- 1
      mut_dist <- input$select_mut_dist

      loci <- matrix(nrow = ng, ncol = 6)
      loci[, 1] <- paste0("rate", 1:ng)
      loci[, 2] <- "1000"
      loci[, 3] <- "1"
      loci[, 4] <- "5e-9"
      loci[, 5] <- "1.5e-8"
      loci[, 6] <- mut_dist
      rv$loci <- loci

      I_mat <- matrix(nrow = ng, ncol = 3 + np)
      I_mat[, 1] <- paste0("locus", 1:ng)
      I_mat[, 2] <- "-I"
      I_mat[, 3] <- as.character(np)
      for (j in 1:ng) for (i in 1:np) I_mat[j, i + 3] <- "10"
      rv$I <- I_mat
    }) |> shiny::bindEvent(
      input$select_data_type, input$num_loci,
      input$select_mut_dist, rv$npops, ignoreNULL = FALSE)

    # From data: load data structure from PHYLIP/VCF
    shiny::observeEvent(input$btn_load_data, {
      ftype <- input$select_file_type

      # Read pop.assign
      pa_file <- input$file_pop_assign
      if (is.null(pa_file)) {
        output$data_load_status <- shiny::renderUI(
          shiny::tags$span(style = "color: red;", "Please upload a population assignment file."))
        return()
      }
      pop_assign <- tryCatch(
        utils::read.table(pa_file$datapath, header = FALSE, sep = "",
                          stringsAsFactors = FALSE),
        error = function(e) NULL)
      if (is.null(pop_assign) || ncol(pop_assign) < 2) {
        output$data_load_status <- shiny::renderUI(
          shiny::tags$span(style = "color: red;",
            "Could not read pop.assign file. Expected 2 columns: sample_name, population_number."))
        return()
      }

      # Build a temporary model to pass to get.data.structure
      tmp_model <- assemble_model(rv)

      tryCatch({
        if (ftype == "phylip") {
          phy_file <- input$file_phylip
          if (is.null(phy_file)) {
            output$data_load_status <- shiny::renderUI(
              shiny::tags$span(style = "color: red;", "Please upload a PHYLIP file."))
            return()
          }
          tmp_model <- get.data.structure(tmp_model,
            path.to.phylip = phy_file$datapath,
            pop.assign = pop_assign, verbose = FALSE)
        } else {
          vcf_file <- input$file_vcf
          cs_file <- input$file_chrom_sizes
          if (is.null(vcf_file)) {
            output$data_load_status <- shiny::renderUI(
              shiny::tags$span(style = "color: red;", "Please upload a VCF file."))
            return()
          }
          if (is.null(cs_file)) {
            output$data_load_status <- shiny::renderUI(
              shiny::tags$span(style = "color: red;", "Please upload a chromosome sizes file."))
            return()
          }
          chrom_sizes <- utils::read.table(cs_file$datapath, header = FALSE, sep = "",
                                            stringsAsFactors = FALSE)
          tmp_model <- get.data.structure(tmp_model,
            path.to.vcf = vcf_file$datapath,
            pop.assign = pop_assign, chrom.sizes = chrom_sizes, verbose = FALSE)
        }

        # Update rv with loaded data
        rv$loci <- revert_dist(tmp_model$loci)
        rv$I <- tmp_model$I

        nloci <- nrow(tmp_model$I)
        npop <- as.integer(tmp_model$I[1, 3])
        pop_sizes <- paste(tmp_model$I[1, 4:(3 + npop)], collapse = "/")

        output$data_load_status <- shiny::renderUI(
          shiny::tags$span(style = "color: green; font-weight: bold;",
            sprintf("Loaded %d loci, %s haploid samples per pop.", nloci, pop_sizes)))

      }, error = function(e) {
        output$data_load_status <- shiny::renderUI(
          shiny::tags$span(style = "color: red;", paste("Error:", e$message)))
      })
    })

    # From data: preview tables
    output$table_loaded_loci <- DT::renderDT({
      dtype <- input$select_data_type
      if (is.null(dtype) || dtype != "fromdata" || is.null(rv$loci)) return(NULL)
      n_show <- min(nrow(rv$loci), 50)
      df <- data.frame(
        Locus = rv$loci[1:n_show, 1],
        Length_bp = as.numeric(rv$loci[1:n_show, 2]),
        stringsAsFactors = FALSE
      )
      DT::datatable(df, selection = "none", rownames = FALSE,
        options = list(dom = "t", paging = FALSE, ordering = FALSE),
        caption = if (nrow(rv$loci) > 50) paste("Showing first 50 of", nrow(rv$loci), "loci") else NULL)
    })

    output$table_loaded_samples <- DT::renderDT({
      dtype <- input$select_data_type
      if (is.null(dtype) || dtype != "fromdata" || is.null(rv$I)) return(NULL)
      np <- as.integer(rv$I[1, 3])
      n_show <- min(nrow(rv$I), 50)
      df <- data.frame(Locus = rv$I[1:n_show, 1], stringsAsFactors = FALSE)
      for (i in 1:np) {
        df[[paste0("Pop", i)]] <- as.numeric(rv$I[1:n_show, i + 3])
      }
      DT::datatable(df, selection = "none", rownames = FALSE,
        options = list(dom = "t", paging = FALSE, ordering = FALSE),
        caption = if (nrow(rv$I) > 50) paste("Showing first 50 of", nrow(rv$I), "loci") else NULL)
    })

    # Manual: loci table
    output$table_loci <- DT::renderDT({
      if (is.null(rv$loci) || input$select_data_type != "manual") return(NULL)
      mut_dist <- if (!is.null(input$select_mut_dist)) input$select_mut_dist else "uniform"

      df <- data.frame(
        Locus = rv$loci[, 1],
        Length_bp = as.numeric(rv$loci[, 2]),
        Inheritance = as.numeric(rv$loci[, 3]),
        Mut_Min = as.numeric(rv$loci[, 4]),
        Mut_Max = as.numeric(rv$loci[, 5]),
        stringsAsFactors = FALSE
      )
      if (mut_dist == "normal") colnames(df)[4:5] <- c("Mut_Mean", "Mut_SD")
      DT::datatable(df, editable = list(target = "cell", disable = list(columns = 0)),
        selection = "none", rownames = FALSE,
        options = list(dom = "t", paging = FALSE, ordering = FALSE))
    })

    shiny::observeEvent(input$table_loci_cell_edit, {
      info <- input$table_loci_cell_edit
      # cols: 0=Locus(ro), 1=Length, 2=Inheritance, 3=Mut_Min, 4=Mut_Max
      col_map <- c("2", "3", "4", "5")  # DT col 1-4 -> matrix col 2-5
      if (info$col >= 1 && info$col <= 4) {
        rv$loci[info$row, as.integer(col_map[info$col])] <- as.character(info$value)
      }
    })

    # Manual: sample sizes table
    output$table_samples <- DT::renderDT({
      if (is.null(rv$I) || input$select_data_type != "manual") return(NULL)
      np <- rv$npops
      I_mat <- rv$I

      df <- data.frame(Locus = I_mat[, 1], stringsAsFactors = FALSE)
      for (i in 1:np) {
        df[[paste0("Pop", i)]] <- as.numeric(I_mat[, i + 3])
      }
      DT::datatable(df, editable = list(target = "cell", disable = list(columns = 0)),
        selection = "none", rownames = FALSE,
        options = list(dom = "t", paging = FALSE, ordering = FALSE))
    })

    shiny::observeEvent(input$table_samples_cell_edit, {
      info <- input$table_samples_cell_edit
      # col 0 = Locus(ro), col 1..npops = Pop1..PopN -> matrix col 4..3+npops
      if (info$col >= 1) {
        rv$I[info$row, info$col + 3] <- as.character(info$value)
      }
    })

    # =======================================================================
    # STATUS BOXES
    # =======================================================================

    output$status_pops <- shinydashboard::renderValueBox({
      shinydashboard::valueBox(
        value = rv$npops, subtitle = "Populations",
        icon = shiny::icon("sitemap"),
        color = if (rv$npops > 0) "green" else "red")
    })
    output$status_nodes <- shinydashboard::renderValueBox({
      nj <- if (is.null(rv$ej)) 0 else nrow(rv$ej)
      shinydashboard::valueBox(
        value = nj, subtitle = "Divergence Nodes",
        icon = shiny::icon("code-branch"),
        color = if (nj > 0) "blue" else "yellow")
    })
    output$status_loci <- shinydashboard::renderValueBox({
      nl <- if (is.null(rv$loci)) 0 else nrow(rv$loci)
      shinydashboard::valueBox(
        value = nl, subtitle = "Loci",
        icon = shiny::icon("dna"),
        color = if (nl > 0) "green" else "red")
    })

    # =======================================================================
    # MODEL SUMMARY & PLOT
    # =======================================================================

    output$txt_model_summary <- shiny::renderPrint({
      cat("Model Configuration:\n")
      cat("====================\n\n")
      cat(sprintf("Tree:                    %s\n", rv$tree_string))
      cat(sprintf("Populations:             %d\n", rv$npops))
      cat(sprintf("Junctions:               %d\n", if (is.null(rv$ej)) 0 else nrow(rv$ej)))
      cat(sprintf("Migration:               %s\n", !is.null(rv$m)))
      cat(sprintf("Ancestral Ne changes:    %s\n", !is.null(rv$en)))
      cat(sprintf("Current Ne parameters:   %d\n", if (is.null(rv$n)) 0 else nrow(rv$n)))
      cat(sprintf("Ancestral Ne parameters: %d\n", if (is.null(rv$en)) 0 else nrow(rv$en$size)))
      cat(sprintf("Migration parameters:    %d\n", if (is.null(rv$m)) 0 else nrow(rv$m)))
      cat(sprintf("Ne distribution:         %s\n", input$select_ne_dist))
      cat(sprintf("Time distribution:       %s\n", input$select_time_dist))
      cat(sprintf("Data source:             %s\n", input$select_data_type))
      cat(sprintf("Loci:                    %d\n", if (is.null(rv$loci)) 0 else nrow(rv$loci)))
    })

    # Plot model
    shiny::observeEvent(input$btn_plot_model, {
      tryCatch({
        avg <- isTRUE(input$check_avg_priors)
        fs <- if (!is.null(input$plot_font_size)) input$plot_font_size else 1.0
        pal <- if (!is.null(input$select_plot_palette)) input$select_plot_palette else "Tableau"
        bg_th <- if (!is.null(input$select_plot_bg)) input$select_plot_bg else "Dark"
        sp <- if (!is.null(input$plot_spacing)) input$plot_spacing else 1.5
        ff <- if (!is.null(input$select_plot_font)) input$select_plot_font else "Palatino"
        pt <- if (!is.null(input$plot_title)) input$plot_title else "Demographic Model"
        bg_cols <- c("Dark" = "#222222", "Light" = "#FFFFFF", "Slate" = "#2F4F4F")
        bg_col <- bg_cols[bg_th]
        output$plot_model_diagram <- shiny::renderPlot({
          plot_model_tree(rv, use_avg = avg, font_scale = fs, palette = pal, bg_theme = bg_th, spacing = sp, font_family = ff, plot_title = pt)
        }, bg = bg_col)
      }, error = function(e) {
        output$plot_model_diagram <- shiny::renderPlot({
          par(bg = "#222222")
          plot.new()
          text(0.5, 0.5, paste("Plot error:", e$message), col = "red", cex = 1.2)
        }, bg = "#222222")
      })
    })

    # Save model diagram as PDF
    output$btn_save_diagram_pdf <- shiny::downloadHandler(
      filename = function() {
        paste0("model_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
      },
      content = function(file) {
        avg <- isTRUE(input$check_avg_priors)
        fs <- if (!is.null(input$plot_font_size)) input$plot_font_size else 1.0
        pal <- if (!is.null(input$select_plot_palette)) input$select_plot_palette else "Tableau"
        bg_th <- if (!is.null(input$select_plot_bg)) input$select_plot_bg else "Dark"
        sp <- if (!is.null(input$plot_spacing)) input$plot_spacing else 1.5
        ff <- if (!is.null(input$select_plot_font)) input$select_plot_font else "Palatino"
        pt <- if (!is.null(input$plot_title)) input$plot_title else "Demographic Model"
        bg_cols <- c("Dark" = "#222222", "Light" = "#FFFFFF", "Slate" = "#2F4F4F")
        pdf(file, width = 10, height = 7, bg = bg_cols[bg_th])
        tryCatch(
          plot_model_tree(rv, use_avg = avg, font_scale = fs, palette = pal, bg_theme = bg_th, spacing = sp, font_family = ff, plot_title = pt),
          error = function(e) {
            plot.new(); text(0.5, 0.5, paste("Error:", e$message), col = "red")
          }
        )
        dev.off()
      }
    )

    # Dynamic plot height for priors
    rv_prior_height <- shiny::reactiveVal("600px")
    output$plot_priors_ui <- shiny::renderUI({
      shiny::plotOutput("plot_priors", height = rv_prior_height())
    })

    # Plot prior distributions (conditioned on constraints)
    shiny::observeEvent(input$btn_plot_priors, {
      tryCatch({
        fs <- if (!is.null(input$plot_font_size)) input$plot_font_size else 1.0
        nsamp <- input$num_prior_samples
        if (is.null(nsamp) || nsamp < 10) nsamp <- 1000

        # Sample via msABC.commander() — same code path as simulations
        model <- tryCatch(assemble_model(rv), error = function(e) NULL)
        if (is.null(model)) {
          output$plot_priors <- shiny::renderPlot({
            par(bg = "#FFFFFF"); plot.new()
            text(0.5, 0.5, "Model not ready", col = "gray40", cex = 1.5)
          }, bg = "#FFFFFF")
          return()
        }
        if (is.null(model$loci) || !is.matrix(model$loci) || nrow(model$loci) == 0)
          model$loci <- matrix(c("loc1", "100", "1", "1e-8", "1e-8", "runif"), nrow = 1)
        if (is.null(model$I) || !is.matrix(model$I)) {
          samp_per_pop <- rep("10", rv$npops)
          model$I <- matrix(c("I", paste(rv$npops), as.character(rv$npops), samp_per_pop), nrow = 1)
        }

        all_samples <- NULL
        for (s in 1:nsamp) {
          cmd <- tryCatch(
            PipeMaster:::msABC.commander(model, use.alpha = FALSE, arg = "/tmp/prior_plot/"),
            error = function(e) NULL)
          if (is.null(cmd)) next
          p <- cmd[[2]]
          vals <- as.numeric(p[2, ])
          names(vals) <- p[1, ]
          all_samples <- rbind(all_samples, vals)
        }

        if (is.null(all_samples) || ncol(all_samples) == 0) {
          output$plot_priors <- shiny::renderPlot({
            bg_c <- c("Dark" = "#222222", "Light" = "#FFFFFF", "Slate" = "#2F4F4F")
            bg_sel <- if (!is.null(input$select_plot_bg)) input$select_plot_bg else "Light"
            fg_c <- if (bg_sel == "Light") "gray40" else "white"
            par(bg = bg_c[bg_sel])
            plot.new()
            text(0.5, 0.5, "No parameters defined yet", col = fg_c, cex = 1.5)
          }, bg = if (!is.null(input$select_plot_bg) && input$select_plot_bg == "Light") "#FFFFFF" else "#222222")
          return()
        }

        pal_name <- if (!is.null(input$select_plot_palette)) input$select_plot_palette else "Colorblind"
        bg_th <- if (!is.null(input$select_plot_bg)) input$select_plot_bg else "Light"
        ff <- if (!is.null(input$select_plot_font)) input$select_plot_font else "Palatino"
        palettes <- list(
          "Tableau"    = c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2",
                           "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7"),
          "Colorblind" = c("#0072B2", "#E69F00", "#009E73", "#CC79A7",
                           "#56B4E9", "#D55E00", "#F0E442", "#999999"),
          "Pastel"     = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
                           "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5"),
          "Vivid"      = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                           "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
        )
        bg_themes <- list(
          "Dark"  = list(bg = "#222222", fg = "white", axis = "#cccccc", grid = "#666666"),
          "Light" = list(bg = "#FFFFFF", fg = "black", axis = "#333333", grid = "#999999"),
          "Slate" = list(bg = "#2F4F4F", fg = "white", axis = "#cccccc", grid = "#5F7F7F")
        )
        pal <- palettes[[pal_name]]
        th <- bg_themes[[bg_th]]

        np <- ncol(all_samples)
        nc <- min(3, np)
        nr <- ceiling(np / nc)
        plot_h <- paste0(max(300, nr * 200), "px")
        rv_prior_height(plot_h)

        output$plot_priors <- shiny::renderPlot({
          par(mfrow = c(nr, nc), mar = c(3, 3, 2, 0.5),
              bg = th$bg, fg = th$fg,
              col.axis = th$axis, col.lab = th$axis, col.main = th$fg,
              family = ff)

          for (i in 1:np) {
            samp <- all_samples[, i]
            nm <- colnames(all_samples)[i]
            col_i <- pal[((i - 1) %% length(pal)) + 1]
            hist(samp, main = nm, col = adjustcolor(col_i, 0.6),
                 border = adjustcolor(col_i, 0.9),
                 xlab = "", ylab = "", axes = FALSE, cex.main = 1.4 * fs,
                 breaks = 30)
            axis(1, col = th$grid, col.ticks = th$grid, col.axis = th$axis, cex.axis = 1.2 * fs)
            axis(2, col = th$grid, col.ticks = th$grid, col.axis = th$axis, las = 1, cex.axis = 1.2 * fs)
          }
        }, bg = th$bg)
      }, error = function(e) {
        output$plot_priors <- shiny::renderPlot({
          par(bg = "#222222")
          plot.new()
          text(0.5, 0.5, paste("Error:", e$message), col = "red", cex = 1.2)
        }, bg = "#222222")
      })
    })

    # Usage code
    output$txt_usage_code <- shiny::renderText({
      "model <- main.menu.gui()\n\n# Simulate reference table\nsim.sumstats(model, nsim.blocks = 100, block.size = 10,\n             output.name = \"my_sims\")\n\n# Or with existing model as template:\nmodel2 <- main.menu.gui(model)"
    })

    # =======================================================================
    # BUILD / DOWNLOAD / CANCEL
    # =======================================================================

    shiny::observeEvent(input$btn_build_model, {
      model <- assemble_model(rv)
      shiny::stopApp(returnValue = model)
    })

    output$btn_download_model <- shiny::downloadHandler(
      filename = function() { "PipeMaster_model.rds" },
      content = function(file) {
        model <- assemble_model(rv)
        saveRDS(model, file)
      }
    )

    shiny::observeEvent(input$btn_cancel, {
      shiny::stopApp(returnValue = NULL)
    })

  } # end server

  # ===========================================================================
  # LAUNCH APP
  # ===========================================================================
  shiny::runApp(shiny::shinyApp(ui = ui, server = server))
}
