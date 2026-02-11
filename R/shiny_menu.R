#' Shiny GUI for Model Builder
#' @description Launches a Shiny web app as a graphical alternative to main.menu().
#'   Uses a shinydashboard layout with sidebar navigation.
#' @param input An optional existing Model object to load as template. Default NULL.
#' @return A Model object when the user clicks "Build Model", or NULL if cancelled.
#' @examples
#' \dontrun{
#' my.model <- main.menu.gui()
#' sim.sumstat(my.model)
#' }
#' @export
main.menu.gui <- function(input = NULL) {

  # Save the template before Shiny shadows 'input'
  template_model <- input

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

  # Build ancestral Ne matrices (size + time)
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

  # Build migration parameter matrix
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
  build_condition_matrices <- function(n, ej, en, m, em) {
    size_names <- n[, 1]
    time_names <- NULL
    mig_names  <- NULL

    if (!is.null(ej) && is.matrix(ej)) time_names <- ej[, 1]
    if (!is.null(m)  && is.matrix(m))  mig_names  <- m[, 1]
    if (!is.null(en)) {
      size_names <- c(size_names, en$size[, 1])
      time_names <- c(time_names, en$time[, 1])
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

    if (!is.null(rv$en)) {
      en <- list(size = convert_dist(rv$en$size), time = convert_dist(rv$en$time))
      if (is.null(nrow(en$size))) en <- NULL
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

    conds <- build_condition_matrices(rv$n, rv$ej, rv$en, rv$m, rv$em)
    if (!is.null(rv$size_matrix)) conds$size.matrix <- rv$size_matrix
    if (!is.null(rv$mig_matrix))  conds$mig.matrix  <- rv$mig_matrix
    if (!is.null(rv$time_matrix)) conds$time.matrix  <- rv$time_matrix
    model$conds <- conds
    model$tree  <- rv$tree_string

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
                if (!is.null(rv$en)) sample_par_vec(rv$en$size))
    time_v <- c(if (!is.null(rv$ej) && is.matrix(rv$ej)) sample_par_vec(rv$ej),
                if (!is.null(rv$en)) sample_par_vec(rv$en$time),
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
  plot_model_tree <- function(rv, use_avg = TRUE) {
    npops <- rv$npops
    if (is.null(npops) || npops < 1) {
      plot.new(); text(0.5, 0.5, "No populations defined", col = "white", cex = 1.5)
      return(invisible(NULL))
    }

    # Get parameter values with conditions applied
    # use_avg = TRUE: compute conditional means by sampling (respects constraints)
    # use_avg = FALSE: single conditioned sample
    cond_vals <- if (use_avg) conditioned_means(rv) else sample_all_conditioned(rv)
    if (is.null(cond_vals) || length(cond_vals) == 0) {
      plot.new(); text(0.5, 0.5, "No parameters to plot", col = "white", cex = 1.5)
      return(invisible(NULL))
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
    events <- list()
    if (!is.null(rv$ej) && is.matrix(rv$ej) && nrow(rv$ej) > 0) {
      for (i in 1:nrow(rv$ej)) {
        pair <- as.integer(strsplit(rv$ej[i, 3], " ")[[1]])
        events[[length(events) + 1]] <- list(
          type = "join", time = lookup(rv$ej[i, 1]),
          src = pair[1], tgt = pair[2])
      }
    }
    if (!is.null(rv$en) && !is.null(rv$en$size) && is.matrix(rv$en$size)) {
      for (i in 1:nrow(rv$en$size)) {
        events[[length(events) + 1]] <- list(
          type = "ne_change", time = lookup(rv$en$time[i, 1]),
          pop = as.integer(rv$en$size[i, 3]),
          ne = lookup(rv$en$size[i, 1]))
      }
    }
    if (length(events) > 0)
      events <- events[order(sapply(events, `[[`, "time"))]

    # Track state and build drawing primitives
    x_pos  <- setNames(as.numeric(1:npops), as.character(1:npops))
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
          segs[[length(segs) + 1]] <- list(
            pop = s, x = x_pos[s], t0 = last_t[s], t1 = ev$time, ne = cur_ne[s])
          segs[[length(segs) + 1]] <- list(
            pop = g, x = x_pos[g], t0 = last_t[g], t1 = ev$time, ne = cur_ne[g])
          new_x <- mean(c(x_pos[s], x_pos[g]))
          merges[[length(merges) + 1]] <- list(
            time = ev$time,
            x_src = x_pos[s], ne_src = cur_ne[s],
            x_tgt = x_pos[g], ne_tgt = cur_ne[g],
            x_new = new_x, ne_new = cur_ne[g])
          alive[s] <- FALSE
          x_pos[g] <- new_x
          last_t[g] <- ev$time
        }
      }
    }

    # Close remaining active lineages (root extension)
    all_t1  <- if (length(segs) > 0) sapply(segs, `[[`, "t1") else 0
    max_t   <- max(c(all_t1, 1))
    root_ext <- max_t * 0.25
    for (p in names(alive)) {
      if (alive[p]) {
        segs[[length(segs) + 1]] <- list(
          pop = p, x = x_pos[p], t0 = last_t[p],
          t1 = max_t + root_ext, ne = cur_ne[p])
      }
    }

    # Scaling
    all_ne <- sapply(segs, `[[`, "ne")
    y_max  <- max(sapply(segs, `[[`, "t1"))
    ne_max <- max(all_ne); ne_min <- min(all_ne)
    hw_max <- 0.38; hw_min <- 0.06
    ne2hw <- function(ne) {
      if (ne_max == ne_min) return((hw_max + hw_min) / 2)
      hw_min + (hw_max - hw_min) * (ne - ne_min) / (ne_max - ne_min)
    }

    pal <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2",
             "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7",
             "#9C755F", "#BAB0AC")
    pcol <- function(p) pal[((as.integer(p) - 1) %% length(pal)) + 1]

    delta <- y_max * 0.015

    # Plot area
    all_x <- unlist(lapply(segs, function(s) c(s$x - ne2hw(s$ne), s$x + ne2hw(s$ne))))
    x_lo <- min(all_x) - 0.3; x_hi <- max(all_x) + 0.3

    par(mar = c(3.5, 5, 2.5, 1.5), bg = "#222222", fg = "white",
        col.axis = "#cccccc", col.lab = "#cccccc", col.main = "white")
    plot(NULL, xlim = c(x_lo, x_hi), ylim = c(-y_max * 0.07, y_max),
         xlab = "", ylab = "Time (generations ago)",
         main = "Demographic Model", axes = FALSE)
    axis(2, las = 1, col = "#666666", col.ticks = "#666666", col.axis = "#cccccc")

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
              col = adjustcolor("#888888", 0.5), border = NA)
      # Right trapezoid
      polygon(x = c(rx - rhw, rx + rhw, nx + nhw, nx),
              y = c(mg$time, mg$time, mg$time + delta, mg$time + delta),
              col = adjustcolor("#888888", 0.5), border = NA)
    }

    # Draw rectangles
    for (seg in segs) {
      hw <- ne2hw(seg$ne)
      col <- pcol(seg$pop)
      rect(seg$x - hw, seg$t0, seg$x + hw, seg$t1,
           col = adjustcolor(col, alpha.f = 0.75),
           border = adjustcolor(col, alpha.f = 0.9), lwd = 1.5)
    }

    # Population labels
    labels <- if (!is.null(rv$pop_labels)) rv$pop_labels else as.character(1:npops)
    for (i in 1:npops) {
      text(i, -y_max * 0.04, paste0("Pop ", labels[i]),
           col = pcol(as.character(i)), cex = 0.9, font = 2)
    }

    # Ne annotations inside rectangles
    for (seg in segs) {
      seg_h <- seg$t1 - seg$t0
      if (seg_h > y_max * 0.05) {
        text(seg$x, (seg$t0 + seg$t1) / 2,
             format(round(seg$ne), big.mark = ","),
             col = "white", cex = 0.55, srt = 90)
      }
    }

    # Time annotations at merge points
    for (mg in merges) {
      text(x_hi - 0.05, mg$time,
           format(round(mg$time), big.mark = ","),
           col = "#aaaaaa", cex = 0.7, adj = c(1, 0.5))
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
             col = "#FFD700", cex = 0.5)
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
      shinydashboard::sidebarMenu(
        id = "sidebar",
        shinydashboard::menuItem("Getting Started",      tabName = "start",      icon = shiny::icon("play-circle")),
        shinydashboard::menuItem("Population Structure",  tabName = "structure",  icon = shiny::icon("sitemap")),
        shinydashboard::menuItem("Demography (Ne)",       tabName = "demography", icon = shiny::icon("chart-line")),
        shinydashboard::menuItem("Migration",             tabName = "migration",  icon = shiny::icon("exchange-alt")),
        shinydashboard::menuItem("Time Priors",           tabName = "time",       icon = shiny::icon("clock")),
        shinydashboard::menuItem("Conditions",            tabName = "conditions", icon = shiny::icon("filter")),
        shinydashboard::menuItem("Gene Setup",            tabName = "genes",      icon = shiny::icon("dna")),
        shinydashboard::menuItem("Build & Export",        tabName = "export",     icon = shiny::icon("download"))
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
                shiny::tags$li(shiny::strong("Gene Setup:"), " Configure loci and mutation rates"),
                shiny::tags$li(shiny::strong("Build & Export:"), " Generate model object for sim.sumstat()")
              ),
              shiny::hr(),
              shiny::h4("Important: Data Structure"),
              shiny::p("Before simulating data, you must run ",
                shiny::code("get.data.structure()"),
                " on your model object. This function reads your observed FASTA alignments",
                " and a population assignment file to extract the sample sizes per population",
                " and sequence length per locus, updating the model accordingly."),
              shiny::tags$pre(
                "model <- main.menu.gui()\n",
                "model <- get.data.structure(model, path.to.fasta = \"path/to/fastas\",\n",
                "                            pop.assign = my_pop_table)\n",
                "sim.sumstat(model)"
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
              shiny::h4("Tree Visualization:"),
              shiny::plotOutput("plot_tree", height = "400px")
            )
          ),
          shiny::fluidRow(
            shinydashboard::box(
              title = "Detected Junctions (Divergence Events)", width = 6,
              status = "info", solidHeader = TRUE,
              DT::DTOutput("table_junctions")
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
                    "Include ancestral Ne changes through time", value = FALSE)
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
                shiny::column(6,
                  shiny::checkboxInput("check_migration",
                    "Enable migration between populations", value = FALSE)
                ),
                shiny::column(6,
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
                DT::DTOutput("table_migration"),
                shiny::br(),
                shiny::helpText("Click a cell to edit.")
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
                    shiny::p("Enable ancestral Ne changes in the Demography tab first.")
                  )
                ),
                shiny::tabPanel("Migration Change Times",
                  shiny::br(),
                  shiny::p("Migration change times will be available in a future version.")
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
                  shiny::selectInput("select_data_type", "Data Type:",
                    choices = c("Sanger sequencing" = "sanger", "Genomic (NGS)" = "genomic"),
                    selected = "sanger")
                ),
                shiny::column(4,
                  shiny::conditionalPanel(
                    condition = "input.select_data_type == 'sanger'",
                    shiny::numericInput("num_loci", "Number of Loci:", value = 1, min = 1, max = 10000)
                  ),
                  shiny::conditionalPanel(
                    condition = "input.select_data_type == 'genomic'",
                    shiny::numericInput("num_genomic_loci",
                      "Number of genomic loci (include invariable):", value = 1000, min = 1)
                  )
                ),
                shiny::column(4,
                  shiny::selectInput("select_mut_dist", "Mutation Rate Distribution:",
                    choices = c("Uniform" = "uniform", "Normal" = "normal"), selected = "uniform")
                )
              )
            )
          ),
          shiny::fluidRow(
            shinydashboard::box(
              title = "Locus Parameters", width = 12,
              status = "info", solidHeader = TRUE,
              DT::DTOutput("table_loci"),
              shiny::br(),
              shiny::helpText("Click a cell to edit mutation rate priors.")
            )
          ),
          shiny::fluidRow(
            shinydashboard::box(
              title = "Sample Sizes per Population", width = 12,
              status = "info", solidHeader = TRUE,
              DT::DTOutput("table_samples"),
              shiny::br(),
              shiny::helpText("Click a cell to edit sample sizes per population per locus.")
            )
          )
        ),

        # =====================================================================
        # TAB: Build & Export
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
              shiny::checkboxInput("check_avg_priors",
                "Use average of priors (uncheck to sample from priors)", value = TRUE),
              shiny::actionButton("btn_plot_model", "Generate Model Plot",
                icon = shiny::icon("chart-area"), class = "btn-info"),
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
              shiny::plotOutput("plot_priors", height = "600px")
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
      loci          = NULL,
      I             = NULL,
      ne_changes    = NULL,
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
        if (!is.null(tmpl$conds$size.matrix)) rv$size_matrix <- tmpl$conds$size.matrix
        if (!is.null(tmpl$conds$mig.matrix))  rv$mig_matrix  <- tmpl$conds$mig.matrix
        if (!is.null(tmpl$conds$time.matrix)) rv$time_matrix <- tmpl$conds$time.matrix

        shiny::updateTextInput(session, "txt_tree",
          value = if (!is.null(tmpl$tree)) tmpl$tree else "1")
        if (!is.null(rv$m))  shiny::updateCheckboxInput(session, "check_migration", value = TRUE)
        if (!is.null(rv$en)) shiny::updateCheckboxInput(session, "check_ancestral_ne", value = TRUE)
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

      # Reset ancestral Ne
      if (isTRUE(input$check_ancestral_ne)) {
        rv$ne_changes <- rep(1, parsed$npops)
        rv$en <- make_anc_ne(rv$ne_changes, cur_dist)
      } else {
        rv$en <- NULL
      }

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

    # Tree plot
    output$plot_tree <- shiny::renderPlot({
      if (!is.null(rv$tree_obj)) {
        ape::plot.phylo(rv$tree_obj, main = "Population Tree",
          cex = 1.5, font = 2, edge.width = 3, edge.color = "#3c8dbc")
      } else {
        plot.new()
        if (rv$npops == 1) {
          text(0.5, 0.5, "Single population (no tree)", cex = 1.5, col = "gray50")
        } else {
          text(0.5, 0.5, "Click 'Validate & Apply Tree' to see visualization", cex = 1.2, col = "gray50")
        }
      }
    })

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

    # =======================================================================
    # DEMOGRAPHY (Ne)
    # =======================================================================

    # Ne distribution change
    shiny::observeEvent(input$select_ne_dist, {
      dist <- input$select_ne_dist
      if (!is.null(rv$n)) rv$n[, 6] <- dist
      if (!is.null(rv$en)) rv$en$size[, 6] <- dist
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

    # =======================================================================
    # MIGRATION
    # =======================================================================

    shiny::observeEvent(input$check_migration, {
      if (isTRUE(input$check_migration) && rv$npops > 1) {
        rv$m <- make_mig_par(rv$npops, input$select_mig_dist)
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

    # =======================================================================
    # TIME PRIORS
    # =======================================================================

    shiny::observeEvent(input$select_time_dist, {
      dist <- input$select_time_dist
      if (!is.null(rv$ej)) rv$ej[, 6] <- dist
      if (!is.null(rv$en)) rv$en$time[, 6] <- dist
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

    # =======================================================================
    # CONDITIONS
    # =======================================================================

    # Reactively build condition parameter lists
    size_params <- shiny::reactive({
      nms <- rv$n[, 1]
      if (!is.null(rv$en)) nms <- c(nms, rv$en$size[, 1])
      nms
    })

    time_params <- shiny::reactive({
      nms <- character()
      if (!is.null(rv$ej)) nms <- c(nms, rv$ej[, 1])
      if (!is.null(rv$en)) nms <- c(nms, rv$en$time[, 1])
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

    shiny::observe({
      dtype <- input$select_data_type
      np <- rv$npops

      if (dtype == "sanger") {
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

      } else {
        nl <- input$num_genomic_loci
        if (is.null(nl) || is.na(nl) || nl < 1) nl <- 1000
        mut_dist <- input$select_mut_dist

        loci <- matrix(nrow = 1, ncol = 6)
        loci[, 1] <- "genomic"
        loci[, 2] <- "bp"
        loci[, 3] <- as.character(nl)
        loci[, 4] <- "1e-11"
        loci[, 5] <- "1e-9"
        loci[, 6] <- mut_dist
        rv$loci <- loci

        I_mat <- matrix(nrow = 1, ncol = 3 + np)
        I_mat[, 1] <- "genomic"
        I_mat[, 2] <- "-I"
        I_mat[, 3] <- as.character(np)
        for (i in 1:np) I_mat[1, i + 3] <- NA
        rv$I <- I_mat
      }
    }) |> shiny::bindEvent(
      input$select_data_type, input$num_loci, input$num_genomic_loci,
      input$select_mut_dist, rv$npops, ignoreNULL = FALSE)

    # Loci table
    output$table_loci <- DT::renderDT({
      if (is.null(rv$loci)) return(NULL)
      mut_dist <- if (!is.null(input$select_mut_dist)) input$select_mut_dist else "uniform"
      dtype <- input$select_data_type

      if (dtype == "sanger") {
        df <- data.frame(
          Locus = rv$loci[, 1],
          Length_bp = as.numeric(rv$loci[, 2]),
          Inheritance = as.numeric(rv$loci[, 3]),
          Mut_Min = as.numeric(rv$loci[, 4]),
          Mut_Max = as.numeric(rv$loci[, 5]),
          stringsAsFactors = FALSE
        )
        if (mut_dist == "normal") colnames(df)[4:5] <- c("Mut_Mean", "Mut_SD")
      } else {
        df <- data.frame(
          Type = rv$loci[, 1],
          N_Loci = as.numeric(rv$loci[, 3]),
          Mut_Min = as.numeric(rv$loci[, 4]),
          Mut_Max = as.numeric(rv$loci[, 5]),
          stringsAsFactors = FALSE
        )
        if (mut_dist == "normal") colnames(df)[3:4] <- c("Mut_Mean", "Mut_SD")
      }
      DT::datatable(df, editable = list(target = "cell", disable = list(columns = 0)),
        selection = "none", rownames = FALSE,
        options = list(dom = "t", paging = FALSE, ordering = FALSE))
    })

    shiny::observeEvent(input$table_loci_cell_edit, {
      info <- input$table_loci_cell_edit
      dtype <- input$select_data_type
      if (dtype == "sanger") {
        # cols: 0=Locus(ro), 1=Length, 2=Inheritance, 3=Mut_Min, 4=Mut_Max
        col_map <- c("2", "3", "4", "5")  # DT col 1-4 -> matrix col 2-5
        if (info$col >= 1 && info$col <= 4) {
          rv$loci[info$row, as.integer(col_map[info$col])] <- as.character(info$value)
        }
      } else {
        # cols: 0=Type(ro), 1=N_Loci, 2=Mut_Min, 3=Mut_Max
        if (info$col == 1) rv$loci[info$row, 3] <- as.character(info$value)
        if (info$col == 2) rv$loci[info$row, 4] <- as.character(info$value)
        if (info$col == 3) rv$loci[info$row, 5] <- as.character(info$value)
      }
    })

    # Sample sizes table
    output$table_samples <- DT::renderDT({
      if (is.null(rv$I)) return(NULL)
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
      cat(sprintf("Data type:               %s\n", input$select_data_type))
      cat(sprintf("Loci:                    %d\n", if (is.null(rv$loci)) 0 else nrow(rv$loci)))
    })

    # Plot model
    shiny::observeEvent(input$btn_plot_model, {
      tryCatch({
        avg <- isTRUE(input$check_avg_priors)
        output$plot_model_diagram <- shiny::renderPlot({
          plot_model_tree(rv, use_avg = avg)
        }, bg = "#222222")
      }, error = function(e) {
        output$plot_model_diagram <- shiny::renderPlot({
          par(bg = "#222222")
          plot.new()
          text(0.5, 0.5, paste("Plot error:", e$message), col = "red", cex = 1.2)
        }, bg = "#222222")
      })
    })

    # Plot prior distributions (conditioned on constraints)
    shiny::observeEvent(input$btn_plot_priors, {
      tryCatch({
        nsamp <- input$num_prior_samples
        if (is.null(nsamp) || nsamp < 10) nsamp <- 1000

        # Sample nsamp times with conditions applied
        all_samples <- NULL
        for (s in 1:nsamp) {
          v <- sample_all_conditioned(rv)
          all_samples <- rbind(all_samples, v)
        }

        if (is.null(all_samples) || ncol(all_samples) == 0) {
          output$plot_priors <- shiny::renderPlot({
            par(bg = "#222222")
            plot.new()
            text(0.5, 0.5, "No parameters defined yet", col = "white", cex = 1.5)
          }, bg = "#222222")
          return()
        }

        output$plot_priors <- shiny::renderPlot({
          np <- ncol(all_samples)
          nc <- min(3, np)
          nr <- ceiling(np / nc)
          par(mfrow = c(nr, nc), mar = c(3, 3, 2.5, 1),
              bg = "#222222", fg = "white",
              col.axis = "#cccccc", col.lab = "#cccccc", col.main = "white")

          for (i in 1:np) {
            samp <- all_samples[, i]
            nm <- colnames(all_samples)[i]
            d <- stats::density(samp)
            plot(d, main = nm, col = "#F28E2B", lwd = 2,
                 xlab = "", ylab = "", axes = FALSE)
            polygon(d$x, d$y, col = adjustcolor("#F28E2B", 0.3), border = NA)
            axis(1, col = "#666666", col.ticks = "#666666", col.axis = "#cccccc")
            axis(2, col = "#666666", col.ticks = "#666666", col.axis = "#cccccc", las = 1)
          }
        }, bg = "#222222")
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
      "model <- main.menu.gui()\n\n# IMPORTANT: Read data structure before simulating\nmodel <- get.data.structure(model, path.to.fasta = \"path/to/fastas\",\n                            pop.assign = my_pop_table)\n\nsim.sumstat(model)\n\n# Or with existing model as template:\nmodel2 <- main.menu.gui(model)"
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
