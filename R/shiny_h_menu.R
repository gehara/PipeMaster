#' Create an hModel object for hierarchical codemographic analysis
#' @description Bundles all priors and settings needed by \code{sim.coexp.ngs()} into
#'   a single S3 object of class \code{"hModel"}.
#' @param Ne.prior Data frame (nspecies x 4): species name, distribution, min Ne, max Ne.
#' @param NeA.prior Data frame (nspecies x 4): species name, distribution, min ancestral Ne ratio, max ancestral Ne ratio.
#' @param time.prior Data frame (nspecies x 5): species name, distribution, min time, max time, generation time.
#' @param gene.prior List of data frames (one per species). Each data frame is a locfile with columns: id, n, pop, length.
#' @param coexp.prior Numeric vector of length 2: c(min, max) uniform prior on coexpansion time.
#' @param var.zeta "FREE" or a numeric value (0-1). Proportion of coexpanding species.
#' @param th Numeric. Threshold for minimum time difference between Ts and population-specific times. Default 0.
#' @param mu.rates List of per-species mutation rate distribution specs. Each element is a list like
#'   \code{list("rnorm", NA, 1e-8, 1.5e-9)}. Must have one element per species.
#' @param alpha Logical. If TRUE, demographic changes are exponential. If FALSE, sudden. Default FALSE.
#' @param phylip_paths Optional named list mapping species names to file paths of multi-locus PHYLIP files used to build gene.prior. Stored for reference and GUI template reload. Default NULL.
#' @return An object of class \code{"hModel"}.
#' @examples
#' \dontrun{
#' Ne   <- data.frame(species=c("sp1","sp2"), dist="uniform", min=10000, max=500000)
#' NeA  <- data.frame(species=c("sp1","sp2"), dist="uniform", min=0.01, max=10)
#' Time <- data.frame(species=c("sp1","sp2"), dist="uniform", min=10000, max=500000, gen=1)
#' gp   <- list(locfile_sp1, locfile_sp2)  # data.frames with columns: id, n, pop, length
#' mu   <- list(list("rnorm", NA, 1e-8, 1.5e-9),
#'              list("rnorm", NA, 1e-8, 1.5e-9))
#' hm   <- hModel(Ne, NeA, Time, gp, coexp.prior=c(10000,500000), mu.rates=mu)
#' sim.coexp.ngs(hm, nsims=100)
#' }
#' @noRd
hModel <- function(Ne.prior, NeA.prior, time.prior, gene.prior,
                   coexp.prior, var.zeta = "FREE", th = 0,
                   mu.rates, alpha = FALSE, phylip_paths = NULL,
                   model = c("single", "PT", "bins")) {
  model <- match.arg(model)

  # --- Input validation ---
  if (!is.data.frame(Ne.prior))
    stop("Ne.prior must be a data.frame")
  if (ncol(Ne.prior) != 4)
    stop("Ne.prior must have 4 columns: species, distribution, min, max")

  if (!is.data.frame(NeA.prior))
    stop("NeA.prior must be a data.frame")
  if (ncol(NeA.prior) != 4)
    stop("NeA.prior must have 4 columns: species, distribution, min, max")

  if (!is.data.frame(time.prior))
    stop("time.prior must be a data.frame")
  if (ncol(time.prior) != 5)
    stop("time.prior must have 5 columns: species, distribution, min, max, generation_time")

  nsp <- nrow(Ne.prior)

  if (nrow(NeA.prior) != nsp)
    stop("NeA.prior must have the same number of rows as Ne.prior")
  if (nrow(time.prior) != nsp)
    stop("time.prior must have the same number of rows as Ne.prior")

  if (!is.list(gene.prior))
    stop("gene.prior must be a list of data.frames (one per species)")
  if (length(gene.prior) != nsp)
    stop("gene.prior must have one element per species (", nsp, " expected)")

  if (!is.numeric(coexp.prior) || length(coexp.prior) != 2)
    stop("coexp.prior must be a numeric vector of length 2: c(min, max)")
  if (coexp.prior[1] > coexp.prior[2])
    stop("coexp.prior[1] (min) must be <= coexp.prior[2] (max)")

  if (!identical(var.zeta, "FREE")) {
    if (!is.numeric(var.zeta) || length(var.zeta) != 1 || var.zeta < 0 || var.zeta > 1)
      stop("var.zeta must be 'FREE' or a single numeric value between 0 and 1")
  }

  if (!is.numeric(th) || length(th) != 1 || th < 0)
    stop("th must be a non-negative numeric value")

  if (!is.list(mu.rates))
    stop("mu.rates must be a list of per-species mutation rate distribution specs")
  if (length(mu.rates) != nsp)
    stop("mu.rates must have one element per species (", nsp, " expected, got ", length(mu.rates), ")")
  for (i in seq_len(nsp)) {
    if (!is.list(mu.rates[[i]]) || length(mu.rates[[i]]) < 3)
      stop("mu.rates[[", i, "]] must be a list with at least 3 elements: distribution name, NA, and distribution parameters")
  }

  if (!is.logical(alpha) || length(alpha) != 1)
    stop("alpha must be a single logical value (TRUE or FALSE)")

  # --- Build object ---
  obj <- list(
    Ne.prior    = Ne.prior,
    NeA.prior   = NeA.prior,
    time.prior  = time.prior,
    gene.prior  = gene.prior,
    coexp.prior = coexp.prior,
    var.zeta    = var.zeta,
    th          = th,
    mu.rates    = mu.rates,
    alpha       = alpha,
    nspecies    = nsp,
    species     = as.character(Ne.prior[, 1]),
    phylip_paths = phylip_paths,
    model       = model
  )

  class(obj) <- "hModel"
  return(obj)
}

#' Print method for hModel objects
#' @param x An hModel object.
#' @param ... Additional arguments (ignored).
#' @noRd
print.hModel <- function(x, ...) {
  cat("Hierarchical codemographic model (hModel)\n")
  cat("==========================================\n\n")
  cat(sprintf("Species (%d):  %s\n", x$nspecies, paste(x$species, collapse = ", ")))
  cat(sprintf("var.zeta:     %s\n", as.character(x$var.zeta)))
  cat(sprintf("coexp.prior:  [%s, %s]\n", x$coexp.prior[1], x$coexp.prior[2]))
  cat(sprintf("th:           %s\n", x$th))
  cat(sprintf("alpha:        %s\n", x$alpha))
  cat("\nPer-species details:\n")
  for (i in seq_along(x$gene.prior)) {
    mr <- x$mu.rates[[i]]
    mu_str <- sprintf("%s(%s)", mr[[1]], paste(mr[3:length(mr)], collapse = ", "))
    phy <- if (!is.null(x$phylip_paths[[x$species[i]]])) paste0(" [", x$phylip_paths[[x$species[i]]], "]") else ""
    cat(sprintf("  %s: %d loci, mu.rates=%s%s\n", x$species[i], nrow(x$gene.prior[[i]]), mu_str, phy))
  }
  invisible(x)
}

#' Shiny GUI for building hModel objects
#' @description Launches a Shiny web app for interactively building hierarchical
#'   codemographic models (\code{hModel} objects) for use with \code{sim.coexp.ngs()}.
#'   Uses a shinydashboard layout with sidebar navigation.
#' @param input An optional existing hModel object to load as template. Default NULL.
#' @return An hModel object when the user clicks "Build Model", or NULL if cancelled.
#' @examples
#' \dontrun{
#' hm <- h.menu.gui()
#' sim.coexp.ngs(hm, nsims = 100)
#' }
#' @noRd
h.menu.gui <- function(input = NULL) {

  # Save the template before Shiny shadows 'input'
  template_hmodel <- input

  # Check for required packages
  for (pkg in c("shinydashboard", "shinyjs", "DT")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(paste0("Package '", pkg, "' is required. Install with: install.packages('", pkg, "')"))
  }

  # ===========================================================================
  # HELPER FUNCTIONS
  # ===========================================================================

  # Build a default locfile data.frame for a species
  make_locfile <- function(nloci, nsamp, bp_length) {
    data.frame(
      id     = paste0("locus", seq_len(nloci)),
      n      = rep(as.integer(nsamp), nloci),
      pop    = rep(1L, nloci),
      length = rep(as.integer(bp_length), nloci),
      stringsAsFactors = FALSE
    )
  }

  # Build a locfile data.frame from a parsed multi-locus phylip file
  # loci: list of character matrices (one per locus) from read.phylip.loci()
  locfile_from_phylip <- function(loci) {
    nloci <- length(loci)
    data.frame(
      id     = paste0("locus", seq_len(nloci)),
      n      = vapply(loci, nrow, integer(1)),
      pop    = rep(1L, nloci),
      length = vapply(loci, ncol, integer(1)),
      stringsAsFactors = FALSE
    )
  }

  # Assemble hModel from reactive values
  assemble_hmodel <- function(rv) {
    nsp <- rv$nspecies
    sp_names <- rv$species_names

    # Build Ne.prior data.frame
    Ne.prior <- data.frame(
      species = sp_names,
      dist    = rep("uniform", nsp),
      min     = rv$ne_min,
      max     = rv$ne_max,
      stringsAsFactors = FALSE
    )

    # Build NeA.prior data.frame
    NeA.prior <- data.frame(
      species = sp_names,
      dist    = rep("uniform", nsp),
      min     = rv$nea_min,
      max     = rv$nea_max,
      stringsAsFactors = FALSE
    )

    # Build time.prior data.frame
    time.prior <- data.frame(
      species = sp_names,
      dist    = rep("uniform", nsp),
      min     = rv$time_min,
      max     = rv$time_max,
      gen     = rv$gen_time,
      stringsAsFactors = FALSE
    )

    # gene.prior (list of locfiles)
    gene.prior <- rv$gene_prior

    # coexp.prior
    coexp.prior <- c(rv$coexp_min, rv$coexp_max)

    # var.zeta
    var.zeta <- rv$var_zeta

    # th
    th <- rv$th_val

    # mu.rates (per-species)
    mu.rates <- lapply(seq_len(nsp), function(i) {
      list(rv$mu_dist[[i]], NA, rv$mu_p1[[i]], rv$mu_p2[[i]])
    })

    # alpha
    alpha <- rv$alpha_val

    # Collect phylip paths (non-NULL entries only)
    pp <- rv$phylip_paths
    if (length(pp) == 0) pp <- NULL

    hModel(Ne.prior = Ne.prior, NeA.prior = NeA.prior, time.prior = time.prior,
           gene.prior = gene.prior, coexp.prior = coexp.prior,
           var.zeta = var.zeta, th = th, mu.rates = mu.rates, alpha = alpha,
           phylip_paths = pp)
  }

  # ===========================================================================
  # UI
  # ===========================================================================
  ui <- shinydashboard::dashboardPage(
    skin = "purple",

    shinydashboard::dashboardHeader(
      title = "PipeMaster - hModel Builder",
      titleWidth = 320
    ),

    shinydashboard::dashboardSidebar(
      width = 280,
      shinydashboard::sidebarMenu(
        id = "sidebar",
        shinydashboard::menuItem("Getting Started",      tabName = "start",       icon = shiny::icon("play-circle")),
        shinydashboard::menuItem("Species Setup",        tabName = "species",     icon = shiny::icon("paw")),
        shinydashboard::menuItem("Ne Priors",            tabName = "ne_priors",   icon = shiny::icon("chart-bar")),
        shinydashboard::menuItem("Ancestral Ne",         tabName = "nea_priors",  icon = shiny::icon("chart-line")),
        shinydashboard::menuItem("Time Priors",          tabName = "time_priors", icon = shiny::icon("clock")),
        shinydashboard::menuItem("Gene/Locus Setup",     tabName = "gene_setup",  icon = shiny::icon("dna")),
        shinydashboard::menuItem("Hierarchical Params",  tabName = "hier_params", icon = shiny::icon("layer-group")),
        shinydashboard::menuItem("Build & Export",       tabName = "export",      icon = shiny::icon("download"))
      ),
      shiny::hr(),
      shinydashboard::box(
        title = "Model Status", width = 12, solidHeader = FALSE,
        shinydashboard::valueBoxOutput("hm_status_species", width = 12),
        shinydashboard::valueBoxOutput("hm_status_loci",    width = 12)
      )
    ),

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
        # TAB 1: Getting Started
        # =====================================================================
        shinydashboard::tabItem(tabName = "start",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Welcome to the hModel Builder", width = 12,
              status = "primary", solidHeader = TRUE,
              shiny::h3("Build a hierarchical codemographic model"),
              shiny::p("This GUI creates an hModel object for use with sim.coexp.ngs()."),
              shiny::hr(),
              shiny::h4("Workflow:"),
              shiny::tags$ol(
                shiny::tags$li(shiny::strong("Species Setup:"), " Define number and names of species"),
                shiny::tags$li(shiny::strong("Ne Priors:"), " Set effective population size priors per species"),
                shiny::tags$li(shiny::strong("Ancestral Ne:"), " Set ancestral Ne ratio priors per species"),
                shiny::tags$li(shiny::strong("Time Priors:"), " Set time of demographic change priors + generation times"),
                shiny::tags$li(shiny::strong("Gene/Locus Setup:"), " Configure loci per species (nloci, samples, bp length)"),
                shiny::tags$li(shiny::strong("Hierarchical Params:"), " Coexpansion prior, var.zeta, threshold, alpha, mutation rates"),
                shiny::tags$li(shiny::strong("Build & Export:"), " Review, build, download, or return to R")
              ),
              shiny::hr(),
              shiny::h4("Usage example:"),
              shiny::tags$pre(
                "hm <- h.menu.gui()\n",
                "sim.coexp.ngs(hm, nsims = 1000)\n\n",
                "# Or override specific parameters:\n",
                "sim.coexp.ngs(hm, nsims = 1000, th = 100)"
              )
            )
          )
        ),

        # =====================================================================
        # TAB 2: Species Setup
        # =====================================================================
        shinydashboard::tabItem(tabName = "species",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Species Configuration", width = 12,
              status = "primary", solidHeader = TRUE,
              shiny::fluidRow(
                shiny::column(4,
                  shiny::numericInput("hm_nspecies", "Number of species:",
                    value = 2, min = 2, max = 100, step = 1)
                ),
                shiny::column(4,
                  shiny::br(),
                  shiny::actionButton("hm_btn_apply_species", "Apply",
                    icon = shiny::icon("check"), class = "btn-success")
                )
              ),
              shiny::hr(),
              shiny::h4("Species Names:"),
              shiny::p("Edit the species names below. Click 'Apply' above first to set the number of species."),
              shiny::uiOutput("hm_species_names_ui"),
              shiny::br(),
              shiny::actionButton("hm_btn_update_names", "Update Names",
                icon = shiny::icon("sync"), class = "btn-primary")
            )
          )
        ),

        # =====================================================================
        # TAB 3: Ne Priors
        # =====================================================================
        shinydashboard::tabItem(tabName = "ne_priors",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Effective Population Size (Ne) Priors", width = 12,
              status = "primary", solidHeader = TRUE,
              shiny::p("Uniform prior bounds for Ne of each species."),
              DT::DTOutput("hm_table_ne"),
              shiny::br(),
              shiny::helpText("Click a cell to edit. Species column is read-only.")
            )
          )
        ),

        # =====================================================================
        # TAB 4: Ancestral Ne
        # =====================================================================
        shinydashboard::tabItem(tabName = "nea_priors",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Ancestral Ne Ratio Priors", width = 12,
              status = "primary", solidHeader = TRUE,
              shiny::p("Uniform prior bounds for the ancestral Ne ratio (NeA/Ne) of each species."),
              DT::DTOutput("hm_table_nea"),
              shiny::br(),
              shiny::helpText("Click a cell to edit. Species column is read-only.")
            )
          )
        ),

        # =====================================================================
        # TAB 5: Time Priors
        # =====================================================================
        shinydashboard::tabItem(tabName = "time_priors",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Time of Demographic Change Priors", width = 12,
              status = "primary", solidHeader = TRUE,
              shiny::p("Uniform prior bounds for time of demographic change (in years) and generation time per species."),
              DT::DTOutput("hm_table_time"),
              shiny::br(),
              shiny::helpText("Click a cell to edit. Species column is read-only.")
            )
          )
        ),

        # =====================================================================
        # TAB 6: Gene/Locus Setup
        # =====================================================================
        shinydashboard::tabItem(tabName = "gene_setup",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Load Data from Phylip File", width = 12,
              status = "success", solidHeader = TRUE,
              shiny::p("Point to a multi-locus sequential PHYLIP file for each species.",
                " The file is parsed to extract locus count, bp lengths, and sample sizes",
                " (same format as obs.sumstat.ngs())."),
              shiny::fluidRow(
                shiny::column(3,
                  shiny::selectInput("hm_gene_species", "Select species:", choices = NULL)
                ),
                shiny::column(6,
                  shiny::textInput("hm_phylip_path", "Path to PHYLIP file:",
                    placeholder = "/path/to/species_data.phy")
                ),
                shiny::column(3,
                  shiny::br(),
                  shiny::actionButton("hm_btn_load_phylip", "Load from Phylip",
                    icon = shiny::icon("file-import"), class = "btn-success")
                )
              ),
              shiny::uiOutput("hm_phylip_status")
            )
          ),
          shiny::fluidRow(
            shinydashboard::box(
              title = "Manual Configuration (alternative)", width = 12,
              status = "warning", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
              shiny::p("Use this if you don't have a phylip file yet, or want to override values."),
              shiny::fluidRow(
                shiny::column(3,
                  shiny::numericInput("hm_gene_nloci", "Number of loci:", value = 10, min = 1, max = 100000)
                ),
                shiny::column(3,
                  shiny::numericInput("hm_gene_nsamp", "Sample size (n):", value = 10, min = 2, max = 1000)
                ),
                shiny::column(3,
                  shiny::numericInput("hm_gene_bp", "BP length:", value = 500, min = 1, max = 100000)
                )
              ),
              shiny::fluidRow(
                shiny::column(3,
                  shiny::actionButton("hm_btn_gene_apply", "Apply to selected species",
                    icon = shiny::icon("check"), class = "btn-info")
                ),
                shiny::column(3,
                  shiny::actionButton("hm_btn_gene_apply_all", "Apply to ALL species",
                    icon = shiny::icon("sync"), class = "btn-warning")
                )
              )
            )
          ),
          shiny::fluidRow(
            shinydashboard::box(
              title = "Locfile Preview", width = 12,
              status = "info", solidHeader = TRUE,
              shiny::h4(shiny::textOutput("hm_locfile_title")),
              DT::DTOutput("hm_table_locfile_preview")
            )
          )
        ),

        # =====================================================================
        # TAB 7: Hierarchical Params
        # =====================================================================
        shinydashboard::tabItem(tabName = "hier_params",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Coexpansion Prior", width = 6,
              status = "primary", solidHeader = TRUE,
              shiny::numericInput("hm_coexp_min", "Coexpansion time min:", value = 10000, min = 0),
              shiny::numericInput("hm_coexp_max", "Coexpansion time max:", value = 500000, min = 0)
            ),
            shinydashboard::box(
              title = "Zeta Parameter", width = 6,
              status = "primary", solidHeader = TRUE,
              shiny::radioButtons("hm_zeta_type", "var.zeta mode:",
                choices = c("FREE (uniform over possible values)" = "FREE",
                            "Fixed value" = "fixed"),
                selected = "FREE"),
              shiny::conditionalPanel(
                condition = "input.hm_zeta_type == 'fixed'",
                shiny::numericInput("hm_zeta_val", "Fixed zeta value (0-1):",
                  value = 0.5, min = 0, max = 1, step = 0.05)
              )
            )
          ),
          shiny::fluidRow(
            shinydashboard::box(
              title = "Threshold & Alpha", width = 6,
              status = "info", solidHeader = TRUE,
              shiny::numericInput("hm_th", "Threshold (th):", value = 0, min = 0),
              shiny::helpText("Minimum time difference between Ts and population-specific times. Set to 0 for the Chan et al. (2014) model."),
              shiny::checkboxInput("hm_alpha", "Exponential demographic changes (alpha=TRUE)", value = FALSE),
              shiny::helpText("If unchecked, sudden demographic changes are used (default).")
            ),
            shinydashboard::box(
              title = "Mutation Rate Distribution (per species)", width = 6,
              status = "info", solidHeader = TRUE,
              shiny::p("Each species has its own mutation rate prior distribution.",
                " For rnorm: param1=mean, param2=SD. For runif: param1=min, param2=max.",
                " For rlnorm: param1=meanlog, param2=sdlog."),
              DT::DTOutput("hm_table_mu"),
              shiny::br(),
              shiny::helpText("Click a cell to edit. Species column is read-only.")
            )
          )
        ),

        # =====================================================================
        # TAB 8: Build & Export
        # =====================================================================
        shinydashboard::tabItem(tabName = "export",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Model Summary", width = 12,
              status = "primary", solidHeader = TRUE,
              shiny::verbatimTextOutput("hm_txt_summary")
            )
          ),
          shiny::fluidRow(
            shinydashboard::box(
              title = "Export", width = 12,
              status = "success", solidHeader = TRUE,
              shiny::fluidRow(
                shiny::column(4,
                  shiny::actionButton("hm_btn_build",
                    "Build hModel & Return to R",
                    icon = shiny::icon("check"), class = "btn-success btn-lg")
                ),
                shiny::column(4,
                  shiny::downloadButton("hm_btn_download",
                    "Download hModel (.rds)", class = "btn-primary btn-lg")
                ),
                shiny::column(4,
                  shiny::actionButton("hm_btn_cancel", "Cancel",
                    icon = shiny::icon("times"), class = "btn-danger btn-lg")
                )
              ),
              shiny::hr(),
              shiny::h5("After building, use in R:"),
              shiny::verbatimTextOutput("hm_txt_usage")
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

    # --- Reactive values ---
    rv <- shiny::reactiveValues(
      nspecies      = 2,
      species_names = c("sp1", "sp2"),
      ne_min        = c(10000, 10000),
      ne_max        = c(500000, 500000),
      nea_min       = c(0.01, 0.01),
      nea_max       = c(10, 10),
      time_min      = c(10000, 10000),
      time_max      = c(500000, 500000),
      gen_time      = c(1, 1),
      gene_prior    = list(
        make_locfile(10, 10, 500),
        make_locfile(10, 10, 500)
      ),
      coexp_min     = 10000,
      coexp_max     = 500000,
      var_zeta      = "FREE",
      th_val        = 0,
      alpha_val     = FALSE,
      mu_dist       = c("rnorm", "rnorm"),
      mu_p1         = c(1e-8, 1e-8),
      mu_p2         = c(1.5e-9, 1.5e-9),
      phylip_paths  = list()
    )

    # --- Load template if provided ---
    shiny::observe({
      if (!is.null(template_hmodel) && inherits(template_hmodel, "hModel")) {
        tmpl <- template_hmodel
        nsp <- tmpl$nspecies

        rv$nspecies      <- nsp
        rv$species_names <- tmpl$species
        rv$ne_min        <- as.numeric(tmpl$Ne.prior[, 3])
        rv$ne_max        <- as.numeric(tmpl$Ne.prior[, 4])
        rv$nea_min       <- as.numeric(tmpl$NeA.prior[, 3])
        rv$nea_max       <- as.numeric(tmpl$NeA.prior[, 4])
        rv$time_min      <- as.numeric(tmpl$time.prior[, 3])
        rv$time_max      <- as.numeric(tmpl$time.prior[, 4])
        rv$gen_time      <- as.numeric(tmpl$time.prior[, 5])
        rv$gene_prior    <- tmpl$gene.prior
        if (!is.null(tmpl$phylip_paths)) rv$phylip_paths <- tmpl$phylip_paths
        rv$coexp_min     <- tmpl$coexp.prior[1]
        rv$coexp_max     <- tmpl$coexp.prior[2]
        rv$var_zeta      <- tmpl$var.zeta
        rv$th_val        <- tmpl$th
        rv$alpha_val     <- tmpl$alpha
        rv$mu_dist       <- vapply(tmpl$mu.rates, function(m) m[[1]], character(1))
        rv$mu_p1         <- vapply(tmpl$mu.rates, function(m) m[[3]], numeric(1))
        rv$mu_p2         <- vapply(tmpl$mu.rates, function(m) m[[4]], numeric(1))

        shiny::updateNumericInput(session, "hm_nspecies", value = nsp)
        shiny::updateNumericInput(session, "hm_coexp_min", value = tmpl$coexp.prior[1])
        shiny::updateNumericInput(session, "hm_coexp_max", value = tmpl$coexp.prior[2])
        shiny::updateNumericInput(session, "hm_th", value = tmpl$th)
        shiny::updateCheckboxInput(session, "hm_alpha", value = tmpl$alpha)

        if (identical(tmpl$var.zeta, "FREE")) {
          shiny::updateRadioButtons(session, "hm_zeta_type", selected = "FREE")
        } else {
          shiny::updateRadioButtons(session, "hm_zeta_type", selected = "fixed")
          shiny::updateNumericInput(session, "hm_zeta_val", value = tmpl$var.zeta)
        }
      }
    }) |> shiny::bindEvent(session$clientData$url_protocol, once = TRUE)

    # =======================================================================
    # SPECIES SETUP
    # =======================================================================

    shiny::observeEvent(input$hm_btn_apply_species, {
      nsp <- input$hm_nspecies
      if (is.null(nsp) || is.na(nsp) || nsp < 2) {
        shiny::showNotification("Need at least 2 species.", type = "error")
        return()
      }

      old_nsp <- rv$nspecies
      rv$nspecies <- nsp

      # Expand or shrink vectors
      expand_vec <- function(old, n, default) {
        if (length(old) >= n) return(old[1:n])
        c(old, rep(default, n - length(old)))
      }

      rv$species_names <- expand_vec(rv$species_names, nsp, paste0("sp", seq_len(nsp))[(length(rv$species_names)+1):nsp])
      # Fix names to ensure proper sp1, sp2 etc.
      if (nsp > old_nsp) {
        rv$species_names <- c(rv$species_names[1:old_nsp], paste0("sp", (old_nsp+1):nsp))
      } else {
        rv$species_names <- rv$species_names[1:nsp]
      }

      rv$ne_min   <- expand_vec(rv$ne_min,   nsp, 10000)
      rv$ne_max   <- expand_vec(rv$ne_max,   nsp, 500000)
      rv$nea_min  <- expand_vec(rv$nea_min,  nsp, 0.01)
      rv$nea_max  <- expand_vec(rv$nea_max,  nsp, 10)
      rv$time_min <- expand_vec(rv$time_min, nsp, 10000)
      rv$time_max <- expand_vec(rv$time_max, nsp, 500000)
      rv$gen_time <- expand_vec(rv$gen_time, nsp, 1)
      rv$mu_dist  <- expand_vec(rv$mu_dist, nsp, "rnorm")
      rv$mu_p1    <- expand_vec(rv$mu_p1,   nsp, 1e-8)
      rv$mu_p2    <- expand_vec(rv$mu_p2,   nsp, 1.5e-9)

      # Expand gene_prior list
      old_gp <- rv$gene_prior
      if (length(old_gp) >= nsp) {
        rv$gene_prior <- old_gp[1:nsp]
      } else {
        for (i in (length(old_gp)+1):nsp) {
          old_gp[[i]] <- make_locfile(10, 10, 500)
        }
        rv$gene_prior <- old_gp
      }

      # Update species dropdown
      shiny::updateSelectInput(session, "hm_gene_species",
        choices = rv$species_names, selected = rv$species_names[1])

      shiny::showNotification(paste(nsp, "species configured."), type = "message")
    })

    # Species names UI
    output$hm_species_names_ui <- shiny::renderUI({
      nsp <- rv$nspecies
      inputs <- lapply(1:nsp, function(i) {
        shiny::textInput(paste0("hm_sp_name_", i),
          paste("Species", i, ":"),
          value = rv$species_names[i])
      })
      do.call(shiny::tagList, inputs)
    })

    shiny::observeEvent(input$hm_btn_update_names, {
      nsp <- rv$nspecies
      new_names <- sapply(1:nsp, function(i) {
        val <- input[[paste0("hm_sp_name_", i)]]
        if (is.null(val) || trimws(val) == "") paste0("sp", i) else trimws(val)
      })
      rv$species_names <- new_names
      shiny::updateSelectInput(session, "hm_gene_species",
        choices = new_names, selected = new_names[1])
      shiny::showNotification("Species names updated.", type = "message")
    })

    # =======================================================================
    # Ne PRIORS TABLE
    # =======================================================================

    output$hm_table_ne <- DT::renderDT({
      df <- data.frame(
        Species = rv$species_names,
        Min_Ne  = rv$ne_min,
        Max_Ne  = rv$ne_max,
        stringsAsFactors = FALSE
      )
      DT::datatable(df, editable = list(target = "cell", disable = list(columns = 0)),
        selection = "none", rownames = FALSE,
        options = list(dom = "t", paging = FALSE, ordering = FALSE))
    })

    shiny::observeEvent(input$hm_table_ne_cell_edit, {
      info <- input$hm_table_ne_cell_edit
      if (info$col == 1) rv$ne_min[info$row] <- as.numeric(info$value)
      if (info$col == 2) rv$ne_max[info$row] <- as.numeric(info$value)
    })

    # =======================================================================
    # NeA PRIORS TABLE
    # =======================================================================

    output$hm_table_nea <- DT::renderDT({
      df <- data.frame(
        Species       = rv$species_names,
        Min_NeA_ratio = rv$nea_min,
        Max_NeA_ratio = rv$nea_max,
        stringsAsFactors = FALSE
      )
      DT::datatable(df, editable = list(target = "cell", disable = list(columns = 0)),
        selection = "none", rownames = FALSE,
        options = list(dom = "t", paging = FALSE, ordering = FALSE))
    })

    shiny::observeEvent(input$hm_table_nea_cell_edit, {
      info <- input$hm_table_nea_cell_edit
      if (info$col == 1) rv$nea_min[info$row] <- as.numeric(info$value)
      if (info$col == 2) rv$nea_max[info$row] <- as.numeric(info$value)
    })

    # =======================================================================
    # TIME PRIORS TABLE
    # =======================================================================

    output$hm_table_time <- DT::renderDT({
      df <- data.frame(
        Species        = rv$species_names,
        Min_Time       = rv$time_min,
        Max_Time       = rv$time_max,
        Generation_Time = rv$gen_time,
        stringsAsFactors = FALSE
      )
      DT::datatable(df, editable = list(target = "cell", disable = list(columns = 0)),
        selection = "none", rownames = FALSE,
        options = list(dom = "t", paging = FALSE, ordering = FALSE))
    })

    shiny::observeEvent(input$hm_table_time_cell_edit, {
      info <- input$hm_table_time_cell_edit
      if (info$col == 1) rv$time_min[info$row]  <- as.numeric(info$value)
      if (info$col == 2) rv$time_max[info$row]  <- as.numeric(info$value)
      if (info$col == 3) rv$gen_time[info$row]   <- as.numeric(info$value)
    })

    # =======================================================================
    # GENE/LOCUS SETUP
    # =======================================================================

    # Update species dropdown when names change
    shiny::observe({
      shiny::updateSelectInput(session, "hm_gene_species",
        choices = rv$species_names, selected = rv$species_names[1])
    })

    # Load from phylip file for the selected species
    shiny::observeEvent(input$hm_btn_load_phylip, {
      sp_name <- input$hm_gene_species
      if (is.null(sp_name)) return()
      idx <- match(sp_name, rv$species_names)
      if (is.na(idx)) return()

      phy_path <- trimws(input$hm_phylip_path)
      if (is.null(phy_path) || phy_path == "") {
        shiny::showNotification("Please enter a path to a PHYLIP file.", type = "error")
        return()
      }

      # Expand ~ in path
      phy_path <- path.expand(phy_path)

      if (!file.exists(phy_path)) {
        shiny::showNotification(paste("File not found:", phy_path), type = "error", duration = 8)
        return()
      }

      tryCatch({
        loci <- read.phylip.loci(phy_path)
        inherit <- input$hm_gene_inherit
        if (is.null(inherit)) inherit <- 1

        lf <- locfile_from_phylip(loci)
        rv$gene_prior[[idx]] <- lf
        rv$phylip_paths[[sp_name]] <- phy_path

        shiny::showNotification(
          paste0("Loaded ", nrow(lf), " loci for ", sp_name,
                 " (", nrow(loci[[1]]), " samples, ",
                 min(vapply(loci, ncol, integer(1))), "-",
                 max(vapply(loci, ncol, integer(1))), " bp)"),
          type = "message", duration = 8)
      }, error = function(e) {
        shiny::showNotification(paste("Error reading phylip file:", e$message),
          type = "error", duration = 10)
      })
    })

    # Phylip load status display
    output$hm_phylip_status <- shiny::renderUI({
      sp_name <- input$hm_gene_species
      if (is.null(sp_name)) return(NULL)
      phy_path <- rv$phylip_paths[[sp_name]]
      if (is.null(phy_path)) {
        shiny::tags$div(class = "error-message",
          shiny::icon("info-circle"), " No phylip file loaded for this species. ",
          "Enter a path and click 'Load from Phylip', or use manual configuration below.")
      } else {
        idx <- match(sp_name, rv$species_names)
        nloci <- if (!is.null(rv$gene_prior[[idx]])) nrow(rv$gene_prior[[idx]]) else 0
        shiny::tags$div(class = "success-message",
          shiny::icon("check-circle"),
          sprintf(" Loaded: %s (%d loci)", basename(phy_path), nloci))
      }
    })

    # Manual: apply to selected species
    shiny::observeEvent(input$hm_btn_gene_apply, {
      sp_name <- input$hm_gene_species
      if (is.null(sp_name)) return()
      idx <- match(sp_name, rv$species_names)
      if (is.na(idx)) return()

      nloci   <- input$hm_gene_nloci
      nsamp   <- input$hm_gene_nsamp
      bp      <- input$hm_gene_bp
      inherit <- input$hm_gene_inherit

      if (is.null(nloci) || is.null(nsamp) || is.null(bp)) return()

      rv$gene_prior[[idx]] <- make_locfile(nloci, nsamp, bp)
      rv$phylip_paths[[sp_name]] <- NULL
      shiny::showNotification(paste("Locfile updated for", sp_name, "(manual)"), type = "message")
    })

    # Manual: apply to all species
    shiny::observeEvent(input$hm_btn_gene_apply_all, {
      nloci   <- input$hm_gene_nloci
      nsamp   <- input$hm_gene_nsamp
      bp      <- input$hm_gene_bp
      inherit <- input$hm_gene_inherit

      if (is.null(nloci) || is.null(nsamp) || is.null(bp)) return()

      for (i in seq_len(rv$nspecies)) {
        rv$gene_prior[[i]] <- make_locfile(nloci, nsamp, bp)
        rv$phylip_paths[[rv$species_names[i]]] <- NULL
      }
      shiny::showNotification("Locfile updated for all species (manual).", type = "message")
    })

    # Locfile preview title
    output$hm_locfile_title <- shiny::renderText({
      sp_name <- input$hm_gene_species
      if (is.null(sp_name)) return("Select a species")
      idx <- match(sp_name, rv$species_names)
      if (is.na(idx) || is.null(rv$gene_prior[[idx]])) return(paste(sp_name, "- no data"))
      nloci <- nrow(rv$gene_prior[[idx]])
      paste0(sp_name, " - ", nloci, " loci (showing first 10)")
    })

    # Locfile preview table
    output$hm_table_locfile_preview <- DT::renderDT({
      sp_name <- input$hm_gene_species
      if (is.null(sp_name)) return(NULL)
      idx <- match(sp_name, rv$species_names)
      if (is.na(idx) || is.null(rv$gene_prior[[idx]])) return(NULL)

      lf <- rv$gene_prior[[idx]]
      preview <- utils::head(lf, 10)
      DT::datatable(preview, rownames = FALSE, selection = "none",
        options = list(dom = "t", paging = FALSE, ordering = FALSE))
    })

    # =======================================================================
    # HIERARCHICAL PARAMS (read from inputs reactively)
    # =======================================================================

    shiny::observe({
      rv$coexp_min <- input$hm_coexp_min
      rv$coexp_max <- input$hm_coexp_max
    })

    shiny::observe({
      if (input$hm_zeta_type == "FREE") {
        rv$var_zeta <- "FREE"
      } else {
        rv$var_zeta <- input$hm_zeta_val
      }
    })

    shiny::observe({
      rv$th_val    <- input$hm_th
      rv$alpha_val <- input$hm_alpha
    })

    # Per-species mutation rate table
    output$hm_table_mu <- DT::renderDT({
      df <- data.frame(
        Species = rv$species_names,
        Distribution = rv$mu_dist,
        Param1 = rv$mu_p1,
        Param2 = rv$mu_p2,
        stringsAsFactors = FALSE
      )
      DT::datatable(df, editable = list(target = "cell", disable = list(columns = 0)),
        selection = "none", rownames = FALSE,
        options = list(dom = "t", paging = FALSE, ordering = FALSE))
    })

    shiny::observeEvent(input$hm_table_mu_cell_edit, {
      info <- input$hm_table_mu_cell_edit
      if (info$col == 1) rv$mu_dist[info$row] <- as.character(info$value)
      if (info$col == 2) rv$mu_p1[info$row]   <- as.numeric(info$value)
      if (info$col == 3) rv$mu_p2[info$row]   <- as.numeric(info$value)
    })

    # =======================================================================
    # STATUS BOXES
    # =======================================================================

    output$hm_status_species <- shinydashboard::renderValueBox({
      shinydashboard::valueBox(
        value = rv$nspecies, subtitle = "Species",
        icon = shiny::icon("paw"),
        color = if (rv$nspecies >= 2) "green" else "red")
    })

    output$hm_status_loci <- shinydashboard::renderValueBox({
      total_loci <- sum(sapply(rv$gene_prior, function(gp) {
        if (is.data.frame(gp)) nrow(gp) else 0
      }))
      shinydashboard::valueBox(
        value = total_loci, subtitle = "Total Loci",
        icon = shiny::icon("dna"),
        color = if (total_loci > 0) "green" else "red")
    })

    # =======================================================================
    # MODEL SUMMARY
    # =======================================================================

    output$hm_txt_summary <- shiny::renderPrint({
      cat("hModel Configuration:\n")
      cat("=====================\n\n")
      cat(sprintf("Species (%d):     %s\n", rv$nspecies,
                  paste(rv$species_names, collapse = ", ")))
      cat(sprintf("var.zeta:         %s\n", as.character(rv$var_zeta)))
      cat(sprintf("coexp.prior:      [%s, %s]\n", rv$coexp_min, rv$coexp_max))
      cat(sprintf("th:               %s\n", rv$th_val))
      cat(sprintf("alpha:            %s\n", rv$alpha_val))
      cat("\nMutation rate distributions:\n")
      for (i in seq_len(rv$nspecies)) {
        cat(sprintf("  %s: %s(%s, %s)\n", rv$species_names[i],
                    rv$mu_dist[i], rv$mu_p1[i], rv$mu_p2[i]))
      }
      cat("\nNe priors:\n")
      for (i in seq_len(rv$nspecies)) {
        cat(sprintf("  %s: [%s, %s]\n", rv$species_names[i], rv$ne_min[i], rv$ne_max[i]))
      }
      cat("\nAncestral Ne ratio priors:\n")
      for (i in seq_len(rv$nspecies)) {
        cat(sprintf("  %s: [%s, %s]\n", rv$species_names[i], rv$nea_min[i], rv$nea_max[i]))
      }
      cat("\nTime priors:\n")
      for (i in seq_len(rv$nspecies)) {
        cat(sprintf("  %s: [%s, %s], gen=%s\n", rv$species_names[i],
                    rv$time_min[i], rv$time_max[i], rv$gen_time[i]))
      }
      cat("\nLoci per species:\n")
      for (i in seq_len(rv$nspecies)) {
        nloci <- if (is.data.frame(rv$gene_prior[[i]])) nrow(rv$gene_prior[[i]]) else 0
        cat(sprintf("  %s: %d loci\n", rv$species_names[i], nloci))
      }
    })

    # Usage code
    output$hm_txt_usage <- shiny::renderText({
      "hm <- h.menu.gui()\nsim.coexp.ngs(hm, nsims = 1000)\n\n# Reload into GUI:\nhm2 <- h.menu.gui(hm)"
    })

    # =======================================================================
    # BUILD / DOWNLOAD / CANCEL
    # =======================================================================

    shiny::observeEvent(input$hm_btn_build, {
      tryCatch({
        model <- assemble_hmodel(rv)
        shiny::stopApp(returnValue = model)
      }, error = function(e) {
        shiny::showNotification(paste("Build error:", e$message),
          type = "error", duration = 10)
      })
    })

    output$hm_btn_download <- shiny::downloadHandler(
      filename = function() { "hModel.rds" },
      content = function(file) {
        tryCatch({
          model <- assemble_hmodel(rv)
          saveRDS(model, file)
        }, error = function(e) {
          shiny::showNotification(paste("Error:", e$message), type = "error")
        })
      }
    )

    shiny::observeEvent(input$hm_btn_cancel, {
      shiny::stopApp(returnValue = NULL)
    })

  } # end server

  # ===========================================================================
  # LAUNCH APP
  # ===========================================================================
  shiny::runApp(shiny::shinyApp(ui = ui, server = server))
}
