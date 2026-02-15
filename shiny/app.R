#' PipeMaster GUI - Shiny Application
#' @description Modern graphical interface for building coalescent models
#' @author Created by Claude AI for PipeMaster package
#' 
#' To run this app:
#'   1. Install required packages (see below)
#'   2. Set working directory to folder containing app_ui.R and app_server.R
#'   3. Run: shiny::runApp()

# ============================================================================
# PACKAGE DEPENDENCIES
# ============================================================================

required_packages <- c(
  "shiny",           # Core Shiny framework
  "shinydashboard",  # Dashboard layout
  "shinyjs",         # JavaScript operations
  "shinyWidgets",    # Additional widgets
  "DT",              # Interactive tables
  "plotly",          # Interactive plots
  "ape",             # Phylogenetic trees
  "phytools"         # Phylogenetic visualization
)

# Check and install missing packages
missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]

if(length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages)
}

# Load packages
lapply(required_packages, library, character.only = TRUE)

# ============================================================================
# LOAD UI AND SERVER
# ============================================================================

source("app_ui.R")
source("app_server.R")

# ============================================================================
# LAUNCH APPLICATION
# ============================================================================

shinyApp(ui = ui, server = server)
