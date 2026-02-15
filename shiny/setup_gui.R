#' PipeMaster GUI - Setup Script
#' @description Automated installation and setup for PipeMaster GUI
#' @author Created by Claude AI
#' 
#' Run this script to automatically:
#'   1. Install required packages
#'   2. Download GUI files
#'   3. Launch the application

cat("\n")
cat("============================================================\n")
cat("   PipeMaster GUI - Automated Setup\n")
cat("============================================================\n\n")

# ============================================================================
# STEP 1: Check R Version
# ============================================================================

cat("Step 1: Checking R version...\n")

r_version <- getRversion()
required_version <- "4.0.0"

if(r_version < required_version) {
  stop(sprintf(
    "R version %s or higher is required. You have %s. Please upgrade R.",
    required_version,
    r_version
  ))
}

cat(sprintf("✓ R version %s detected (OK)\n\n", r_version))

# ============================================================================
# STEP 2: Install Required Packages
# ============================================================================

cat("Step 2: Installing required packages...\n")

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

# Check which packages are missing
missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]

if(length(missing_packages) > 0) {
  cat(sprintf("Installing %d missing packages...\n", length(missing_packages)))
  cat(paste("  -", missing_packages, collapse = "\n"), "\n")
  
  install.packages(missing_packages, quiet = TRUE)
  
  cat("✓ Packages installed successfully\n\n")
} else {
  cat("✓ All required packages already installed\n\n")
}

# Verify installation
failed_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]

if(length(failed_packages) > 0) {
  stop(sprintf(
    "Failed to install: %s\nPlease install manually with: install.packages(c(%s))",
    paste(failed_packages, collapse = ", "),
    paste(sprintf('"%s"', failed_packages), collapse = ", ")
  ))
}

# ============================================================================
# STEP 3: Check for PipeMaster Package
# ============================================================================

cat("Step 3: Checking for PipeMaster package...\n")

if(!"PipeMaster" %in% installed.packages()[,"Package"]) {
  cat("⚠ PipeMaster package not found\n")
  cat("Installing from GitHub...\n")
  
  if(!"devtools" %in% installed.packages()[,"Package"]) {
    install.packages("devtools", quiet = TRUE)
  }
  
  devtools::install_github("gehara/PipeMaster", quiet = TRUE)
  
  if("PipeMaster" %in% installed.packages()[,"Package"]) {
    cat("✓ PipeMaster installed successfully\n\n")
  } else {
    cat("⚠ Could not install PipeMaster automatically\n")
    cat("Please install manually:\n")
    cat('  devtools::install_github("gehara/PipeMaster")\n\n')
  }
} else {
  cat("✓ PipeMaster package found\n\n")
}

# ============================================================================
# STEP 4: Setup GUI Files
# ============================================================================

cat("Step 4: Setting up GUI files...\n")

# Create GUI directory
gui_dir <- file.path(getwd(), "PipeMaster_GUI")

if(!dir.exists(gui_dir)) {
  dir.create(gui_dir)
  cat(sprintf("✓ Created directory: %s\n", gui_dir))
} else {
  cat(sprintf("✓ Directory exists: %s\n", gui_dir))
}

# Check if GUI files exist
gui_files <- c("app.R", "app_ui.R", "app_server.R")
files_exist <- file.exists(file.path(gui_dir, gui_files))

if(all(files_exist)) {
  cat("✓ All GUI files present\n\n")
} else {
  cat("⚠ Some GUI files missing:\n")
  missing_files <- gui_files[!files_exist]
  cat(paste("  -", missing_files, collapse = "\n"), "\n")
  cat("\nPlease ensure these files are in the PipeMaster_GUI directory:\n")
  cat("  - app.R\n")
  cat("  - app_ui.R\n")
  cat("  - app_server.R\n\n")
  
  response <- readline("Copy files to this directory now? (Y/N): ")
  
  if(toupper(response) %in% c("Y", "YES")) {
    cat("\nWaiting for file copy...\n")
    cat("Press Enter when files are in place: ")
    readline()
  } else {
    stop("Setup incomplete. Please add GUI files and run setup again.")
  }
}

# ============================================================================
# STEP 5: Test Installation
# ============================================================================

cat("Step 5: Testing installation...\n")

test_passed <- TRUE

# Test package loading
for(pkg in required_packages) {
  test_result <- tryCatch({
    library(pkg, character.only = TRUE)
    TRUE
  }, error = function(e) {
    cat(sprintf("✗ Error loading %s: %s\n", pkg, e$message))
    FALSE
  })
  
  if(!test_result) {
    test_passed <- FALSE
  }
}

if(test_passed) {
  cat("✓ All packages load successfully\n\n")
} else {
  stop("Package loading test failed. Please check error messages above.")
}

# ============================================================================
# STEP 6: Launch GUI
# ============================================================================

cat("============================================================\n")
cat("   Setup Complete!\n")
cat("============================================================\n\n")

cat("PipeMaster GUI is ready to use.\n\n")

cat("To launch the GUI:\n")
cat(sprintf('  1. setwd("%s")\n', gui_dir))
cat('  2. shiny::runApp()\n\n')

cat("Or run now automatically:\n\n")

response <- readline("Launch PipeMaster GUI now? (Y/N): ")

if(toupper(response) %in% c("Y", "YES")) {
  cat("\nLaunching GUI...\n\n")
  
  setwd(gui_dir)
  shiny::runApp(launch.browser = TRUE)
  
} else {
  cat("\nYou can launch the GUI later with:\n")
  cat(sprintf('  setwd("%s")\n', gui_dir))
  cat('  shiny::runApp()\n\n')
}

cat("============================================================\n")
cat("   For help, see GUI_USER_GUIDE.md\n")
cat("============================================================\n\n")
