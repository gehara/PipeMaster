#' PipeMaster GUI - Shiny Application Server
#' @description Server logic for PipeMaster graphical interface

library(shiny)
library(ape)
library(phytools)

# Define server logic
server <- function(input, output, session) {
  
  # ============================================================================
  # REACTIVE VALUES
  # ============================================================================
  
  rv <- reactiveValues(
    model = NULL,
    tree = NULL,
    n_pops = 0,
    nodes = NULL,
    current_ne = NULL,
    ancestral_ne = NULL,
    migration = NULL,
    loci = NULL,
    samples = NULL,
    conditions = list(
      size = NULL,
      mig = NULL,
      time = NULL
    ),
    recent_changes = list(),
    validation_errors = list()
  )
  
  # ============================================================================
  # HELPER FUNCTIONS
  # ============================================================================
  
  # Add change to recent changes log
  add_change <- function(message) {
    timestamp <- format(Sys.time(), "%H:%M:%S")
    rv$recent_changes <- c(
      list(list(time = timestamp, msg = message)),
      rv$recent_changes
    )
    # Keep only last 10 changes
    if(length(rv$recent_changes) > 10) {
      rv$recent_changes <- rv$recent_changes[1:10]
    }
  }
  
  # Validate tree string
  validate_tree_string <- function(tree_str) {
    if(tree_str == "1") {
      return(list(valid = TRUE, message = "Single population model"))
    }
    
    n_open <- str_count(tree_str, fixed("("))
    n_close <- str_count(tree_str, fixed(")"))
    n_comma <- str_count(tree_str, fixed(","))
    
    if(n_open != n_close) {
      return(list(
        valid = FALSE,
        message = sprintf("Parentheses mismatch: %d '(' vs %d ')'", n_open, n_close)
      ))
    }
    
    if(n_open != n_comma) {
      return(list(
        valid = FALSE,
        message = sprintf("Not bifurcating: %d nodes vs %d commas", n_open, n_comma)
      ))
    }
    
    # Try to parse with ape
    tree_obj <- tryCatch(
      read.tree(text = tree_str),
      error = function(e) NULL
    )
    
    if(is.null(tree_obj)) {
      return(list(valid = FALSE, message = "Invalid Newick format"))
    }
    
    return(list(valid = TRUE, message = "Valid tree!", tree = tree_obj))
  }
  
  # Initialize default Ne table
  initialize_ne_table <- function(n_pops) {
    data.frame(
      Population = paste0("pop", 1:n_pops),
      Parameter = paste0("Ne0.pop", 1:n_pops),
      Min = rep(100000, n_pops),
      Max = rep(500000, n_pops),
      stringsAsFactors = FALSE
    )
  }
  
  # ============================================================================
  # GETTING STARTED TAB
  # ============================================================================
  
  # Start new model
  observeEvent(input$btn_new_model, {
    showModal(modalDialog(
      title = "Start New Model",
      "This will clear any existing model. Continue?",
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_new_model", "Yes, Start New", class = "btn-primary")
      )
    ))
  })
  
  observeEvent(input$confirm_new_model, {
    rv$model <- NULL
    rv$tree <- NULL
    rv$n_pops <- 0
    rv$nodes <- NULL
    rv$current_ne <- NULL
    rv$ancestral_ne <- NULL
    rv$migration <- NULL
    rv$loci <- NULL
    rv$samples <- NULL
    rv$conditions <- list(size = NULL, mig = NULL, time = NULL)
    rv$recent_changes <- list()
    
    add_change("New model started")
    removeModal()
    
    showNotification("New model initialized!", type = "message")
  })
  
  # Load existing model
  observeEvent(input$file_load_model, {
    req(input$file_load_model)
    
    tryCatch({
      ext <- tools::file_ext(input$file_load_model$name)
      
      if(ext == "rds") {
        model <- readRDS(input$file_load_model$datapath)
      } else if(ext == "RData") {
        load(input$file_load_model$datapath)
        # Assume model is named 'model'
        if(!exists("model")) {
          stop("No 'model' object found in RData file")
        }
      }
      
      # Load model data into reactive values
      rv$model <- model
      # ... (parse model object and populate rv)
      
      add_change(sprintf("Loaded model from %s", input$file_load_model$name))
      showNotification("Model loaded successfully!", type = "message")
      
    }, error = function(e) {
      showNotification(
        paste("Error loading model:", e$message),
        type = "error",
        duration = 10
      )
    })
  })
  
  # Recent changes display
  output$recent_changes <- renderUI({
    req(length(rv$recent_changes) > 0)
    
    changes_html <- lapply(rv$recent_changes, function(change) {
      tags$div(
        style = "margin-bottom: 5px;",
        tags$span(
          style = "color: #999; font-size: 11px;",
          change$time
        ),
        " - ",
        change$msg
      )
    })
    
    tagList(changes_html)
  })
  
  # ============================================================================
  # POPULATION STRUCTURE TAB
  # ============================================================================
  
  # Validate tree button
  observeEvent(input$btn_validate_tree, {
    req(input$txt_tree)
    
    result <- validate_tree_string(input$txt_tree)
    
    if(result$valid) {
      rv$tree <- result$tree
      
      if(input$txt_tree == "1") {
        rv$n_pops <- 1
        rv$nodes <- NULL
      } else {
        rv$n_pops <- length(rv$tree$tip.label)
        # Extract nodes
        rv$nodes <- data.frame(
          Node = 1:(rv$tree$Nnode),
          Tips = sapply(1:(rv$tree$Nnode), function(i) {
            paste(rv$tree$tip.label[which(rv$tree$edge[,1] == (i + rv$n_pops))], collapse = ", ")
          }),
          stringsAsFactors = FALSE
        )
      }
      
      # Initialize Ne table
      rv$current_ne <- initialize_ne_table(rv$n_pops)
      
      add_change(sprintf("Tree validated: %d populations", rv$n_pops))
      
      showNotification(
        paste("Tree validated!", result$message),
        type = "message"
      )
    } else {
      showNotification(
        paste("Tree validation failed:", result$message),
        type = "error",
        duration = 10
      )
    }
  })
  
  # Tree validation message
  output$tree_validation_message <- renderUI({
    req(input$txt_tree)
    
    result <- validate_tree_string(input$txt_tree)
    
    if(result$valid) {
      tags$div(
        class = "success-message",
        icon("check-circle"),
        " Valid tree! ",
        sprintf("Detected %d populations.", rv$n_pops)
      )
    } else {
      tags$div(
        class = "error-message",
        icon("exclamation-triangle"),
        " Error: ",
        result$message
      )
    }
  })
  
  # Plot tree
  output$plot_tree <- renderPlot({
    req(rv$tree)
    
    if(!is.null(rv$tree)) {
      plot(rv$tree, 
           main = "Population Tree",
           cex = 1.2,
           font = 2,
           edge.width = 2,
           edge.color = "#3c8dbc")
      axisPhylo()
    }
  })
  
  # Nodes table
  output$table_nodes <- renderDT({
    req(rv$nodes)
    
    datatable(
      rv$nodes,
      options = list(
        pageLength = 10,
        dom = 't'
      ),
      rownames = FALSE
    )
  })
  
  # Example tree button
  observeEvent(input$btn_example_tree, {
    updateTextInput(session, "txt_tree", value = "((A,B),(C,D));")
  })
  
  # ============================================================================
  # DEMOGRAPHY TAB
  # ============================================================================
  
  # Current Ne table
  output$table_current_ne <- renderDT({
    req(rv$current_ne)
    
    datatable(
      rv$current_ne,
      selection = 'single',
      options = list(
        pageLength = 10,
        dom = 'tp'
      ),
      rownames = FALSE,
      editable = TRUE
    )
  })
  
  # Edit current Ne
  observeEvent(input$btn_edit_current_ne, {
    req(input$table_current_ne_rows_selected)
    
    selected_row <- input$table_current_ne_rows_selected
    
    showModal(modalDialog(
      title = sprintf("Edit %s", rv$current_ne$Parameter[selected_row]),
      
      if(input$select_ne_dist == "uniform") {
        tagList(
          numericInput("modal_ne_min",
                       "Minimum:",
                       value = rv$current_ne$Min[selected_row],
                       min = 0),
          numericInput("modal_ne_max",
                       "Maximum:",
                       value = rv$current_ne$Max[selected_row],
                       min = 0)
        )
      } else {
        tagList(
          numericInput("modal_ne_mean",
                       "Mean:",
                       value = (rv$current_ne$Min[selected_row] + rv$current_ne$Max[selected_row]) / 2,
                       min = 0),
          numericInput("modal_ne_sd",
                       "Standard Deviation:",
                       value = (rv$current_ne$Max[selected_row] - rv$current_ne$Min[selected_row]) / 4,
                       min = 0)
        )
      },
      
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_edit_ne", "Save", class = "btn-primary")
      )
    ))
  })
  
  observeEvent(input$confirm_edit_ne, {
    selected_row <- input$table_current_ne_rows_selected
    
    if(input$select_ne_dist == "uniform") {
      rv$current_ne$Min[selected_row] <- input$modal_ne_min
      rv$current_ne$Max[selected_row] <- input$modal_ne_max
    } else {
      # Store as mean ± 2*SD for uniform approximation
      rv$current_ne$Min[selected_row] <- input$modal_ne_mean - 2 * input$modal_ne_sd
      rv$current_ne$Max[selected_row] <- input$modal_ne_mean + 2 * input$modal_ne_sd
    }
    
    add_change(sprintf("Updated %s", rv$current_ne$Parameter[selected_row]))
    removeModal()
  })
  
  # ============================================================================
  # GENE SETUP TAB
  # ============================================================================
  
  # Initialize loci table when number of loci changes
  observeEvent(input$num_loci, {
    req(input$num_loci > 0)
    
    rv$loci <- data.frame(
      Locus = paste0("locus", 1:input$num_loci),
      Length = rep(1000, input$num_loci),
      Inheritance = rep(1, input$num_loci),
      Mut_Rate_Min = rep(5e-9, input$num_loci),
      Mut_Rate_Max = rep(1.5e-8, input$num_loci),
      stringsAsFactors = FALSE
    )
    
    # Initialize samples table
    if(!is.null(rv$n_pops) && rv$n_pops > 0) {
      sample_matrix <- matrix(10, nrow = input$num_loci, ncol = rv$n_pops)
      colnames(sample_matrix) <- paste0("Pop", 1:rv$n_pops)
      
      rv$samples <- data.frame(
        Locus = paste0("locus", 1:input$num_loci),
        sample_matrix,
        stringsAsFactors = FALSE
      )
    }
  })
  
  # Loci table
  output$table_loci <- renderDT({
    req(rv$loci)
    
    datatable(
      rv$loci,
      selection = 'single',
      options = list(
        pageLength = 10,
        dom = 'tp',
        scrollX = TRUE
      ),
      rownames = FALSE,
      editable = TRUE
    )
  })
  
  # Samples table
  output$table_samples <- renderDT({
    req(rv$samples)
    
    datatable(
      rv$samples,
      options = list(
        pageLength = 10,
        dom = 'tp',
        scrollX = TRUE
      ),
      rownames = FALSE,
      editable = TRUE
    )
  })
  
  # Auto-fill loci
  observeEvent(input$btn_auto_fill_loci, {
    req(rv$loci)
    
    showModal(modalDialog(
      title = "Auto-fill Locus Parameters",
      
      numericInput("modal_locus_length",
                   "Sequence Length (bp):",
                   value = 1000,
                   min = 1),
      numericInput("modal_mut_min",
                   "Min Mutation Rate:",
                   value = 5e-9,
                   min = 0),
      numericInput("modal_mut_max",
                   "Max Mutation Rate:",
                   value = 1.5e-8,
                   min = 0),
      
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_auto_fill", "Apply to All", class = "btn-primary")
      )
    ))
  })
  
  observeEvent(input$confirm_auto_fill, {
    rv$loci$Length <- input$modal_locus_length
    rv$loci$Mut_Rate_Min <- input$modal_mut_min
    rv$loci$Mut_Rate_Max <- input$modal_mut_max
    
    add_change("Auto-filled locus parameters")
    removeModal()
  })
  
  # ============================================================================
  # MODEL SUMMARY TAB
  # ============================================================================
  
  output$txt_model_summary <- renderPrint({
    cat("Model Configuration:\n")
    cat("===================\n\n")
    
    cat(sprintf("Tree: %s\n", ifelse(is.null(rv$tree), "Not set", input$txt_tree)))
    cat(sprintf("Populations: %d\n", rv$n_pops))
    
    if(!is.null(rv$nodes)) {
      cat(sprintf("Nodes: %d\n", nrow(rv$nodes)))
    }
    
    cat(sprintf("\nNe Distribution: %s\n", input$select_ne_dist))
    cat(sprintf("Migration: %s\n", ifelse(input$check_migration, "Enabled", "Disabled")))
    cat(sprintf("Ancestral Ne: %s\n", ifelse(input$check_ancestral_ne, "Enabled", "Disabled")))
    
    cat(sprintf("\nData Type: %s\n", input$select_data_type))
    cat(sprintf("Number of Loci: %d\n", input$num_loci))
    cat(sprintf("Mutation Rate Distribution: %s\n", input$select_mut_dist))
    
    cat("\n")
  })
  
  # ============================================================================
  # STATUS BOXES
  # ============================================================================
  
  output$status_populations <- renderValueBox({
    valueBox(
      value = rv$n_pops,
      subtitle = "Populations",
      icon = icon("sitemap"),
      color = if(rv$n_pops > 0) "green" else "red"
    )
  })
  
  output$status_nodes <- renderValueBox({
    n_nodes <- if(is.null(rv$nodes)) 0 else nrow(rv$nodes)
    valueBox(
      value = n_nodes,
      subtitle = "Divergence Nodes",
      icon = icon("code-branch"),
      color = if(n_nodes > 0) "blue" else "yellow"
    )
  })
  
  output$status_loci <- renderValueBox({
    valueBox(
      value = input$num_loci,
      subtitle = "Loci",
      icon = icon("dna"),
      color = if(input$num_loci > 0) "green" else "red"
    )
  })
  
  # ============================================================================
  # EXPORT TAB
  # ============================================================================
  
  # Simulation code display
  output$txt_simulation_code <- renderText({
    sprintf('# Load your model
model <- readRDS("%s.rds")

# Run simulations
library(PipeMaster)

sim.ms.sumstat(
  model = model,
  nsim.blocks = 10,
  sim.block.size = 1000,
  output.name = "%s",
  perpop.SS = TRUE,
  overall.SS = TRUE
)
', input$txt_model_name, input$txt_model_name)
  })
  
  # Download model
  output$btn_download_model <- downloadHandler(
    filename = function() {
      paste0(input$txt_model_name, ".", input$select_export_format)
    },
    content = function(file) {
      # Build model object from reactive values
      model <- list(
        tree = input$txt_tree,
        n_pops = rv$n_pops,
        flags = list(
          n = rv$current_ne,
          ej = rv$nodes,
          en = rv$ancestral_ne,
          m = rv$migration
        ),
        loci = rv$loci,
        I = rv$samples,
        conds = rv$conditions
      )
      
      if(input$select_export_format == "rds") {
        saveRDS(model, file)
      } else {
        save(model, file = file)
      }
    }
  )
  
  # Validate model
  observeEvent(input$btn_validate_model, {
    errors <- character()
    warnings <- character()
    
    # Check required fields
    if(is.null(rv$tree) || rv$n_pops == 0) {
      errors <- c(errors, "Population structure not defined")
    }
    
    if(is.null(rv$current_ne)) {
      errors <- c(errors, "Ne priors not set")
    }
    
    if(is.null(rv$loci) || nrow(rv$loci) == 0) {
      errors <- c(errors, "No loci configured")
    }
    
    if(is.null(rv$samples)) {
      errors <- c(errors, "Sample sizes not set")
    }
    
    # Check for warnings
    if(input$check_migration && is.null(rv$migration)) {
      warnings <- c(warnings, "Migration enabled but not configured")
    }
    
    rv$validation_errors <- list(errors = errors, warnings = warnings)
  })
  
  # Validation results
  output$model_validation_results <- renderUI({
    req(length(rv$validation_errors) > 0)
    
    errors <- rv$validation_errors$errors
    warnings <- rv$validation_errors$warnings
    
    tagList(
      if(length(errors) > 0) {
        tags$div(
          class = "error-message",
          h5(icon("exclamation-circle"), " Errors:"),
          tags$ul(
            lapply(errors, function(e) tags$li(e))
          )
        )
      },
      
      if(length(warnings) > 0) {
        tags$div(
          style = "color: #f39c12; padding: 10px; background-color: #fcf8e3; border-radius: 5px; margin-top: 10px;",
          h5(icon("exclamation-triangle"), " Warnings:"),
          tags$ul(
            lapply(warnings, function(w) tags$li(w))
          )
        )
      },
      
      if(length(errors) == 0 && length(warnings) == 0) {
        tags$div(
          class = "success-message",
          h5(icon("check-circle"), " Model is valid!"),
          p("All required parameters are set. You can download and use this model.")
        )
      }
    )
  })
}
