#' PipeMaster GUI - Shiny Application UI
#' @description Modern graphical interface for PipeMaster model building
#' @author Created by Claude AI

library(shiny)
library(shinydashboard)
library(shinyjs)
library(DT)
library(plotly)
library(shinyWidgets)

# Define UI
ui <- dashboardPage(
  
  # ============================================================================
  # HEADER
  # ============================================================================
  dashboardHeader(
    title = "PipeMaster - ABC Model Builder",
    titleWidth = 350,
    tags$li(class = "dropdown",
            tags$a(href = "https://github.com/gehara/PipeMaster",
                   target = "_blank",
                   icon("github"),
                   "Documentation"
            )
    )
  ),
  
  # ============================================================================
  # SIDEBAR
  # ============================================================================
  dashboardSidebar(
    width = 300,
    sidebarMenu(
      id = "sidebar",
      menuItem("Getting Started", tabName = "start", icon = icon("play-circle")),
      menuItem("Population Structure", tabName = "structure", icon = icon("sitemap")),
      menuItem("Demography (Ne)", tabName = "demography", icon = icon("chart-line")),
      menuItem("Migration", tabName = "migration", icon = icon("exchange-alt")),
      menuItem("Time Priors", tabName = "time", icon = icon("clock")),
      menuItem("Conditions", tabName = "conditions", icon = icon("filter")),
      menuItem("Gene Setup", tabName = "genes", icon = icon("dna")),
      menuItem("Model Summary", tabName = "summary", icon = icon("eye")),
      menuItem("Export Model", tabName = "export", icon = icon("download"))
    ),
    
    hr(),
    
    # Model status indicator
    box(
      title = "Model Status",
      width = 12,
      status = "info",
      solidHeader = FALSE,
      collapsible = FALSE,
      
      valueBoxOutput("status_populations", width = 12),
      valueBoxOutput("status_nodes", width = 12),
      valueBoxOutput("status_loci", width = 12)
    )
  ),
  
  # ============================================================================
  # BODY
  # ============================================================================
  dashboardBody(
    
    # Custom CSS
    tags$head(
      tags$style(HTML("
        .main-header .logo { font-weight: bold; }
        .content-wrapper, .right-side { background-color: #f4f6f9; }
        .box { border-radius: 5px; }
        .small-box { border-radius: 5px; }
        .btn-primary { background-color: #3c8dbc; }
        .tab-content { padding: 20px; }
        .dataTable { font-size: 12px; }
        .progress-bar { background-color: #3c8dbc; }
        .info-box { min-height: 80px; }
        .error-message { 
          color: #dd4b39; 
          font-weight: bold; 
          padding: 10px;
          background-color: #f2dede;
          border-radius: 5px;
        }
        .success-message {
          color: #00a65a;
          font-weight: bold;
          padding: 10px;
          background-color: #d4edda;
          border-radius: 5px;
        }
        .tree-display {
          font-family: 'Courier New', monospace;
          background-color: #f5f5f5;
          padding: 15px;
          border-radius: 5px;
          border: 1px solid #ddd;
        }
      "))
    ),
    
    # Enable shinyjs
    useShinyjs(),
    
    # ============================================================================
    # TAB ITEMS
    # ============================================================================
    tabItems(
      
      # --------------------------------------------------------------------------
      # GETTING STARTED TAB
      # --------------------------------------------------------------------------
      tabItem(
        tabName = "start",
        
        fluidRow(
          box(
            title = "Welcome to PipeMaster GUI",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            
            h3("Build your coalescent model visually"),
            p("This graphical interface makes it easy to set up complex population genetic models for ABC analysis."),
            
            hr(),
            
            h4("Quick Start Options:"),
            
            fluidRow(
              column(6,
                     actionButton("btn_new_model", 
                                  "Start New Model",
                                  icon = icon("plus-circle"),
                                  class = "btn-primary btn-lg btn-block",
                                  style = "margin-bottom: 10px;")
              ),
              column(6,
                     fileInput("file_load_model",
                               "Load Existing Model",
                               accept = c(".rds", ".RData"),
                               buttonLabel = "Browse...",
                               placeholder = "Choose model file")
              )
            ),
            
            hr(),
            
            h4("Model Building Workflow:"),
            tags$ol(
              tags$li(strong("Population Structure:"), " Define your population tree topology"),
              tags$li(strong("Demography:"), " Set effective population size priors"),
              tags$li(strong("Migration:"), " Configure gene flow between populations (optional)"),
              tags$li(strong("Time Priors:"), " Set temporal parameters"),
              tags$li(strong("Conditions:"), " Add parameter constraints"),
              tags$li(strong("Gene Setup:"), " Configure loci and mutation rates"),
              tags$li(strong("Export:"), " Save and run simulations")
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Recent Changes",
            width = 6,
            status = "info",
            solidHeader = TRUE,
            
            uiOutput("recent_changes")
          ),
          
          box(
            title = "Tips & Tricks",
            width = 6,
            status = "warning",
            solidHeader = TRUE,
            
            tags$ul(
              tags$li("Use the sidebar to navigate between sections"),
              tags$li("Green checkmarks indicate completed sections"),
              tags$li("Red warnings show what still needs attention"),
              tags$li("Hover over info icons for helpful tooltips"),
              tags$li("Your progress is automatically saved")
            )
          )
        )
      ),
      
      # --------------------------------------------------------------------------
      # POPULATION STRUCTURE TAB
      # --------------------------------------------------------------------------
      tabItem(
        tabName = "structure",
        
        fluidRow(
          box(
            title = "Population Structure",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            
            p("Define the topology of your population tree in Newick format, or specify a single population."),
            
            fluidRow(
              column(8,
                     textInput("txt_tree",
                               "Tree Topology (Newick format):",
                               value = "",
                               placeholder = "e.g., ((A,B),(C,D)); or 1 for single population")
              ),
              column(4,
                     br(),
                     actionButton("btn_validate_tree",
                                  "Validate Tree",
                                  icon = icon("check-circle"),
                                  class = "btn-success")
              )
            ),
            
            uiOutput("tree_validation_message"),
            
            hr(),
            
            h4("Tree Visualization:"),
            plotOutput("plot_tree", height = "400px")
          )
        ),
        
        fluidRow(
          box(
            title = "Detected Nodes",
            width = 6,
            status = "info",
            solidHeader = TRUE,
            
            DTOutput("table_nodes")
          ),
          
          box(
            title = "Tree Examples",
            width = 6,
            status = "warning",
            solidHeader = TRUE,
            
            h5("Common Tree Topologies:"),
            tags$ul(
              tags$li(strong("Single population:"), code("1")),
              tags$li(strong("Two populations:"), code("(A,B);")),
              tags$li(strong("Three populations (bifurcating):"), code("((A,B),C);")),
              tags$li(strong("Four populations:"), code("((A,B),(C,D));")),
              tags$li(strong("With branch lengths:"), code("((A:0.1,B:0.2):0.3,C:0.4);"))
            ),
            
            actionButton("btn_example_tree",
                         "Load Example Tree",
                         icon = icon("magic"),
                         class = "btn-info")
          )
        )
      ),
      
      # --------------------------------------------------------------------------
      # DEMOGRAPHY TAB
      # --------------------------------------------------------------------------
      tabItem(
        tabName = "demography",
        
        fluidRow(
          box(
            title = "Effective Population Size (Ne) Priors",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            
            fluidRow(
              column(6,
                     selectInput("select_ne_dist",
                                 "Prior Distribution:",
                                 choices = c("Uniform" = "uniform",
                                             "Normal" = "normal"),
                                 selected = "uniform")
              ),
              column(6,
                     checkboxInput("check_ancestral_ne",
                                   "Include ancestral Ne changes",
                                   value = FALSE)
              )
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Current Ne Priors",
            width = 12,
            status = "info",
            solidHeader = TRUE,
            
            DTOutput("table_current_ne"),
            
            br(),
            
            actionButton("btn_edit_current_ne",
                         "Edit Selected",
                         icon = icon("edit"),
                         class = "btn-primary")
          )
        ),
        
        conditionalPanel(
          condition = "input.check_ancestral_ne == true",
          
          fluidRow(
            box(
              title = "Ancestral Ne Priors",
              width = 12,
              status = "info",
              solidHeader = TRUE,
              
              DTOutput("table_ancestral_ne"),
              
              br(),
              
              fluidRow(
                column(6,
                       actionButton("btn_add_ancestral_ne",
                                    "Add Ancestral Ne",
                                    icon = icon("plus"),
                                    class = "btn-success")
                ),
                column(6,
                       actionButton("btn_edit_ancestral_ne",
                                    "Edit Selected",
                                    icon = icon("edit"),
                                    class = "btn-primary")
                )
              )
            )
          )
        )
      ),
      
      # --------------------------------------------------------------------------
      # MIGRATION TAB
      # --------------------------------------------------------------------------
      tabItem(
        tabName = "migration",
        
        fluidRow(
          box(
            title = "Gene Flow Configuration",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            
            fluidRow(
              column(6,
                     checkboxInput("check_migration",
                                   "Enable migration between populations",
                                   value = FALSE)
              ),
              column(6,
                     conditionalPanel(
                       condition = "input.check_migration == true",
                       selectInput("select_mig_dist",
                                   "Migration Prior Distribution:",
                                   choices = c("Uniform" = "uniform",
                                               "Normal" = "normal"),
                                   selected = "uniform")
                     )
              )
            )
          )
        ),
        
        conditionalPanel(
          condition = "input.check_migration == true",
          
          fluidRow(
            box(
              title = "Current Migration Rates (4Nm)",
              width = 12,
              status = "info",
              solidHeader = TRUE,
              
              DTOutput("table_migration"),
              
              br(),
              
              actionButton("btn_edit_migration",
                           "Edit Selected",
                           icon = icon("edit"),
                           class = "btn-primary")
            )
          ),
          
          fluidRow(
            box(
              title = "Migration Matrix Visualization",
              width = 12,
              status = "info",
              solidHeader = TRUE,
              
              plotlyOutput("plot_migration_matrix", height = "400px")
            )
          )
        )
      ),
      
      # --------------------------------------------------------------------------
      # TIME PRIORS TAB
      # --------------------------------------------------------------------------
      tabItem(
        tabName = "time",
        
        fluidRow(
          box(
            title = "Temporal Parameters",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            
            p("Configure time priors for population divergence, Ne changes, and migration changes."),
            
            tabsetPanel(
              id = "time_tabs",
              
              tabPanel(
                "Divergence Times",
                br(),
                DTOutput("table_divergence_times"),
                br(),
                actionButton("btn_edit_div_time",
                             "Edit Selected",
                             icon = icon("edit"),
                             class = "btn-primary")
              ),
              
              tabPanel(
                "Ne Change Times",
                br(),
                conditionalPanel(
                  condition = "input.check_ancestral_ne == true",
                  DTOutput("table_ne_change_times"),
                  br(),
                  actionButton("btn_edit_ne_time",
                               "Edit Selected",
                               icon = icon("edit"),
                               class = "btn-primary")
                ),
                conditionalPanel(
                  condition = "input.check_ancestral_ne == false",
                  p("Enable ancestral Ne changes in the Demography tab to configure timing.")
                )
              ),
              
              tabPanel(
                "Migration Change Times",
                br(),
                conditionalPanel(
                  condition = "input.check_migration == true",
                  DTOutput("table_mig_change_times"),
                  br(),
                  actionButton("btn_edit_mig_time",
                               "Edit Selected",
                               icon = icon("edit"),
                               class = "btn-primary")
                ),
                conditionalPanel(
                  condition = "input.check_migration == false",
                  p("Enable migration in the Migration tab to configure timing.")
                )
              )
            )
          )
        )
      ),
      
      # --------------------------------------------------------------------------
      # CONDITIONS TAB
      # --------------------------------------------------------------------------
      tabItem(
        tabName = "conditions",
        
        fluidRow(
          box(
            title = "Parameter Constraints",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            
            p("Define relationships between parameters (e.g., Ne0.pop1 > Ne0.pop2)."),
            
            tabsetPanel(
              id = "condition_tabs",
              
              tabPanel(
                "Size Conditions",
                br(),
                DTOutput("table_size_conditions"),
                br(),
                actionButton("btn_add_size_condition",
                             "Add Condition",
                             icon = icon("plus"),
                             class = "btn-success")
              ),
              
              tabPanel(
                "Migration Conditions",
                br(),
                conditionalPanel(
                  condition = "input.check_migration == true",
                  DTOutput("table_mig_conditions"),
                  br(),
                  actionButton("btn_add_mig_condition",
                               "Add Condition",
                               icon = icon("plus"),
                               class = "btn-success")
                ),
                conditionalPanel(
                  condition = "input.check_migration == false",
                  p("Migration is not enabled.")
                )
              ),
              
              tabPanel(
                "Time Conditions",
                br(),
                DTOutput("table_time_conditions"),
                br(),
                actionButton("btn_add_time_condition",
                             "Add Condition",
                             icon = icon("plus"),
                             class = "btn-success")
              )
            )
          )
        )
      ),
      
      # --------------------------------------------------------------------------
      # GENE SETUP TAB
      # --------------------------------------------------------------------------
      tabItem(
        tabName = "genes",
        
        fluidRow(
          box(
            title = "Genetic Data Configuration",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            
            fluidRow(
              column(4,
                     selectInput("select_data_type",
                                 "Data Type:",
                                 choices = c("Sanger sequencing" = "sanger",
                                             "Genomic (NGS)" = "genomic"),
                                 selected = "sanger")
              ),
              column(4,
                     numericInput("num_loci",
                                  "Number of Loci:",
                                  value = 10,
                                  min = 1,
                                  max = 10000)
              ),
              column(4,
                     selectInput("select_mut_dist",
                                 "Mutation Rate Distribution:",
                                 choices = c("Uniform" = "uniform",
                                             "Normal" = "normal"),
                                 selected = "uniform")
              )
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Locus Parameters",
            width = 12,
            status = "info",
            solidHeader = TRUE,
            
            DTOutput("table_loci"),
            
            br(),
            
            fluidRow(
              column(6,
                     actionButton("btn_edit_locus",
                                  "Edit Selected",
                                  icon = icon("edit"),
                                  class = "btn-primary")
              ),
              column(6,
                     actionButton("btn_auto_fill_loci",
                                  "Auto-fill Default Values",
                                  icon = icon("magic"),
                                  class = "btn-info")
              )
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Sample Sizes per Population",
            width = 12,
            status = "info",
            solidHeader = TRUE,
            
            DTOutput("table_samples"),
            
            br(),
            
            actionButton("btn_edit_samples",
                         "Edit Sample Sizes",
                         icon = icon("users"),
                         class = "btn-primary")
          )
        )
      ),
      
      # --------------------------------------------------------------------------
      # MODEL SUMMARY TAB
      # --------------------------------------------------------------------------
      tabItem(
        tabName = "summary",
        
        fluidRow(
          box(
            title = "Model Overview",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            
            actionButton("btn_refresh_summary",
                         "Refresh Summary",
                         icon = icon("sync"),
                         class = "btn-info"),
            
            hr(),
            
            h4("Model Configuration:"),
            verbatimTextOutput("txt_model_summary")
          )
        ),
        
        fluidRow(
          box(
            title = "Model Visualization",
            width = 12,
            status = "info",
            solidHeader = TRUE,
            
            checkboxInput("check_exponential",
                          "Use exponential size change (alpha)",
                          value = FALSE),
            
            conditionalPanel(
              condition = "input.check_exponential == true",
              textInput("txt_exp_pops",
                        "Populations with exponential change (comma-separated numbers):",
                        value = "1")
            ),
            
            actionButton("btn_plot_model",
                         "Generate Model Plot",
                         icon = icon("chart-area"),
                         class = "btn-success"),
            
            hr(),
            
            plotOutput("plot_model_diagram", height = "600px")
          )
        )
      ),
      
      # --------------------------------------------------------------------------
      # EXPORT TAB
      # --------------------------------------------------------------------------
      tabItem(
        tabName = "export",
        
        fluidRow(
          box(
            title = "Export Model",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            
            h4("Save your model for simulations:"),
            
            fluidRow(
              column(6,
                     textInput("txt_model_name",
                               "Model Name:",
                               value = "my_model",
                               placeholder = "Enter model name")
              ),
              column(6,
                     selectInput("select_export_format",
                                 "Export Format:",
                                 choices = c("R object (.rds)" = "rds",
                                             "R workspace (.RData)" = "RData"),
                                 selected = "rds")
              )
            ),
            
            br(),
            
            downloadButton("btn_download_model",
                           "Download Model File",
                           class = "btn-success btn-lg"),
            
            hr(),
            
            h4("Ready to run simulations?"),
            p("After downloading your model, use the following R commands:"),
            
            verbatimTextOutput("txt_simulation_code")
          )
        ),
        
        fluidRow(
          box(
            title = "Model Validation",
            width = 12,
            status = "warning",
            solidHeader = TRUE,
            
            actionButton("btn_validate_model",
                         "Validate Model",
                         icon = icon("check-double"),
                         class = "btn-warning"),
            
            br(), br(),
            
            uiOutput("model_validation_results")
          )
        )
      )
    )
  )
)
