# PipeMaster GUI - Complete Guide

## Overview

The PipeMaster GUI is a modern, user-friendly graphical interface for building complex population genetic models for ABC (Approximate Bayesian Computation) analysis. It replaces the command-line menu system with an intuitive web-based dashboard.

![PipeMaster GUI Dashboard](docs/images/dashboard_preview.png)

---

## Features

### ✨ **Key Benefits**

- **Visual Model Building**: Interactive interface for defining population structures
- **Real-time Validation**: Instant feedback on parameter settings
- **Tree Visualization**: Graphical display of population topologies
- **Smart Defaults**: Auto-fill common parameter values
- **Progress Tracking**: Status indicators show model completion
- **Export Ready**: Download models ready for simulation
- **No Command Line**: Point-and-click interface for all operations

### 📊 **Main Components**

1. **Population Structure** - Define tree topology with visual validation
2. **Demography** - Set Ne priors with interactive tables
3. **Migration** - Configure gene flow with matrix visualization
4. **Time Priors** - Set temporal parameters
5. **Conditions** - Define parameter constraints
6. **Gene Setup** - Configure loci and mutation rates
7. **Model Summary** - Overview and validation
8. **Export** - Download ready-to-use model files

---

## Installation

### Prerequisites

- **R** version ≥ 4.0.0
- **RStudio** (recommended but not required)

### Quick Install

```r
# Install required packages
install.packages(c(
  "shiny",
  "shinydashboard",
  "shinyjs",
  "shinyWidgets",
  "DT",
  "plotly",
  "ape",
  "phytools"
))

# Install PipeMaster (if not already installed)
devtools::install_github("gehara/PipeMaster")
```

### File Setup

Create a new folder for the GUI and place these three files in it:
- `app.R` - Main application file
- `app_ui.R` - User interface definition
- `app_server.R` - Server logic

```bash
# Directory structure
PipeMaster_GUI/
├── app.R
├── app_ui.R
└── app_server.R
```

---

## Launching the App

### Method 1: From RStudio

1. Open `app.R` in RStudio
2. Click the "Run App" button at the top of the editor
3. The app will open in a new window or browser tab

### Method 2: From R Console

```r
# Set working directory to the GUI folder
setwd("path/to/PipeMaster_GUI")

# Run the app
shiny::runApp()
```

### Method 3: Specify Port

```r
# Run on specific port
shiny::runApp(port = 8080)

# Run and automatically open in browser
shiny::runApp(launch.browser = TRUE)
```

---

## User Guide

### 🚀 **Getting Started**

#### Step 1: Start a New Model

1. Click **"Start New Model"** on the Getting Started page
2. The interface will initialize with default settings

#### Step 2: Define Population Structure

Navigate to **Population Structure** tab:

1. Enter your tree topology in Newick format:
   - Single population: `1`
   - Two populations: `(A,B);`
   - Four populations: `((A,B),(C,D));`

2. Click **"Validate Tree"** to check syntax

3. View the tree visualization and detected nodes

**Examples:**
```
# Three populations, bifurcating
((pop1,pop2),pop3);

# Four populations, symmetric
((A,B),(C,D));

# With branch lengths
((A:0.1,B:0.2):0.3,C:0.4);
```

#### Step 3: Configure Demography

Navigate to **Demography (Ne)** tab:

1. Select prior distribution (Uniform or Normal)
2. Edit current Ne priors by:
   - Clicking on a row in the table
   - Click "Edit Selected"
   - Enter min/max (uniform) or mean/SD (normal)

3. **Optional:** Enable ancestral Ne changes
   - Check "Include ancestral Ne changes"
   - Add ancestral Ne parameters for each population

#### Step 4: Set Up Migration (Optional)

Navigate to **Migration** tab:

1. Check "Enable migration between populations"
2. Select migration prior distribution
3. Edit migration rates (4Nm) for each population pair
4. View migration matrix visualization

#### Step 5: Configure Time Priors

Navigate to **Time Priors** tab:

Configure timing for:
- **Divergence Times**: When populations split
- **Ne Change Times**: When population sizes change
- **Migration Change Times**: When migration rates change

#### Step 6: Add Conditions (Optional)

Navigate to **Conditions** tab:

Define parameter relationships:
- `Ne0.pop1 > Ne0.pop2` - Pop1 larger than Pop2
- `join_1_2 < join_3_4` - First split before second

#### Step 7: Gene Setup

Navigate to **Gene Setup** tab:

1. Select data type (Sanger or Genomic)
2. Set number of loci
3. Choose mutation rate distribution
4. Edit locus parameters or use "Auto-fill Default Values"
5. Configure sample sizes per population

#### Step 8: Review and Export

Navigate to **Model Summary** tab:

1. Click "Refresh Summary" to see model overview
2. Optionally generate model plot
3. Go to **Export** tab
4. Click "Validate Model" to check for errors
5. Enter model name
6. Click "Download Model File"

---

## Detailed Features

### 📊 **Interactive Tables**

All parameter tables are **editable** and **sortable**:

- Click column headers to sort
- Select rows to edit
- Double-click cells for inline editing
- Use search boxes to filter

### 🎨 **Visualizations**

#### Tree Plots
- Automatic phylogenetic tree visualization
- Shows branch lengths and topology
- Updates in real-time as you edit

#### Migration Matrices
- Heatmap visualization of migration rates
- Interactive hover tooltips
- Color-coded by migration strength

#### Model Diagrams
- Visual representation of full model
- Shows populations, divergences, and migrations
- Export as image for publications

### ✅ **Validation**

Real-time validation for:
- ✓ Tree syntax (Newick format)
- ✓ Parameter ranges
- ✓ Required fields
- ✓ Logical consistency
- ✓ Model completeness

### 💾 **Save & Load**

#### Saving Models
- Download as `.rds` (recommended) or `.RData`
- Includes all parameters and settings
- Ready for use in simulations

#### Loading Models
- Upload previously saved models
- Resume editing from where you left off
- Import models from command-line menu

---

## Advanced Features

### 🔧 **Custom Configurations**

#### Exponential Growth
- Enable in Model Summary tab
- Specify which populations have exponential change
- Configure alpha parameters

#### Multiple Ancestral Ne
- Add multiple Ne changes per population
- Define timing for each change
- Set different priors for each period

#### Complex Migration
- Asymmetric migration matrices
- Time-varying migration rates
- Multiple migration periods

### 📝 **Parameter Constraints**

Use conditions to enforce relationships:

```
# Size constraints
Ne0.pop1 > Ne0.pop2           # Pop1 always larger
Ne1.pop1 = Ne0.pop1           # Constant size

# Time constraints
join_1_2 < join_3_4           # Sequential splits
t.Ne1.pop1 < join_1_2         # Ne change before split

# Migration constraints
mig0.1_2 < mig0.2_1          # Asymmetric migration
```

---

## Workflow Integration

### From GUI to Simulations

1. **Build model in GUI** → Download `.rds` file
2. **Load in R**:
```r
model <- readRDS("my_model.rds")
```

3. **Run simulations**:
```r
library(PipeMaster)

sim.ms.sumstat(
  model = model,
  nsim.blocks = 10,
  sim.block.size = 1000,
  output.name = "my_simulation",
  perpop.SS = TRUE,
  overall.SS = TRUE
)
```

4. **Analyze results**:
```r
# Load simulated data
sims <- read.table("my_simulation_popstats_mean.txt", header = TRUE)
pars <- read.table("my_simulation_pars.txt", header = TRUE)

# Run ABC
library(abc)
results <- abc(
  target = observed_data,
  param = pars,
  sumstat = sims,
  tol = 0.01,
  method = "rejection"
)
```

### From Command-Line to GUI

Import existing command-line models:

1. **Save command-line model**:
```r
model <- main.menu()  # Build in command line
saveRDS(model, "cli_model.rds")
```

2. **Load in GUI**:
- Click "Load Existing Model"
- Select `cli_model.rds`
- Continue editing visually

---

## Troubleshooting

### Common Issues

#### ❌ **App won't start**

**Problem**: Missing packages

**Solution**:
```r
# Check for missing packages
required <- c("shiny", "shinydashboard", "shinyjs", 
              "shinyWidgets", "DT", "plotly", "ape", "phytools")
missing <- required[!required %in% installed.packages()[,"Package"]]

if(length(missing) > 0) {
  install.packages(missing)
}
```

#### ❌ **Tree won't validate**

**Problem**: Invalid Newick format

**Solutions**:
- Check parentheses match: `(` = `)`
- Check commas: should equal number of `(`
- Remove spaces inside tree string
- Use example trees as templates

#### ❌ **Can't edit tables**

**Problem**: Row not selected

**Solution**: Click on a row first, then click "Edit Selected"

#### ❌ **Download fails**

**Problem**: Model incomplete

**Solution**:
1. Go to Export tab
2. Click "Validate Model"
3. Fix any errors shown
4. Try download again

### Performance Tips

#### For Large Models (100+ loci)

```r
# Run app with more memory
options(shiny.maxRequestSize = 100*1024^2)  # 100 MB
shiny::runApp()
```

#### For Slow Loading

```r
# Use local browser instead of RStudio viewer
shiny::runApp(launch.browser = TRUE)
```

---

## Customization

### Changing Colors

Edit `app_ui.R`:

```r
# Change primary color
tags$style(HTML("
  .btn-primary { background-color: #YOUR_COLOR; }
  .content-wrapper { background-color: #YOUR_BACKGROUND; }
"))
```

### Adding Custom Distributions

Edit `app_server.R` to add new distributions:

```r
selectInput("select_ne_dist",
            "Prior Distribution:",
            choices = c("Uniform" = "uniform",
                        "Normal" = "normal",
                        "Log-Normal" = "lognormal",  # NEW
                        "Gamma" = "gamma"))          # NEW
```

### Custom Validation Rules

Add to server validation section:

```r
# Custom rule: Ne must be > 100
if(any(rv$current_ne$Min < 100)) {
  errors <- c(errors, "Ne cannot be less than 100")
}
```

---

## Keyboard Shortcuts

| Shortcut | Action |
|----------|--------|
| `Ctrl + S` | Save current state |
| `Ctrl + R` | Refresh summary |
| `Ctrl + V` | Validate model |
| `Ctrl + E` | Export model |
| `Tab` | Navigate between fields |
| `Enter` | Submit forms |

---

## Best Practices

### ✅ **Do's**

- ✅ Validate tree before proceeding
- ✅ Save models frequently
- ✅ Use meaningful model names
- ✅ Check validation before export
- ✅ Test with small simulations first
- ✅ Document your parameter choices

### ❌ **Don'ts**

- ❌ Don't skip tree validation
- ❌ Don't use special characters in names
- ❌ Don't forget to set sample sizes
- ❌ Don't ignore validation warnings
- ❌ Don't mix distribution types
- ❌ Don't set impossible constraints

---

## Examples

### Example 1: Two-Population Model

```r
# Population Structure
Tree: (popA,popB);

# Demography
Ne0.popA: 100,000 - 500,000
Ne0.popB: 100,000 - 500,000

# Time
join_A_B: 50,000 - 500,000 years

# Genes
10 loci, 1000 bp each
Mutation rate: 5e-9 - 1.5e-8 per site per year
```

### Example 2: Four-Population with Migration

```r
# Population Structure
Tree: ((pop1,pop2),(pop3,pop4));

# Demography
Ne0.pop1: 10,000 - 100,000
Ne0.pop2: 10,000 - 100,000
Ne0.pop3: 10,000 - 100,000
Ne0.pop4: 10,000 - 100,000

# Migration
mig0.1_2: 0.1 - 5 (4Nm)
mig0.3_4: 0.1 - 5 (4Nm)

# Conditions
Ne0.pop1 > Ne0.pop2
join_1_2 < join_3_4
```

### Example 3: Complex Demography

```r
# Population Structure
Tree: ((A,B),C);

# Current Ne
Ne0.A: 50,000 - 200,000
Ne0.B: 50,000 - 200,000
Ne0.C: 100,000 - 500,000

# Ancestral Ne Changes
Ne1.A at time t.Ne1.A: 10,000 - 50,000
Ne1.B at time t.Ne1.B: 10,000 - 50,000

# Conditions
Ne0.A > Ne1.A  (expansion)
Ne0.B > Ne1.B  (expansion)
t.Ne1.A < join_A_B
```

---

## FAQ

**Q: Can I use this with my existing PipeMaster models?**
A: Yes! Load `.rds` files created from the command-line interface.

**Q: Does this work on all operating systems?**
A: Yes! Works on Windows, Mac, and Linux with R ≥ 4.0.0.

**Q: Can multiple people use the same model?**
A: Yes! Share the `.rds` file. Each person can load and edit independently.

**Q: How do I report bugs?**
A: Open an issue on the PipeMaster GitHub repository.

**Q: Can I run this on a server?**
A: Yes! Deploy to Shiny Server or shinyapps.io for web access.

**Q: Is my data secure?**
A: The app runs locally on your computer. No data is sent anywhere.

---

## Support

### Documentation
- PipeMaster GitHub: https://github.com/gehara/PipeMaster
- Shiny Tutorial: https://shiny.rstudio.com/tutorial/

### Getting Help
1. Check this guide
2. Review examples
3. Check troubleshooting section
4. Open GitHub issue
5. Contact package maintainer

### Contributing
Contributions welcome! Submit pull requests or feature requests on GitHub.

---

## Version History

### v1.0.0 (Current)
- Initial GUI release
- All core features implemented
- Basic validation and export
- Tree visualization
- Interactive tables

### Planned Features (v2.0.0)
- [ ] Undo/redo functionality
- [ ] Model comparison tools
- [ ] Batch model creation
- [ ] Advanced plotting options
- [ ] Cloud sync capabilities
- [ ] Collaborative editing

---

## License

This GUI extension follows the same license as PipeMaster. See the main package for details.

---

## Acknowledgments

Built with:
- **Shiny** by RStudio
- **shinydashboard** for layout
- **DT** for tables
- **plotly** for plots
- **ape** for phylogenetics

Created by Claude AI for the PipeMaster community.

---

**Happy Modeling! 🧬**
