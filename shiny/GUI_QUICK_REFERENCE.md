# PipeMaster GUI - Quick Reference Card

## 🚀 Quick Start

```r
# 1. Install packages
source("setup_gui.R")

# 2. Launch app
setwd("PipeMaster_GUI")
shiny::runApp()
```

---

## 📋 Model Building Checklist

- [ ] **Population Structure**: Define tree topology
- [ ] **Demography**: Set Ne priors  
- [ ] **Migration**: Configure gene flow (optional)
- [ ] **Time Priors**: Set temporal parameters
- [ ] **Conditions**: Add constraints (optional)
- [ ] **Gene Setup**: Configure loci and mutation
- [ ] **Validate**: Check model completeness
- [ ] **Export**: Download model file

---

## 🌳 Tree Format Examples

| Populations | Tree String |
|-------------|-------------|
| 1 | `1` |
| 2 | `(A,B);` |
| 3 | `((A,B),C);` |
| 4 | `((A,B),(C,D));` |
| With lengths | `((A:0.1,B:0.2):0.3,C:0.4);` |

---

## 🎛️ Parameter Types

### Ne (Effective Population Size)
- **Current Ne**: Ne at present time
- **Ancestral Ne**: Ne at past times
- **Units**: 4Nm (diploid)

### Migration (4Nm)
- **Current**: Gene flow at present
- **Ancestral**: Gene flow in past
- **Range**: Typically 0.1 - 10

### Time (years)
- **Divergence**: When populations split
- **Ne change**: When size changes
- **Migration change**: When flow changes

### Mutation (per site per year)
- **Typical range**: 1e-9 to 1e-7
- **Sanger**: ~1e-8
- **Genomic**: Variable

---

## 🔧 Common Tasks

### Edit Parameter
1. Click row in table
2. Click "Edit Selected"
3. Enter new values
4. Click "Save"

### Add Condition
1. Go to Conditions tab
2. Click "Add Condition"
3. Format: `param1 < param2`
4. Options: `<`, `>`, `=`

### Auto-fill Loci
1. Gene Setup tab
2. Click "Auto-fill Default Values"
3. Set common parameters
4. Apply to all loci

### Validate Model
1. Export tab
2. Click "Validate Model"
3. Fix any errors
4. Download when valid

---

## 📊 Distribution Types

### Uniform
- **Parameters**: Min, Max
- **Use**: Wide uncertainty
- **Example**: Ne 10,000 - 100,000

### Normal
- **Parameters**: Mean, SD
- **Use**: Known approximate value
- **Example**: μ=50,000, σ=10,000

---

## ⚡ Keyboard Shortcuts

| Key | Action |
|-----|--------|
| `Ctrl+S` | Save state |
| `Ctrl+R` | Refresh |
| `Ctrl+V` | Validate |
| `Ctrl+E` | Export |
| `Tab` | Next field |
| `Enter` | Submit |

---

## 🐛 Troubleshooting

### Tree won't validate
- ✓ Check parentheses match
- ✓ Check commas = ( count
- ✓ Remove spaces
- ✓ Use examples

### Can't edit table
- ✓ Click row first
- ✓ Then click "Edit"

### Download fails
- ✓ Validate model first
- ✓ Fix all errors
- ✓ Try again

### App won't start
- ✓ Run `setup_gui.R`
- ✓ Check R version ≥ 4.0
- ✓ Install packages

---

## 💡 Best Practices

### ✅ Do
- Validate tree early
- Save models often
- Use meaningful names
- Check validation
- Test small first

### ❌ Don't
- Skip validation
- Use special chars
- Forget samples
- Ignore warnings
- Mix distributions

---

## 📁 File Formats

### Model Files
- **`.rds`**: Recommended
- **`.RData`**: Alternative
- Contains all parameters

### Output Files
- `*_popstats_mean.txt`
- `*_overallstats_mean.txt`
- `*_pars.txt`

---

## 🔗 Workflow

```
GUI → Build Model → Download .rds
  ↓
Load in R → Run Simulations → ABC Analysis
  ↓
Results → Interpret → Publish
```

---

## 💻 Simulation Code

```r
# Load model
model <- readRDS("my_model.rds")

# Run simulations
library(PipeMaster)

sim.ms.sumstat(
  model = model,
  nsim.blocks = 10,
  sim.block.size = 1000,
  output.name = "sim"
)

# Analyze with ABC
library(abc)
results <- abc(
  target = obs,
  param = pars,
  sumstat = sims,
  tol = 0.01
)
```

---

## 📞 Support

- **Docs**: GUI_USER_GUIDE.md
- **GitHub**: github.com/gehara/PipeMaster
- **Issues**: Submit via GitHub
- **Email**: Package maintainer

---

## 🎯 Example Models

### Two-Pop Simple
```
Tree: (A,B);
Ne: 100k-500k each
Divergence: 50k-500k ya
10 loci, 1000bp
```

### Four-Pop Migration
```
Tree: ((1,2),(3,4));
Ne: 10k-100k each
Migration: 0.1-5 (4Nm)
join_1_2 < join_3_4
```

### Complex Demography
```
Tree: ((A,B),C);
Current Ne: 50k-200k
Ancestral Ne: 10k-50k
Expansion in A,B
```

---

**Version 1.0** | Created with ❤️ for PipeMaster
