# THIS SITE IS UNDER CONSTRUCTION. PLEASE BE PATIENT UNTIL A STABLE VERSION OF THE SOFTWARE IS AVAILABLE.
# jpepR

**jpepR** implements the *Joint Pleiotropic and Epigenomic Partitioning* (J-PEP) method, a framework for identifying interpretable SNP clusters that jointly reflect trait-specific pleiotropy and tissue-specific epigenomic enrichment.

---

## ✨ Features

- Fine-maps GWAS traits using P-values and LD blocks
- Builds SNP × trait (`V_trait`) and SNP × tissue (`V_tissue`) matrices
- Performs joint Bayesian non-negative matrix factorization (bNMF)
- Outputs interpretable SNP clusters, trait and tissue loadings

---

## 📦 Installation

You can install the development version of **jpepR** from GitHub:

```r
# Install from GitHub
install.packages("devtools")  # if not already installed
devtools::install_github("kernerg/jpepR")
```

## 🚀 Quick Start

```r
# Load the package
library(jpepR)

# Step 1: Fine-map focal and auxiliary traits
finemap_to_sqlite(
  trait = "T2D",
  sumstat_file = "path/to/T2D.sumstats.gz",
  ld_blocks = "path/to/ld_blocks.txt",
  output_dir = "results/"
)

# Step 2: Construct V_trait matrix
make_vtrait(
  focal_trait = "T2D",
  aux_traits = c("WHR", "BMI", "TG"),
  output_dir = "results/"
)

# Step 3: Construct V_tissue matrix
make_vtissue(
  focal_trait = "T2D",
  output_dir = "results/",
  epi_dir = "default"
)

# Step 4: Intersect SNPs across matrices
vmats <- intersect_vmatrices(
  focal_trait = "T2D",
  output_dir = "results/"
)

# Step 5: Run J-PEP model
results <- run_jpep_model(
  V_trait = vmats$V_trait,
  V_tissue = vmats$V_tissue
)

# Step 6: Save results to disk
save_jpep_results(
  consensus = results,
  focal_trait = "T2D",
  output_dir = "results/"
)
```

## 📘 Vignettes

You can explore the full J-PEP pipeline using built-in example data in the following vignette:

```r
browseVignettes("jpepR")
```

---

### 🧪 Dependencies

```markdown
## 🧪 Dependencies

The following packages are required or suggested for using `jpepR`:

### Imports
- `data.table`
- `DBI`
- `RSQLite`
- `glue`
- `Rmpfr`
- `igraph`
- `clue`

### Suggested (for development and vignettes)
- `devtools`
- `rmarkdown`
- `knitr`
```
---

## 📄 License

This package is released under the **MIT License**.  
See the [LICENSE](LICENSE) file for details.

---

## 👤 Author

**Gaspard Kerner**  
[gkerner@hsph.harvard.edu](mailto:gkerner@hsph.harvard.edu)

---

## 🔗 Related Projects

- [Original publication (Kerner et al., 2025, *medRxiv*)]
- [EpiMap (Boix et al., 2021, *Nature*)] 
- [Pleiotropic Decomposition Models(e.g. Ulder et al. 2018 *Plos Med*, Suzuki et al. 2024 *Nature* and Smith et al. 2024 *Nat Med*)]
