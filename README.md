
# jpepR

**jpepR** implements the *Joint Pleiotropic and Epigenomic Partitioning* (J-PEP) method, a framework for identifying interpretable SNP clusters that jointly reflect trait-specific pleiotropy and tissue-specific epigenomic enrichment.

> ⚠️ **THIS SITE IS UNDER CONSTRUCTION. PLEASE BE PATIENT UNTIL A STABLE VERSION OF THE SOFTWARE IS AVAILABLE.**

---

## 📦 Installation

You can install the development version of **jpepR** from GitHub:

```r
# Install from GitHub
install.packages("devtools")  # if not already installed
devtools::install_github("kernerg/jpepR")
```

---

## 🧪 Dependencies

### Imports
- `data.table`
- `DBI`
- `RSQLite`
- `glue`
- `Rmpfr`
- `igraph`
- `clue`

### Suggested (for development, vignettes, and plotting)
- `devtools`
- `rmarkdown`
- `knitr`
- `pheatmap`
- `RColorBrewer`

---

## 📘 How to Use

To learn how to use **jpepR**, please refer to the detailed vignette:

```r
browseVignettes("jpepR")
```

or view the source directly at:

➡️ [`vignettes/jpepR-test.Rmd`](vignettes/jpepR-test.Rmd)

This vignette walks through the full J-PEP workflow using example data, from fine-mapping and V matrix construction to joint factorization and visualization.

---

## 🚀 Quick Test Run (Minimal Example)

This **minimal test** downloads example data and runs J-PEP from start to finish.  
For a complete walkthrough, see the full vignette (`vignettes/jpepR-test.Rmd`).

```r
# Load jpepR
library(jpepR)

# Set up a working directory
setwd("/n/groups/price/gaspard/jpepR_clean_env")
dir.create("data-raw", showWarnings = FALSE)

# Download example test data (~300MB)
jpepR::download_test_data(dest_dir = "data-raw")

# Define inputs
trait <- "T2D"
focal_loc <- list(T2D = "data-raw/T2D_postfinemap.tsv")
load("data-raw/T2D_auxtraits.rda")  # loads 'subset_aux'
output_dir <- file.path(tempdir(), "jpep_test_output")

# Construct V matrices
make_vtrait(trait, NULL, output_dir, fine_mapped_files = focal_loc,
            full_matrix_path = "data-raw/big_pleio_matrix.tsv.gz",
            subset = subset_aux)
make_vtissue(trait, output_dir, epi_dir = "default")
vmats <- intersect_vmatrices(trait, output_dir)

# Run J-PEP model
runs <- run_jpep_model("JPEP", vmats$V_trait, vmats$V_tissue,
                       Tolerance = 1e-7, k = 15, NbrReps = 10, attempts = 5)

# Process results
results <- concensus_run("JPEP", runs, vmats$V_trait, vmats$V_tissue)

# Retrieve matrices
W <- results$W
H_trait <- results$H_trait
H_tissue <- results$H_tissue
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
- [Pleiotropic decomposition models (Ulder et al., 2018 *PLoS Med*; Suzuki et al., 2024 *Nature*; Smith et al., 2024 *Nat Med*)]
