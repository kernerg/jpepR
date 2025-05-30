
# jpepR

**jpepR** implements the *Joint Pleiotropic and Epigenomic Partitioning* (J-PEP) method, a framework for identifying interpretable SNP clusters that jointly reflect trait-specific pleiotropy and tissue-specific epigenomic enrichment.

> ⚠️ **THIS SITE IS UNDER CONSTRUCTION. PLEASE BE PATIENT UNTIL A STABLE VERSION OF THE SOFTWARE IS AVAILABLE.**

> 🚧 **In the meantime, feel free to try it — the example data and pipeline should run smoothly!**

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

## 🚀 Quick Test Run (Not so quick though)

This **minimal test** downloads example data and runs J-PEP from start to finish.  
For a complete walkthrough, see the full vignette (`vignettes/jpepR-test.Rmd`).

```r
# Load jpepR
library(jpepR)

# Set up a working directory
setwd("/path/to/your/folder/jpepR_clean_env")
dir.create("data-raw", showWarnings = FALSE)

# Download example test data (~900MB)
jpepR::download_test_data(dest_dir = "data-raw")

# Define inputs (T2D Fine-map [single-causal] from Smith et al. 2024 Nat Med)
trait             <- "T2D"
trait_loc         <- list(T2D = "data-raw/T2D_postfinemap.tsv")
aux_trait_mat_loc <- "data-raw/big_pleio_matrix.tsv.gz"
output_dir        <- file.path(tempdir(), "jpep_test_output")
load("data-raw/T2D_auxtraits.rda")  # loads default auxiliary trait subset for T2D

# Build V matrices (⚠️ make_vtissue may take a few minutes)
make_vtrait(
  trait,
  NULL,
  output_dir,
  fine_mapped_files = trait_loc,
  full_matrix_path = aux_trait_mat_loc,
  subset = subset_aux
)
make_vtissue(
  trait,
  output_dir,
  epi_dir = "default"
)

# Aligns SNPs in V_trait and V_tissue
vmats <- intersect_vmatrices(trait, output_dir)

# Run J-PEP (👀 use "Pleiotropic" or "Epigenomic" for other models)
runs <- run_jpep_model(
  "JPEP",
  vmats$V_trait,
  vmats$V_tissue,
  Tolerance = 1e-7,
  k = 15,
  NbrReps = 10,
  attempts = 5
)

# Process results
results <- concensus_run(
  "JPEP",
  runs,
  vmats$V_trait,
  vmats$V_tissue
)

# Retrieve W and H matrices
W <- results$W
H_trait <- results$H_trait
H_tissue <- results$H_tissue
```

---

## 🧬 How to Run J-PEP on Your Own Trait

To run J-PEP on a trait of your choice, simply adjust the inputs shown in the **🚀 Quick Test Run** section:

- **Trait name**  
  Replace `"T2D"` with your desired trait name in `trait`.
  Replace `"T2D"` with your desired trait name in `trait_loc`'s list name.

- **Fine-mapped summary statistics**  
  Update `trait_loc` to point to your fine-mapped file (TSV or similar), which must contain the following columns:

  | CHR | BP | RSID | A1 | A2 | BETA | PIP |
  |-----|----|------|----|----|------|-----|
  | 1–22 | hg19 position | SNP ID | Effect allele | Non-effect allele | Effect size for A1 | Posterior inclusion probability |

- **SNP matching**  
  By default, SNPs are matched using the `RSID` column.
  Make sure RSIDs in your focal trait file are consistent with those used in the auxiliary trait matrix.

- **Auxiliary trait matrix**  
  If not specified, J-PEP will use a default matrix of ~1 million SNPs × 164 traits, fine-mapped under a *single causal variant* model.

- **Subset of auxiliary traits** *(optional)*  
  You can restrict analysis to a subset of traits by passing a character vector `subset_aux`, listing preferred auxiliary traits from the full matrix (check column names of **Auxiliary trait matrix**).

> 🔍 For a complete example, refer to the vignette:  
> [`vignettes/jpepR-test.Rmd`](vignettes/jpepR-test.Rmd)

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
