
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
