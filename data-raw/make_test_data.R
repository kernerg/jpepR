library(data.table)
library(usethis)

# Example data paths
ld_file <- "data-raw/ld_blocks.txt"
exon_file <- "data-raw/Homo_sapiens.GRCh37.87_exon.bed"
meta_file <- "data-raw/Epimap_metadata.txt"
sampleID_file <- "data-raw/Epimap_sampleID.txt"
bckg_file <- "data-raw/Epimap_bckg_expectation.rds"

# Load files
ld_blocks <- fread(ld_file)
exons <- fread(exon_file)
metadata <- fread(meta_file)
sampleID <- fread(sampleID_file, header = F)
bckg_expectation <- readRDS(bckg_file)

# Save as .rda for use in examples or vignettes
use_data(ld_blocks, exons, metadata, sampleID,
         bckg_expectation, overwrite = TRUE)

