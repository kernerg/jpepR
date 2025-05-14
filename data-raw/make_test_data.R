library(data.table)
library(usethis)

# Example data paths
focal_file <- "data-raw/focal_sumstats.gz"
aux_files <- paste0("data-raw/auxiliary_0", 1:5, "_sumstats.gz")
ld_file <- "data-raw/ld_blocks.txt"
exon_file <- "data-raw/Homo_sapiens.GRCh37.87_exon.bed"
meta_file <- "data-raw/Epimap_metadata.txt"
sampleID_file <- "data-raw/Epimap_sampleID.txt"
track_file <- "data-raw/Epimap_tracks_overlapping_focal.tsv.gz"
bckg_file <- "data-raw/Epimap_bckg_expectation.rds"

# Load files
focal_sumstats <- fread(focal_file)
auxiliary_01 <- fread(aux_files[1])
auxiliary_02 <- fread(aux_files[2])
auxiliary_03 <- fread(aux_files[3])
auxiliary_04 <- fread(aux_files[4])
auxiliary_05 <- fread(aux_files[5])
ld_blocks <- fread(ld_file)
exons <- fread(exon_file)
metadata <- fread(meta_file)
sampleID <- fread(sampleID_file, header = F)
tracks <- fread(track_file)
bckg_expectation <- readRDS(bckg_file)

# Save as .rda for use in examples or vignettes
use_data(focal_sumstats, auxiliary_01, auxiliary_02, auxiliary_03,
         auxiliary_04, auxiliary_05, ld_blocks,
         exons, metadata, sampleID, tracks,
         bckg_expectation, overwrite = TRUE)

