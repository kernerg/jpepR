

#' Example summary statistics for a focal trait (simulated)
#'
#' @format A data.table with columns CHR, BP, RSID, A1, A2, P, etc.
#' @source Original data from [source]
"focal_sumstats"

#' Example summary statistics for auxiliary trait 1 (simulated)
#'
#' @format A data.table with columns CHR, BP, RSID, A1, A2, P, etc.
#' @source Original data from [source]
"auxiliary_01_sumstats"

#' Example summary statistics for auxiliary trait 2 (simulated)
#'
#' @format A data.table with columns CHR, BP, RSID, A1, A2, P, etc.
#' @source Original data from [source]
"auxiliary_02_sumstats"

#' Example summary statistics for auxiliary trait 3 (simulated)
#'
#' @format A data.table with columns CHR, BP, RSID, A1, A2, P, etc.
#' @source Original data from [source]
"auxiliary_03_sumstats"

#' Example summary statistics for auxiliary trait 4 (simulated)
#'
#' @format A data.table with columns CHR, BP, RSID, A1, A2, P, etc.
#' @source Original data from [source]
"auxiliary_04_sumstats"

#' Example summary statistics for auxiliary trait 5 (simulated)
#'
#' @format A data.table with columns CHR, BP, RSID, A1, A2, P, etc.
#' @source Original data from [source]
"auxiliary_05_sumstats"

#' Example LD block file
#'
#' Used for defining fine-mapping regions in J-PEP.
#'
#' @format A data.table with columns: LOCUS, CHR, START, STOP
#' @source Berisa & Pickrell (2016)
"ld_blocks"

#' Exon data
#'
#' Used to focus on non-exonic regulatory regions.
#'
#' @format A data.table with columns: CHR, START, STOP, gene, b, strand
#' @source Missing
"Homo_sapiens.GRCh37.87_exon.bed"

#' EpiMap epigenomic tracks - sampled to overlaping focal SNPs (in consideration of memory usage)
#'
#' Tracks from the original EpiMap publication.
#'
#' @format A data.table with columns: CHR, START, STOP, track, LENGTH
#' @source EpiMap publication (Boix et al. 2021 Nature)
"Epimap_tracks_overlapping_focal.tsv.gz"

#' EpiMap track background expectations
#'
#' Computed on 1KG Phase 3 chr 1
#'
#' @format A data.table Numeric vector with track names
#' @source EpiMap and 1KG Phase 3
"Epimap_bckg_expectation.rds"

#' Epimap metadata
#'
#' Used to link tracks to tissues.
#'
#' @format A data.table with tissue and other metadata information on tracks.
#' @source EpiMap
"Epimap_metadata.txt"


#' Epimap sample IDs
#'
#' Used to link Epimap tracks to sample IDs.
#'
#' @format A data.table with one column with sample IDs.
#' @source EpiMap
"Epimap_sampleID.txt"
