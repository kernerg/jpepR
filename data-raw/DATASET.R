

#' Example T2D finemapping file
#'
#' Used for building V_trait together with auxiliary traits.
#'
#' @format A data.table with columns: CHR, BP, A1, A2, RSID, PIP, BETA
#' @source Smith et al. (2024) Nat Med
"T2D_postfinemap.tsv"

#' Example auxiliary traits for T2D
#'
#' Used for building V_trait together with focal trait T2D.
#'
#' @format A character vector with auxiliary trait names.
#' @source Kerner et al. (2025)
"T2D_auxtraits.rda"

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
