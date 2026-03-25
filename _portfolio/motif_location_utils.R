# =============================================================================
# Motif Location Utilities
# Author: Yan Jin
#
# Functions for mapping motif hits from HOMER/FIMO output back to genomic
# coordinates, and for resolving overlapping motifs within the same peak.
#
# Dependencies: IRanges, S4Vectors (Bioconductor)
# =============================================================================

# -----------------------------------------------------------------------------
# compute_motif_locations()
#
# Given a motif hit table (from HOMER/FIMO) and a peak coordinate table,
# computes the genomic start/end of each motif hit based on peak center
# and motif offset.
#
# Args:
#   target_df    : data.frame with columns: PositionID, Offset, Sequence,
#                  Motif.Name, Strand, MotifScore
#   coord_df     : data.frame with columns: Chrom, Start, End, PositionID
#   motif_len    : integer, length of the motif (default 15)
#   enforce_order: logical, ensure MotifStart <= MotifEnd (default TRUE)
#   round_fun    : function used to round peak center (default round)
#
# Returns:
#   data.frame with genomic coordinates of each motif hit
# -----------------------------------------------------------------------------
compute_motif_locations <- function(target_df, coord_df,
                                    motif_len    = 15,
                                    enforce_order = TRUE,
                                    round_fun    = round) {
  req_target <- c("PositionID", "Offset", "Sequence", "Motif.Name", "Strand", "MotifScore")
  req_coord  <- c("Chrom", "Start", "End", "PositionID")

  if (!all(req_target %in% names(target_df))) {
    stop("target_df missing columns: ", paste(setdiff(req_target, names(target_df)), collapse = ", "))
  }
  if (!all(req_coord %in% names(coord_df))) {
    stop("coord_df missing columns: ", paste(setdiff(req_coord, names(coord_df)), collapse = ", "))
  }

  # Left join: keep all rows from target_df
  df <- merge(target_df, coord_df, by = "PositionID", all.x = TRUE, sort = FALSE)

  # Compute peak center
  center <- round_fun((df$Start + df$End) / 2)

  # Motif start based on offset from center
  motif_start <- center + as.integer(df$Offset)

  # Compute motif end based on strand
  is_plus <- df$Strand %in% c("+", "plus", "POS", "pos")
  motif_end <- ifelse(is_plus,
                      motif_start + motif_len,
                      motif_start - motif_len)

  # Ensure start <= end if requested
  if (enforce_order) {
    start_out <- pmin(motif_start, motif_end, na.rm = TRUE)
    end_out   <- pmax(motif_start, motif_end, na.rm = TRUE)
  } else {
    start_out <- motif_start
    end_out   <- motif_end
  }

  out <- data.frame(
    Chrom      = df$Chrom,
    MotifStart = as.integer(start_out),
    MotifEnd   = as.integer(end_out),
    PositionID = df$PositionID,
    Offset     = as.integer(df$Offset),
    Sequence   = df$Sequence,
    MotifName  = df[["Motif.Name"]],
    Strand     = df$Strand,
    MotifScore = df$MotifScore,
    PeakStart  = df$Start,
    PeakEnd    = df$End,
    PeakCenter = as.integer(center),
    stringsAsFactors = FALSE
  )

  # Adjust endpoints by strand (half-open interval convention):
  # (+) strand: MotifEnd - 1
  idx_plus  <- which(out$Strand == "+" & !is.na(out$MotifEnd))
  if (length(idx_plus))  out$MotifEnd[idx_plus]   <- out$MotifEnd[idx_plus]   - 1L

  # (-) strand: MotifStart + 1
  idx_minus <- which(out$Strand == "-" & !is.na(out$MotifStart))
  if (length(idx_minus)) out$MotifStart[idx_minus] <- out$MotifStart[idx_minus] + 1L

  out
}


# -----------------------------------------------------------------------------
# select_best_by_overlap()
#
# Within each peak (PositionID), identifies clusters of overlapping motif hits
# and retains only the highest-scoring hit per cluster.
#
# Args:
#   df        : data.frame output from compute_motif_locations()
#   id_col    : column name for peak ID (default "PositionID")
#   start_col : column name for motif start (default "MotifStart")
#   end_col   : column name for motif end (default "MotifEnd")
#   score_col : column name for motif score (default "MotifScore")
#
# Returns:
#   data.frame with one representative (best-scoring) motif per overlap cluster
# -----------------------------------------------------------------------------
select_best_by_overlap <- function(df,
                                   id_col    = "PositionID",
                                   start_col = "MotifStart",
                                   end_col   = "MotifEnd",
                                   score_col = "MotifScore") {
  library(IRanges)

  need <- c(id_col, start_col, end_col, score_col)
  miss <- setdiff(need, names(df))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))

  groups <- split(df, df[[id_col]], drop = TRUE)

  picked_list <- lapply(groups, function(g) {
    st <- as.integer(g[[start_col]])
    en <- as.integer(g[[end_col]])
    sc <- as.numeric(g[[score_col]])

    ir  <- IRanges(start = st, end = en)
    red <- reduce(ir, with.revmap = TRUE)
    revmap <- S4Vectors::mcols(red)$revmap

    keep_idx <- integer(length(revmap))
    for (i in seq_along(revmap)) {
      idxs     <- as.integer(revmap[[i]])
      sc_local <- sc[idxs]
      sc_local[is.na(sc_local)] <- -Inf
      # If all scores are NA, keep the first entry in the cluster
      best <- if (all(is.infinite(sc_local) & sc_local < 0)) idxs[1L] else idxs[which.max(sc_local)]
      keep_idx[i] <- best
    }

    g[sort(unique(keep_idx)), , drop = FALSE]
  })

  out <- do.call(rbind, picked_list)
  rownames(out) <- NULL

  # Sort by peak ID then genomic position
  if (start_col %in% names(out)) {
    out <- out[order(out[[id_col]], out[[start_col]], out[[end_col]]), , drop = FALSE]
  }

  out
}


# =============================================================================
# Example usage
# =============================================================================
#
# motif_hits <- read.table("AR_full_motif.txt", sep = "\t", header = TRUE)
# peak_coords <- read.table("MergedChipseq.bed", sep = "\t", header = FALSE,
#                           col.names = c("Chrom", "Start", "End", "PositionID"))
#
# # Step 1: map motif hits to genomic coordinates
# motif_locs <- compute_motif_locations(motif_hits, peak_coords, motif_len = 16)
#
# # Step 2: keep only best non-overlapping hit per peak
# motif_locs_filtered <- select_best_by_overlap(motif_locs)
#
# write.csv(motif_locs_filtered, "motif_locations_filtered.csv", row.names = FALSE)
