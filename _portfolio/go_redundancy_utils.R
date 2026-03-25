# =============================================================================
# GO Term Redundancy Reduction Utilities
# Author: Yan Jin
#
# Functions for reducing redundant GO terms from enrichment analysis results
# using semantic similarity (rrvgo), and for mapping parent GO terms back to
# their associated gene lists.
#
# Dependencies: rrvgo, org.Mm.eg.db (or other OrgDb), dplyr, tidyr
# =============================================================================

# -----------------------------------------------------------------------------
# run_rrvgo()
#
# Computes GO semantic similarity and reduces redundant terms using rrvgo.
# Scores are derived from -log10(q-value), so more significant terms are
# weighted higher during parent term selection.
#
# Args:
#   df        : data.frame with GO enrichment results (e.g., from GREAT/clusterProfiler)
#   go_col    : column name containing GO IDs (default "ID")
#   qval_col  : column name containing FDR q-values (default "q-value FDR B&H")
#   orgdb     : OrgDb identifier string (default "org.Mm.eg.db" for mouse)
#   ont       : GO ontology: "BP", "MF", or "CC" (default "BP")
#   method    : similarity method passed to calculateSimMatrix (default "Rel")
#   threshold : similarity threshold for merging terms (default 0.7)
#
# Returns:
#   list with:
#     $scores       — named numeric vector of -log10(q-value) per GO ID
#     $simMatrix    — semantic similarity matrix
#     $reducedTerms — data.frame with parent term assignments
# -----------------------------------------------------------------------------
run_rrvgo <- function(df,
                      go_col    = "ID",
                      qval_col  = "q-value FDR B&H",
                      orgdb     = "org.Mm.eg.db",
                      ont       = "BP",
                      method    = "Rel",
                      threshold = 0.7) {
  library(rrvgo)

  stopifnot(go_col   %in% colnames(df))
  stopifnot(qval_col %in% colnames(df))

  # Build score vector: -log10(q-value), higher = more significant
  scores <- setNames(-log10(df[[qval_col]]), df[[go_col]])
  scores <- scores[is.finite(scores)]   # remove NA / Inf

  # Compute pairwise semantic similarity
  simMatrix <- calculateSimMatrix(
    names(scores),
    orgdb  = orgdb,
    ont    = ont,
    method = method
  )

  # Reduce redundant terms
  reducedTerms <- suppressMessages(
    reduceSimMatrix(
      simMatrix,
      scores,
      threshold = threshold,
      orgdb     = orgdb
    )
  )

  list(
    scores       = scores,
    simMatrix    = simMatrix,
    reducedTerms = reducedTerms
  )
}


# -----------------------------------------------------------------------------
# extract_parent_go_genes()
#
# After redundancy reduction, maps each representative (parent) GO term back
# to its member GO terms and extracts the union of associated genes.
#
# Args:
#   rrvgo_obj       : output list from run_rrvgo()
#   enrich_df       : original enrichment data.frame (same as used in run_rrvgo)
#   go_col_enrich   : column name for GO IDs in enrich_df (default "ID")
#   gene_col_enrich : column name for gene lists in enrich_df
#                     (default "Hit in Query List", comma-separated)
#   sep             : separator used in the gene list column (default ",")
#
# Returns:
#   list with:
#     $parent_go_map       — data.frame mapping each GO term to its parent
#     $parent_go_gene_long — long-format data.frame (parent, go, gene)
#     $parent_gene_summary — one row per parent: gene list + count
#     $all_parent_genes    — character vector of all unique genes
# -----------------------------------------------------------------------------
extract_parent_go_genes <- function(rrvgo_obj,
                                    enrich_df,
                                    go_col_enrich   = "ID",
                                    gene_col_enrich = "Hit in Query List",
                                    sep             = ",") {
  library(dplyr)
  library(tidyr)

  stopifnot("reducedTerms" %in% names(rrvgo_obj))
  stopifnot(all(c("parent", "go") %in% colnames(rrvgo_obj$reducedTerms)))
  stopifnot(go_col_enrich   %in% colnames(enrich_df))
  stopifnot(gene_col_enrich %in% colnames(enrich_df))

  # Step 1: parent → member GO mapping from rrvgo output
  parent_go_map <- rrvgo_obj$reducedTerms[, c("parent", "go")]

  # Step 2: subset enrichment table to GO terms present in rrvgo output
  go_enrichment_subset <- enrich_df[enrich_df[[go_col_enrich]] %in% parent_go_map$go, ]

  # Step 3: join parent-GO map with enrichment results
  parent_go_gene_df <- merge(
    parent_go_map,
    go_enrichment_subset,
    by.x = "go",
    by.y = go_col_enrich,
    all.x = TRUE
  )
  parent_go_gene_df <- parent_go_gene_df[, c("parent", "go", gene_col_enrich)]
  colnames(parent_go_gene_df)[3] <- "genes"

  # Step 4: expand comma-separated gene lists to one gene per row
  parent_go_gene_df <- tidyr::separate_rows(parent_go_gene_df, genes, sep = sep)

  # Step 5: summarise unique genes per parent term
  parent_gene_summary <- parent_go_gene_df |>
    dplyr::group_by(parent) |>
    dplyr::summarise(
      genes      = list(unique(genes)),
      gene_count = dplyr::n_distinct(genes),
      .groups    = "drop"
    )

  list(
    parent_go_map       = parent_go_map,
    parent_go_gene_long = parent_go_gene_df,
    parent_gene_summary = parent_gene_summary,
    all_parent_genes    = unique(parent_go_gene_df$genes)
  )
}


# -----------------------------------------------------------------------------
# batch_export_rrvgo()
#
# Loops over a named list of run_rrvgo() results and saves heatmap, scatter,
# treemap plots and reducedTerms CSV for each entry.
#
# Args:
#   rrvgo_list : named list of run_rrvgo() outputs
#   outdir     : output directory path (created if it does not exist)
#   prefix     : filename prefix for outputs (default "motif")
#   width_*    : plot widths in inches
#   height_*   : plot heights in inches
# -----------------------------------------------------------------------------
batch_export_rrvgo <- function(rrvgo_list,
                               outdir,
                               prefix     = "motif",
                               width_heat  = 8, height_heat  = 8,
                               width_scat  = 7, height_scat  = 6,
                               width_tree  = 8, height_tree  = 6) {
  library(rrvgo)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  for (name in names(rrvgo_list)) {
    message("Exporting: ", name)
    obj <- rrvgo_list[[name]]
    base <- file.path(outdir, paste0(prefix, "_", name))

    # CSV
    write.csv(obj$reducedTerms, paste0(base, "_reducedTerms.csv"), row.names = FALSE)

    # Heatmap
    pdf(paste0(base, "_heatmap.pdf"), width = width_heat, height = height_heat)
    heatmapPlot(obj$simMatrix, obj$reducedTerms,
                annotateParent = TRUE, annotationLabel = "parentTerm", fontsize = 6)
    dev.off()

    # Scatter plot
    pdf(paste0(base, "_scatter.pdf"), width = width_scat, height = height_scat)
    scatterPlot(obj$simMatrix, obj$reducedTerms)
    dev.off()

    # Treemap
    pdf(paste0(base, "_treemap.pdf"), width = width_tree, height = height_tree)
    treemapPlot(obj$reducedTerms)
    dev.off()
  }

  invisible(NULL)
}


# =============================================================================
# Example usage
# =============================================================================
#
# library(org.Mm.eg.db)
#
# # Run GO redundancy reduction for one motif's enrichment result
# rrvgo_motif1 <- run_rrvgo(
#   df        = motif1_enrichment,   # data.frame with GO IDs and q-values
#   go_col    = "ID",
#   qval_col  = "q-value FDR B&H",
#   orgdb     = "org.Mm.eg.db",
#   threshold = 0.7
# )
#
# # Extract genes associated with each parent GO term
# res_motif1 <- extract_parent_go_genes(
#   rrvgo_obj = rrvgo_motif1,
#   enrich_df = motif1_enrichment
# )
# res_motif1$parent_gene_summary
#
# # Batch export plots and tables for multiple motifs
# rrvgo_list <- list(motif1 = rrvgo_motif1, motif2 = rrvgo_motif2)
# batch_export_rrvgo(rrvgo_list, outdir = "results/GO_plots", prefix = "TCDD_F")
