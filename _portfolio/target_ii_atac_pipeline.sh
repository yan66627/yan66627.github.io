#!/bin/bash
# =============================================================================
# TaRGET II ATAC-seq: Batch HOMER Motif Scanning → Peak-to-Gene Pipeline
# =============================================================================
# Project : TaRGET II Consortium — Environmental Epigenomics
# Dataset : Male adult liver (Li, adt, M), mm10
# Purpose : For each transcription factor motif of interest, scan HOMER
#           ATAC-seq peak files, convert hits to BED, and map each peak
#           to the nearest TSS within 1 Mb (protein-coding + lncRNA).
# Input   : motif_list.txt  — tab-separated: <homer_dir> <motif_id> <gene_name>
# Output  : nearest_TSS_1Mb_pc_lncRNA_mm10/<sample>_<gene>_peak_gene_1Mb.txt
# Author  : Yan Jin
# =============================================================================

# -----------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------
INPUT_LIST="motif_list.txt"          # columns: homer_dir  motif_id  gene_name
TSS_BED="gencode_mm10_gene_TSS_pc_lncRNA.sorted.bed"   # pre-sorted TSS reference
TEMP_DIR="homer_preparsed_cache"     # HOMER genome pre-parse cache
SORTED_DIR="sorted_peaks"            # intermediate sorted BED files
OUT_DIR="nearest_TSS_1Mb_pc_lncRNA_mm10"

mkdir -p "$TEMP_DIR" "$SORTED_DIR" "$OUT_DIR"

echo "Input list  : $INPUT_LIST"
echo "TSS reference: $TSS_BED"
echo "============================================================"

# -----------------------------------------------------------------------
# Steps 1-5: For each motif entry — scan peaks, extract hits, make BED
# -----------------------------------------------------------------------
while read -r file motif gene; do

    echo ""
    echo ">>> FILE=${file}  MOTIF=${motif}  GENE=${gene}"

    OUT_RAW="${file}_${gene}.txt"
    OUT_PEAKS="${file}_${gene}_peak.txt"
    OUT_BED="${file}_${gene}_peak.bed"
    OUT_SORTED="${file}_${gene}_peak.sorted.bed"
    MOTIF_FILE="${file}/homerResults/${motif}.motif"

    # -- Step 1: Validate motif file -------------------------------------
    if [[ ! -f "$MOTIF_FILE" ]]; then
        echo "WARNING: motif file not found — $MOTIF_FILE  (skipping)"
        continue
    fi

    # -- Step 2: HOMER motif scan ----------------------------------------
    # -size given  : use peak sizes as defined in the input BED
    # -find        : score each peak against the specified motif matrix
    # -preparsedDir: reuse cached genome pre-parse to speed up batch runs
    findMotifsGenome.pl "${file}.txt" mm10 "${file}" \
        -size given \
        -find "$MOTIF_FILE" \
        -preparsedDir "$TEMP_DIR" \
        > "$OUT_RAW"
    echo "HOMER output : $OUT_RAW"

    # -- Step 3: Extract unique peak IDs ---------------------------------
    cut -f1 "$OUT_RAW" | sed '1d' | sort -u > "$OUT_PEAKS"
    echo "Peak list    : $OUT_PEAKS"

    # -- Step 4: Convert peak IDs to BED ---------------------------------
    # Peak ID format: chr_start_end  →  chr  start  end  chr_start_end
    awk -F'_' 'BEGIN{OFS="\t"} {print $1, $2, $3, $0}' \
        "$OUT_PEAKS" > "$OUT_BED"
    echo "BED file     : $OUT_BED"

    # -- Step 5: Sort BED by coordinate ----------------------------------
    sort -k1,1 -k2,2n "$OUT_BED" > "$OUT_SORTED"
    echo "Sorted BED   : $OUT_SORTED"
    echo "------------------------------------------------------------"

done < "$INPUT_LIST"

# -----------------------------------------------------------------------
# Step 6: Collect all sorted BED files
# -----------------------------------------------------------------------
mv *_peak.sorted.bed "$SORTED_DIR"/ 2>/dev/null
echo ""
echo "Moved sorted BED files → $SORTED_DIR/"

# -----------------------------------------------------------------------
# Step 7: bedtools closest — map each peak to nearest TSS within 1 Mb
# -----------------------------------------------------------------------
# -D a  : report distance signed relative to feature A
# -k 1  : report only the single nearest feature
# Filter: distance != "." AND |distance| <= 1,000,000 bp
for PEAK in "$SORTED_DIR"/*.bed; do
    BASENAME=$(basename "$PEAK" .bed)

    bedtools closest \
        -a "$PEAK" \
        -b "$TSS_BED" \
        -D a \
        -k 1 \
        -sorted \
    | awk 'BEGIN{OFS="\t"} $NF != "." && ($NF <= 1000000 && $NF >= -1000000)' \
        > "$OUT_DIR/${BASENAME}_peak_gene_1Mb.txt"

    echo "Nearest-TSS output: $OUT_DIR/${BASENAME}_peak_gene_1Mb.txt"
done

echo ""
echo "============================================================"
echo "Pipeline complete."
echo "============================================================"
