#!/bin/bash
# ============================================================
# Novel loci: LD-based clumping, run LOCALLY on your Mac with
# plink2 against the downloaded UKB sample (fileforclumping.bed/.bim/.fam).
# ============================================================
# This replaces the cluster-based run_plink_clump.sh -- since
# fileforclumping is a single bed/bim/fam covering (presumably) all
# chromosomes, rather than a UKB deployment split per-chromosome, we
# run ONE combined clump job instead of looping per chromosome.
#
# EDIT THESE THREE PATHS if your layout differs:
PLINK2="./plink2"                                   # the binary you downloaded, in Imaging_gen_v2
BFILE="$HOME/Downloads/fileforclumping"              # expects .bed/.bim/.fam with this prefix
INPUT="Outputs/Novel_loci/clump_input/all_chr_clump_input.txt"   # combined SNP+P file, already built
OUT="Outputs/Novel_loci/clump_output/all_chr_clumped"
CLUMP_KB=250     # PLINK default; widen to 500-1000 to be more conservative
CLUMP_R2=0.1     # PLINK default LD threshold for same-clump membership
# ============================================================

mkdir -p "$(dirname "$OUT")"

# Make sure the binary is executable (macOS sometimes strips this on download)
chmod +x "$PLINK2" 2>/dev/null

"$PLINK2" \
  --bfile "$BFILE" \
  --clump "$INPUT" \
  --clump-snp-field SNP \
  --clump-field P \
  --clump-p1 1 \
  --clump-p2 1 \
  --clump-kb ${CLUMP_KB} \
  --clump-r2 ${CLUMP_R2} \
  --out "$OUT"

echo ""
echo "Done. Output: ${OUT}.clumps (plink2) or ${OUT}.clumped (plink1.9-style, depending on your build)."
echo ""
echo "NOTE: plink2's --clump flag set has changed across versions/builds. If any"
echo "flag above errors with 'unrecognized', run:"
echo "    $PLINK2 --help clump"
echo "and drop/rename the offending flag (most likely candidates: --clump-p2 may"
echo "not exist in newer builds since plink2 simplified to a single p-value"
echo "threshold, and the output file extension may be .clumps instead of .clumped)."
echo ""
echo "Once this finishes, send back the output file(s) and I'll cross-reference"
echo "against snp_origin_lookup.csv to classify novel vs. already-known loci."
