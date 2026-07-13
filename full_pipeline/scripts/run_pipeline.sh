#!/usr/bin/env bash
#
# End-to-end CAFE5 pipeline driver:
#   proteomes  --OrthoFinder-->  orthogroup counts + species tree
#              --python------->  CAFE input table (filtered)
#              --phytools----->  ultrametric species tree
#              --CAFE5-------->  gene-family evolution inference
#
# Prerequisites (installed once into a micromamba env named "cafepipe"):
#   orthofinder, diamond, mcl, fastme, r-base, r-phytools, r-ape
# and a compiled CAFE5 binary (built from this repository with cmake).
#
# Edit the three paths below for your machine, then: bash run_pipeline.sh
set -euo pipefail

# ---- configuration ---------------------------------------------------------
HERE="$(cd "$(dirname "$0")/.." && pwd)"          # full_pipeline/
ENV="${CAFE_ENV:-$HERE/mamba/envs/cafepipe}"      # micromamba env prefix
CAFE="${CAFE_BIN:-$HERE/../build/cafe5}"          # compiled cafe5 binary
PROTEOMES="${PROTEOMES:-$HERE/proteomes}"         # input .faa files (1/species)
THREADS="${THREADS:-4}"

export PATH="$ENV/bin:$PATH"                       # so orthofinder's python has numpy
mkdir -p "$HERE/results"

# ---- 1. OrthoFinder: cluster genes into orthogroups + infer species tree ---
echo "== [1/4] OrthoFinder =="
rm -rf "$HERE/orthofinder_out"
orthofinder -f "$PROTEOMES" -t "$THREADS" -a "$THREADS" -o "$HERE/orthofinder_out"
RES="$(dirname "$(find "$HERE/orthofinder_out" -name Orthogroups.GeneCount.tsv | head -1)")/.."

# ---- 2. Build the CAFE count table (filtered) ------------------------------
echo "== [2/4] format + filter gene counts =="
python3 "$HERE/scripts/orthogroups_to_cafe.py" \
    -i "$RES/Orthogroups/Orthogroups.GeneCount.tsv" \
    -o "$HERE/results/cafe_input.tsv"

# ---- 3. phytools: make the species tree ultrametric ------------------------
echo "== [3/4] phytools ultrametric tree =="
Rscript "$HERE/scripts/make_ultrametric.R" \
    "$RES/Species_Tree/SpeciesTree_rooted.txt" \
    "$HERE/results/species_tree.ultrametric.nwk"

# ---- 4. CAFE5: three analyses ----------------------------------------------
echo "== [4/4] CAFE5 =="
T="$HERE/results/species_tree.ultrametric.nwk"
I="$HERE/results/cafe_input.tsv"
"$CAFE" -i "$I" -t "$T" -p        -o "$HERE/results/run1_base_single"       # global lambda
"$CAFE" -i "$I" -t "$T" -p -k 3   -o "$HERE/results/run2_gamma_k3"          # among-family rate variation
"$CAFE" -i "$I" -t "$T" -p -e     -o "$HERE/results/run3_base_errormodel"   # + estimated error model

echo "Done. Results in $HERE/results/"
