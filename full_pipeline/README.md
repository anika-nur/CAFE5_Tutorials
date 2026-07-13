# A complete CAFE5 pipeline: OrthoFinder → phytools → CAFE5

This directory contains a full, reproducible worked example of a CAFE5 gene
family evolution analysis, starting from raw protein FASTA files and ending
with per-branch expansion/contraction inferences. It exercises every
preparatory step CAFE5 requires and then runs three CAFE5 models.

```
proteomes/*.faa
      │  OrthoFinder  (DIAMOND all-vs-all + MCL clustering + FastME species tree)
      ▼
Orthogroups.GeneCount.tsv          SpeciesTree_rooted.txt   (branch lengths = substitutions)
      │  orthogroups_to_cafe.py            │  make_ultrametric.R  (phytools force.ultrametric)
      ▼                                     ▼
results/cafe_input.tsv             results/species_tree.ultrametric.nwk
      │                                     │
      └──────────────── CAFE5 ──────────────┘
                          ▼
        λ estimate, ancestral counts, per-branch expansions/contractions,
        rapidly-evolving families
```

## The dataset

Six *Mycoplasma* proteomes (small bacterial genomes, 476–1130 proteins each)
shipped as OrthoFinder's bundled `ExampleData`. They were used here because the
session's network policy blocks the usual proteome sources (UniProt, Ensembl,
NCBI); this bundled set is fetched with the OrthoFinder package and is
purpose-built for exactly this clustering+tree workflow. To run the pipeline on
your own organisms, just drop one protein FASTA per species into `proteomes/`
and re-run — nothing else changes.

| Species | Proteins |
|---|---|
| Mycoplasma_agalactiae | 820 |
| Mycoplasma_arthritidis | 618 |
| Mycoplasma_gallisepticum | 763 |
| Mycoplasma_genitalium | 476 |
| Mycoplasma_haemocanis | 1130 |
| Mycoplasma_hyopneumoniae | 674 |

## Tools

| Step | Tool | How it was obtained |
|---|---|---|
| Orthogroup inference + species tree | OrthoFinder 3 (+ DIAMOND, MCL, FastME) | bioconda, via micromamba |
| Ultrametric tree | R + **phytools** (`force.ultrametric`) + ape | conda-forge, via micromamba |
| Gene family evolution | **CAFE5** | compiled from this repository with cmake/g++ |

Reproduce the toolchain:

```bash
# micromamba bootstrapped from conda-forge
micromamba create -y -p ./mamba/envs/cafepipe -c conda-forge -c bioconda \
    orthofinder r-base r-phytools r-ape
# CAFE5 built from the repo root:
cmake -S .. -B ../build -DCMAKE_BUILD_TYPE=Release && make -C ../build cafe5
```

## Running it

```bash
bash scripts/run_pipeline.sh
```

This regenerates everything under `results/`. Heavy/regenerable artifacts
(`proteomes/`, `orthofinder_out/`, the conda env) are git-ignored; the tracked
outputs are the scripts, the formatted CAFE input, the ultrametric tree, and
the CAFE5 result files.

## Preparatory steps in detail

### 1. OrthoFinder — defining gene families
OrthoFinder runs an all-vs-all DIAMOND search, clusters genes into orthogroups
with MCL, and infers a rooted species tree (STAG/STRIDE + FastME). Of 4481
genes, 84.4% fell into 686 orthogroups; 192 orthogroups contained all six
species. The two files we carry forward are
`Orthogroups/Orthogroups.GeneCount.tsv` (the family × species count matrix) and
`Species_Tree/SpeciesTree_rooted.txt`.

### 2. `orthogroups_to_cafe.py` — building the CAFE count table
Converts the OrthoFinder count matrix into CAFE's `Desc / Family ID / <species…>`
format and applies the tutorial's two filters:
* **presence filter** — drop families found in fewer than 2 species (CAFE assumes
  a family was present at the root; single-species families break ancestral
  reconstruction). **129 dropped.**
* **size filter** — divert families with ≥100 copies in any species to a
  `*.large.tsv` side file (high-variance families destabilise the λ search).
  **0 here** (the large *Mycoplasma_haemocanis*-specific expansions were already
  removed by the presence filter).

Result: **557 families** in `results/cafe_input.tsv`.

### 3. `make_ultrametric.R` — phytools time-tree conversion
OrthoFinder's tree has branch lengths in substitutions/site; CAFE requires an
**ultrametric** tree. `phytools::force.ultrametric(method="extend")` equalises
all root-to-tip distances (we then rescale root height to 1).

> **Caveat (see the main CAFE README "Time Trees" section):** `force.ultrametric`
> is *not* divergence-time estimation. The branch lengths become ultrametric but
> are not calibrated to real time, so λ is in "events per unit smoothed
> distance", not per million years. For a publication-grade analysis you would
> date the tree with fossil calibrations (e.g. r8s / `ape::chronos`). We use
> phytools here because the task called for it and it produces a valid,
> CAFE-acceptable tree. `method="extend"` is used rather than `"nnls"` because
> nnls collapsed the ingroup root branch to length 0, which CAFE rejects.

## CAFE5 results

### Run 1 — Base model, single global λ
```
Model Base Final Likelihood (-lnL): 1920.44
Lambda: 0.17623
```
* **32 / 557** families are significantly changing (p < 0.05,
  `Base_family_results.txt`).
* Per-branch expansions/contractions (`Base_clade_results.txt`) — the
  genome-reduced *M. genitalium* shows strong net **contraction** (10 gains vs
  54 losses), consistent with its known reductive evolution:

| Branch | Increase | Decrease |
|---|---|---|
| Mycoplasma_hyopneumoniae | 10 | 0 |
| Mycoplasma_haemocanis | 7 | 138 |
| Mycoplasma_agalactiae | 30 | 34 |
| Mycoplasma_arthritidis | 17 | 95 |
| Mycoplasma_gallisepticum | 19 | 7 |
| Mycoplasma_genitalium | 10 | 54 |

### Run 2 — Gamma model, among-family rate variation (k=3)
```
Model Gamma Final Likelihood (-lnL): 1892.12
Lambda: 0.17273
Alpha:  6.46143
```
The Gamma model lets each family fall into one of 3 discrete rate categories
drawn from a Gamma(α) distribution, with α estimated jointly with λ. Its
likelihood (1892.12) is better than the Base model's (1920.44), as expected
with the extra parameter. The high α (≈6.5) means the among-family rate
variation is modest — rates are fairly concentrated around the mean. The same
**32** families come out significant.

### Run 3 — Base model with an estimated error model
```
Model Base Final Likelihood (-lnL): 1903.21
Lambda:  0.14552
Epsilon: 0.02210
```
CAFE estimates a genome-assembly/annotation error rate (ε ≈ 0.022) directly
from the data (`Base_error_model.txt`) and re-estimates λ accounting for it.
**λ drops from 0.176 to 0.146**: uncorrected error artificially inflates the
apparent rate of gene-family change, so controlling for it lowers the estimate —
exactly the behaviour described in the CAFE tutorial.

### Model comparison summary
| Run | Model | −lnL | λ | other |
|---|---|---|---|---|
| 1 | Base, global λ | 1920.44 | 0.176 | 32 sig. families |
| 2 | Gamma, k=3 | 1892.12 | 0.173 | α = 6.46 |
| 3 | Base + error model | 1903.21 | 0.146 | ε = 0.022 |

## Output file guide (per model)
* `*_results.txt` — model, final likelihood, estimated λ (and α / ε).
* `*_family_results.txt` — per-family p-value and significance flag.
* `*_clade_results.txt` — expansions/contractions per branch.
* `*_asr.tre` — ancestral state reconstruction (Nexus), one tree per family.
* `*_count.tab` / `*_change.tab` — reconstructed counts and parent-child deltas.
* `*_branch_probabilities.tab` — per-node probability of the observed size.

Extract just the significantly-changing family trees for viewing in
Dendroscope/FigTree:
```bash
echo $'#nexus\nbegin trees;' > Significant_trees.tre
grep "*" results/run1_base_single/Base_asr.tre >> Significant_trees.tre
echo "end;" >> Significant_trees.tre
```
