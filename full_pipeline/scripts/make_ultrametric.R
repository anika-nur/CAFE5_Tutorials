#!/usr/bin/env Rscript
# Turn OrthoFinder's rooted species tree (branch lengths in substitutions/site)
# into an ultrametric tree that CAFE5 will accept, using phytools.
#
# IMPORTANT CAVEAT: force.ultrametric() equalises root-to-tip distances by
# smoothing the existing (substitution) branch lengths. It does NOT perform
# divergence-time estimation, so the resulting branch lengths are ultrametric
# but NOT calibrated to real time. CAFE's lambda is therefore expressed in
# "events per unit of this smoothed distance". For a properly time-calibrated
# analysis you would instead date the tree (e.g. r8s / ape::chronos with fossil
# calibrations). We additionally rescale so the root height = 1 to give a clean
# relative time axis.
#
# Usage: Rscript make_ultrametric.R <in.tree> <out.tree>

suppressMessages({library(ape); library(phytools)})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: make_ultrametric.R <in.tree> <out.tree>")
infile <- args[1]; outfile <- args[2]

tr <- read.tree(infile)
cat("Input tree:", length(tr$tip.label), "tips\n")
cat("  rooted:      ", is.rooted(tr), "\n")
cat("  binary:      ", is.binary(tr), "\n")
cat("  ultrametric: ", is.ultrametric(tr), "\n")

if (!is.rooted(tr))  stop("tree is not rooted; CAFE requires a rooted tree")
if (!is.binary(tr))  { cat("  -> resolving polytomies with multi2di\n"); tr <- multi2di(tr) }

# phytools: make the tree ultrametric.
#   method="extend" lengthens only the terminal edges so every tip lines up at
#   the same depth. Unlike method="nnls" it never collapses an internal edge to
#   length 0 (which CAFE rejects as an "Invalid branch length"), so it is the
#   safe choice for feeding CAFE.
um <- force.ultrametric(tr, method = "extend")

# drop node labels / support values (CAFE wants a clean Newick) and rescale
um$node.label <- NULL
root_height <- max(node.depth.edgelength(um))
um$edge.length <- um$edge.length / root_height   # root-to-tip depth = 1

# safety guard: CAFE cannot have any zero-length branch
if (any(um$edge.length <= 0)) {
  eps <- 1e-6
  um$edge.length[um$edge.length <= 0] <- eps
  cat("  (bumped", sum(um$edge.length <= eps), "zero-length edge(s) to", eps, ")\n")
}

cat("\nOutput tree:\n")
cat("  ultrametric: ", is.ultrametric(um), "\n")
cat("  binary:      ", is.binary(um), "\n")
cat("  root->tip depth (rescaled): ", max(node.depth.edgelength(um)), "\n")

write.tree(um, file = outfile)
cat("\nWrote", outfile, "\n")
