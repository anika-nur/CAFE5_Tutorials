#!/usr/bin/env python3
"""
Convert an OrthoFinder `Orthogroups.GeneCount.tsv` table into a CAFE5 input
file, applying the two filters recommended in the CAFE tutorial:

  1. clade/presence filter -- keep only families with >=1 gene copy in at
     least two species. Families present in a single species cannot be used
     for ancestral reconstruction and violate CAFE's assumption that a family
     was present in the most recent common ancestor at the root.

  2. size filter -- families in which any single species has >= 100 gene
     copies are written to a separate `*.large.tsv` file. Such high-variance
     families destabilise the maximum-likelihood lambda search, so they are
     analysed afterwards with lambda fixed to the value estimated on the
     well-behaved families.

Output columns:  Desc <tab> "Family ID" <tab> <one column per species>
"""
import argparse
import sys

SIZE_CUTOFF = 100
MIN_SPECIES_PRESENT = 2


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--infile", required=True,
                    help="OrthoFinder Orthogroups.GeneCount.tsv")
    ap.add_argument("-o", "--outfile", required=True,
                    help="filtered CAFE input to write")
    args = ap.parse_args()

    with open(args.infile) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        # OrthoFinder header: Orthogroup <sp1> ... <spN> Total
        if header[0] != "Orthogroup":
            sys.exit("ERROR: does not look like an Orthogroups.GeneCount.tsv "
                     "(first column is %r, expected 'Orthogroup')" % header[0])
        species = header[1:]
        if species and species[-1] == "Total":
            species = species[:-1]

        kept, large, dropped = [], [], 0
        for line in fh:
            tok = line.rstrip("\n").split("\t")
            fam = tok[0]
            counts = [int(x) for x in tok[1:1 + len(species)]]

            n_present = sum(1 for c in counts if c >= 1)
            if n_present < MIN_SPECIES_PRESENT:
                dropped += 1
                continue

            row = ["(null)", fam] + [str(c) for c in counts]
            if any(c >= SIZE_CUTOFF for c in counts):
                large.append(row)
            else:
                kept.append(row)

    out_header = "\t".join(["Desc", "Family ID"] + species) + "\n"

    with open(args.outfile, "w") as out:
        out.write(out_header)
        for row in kept:
            out.write("\t".join(row) + "\n")

    large_path = None
    if large:
        large_path = args.outfile.rsplit(".", 1)[0] + ".large.tsv"
        with open(large_path, "w") as out:
            out.write(out_header)
            for row in large:
                out.write("\t".join(row) + "\n")

    print("species (%d):        %s" % (len(species), ", ".join(species)))
    print("families kept:       %d  -> %s" % (len(kept), args.outfile))
    print("families set aside:  %d  (>=%d copies in a species)%s"
          % (len(large), SIZE_CUTOFF,
             "  -> " + large_path if large_path else ""))
    print("families dropped:    %d  (present in <%d species)"
          % (dropped, MIN_SPECIES_PRESENT))


if __name__ == "__main__":
    main()
