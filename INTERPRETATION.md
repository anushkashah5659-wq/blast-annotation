# Interpretation of BLAST results

**Query:** positions 1-800 of the phiX174 genome (NC_001422.1), searched with
blastn (megablast) against NCBI nt, E-value threshold 1e-5, up to 25 hits.

## What the hits are

All 10 hits in the annotation table are complete phiX174 genomes (5,386 bp):
three records described as "evolved" J1, J2 and J3 lines (MH378441-MH378443)
and seven isolates from the HM7536xx-HM7537xx series. Every one matches the
query with 100% identity over the full 800 bp (100% query cover), with an
E-value of 0 and an identical bit score of 1478.44.

## How reliable the matches are

The matches are extremely strong. An E-value of 0 means a match this good would
not occur by chance, and 100% identity over 100% of the query leaves no doubt
that the query comes from phiX174. This was expected, because the query was cut
from a phiX174 reference genome, so the search effectively confirms the query's
identity.

## Limits of the interpretation

- The hits all have the same score, so their order is arbitrary. Rank 1 is not
  a better match than rank 10.
- The script labels these hits "near-identical (same species/strain)". The
  species call is safe, but an 800 bp region cannot distinguish strains, since
  the evolved lines and isolates are identical in this region. Telling them
  apart would need the full genome or a more variable region.
- Row 6 combines two identical records (HM753673 and HM753674), because BLAST
  merges identical sequences into one hit.
- Only the top 10 of 25 requested hits are in the table. They show no variety,
  so this search cannot show how phiX174 differs from related phages. A
  longer or more variable query, or a search that excludes phiX174, would be
  needed for that.

## Conclusion

The query is phiX174 DNA. The best matches are many sequenced phiX174 genomes,
including laboratory-evolved lines, and the 800 bp region is fully conserved
across them.
