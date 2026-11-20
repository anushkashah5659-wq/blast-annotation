# blast-annotation

Search NCBI's public nucleotide database (nt) with a query sequence using BLAST,
record the search parameters, interpret the hits, and produce an annotation table
of the best matches.

## What is in this repo

| File | What it is |
|---|---|
| `blast_search.py` | Runs the BLAST search and builds the tables |
| `.github/workflows/blast.yml` | Runs the script on GitHub (manual "Run workflow" button) |
| `data/phix174.fasta` | Source genome (phiX174, NCBI NC_001422.1), downloaded by the workflow |
| `output/blast_raw.xml` | Raw BLAST output, unedited |
| `output/search_parameters.json` | Parameters I set (program, database, E-value, query, date) |
| `output/blast_reported_settings.json` | Settings BLAST reported back (version, scoring, filter) |
| `output/annotation_table.tsv` | Annotation table of the top 10 hits |
| `output/interpretation_auto.md` | Rule-based summary of the hits |
| `INTERPRETATION.md` | My written interpretation of the results |

## How to run

On GitHub: **Actions** tab -> **BLAST search** -> **Run workflow**. It takes a few
minutes, then commits the results into `output/`.

On your own computer:

    pip install biopython
    python3 blast_search.py data/phix174.fasta --start 1 --end 800 --outdir output

Any single- or multi-record FASTA works as input (the first record is used).
Re-process a saved result without re-running the search: `python3 blast_search.py --xml output/blast_raw.xml`

## Search parameters

- Query: positions 1-800 of phiX174 (NC_001422.1), 800 bp
- Program: blastn (megablast), remote search at NCBI
- Database: nt
- E-value threshold: 1e-5
- Max hits requested: 25
- Low-complexity filter: BLAST default
- The run date is stored in `output/search_parameters.json` because nt changes over time.

## How hits are interpreted

For each hit I report percent identity, query cover (share of the query covered by
all matching segments), E-value and bit score. Hits are put into categories:

- **Near-identical**: identity >= 97% and query cover >= 80%
- **Related**: identity >= 80% and query cover >= 50%
- **Partial/local only**: anything else that is significant
- **Not significant**: E-value above the threshold

These cut-offs are rules of thumb I chose, not official BLAST rules.
Identity is taken from the best-scoring segment of each hit.

## Decisions

- Used an 800 bp slice of a known genome so the search is fast and the query is real, public data.
- Kept the table to 10 hits to keep it readable.
- The automatic summary is only a starting point; my own reading is in `INTERPRETATION.md`.
