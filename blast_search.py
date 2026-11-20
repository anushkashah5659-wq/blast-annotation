#!/usr/bin/env python3
"""Search NCBI BLAST with a query sequence, record the search parameters,
and turn the hits into an annotation table plus a short interpretation.

Search needs Biopython (pip install biopython) and internet access.
Parsing/annotation uses the standard library only, so --xml can re-process
a saved BLAST XML file offline.
"""

import argparse
import json
import sys
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path


def read_first_fasta(path):
    """Return (id, sequence) of the first record in a FASTA file."""
    rec_id, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if rec_id is not None:
                    break
                rec_id = line[1:].split()[0]
            elif rec_id is not None:
                chunks.append(line)
    if rec_id is None:
        sys.exit("No FASTA record found in " + path)
    return rec_id, "".join(chunks).upper()


def run_blast(seq, params, xml_path):
    """Send the query to NCBI BLAST (remote) and save the raw XML."""
    from Bio.Blast import NCBIWWW  # imported here so offline parsing needs no Biopython

    handle = NCBIWWW.qblast(
        params["program"], params["database"], seq,
        expect=params["evalue_threshold"],
        hitlist_size=params["max_hits"],
        megablast=params["megablast"],
        format_type="XML",
    )
    Path(xml_path).write_text(handle.read())


def text(node, tag, default=""):
    el = node.find(tag)
    return el.text if el is not None and el.text is not None else default


def union_length(intervals):
    """Number of query positions covered by a set of (start, end) intervals."""
    total, cur_s, cur_e = 0, None, None
    for s, e in sorted((min(a, b), max(a, b)) for a, b in intervals):
        if cur_e is None or s > cur_e + 1:
            if cur_e is not None:
                total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        total += cur_e - cur_s + 1
    return total


def parse_blast_xml(xml_path):
    root = ET.parse(xml_path).getroot()
    info = {
        "blast_program": text(root, "BlastOutput_program"),
        "blast_version": text(root, "BlastOutput_version"),
        "database_searched": text(root, "BlastOutput_db"),
        "query_length": int(text(root, "BlastOutput_query-len", "0")),
        "evalue_threshold_used": text(root, "BlastOutput_param/Parameters/Parameters_expect"),
        "gap_open": text(root, "BlastOutput_param/Parameters/Parameters_gap-open"),
        "gap_extend": text(root, "BlastOutput_param/Parameters/Parameters_gap-extend"),
        "filter": text(root, "BlastOutput_param/Parameters/Parameters_filter"),
        "match_score": text(root, "BlastOutput_param/Parameters/Parameters_sc-match"),
        "mismatch_score": text(root, "BlastOutput_param/Parameters/Parameters_sc-mismatch"),
    }
    qlen = info["query_length"]
    hits = []
    for hit in root.iter("Hit"):
        hsps = list(hit.iter("Hsp"))
        if not hsps:
            continue
        best = max(hsps, key=lambda h: float(text(h, "Hsp_bit-score", "0")))
        intervals = [(int(text(h, "Hsp_query-from")), int(text(h, "Hsp_query-to"))) for h in hsps]
        aln = int(text(best, "Hsp_align-len", "0"))
        ident = int(text(best, "Hsp_identity", "0"))
        hits.append({
            "accession": text(hit, "Hit_accession"),
            "description": text(hit, "Hit_def"),
            "subject_length": int(text(hit, "Hit_len", "0")),
            "query_cover_pct": round(100 * union_length(intervals) / qlen, 1) if qlen else 0.0,
            "pct_identity": round(100 * ident / aln, 2) if aln else 0.0,
            "aln_length": aln,
            "n_hsps": len(hsps),
            "evalue": float(text(best, "Hsp_evalue", "nan")),
            "bit_score": float(text(best, "Hsp_bit-score", "0")),
        })
    return info, hits


def classify(h, evalue_threshold=1e-5):
    """Simple, documented rules for reading a hit (see README)."""
    if h["evalue"] > evalue_threshold:
        return "Not significant"
    if h["pct_identity"] >= 97 and h["query_cover_pct"] >= 80:
        return "Near-identical (same species/strain)"
    if h["pct_identity"] >= 80 and h["query_cover_pct"] >= 50:
        return "Related (same genus/family level)"
    return "Partial/local similarity only"


def write_outputs(info, hits, qid, outdir, top_n):
    annotated = []
    for i, h in enumerate(hits, 1):
        annotated.append({"rank": i, **h, "interpretation": classify(h)})

    cols = ["rank", "accession", "description", "subject_length", "query_cover_pct",
            "pct_identity", "aln_length", "n_hsps", "evalue", "bit_score", "interpretation"]
    top = annotated[:top_n]
    with open(Path(outdir) / "annotation_table.tsv", "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in top:
            out.write("\t".join(f"{r[c]:.3g}" if c in ("evalue",) else str(r[c]) for c in cols) + "\n")

    counts = {}
    for r in annotated:
        counts[r["interpretation"]] = counts.get(r["interpretation"], 0) + 1

    lines = [f"# Automatic interpretation of BLAST hits for `{qid}`", "",
             f"- Query length: {info['query_length']} bp",
             f"- Program / database: {info['blast_program']} / {info['database_searched']}",
             f"- Hits returned: {len(annotated)}", ""]
    if annotated:
        b = annotated[0]
        lines += [f"**Best hit:** {b['accession']} - {b['description']}",
                  f"({b['pct_identity']}% identity, {b['query_cover_pct']}% query cover, E = {b['evalue']:.2g}) "
                  f"-> {b['interpretation']}.", "", "**Hits by category:**", ""]
        lines += [f"- {k}: {v}" for k, v in sorted(counts.items(), key=lambda kv: -kv[1])]
    else:
        lines.append("No hits were returned.")
    lines += ["", "How categories are decided (E-value <= 1e-5 required for all except 'Not significant'):",
              "- Near-identical: identity >= 97% and query cover >= 80%",
              "- Related: identity >= 80% and query cover >= 50%",
              "- Otherwise: partial/local similarity only", "",
              "These are rules of thumb. Read the description column yourself before drawing conclusions.", ""]
    (Path(outdir) / "interpretation_auto.md").write_text("\n".join(lines))


def main():
    p = argparse.ArgumentParser(description="BLAST search + annotation table.")
    p.add_argument("query", nargs="?", help="query FASTA (first record is used)")
    p.add_argument("--start", type=int, help="optional 1-based start of region to use as query")
    p.add_argument("--end", type=int, help="optional 1-based end of region to use as query")
    p.add_argument("--program", default="blastn")
    p.add_argument("--database", default="nt")
    p.add_argument("--evalue", type=float, default=1e-5)
    p.add_argument("--max-hits", type=int, default=25)
    p.add_argument("--no-megablast", action="store_true", help="use plain blastn instead of megablast")
    p.add_argument("--top", type=int, default=10, help="rows in annotation table")
    p.add_argument("--outdir", default="output")
    p.add_argument("--xml", help="skip the search and process this saved BLAST XML file")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    xml_path = Path(args.xml) if args.xml else outdir / "blast_raw.xml"

    qid = "query"
    if not args.xml:
        if not args.query:
            sys.exit("Give a query FASTA file (or --xml to re-process a saved result).")
        qid, seq = read_first_fasta(args.query)
        if args.start or args.end:
            s, e = (args.start or 1), (args.end or len(seq))
            seq = seq[s - 1:e]
            qid = f"{qid}:{s}-{e}"
        params = {
            "query_id": qid, "query_length": len(seq),
            "program": args.program, "database": args.database,
            "evalue_threshold": args.evalue, "max_hits": args.max_hits,
            "megablast": not args.no_megablast,
            "run_date_utc": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S"),
        }
        run_blast(seq, params, xml_path)
        (outdir / "search_parameters.json").write_text(json.dumps(params, indent=2) + "\n")

    info, hits = parse_blast_xml(xml_path)
    # Parameters reported back by BLAST itself (version, scoring, filter) are kept too.
    (outdir / "blast_reported_settings.json").write_text(json.dumps(info, indent=2) + "\n")
    write_outputs(info, hits, qid, outdir, args.top)
    print(f"{len(hits)} hits parsed; wrote annotation_table.tsv and interpretation_auto.md to {outdir}/")


if __name__ == "__main__":
    main()
