#!/usr/bin/env bash
# ExpansionHunter v5 on the 126 bp tandem copy upstream of BRCA1 exon 1. Research use only.
#
# usage: eh_run.sh SAMPLES.tsv REFERENCE.fa OUTDIR [TRUTH.tsv]
#   SAMPLES.tsv  sample<TAB>CRAM URL or path (1000 Genomes 30x CRAM, GRCh38)
#   REFERENCE.fa GRCh38 FASTA with chr17 (indexed); also used to decode the CRAM
#   OUTDIR       per-sample region BAMs and ExpansionHunter JSON go here
#   TRUTH.tsv    optional: sample<TAB>carrier|reference; adds the truth column to the summary
# Environment: EXPANSIONHUNTER (default: ExpansionHunter), SAMTOOLS (default: samtools)
set -euo pipefail
SAMPLES=$1; REF=$2; OUT=$3; TRUTH=${4:-}
HERE=$(cd "$(dirname "$0")" && pwd)
EH=${EXPANSIONHUNTER:-ExpansionHunter}
SAMTOOLS=${SAMTOOLS:-samtools}
REGION=chr17:43127600-43128900
mkdir -p "$OUT"

grep -v '^#' "$SAMPLES" | while IFS=$'\t' read -r sample cram _; do
  [[ -z $sample ]] && continue
  bam="$OUT/$sample.region.bam"
  if [[ ! -s $bam ]]; then
    "$SAMTOOLS" view -b -T "$REF" -o "$bam.tmp" "$cram" "$REGION"
    mv "$bam.tmp" "$bam"
    "$SAMTOOLS" index "$bam"
  fi
  [[ -s $OUT/$sample.json ]] || "$EH" --reads "$bam" --reference "$REF" \
    --variant-catalog "$HERE/eh_catalog.json" --output-prefix "$OUT/$sample" >/dev/null 2>&1
done

python3 - "$OUT" "$TRUTH" <<'PY'
import csv, glob, json, os, sys
out, truth_path = sys.argv[1], sys.argv[2]
truth = {}
if truth_path:
    truth = {r[0]: r[1] for r in csv.reader(open(truth_path), delimiter="\t") if r and not r[0].startswith("#")}
print("sample\ttruth\teh_genotype\teh_ci\teh_spanning\teh_inrepeat")
for f in sorted(glob.glob(os.path.join(out, "*.json"))):
    s = os.path.basename(f)[:-5]
    v = json.load(open(f))["LocusResults"]["BRCA1_TC126"]["Variants"]["BRCA1_TC126"]
    print(s, truth.get(s, ""), v.get("Genotype", ""), v.get("GenotypeConfidenceInterval", ""),
          v.get("CountsOfSpanningReads", ""), v.get("CountsOfInrepeatReads", ""), sep="\t")
PY
