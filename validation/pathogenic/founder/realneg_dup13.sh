#!/usr/bin/env bash
# Real genomes as negatives for the BRCA1 exon 13 duplication: the 200 samples of validation/longread/frozen/samples.tsv,
# 1000 Genomes 30x short reads of chr17:43076000-43087000, analysed as in run.sh (all reads only).
# Usage: realneg.sh REFERENCE OUTDIR   (OUTDIR/regions.bed and OUTDIR/genes.fa from run.sh)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(cd "$HERE/../../.." && pwd)"
REF="${1:?}"; OUT="${2:?}"; mkdir -p "$OUT/real"
grep -v '^#' "$ROOT/validation/longread/frozen/samples.tsv" | while IFS=$'\t' read -r s url; do
  [[ -s $OUT/real/$s.json ]] && continue
  ok=0
  for i in 1 2 3; do
    (cd "$OUT/real" && samtools view -h -T "$REF" "$url" chr17:43076000-43087000 > "$s.sam.tmp" 2>/dev/null) && ok=1 && break
    sleep 5
  done
  [[ $ok == 1 ]] || { echo "$s fetch failed" >&2; rm -f "$OUT/real/$s.sam.tmp"; continue; }
  mv "$OUT/real/$s.sam.tmp" "$OUT/real/$s.sam"
  samtools depth "$OUT/real/$s.sam" | awk '$1=="chr17" && $2==43084362 {print $3}' > "$OUT/real/$s.depth"
  python3 "$ROOT/bin/brca12_insertion_signals.py" --regions "$OUT/regions.bed" --gene-fasta "$OUT/genes.fa" \
    --out "$OUT/real/$s.json" < "$OUT/real/$s.sam"
  rm -f "$OUT/real/$s.sam" "$OUT/real"/*.crai
done
