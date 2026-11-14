#!/usr/bin/env bash
# Real genomes: the 200 samples of ../longread/frozen/samples.tsv, 1000 Genomes 30x reads of both gene regions.
# Usage: realneg.sh REFERENCE OUTDIR   (OUTDIR/regions.bed and genes.fa from run.sh)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(cd "$HERE/../.." && pwd)"
REF="${1:?}"; OUT="${2:?}"; mkdir -p "$OUT/real"
grep -v '^#' "$ROOT/validation/longread/frozen/samples.tsv" | while IFS=$'\t' read -r s url; do
  [[ -s $OUT/real/$s.json ]] && continue
  ok=0
  for i in 1 2 3; do
    (cd "$OUT/real" && samtools view -h -T "$REF" "$url" chr17:43044295-43170327 chr13:32315077-32400268 > "$s.sam.tmp" 2>/dev/null) && ok=1 && break
    sleep 5
  done
  [[ $ok == 1 ]] || { echo "$s fetch failed" >&2; rm -f "$OUT/real/$s.sam.tmp"; continue; }
  mv "$OUT/real/$s.sam.tmp" "$OUT/real/$s.sam"
  python3 "$ROOT/bin/brca12_insertion_signals.py" --regions "$OUT/regions.bed" --gene-fasta "$OUT/genes.fa" \
    --out "$OUT/real/$s.json" < "$OUT/real/$s.sam"
  rm -f "$OUT/real/$s.sam" "$OUT/real"/*.crai
done
