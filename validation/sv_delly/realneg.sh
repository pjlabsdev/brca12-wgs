#!/usr/bin/env bash
# Real genomes: the 200 samples of ../longread/frozen/samples.tsv plus real deletion carriers
# (real_carriers.tsv), 1000 Genomes 30x reads of both gene regions.
# Delly (pipeline rule) and the read-pair research arm on each sample.
# Usage: realneg.sh REFERENCE OUTDIR   (OUTDIR/regions.bed, genes.fa and exclude.bed from run.sh)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(cd "$HERE/../.." && pwd)"
REF="${1:?}"; OUT="${2:?}"; DELLY="${DELLY:-delly}"; mkdir -p "$OUT/real"
cat <(grep -v '^#' "$ROOT/validation/longread/frozen/samples.tsv") <(grep -v '^#' "$HERE/real_carriers.tsv" | cut -f1,2) |
while IFS=$'\t' read -r s url; do
  [[ -s $OUT/real/$s.json && -s $OUT/real/$s.delly.vcf ]] && continue
  ok=0
  for i in 1 2 3; do
    (cd "$OUT/real" && samtools view -h -T "$REF" "$url" chr17:43044295-43170327 chr13:32315077-32400268 > "$s.sam.tmp" 2>/dev/null) && ok=1 && break
    sleep 5
  done
  [[ $ok == 1 ]] || { echo "$s fetch failed" >&2; rm -f "$OUT/real/$s.sam.tmp"; continue; }
  mv "$OUT/real/$s.sam.tmp" "$OUT/real/$s.sam"
  python3 "$ROOT/bin/brca12_insertion_signals.py" --regions "$OUT/regions.bed" --gene-fasta "$OUT/genes.fa" \
    --read-pair-svs --out "$OUT/real/$s.json" < "$OUT/real/$s.sam"
  # keep only the @SQ lines of chr13 and chr17, so that the BAM matches the chr13/chr17 test reference
  awk '!/^@SQ/ || /SN:chr1[37]\t/' "$OUT/real/$s.sam" | samtools sort -o "$OUT/real/$s.bam" - 2>/dev/null
  samtools index "$OUT/real/$s.bam"
  "$DELLY" call -g "$REF" $( [[ -s $OUT/exclude.bed ]] && echo "-x $OUT/exclude.bed" ) -o "$OUT/real/$s.delly.bcf" "$OUT/real/$s.bam" > /dev/null 2>&1 || true
  if [[ -s $OUT/real/$s.delly.bcf ]]; then bcftools view -f PASS -i 'GT!="0/0" && GT!="./." && (SVTYPE=="DEL" || SVTYPE=="DUP" || SVTYPE=="INV" || SVTYPE=="BND")' "$OUT/real/$s.delly.bcf" > "$OUT/real/$s.delly.vcf"; else echo "##fileformat=VCFv4.2 (no Delly output)" > "$OUT/real/$s.delly.vcf"; fi
  rm -f "$OUT/real/$s.sam" "$OUT/real/$s.bam" "$OUT/real/$s.bam.bai" "$OUT/real/$s.delly.bcf"* "$OUT/real"/*.crai
done
