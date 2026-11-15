#!/usr/bin/env bash
# Frozen Delly-versus-Manta test: simulated deletions and inversions (events.tsv) in real HG002 background reads,
# three depths; Delly 1.3.3 (rule of the pipeline, see PREREGISTRATION.md), Manta 1.6 and the read-pair research
# arm of bin/brca12_insertion_signals.py on the same reads. Also the earlier
# spike-in scenarios and the 40 simulated tandem duplications (../tandem_dup) as samples without a deletion
# or inversion.
# Usage: run.sh SPIKEDIR TDUPDIR REFERENCE OUTDIR
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(cd "$HERE/../.." && pwd)"
SP="${1:?}"; TD="${2:?}"; REF="${3:?}"; OUT="${4:?}"; THREADS="${THREADS:-4}"; MANTA="${MANTA:-configManta.py}"; DELLY="${DELLY:-delly}"
mkdir -p "$OUT"/{sim,res}
awk 'BEGIN{OFS="\t"} !/^#/ {print $1,$2,$3,$4}' "$ROOT/assets/brca12_regions.bed" > "$OUT/regions.bed"
[[ -s $OUT/genes.fa ]] || samtools faidx "$REF" chr17:43044295-43170327 chr13:32315077-32400268 > "$OUT/genes.fa"
awk 'BEGIN{OFS="\t"} $1!="chr13" && $1!="chr17" {print $1,0,$2}' "$REF.fai" > "$OUT/exclude.bed"
printf 'chr13\t32295077\t32420268\nchr17\t43024295\t43190327\n' | bgzip -c > "$OUT/manta_region.bed.gz"; tabix -f -p bed "$OUT/manta_region.bed.gz"
analyse() {
  local n=$1 bam=$2 depths=$3
  for f in $depths; do
    local t="$n.f$f"; [[ -s $OUT/res/$t.json && -s $OUT/res/$t.manta.vcf.gz && -s $OUT/res/$t.delly.vcf ]] && continue
    if [[ $f == 1 ]]; then cp "$bam" "$OUT/res/$t.bam"; else samtools view -b -s "$((7 + f)).$( [[ $f == 2 ]] && echo 5 || echo 333)" -o "$OUT/res/$t.bam" "$bam"; fi
    samtools index "$OUT/res/$t.bam"
    samtools view -h "$OUT/res/$t.bam" chr17:43044295-43170327 chr13:32315077-32400268 | python3 "$ROOT/bin/brca12_insertion_signals.py" \
      --regions "$OUT/regions.bed" --gene-fasta "$OUT/genes.fa" --read-pair-svs --out "$OUT/res/$t.json"
    "$DELLY" call -g "$REF" $( [[ -s $OUT/exclude.bed ]] && echo "-x $OUT/exclude.bed" ) -o "$OUT/res/$t.delly.bcf" "$OUT/res/$t.bam" > /dev/null 2>&1 || true
    if [[ -s $OUT/res/$t.delly.bcf ]]; then bcftools view -f PASS -i 'GT!="0/0" && GT!="./." && (SVTYPE=="DEL" || SVTYPE=="DUP" || SVTYPE=="INV" || SVTYPE=="BND")' "$OUT/res/$t.delly.bcf" > "$OUT/res/$t.delly.vcf"; else echo "##fileformat=VCFv4.2 (no Delly output)" > "$OUT/res/$t.delly.vcf"; fi
    rm -f "$OUT/res/$t.delly.bcf" "$OUT/res/$t.delly.bcf.csi"
    rm -rf "$OUT/res/$t.manta"
    "$MANTA" --bam "$OUT/res/$t.bam" --referenceFasta "$REF" --callRegions "$OUT/manta_region.bed.gz" --runDir "$OUT/res/$t.manta" > /dev/null
    "$OUT/res/$t.manta/runWorkflow.py" -j "$THREADS" > /dev/null 2>&1
    cp "$OUT/res/$t.manta/results/variants/diploidSV.vcf.gz" "$OUT/res/$t.manta.vcf.gz"
    rm -rf "$OUT/res/$t.manta" "$OUT/res/$t.bam" "$OUT/res/$t.bam.bai"
  done
}
grep -v '^#' "$HERE/events.tsv" | cut -f1 | while read -r n; do
  d="$OUT/sim/$n"
  if [[ ! -s $d/window.bam ]]; then
    seed=$((10#${n#*_})); [[ $n == INV_* ]] && seed=$((seed + 100))
    python3 "$ROOT/validation/pathogenic/spike.py" simulate --sample "$n" --events "$HERE/events.tsv" \
      --seed "$seed" --fragments "$SP/fragments.tsv.gz" --reference "$REF" --outdir "$d"
    bwa mem -t "$THREADS" -R '@RG\tID:spike\tSM:spike\tPL:ILLUMINA\tLB:spike' "$REF" \
      "$d/${n}_R1.fastq.gz" "$d/${n}_R2.fastq.gz" 2> "$d/bwa.log" | samtools sort -o "$d/window.bam" -
    samtools index "$d/window.bam"; rm -f "$d"/*.fastq.gz
  fi
  analyse "$n" "$d/window.bam" "1 2 3"
done
for n in SPIKE_NEG SPIKE_S1 SPIKE_S2 SPIKE_L1 SPIKE_L2 SPIKE_L3 SPIKE_L4 SPIKE_L5 SPIKE_L6; do
  [[ -s $SP/sim/$n/window.bam.bai ]] || samtools index "$SP/sim/$n/window.bam"
  analyse "$n" "$SP/sim/$n/window.bam" "1 2 3"
done
for d in "$TD"/sim/TDUP_*; do analyse "$(basename "$d")" "$d/window.bam" "1 2 3"; done
