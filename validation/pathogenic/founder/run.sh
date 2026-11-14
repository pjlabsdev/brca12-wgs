#!/usr/bin/env bash
# Simulated BRCA2 c.156_157insAlu carriers and spike-in negatives in real HG002 background reads.
# Usage: run.sh SPIKEDIR REFERENCE OUTDIR
#   SPIKEDIR: work dir of ../run_spikeins.sh (fragments.tsv.gz, sim/<scenario>/window.bam)
# Every window BAM is analysed at three depths (all reads, 1/2, 1/3) with
# bin/brca12_insertion_signals.py --gene-fasta (catalog realignment and read-end signature) and Manta.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(cd "$HERE/../../.." && pwd)"
SP="${1:?}"; REF="${2:?}"; OUT="${3:?}"; THREADS="${THREADS:-4}"; MANTA="${MANTA:-configManta.py}"
mkdir -p "$OUT"/{sim,res}
awk 'BEGIN{OFS="\t"} !/^#/ {print $1,$2,$3,$4}' "$ROOT/assets/brca12_regions.bed" > "$OUT/regions.bed"
[[ -s $OUT/genes.fa ]] || samtools faidx "$REF" chr17:43044295-43170327 chr13:32315077-32400268 > "$OUT/genes.fa"
printf 'chr13\t32300000\t32340000\n' | bgzip -c > "$OUT/manta_region.bed.gz"; tabix -f -p bed "$OUT/manta_region.bed.gz"

analyse() {   # name bam
  local n=$1 bam=$2
  for f in 1 2 3; do
    local t="$n.f$f"; [[ -s $OUT/res/$t.json && -s $OUT/res/$t.manta.vcf.gz ]] && continue
    if [[ $f == 1 ]]; then cp "$bam" "$OUT/res/$t.bam"; else samtools view -b -s "$((7 + f)).$( [[ $f == 2 ]] && echo 5 || echo 333)" -o "$OUT/res/$t.bam" "$bam"; fi
    samtools index "$OUT/res/$t.bam"
    samtools depth -a -r chr13:32319165-32319165 "$OUT/res/$t.bam" | cut -f3 > "$OUT/res/$t.depth"
    samtools view -h "$OUT/res/$t.bam" chr13:32315077-32400268 | python3 "$ROOT/bin/brca12_insertion_signals.py" \
      --regions "$OUT/regions.bed" --gene-fasta "$OUT/genes.fa" --out "$OUT/res/$t.json"
    rm -rf "$OUT/res/$t.manta"
    "$MANTA" --bam "$OUT/res/$t.bam" --referenceFasta "$REF" --callRegions "$OUT/manta_region.bed.gz" \
      --runDir "$OUT/res/$t.manta" > /dev/null
    "$OUT/res/$t.manta/runWorkflow.py" -j "$THREADS" > /dev/null 2>&1
    cp "$OUT/res/$t.manta/results/variants/diploidSV.vcf.gz" "$OUT/res/$t.manta.vcf.gz"
    rm -rf "$OUT/res/$t.manta" "$OUT/res/$t.bam" "$OUT/res/$t.bam.bai"
  done
}

for fa in "$HERE"/elements/*.fa; do
  n=$(basename "$fa" .fa); d="$OUT/sim/$n"
  if [[ ! -s $d/window.bam ]]; then
    python3 "$ROOT/validation/pathogenic/spike.py" simulate --sample "$n" --events "$HERE/events.tsv" \
      --alu-sequence "$fa" --fragments "$SP/fragments.tsv.gz" --reference "$REF" --outdir "$d"
    bwa mem -t "$THREADS" -R '@RG\tID:spike\tSM:spike\tPL:ILLUMINA\tLB:spike' "$REF" \
      "$d/${n}_R1.fastq.gz" "$d/${n}_R2.fastq.gz" 2> "$d/bwa.log" | samtools sort -o "$d/window.bam" -
    samtools index "$d/window.bam"; rm -f "$d"/*.fastq.gz
  fi
  analyse "$n" "$d/window.bam"
done
for n in SPIKE_NEG SPIKE_S1 SPIKE_S2 SPIKE_L1 SPIKE_L2 SPIKE_L3 SPIKE_L4 SPIKE_L5 SPIKE_L6; do
  [[ -s $SP/sim/$n/window.bam.bai ]] || samtools index "$SP/sim/$n/window.bam"
  analyse "$n" "$SP/sim/$n/window.bam"
done
