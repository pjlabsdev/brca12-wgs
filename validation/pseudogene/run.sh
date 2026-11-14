#!/usr/bin/env bash
# BRCA1 exon 1-2 rearrangements through the pseudogene BRCA1P1, as SIMULATED spike-ins in the real
# HG002 background of ../pathogenic/run_spikeins.sh, run through the whole pipeline.
# Usage: run.sh WORKDIR   (WORKDIR as for ../pathogenic/run_spikeins.sh; needs cnvkit.py and Rscript on PATH)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(cd "$HERE/../.." && pwd)"
WORK="${1:?}"; SP="$WORK/spike"; REF="$WORK/data/ref.fa"; OUT="$WORK/pseudogene"; THREADS="${THREADS:-4}"
mkdir -p "$OUT" && cd "$OUT"
for s in SPIKE_PSI SPIKE_PSIDUP; do
  [[ -s out/$s/brca12_results.json ]] && continue
  python3 "$ROOT/validation/pathogenic/spike.py" simulate --sample "$s" --events "$HERE/events.tsv" \
    --fragments "$SP/fragments.tsv.gz" --reference "$REF" --outdir "sim/$s"
  bwa mem -t "$THREADS" -R '@RG\tID:spike\tSM:spike\tPL:ILLUMINA\tLB:spike' "$REF" \
    "sim/$s/${s}_R1.fastq.gz" "sim/$s/${s}_R2.fastq.gz" 2> "sim/$s/bwa.log" | samtools sort -@2 -o "sim/$s/window.bam" -
  samtools merge -@"$THREADS" -f "$s.bam" "$SP/bg_clean.bam" "sim/$s/window.bam" && samtools index "$s.bam"
  nextflow run "$ROOT/main.nf" -profile standard -work-dir "$OUT/work_$s" --sample "$s" --bam "$OUT/$s.bam" \
    --reference "$REF" --vep_cache "$WORK/data/vep_cache" --clinvar_vcf "$WORK/data/clinvar.vcf.gz" \
    --vep "$ROOT/validation/vep_passthrough.sh" --cnvkit "${CNVKIT:-cnvkit.py}" --threads "$THREADS" --outdir "$OUT/out/$s" < /dev/null
  rm -rf "$OUT/work_$s" "$s.bam" "$s.bam.bai"
done
