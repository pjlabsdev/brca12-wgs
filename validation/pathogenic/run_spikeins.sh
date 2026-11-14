#!/usr/bin/env bash
# Spike-in validation: synthetic pathogenic-type variants in real HG002 background reads.
#
# Background: GIAB/NHGRI HG002 HiSeq 2x148 bp, chr13 + chr17 at ~30x (so that CNVkit has a real
# chromosome-wide baseline). Inside BRCA1/BRCA2 +/- 20 kb the real pairs are replaced by pairs with
# the same positions, fragment sizes and base qualities, read from a reference haplotype and a
# haplotype carrying the spiked events (see spike.py). The replaced pairs are aligned with BWA-MEM,
# merged with the untouched background and the pipeline is run with --bam.
# Real tools: bwa, samtools, GATK HaplotypeCaller, bcftools, mosdepth, Manta, CNVkit. Stand-in: VEP.
#
# Usage: validation/pathogenic/run_spikeins.sh WORKDIR [SAMPLE ...]
#        (WORKDIR/data as prepared by validation/run_giab.sh; needs cnvkit.py and Rscript on PATH)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
HERE="$ROOT/validation/pathogenic"
WORK="${1:?usage: run_spikeins.sh WORKDIR [SAMPLE ...]}"; shift
DATA="$WORK/data"
SP="$WORK/spike"
THREADS="${THREADS:-4}"
CNVKIT="${CNVKIT:-cnvkit.py}"
S=https://s3.amazonaws.com/giab
BAM=$S/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.60x.1.bam
SAMPLES=("$@")
[[ ${#SAMPLES[@]} -gt 0 ]] || SAMPLES=(SPIKE_NEG SPIKE_S1 SPIKE_S2 SPIKE_L1 SPIKE_L2 SPIKE_L3 SPIKE_L4 SPIKE_L5 SPIKE_L6)
mkdir -p "$SP" && cd "$SP"

[[ -s bg.bam ]] || { samtools view -@3 -b -s 42.5 -o bg.bam "$BAM" chr13 chr17 && samtools index -@3 bg.bam; }
[[ -s bg_clean.bam.bai ]] || python3 "$HERE/spike.py" prepare --bam bg.bam --fragments fragments.tsv.gz --out bg_clean.bam
awk 'BEGIN{OFS="\t"} !/^#/ {print $1,$2,$3}' "$ROOT/assets/brca12_regions.bed" > genes.bed

for s in "${SAMPLES[@]}"; do
  out="$SP/out/$s"
  if [[ ! -s $out/brca12_results.json ]]; then
    python3 "$HERE/spike.py" simulate --sample "$s" --events "$HERE/small_events.tsv" "$HERE/large_events.tsv" \
      --fragments fragments.tsv.gz --reference "$DATA/ref.fa" --outdir "sim/$s"
    bwa mem -t "$THREADS" -R '@RG\tID:spike\tSM:spike\tPL:ILLUMINA\tLB:spike' "$DATA/ref.fa" \
      "sim/$s/${s}_R1.fastq.gz" "sim/$s/${s}_R2.fastq.gz" 2> "sim/$s/bwa.log" | samtools sort -@2 -o "sim/$s/window.bam" -
    samtools merge -@"$THREADS" -f "$s.bam" bg_clean.bam "sim/$s/window.bam" && samtools index -@"$THREADS" "$s.bam"
    (cd "$SP" && nextflow run "$ROOT/main.nf" -profile standard -work-dir "$SP/work_$s" \
       --sample "$s" --bam "$SP/$s.bam" \
       --reference "$DATA/ref.fa" --vep_cache "$DATA/vep_cache" --clinvar_vcf "$DATA/clinvar.vcf.gz" \
       --vep "$ROOT/validation/vep_passthrough.sh" --cnvkit "$CNVKIT" \
       --threads "$THREADS" --outdir "$out" < /dev/null)
    rm -rf "$SP/work_$s" "$s.bam" "$s.bam.bai"
  fi
  bgzip -f -k "sim/$s/$s.truth.vcf" && tabix -f "sim/$s/$s.truth.vcf.gz"
  python3 "$ROOT/validation/giab_benchmark.py" --sample "$s" --calls "$out/calls/$s.small.vcf.gz" \
    --truth "sim/$s/$s.truth.vcf.gz" --truth-bed genes.bed \
    --regions "$ROOT/assets/brca12_regions.bed" --exons "$ROOT/assets/brca12_exons.bed" \
    --sdf "$DATA/ref.sdf" --outdir "$out/benchmark" --threads "$THREADS" > /dev/null
done
python3 "$HERE/spike_eval.py" --simdir "$SP/sim" --outdir "$SP/out" --small "$HERE/small_events.tsv" --out "$SP/spike_summary"
