#!/usr/bin/env bash
# Mosaic (low allele fraction) validation of --mosaic: the 94 small events of SPIKE_S1 (SNVs and
# indels up to 40 bp, 5 modelled on known pathogenic variants) are spiked into real HG002 reads at
# allele fractions of 5, 10, 15, 20 and 30 % (share of read pairs drawn from the event haplotype),
# plus a negative control without events. The window reads are run from FASTQ against chr13+chr17
# with GATK HaplotypeCaller and, with --mosaic, GATK Mutect2. The background sample's own variants
# are not carried into the windows, so every call that is not a spiked event is a false positive.
# CNVkit is a stand-in here (the reads cover only the gene windows); VEP is a passthrough.
#
# Usage: validation/pathogenic/run_mosaic.sh WORKDIR   (WORKDIR as for run_spikeins.sh, with
#        WORKDIR/spike/fragments.tsv.gz from 'spike.py prepare')
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
HERE="$ROOT/validation/pathogenic"
WORK="${1:?usage: run_mosaic.sh WORKDIR}"
DATA="$WORK/data"
SP="$WORK/spike"
MO="$WORK/mosaic"
THREADS="${THREADS:-2}"
mkdir -p "$MO" && cd "$MO"
: > clinvar.vcf && printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > clinvar.vcf
bgzip -f clinvar.vcf && tabix -f -p vcf clinvar.vcf.gz
mkdir -p vep_cache stand_ins && ln -sf "$ROOT/tests/fake_tools/faketool.py" stand_ins/cnvkit.py

for af in 05 10 15 20 30 NEG; do
  s="MOSAIC_$af"
  out="$MO/out/$s"
  if [[ ! -s $out/brca12_results.json ]]; then
    if [[ $af == NEG ]]; then set_args=(--event-set SPIKE_NEG --fraction 0.5); else set_args=(--event-set SPIKE_S1 --fraction "0.$af"); fi
    python3 "$HERE/spike.py" simulate --sample "$s" "${set_args[@]}" --events "$HERE/small_events.tsv" \
      --fragments "$SP/fragments.tsv.gz" --reference "$DATA/ref.fa" --outdir "sim/$s"
    nextflow run "$ROOT/main.nf" -profile standard -work-dir "$MO/work_$s" \
      --sample "$s" --reads_r1 "$MO/sim/$s/${s}_R1.fastq.gz" --reads_r2 "$MO/sim/$s/${s}_R2.fastq.gz" --mosaic \
      --reference "$DATA/ref.fa" --vep_cache "$MO/vep_cache" --clinvar_vcf "$MO/clinvar.vcf.gz" \
      --vep "$ROOT/validation/vep_passthrough.sh" --cnvkit "$MO/stand_ins/cnvkit.py" \
      --threads "$THREADS" --outdir "$out"
    rm -rf "$MO/work_$s"
  fi
done
python3 "$HERE/mosaic_eval.py" --workdir "$MO" --out "$MO/mosaic_summary"
