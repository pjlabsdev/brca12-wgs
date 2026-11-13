#!/usr/bin/env bash
# Real carriers of pathogenic BRCA1/BRCA2 variants from the 1000 Genomes Project.
#
# Samples: carriers_1000g.tsv (found in the NYGC 30x joint call set of the 3,202 samples, open data).
# Reads:   NYGC 30x CRAMs on the public 1000 Genomes S3 bucket; read pairs over BRCA1/BRCA2 +/- 10 kb
#          are extracted and turned back into FASTQ, so the pipeline runs from FASTQ (BWA-MEM included).
# Checks:  (1) is the known pathogenic variant called, heterozygous;
#          (2) concordance of all small-variant calls with the NYGC call set (PASS) for that sample.
# Stand-ins as in the GIAB benchmark: VEP (passthrough) and CNVkit (neutral segments).
#
# Usage: validation/pathogenic/run_carriers.sh WORKDIR   (WORKDIR/data must hold ref.fa, ref.sdf,
#        clinvar.vcf.gz, vep_cache and stand_ins as prepared by validation/run_giab.sh)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
HERE="$ROOT/validation/pathogenic"
WORK="${1:?usage: run_carriers.sh WORKDIR}"
DATA="$WORK/data"
KG="$WORK/kg"
THREADS="${THREADS:-4}"
S3=https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage
CALLS=$S3/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05
REG="chr17:43034295-43180327 chr13:32305077-32410268"
mkdir -p "$KG" && cd "$KG"

samples=$(grep -v '^#' "$HERE/carriers_1000g.tsv" | cut -f1 | paste -sd,)
if [[ ! -s nygc.vcf.gz ]]; then
  bcftools view -s "$samples" -r chr13:32315077-32400268 "${CALLS}_chr13.recalibrated_variants.vcf.gz" -Oz -o nygc13.vcf.gz
  bcftools view -s "$samples" -r chr17:43044295-43170327 "${CALLS}_chr17.recalibrated_variants.vcf.gz" -Oz -o nygc17.vcf.gz
  bcftools concat nygc13.vcf.gz nygc17.vcf.gz -Oz -o nygc.vcf.gz && tabix -f nygc.vcf.gz
fi
awk 'BEGIN{OFS="\t"} !/^#/ {print $1,$2,$3}' "$ROOT/assets/brca12_regions.bed" > genes.bed

grep -v '^#' "$HERE/carriers_1000g.tsv" | while IFS=$'\t' read -r s run prefix _; do
  out="$WORK/carriers/$s"
  if [[ ! -s ${s}_R1.fastq.gz ]]; then
    for attempt in 1 2 3 4; do   # the remote read occasionally stalls
      timeout 300 samtools view -@2 -b -T "$DATA/ref.fa" -o "$s.region.bam" "$S3/$prefix/$run/$s.final.cram" $REG && break
      [[ $attempt -lt 4 ]] || exit 1
    done
    samtools collate -@2 -u -O "$s.region.bam" | samtools fastq -@2 -1 "${s}_R1.fastq.gz" -2 "${s}_R2.fastq.gz" -0 /dev/null -s /dev/null -n -
  fi
  if [[ ! -s $out/brca12_results.json ]]; then
    (cd "$WORK" && nextflow run "$ROOT/main.nf" -profile standard -work-dir "$WORK/work_$s" \
       --sample "$s" --reads_r1 "$KG/${s}_R1.fastq.gz" --reads_r2 "$KG/${s}_R2.fastq.gz" \
       --reference "$DATA/ref.fa" --vep_cache "$DATA/vep_cache" --clinvar_vcf "$DATA/clinvar.vcf.gz" \
       --vep "$ROOT/validation/vep_passthrough.sh" --cnvkit "$DATA/stand_ins/cnvkit.py" \
       --threads "$THREADS" --outdir "$out" < /dev/null)
    rm -rf "$WORK/work_$s"
  fi
  # NYGC calls of this sample: PASS, non-reference, unused alleles removed
  bcftools view -s "$s" -f PASS,. nygc.vcf.gz | bcftools view -c 1 -a | bcftools norm -m -any -f "$DATA/ref.fa" \
    | bcftools view -e 'ALT="*"' -Oz -o "$s.nygc.vcf.gz" && tabix -f "$s.nygc.vcf.gz"
  python3 "$ROOT/validation/giab_benchmark.py" --sample "$s" --calls "$out/calls/$s.small.vcf.gz" \
    --truth "$s.nygc.vcf.gz" --truth-bed genes.bed \
    --regions "$ROOT/assets/brca12_regions.bed" --exons "$ROOT/assets/brca12_exons.bed" \
    --sdf "$DATA/ref.sdf" --outdir "$out/benchmark" --threads "$THREADS"
done
python3 "$HERE/carriers_eval.py" --table "$HERE/carriers_1000g.tsv" --workdir "$WORK/carriers" \
  --reference "$DATA/ref.fa" --out "$WORK/carriers/carriers_summary"
