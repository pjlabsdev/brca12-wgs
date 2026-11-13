#!/usr/bin/env bash
# GIAB benchmark of brca12-wgs for HG001 and HG002 (GRCh38, NIST v4.2.1).
#
# Reads: GIAB/NHGRI Illumina HiSeq 2x148 bp novoalign BAMs on the public GIAB S3 mirror; read pairs
# over BRCA1/BRCA2 +/- 10 kb are extracted, downsampled to ~30-40x and turned back into FASTQ, so the
# pipeline runs from FASTQ (BWA-MEM alignment included). Reference: chr13 + chr17 of the GIAB
# GRCh38 no-alt analysis set. Truth: GIAB v4.2.1 benchmark VCF and high-confidence BED.
#
# Real tools: bwa, samtools, GATK HaplotypeCaller, bcftools, mosdepth, Manta.
# Stand-ins:  VEP (passthrough, no cache available) and CNVkit (neutral segments); annotation and
#             CNV calling are therefore NOT part of this benchmark.
#
# Usage: validation/run_giab.sh WORKDIR   (needs nextflow, bwa, samtools, gatk, bcftools, tabix,
#        bgzip, mosdepth, rtg, configManta.py on PATH; ~2 GB disk)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
WORK="${1:?usage: run_giab.sh WORKDIR}"
mkdir -p "$WORK/data" && cd "$WORK/data"
S=https://s3.amazonaws.com/giab
REF=$S/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz
REG="chr17:43034295-43180327 chr13:32305077-32410268"
THREADS="${THREADS:-4}"

if [[ ! -s ref.fa.bwt ]]; then
  samtools faidx "$REF" chr13 chr17 > ref.fa
  samtools faidx ref.fa
  bwa index ref.fa
  gatk CreateSequenceDictionary -R ref.fa
  rtg format -o ref.sdf ref.fa
fi

declare -A BAM=(
  [HG001]=$S/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/NHGRI_Illumina300X_novoalign_bams/HG001.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.300x.bam
  [HG002]=$S/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.60x.1.bam
)
# samtools view -s SEED.FRACTION: HG001 300x -> ~36x, HG002 60x -> ~30x
declare -A FRAC=([HG001]=42.12 [HG002]=42.5)
declare -A TRUTHDIR=([HG001]=NA12878_HG001 [HG002]=AshkenazimTrio/HG002_NA24385_son)
declare -A TRUTHBED=([HG001]=HG001_GRCh38_1_22_v4.2.1_benchmark.bed [HG002]=HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed)

: > clinvar.vcf && printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > clinvar.vcf
bgzip -f clinvar.vcf && tabix -f -p vcf clinvar.vcf.gz
mkdir -p vep_cache stand_ins
ln -sf "$ROOT/tests/fake_tools/faketool.py" stand_ins/cnvkit.py

for s in HG001 HG002; do
  if [[ ! -s ${s}_R1.fastq.gz ]]; then
    samtools view -@2 -b -s "${FRAC[$s]}" -o "$s.region.bam" "${BAM[$s]}" $REG
    samtools collate -@2 -u -O "$s.region.bam" | samtools fastq -@2 -1 "${s}_R1.fastq.gz" -2 "${s}_R2.fastq.gz" -0 /dev/null -s /dev/null -n -
  fi
  base=$S/release/${TRUTHDIR[$s]}/NISTv4.2.1/GRCh38
  [[ -s $s.truth.vcf.gz ]] || curl -sfo $s.truth.vcf.gz $base/${s}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
  [[ -s $s.truth.vcf.gz.tbi ]] || curl -sfo $s.truth.vcf.gz.tbi $base/${s}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
  [[ -s $s.bench.bed ]] || curl -sfo $s.bench.bed $base/${TRUTHBED[$s]}

  out="$WORK/$s"
  (cd "$WORK" && nextflow run "$ROOT/main.nf" -profile standard -work-dir "$WORK/work_$s" \
     --sample "$s" --reads_r1 "$WORK/data/${s}_R1.fastq.gz" --reads_r2 "$WORK/data/${s}_R2.fastq.gz" \
     --reference "$WORK/data/ref.fa" --vep_cache "$WORK/data/vep_cache" --clinvar_vcf "$WORK/data/clinvar.vcf.gz" \
     --vep "$ROOT/validation/vep_passthrough.sh" --cnvkit "$WORK/data/stand_ins/cnvkit.py" \
     --threads "$THREADS" --outdir "$out")
  python3 "$ROOT/validation/giab_benchmark.py" --sample "$s" --calls "$out/calls/$s.small.vcf.gz" \
    --truth "$WORK/data/$s.truth.vcf.gz" --truth-bed "$WORK/data/$s.bench.bed" \
    --regions "$ROOT/assets/brca12_regions.bed" --exons "$ROOT/assets/brca12_exons.bed" \
    --sdf "$WORK/data/ref.sdf" --outdir "$out/benchmark" --threads "$THREADS"
done
