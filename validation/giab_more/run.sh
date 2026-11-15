#!/usr/bin/env bash
# GIAB benchmark of brca12-wgs for HG003 to HG007 (GRCh38), frozen in PREREGISTRATION.md.
# Same method as validation/run_giab.sh (HG001, HG002); only the samples differ.
#
# Reads: GIAB/NHGRI Illumina HiSeq 2x148 bp novoalign BAMs on the public GIAB S3 mirror; read pairs
# over BRCA1/BRCA2 +/- 10 kb are extracted, downsampled to ~30-40x and turned back into FASTQ, so the
# pipeline runs from FASTQ (BWA-MEM alignment included). Reference: chr13 + chr17 of the GIAB
# GRCh38 no-alt analysis set. Truth: GIAB v4.2.1 benchmark VCF and high-confidence BED.
#
# Real tools: bwa, samtools, GATK HaplotypeCaller, bcftools, mosdepth, Delly (pipeline default).
# Stand-ins:  VEP (passthrough, no cache available) and CNVkit (neutral segments); annotation and
#             CNV calling are therefore NOT part of this benchmark.
#
# Usage: validation/run_giab.sh WORKDIR   (needs nextflow, bwa, samtools, gatk, bcftools, tabix,
#        bgzip, mosdepth, rtg, delly on PATH; ~3 GB disk)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
WORK="${1:?usage: run.sh WORKDIR}"
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

A=$S/data/AshkenazimTrio
C=$S/data/ChineseTrio
declare -A BAM=(
  [HG003]=$A/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG003.GRCh38.60x.1.bam
  [HG004]=$A/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG004.GRCh38.60x.1.bam
  [HG005]=$C/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/NHGRI_Illumina300X_Chinesetrio_novoalign_bams/HG005.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.300x.bam
  [HG006]=$C/HG006_NA24694-huCA017E_father/NA24694_Father_HiSeq100x/NHGRI_Illumina100X_Chinesetrio_novoalign_bams/HG006.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.100x.bam
  [HG007]=$C/HG007_NA24695-hu38168_mother/NA24695_Mother_HiSeq100x/NHGRI_Illumina100X_Chinesetrio_novoalign_bams/HG007.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.100x.bam
)
# samtools view -s SEED.FRACTION to about 30x
declare -A FRAC=([HG003]=42.5 [HG004]=42.5 [HG005]=42.1 [HG006]=42.3 [HG007]=42.3)
R=$S/release
declare -A TRUTHVCF=(
  [HG003]=$R/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
  [HG004]=$R/AshkenazimTrio/HG004_NA24143_mother/NISTv4.2.1/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
  [HG005]=$R/ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh38/HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
  [HG006]=$R/ChineseTrio/HG006_NA24694_father/NISTv3.3.2/GRCh38/HG006_GIAB_GRCh38_highconf_CG-IllFB-IllSNT-10X_CHROM1-22_v.3.3.2_highconf.vcf.gz
  [HG007]=$R/ChineseTrio/HG007_NA24695_mother/NISTv3.3.2/GRCh38/HG007_GIAB_GRCh38_highconf_CG-IllFB-IllSNT-10X_CHROM1-22_v.3.3.2_highconf.vcf.gz
)
declare -A TRUTHBED=(
  [HG003]=$R/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
  [HG004]=$R/AshkenazimTrio/HG004_NA24143_mother/NISTv4.2.1/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
  [HG005]=$R/ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh38/HG005_GRCh38_1_22_v4.2.1_benchmark.bed
  [HG006]=$R/ChineseTrio/HG006_NA24694_father/NISTv3.3.2/GRCh38/HG006_GIAB_GRCh38_highconf_CG-IllFB-IllSNT-10X_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed
  [HG007]=$R/ChineseTrio/HG007_NA24695_mother/NISTv3.3.2/GRCh38/HG007_GIAB_GRCh38_highconf_CG-IllFB-IllSNT-10X_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed
)

: > clinvar.vcf && printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > clinvar.vcf
bgzip -f clinvar.vcf && tabix -f -p vcf clinvar.vcf.gz
mkdir -p vep_cache stand_ins
ln -sf "$ROOT/tests/fake_tools/faketool.py" stand_ins/cnvkit.py

for s in ${SAMPLES:-HG003 HG004 HG005 HG006 HG007}; do
  if [[ ! -s ${s}_R1.fastq.gz ]]; then
    samtools view -@2 -b -s "${FRAC[$s]}" -o "$s.region.bam" "${BAM[$s]}" $REG
    samtools collate -@2 -u -O "$s.region.bam" | samtools fastq -@2 -1 "${s}_R1.fastq.gz" -2 "${s}_R2.fastq.gz" -0 /dev/null -s /dev/null -n -
  fi
  [[ -s $s.truth.vcf.gz ]] || curl -sfo $s.truth.vcf.gz "${TRUTHVCF[$s]}"
  [[ -s $s.truth.vcf.gz.tbi ]] || curl -sfo $s.truth.vcf.gz.tbi "${TRUTHVCF[$s]}.tbi"
  [[ -s $s.bench.bed ]] || curl -sfo $s.bench.bed "${TRUTHBED[$s]}"

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
