#!/usr/bin/env bash
# Arm B in the cloud (deviation from validation/giab_more/ADDENDUM_arm_b.md, see README.md here):
# the same Google NovaSeq 6000 PCR-free 30x reads of HG003 to HG007 that arm B planned to align as
# whole genomes, but only the read pairs over BRCA1/BRCA2 +/- 10 kb. They are taken from Google's
# public GRCh38 BAMs of those reads (index range requests, a few MB per sample), turned back into
# FASTQ and run through the pipeline from FASTQ like arm A. Truth: GIAB v4.2.1 for all five.
#
# Real tools: bwa, samtools, GATK HaplotypeCaller, bcftools, mosdepth, Delly.
# Stand-ins:  VEP (passthrough) and CNVkit (neutral segments, with a stand-in R), as in arm A.
#
# Usage: CODE=/path/to/brca12-wgs-checkout validation/giab_novaseq/run.sh WORKDIR [REFDIR]
#   results go to WORKDIR/<code revision>/<sample>
#   CODE    pipeline checkout to test (default: this repository); arm B froze main a705662
#   REFDIR  folder with ref.fa (chr13 + chr17 of the GIAB GRCh38 no-alt set) + bwa/GATK/rtg indexes,
#           as built by validation/giab_more/run.sh (default: WORKDIR/data)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
CODE="$(cd "${CODE:-$ROOT}" && pwd)"
WORK="${1:?usage: run.sh WORKDIR [REFDIR]}"
mkdir -p "$WORK/data"
WORK="$(cd "$WORK" && pwd)"            # absolute before any cd, so relative arguments work
REFDIR="$(cd "${2:-$WORK/data}" && pwd)"
REV="$(git -C "$CODE" rev-parse --short HEAD 2>/dev/null || echo unknown)"
cd "$WORK/data"
G=https://storage.googleapis.com/brain-genomics-public/research/sequencing/grch38/bam/novaseq/wgs_pcr_free/30x
S=https://s3.amazonaws.com/giab/release
REG="chr17:43034295-43180327 chr13:32305077-32410268"
THREADS="${THREADS:-4}"
test -s "$REFDIR/ref.fa.bwt" || { echo "no reference index in $REFDIR (build it with validation/giab_more/run.sh)"; exit 1; }

declare -A TRUTH=(
  [HG003]=$S/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark
  [HG004]=$S/AshkenazimTrio/HG004_NA24143_mother/NISTv4.2.1/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark
  [HG005]=$S/ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh38/HG005_GRCh38_1_22_v4.2.1_benchmark
  [HG006]=$S/ChineseTrio/HG006_NA24694_father/NISTv4.2.1/GRCh38/HG006_GRCh38_1_22_v4.2.1_benchmark
  [HG007]=$S/ChineseTrio/HG007_NA24695_mother/NISTv4.2.1/GRCh38/HG007_GRCh38_1_22_v4.2.1_benchmark
)
# the Ashkenazi trio has a BED without the regions that are inconsistent in the trio, as in arm A
declare -A BEDSUFFIX=([HG003]=_noinconsistent.bed [HG004]=_noinconsistent.bed [HG005]=.bed [HG006]=.bed [HG007]=.bed)

[[ -s clinvar.vcf.gz.tbi ]] || {
  printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > clinvar.vcf
  bgzip -f clinvar.vcf && tabix -f -p vcf clinvar.vcf.gz; }
mkdir -p vep_cache stand_ins
ln -sf "$CODE/tests/fake_tools/faketool.py" stand_ins/cnvkit.py
ln -sf "$CODE/tests/fake_tools/faketool.py" stand_ins/Rscript   # CNVkit stand-in: R/DNAcopy check

for s in ${SAMPLES:-HG003 HG004 HG005 HG006 HG007}; do
  if [[ ! -s ${s}_R1.fastq.gz ]]; then
    samtools view -@2 -b -o "$s.region.bam" "$G/$s.novaseq.pcr-free.30x.dedup.grch38.bam" $REG
    samtools collate -@2 -u -O -T "$s.col" "$s.region.bam" \
      | samtools fastq -@2 -1 "${s}_R1.fastq.gz" -2 "${s}_R2.fastq.gz" -0 /dev/null -s /dev/null -n -
  fi
  [[ -s $s.truth.vcf.gz ]] || curl -sfo $s.truth.vcf.gz "${TRUTH[$s]}.vcf.gz"
  [[ -s $s.truth.vcf.gz.tbi ]] || curl -sfo $s.truth.vcf.gz.tbi "${TRUTH[$s]}.vcf.gz.tbi"
  [[ -s $s.bench.bed ]] || curl -sfo $s.bench.bed "${TRUTH[$s]}${BEDSUFFIX[$s]}"

  # results are kept per code revision, so one WORKDIR can hold runs of several versions
  out="$WORK/$REV/$s"
  [[ -s $out/benchmark/summary.json ]] && continue
  (cd "$WORK" && nextflow run "$CODE/main.nf" -profile standard -work-dir "$WORK/work_${REV}_$s" \
     --sample "$s" --reads_r1 "$WORK/data/${s}_R1.fastq.gz" --reads_r2 "$WORK/data/${s}_R2.fastq.gz" \
     --reference "$REFDIR/ref.fa" --vep_cache "$WORK/data/vep_cache" --clinvar_vcf "$WORK/data/clinvar.vcf.gz" \
     --vep "$CODE/validation/vep_passthrough.sh" --cnvkit "$WORK/data/stand_ins/cnvkit.py" --rscript "$WORK/data/stand_ins/Rscript" \
     --threads "$THREADS" --outdir "$out")
  python3 "$CODE/validation/giab_benchmark.py" --sample "$s" --calls "$out/calls/$s.small.vcf.gz" \
    --truth "$WORK/data/$s.truth.vcf.gz" --truth-bed "$WORK/data/$s.bench.bed" \
    --regions "$CODE/assets/brca12_regions.bed" --exons "$CODE/assets/brca12_exons.bed" \
    --sdf "$REFDIR/ref.sdf" --outdir "$out/benchmark" --threads "$THREADS"
  rm -rf "$WORK/work_${REV}_$s"
done
