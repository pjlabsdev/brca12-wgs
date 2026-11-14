#!/usr/bin/env bash
# Exome and amplicon benchmark of brca12-wgs (--assay exome / amplicon) against GIAB v4.2.1.
#
# Three public GIAB data sets, all aligned to hg19/GRCh37 by their producers. Reads over
# BRCA1/BRCA2 +/- 10 kb are extracted, turned back into FASTQ and run through the pipeline from
# FASTQ (BWA-MEM to GRCh38 included):
#   HG002_exome_agilent  Oslo University Hospital, Agilent SureSelect v5, HiSeq 2x125 (paired)
#   HG001_exome_truseq   University of Nebraska, Illumina TruSeq exome (paired), kit BED public
#   HG002_ampliseq_ion   Thermo Fisher Ion AmpliSeq Exome, Ion Proton (single-end amplicons),
#                        one of four runs; the amplicon "effective regions" BED is public
# Targets: the vendor BED where it is public (TruSeq, AmpliSeq), otherwise the MANE exons +/- 50 bp.
# hg19 target intervals are moved to GRCh38 with the single offset per locus (Ensembl assembly
# mapping: BRCA1 +1,847,983, BRCA2 -574,137; both loci map as one block).
# Truth: GIAB v4.2.1 benchmark VCF/BED (GRCh38), as in validation/run_giab.sh.
#
# Real tools: bwa, samtools, GATK, bcftools, mosdepth, CNVkit, Manta. VEP is a passthrough (no cache).
# Usage: validation/targeted/run_targeted.sh WORKDIR GIAB_DATA_DIR
#   GIAB_DATA_DIR holds ref.fa (+ BWA index, .dict), ref.sdf and HG00{1,2}.truth.vcf.gz/.bench.bed
#   from validation/run_giab.sh.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
WORK="${1:?usage: run_targeted.sh WORKDIR GIAB_DATA_DIR}"
DATA="$(cd "${2:?usage: run_targeted.sh WORKDIR GIAB_DATA_DIR}" && pwd)"
THREADS="${THREADS:-2}"
S=https://giab.s3.amazonaws.com
mkdir -p "$WORK/input" && cd "$WORK/input"

# hg19 windows: BRCA1 GRCh38 chr17:43034295-43180327, BRCA2 chr13:32305077-32410268 (+/- 10 kb)
REG19_NUM="17:41186312-41332344 13:32879214-32984405"
REG19_CHR="chr17:41186312-41332344 chr13:32879214-32984405"

lift_bed() {  # hg19 BED (chr or numeric) -> GRCh38 BED restricted to the BRCA loci
  awk 'BEGIN{OFS="\t"} !/^(#|track|browser)/ {
         c=$1; sub(/^chr/,"",c)
         if (c=="17" && $3>41186312 && $2<41332344) { off=1847983 }
         else if (c=="13" && $3>32879214 && $2<32984405) { off=-574137 }
         else next
         name = (NF>=4 ? $4 : ".")
         print "chr" c, $2+off, $3+off, name }' "$1" | sort -k1,1 -k2,2n
}

paired_fastq() {  # $1 sample, $2 BAM URL, $3 regions
  [[ -s ${1}_R1.fastq.gz ]] && return
  samtools view -@2 -b -F 0x900 -o "$1.region.bam" "$2" $3
  samtools collate -@2 -u -O "$1.region.bam" | samtools fastq -@2 -1 "${1}_R1.fastq.gz" -2 "${1}_R2.fastq.gz" -0 /dev/null -s /dev/null -n -
  rm -f "$1.region.bam"
}

# 1. HG002 Agilent SureSelect v5 exome (no public kit BED: MANE exons +/- 50 bp as targets)
paired_fastq HG002_exome_agilent \
  "$S/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam" \
  "$REG19_NUM"
[[ -s HG002_exome_agilent.targets.bed ]] || awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $2-50, $3+50, $4}' "$ROOT/assets/brca12_exons.bed" > HG002_exome_agilent.targets.bed

# 2. HG001 Illumina TruSeq exome with the vendor target BED
paired_fastq HG001_exome_truseq "$S/data/NA12878/Nebraska_NA12878_HG001_TruSeq_Exome/NIST-hg001-7001-ready.bam" "$REG19_NUM"
if [[ ! -s HG001_exome_truseq.targets.bed ]]; then
  curl -sf "$S/data/NA12878/Nebraska_NA12878_HG001_TruSeq_Exome/TruSeq_exome_targeted_regions.hg19.bed" -o truseq.hg19.bed
  lift_bed truseq.hg19.bed > HG001_exome_truseq.targets.bed
fi

# 3. HG002 Ion AmpliSeq Exome, single-end amplicon reads (primers not public: no primer clipping)
if [[ ! -s HG002_ampliseq_ion_R1.fastq.gz ]]; then
  samtools view -@2 -b -F 0x900 -o ion.region.bam \
    "$S/data/AshkenazimTrio/HG002_NA24385_son/ion_exome/HG002_NA24385_SRR1767406_IonXpress_020_rawlib_24028.bam" $REG19_CHR
  samtools fastq -@2 ion.region.bam | gzip -c > HG002_ampliseq_ion_R1.fastq.gz
  rm -f ion.region.bam
fi
if [[ ! -s HG002_ampliseq_ion.targets.bed ]]; then
  curl -sf "$S/data/AshkenazimTrio/analysis/IonTorrent_TVC_03162015/AmpliseqExome.20141120_effective_regions.bed" -o ampliseq.hg19.bed
  lift_bed ampliseq.hg19.bed > HG002_ampliseq_ion.targets.bed
fi

: > clinvar.vcf && printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > clinvar.vcf
bgzip -f clinvar.vcf && tabix -f -p vcf clinvar.vcf.gz
mkdir -p vep_cache

run() {  # $1 name, $2 truth sample, $3 assay, $4.. extra pipeline arguments
  local name=$1 truth=$2 assay=$3; shift 3
  local out="$WORK/$name"
  if [[ ! -s $out/brca12_results.json ]]; then
    (cd "$WORK" && nextflow run "$ROOT/main.nf" -profile standard -work-dir "$WORK/work_$name" \
       --sample "$name" --assay "$assay" --targets_bed "$WORK/input/$name.targets.bed" --mosaic "$@" \
       --reference "$DATA/ref.fa" --vep_cache "$WORK/input/vep_cache" --clinvar_vcf "$WORK/input/clinvar.vcf.gz" \
       --vep "$ROOT/validation/vep_passthrough.sh" --threads "$THREADS" --outdir "$out")
  fi
  python3 "$ROOT/validation/giab_benchmark.py" --sample "$name" --calls "$out/calls/$name.small.vcf.gz" \
    --truth "$DATA/$truth.truth.vcf.gz" --truth-bed "$DATA/$truth.bench.bed" \
    --regions "$ROOT/assets/brca12_regions.bed" --exons "$ROOT/assets/brca12_exons.bed" \
    --targets "$WORK/input/$name.targets.bed" \
    --sdf "$DATA/ref.sdf" --outdir "$out/benchmark" --threads "$THREADS"
}

run HG002_exome_agilent HG002 exome \
  --reads_r1 "$WORK/input/HG002_exome_agilent_R1.fastq.gz" --reads_r2 "$WORK/input/HG002_exome_agilent_R2.fastq.gz"
run HG001_exome_truseq HG001 exome \
  --reads_r1 "$WORK/input/HG001_exome_truseq_R1.fastq.gz" --reads_r2 "$WORK/input/HG001_exome_truseq_R2.fastq.gz"
run HG002_ampliseq_ion HG002 amplicon --platform iontorrent \
  --reads_r1 "$WORK/input/HG002_ampliseq_ion_R1.fastq.gz"

python3 - "$WORK" <<'PY'
import json, os, sys
work = sys.argv[1]
rows = []
for name in ("HG002_exome_agilent", "HG001_exome_truseq", "HG002_ampliseq_ion"):
    bench = json.load(open(os.path.join(work, name, "benchmark", "summary.json")))
    qc = json.load(open(os.path.join(work, name, "brca12_results.json"))).get("readQC") or {}
    rows.append({"dataset": name, "benchmark": bench["results"], "targetSummary": qc.get("targetSummary"),
                 "regionStats": qc.get("regionStats"), "flags": qc.get("flags")})
json.dump(rows, open(os.path.join(work, "targeted_summary.json"), "w"), indent=1)
print(json.dumps([{"dataset": r["dataset"], **{f"{b['region']}:{b['type']}": f"{b['TP']}/{b['truthTotal']} FP {b['FP']}"
                   for b in r["benchmark"]}} for r in rows], indent=1))
PY
