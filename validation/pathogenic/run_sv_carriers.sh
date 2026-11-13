#!/usr/bin/env bash
# Real structural variants in BRCA1/BRCA2 from the 1000 Genomes Project (not disease-causing).
#
# Sites:   sv_sites_1000g.tsv, two Alu insertions in BRCA2 and two deletions in BRCA1 from the phased
#          SNV/INDEL/SV panel of the 3,202 samples (short-read calls: GATK-SV/MELT and svtools).
# Samples: sv_carriers_1000g.tsv, carriers according to the panel plus six samples that are 0|0 at
#          all four sites. Reads as in run_carriers.sh (NYGC 30x CRAM, gene +/- 10 kb, back to FASTQ).
#          Mates that the original alignment placed outside the two windows are only included
#          with MATES=1 (output then goes to WORKDIR/sv_carriers_mates).
# Stand-ins: VEP (passthrough) and CNVkit (neutral segments).
#
# Usage: validation/pathogenic/run_sv_carriers.sh WORKDIR
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
HERE="$ROOT/validation/pathogenic"
WORK="${1:?usage: run_sv_carriers.sh WORKDIR}"
DATA="$WORK/data"
SUFFIX="${MATES:+_mates}"
KG="$WORK/svkg$SUFFIX"
FULLREF=https://1000genomes.s3.amazonaws.com/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
TABLE="${TABLE:-$HERE/sv_carriers_1000g.tsv}"
THREADS="${THREADS:-4}"
S3=https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage
REG="chr17:43034295-43180327 chr13:32305077-32410268"
mkdir -p "$KG" && cd "$KG"

grep -v '^#' "$TABLE" | while IFS=$'\t' read -r s run prefix _; do
  out="$WORK/sv_carriers$SUFFIX/$s"
  [[ -s $out/brca12_results.json ]] && continue
  if [[ ! -s ${s}_R1.fastq.gz ]]; then
    for attempt in 1 2 3 4; do   # the remote read occasionally stalls
      timeout 300 samtools view -@2 -b -T "$DATA/ref.fa" -o "$s.region.bam" "$S3/$prefix/$run/$s.final.cram" $REG && break
      [[ $attempt -lt 4 ]] || exit 1
    done
    if [[ ${MATES:-0} == 1 ]]; then
      # also fetch mates that the original alignment placed outside the windows (needs the full
      # reference to decode other contigs; slow, several minutes per sample)
      samtools view "$s.region.bam" | awk -v OFS='\t' '{c=($7=="="?$3:$7); out=($7!="=") || ($3=="chr13" && ($8<32305077 || $8>32410268)) || ($3=="chr17" && ($8<43034295 || $8>43180327)); if (out && $7!="*") {print c,$8-1,$8 > "'"$s"'.mates.bed"; print $1 > "'"$s"'.names.txt"}}'
      sort -k1,1 -k2,2n -u "$s.mates.bed" -o "$s.mates.bed"
      timeout 1500 samtools view -@2 -b -M -L "$s.mates.bed" -T "$FULLREF" -o "$s.matesall.bam" "$S3/$prefix/$run/$s.final.cram"
      samtools view -b -N "$s.names.txt" -o "$s.mates.bam" "$s.matesall.bam" && rm -f "$s.matesall.bam"
      samtools cat -o "$s.both.bam" "$s.region.bam" "$s.mates.bam" && mv "$s.both.bam" "$s.region.bam"
    fi
    samtools collate -@2 -u -O "$s.region.bam" | samtools fastq -@2 -1 "${s}_R1.fastq.gz" -2 "${s}_R2.fastq.gz" -0 /dev/null -s /dev/null -n -
  fi
  (cd "$WORK" && nextflow run "$ROOT/main.nf" -profile standard -work-dir "$WORK/work_$s" \
     --sample "$s" --reads_r1 "$KG/${s}_R1.fastq.gz" --reads_r2 "$KG/${s}_R2.fastq.gz" \
     --reference "$DATA/ref.fa" --vep_cache "$DATA/vep_cache" --clinvar_vcf "$DATA/clinvar.vcf.gz" \
     --vep "$ROOT/validation/vep_passthrough.sh" --cnvkit "$DATA/stand_ins/cnvkit.py" \
     --threads "$THREADS" --outdir "$out" < /dev/null)
  rm -rf "$WORK/work_$s"
done
python3 "$HERE/sv_carriers_eval.py" --sites "$HERE/sv_sites_1000g.tsv" --carriers "$TABLE" \
  --workdir "$WORK/sv_carriers$SUFFIX" --out "$WORK/sv_carriers$SUFFIX/sv_carriers_summary"
