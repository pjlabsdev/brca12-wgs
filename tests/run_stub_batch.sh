#!/usr/bin/env bash
# Wiring test for the multi-sample run: two samples from a sample sheet through bin/brca12_batch.py
# with stand-in tools (tests/fake_tools), then the cohort summary. Does NOT test the real tools.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
FAKE="$WORK/fake"; mkdir -p "$FAKE"
for tool in bwa samtools gatk bcftools tabix bgzip vep configManta.py cnvkit.py mosdepth; do
  ln -s "$ROOT/tests/fake_tools/faketool.py" "$FAKE/$tool"
done
cd "$WORK"
printf ">chr13\nN\n>chr17\nN\n" > ref.fa
printf "chr13\t114364328\t0\t60\t61\nchr17\t83257441\t0\t60\t61\n" > ref.fa.fai
touch ref.fa.amb ref.fa.ann ref.fa.bwt ref.fa.pac ref.fa.sa ref.dict R1.fq.gz R2.fq.gz clinvar.vcf.gz clinvar.vcf.gz.tbi
mkdir -p vep
printf "sample,reads_r1,reads_r2,bam\nS01,%s,%s,\nS02,%s,%s,\n" "$WORK/R1.fq.gz" "$WORK/R2.fq.gz" "$WORK/R1.fq.gz" "$WORK/R2.fq.gz" > samples.csv
PATH="$FAKE:$PATH" python3 "$ROOT/bin/brca12_batch.py" --samplesheet samples.csv --outdir "$WORK/batch" --run-name STUB --jobs 2 -- \
  -profile standard --reference "$WORK/ref.fa" --vep_cache "$WORK/vep" --clinvar_vcf "$WORK/clinvar.vcf.gz" \
  --exons_bed "$ROOT/tests/data/exons.test.bed" --threads 1
for f in S01/brca12_results.json S02/brca12_report.xlsx cohort/cohort_samples.tsv cohort/cohort_variants.tsv cohort/cohort_contamination.tsv cohort/cohort.xlsx; do
  test -s "$WORK/batch/$f" || { echo "missing $f"; exit 1; }
done
test "$(tail -n +2 "$WORK/batch/cohort/cohort_samples.tsv" | wc -l | tr -d ' ')" = 2
echo "stub batch OK"
