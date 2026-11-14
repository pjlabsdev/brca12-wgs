#!/usr/bin/env bash
# End-to-end wiring test: runs the whole Nextflow pipeline with stand-in tools
# (tests/fake_tools) on a dummy GRCh38 index. Checks that every process connects
# and that the report files are produced. It does NOT test the real bioinformatics tools.
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
PATH="$FAKE:$PATH" nextflow run "$ROOT/main.nf" -profile standard \
  --sample TEST --reference "$WORK/ref.fa" --reads_r1 "$WORK/R1.fq.gz" --reads_r2 "$WORK/R2.fq.gz" \
  --vep_cache "$WORK/vep" --clinvar_vcf "$WORK/clinvar.vcf.gz" \
  --exons_bed "$ROOT/tests/data/exons.test.bed" --threads 2 --outdir "$WORK/out"
for f in brca12_report.html brca12_report.md brca12_results.json brca12_commands.sh brca12_report.xlsx; do
  test -s "$WORK/out/$f" || { echo "missing $f"; exit 1; }
done
python3 -c "import json,sys; r=json.load(open(sys.argv[1])); assert r['schema']=='genome.brca12.v1'; assert r['genes']['BRCA1']['smallVariants']; print('stub pipeline OK')" "$WORK/out/brca12_results.json"
