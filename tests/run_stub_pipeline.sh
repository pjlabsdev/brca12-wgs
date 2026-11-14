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
for f in brca12_report.html brca12_report.md brca12_results.json brca12_commands.sh; do
  test -s "$WORK/out/$f" || { echo "missing $f"; exit 1; }
done
python3 -c "import json,sys; r=json.load(open(sys.argv[1])); assert r['schema']=='genome.brca12.v1'; assert r['genes']['BRCA1']['smallVariants']; assert r['readQC']['regionStats']['readsInRegions']==1000; print('stub pipeline (WGS) OK')" "$WORK/out/brca12_results.json"

# Single-end amplicon panel with primer clipping, target QC and mosaic candidates.
printf "chr17\t43044900\t43045200\tBRCA1_amp1\nchr13\t32315400\t32315700\tBRCA2_amp1\n" > targets.bed
printf "chr17\t43044900\t43044922\tBRCA1_amp1_F\t0\t+\nchr17\t43045178\t43045200\tBRCA1_amp1_R\t0\t-\n" > primers.bed
PATH="$FAKE:$PATH" nextflow run "$ROOT/main.nf" -profile standard \
  --sample AMP --reference "$WORK/ref.fa" --reads_r1 "$WORK/R1.fq.gz" \
  --vep_cache "$WORK/vep" --clinvar_vcf "$WORK/clinvar.vcf.gz" \
  --exons_bed "$ROOT/tests/data/exons.test.bed" --threads 2 --outdir "$WORK/out_amp" \
  --assay amplicon --targets_bed "$WORK/targets.bed" --primers_bed "$WORK/primers.bed" --mosaic
python3 - "$WORK/out_amp/brca12_results.json" <<'PY'
import json, sys
r = json.load(open(sys.argv[1]))
assert r["assay"] == "amplicon"
qc = r["readQC"]
assert qc["primerClipping"]["clippedPct"] == 93.9, qc["primerClipping"]
assert qc["targetSummary"]["targetCount"] == 2 and "uniformityPct" in qc["targetSummary"], qc["targetSummary"]
mosaic = r["mosaicCandidates"]
assert [m["pos"] for m in mosaic] == [43094000], mosaic
assert any("Amplikon" in item for item in r["limitations"])
assert r["input"]["singleEnd"] and r["svs"] == []
print("stub pipeline (amplicon, mosaic) OK")
PY
