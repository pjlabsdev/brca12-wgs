#!/usr/bin/env bash
# End-to-end wiring test: runs the whole Nextflow pipeline with stand-in tools
# (tests/fake_tools) on a dummy GRCh38 index. Checks that every process connects
# and that the report files are produced. It does NOT test the real bioinformatics tools.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
FAKE="$WORK/fake"; mkdir -p "$FAKE"
for tool in bwa samtools gatk bcftools tabix bgzip vep configManta.py delly cnvkit.py Rscript mosdepth; do
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
python3 -c "import json,sys; r=json.load(open(sys.argv[1])); assert r['schema']=='genome.brca12.v1'; assert r['genes']['BRCA1']['smallVariants']; assert r['readQC']['regionStats']['readsInRegions']==1000; assert 'exonCopyNumber' in r and 'exonCopyNumberCallCount' in r['summary']; assert any(s['caller']=='Delly' for s in r['svs']), r['svs']; assert r['toolVersions'].get('delly')=='v1.3.3', r['toolVersions']; assert r['toolVersions'].get('dnacopy')=='1.76.0', r['toolVersions']; print('stub pipeline (WGS) OK')" "$WORK/out/brca12_results.json"

# The coverage step must read only the analysed regions through the index, not the whole alignment.
grep -q -- "view -b -M -L coverage_regions.bed" "$WORK/out/coverage/coverage.command.sh" || { echo "coverage step does not fetch regions through the index"; exit 1; }
if grep -E "mosdepth .*--by [^ ]*exons" "$WORK/out/coverage/coverage.command.sh" | grep -v -q "coverage_regions.bam"; then echo "mosdepth --by still reads the whole alignment"; exit 1; fi
echo "stub pipeline (coverage reads regions only) OK"

# An R without DNAcopy must stop the run with a clear message instead of a silent CNVkit result
# without segments (CNVkit batch -p does not pass on errors of its worker processes).
if FAKE_NO_DNACOPY=1 PATH="$FAKE:$PATH" nextflow run "$ROOT/main.nf" -profile standard \
  --sample NODNA --reference "$WORK/ref.fa" --reads_r1 "$WORK/R1.fq.gz" --reads_r2 "$WORK/R2.fq.gz" \
  --vep_cache "$WORK/vep" --clinvar_vcf "$WORK/clinvar.vcf.gz" \
  --exons_bed "$ROOT/tests/data/exons.test.bed" --threads 2 --outdir "$WORK/out_nodna" > nodna.log 2>&1; then
  echo "run without DNAcopy did not fail"; exit 1
fi
grep -q "CNVkit needs R with the Bioconductor package DNAcopy" nodna.log || { tail -40 nodna.log; echo "DNAcopy message missing"; exit 1; }
echo "stub pipeline (R without DNAcopy stops with a clear message) OK"

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

# Aligned input as CRAM with a .crai index (regression: the CRAM itself was taken as its own index).
printf 'cram' > "$WORK/IN.cram"; printf 'crai' > "$WORK/IN.cram.crai"
PATH="$FAKE:$PATH" nextflow run "$ROOT/main.nf" -profile standard \
  --sample CRAMIN --reference "$WORK/ref.fa" --bam "$WORK/IN.cram" \
  --vep_cache "$WORK/vep" --clinvar_vcf "$WORK/clinvar.vcf.gz" \
  --exons_bed "$ROOT/tests/data/exons.test.bed" --threads 2 --outdir "$WORK/out_cram" --mosaic
python3 -c "import json,sys; r=json.load(open(sys.argv[1])); assert r['input']['mode']=='bam' and r['input']['files']==['IN.cram'], r['input']; assert r['genes']['BRCA1']['smallVariants'] and r['mosaicCalling']['enabled']; print('stub pipeline (CRAM input) OK')" "$WORK/out_cram/brca12_results.json"

# Spaces in paths: the pipeline folder, the input folder and the output folder each contain a space
# (for example a folder called "brca12-wgs Start 2" in Downloads).
SP="$WORK/with space"; mkdir -p "$SP/pipeline copy" "$SP/in put"
cp -R "$ROOT/main.nf" "$ROOT/nextflow.config" "$ROOT/bin" "$ROOT/assets" "$ROOT/envs" "$SP/pipeline copy/"
cp "$WORK"/ref.* "$WORK"/R1.fq.gz "$WORK"/R2.fq.gz "$WORK"/clinvar.vcf.gz "$WORK"/clinvar.vcf.gz.tbi "$SP/in put/"
PATH="$FAKE:$PATH" nextflow run "$SP/pipeline copy/main.nf" -profile standard \
  --sample SPACE --reference "$SP/in put/ref.fa" --reads_r1 "$SP/in put/R1.fq.gz" --reads_r2 "$SP/in put/R2.fq.gz" \
  --vep_cache "$WORK/vep" --clinvar_vcf "$SP/in put/clinvar.vcf.gz" \
  --exons_bed "$ROOT/tests/data/exons.test.bed" --threads 2 --outdir "$SP/out put"
for f in brca12_report.html brca12_report.md brca12_results.json brca12_commands.sh brca12_report.xlsx; do
  test -s "$SP/out put/$f" || { echo "missing $f (paths with spaces)"; exit 1; }
done
echo "stub pipeline (paths with spaces) OK"
