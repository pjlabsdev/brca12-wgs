#!/usr/bin/env bash
# Assembly-based truth: samples of the Human Pangenome Reference Consortium (HPRC) year 1.
#
# Truth:   the HPRC v1.1 minigraph-cactus pangenome, GRCh38-based VCF (decomposed with vcfwave),
#          built from haplotype-resolved long-read assemblies. Independent of short-read callers.
# Reads:   the same individuals in the 1000 Genomes NYGC 30x data (hprc_samples_1000g.tsv),
#          prepared as in pathogenic/run_sv_carriers.sh with complete read pairs.
# Stand-ins: VEP (passthrough) and CNVkit (neutral segments).
#
# Usage: validation/hprc/run_hprc.sh WORKDIR     (WORKDIR as prepared by validation/run_giab.sh)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
HERE="$ROOT/validation/hprc"
WORK="${1:?usage: run_hprc.sh WORKDIR}"
VCF=https://human-pangenomics.s3.amazonaws.com/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz
mkdir -p "$WORK/hprc"
for r in chr13:32315077-32400268 chr17:43044295-43170327; do
  out="$WORK/hprc/hprc.${r%%:*}.vcf.gz"
  [[ -s $out.tbi ]] || { bcftools view "$VCF" "$r" -Oz -o "$out" && tabix -f "$out"; }
done
MATES=1 OUTDIR=hprc_runs TABLE_ONLY=1 TABLE="$HERE/hprc_samples_1000g.tsv" \
  bash "$ROOT/validation/pathogenic/run_sv_carriers.sh" "$WORK"
python3 "$HERE/hprc_eval.py" --samples "$HERE/hprc_samples_1000g.tsv" --truth-dir "$WORK/hprc" \
  --workdir "$WORK/hprc_runs" --sdf "$WORK/data/ref.sdf" --out "$WORK/hprc_runs/hprc_summary"
