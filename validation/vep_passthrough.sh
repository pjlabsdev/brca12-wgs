#!/usr/bin/env bash
# Stand-in for Ensembl VEP during the GIAB benchmark when no VEP cache is available.
# It copies the variants unchanged and only declares an empty CSQ field, so the
# benchmark measures calling, not annotation. Never use it for real analyses.
set -euo pipefail
if [[ "${1:-}" == "--help" ]]; then echo "ensembl-vep : passthrough-benchmark-only"; exit 0; fi
in=""; out=""
while [[ $# -gt 0 ]]; do
  case "$1" in -i) in="$2"; shift 2;; -o) out="$2"; shift 2;; *) shift;; esac
done
fields="Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|VARIANT_CLASS|CANONICAL|MANE_SELECT|PICK|MAX_AF|gnomADe_AF|gnomADg_AF|ClinVar|ClinVar_CLNSIG|ClinVar_CLNREVSTAT|ClinVar_CLNDN|ClinVar_CLNSIGCONF"
hdr="$(mktemp)"
echo "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: $fields\">" > "$hdr"
bcftools annotate --header-lines "$hdr" "$in" -Oz -o "$out"
rm -f "$hdr"
