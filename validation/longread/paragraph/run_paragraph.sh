#!/bin/bash
# Paragraph 2.3 (bioconda) on the frozen samples, same allele sequences as the realignment catalogue.
# Usage: run_paragraph.sh MANIFEST REFERENCE OUTDIR
# MANIFEST: id, path (BAM of the three windows), depth (30), read length (150), tab separated with header.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
multigrmpy.py -i "$here/variants.vcf" -m "$1" -r "$2" -o "$3" -t 4 -M 400
