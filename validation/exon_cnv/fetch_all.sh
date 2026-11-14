#!/bin/bash
# usage: fetch_all.sh PY REF BED OUTDIR LIST.tsv
PY=$1; REF=$2; BED=$3; OUT=$4; LIST=$5
cut -f1-3 "$LIST" | xargs -P 6 -L 1 bash -c 'for i in 1 2 3; do timeout 900 "$0" -I fetch.py "$1" "$2" "$3" "$4" "$5" "$6" 2>/dev/null && break; done' "$PY" "$REF" "$BED" "$OUT"
ls "$OUT"/*.json | wc -l
