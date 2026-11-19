#!/bin/bash
# brca12-wgs Start: double-click (macOS: .command opens in Terminal) or run `bash brca12-wgs-start.command`.
# Asks for a BAM or CRAM (drag the file into the window, Terminal inserts its path, or type the path),
# shows input and output folder, starts the Nextflow pipeline after Enter and opens the report.
# Reference, VEP cache, ClinVar VCF and Nextflow options are asked once and kept in
# ~/.config/brca12-wgs/start.conf (plain shell variables, editable).
# BRCA12_DRYRUN=1 stops after the overview without starting the pipeline.
set -u

HERE="$(cd "$(dirname "$0")" && pwd)"
CONF="$HOME/.config/brca12-wgs/start.conf"
mkdir -p "$(dirname "$CONF")"
REPO="" REFERENCE="" VEP_CACHE="" CLINVAR_VCF="" PROFILE="conda" NEXTFLOW="" EXTRA_ARGS=""
[ -f "$CONF" ] && . "$CONF"

if [ "${LANG:-}" != "${LANG#de}" ] || [ "$(defaults read -g AppleLocale 2>/dev/null | cut -c1-2)" = "de" ]; then DE=1; else DE=0; fi
say() { if [ "$DE" = 1 ]; then printf '%s\n' "$1"; else printf '%s\n' "$2"; fi; }
ask() { if [ "$DE" = 1 ]; then printf '\n%s\n> ' "$1"; else printf '\n%s\n> ' "$2"; fi; }
finish() { printf '\n'; say "Enter zum Schließen." "Press Enter to close."; read -r _; exit "${1:-0}"; }

# a path as dragged into Terminal: quotes or backslash escapes, trailing space
clean() {
  local p="$1"
  p="${p#"${p%%[![:space:]]*}"}"; p="${p%"${p##*[![:space:]]}"}"
  case "$p" in \'*\') p="${p#\'}"; p="${p%\'}" ;; \"*\") p="${p#\"}"; p="${p%\"}" ;; *) p="$(printf '%s' "$p" | sed 's/\\\(.\)/\1/g')" ;; esac
  case "$p" in "~/"*) p="$HOME/${p#\~/}" ;; esac
  printf '%s' "$p"
}

# askpath DE EN KIND [DEFAULT]: KIND f = file, d = directory; Enter keeps DEFAULT
askpath() {
  local p
  while true; do
    if [ -n "${4:-}" ]; then ask "$1 [Enter: $4]" "$2 [Enter: $4]"; else ask "$1" "$2"; fi
    IFS= read -r p || finish 1
    p="$(clean "$p")"
    [ -z "$p" ] && p="${4:-}"
    if [ "$3" = d ] && [ -d "$p" ]; then REPLY_PATH="$p"; return; fi
    if [ "$3" = f ] && [ -f "$p" ]; then REPLY_PATH="$p"; return; fi
    say "Nicht gefunden: $p" "Not found: $p"
  done
}

save() {
  {
    printf '# brca12-wgs Start settings (written by brca12-wgs-start.command)\n'
    for v in REPO REFERENCE VEP_CACHE CLINVAR_VCF PROFILE NEXTFLOW EXTRA_ARGS; do printf '%s=%q\n' "$v" "${!v}"; done
  } > "$CONF"
}

clear 2>/dev/null
say "brca12-wgs · BRCA1/BRCA2 aus Genom-, Exom- oder Paneldaten" "brca12-wgs · BRCA1/BRCA2 from genome, exome or panel data"
say "Nur für Forschung (RUO). Ergebnisse brauchen eine klinisch validierte Bestätigung." "Research use only (RUO). Results need clinically validated confirmation."

# pipeline folder (main.nf): next to this file, one level up, or as remembered
if [ ! -f "$REPO/main.nf" ]; then
  REPO=""
  for c in "$HERE/.." "$HERE"; do [ -f "$c/main.nf" ] && REPO="$(cd "$c" && pwd)" && break; done
fi
if [ -z "$REPO" ]; then askpath "Ordner von brca12-wgs (mit main.nf) hineinziehen:" "Drop the brca12-wgs folder (with main.nf):" d; REPO="$REPLY_PATH"; fi
[ -f "$REPO/main.nf" ] || { say "In $REPO liegt keine main.nf." "No main.nf in $REPO."; finish 1; }
[ -z "$NEXTFLOW" ] && NEXTFLOW="$(command -v nextflow 2>/dev/null)"
if [ -z "$NEXTFLOW" ] || [ ! -x "$NEXTFLOW" ]; then
  askpath "Nextflow wurde nicht gefunden. Programm nextflow hineinziehen:" "Nextflow was not found. Drop the nextflow program:" f; NEXTFLOW="$REPLY_PATH"
fi

# input
askpath "BAM oder CRAM hier hineinziehen oder Pfad eingeben, dann Enter:" "Drop a BAM or CRAM here or enter its path, then press Enter:" f
IN="$REPLY_PATH"
lower="$(printf '%s' "$IN" | tr '[:upper:]' '[:lower:]')"
case "$lower" in
  *.bam)
    if [ ! -f "$IN.bai" ] && [ ! -f "${IN%.*}.bai" ] && [ ! -f "$IN.csi" ]; then say "Zur BAM fehlt der Index (.bai daneben)." "The BAM index (.bai) is missing."; finish 1; fi
    askpath "Referenz (FASTA mit .fai), mit der die BAM ausgerichtet wurde:" "Reference (FASTA with .fai) the BAM was aligned to:" f "$REFERENCE" ;;
  *.cram)
    if [ ! -f "$IN.crai" ] && [ ! -f "${IN%.*}.crai" ]; then say "Zur CRAM fehlt der Index (.crai daneben)." "The CRAM index (.crai) is missing."; finish 1; fi
    askpath "Referenz der CRAM (FASTA mit .fai) hineinziehen:" "Drop the reference of the CRAM (FASTA with .fai):" f "$REFERENCE" ;;
  *) say "Unbekanntes Format. Erlaubt sind BAM oder CRAM." "Unknown format. Allowed are BAM or CRAM."; finish 1 ;;
esac
REFERENCE="$REPLY_PATH"
[ -f "$REFERENCE.fai" ] || { say "Zur Referenz fehlt der Index (.fai daneben)." "The reference index (.fai) is missing."; finish 1; }
[ -d "$VEP_CACHE" ] || { askpath "VEP-Cache-Ordner (einmalig) hineinziehen:" "Drop the VEP cache folder (once):" d; VEP_CACHE="$REPLY_PATH"; }
[ -f "$CLINVAR_VCF" ] || { askpath "ClinVar-VCF (clinvar.vcf.gz, einmalig) hineinziehen:" "Drop the ClinVar VCF (clinvar.vcf.gz, once):" f; CLINVAR_VCF="$REPLY_PATH"; }
save

# sample name and a new output folder
S="$(basename "$IN")"
for x in .bam .cram .BAM .CRAM; do S="${S%"$x"}"; done
[ -z "$S" ] && S=sample
OUT="$HOME/Downloads/brca12-$S"; n=2
while [ -e "$OUT" ]; do OUT="$HOME/Downloads/brca12-$S-$n"; n=$((n + 1)); done
T="$( (sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4) | head -1)"
if [ "$DE" = 1 ]; then RL=de; else RL=en; fi

say "" ""
say "Eingabe:    $IN" "Input:      $IN"
say "Referenz:   $REFERENCE" "Reference:  $REFERENCE"
say "Ausgabe:    $OUT" "Output:     $OUT"
say "Threads:    $T · Profil $PROFILE · Bericht $RL" "Threads:    $T · profile $PROFILE · report $RL"
say "Einstellungen: $CONF" "Settings:   $CONF"
if [ "${BRCA12_DRYRUN:-0}" = 1 ]; then say "" ""; say "Trockenlauf: hier würde der Lauf starten." "Dry run: the run would start here."; exit 0; fi
ask "Start mit Enter (oder Strg-C zum Abbrechen). Ein Genom dauert je nach Rechner mehrere Stunden." "Press Enter to start (or Ctrl-C to cancel). A genome takes several hours depending on the machine."; read -r _

mkdir -p "$OUT"
cd "$OUT" || finish 1
start=$(date +%s)
# shellcheck disable=SC2086  # EXTRA_ARGS holds further Nextflow options, split on purpose
"$NEXTFLOW" run "$REPO/main.nf" -profile "$PROFILE" -work-dir "$OUT/work" \
  --sample "$S" --bam "$IN" --reference "$REFERENCE" --vep_cache "$VEP_CACHE" --clinvar_vcf "$CLINVAR_VCF" \
  --report_lang "$RL" --threads "$T" --outdir "$OUT/results" $EXTRA_ARGS 2>&1 | tee "$OUT/brca12-start.log"
rc=${PIPESTATUS[0]}
mins=$(( ($(date +%s) - start) / 60 ))

printf '\n'
if [ "$rc" != 0 ]; then say "Der Lauf ist mit Fehler $rc beendet worden. Protokoll: $OUT/brca12-start.log" "The run stopped with error $rc. Log: $OUT/brca12-start.log"; finish "$rc"; fi
say "Fertig nach $mins min. Ergebnisse in $OUT/results:" "Done after $mins min. Results in $OUT/results:"
find "$OUT/results" -maxdepth 2 -name 'brca12_report*' -print 2>/dev/null
REPORT="$(find "$OUT/results" -maxdepth 2 -name 'brca12_report.html' -print 2>/dev/null | head -1)"
if command -v open >/dev/null 2>&1; then open "$OUT/results"; [ -n "$REPORT" ] && open "$REPORT"; fi
finish 0
