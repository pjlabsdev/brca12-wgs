#!/bin/bash
# brca12-wgs Start: double-click (macOS: .command opens in Terminal) or run `bash brca12-wgs-start.command`.
# Asks for a BAM or CRAM (drag the file into the window, Terminal inserts its path, or type the path),
# shows input and output folder, starts the Nextflow pipeline after Enter with an animated progress
# display (one line per finished step, read from the Nextflow trace) and opens the report.
# Reference, VEP cache, ClinVar VCF and Nextflow options are asked once and kept in
# ~/.config/brca12-wgs/start.conf (plain shell variables, editable); an optional env.sh next to it
# can set the environment (for example JAVA_HOME and PATH for a local tool installation).
# BRCA12_DRYRUN=1 stops after the overview without starting the pipeline.
set -u
case "${LC_ALL:-${LC_CTYPE:-${LANG:-}}}" in *UTF-8*) ;; *) export LC_CTYPE=en_US.UTF-8 ;; esac

HERE="$(cd "$(dirname "$0")" && pwd)"
CONF="$HOME/.config/brca12-wgs/start.conf"
mkdir -p "$(dirname "$CONF")"
REPO="" REFERENCE="" VEP_CACHE="" CLINVAR_VCF="" PROFILE="conda" NEXTFLOW="" EXTRA_ARGS=""
[ -f "$CONF" ] && . "$CONF"
# optional environment for Nextflow and the tools (JAVA_HOME, PATH), kept apart from start.conf
[ -f "$(dirname "$CONF")/env.sh" ] && . "$(dirname "$CONF")/env.sh"

if [ "${LANG:-}" != "${LANG#de}" ] || [ "$(defaults read -g AppleLocale 2>/dev/null | cut -c1-2)" = "de" ]; then DE=1; else DE=0; fi
tr_() { if [ "$DE" = 1 ]; then printf '%s' "$1"; else printf '%s' "$2"; fi; }

# ── display (PJ Labs: gold accent, no bold) ──────────────────────────────────
TTY=0; [ -t 1 ] && TTY=1
if [ "$TTY" = 1 ] && [ -z "${NO_COLOR:-}" ]; then
  RST=$'\e[0m'; ACC=$'\e[38;5;178m'; OK=$'\e[38;5;71m'; ERR=$'\e[38;5;167m'; MUT=$'\e[38;5;244m'; FNT=$'\e[38;5;238m'
else
  RST=""; ACC=""; OK=""; ERR=""; MUT=""; FNT=""
fi
W=64
FRAMES=(⠋ ⠙ ⠹ ⠸ ⠼ ⠴ ⠦ ⠧ ⠇ ⠏)
PID=""
rep() { local o="" i=0; while [ "$i" -lt "$2" ]; do o="$o$1"; i=$((i + 1)); done; printf '%s' "$o"; }
cursor() { [ "$TTY" = 1 ] && printf '\e[?25%s' "$1"; }
title() { [ "$TTY" = 1 ] && printf '\e]0;%s\a' "$1"; }
now() { perl -MTime::HiRes=time -e 'printf "%.3f", time' 2>/dev/null || date +%s; }
dur() { awk -v s="$1" 'BEGIN { s = int(s); if (s < 60) printf "%ds", s; else if (s < 3600) printf "%dm %02ds", s / 60, s % 60; else printf "%dh %02dm", s / 3600, (s % 3600) / 60 }'; }
since() { dur "$(awk -v a="$1" -v b="$(now)" 'BEGIN { print b - a }')"; }
row() {  # row INDENT LEFT RIGHT [COLOUR_LEFT] [COLOUR_RIGHT]
  local gap=$((W - $1 - ${#2} - ${#3})); [ "$gap" -lt 1 ] && gap=1
  printf '%s%s%s%s%s%s%s%s\n' "$(rep ' ' "$1")" "${4:-}" "$2" "$RST" "$(rep ' ' "$gap")" "${5:-$MUT}" "$3" "$RST"
}
divider() { printf '  %s%s%s\n' "$FNT" "$(rep ─ $((W - 2)))" "$RST"; }
result() {  # result ok|fail LABEL DETAIL TIME
  local sym="✓" col="$OK" l="$2"; [ "$1" = fail ] && sym="✗" && col="$ERR"
  [ -n "$3" ] && l="$l  $3"
  local gap=$((W - 6 - ${#l} - ${#4})); [ "$gap" -lt 1 ] && gap=1
  local text="$2"; [ -n "$3" ] && text="$text  $MUT$3$RST"; [ "$1" = fail ] && text="$ERR$2${3:+  $3}$RST"
  printf '    %s%s%s %s%s%s%s%s\n' "$col" "$sym" "$RST" "$text" "$(rep ' ' "$gap")" "$MUT" "$4" "$RST"
}
spin() { [ "$TTY" = 1 ] && printf '\r\e[2K    %s%s%s %s' "$ACC" "${FRAMES[$(($1 % 10))]}" "$RST" "$2"; }
clearline() { [ "$TTY" = 1 ] && printf '\r\e[2K'; }
bar() {  # bar DONE TOTAL WIDTH
  local n=0; [ "$2" -gt 0 ] && n=$(($1 * $3 / $2)); [ "$n" -gt "$3" ] && n=$3
  printf '%s%s%s%s%s' "$ACC" "$(rep ━ "$n")" "$FNT" "$(rep ━ $(($3 - n)))" "$RST"
}
ask() { printf '\n  %s\n  %s›%s ' "$(tr_ "$1" "$2")" "$ACC" "$RST"; }
finish() { cursor h; printf '\n'; printf '  %s%s%s' "$MUT" "$(tr_ "Enter zum Schließen" "Press Enter to close")" "$RST"; read -r _; exit "${1:-0}"; }
stop() { printf '  %s%s%s\n' "$ERR" "$1" "$RST"; finish 1; }
trap 'cursor h' EXIT
trap 'cursor h; [ -n "$PID" ] && kill "$PID" 2>/dev/null; exit 130' INT TERM HUP

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
    printf '  %s%s %s%s\n' "$ERR" "$(tr_ "Nicht gefunden:" "Not found:")" "$p" "$RST"
  done
}

save() {
  {
    printf '# brca12-wgs Start settings (written by brca12-wgs-start.command)\n'
    for v in REPO REFERENCE VEP_CACHE CLINVAR_VCF PROFILE NEXTFLOW EXTRA_ARGS; do printf '%s=%q\n' "$v" "${!v}"; done
  } > "$CONF"
}
short() { case "$1" in "$HOME"/*) printf '~/%s' "${1#"$HOME"/}" ;; *) printf '%s' "$1" ;; esac; }

# ── header ───────────────────────────────────────────────────────────────────
clear 2>/dev/null
title "brca12-wgs"
printf '\n'
row 2 "brca12-wgs" "$(date '+%d.%m.%Y  %H:%M')" "$ACC" "$MUT"
printf '  %s%s%s\n' "$MUT" "$(tr_ "BRCA1/BRCA2 aus Genom-, Exom- oder Paneldaten · nur für Forschung (RUO)" "BRCA1/BRCA2 from genome, exome or panel data · research use only (RUO)")" "$RST"
divider

# pipeline folder (main.nf): next to this file, one level up, or as remembered
if [ ! -f "$REPO/main.nf" ]; then
  REPO=""
  for c in "$HERE/.." "$HERE"; do [ -f "$c/main.nf" ] && REPO="$(cd "$c" && pwd)" && break; done
fi
if [ -z "$REPO" ]; then askpath "Ordner von brca12-wgs (mit main.nf) hineinziehen:" "Drop the brca12-wgs folder (with main.nf):" d; REPO="$REPLY_PATH"; fi
[ -f "$REPO/main.nf" ] || stop "$(tr_ "In $REPO liegt keine main.nf." "No main.nf in $REPO.")"
[ -z "$NEXTFLOW" ] && NEXTFLOW="$(command -v nextflow 2>/dev/null)"
if [ -z "$NEXTFLOW" ] || [ ! -x "$NEXTFLOW" ]; then
  askpath "Nextflow wurde nicht gefunden. Programm nextflow hineinziehen:" "Nextflow was not found. Drop the nextflow program:" f; NEXTFLOW="$REPLY_PATH"
fi

# ── input ────────────────────────────────────────────────────────────────────
askpath "BAM oder CRAM hier hineinziehen oder Pfad eingeben, dann Enter:" "Drop a BAM or CRAM here or enter its path, then press Enter:" f
IN="$REPLY_PATH"
lower="$(printf '%s' "$IN" | tr '[:upper:]' '[:lower:]')"
case "$lower" in
  *.bam)
    if [ ! -f "$IN.bai" ] && [ ! -f "${IN%.*}.bai" ] && [ ! -f "$IN.csi" ]; then stop "$(tr_ "Zur BAM fehlt der Index (.bai daneben)." "The BAM index (.bai) is missing.")"; fi
    askpath "Referenz (FASTA mit .fai), mit der die BAM ausgerichtet wurde:" "Reference (FASTA with .fai) the BAM was aligned to:" f "$REFERENCE" ;;
  *.cram)
    if [ ! -f "$IN.crai" ] && [ ! -f "${IN%.*}.crai" ]; then stop "$(tr_ "Zur CRAM fehlt der Index (.crai daneben)." "The CRAM index (.crai) is missing.")"; fi
    askpath "Referenz der CRAM (FASTA mit .fai) hineinziehen:" "Drop the reference of the CRAM (FASTA with .fai):" f "$REFERENCE" ;;
  *) stop "$(tr_ "Unbekanntes Format. Erlaubt sind BAM oder CRAM." "Unknown format. Allowed are BAM or CRAM.")" ;;
esac
REFERENCE="$REPLY_PATH"
[ -f "$REFERENCE.fai" ] || stop "$(tr_ "Zur Referenz fehlt der Index (.fai daneben)." "The reference index (.fai) is missing.")"
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

printf '\n'
divider
row 2 "$(tr_ "Probe" "Sample")" "$S" "$MUT" ""
row 2 "$(tr_ "Eingabe" "Input")" "$(basename "$IN")" "$MUT" ""
row 2 "$(tr_ "Referenz" "Reference")" "$(basename "$REFERENCE")" "$MUT" ""
row 2 "$(tr_ "Ausgabe" "Output")" "$(short "$OUT")" "$MUT" ""
row 2 "Threads" "$T  ·  $PROFILE  ·  $(tr_ "Bericht" "report") $RL" "$MUT" ""
divider
if [ "${BRCA12_DRYRUN:-0}" = 1 ]; then printf '\n  %s%s%s\n' "$MUT" "$(tr_ "Trockenlauf: hier würde der Lauf starten." "Dry run: the run would start here.")" "$RST"; exit 0; fi
ask "Start mit Enter (oder Strg-C zum Abbrechen). Ein Genom dauert je nach Rechner mehrere Stunden." "Press Enter to start (or Ctrl-C to cancel). A genome takes several hours depending on the machine."; read -r _ || exit 1

mkdir -p "$OUT" || finish 1
cd "$OUT" || finish 1
LOG="$OUT/brca12-start.log"
TRACE="$OUT/results/pipeline_info/trace.tsv"

# ── run with progress ────────────────────────────────────────────────────────
printf '\n'
title "brca12-wgs · $S"
START=$(now)
# shellcheck disable=SC2086  # EXTRA_ARGS holds further Nextflow options, split on purpose
"$NEXTFLOW" -log "$OUT/.nextflow.log" run "$REPO/main.nf" -ansi-log false -profile "$PROFILE" -work-dir "$OUT/work" \
  --sample "$S" --bam "$IN" --reference "$REFERENCE" --vep_cache "$VEP_CACHE" --clinvar_vcf "$CLINVAR_VCF" \
  --report_lang "$RL" --threads "$T" --outdir "$OUT/results" $EXTRA_ARGS >"$LOG" 2>&1 & PID=$!

# progress: submitted steps from the Nextflow log, finished steps from the trace file (one tick line
# each with its run time), running steps next to the spinner
shown=0; i=0; cursor l
show_done() {  # tick lines for trace rows not shown yet
  [ -f "$TRACE" ] || return 0
  local rows line name status rt
  rows=$(($(wc -l < "$TRACE") - 1))
  while [ "$shown" -lt "$rows" ]; do
    shown=$((shown + 1))
    line=$(awk -F'\t' -v k=$((shown + 1)) 'NR == k { print $4 "\t" $5 "\t" $9 }' "$TRACE")
    name=$(printf '%s' "$line" | cut -f1 | sed 's/ (.*//'); status=$(printf '%s' "$line" | cut -f2); rt=$(printf '%s' "$line" | cut -f3)
    clearline
    case "$status" in
      COMPLETED) result ok "$name" "" "$rt" ;;
      CACHED) result ok "$name" "$(tr_ "aus Zwischenspeicher" "cached")" "$rt" ;;
      *) result fail "$name" "$status" "$rt" ;;
    esac
  done
}
running() {  # names of submitted steps that are not in the trace yet
  { [ -f "$TRACE" ] && sed 1d "$TRACE" | cut -f4 | sed 's/^/D	/'
    grep 'Submitted process >' "$LOG" 2>/dev/null | sed 's/.*Submitted process > /S	/'; } |
    awk -F'\t' '$1 == "D" { d[$2] = 1; next } $1 == "S" && !($2 in d) { n = $2; sub(/ \(.*/, "", n); o = o (c++ ? ", " : "") n } END { print o }'
}
while kill -0 "$PID" 2>/dev/null; do
  show_done
  sub=$(grep -c 'Submitted process >' "$LOG" 2>/dev/null)
  run_now="$(running)"
  spin "$i" "$(bar "$shown" "${sub:-0}" 20)  $MUT$shown / ${sub:-0}  ${run_now:-$(tr_ "startet" "starting")}  $(since "$START")$RST"
  i=$((i + 1)); sleep 0.25
done
wait "$PID"; rc=$?; PID=""
show_done; clearline; cursor h

printf '\n'
divider
if [ "$rc" != 0 ]; then
  row 2 "$(tr_ "Fehler $rc" "Error $rc")" "$(since "$START")" "$ERR"
  grep -i -e 'error' -e 'caused by' "$LOG" | tail -4 | sed "s/^/    $MUT/; s/\$/$RST/"
  printf '    %s%s %s%s\n' "$MUT" "$(tr_ "Protokoll:" "Log:")" "$LOG" "$RST"
  finish "$rc"
fi
row 2 "$(tr_ "Fertig" "Done")" "$(since "$START")" "$OK"
printf '\n'
find "$OUT/results" -maxdepth 2 -name 'brca12_report*' -print 2>/dev/null | while IFS= read -r f; do printf '  %s%s%s\n' "$MUT" "$(short "$f")" "$RST"; done
REPORT="$(find "$OUT/results" -maxdepth 2 -name 'brca12_report.html' -print 2>/dev/null | head -1)"
title "brca12-wgs · $S · $(tr_ "fertig" "done")"
if command -v open >/dev/null 2>&1; then open "$OUT/results"; [ -n "$REPORT" ] && open "$REPORT"; fi
finish 0
