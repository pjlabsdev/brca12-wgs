#!/bin/bash
# brca12-wgs Start: double-click (macOS: .command opens in Terminal), on Linux run
# `bash brca12-wgs-start.sh` (the same file), on Windows double-click brca12-wgs-start.bat (runs it in WSL 2).
# Language: English; German when the system language is German. BRCA12_LANG=en or =de sets it.
# Asks for a BAM or CRAM (drag the file into the window, Terminal inserts its path, or type the path),
# shows input and output folder, starts the Nextflow pipeline after Enter with an animated progress
# display (one line per finished step, read from the Nextflow trace) and opens the report.
# Reference, VEP cache, ClinVar VCF and Nextflow options are asked once and kept in
# ~/.config/brca12-wgs/start.conf (plain shell variables, editable); an optional env.sh next to it
# can set the environment (for example JAVA_HOME and PATH for a local tool installation).
# BRCA12_DRYRUN=1 stops after the overview without starting the pipeline or saving settings.
# At every start it asks GitHub for the latest release (one request, 5 s limit, silent on errors; nothing
# about you or your data is sent) and offers the update (default no). The new package replaces this
# folder's program files in place and is started; results and settings stay. BRCA12_NO_UPDATE=1 switches
# the check off.
# If the input or reference disk drops out (external USB disk), the start file waits until both are
# readable again and continues with -resume (finished steps are kept), at most BRCA12_RESUME times (3).
# If a task hangs on the lost disk, Nextflow is stopped after BRCA12_STALL seconds (60) without the
# input or reference, so that the resume can start.
set -u
case "${LC_ALL:-${LC_CTYPE:-${LANG:-}}}" in *UTF-8*) ;; *) export LC_CTYPE=en_US.UTF-8 ;; esac

HERE="$(cd "$(dirname "$0")" && pwd)"
# macOS: a downloaded start folder carries the quarantine flag and its programs would be blocked on
# other Macs; drop it for this folder only (start file and the pipeline shipped next to it).
if [ "$(uname)" = Darwin ]; then xattr -dr com.apple.quarantine "$HERE" 2>/dev/null || true; fi
CONF="$HOME/.config/brca12-wgs/start.conf"
mkdir -p "$(dirname "$CONF")"
REPO="" REFERENCE="" VEP_CACHE="" CLINVAR_VCF="" PROFILE="conda" NEXTFLOW="" EXTRA_ARGS=""
[ -f "$CONF" ] && . "$CONF"
# optional environment for Nextflow and the tools (JAVA_HOME, PATH), kept apart from start.conf
[ -f "$(dirname "$CONF")/env.sh" ] && . "$(dirname "$CONF")/env.sh"

# English is the main language; German is used when the system language is German or with BRCA12_LANG=de
WSL=0; { grep -qi microsoft /proc/version 2>/dev/null || [ "${BRCA12_WSL:-0}" = 1 ]; } && command -v wslpath >/dev/null 2>&1 && WSL=1
case "${BRCA12_LANG:-}" in
  de*|DE*) DE=1 ;;
  ?*) DE=0 ;;
  *) syslang="${LC_ALL:-${LC_MESSAGES:-${LANG:-}}}"
     [ "$(uname)" = Darwin ] && syslang="$(defaults read -g AppleLocale 2>/dev/null || printf '%s' "$syslang")"
     [ "$WSL" = 1 ] && syslang="$(powershell.exe -NoProfile -Command '(Get-UICulture).Name' 2>/dev/null | tr -d '\r' || printf '%s' "$syslang")"
     case "$syslang" in de*) DE=1 ;; *) DE=0 ;; esac ;;
esac
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
  # Windows path (C:\Users\..., dragged into the window with or without quotes): translate for WSL
  if [ "$WSL" = 1 ]; then
    local w="$p"; w="${w#\"}"; w="${w%\"}"
    case "$w" in [A-Za-z]:[\\/]*) wslpath -u "$w"; return ;; esac
  fi
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
row 2 "brca12-wgs" "$(if [ "$DE" = 1 ]; then date '+%d.%m.%Y  %H:%M'; else date '+%Y-%m-%d  %H:%M'; fi)" "$ACC" "$MUT"
printf '  %s%s%s\n' "$MUT" "$(tr_ "BRCA1/BRCA2 aus Genom-, Exom- oder Paneldaten · nur für Forschung (RUO)" "BRCA1/BRCA2 from genome, exome or panel data · research use only (RUO)")" "$RST"
divider

# pipeline folder (main.nf): the one shipped with this file (pipeline/ next to it, as in the start
# ZIP, or the repository one level up) comes first, then the remembered one
for c in "$HERE/pipeline" "$HERE/.." "$HERE"; do
  if [ -f "$c/main.nf" ]; then REPO="$(cd "$c" && pwd)"; break; fi
done
[ -f "$REPO/main.nf" ] || REPO=""
if [ -z "$REPO" ]; then askpath "Ordner von brca12-wgs (mit main.nf) hineinziehen:" "Drop the brca12-wgs folder (with main.nf):" d; REPO="$REPLY_PATH"; fi
[ -f "$REPO/main.nf" ] || stop "$(tr_ "In $REPO liegt keine main.nf." "No main.nf in $REPO.")"

# ── update check ─────────────────────────────────────────────────────────────
# newer A B: true when version A (1.2.3 or v1.2.3) is newer than B
newer() { awk -v a="${1#v}" -v b="${2#v}" 'BEGIN { n = split(a, x, "."); m = split(b, y, "."); for (i = 1; i <= (n > m ? n : m); i++) { if (x[i] + 0 > y[i] + 0) exit 0; if (x[i] + 0 < y[i] + 0) exit 1 } exit 1 }'; }
VERSION="$(sed -n "s/^ *version *= *'\([^']*\)'.*/\1/p" "$REPO/nextflow.config" 2>/dev/null | head -1)"
case "$(uname)" in Darwin) PKG_OS=macos; PKG_START=brca12-wgs-start.command ;; *) PKG_OS=linux; PKG_START=brca12-wgs-start.sh ;; esac
[ "$WSL" = 1 ] && PKG_OS=windows
update_check() {
  [ "${BRCA12_NO_UPDATE:-0}" = 1 ] && return 0
  [ -n "$VERSION" ] && command -v curl >/dev/null 2>&1 || return 0
  local api="${BRCA12_UPDATE_URL:-https://api.github.com/repos/pjlabsdev/brca12-wgs/releases/latest}" json latest url ans tmp
  json="$(curl -fsSL --max-time 5 -H 'Accept: application/vnd.github+json' "$api" 2>/dev/null)" || return 0
  latest="$(printf '%s' "$json" | tr ',' '\n' | sed -n 's/.*"tag_name" *: *"v\{0,1\}\([0-9][0-9.]*\)".*/\1/p' | head -1)"
  [ -n "$latest" ] && newer "$latest" "$VERSION" || return 0
  row 2 "$(tr_ "Neue Version $latest verfügbar" "New release $latest available")" "$(tr_ "installiert" "installed") $VERSION" "$ACC" "$MUT"
  url="$(printf '%s' "$json" | tr ',' '\n' | sed -n 's/.*"browser_download_url" *: *"\([^"]*\/brca12-wgs-'"$PKG_OS"'\.zip\)".*/\1/p' | head -1)"
  # only a start package (pipeline/ next to this file) can update itself; a git checkout is updated with git
  if [ -z "$url" ] || [ ! -f "$HERE/pipeline/main.nf" ] || ! command -v unzip >/dev/null 2>&1; then
    printf '  %sDownload: https://github.com/pjlabsdev/brca12-wgs/releases/latest%s\n' "$MUT" "$RST"
    return 0
  fi
  ask "Jetzt aktualisieren? [j/N]" "Update now? [y/N]"
  IFS= read -r ans || return 0
  case "$ans" in [yYjJ]*) ;; *) return 0 ;; esac
  tmp="$(mktemp -d)" || return 0
  if ! curl -fsSL --max-time 300 -o "$tmp/update.zip" "$url" 2>/dev/null || ! unzip -q "$tmp/update.zip" -d "$tmp/x" 2>/dev/null \
     || [ ! -f "$tmp/x/brca12-wgs/pipeline/main.nf" ] || [ ! -f "$tmp/x/brca12-wgs/$PKG_START" ]; then
    rm -rf "$tmp"
    printf '  %s%s%s\n' "$ERR" "$(tr_ "Update fehlgeschlagen, weiter mit $VERSION." "Update failed, continuing with $VERSION.")" "$RST"
    return 0
  fi
  # in place: the pipeline is replaced as a whole, the files next to it one by one (each renamed into
  # place, so this running file is not overwritten under bash); results and settings stay
  if ! { rm -rf "$HERE/pipeline.update" && mv "$tmp/x/brca12-wgs/pipeline" "$HERE/pipeline.update" \
         && rm -rf "$HERE/pipeline.old" && mv "$HERE/pipeline" "$HERE/pipeline.old" && mv "$HERE/pipeline.update" "$HERE/pipeline"; }; then
    [ -d "$HERE/pipeline" ] || mv "$HERE/pipeline.old" "$HERE/pipeline" 2>/dev/null
    rm -rf "$tmp" "$HERE/pipeline.update"
    printf '  %s%s%s\n' "$ERR" "$(tr_ "Update fehlgeschlagen (kein Schreibrecht in $(short "$HERE")), weiter mit $VERSION." "Update failed (cannot write to $(short "$HERE")), continuing with $VERSION.")" "$RST"
    return 0
  fi
  rm -rf "$HERE/pipeline.old"
  local f
  for f in "$tmp/x/brca12-wgs"/*; do
    [ -f "$f" ] || continue
    cp -p "$f" "$HERE/.$(basename "$f").new" && mv -f "$HERE/.$(basename "$f").new" "$HERE/$(basename "$f")"
  done
  rm -rf "$tmp"
  printf '  %s%s%s\n' "$OK" "$(tr_ "Auf $latest aktualisiert." "Updated to $latest.")" "$RST"
  BRCA12_NO_UPDATE=1 exec bash "$HERE/$PKG_START"
}
update_check

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
# a dry run only shows what would happen and leaves the settings as they are
[ "${BRCA12_DRYRUN:-0}" = 1 ] || save

# sample name and a new output folder
S="$(basename "$IN")"
for x in .bam .cram .BAM .CRAM; do S="${S%"$x"}"; done
# the sample name becomes part of file names: only letters, digits, dot, underscore and hyphen
S="$(printf '%s' "$S" | LC_ALL=C sed 's/[^A-Za-z0-9._-]/_/g')"
[ -z "$S" ] && S=sample
# output next to the user's downloads; under WSL the Windows Downloads folder, so that Explorer finds it
BASE="$HOME/Downloads"
if [ "$WSL" = 1 ]; then
  winhome="$(cmd.exe /c 'echo %USERPROFILE%' 2>/dev/null | tr -d '\r')"
  [ -n "$winhome" ] && [ -d "$(wslpath -u "$winhome" 2>/dev/null)/Downloads" ] && BASE="$(wslpath -u "$winhome")/Downloads"
fi
[ -d "$BASE" ] || BASE="$HOME"
[ -n "${BRCA12_OUT:-}" ] && BASE="$BRCA12_OUT"
OUT="$BASE/brca12-$S"; n=2
while [ -e "$OUT" ]; do OUT="$BASE/brca12-$S-$n"; n=$((n + 1)); done
T="$( (sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4) | head -1)"
# results in English first (the HTML page carries German as well, one click away)
RL=en

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
launch() {
  "$NEXTFLOW" -log "$OUT/.nextflow.log" run "$REPO/main.nf" -ansi-log false -resume -profile "$PROFILE" -work-dir "$OUT/work" \
    --sample "$S" --bam "$IN" --reference "$REFERENCE" --vep_cache "$VEP_CACHE" --clinvar_vcf "$CLINVAR_VCF" \
    --report_lang "$RL" --threads "$T" --outdir "$OUT/results" $EXTRA_ARGS >>"$LOG" 2>&1 & PID=$!
}
# a run that ended because a disk dropped out: input or reference unreadable now, or an I/O error or
# a task killed by a signal (bus error, kill) in the log of this attempt, or the watchdog below stopped it
disk_dropout() {
  [ -r "$IN" ] && [ -r "$REFERENCE" ] || return 0
  tail -n +"$1" "$LOG" | grep -qiE 'Input/output error|Device not configured|Bus error|exit status \((13[4-9]|14[0-3])\)|unreadable for [0-9]+ s, Nextflow stopped'
}
# A disk that drops out while Nextflow runs can leave a task hanging on I/O, so Nextflow never ends and
# the resume below never starts. If the input or the reference stays unreadable for BRCA12_STALL
# seconds (60), Nextflow is stopped; the resume loop then waits for the disk and continues.
GONE=""
watch_disk() {
  if [ -r "$IN" ] && [ -r "$REFERENCE" ]; then GONE=""; return; fi
  [ -n "$GONE" ] || GONE=$(date +%s)
  [ "$GONE" = stopped ] && return
  if [ $(($(date +%s) - GONE)) -ge "${BRCA12_STALL:-60}" ]; then
    GONE=stopped
    kill "$PID" 2>/dev/null; sleep 5; kill -9 "$PID" 2>/dev/null
    echo "brca12-wgs start: input or reference unreadable for ${BRCA12_STALL:-60} s, Nextflow stopped" >> "$LOG"
  fi
}
launch

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
  watch_disk
  i=$((i + 1)); sleep 0.25
done
wait "$PID"; rc=$?; PID=""
show_done
attempt=1
while [ "$rc" != 0 ] && [ "$attempt" -le "${BRCA12_RESUME:-3}" ] && disk_dropout "${LOGLINE:-1}"; do
  clearline
  row 4 "$(tr_ "Datenträger nicht erreichbar, warte und setze fort ($attempt)" "Disk not reachable, waiting to resume ($attempt)")" "$(since "$START")" "$ERR"
  until [ -r "$IN" ] && [ -r "$REFERENCE" ] && [ -w "$OUT" ]; do sleep 15; done
  sleep 30
  attempt=$((attempt + 1)); LOGLINE=$(($(wc -l < "$LOG") + 1)); shown=0; GONE=""
  launch
  while kill -0 "$PID" 2>/dev/null; do
    show_done
    sub=$(grep -c 'Submitted process >' "$LOG" 2>/dev/null)
    spin "$i" "$(bar "$shown" "${sub:-0}" 20)  $MUT$(tr_ "fortgesetzt" "resumed")  $(since "$START")$RST"
    watch_disk
    i=$((i + 1)); sleep 0.25
  done
  wait "$PID"; rc=$?; PID=""
  show_done
done
clearline; cursor h

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
if [ "$WSL" = 1 ]; then
  [ -n "$REPORT" ] && explorer.exe "$(wslpath -w "$REPORT")" >/dev/null 2>&1
elif [ "$(uname)" = Darwin ]; then
  open "$OUT/results"; [ -n "$REPORT" ] && open "$REPORT"
elif command -v xdg-open >/dev/null 2>&1 && [ -n "${DISPLAY:-}${WAYLAND_DISPLAY:-}" ]; then
  [ -n "$REPORT" ] && xdg-open "$REPORT" >/dev/null 2>&1
fi
finish 0
