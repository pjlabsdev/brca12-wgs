#!/bin/bash
# brca12-wgs 2.0 Start: double-click (macOS: .command opens in Terminal) or run `bash brca12-wgs-start.command`.
# Asks for the input (BAM, CRAM or FASTQ pair) and, once, for the GRCh38 reference; drag the files into
# the window (Terminal inserts their paths) or type them. Runs brca12-wgs with a progress display and
# opens the result page. Built on the hla-wgs start file (same display and input handling).
# Disk dropout (an external disk that disconnects under load): if the input or the reference is gone
# for BRCA12_STALL seconds (default 120) the run is stopped; once the disk is back it starts again
# (up to BRCA12_RESTARTS times, default 3). Short I/O errors on FASTQ input the program bridges itself.
set -u
case "${LC_ALL:-${LC_CTYPE:-${LANG:-}}}" in *UTF-8*) ;; *) export LC_CTYPE=en_US.UTF-8 ;; esac
HERE="$(cd "$(dirname "$0")" && pwd)"
CONF="$HOME/.brca12-wgs"

if [ -n "${BRCA12_LANG:-}" ]; then [ "$BRCA12_LANG" = de ] && DE=1 || DE=0
elif [ "${LANG:-}" != "${LANG#de}" ] || [ "$(defaults read -g AppleLocale 2>/dev/null | cut -c1-2)" = "de" ]; then DE=1; else DE=0; fi
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
secs() { awk -v a="$1" -v b="$(now)" 'BEGIN { printf "%.1fs", b - a }'; }
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
group() { local sep=","; [ "$DE" = 1 ] && sep="."; printf '%s' "$1" | perl -pe "1 while s/^(\\d+)(\\d{3})/\$1$sep\$2/"; }
ask() { printf '\n  %s\n  %s>%s ' "$(tr_ "$1" "$2")" "$ACC" "$RST"; }
finish() { cursor h; printf '\n'; printf '  %s%s%s' "$MUT" "$(tr_ "Enter zum Schließen" "Press Enter to close")" "$RST"; read -r _; exit "${1:-0}"; }
trap 'cursor h' EXIT
trap 'cursor h; [ -n "$PID" ] && kill "$PID" 2>/dev/null; exit 130' INT TERM HUP

# runs a command in the background with a spinner: step LABEL COMMAND …  (sets T0 for the result line)
step() {
  local label=$1; shift
  T0=$(now)
  if [ "$TTY" = 1 ]; then
    "$@" >"$WORKLOG" 2>&1 & PID=$!; local i=0
    cursor l
    while kill -0 "$PID" 2>/dev/null; do spin "$i" "$label"; i=$((i + 1)); sleep 0.08; done
    wait "$PID"; local rc=$?; PID=""; clearline; return $rc
  else
    "$@" >"$WORKLOG" 2>&1
  fi
}

# download with progress bar: download URL TARGET
download() {
  T0=$(now)
  if [ "$TTY" = 1 ]; then
    local total; total=$(curl -fsSLI "$1" 2>/dev/null | awk 'tolower($1)=="content-length:"{v=$2} END{gsub("\r","",v); print v+0}')
    curl -fsSL -o "$2" "$1" 2>/dev/null & PID=$!; local i=0 size
    cursor l
    while kill -0 "$PID" 2>/dev/null; do
      size=$(stat -f %z "$2" 2>/dev/null || stat -c %s "$2" 2>/dev/null || echo 0)
      if [ "${total:-0}" -gt 0 ]; then
        spin "$i" "Download  $(bar "$size" "$total" 24)  $MUT$(awk -v s="$size" -v t="$total" 'BEGIN{printf "%5.1f / %.1f MB", s/1048576, t/1048576}')$RST"
      else
        spin "$i" "Download  $MUT$(awk -v s="$size" 'BEGIN{printf "%.1f MB", s/1048576}')$RST"
      fi
      i=$((i + 1)); sleep 0.1
    done
    wait "$PID"; local rc=$?; PID=""; clearline; return $rc
  else
    curl -fsSL -o "$2" "$1"
  fi
}

# a path as dragged into Terminal: quotes or backslash escapes, trailing space
clean() {
  local p="$1"
  p="${p#"${p%%[![:space:]]*}"}"; p="${p%"${p##*[![:space:]]}"}"
  case "$p" in \'*\') p="${p#\'}"; p="${p%\'}" ;; \"*\") p="${p#\"}"; p="${p%\"}" ;; *) p="$(printf '%s' "$p" | sed 's/\\\(.\)/\1/g')" ;; esac
  case "$p" in "~/"*) p="$HOME/${p#\~/}" ;; esac
  printf '%s' "$p"
}

askpath() {  # $1 de, $2 en, $3 allow empty
  local p
  while true; do
    ask "$1" "$2"; IFS= read -r p || finish 1
    p="$(clean "$p")"
    if [ -z "$p" ] && [ "${3:-}" = empty ]; then REPLY_PATH=""; return; fi
    if [ -f "$p" ]; then REPLY_PATH="$p"; return; fi
    printf '  %s%s %s%s\n' "$ERR" "$(tr_ "Datei nicht gefunden:" "File not found:")" "$p" "$RST"
  done
}

# brca12-wgs program: next to this file, in the repository build, or in PATH
BIN=""
for c in "$HERE/brca12-wgs" "$HERE/bin/brca12-wgs" "$HERE/../target/release/brca12-wgs" "$(command -v brca12-wgs 2>/dev/null)"; do
  if [ -n "$c" ] && [ -x "$c" ]; then BIN="$c"; break; fi
done
WORKLOG="$(mktemp -t brca12-wgs-start.XXXXXX)"
# a ZIP unpacked after a browser download marks the program as quarantined, and macOS then stops it
# at once; this file was started on purpose, so the mark is removed from the program only
[ -n "$BIN" ] && command -v xattr >/dev/null 2>&1 && xattr -d com.apple.quarantine "$BIN" 2>/dev/null

# ── header ───────────────────────────────────────────────────────────────────
clear 2>/dev/null
title "brca12-wgs"
printf '\n'
row 2 "brca12-wgs" "$(date '+%d.%m.%Y  %H:%M')" "$ACC" "$MUT"
printf '  %s%s%s\n' "$MUT" "$(tr_ "BRCA1 und BRCA2 aus dem Genom · nur für Forschung (RUO)" "BRCA1 and BRCA2 from the genome · research use only (RUO)")" "$RST"
divider
if [ -z "$BIN" ]; then
  printf '\n  %s%s%s\n' "$ERR" "$(tr_ "Das Programm brca12-wgs wurde nicht gefunden (neben dieser Datei, in ../target/release oder im PATH)." "The brca12-wgs program was not found (next to this file, in ../target/release or in PATH).")" "$RST"
  finish 1
fi

# ── reference: BRCA12_WGS_REF, else the one remembered from the last start, else asked once ──
REF="${BRCA12_WGS_REF:-$(cat "$CONF/reference" 2>/dev/null)}"
if [ -z "$REF" ] || [ ! -f "$REF" ] || [ ! -f "$REF.fai" ]; then
  printf '\n'
  row 2 "$(tr_ "Einmalige Einrichtung" "One-time setup")" "GRCh38" "$ACC"
  printf '  %s%s%s\n' "$MUT" "$(tr_ "Gebraucht wird das Referenzgenom GRCh38 als FASTA mit Index (.fai daneben), mit dem die Daten ausgerichtet wurden." "Needed is the GRCh38 reference genome as FASTA with its index (.fai next to it), the one the data were aligned to.")" "$RST"
  while true; do
    askpath "GRCh38-FASTA hineinziehen oder Pfad eingeben:" "Drop the GRCh38 FASTA or enter its path:"
    REF="$REPLY_PATH"
    [ -f "$REF.fai" ] && break
    printf '  %s%s %s.fai%s\n' "$ERR" "$(tr_ "Der Index fehlt:" "The index is missing:")" "$REF" "$RST"
  done
  mkdir -p "$CONF" && printf '%s\n' "$REF" > "$CONF/reference"
fi

# ── input ────────────────────────────────────────────────────────────────────
askpath "BAM, CRAM oder FASTQ (R1) hier hineinziehen oder Pfad eingeben, dann Enter:" "Drop a BAM, CRAM or FASTQ (R1) here or enter its path, then press Enter:"
IN="$REPLY_PATH"
lower="$(printf '%s' "$IN" | tr '[:upper:]' '[:lower:]')"
IN2=""
stop() { printf '  %s%s%s\n' "$ERR" "$1" "$RST"; finish 1; }
case "$lower" in
  *.bam)
    if [ ! -f "$IN.bai" ] && [ ! -f "${IN%.*}.bai" ] && [ ! -f "$IN.csi" ]; then stop "$(tr_ "Zur BAM fehlt der Index (.bai daneben)." "The BAM index (.bai) is missing.")"; fi ;;
  *.cram)
    [ -f "$IN.crai" ] || stop "$(tr_ "Zur CRAM fehlt der Index (.crai daneben)." "The CRAM index (.crai) is missing.")" ;;
  *.fq|*.fastq|*.fq.gz|*.fastq.gz)
    askpath "Zweite FASTQ-Datei (R2) hineinziehen:" "Drop the second FASTQ file (R2):"; IN2="$REPLY_PATH"
    printf '  %s%s%s\n' "$MUT" "$(tr_ "FASTQ: Das Programm liest alle Reads des Genoms und richtet die von BRCA1/BRCA2 selbst aus (bei 30x etwa 6 Minuten). Weitere Programme braucht es nicht." "FASTQ: the program reads every read of the genome and aligns those of BRCA1/BRCA2 itself (about 6 minutes at 30x). No other programs are needed.")" "$RST" ;;
  *) stop "$(tr_ "Unbekanntes Format. Erlaubt sind BAM, CRAM oder FASTQ." "Unknown format. Allowed are BAM, CRAM or FASTQ.")" ;;
esac

# sample name and an output folder that holds no earlier result
S="$(basename "$IN")"
for x in .fastq.gz .fq.gz .fastq .fq .bam .cram .BAM .CRAM; do S="${S%"$x"}"; done
case "$lower" in *.fq|*.fastq|*.fq.gz|*.fastq.gz) S="$(printf '%s' "$S" | sed -E 's/_S[0-9]+_L[0-9]{3}_R1_001$//; s/_L[0-9]{3}_R1_001$//; s/_R1_001$//; s/([._])[Rr]?1$//')" ;; esac
[ -z "$S" ] && S=sample
OUT="$HOME/Downloads/brca12-wgs-$S"; n=2
while [ -e "$OUT/$S.tsv" ]; do OUT="$HOME/Downloads/brca12-wgs-$S-$n"; n=$((n + 1)); done
T="$( (sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4) | head -1)"
SHOWOUT="$OUT"; case "$OUT" in "$HOME"/*) SHOWOUT="~/${OUT#"$HOME"/}" ;; esac
SHOWREF="$REF"; case "$REF" in "$HOME"/*) SHOWREF="~/${REF#"$HOME"/}" ;; esac

printf '\n'
divider
row 2 "$(tr_ "Probe" "Sample")" "$S" "$MUT" ""
row 2 "$(tr_ "Eingabe" "Input")" "$(basename "$IN")${IN2:+ + $(basename "$IN2")}" "$MUT" ""
row 2 "$(tr_ "Referenz" "Reference")" "$(basename "$REF")" "$MUT" ""
row 2 "$(tr_ "Ausgabe" "Output")" "$SHOWOUT" "$MUT" ""
row 2 "Threads" "$T" "$MUT" ""
divider
ask "Start mit Enter (oder Strg-C zum Abbrechen)" "Press Enter to start (or Ctrl-C to cancel)"; read -r _ || exit 1
mkdir -p "$OUT" || finish 1
LOG="$OUT/$S.log"

# ── run with progress ────────────────────────────────────────────────────────
printf '\n'
title "brca12-wgs · $S"
[ "$DE" = 1 ] && export BRCA12_WGS_LANG=de   # result page PREFIX.html in German first
STALL="${BRCA12_STALL:-120}"; RESTARTS="${BRCA12_RESTARTS:-3}"
reachable() { [ -f "$IN" ] && [ -f "$REF" ] && { [ -z "$IN2" ] || [ -f "$IN2" ]; }; }
L_READ="$(tr_ "Reads lesen und ausrichten" "Reading and aligning reads")"; L_READ_OK="$(tr_ "Reads ausgerichtet" "Reads aligned")"
L_RUN="$(tr_ "Analyse" "Analysis")"
START=$(now); attempt=0
while true; do
  if [ -n "$IN2" ]; then
    "$BIN" "$IN" "$IN2" "$OUT/$S" --reference "$REF" --threads "$T" >"$LOG" 2>&1 & PID=$!
  else
    "$BIN" "$IN" "$OUT/$S" --reference "$REF" >"$LOG" 2>&1 & PID=$!
  fi
  # phases from the log: FASTQ reading and alignment ("read pairs read"), then the analysis
  phase=2; [ -n "$IN2" ] && phase=1
  i=0; T0=$(now); gone=""; stopped=0; cursor l
  aligned() {  # the alignment line once it is in the log
    [ "$phase" = 1 ] && grep -q 'read pairs read' "$LOG" 2>/dev/null || return 0
    clearline
    pairs="$(sed -n 's/.*: \([0-9][0-9]*\) read pairs read, \([0-9][0-9]*\) aligned.*/\2/p' "$LOG" | head -1)"
    result ok "$L_READ_OK" "${pairs:+$(group "$pairs") $(tr_ "Paare" "pairs")}" "$(secs "$T0")"
    phase=2; T0=$(now)
  }
  while kill -0 "$PID" 2>/dev/null; do
    aligned
    # watchdog: input or reference unreachable for STALL seconds stops the run
    if reachable; then gone=""
    elif [ -z "$gone" ]; then gone=$(date +%s)
    elif [ $(( $(date +%s) - gone )) -ge "$STALL" ]; then kill "$PID" 2>/dev/null; stopped=1; break
    fi
    if [ "$phase" = 1 ]; then spin "$i" "$L_READ  $MUT$(secs "$T0")$RST"
    else spin "$i" "$L_RUN  $MUT$(tr_ "Abdeckung, Varianten, Kopienzahl, Strukturvarianten" "coverage, variants, copy number, structural variants")  $(secs "$T0")$RST"; fi
    i=$((i + 1)); sleep 0.1
  done
  wait "$PID"; rc=$?; PID=""; clearline; cursor h
  aligned
  [ "$rc" = 0 ] && break
  # a failed or stopped run with the disk gone: wait for it, then start again
  if [ "$stopped" = 1 ] || ! reachable; then
    attempt=$((attempt + 1))
    [ "$attempt" -gt "$RESTARTS" ] && break
    result fail "$(tr_ "Platte nicht erreichbar, warte auf Fortsetzung" "Disk not reachable, waiting to resume")" "" "$(secs "$T0")"
    while ! reachable; do sleep 2; done
    result ok "$(tr_ "Platte wieder da, Lauf startet neu" "Disk back, the run starts again")" "$attempt/$RESTARTS" ""
    continue
  fi
  break
done
if [ "$rc" = 0 ]; then
  result ok "$L_RUN" "" "$(secs "$T0")"
  [ -f "$OUT/$S.html" ] && result ok "$(tr_ "Ergebnisseite" "Result page")" "$S.html" ""
fi

printf '\n'
divider
if [ "$rc" != 0 ]; then
  row 2 "$(tr_ "Fehler $rc" "Error $rc")" "$(secs "$START")" "$ERR"
  tail -3 "$LOG" | sed "s/^/    $MUT/; s/\$/$RST/"
  printf '    %s%s %s%s\n' "$MUT" "$(tr_ "Protokoll:" "Log:")" "$LOG" "$RST"
  finish "$rc"
fi
row 2 "$(tr_ "Fertig" "Done")" "$(secs "$START")" "$OK"
printf '\n'
# one line per gene: result and the statement (English in PREFIX.tsv; German on the page)
awk -F'\t' -v g="$OK" -v a="$ACC" -v e="$ERR" -v r="$RST" 'NR > 1 { c = g; if ($3 == "pathogenic") c = e; else if ($3 == "candidate" || $3 == "not_assessable") c = a
  printf "  %s  %s%s%s\n      %s\n\n", $1, c, $3, r, $14 }' "$OUT/$S.tsv"
printf '  %s%s%s\n' "$MUT" "$(tr_ "Nur für Forschung (RUO). Jeder Befund braucht eine Bestätigung in einem akkreditierten Labor." "Research use only (RUO). Every finding needs confirmation in an accredited laboratory.")" "$RST"
printf '  %s%s%s\n' "$MUT" "$SHOWOUT" "$RST"
title "brca12-wgs · $S · $(tr_ "fertig" "done")"
if [ -z "${BRCA12_NO_OPEN:-}" ] && command -v open >/dev/null 2>&1; then open "$OUT"; [ -f "$OUT/$S.html" ] && open "$OUT/$S.html"; fi
finish 0
