#!/usr/bin/env bash
# Start file of brca12-wgs 2.0 (start/brca12-wgs-2-start.command) with a stand-in program: English and
# German, paths with spaces as dragged into Terminal, reference asked once and remembered, FASTQ pair,
# result per gene; disk dropout: the watchdog stops the run, the start file waits for the disk and
# starts again.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT
PKG="$WORK/my pkg"; mkdir -p "$PKG" "$WORK/home/Downloads" "$WORK/my data"
cp "$ROOT/start/brca12-wgs-2-start.command" "$PKG/start.command"
DATA="$WORK/my data"
printf '>chr13\nN\n' > "$DATA/ref.fa"; printf 'chr13\t1\t7\t1\t2\n' > "$DATA/ref.fa.fai"
touch "$DATA/S 1.bam" "$DATA/S 1.bam.bai" "$DATA/P_R1.fastq.gz" "$DATA/P_R2.fastq.gz"
# stand-in program: logs like brca12-wgs, writes PREFIX.tsv and PREFIX.html; with FAKE_HANG it hangs
# on the first run (as on a lost disk) until it is stopped
cat > "$PKG/brca12-wgs" <<'F'
#!/usr/bin/env bash
args=("$@"); n=${#args[@]}
for ((i = 0; i < n; i++)); do [ "${args[$i]}" = --reference ] && { pre="${args[$((i - 1))]}"; break; }; done
[ -n "${FAKE_HANG:-}" ] && [ ! -e "$FAKE_HANG" ] && { touch "$FAKE_HANG"; sleep 60; exit 1; }
case "$1" in *.fastq.gz) echo "brca12-wgs: 1000 read pairs read, 42 aligned to BRCA1/BRCA2 and the control windows (1 s)" ;; esac
echo "args: $*"
printf 'gene\ttranscript\tresult\tsmall_variants\tclinvar_plp\tloss_of_function\texon_copy_number\tcnv_candidates\tstructural_variants\tcoding_depth\tpct20x\tuncertain_exons\tlow_coverage_bp\tstatement\n' > "$pre.tsv"
printf 'BRCA1\tNM_007294.4\tnone\t0\t0\t0\tnone\t-\tnone\t30\t99\t-\t0\tBRCA1: no pathogenic variant found\n' >> "$pre.tsv"
printf 'BRCA2\tNM_000059.4\tnone\t0\t0\t0\tnone\t-\tnone\t30\t99\t-\t0\tBRCA2: no pathogenic variant found\n' >> "$pre.tsv"
echo "<html></html>" > "$pre.html"
F
chmod +x "$PKG/brca12-wgs"
run() { env HOME="$WORK/home" LANG=C LC_ALL=C BRCA12_NO_OPEN=1 "$@" bash "$PKG/start.command"; }

# first start: reference asked (dragged, quoted), then the BAM (backslash escapes), English
out="$(printf '%s\n%s\n\n\n' "'$DATA/ref.fa'" "$(printf '%s' "$DATA/S 1.bam" | sed 's/ /\\ /g')" | run BRCA12_LANG=en)"
grep -q "research use only (RUO)" <<<"$out" && grep -q "Sample" <<<"$out" && grep -q "Done" <<<"$out" \
  && grep -q "BRCA1  none" <<<"$out" && ! grep -q "Probe\|Fertig" <<<"$out" || { echo "$out"; echo "English run failed"; exit 1; }
[ "$(cat "$WORK/home/.brca12-wgs/reference")" = "$DATA/ref.fa" ] || { echo "reference not remembered"; exit 1; }
[ -s "$WORK/home/Downloads/brca12-wgs-S 1/S 1.tsv" ] && [ -s "$WORK/home/Downloads/brca12-wgs-S 1/S 1.html" ] || { ls -R "$WORK/home/Downloads"; echo "results missing"; exit 1; }
echo "2.0 start file (English, reference asked once, paths with spaces) OK"

# second start: reference remembered, FASTQ pair, German; a new output folder next to the first
out="$(printf '%s\n%s\n\n\n' "$DATA/P_R1.fastq.gz" "$DATA/P_R2.fastq.gz" | run BRCA12_LANG=de)"
grep -q "nur für Forschung (RUO)" <<<"$out" && grep -q "Reads ausgerichtet" <<<"$out" && grep -q "Fertig" <<<"$out" \
  && ! grep -q "GRCh38-FASTA" <<<"$out" || { echo "$out"; echo "German FASTQ run failed"; exit 1; }
grep -q "args: $DATA/P_R1.fastq.gz $DATA/P_R2.fastq.gz $WORK/home/Downloads/brca12-wgs-P/P --reference" "$WORK/home/Downloads/brca12-wgs-P/P.log" || { cat "$WORK/home/Downloads/brca12-wgs-P/P.log"; echo "FASTQ pair not passed"; exit 1; }
out="$(printf '%s\n%s\n\n\n' "$DATA/P_R1.fastq.gz" "$DATA/P_R2.fastq.gz" | run BRCA12_LANG=de)"
[ -s "$WORK/home/Downloads/brca12-wgs-P-2/P.tsv" ] || { echo "earlier result overwritten"; exit 1; }
echo "2.0 start file (German, FASTQ pair, earlier results kept) OK"

# disk dropout: the program hangs, the input disappears, the watchdog stops the run, the disk comes
# back, the run starts again and finishes
printf '%s\n\n\n' "$DATA/S 1.bam" | run BRCA12_LANG=en BRCA12_STALL=1 FAKE_HANG="$WORK/hung" > "$WORK/drop.out" 2>&1 &
RUNPID=$!
for _ in $(seq 100); do [ -e "$WORK/hung" ] && break; sleep 0.1; done
mv "$DATA/S 1.bam" "$DATA/S 1.bam.away"
for _ in $(seq 100); do grep -q "Disk not reachable" "$WORK/drop.out" && break; sleep 0.1; done
grep -q "Disk not reachable, waiting to resume" "$WORK/drop.out" || { kill $RUNPID; cat "$WORK/drop.out"; echo "watchdog did not stop the hanging run"; exit 1; }
mv "$DATA/S 1.bam.away" "$DATA/S 1.bam"
wait $RUNPID || { cat "$WORK/drop.out"; echo "run after disk dropout failed"; exit 1; }
grep -q "Disk back, the run starts again" "$WORK/drop.out" && grep -q "Done" "$WORK/drop.out" || { cat "$WORK/drop.out"; echo "no restart after disk dropout"; exit 1; }
echo "2.0 start file (disk dropout: watchdog stops the run, restart once the disk is back) OK"
