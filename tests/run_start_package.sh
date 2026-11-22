#!/usr/bin/env bash
# Start packages: builds the three ZIPs (start/build_zips.sh), unpacks them into a folder with a
# space and runs the start file from the Linux package: dry run in English and German, a Windows
# path under a stand-in WSL, and one whole run with the stand-in tools (tests/fake_tools).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
bash "$ROOT/start/build_zips.sh" "$WORK/dist" > "$WORK/zips.txt"
[ "$(wc -l < "$WORK/zips.txt")" = 3 ] || { echo "expected three ZIPs"; exit 1; }
for os in macos linux windows; do
  mkdir -p "$WORK/un pack/$os"
  unzip -q "$WORK/dist/brca12-wgs-$os.zip" -d "$WORK/un pack/$os"
  pkg="$WORK/un pack/$os/brca12-wgs"
  for f in README.txt pipeline/main.nf pipeline/nextflow.config pipeline/bin/brca12_report.py \
           pipeline/assets/fonts/OFL.txt pipeline/envs/brca12.yml pipeline/LICENSE; do
    [ -s "$pkg/$f" ] || { echo "$os: missing $f"; exit 1; }
  done
  [ ! -e "$pkg/pipeline/tests" ] && [ ! -e "$pkg/pipeline/validation" ] || { echo "$os: tests or validation data in the package"; exit 1; }
  grep -q "^Deutsch$" "$pkg/README.txt" && [ "$(grep -n -m1 "^Deutsch$" "$pkg/README.txt" | cut -d: -f1)" -gt 20 ] || { echo "$os: README must be English first, then German"; exit 1; }
done
[ -x "$WORK/un pack/macos/brca12-wgs/brca12-wgs-start.command" ] || { echo "macos: start file not executable"; exit 1; }
[ -x "$WORK/un pack/linux/brca12-wgs/brca12-wgs-start.sh" ] || { echo "linux: start file not executable"; exit 1; }
BAT="$WORK/un pack/windows/brca12-wgs/brca12-wgs-start.bat"
[ "$(grep -c $'\r$' "$BAT")" = "$(wc -l < "$BAT")" ] || { echo "windows: .bat needs CRLF line endings"; exit 1; }
grep -q "brca12-wgs-start.sh" "$BAT" && [ -s "$(dirname "$BAT")/brca12-wgs-start.sh" ] || { echo "windows: .bat does not start the shell file"; exit 1; }
if grep -q $'\r' "$(dirname "$BAT")/brca12-wgs-start.sh"; then echo "windows: shell file must keep LF line endings"; exit 1; fi
echo "start packages (three ZIPs, contents) OK"

START="$WORK/un pack/linux/brca12-wgs/brca12-wgs-start.sh"
DATA="$WORK/my data"; mkdir -p "$DATA/vep" "$WORK/home/.config/brca12-wgs" "$WORK/home/Downloads"
printf ">chr13\nN\n>chr17\nN\n" > "$DATA/ref.fa"
printf "chr13\t114364328\t0\t60\t61\nchr17\t83257441\t0\t60\t61\n" > "$DATA/ref.fa.fai"
touch "$DATA/ref.dict" "$DATA/S 1.bam" "$DATA/S 1.bam.bai" "$DATA/clinvar.vcf.gz" "$DATA/clinvar.vcf.gz.tbi"
FAKE="$WORK/fake"; mkdir -p "$FAKE"
for tool in bwa samtools gatk bcftools tabix bgzip vep configManta.py delly cnvkit.py Rscript mosdepth; do
  ln -s "$ROOT/tests/fake_tools/faketool.py" "$FAKE/$tool"
done
{
  printf 'REFERENCE=%q\nVEP_CACHE=%q\nCLINVAR_VCF=%q\nPROFILE=standard\n' "$DATA/ref.fa" "$DATA/vep" "$DATA/clinvar.vcf.gz"
  printf 'EXTRA_ARGS=%q\n' "--exons_bed $ROOT/tests/data/exons.test.bed"
} > "$WORK/home/.config/brca12-wgs/start.conf"
run() { env HOME="$WORK/home" PATH="$FAKE:$PATH" LANG=C LC_ALL=C "$@" bash "$START"; }

out="$(printf '%s\n\n' "$DATA/S 1.bam" | run BRCA12_DRYRUN=1)"
grep -q "Dry run: the run would start here." <<<"$out" && grep -q "research use only (RUO)" <<<"$out" \
  && grep -q "Sample" <<<"$out" && ! grep -q "Probe\|Trockenlauf" <<<"$out" || { echo "$out"; echo "English dry run failed"; exit 1; }
out="$(printf '%s\n\n' "'$DATA/S 1.bam'" | run BRCA12_DRYRUN=1 BRCA12_LANG=de)"
grep -q "Trockenlauf: hier würde der Lauf starten." <<<"$out" && grep -q "nur für Forschung (RUO)" <<<"$out" || { echo "$out"; echo "German dry run failed"; exit 1; }
out="$(printf '%s\n\n' "$DATA/S 1.bam" | run BRCA12_DRYRUN=1 LANG=de_AT.UTF-8 LC_ALL=de_AT.UTF-8)"
grep -q "Trockenlauf" <<<"$out" || { echo "$out"; echo "German system language not used"; exit 1; }
echo "start file (dry run, English first, German second) OK"

# Windows path dragged into the window, translated by a stand-in wslpath; output goes to the Windows Downloads folder
WIN="$WORK/win"; mkdir -p "$WIN" "$WORK/winhome/Downloads"
cat > "$WIN/wslpath" <<W
#!/bin/sh
case "\$1" in -u) case "\$2" in 'C:\\Users\\me') echo "$WORK/winhome" ;; 'D:\\my data\\S 1.bam') echo "$DATA/S 1.bam" ;; *) echo "/mnt/unknown" ;; esac ;; *) echo "\$2" ;; esac
W
printf '#!/bin/sh\nprintf "C:\\\\\\\\Users\\\\\\\\me\\r\\n"\n' > "$WIN/cmd.exe"
printf '#!/bin/sh\nprintf "en-US\\r\\n"\n' > "$WIN/powershell.exe"
chmod +x "$WIN"/*
out="$(printf '%s\n\n' '"D:\my data\S 1.bam"' | PATH="$WIN:$PATH" run BRCA12_DRYRUN=1 BRCA12_WSL=1)"
grep -q "S 1.bam" <<<"$out" && grep -q "$WORK/winhome/Downloads/brca12-S_1" <<<"$out" && grep -q "Dry run" <<<"$out" || { echo "$out"; echo "Windows path under WSL failed"; exit 1; }
echo "start file (Windows path under WSL) OK"

# Update check: a stand-in "latest release" (file:// URLs) that is newer offers the update and replaces
# the program files in place (results next to them stay); an equal version, a failed request and
# BRCA12_NO_UPDATE=1 stay silent; Enter (default no) continues with the installed version.
CONFIG="$(dirname "$START")/pipeline/nextflow.config"
VER=$(awk -F"'" '/^ *version *=/ { print $2; exit }' "$CONFIG")
UP="$WORK/update"; mkdir -p "$UP"
cp -R "$(dirname "$START")" "$UP/brca12-wgs"
sed -i.bak "s/^\( *version *= *\)'[^']*'/\1'99.0.0'/" "$UP/brca12-wgs/pipeline/nextflow.config" && rm "$UP/brca12-wgs/pipeline/nextflow.config.bak"
echo "new" > "$UP/brca12-wgs/pipeline/NEW_FILE"
( cd "$UP" && zip -q -r brca12-wgs-linux.zip brca12-wgs && rm -rf brca12-wgs )
printf '{"tag_name": "v99.0.0", "assets": [{"browser_download_url": "file://%s/brca12-wgs-macos.zip"}, {"browser_download_url": "file://%s/brca12-wgs-linux.zip"}]}\n' "$UP" "$UP" > "$UP/latest.json"
printf '{"tag_name": "v%s", "assets": []}\n' "$VER" > "$UP/same.json"
out="$(printf '%s\n\n' "$DATA/S 1.bam" | run BRCA12_DRYRUN=1 BRCA12_UPDATE_URL="file://$UP/same.json")"
! grep -q "New release" <<<"$out" && grep -q "Dry run" <<<"$out" || { echo "$out"; echo "update offered for the same version"; exit 1; }
out="$(printf '%s\n\n' "$DATA/S 1.bam" | run BRCA12_DRYRUN=1 BRCA12_UPDATE_URL="file://$UP/missing.json")"
! grep -q "New release" <<<"$out" && grep -q "Dry run" <<<"$out" || { echo "$out"; echo "failed update request must stay silent"; exit 1; }
out="$(printf '%s\n\n' "$DATA/S 1.bam" | run BRCA12_DRYRUN=1 BRCA12_NO_UPDATE=1 BRCA12_UPDATE_URL="file://$UP/latest.json")"
! grep -q "New release" <<<"$out" || { echo "$out"; echo "BRCA12_NO_UPDATE=1 ignored"; exit 1; }
out="$(printf '\n%s\n\n' "$DATA/S 1.bam" | run BRCA12_DRYRUN=1 BRCA12_UPDATE_URL="file://$UP/latest.json")"
grep -q "New release 99.0.0 available" <<<"$out" && grep -q "Dry run" <<<"$out" && ! grep -q "99.0.0" "$(dirname "$START")/pipeline/nextflow.config" \
  || { echo "$out"; echo "declined update must continue with the installed version"; exit 1; }
mkdir -p "$(dirname "$START")/brca12-results-kept"
out="$(printf 'y\n%s\n\n' "$DATA/S 1.bam" | run BRCA12_DRYRUN=1 BRCA12_UPDATE_URL="file://$UP/latest.json")"
D="$(dirname "$START")"
grep -q "Updated to 99.0.0" <<<"$out" && grep -q "Dry run" <<<"$out" && grep -q "99.0.0" "$D/pipeline/nextflow.config" \
  && [ -f "$D/pipeline/NEW_FILE" ] && [ -x "$START" ] && [ -d "$D/brca12-results-kept" ] && [ ! -e "$D/pipeline.old" ] && [ ! -e "$D/pipeline.update" ] \
  || { echo "$out"; ls -la "$D"; echo "in-place update failed"; exit 1; }
out="$(printf '%s\n\n' "$DATA/S 1.bam" | run BRCA12_DRYRUN=1 BRCA12_UPDATE_URL="file://$UP/latest.json")"
! grep -q "New release" <<<"$out" || { echo "$out"; echo "updated package still offers the update"; exit 1; }
rm -rf "$UP"
# fresh copy of the package for the whole run below
rm -rf "$WORK/un pack/linux"; mkdir -p "$WORK/un pack/linux"; unzip -q "$WORK/dist/brca12-wgs-linux.zip" -d "$WORK/un pack/linux"
echo "start file (update check: offer, decline, install, silent cases) OK"

if ! command -v nextflow >/dev/null 2>&1; then echo "start file (whole run) skipped: no nextflow"; exit 0; fi
out="$(printf '%s\n\n\n\n' "$DATA/S 1.bam" | run)" || { echo "$out"; tail -30 "$WORK/home/Downloads/brca12-S_1/brca12-start.log" 2>/dev/null; echo "run from the start package failed"; exit 1; }
grep -q "Done" <<<"$out" || { echo "$out"; exit 1; }
for f in brca12_report.html brca12_results.json; do
  test -s "$WORK/home/Downloads/brca12-S_1/results/$f" || { echo "missing $f"; exit 1; }
done
grep -q "<body data-lang='en' data-tab='r'>" "$WORK/home/Downloads/brca12-S_1/results/brca12_report.html" || { echo "report does not show English first"; exit 1; }
echo "start file (whole run from the unpacked package, paths with spaces) OK"
