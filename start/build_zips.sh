#!/usr/bin/env bash
# Build the start packages for macOS, Linux and Windows from the committed state (git archive):
#   start/build_zips.sh [OUTDIR]   ->  OUTDIR/brca12-wgs-{macos,linux,windows}.zip
# The names carry no version, so that .../releases/latest/download/brca12-wgs-macos.zip always points
# to the newest release (the version is the release tag). Each ZIP holds one folder "brca12-wgs" with the start file of that system, README.txt
# (English, then German) and the pipeline in "pipeline/" (without tests, validation data and paper).
# The version comes from nextflow.config. Line endings: .bat CRLF, everything else as committed (LF).
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
out="$(mkdir -p "${1:-$root/dist}" && cd "${1:-$root/dist}" && pwd)"
version="$(sed -n "s/^ *version *= *'\([^']*\)'.*/\1/p" "$root/nextflow.config" | head -1)"
[ -n "$version" ] || { echo "version not found in nextflow.config" >&2; exit 1; }
work="$(mktemp -d)"; trap 'rm -rf "$work"' EXIT
git -C "$root" archive --format=tar HEAD -- main.nf nextflow.config bin assets envs \
    LICENSE README.md CHANGELOG.md CITATION.cff Dockerfile docs/output.md start | tar -x -C "$work"
for os in macos linux windows; do
  name="brca12-wgs"
  dir="$work/$os/$name"
  mkdir -p "$dir/pipeline"
  ( cd "$work" && tar -c --exclude=start --exclude="./macos" --exclude="./linux" --exclude="./windows" . ) | tar -x -C "$dir/pipeline"
  cp "$work/start/README.txt" "$dir/README.txt"
  case "$os" in
    macos)   cp "$work/start/brca12-wgs-start.command" "$dir/brca12-wgs-start.command"; chmod 755 "$dir/brca12-wgs-start.command" ;;
    linux)   cp "$work/start/brca12-wgs-start.command" "$dir/brca12-wgs-start.sh"; chmod 755 "$dir/brca12-wgs-start.sh" ;;
    windows) cp "$work/start/brca12-wgs-start.command" "$dir/brca12-wgs-start.sh"; chmod 755 "$dir/brca12-wgs-start.sh"
             cp "$work/start/brca12-wgs-start.bat" "$dir/brca12-wgs-start.bat" ;;
  esac
  rm -f "$out/$name-$os.zip"
  ( cd "$work/$os" && zip -q -r -X "$out/$name-$os.zip" "$name" )
  echo "$out/$name-$os.zip ($version)"
done
