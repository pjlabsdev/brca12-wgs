#!/usr/bin/env bash
# Print the HTML report to PDF with headless Chrome/Chromium.
#   bin/brca12_pdf.sh brca12_report.html brca12_report.pdf [CHROME]
# Without a usable Chrome/Chromium the step is skipped with a note (exit 0, no PDF).
# Chrome's own header and footer are kept: title (sample and RUO note), date and page numbers.
set -euo pipefail
html="$1"; pdf="$2"; chrome="${3:-}"
if [ -z "$chrome" ]; then
  for candidate in chromium chromium-browser google-chrome google-chrome-stable chrome-headless-shell \
      "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" \
      "/Applications/Chromium.app/Contents/MacOS/Chromium"; do
    if command -v "$candidate" >/dev/null 2>&1; then chrome="$(command -v "$candidate")"; break; fi
    if [ -x "$candidate" ]; then chrome="$candidate"; break; fi
  done
fi
if [ -z "$chrome" ] || ! [ -x "$chrome" ] && ! command -v "$chrome" >/dev/null 2>&1; then
  echo "[brca12] PDF report skipped: no Chrome/Chromium found (set --chrome /path/to/chrome)" >&2
  exit 0
fi
abs="$(cd "$(dirname "$html")" && pwd)/$(basename "$html")"
profile="$(mktemp -d)"
trap 'rm -rf "$profile"' EXIT
"$chrome" --headless --disable-gpu --no-sandbox --no-first-run --user-data-dir="$profile" \
  --print-to-pdf="$pdf" "file://$abs" >/dev/null 2>&1 || true
if [ -s "$pdf" ]; then
  echo "[brca12] PDF report written: $pdf" >&2
else
  echo "[brca12] PDF report skipped: Chrome/Chromium could not print ($chrome)" >&2
  rm -f "$pdf"
fi
