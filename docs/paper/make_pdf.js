// HTML report to PDF with page numbers: node make_pdf.js report.html report.pdf
const { chromium } = require('playwright');
(async () => {
  const b = await chromium.launch();
  const p = await b.newPage();
  const file = require('path').resolve(process.argv[2]);
  await p.goto('file://' + file);
  // the footer is a document of its own: give it the report's embedded IBM Plex Mono (weight 400)
  const html = require('fs').readFileSync(file, 'utf8');
  const face = (html.match(/@font-face\s*{[^}]*IBM Plex Mono[^}]*font-weight:\s*400[^}]*font-style:\s*normal[^}]*}/g) || []).join('');
  await p.pdf({
    path: process.argv[3],
    format: 'A4',
    printBackground: true,
    displayHeaderFooter: true,
    headerTemplate: '<span></span>',
    footerTemplate:
      '<style>' + face + '</style>' +
      '<div style="width:100%;text-align:center;font-family:\'IBM Plex Mono\',monospace;font-size:8px;color:#555;">' +
      '<span class="pageNumber"></span> / <span class="totalPages"></span></div>',
    margin: { top: '17mm', bottom: '18mm', left: '16mm', right: '16mm' },
  });
  await b.close();
})();
