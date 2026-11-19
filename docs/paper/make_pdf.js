// HTML application note to PDF with page numbers: node make_pdf.js note.html note.pdf
// Page margins follow the @page rule of the note (17/16/18/16 mm); the page number "n / N" sits
// centred in the bottom margin and never overlaps the text. Chromium draws the footer outside the
// page document, so the IBM Plex Mono face embedded in the note is passed to the footer as well.
const fs = require('fs');
const path = require('path');
const { chromium } = require('playwright');
(async () => {
  const src = path.resolve(process.argv[2]);
  const html = fs.readFileSync(src, 'utf8');
  const face = (html.match(/@font-face\{font-family:'IBM Plex Mono';font-weight:400;font-style:normal;[^}]*U\+0000-00FF[^}]*\}/) || [''])[0];
  const b = await chromium.launch();
  const p = await b.newPage();
  await p.goto('file://' + src);
  await p.pdf({
    path: process.argv[3],
    format: 'A4',
    printBackground: true,
    displayHeaderFooter: true,
    headerTemplate: '<span></span>',
    footerTemplate:
      `<style>${face}</style>` +
      '<div style="width:100%;text-align:center;font-family:\'IBM Plex Mono\',monospace;font-size:8px;color:#56544C;">' +
      '<span class="pageNumber"></span> / <span class="totalPages"></span></div>',
    margin: { top: '17mm', bottom: '18mm', left: '16mm', right: '16mm' },
  });
  await b.close();
})();
