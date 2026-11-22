"""HTML result page, built after the hla-wgs result page: one interactive file with English and German
(switch at the top), two tabs (result, details per gene), explanations on hover or tap, report number,
overall statement, and every value of brca12_results.json. Printing from the browser prints the open tab
with all genes expanded. The font (IBM Plex Mono, SIL OFL 1.1, see assets/fonts/OFL.txt) is embedded."""
import base64
import html
import json
import os

from . import figures
from .i18n import tr, ui

DOI = "10.5281/zenodo.24173638"
RESEARCH_URL = "https://pjlabs.dev/forschung/brca12-wgs"
SOURCE_URL = "https://github.com/pjlabsdev/brca12-wgs"
FONT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "assets", "fonts")
GENES = ("BRCA1", "BRCA2")

TEXT = {
    "ruo": ("Nur für Forschung (RUO). Kein medizinischer Befund und kein Medizinprodukt.",
            "Research use only (RUO). Not a medical report and not a medical device."),
    "eyebrow": ("BRCA1/BRCA2 aus dem Genom", "BRCA1/BRCA2 from the genome"),
    "tab_r": ("Ergebnis", "Result"),
    "tab_d": ("Details je Gen", "Details per gene"),
    "copy": ("Kopieren", "Copy"),
    "sample": ("Probe", "Sample"),
    "sample_tip": ("Name der Probe (--sample).", "Sample name (--sample)."),
    "report_no": ("Berichtsnummer", "Report number"),
    "report_no_tip": ("Prüfsumme (FNV-1a) über alle Ergebnisse und die Namen der Eingabedateien: gleiches Ergebnis, gleiche Nummer; jede Änderung ergibt eine andere.",
                      "Checksum (FNV-1a) over all results and the names of the input files: same result, same number; any change gives another one."),
    "input": ("Eingabe", "Input"),
    "input_tip": ("Eingabedateien ohne Ordner. Art: {mode}. Referenz: {ref} ({assembly}).",
                  "Input files without folder. Kind: {mode}. Reference: {ref} ({assembly})."),
    "analysis": ("Analyse", "Analysis"),
    "analysis_tip": ("Zeitpunkt des Berichts (UTC): {ts}", "Time of the report (UTC): {ts}"),
    "program": ("Programm", "Program"),
    "program_tip": ("brca12-wgs {version} · DOI {doi} · SNV/Indel-Rufer {caller} · Nextflow-Sitzung {session}",
                    "brca12-wgs {version} · DOI {doi} · SNV/indel caller {caller} · Nextflow session {session}"),
    "database": ("Datenbanken", "Databases"),
    "database_tip": ("Annotation mit Ensembl VEP {vep} und ClinVar (Version der verwendeten ClinVar-VCF). ACMG/AMP-Hinweise nach {spec}. gnomAD: {gnomad}. SpliceAI: {spliceai}. BayesDel: {bayesdel}.",
                     "Annotation with Ensembl VEP {vep} and ClinVar (version of the ClinVar VCF used). ACMG/AMP hints after {spec}. gnomAD: {gnomad}. SpliceAI: {spliceai}. BayesDel: {bayesdel}."),
    "yes": ("ja", "yes"), "no": ("nein", "no"),
    "depth": ("Tiefe über den Exons", "Depth over the exons"),
    "depth_tip": ("Mittlere Tiefe aller Exons je Gen (Reads mit MAPQ >= 20). Ein Exon gilt unter {min:g}x oder unter {pct:g} % der Basen ab 20x als unsicher.",
                  "Mean depth of all exons per gene (reads with MAPQ >= 20). An exon counts as uncertain below {min:g}x or below {pct:g}% of bases at 20x."),
    "quality": ("Qualität", "Quality"),
    "quality_tip": ("Bestanden, wenn kein Exon unsicher ist, die mittlere Exon-Tiefe beider Gene mindestens {min:g}x beträgt und die Read-Prüfung keinen Hinweis meldet. Sonst eingeschränkt, mit Grund.",
                    "Passed if no exon is uncertain, the mean exon depth of both genes is at least {min:g}x and the read check reports no note. Otherwise limited, with the reason."),
    "passed": ("bestanden", "passed"),
    "limited": ("eingeschränkt", "limited"),
    "q_exons": ("unsichere Exons: {n}", "uncertain exons: {n}"),
    "q_depth": ("{gene}: mittlere Exon-Tiefe {depth:g}x", "{gene}: mean exon depth {depth:g}x"),
    "q_reads": ("Hinweise der Read-Prüfung: {n}", "read-check notes: {n}"),
    "verdict_h": ("Gesamtaussage", "Overall statement"),
    "v_none": ("BRCA1 und BRCA2: In den geprüften Kategorien gibt es keinen Hinweis auf eine erbliche Veranlagung.",
               "BRCA1 and BRCA2: no sign of a hereditary predisposition in the categories checked."),
    "v_unclear": ("BRCA1 und BRCA2: kein Befund, der in ClinVar als krankheitsverursachend gilt; {n} Einträge mit unklarer Bedeutung.",
                  "BRCA1 and BRCA2: no finding that ClinVar lists as disease-causing; {n} entries of uncertain significance."),
    "v_attention": ("BRCA1 und BRCA2: {n} Befunde, die für eine erbliche Veranlagung sprechen können und bestätigt werden müssen.",
                    "BRCA1 and BRCA2: {n} findings that may indicate a hereditary predisposition and need confirmation."),
    "v_checked": ("Geprüft: kleine Varianten mit ClinVar und Hinweis auf Funktionsverlust, Kopienzahl (CNVkit und exonweise), strukturelle Varianten, Einschub-Signale{mosaic}. Ein unauffälliges Ergebnis schließt eine Veranlagung nicht aus; eine medizinische Aussage gibt nur ein akkreditiertes Labor.",
                  "Checked: small variants with ClinVar and loss-of-function hints, copy number (CNVkit and per exon), structural variants, insertion signals{mosaic}. An unremarkable result does not rule out a predisposition; only an accredited laboratory gives a medical statement."),
    "v_mosaic": (", Mosaike", ", mosaics"),
    "v_limited": ("Die Aussagekraft ist eingeschränkt: {reasons}.", "The result is limited: {reasons}."),
    "cancer_h": ("Krebsarten und BRCA1/BRCA2", "Cancer types and BRCA1/BRCA2"),
    "cancer_headers": (["Krebsart", "Zusammenhang mit BRCA1", "Zusammenhang mit BRCA2", "Risiko bei krankheitsverursachender Variante", "Dein Ergebnis"],
                       ["Cancer type", "Link with BRCA1", "Link with BRCA2", "Risk with a disease-causing variant", "Your result"]),
    "c_breast_f": ("Brustkrebs (Frauen)", "Breast cancer (women)"),
    "c_breast_m": ("Brustkrebs (Männer)", "Breast cancer (men)"),
    "c_ovarian": ("Eierstockkrebs (mit Eileiter und Bauchfell)", "Ovarian cancer (with fallopian tube and peritoneum)"),
    "c_prostate": ("Prostatakrebs", "Prostate cancer"),
    "c_pancreas": ("Bauchspeicheldrüsenkrebs", "Pancreatic cancer"),
    "c_testis": ("Hodenkrebs", "Testicular cancer"),
    "l_high": ("stark erhöht", "strongly increased"),
    "l_up": ("erhöht", "increased"),
    "l_some": ("leicht erhöht", "modestly increased"),
    "l_none": ("kein gesicherter Zusammenhang", "no established link"),
    "r_none": ("kein Befund in BRCA1/BRCA2, der dieses Risiko erhöht", "no finding in BRCA1/BRCA2 that raises this risk"),
    "r_unclear": ("nur Einträge mit unklarer Bedeutung in {genes}, siehe „Zu prüfen“", "only entries of uncertain significance in {genes}, see \u201cTo review\u201d"),
    "r_attention": ("Befund in {genes}, der dieses Risiko erhöhen kann; Bestätigung nötig, siehe „Zu prüfen“", "finding in {genes} that may raise this risk; needs confirmation, see \u201cTo review\u201d"),
    "r_none_limited": ("kein Befund in BRCA1/BRCA2, der dieses Risiko erhöht; Abdeckung von {genes} eingeschränkt, siehe Qualität",
                       "no finding in BRCA1/BRCA2 that raises this risk; coverage of {genes} limited, see Quality"),
    "r_not_linked": ("BRCA1/BRCA2 sagen hierzu nichts aus", "BRCA1/BRCA2 say nothing about this"),
    "k_breast_f": ("BRCA1 72 %, BRCA2 69 % bis 80 Jahre [1]", "BRCA1 72%, BRCA2 69% by age 80 [1]"),
    "k_breast_m": ("BRCA1 1,2 %, BRCA2 6,8 % bis 70 Jahre [2]", "BRCA1 1.2%, BRCA2 6.8% by age 70 [2]"),
    "k_ovarian": ("BRCA1 44 %, BRCA2 17 % bis 80 Jahre [1]", "BRCA1 44%, BRCA2 17% by age 80 [1]"),
    "k_prostate": ("BRCA2 27 % bis 75 Jahre, oft aggressiv; BRCA1 geringer [3]", "BRCA2 27% by age 75, often aggressive; BRCA1 lower [3]"),
    "k_pancreas": ("erhöht für beide Gene [4]", "increased for both genes [4]"),
    "k_testis": ("kein gesicherter Zusammenhang; am besten belegt ist CHEK2 (nicht geprüft) [5]", "no established link; best supported is CHEK2 (not checked) [5]"),
    "cancer_refs": ("[1] Kuchenbaecker 2017, JAMA, doi:10.1001/jama.2017.7112 · [2] Tai 2007, JNCI, doi:10.1093/jnci/djm203 · [3] Nyberg 2020, Eur Urol, doi:10.1016/j.eururo.2019.08.025 · [4] Li 2022, J Clin Oncol, doi:10.1200/JCO.21.02112 · [5] AlDubayan 2019, JAMA Oncol, doi:10.1001/jamaoncol.2018.6477. Die Risiken gelten für Trägerinnen und Träger einer krankheitsverursachenden Variante, nicht für dich, wenn keine gefunden wurde.",
                    "[1] Kuchenbaecker 2017, JAMA, doi:10.1001/jama.2017.7112 · [2] Tai 2007, JNCI, doi:10.1093/jnci/djm203 · [3] Nyberg 2020, Eur Urol, doi:10.1016/j.eururo.2019.08.025 · [4] Li 2022, J Clin Oncol, doi:10.1200/JCO.21.02112 · [5] AlDubayan 2019, JAMA Oncol, doi:10.1001/jamaoncol.2018.6477. The risks apply to carriers of a disease-causing variant, not to you if none was found."),
    "cancer_note": ("Zusammenhang nach NCCN Guidelines Genetic/Familial High-Risk Assessment: Breast, Ovarian, Pancreatic, and Prostate und ClinGen-Gen-Krankheits-Bewertungen; für Hodenkrebs ist kein Zusammenhang mit BRCA1/BRCA2 gesichert. Die Tabelle zeigt nur, ob in BRCA1 und BRCA2 etwas gefunden wurde. Andere Gene für diese Krebsarten (zum Beispiel PALB2, CHEK2, ATM, TP53, RAD51C, RAD51D, HOXB13, Mismatch-Reparatur-Gene) prüft brca12-wgs nicht. Kein Befund heißt nicht kein Risiko.",
                    "Link after the NCCN Guidelines Genetic/Familial High-Risk Assessment: Breast, Ovarian, Pancreatic, and Prostate and ClinGen gene-disease curations; for testicular cancer no link with BRCA1/BRCA2 is established. The table only shows whether something was found in BRCA1 and BRCA2. Other genes for these cancers (for example PALB2, CHEK2, ATM, TP53, RAD51C, RAD51D, HOXB13, mismatch-repair genes) are not checked by brca12-wgs. No finding does not mean no risk."),
    "how_h": ("Wie die Aussagen zustande kommen", "How the statements are made"),
    "how_sources": ("Datenquellen dieses Laufs: kleine Varianten mit {caller}; Annotation mit Ensembl VEP {vep} auf den MANE-Select-Transkripten (BRCA1 NM_007294.4, BRCA2 NM_000059.4); ClinVar-Ausgabe {clinvar}; Populationshäufigkeit gnomAD ({gnomad}); ACMG/AMP-Hinweise nach ClinGen ENIGMA BRCA1/BRCA2 VCEP v1.2; Exon-Kopienzahl gegen eine Kohorte von 1000-Genomes-Proben (30x).",
                    "Data sources of this run: small variants with {caller}; annotation with Ensembl VEP {vep} on the MANE Select transcripts (BRCA1 NM_007294.4, BRCA2 NM_000059.4); ClinVar release {clinvar}; population frequency gnomAD ({gnomad}); ACMG/AMP hints after ClinGen ENIGMA BRCA1/BRCA2 VCEP v1.2; exon copy number against a cohort of 1000 Genomes samples (30x)."),
    "how_attention": ("Als Befund, der bestätigt werden muss, zählt: eine Variante, die ClinVar ohne widersprüchliche Einträge als pathogen oder wahrscheinlich pathogen führt; ein Kandidat für Funktionsverlust (vorzeitiger Stopp, Leserasterverschiebung, Spleißstelle ±1/2 oder Startcodon, PVS1 nach ENIGMA Tabelle 4); ein Kopienzahl-Kandidat aus CNVkit; eine Exon-Kopienzahl-Abweichung gegenüber der Kohorte; eine strukturelle Variante mit FILTER PASS.",
                      "Counted as a finding that needs confirmation: a variant that ClinVar lists as pathogenic or likely pathogenic without conflicting entries; a loss-of-function candidate (premature stop, frameshift, splice site ±1/2 or start codon, PVS1 after ENIGMA Table 4); a copy-number candidate from CNVkit; an exon copy-number deviation against the cohort; a structural variant with FILTER PASS."),
    "how_unclear": ("Als unklar zählt: ClinVar unklare Bedeutung, Mosaik-Kandidaten, strukturelle Varianten ohne PASS und Einschub-Signale. Alle übrigen Varianten (gutartig, wahrscheinlich gutartig, häufig in der Bevölkerung, ohne ClinVar-Eintrag und ohne Funktionsverlust) zählen nicht als Befund; sie stehen in den Details je Gen.",
                    "Counted as uncertain: ClinVar uncertain significance, mosaic candidates, structural variants without PASS and insertion signals. All other variants (benign, likely benign, common in the population, without a ClinVar entry and without loss of function) do not count as a finding; they are listed in the details per gene."),
    "how_rule": ("Gesamtaussage: mindestens ein Befund ergibt „bestätigen“, sonst ergeben unklare Einträge „unklar“, sonst „kein Hinweis“. Die Krebsarten-Tabelle wendet dieselbe Regel je Gen an und überträgt sie auf die Krebsarten, die mit diesem Gen verbunden sind. Die Qualität ist eingeschränkt, wenn ein Exon unsicher abgedeckt ist, die mittlere Tiefe unter {min:g}x liegt oder die Read-Prüfung einen Hinweis meldet.",
                 "Overall statement: at least one finding gives \u201cneeds confirmation\u201d, otherwise uncertain entries give \u201cuncertain\u201d, otherwise \u201cno sign\u201d. The cancer types table applies the same rule per gene and carries it over to the cancers linked to that gene. Quality is limited if an exon is covered uncertainly, the mean depth is below {min:g}x or the read check reports a note."),
    "how_trace": ("Nachvollziehbar: jeder Wert steht in brca12_results.json, jeder ausgeführte Befehl in brca12_commands.sh, alle Programmversionen unter „Details je Gen“ → „Lauf“. Methode und Prüfung an GIAB- und 1000-Genomes-Proben: Forschungsbeitrag {url}.",
                  "Traceable: every value is in brca12_results.json, every executed command in brca12_commands.sh, all program versions under \u201cDetails per gene\u201d \u2192 \u201cRun\u201d. Method and validation on GIAB and 1000 Genomes samples: research article {url}."),
    "card_small": ("SNVs und Indels gegenüber der Referenz", "SNVs and indels relative to the reference"),
    "card_small_tip": ("Alle kleinen Varianten in BRCA1 und BRCA2, auch häufige und harmlose. Einzeln in den Details je Gen.",
                       "All small variants in BRCA1 and BRCA2, including common and benign ones. Listed in the details per gene."),
    "card_plp": ("ClinVar pathogen oder wahrscheinlich pathogen", "ClinVar pathogenic or likely pathogenic"),
    "card_plp_tip": ("Varianten mit ClinVar-Angabe pathogen oder wahrscheinlich pathogen, ohne widersprüchliche Einträge. Dazu {lof} Kandidaten für Funktionsverlust (PVS1-Hinweis).",
                     "Variants with a ClinVar entry pathogenic or likely pathogenic, without conflicting entries. In addition {lof} loss-of-function candidates (PVS1 hint)."),
    "card_struct": ("Kandidaten für CNV und strukturelle Varianten", "Candidates for CNVs and structural variants"),
    "card_struct_tip": ("CNVkit-Kandidaten {cnv}, Exon-Kopienzahl (Kohorte) {exon}, SV-Rufer {sv}, Einschub-Signale {ins}. Schwächere Hinweise stehen in den Details je Gen.",
                        "CNVkit candidates {cnv}, exon copy number (cohort) {exon}, SV caller {sv}, insertion signals {ins}. Weaker hints are in the details per gene."),
    "card_exons": ("unsichere Exons", "uncertain exons"),
    "card_exons_tip": ("Exons mit zu geringer Tiefe, zu wenig Basen ab 20x oder vielen mehrdeutig kartierten Reads, von {total} Exons.",
                       "Exons with too little depth, too few bases at 20x or many ambiguously mapped reads, out of {total} exons."),
    "low_depth": ("Geringe Tiefe: {gene} hat im Mittel {depth:g}x über den Exons (Schwelle {min:g}x). Kleine Varianten, Kopienzahl und strukturelle Varianten sind dort unsicherer; ein leeres Ergebnis sagt wenig.",
                  "Low depth: {gene} has a mean of {depth:g}x over the exons (threshold {min:g}x). Small variants, copy number and structural variants are less certain there; an empty result says little."),
    "notable": ("Zu prüfen", "To review"),
    "notable_headers": (["Gen", "Art", "Angabe", "Genotyp", "ClinVar", "Hinweise"],
                        ["Gene", "Kind", "Finding", "Genotype", "ClinVar", "Notes"]),
    "notable_none": ("Keine Einträge in den geprüften Kategorien.", "No entries in the checked categories."),
    "notable_note": ("Maschinell zusammengestellt, keine Klassifikation. ACMG/AMP-Angaben sind Hinweise auf einzelne Kriterien. Jede Angabe braucht eine Bestätigung in einem akkreditierten Labor.",
                     "Compiled automatically, not a classification. ACMG/AMP entries are hints on single criteria. Every entry needs confirmation in an accredited laboratory."),
    "kind_plp": ("ClinVar pathogen / wahrscheinlich pathogen", "ClinVar pathogenic / likely pathogenic"),
    "kind_lof": ("Kandidat Funktionsverlust", "Loss-of-function candidate"),
    "kind_vus": ("ClinVar unklare Bedeutung", "ClinVar uncertain significance"),
    "kind_mosaic": ("Mosaik-Kandidat", "Mosaic candidate"),
    "kind_cnv": ("CNV-Kandidat (CNVkit)", "CNV candidate (CNVkit)"),
    "kind_exon": ("Exon-Kopienzahl (Kohorte)", "Exon copy number (cohort)"),
    "kind_sv": ("Strukturelle Variante", "Structural variant"),
    "kind_ins": ("Einschub-Signal", "Insertion signal"),
    "coverage": ("Abdeckung", "Coverage"),
    "cov_headers": (["Gen", "Exons", "mittlere Tiefe", "geringste Exon-Tiefe", "unsichere Exons"],
                    ["Gene", "Exons", "Mean depth", "Lowest exon depth", "Uncertain exons"]),
    "method": ("Methode und Grenzen", "Method and limitations"),
    "method_text": ("Die Reads über BRCA1 und BRCA2 werden mit {caller} auf kleine Varianten untersucht, mit Ensembl VEP und ClinVar annotiert und nach den ENIGMA-Schwellen (ClinGen ENIGMA BRCA1/BRCA2 VCEP v1.2) mit Hinweisen auf einzelne ACMG/AMP-Kriterien versehen. Kopienzahl kommt aus CNVkit und einem exonweisen Vergleich mit einer Kohorte, strukturelle Varianten aus dem SV-Rufer und aus Einschub-Signalen. Die Abdeckung wird je Exon geprüft. Die Prüfung an GIAB- und 1000-Genomes-Proben steht im Forschungsbeitrag.",
                    "The reads over BRCA1 and BRCA2 are searched for small variants with {caller}, annotated with Ensembl VEP and ClinVar and given hints on single ACMG/AMP criteria after the ENIGMA thresholds (ClinGen ENIGMA BRCA1/BRCA2 VCEP v1.2). Copy number comes from CNVkit and a per-exon comparison with a cohort, structural variants from the SV caller and from insertion signals. Coverage is checked per exon. The validation on GIAB and 1000 Genomes samples is in the research article."),
    "details_note": ("Ein Klick öffnet ein Gen. Alle Werte aus brca12_results.json. Gedruckt wird der Reiter, der gerade offen ist; dieser hier mit allen Genen aufgeklappt. Orange hinterlegt: unsichere Exons und Varianten in unsicheren Exons.",
                     "A click opens a gene. All values from brca12_results.json. Printing prints the open tab; this one with all genes expanded. Orange background: uncertain exons and variants in uncertain exons."),
    "gene_kv": ("{chrom}, {strand}, Transkript {transcript} · mittlere Exon-Tiefe {depth} · unsichere Exons {uncertain} · kleine Varianten {small}",
                "{chrom}, {strand}, transcript {transcript} · mean exon depth {depth} · uncertain exons {uncertain} · small variants {small}"),
    "gene_sum": ("{small} Varianten · {notable} zu prüfen · {uncertain} unsichere Exons", "{small} variants · {notable} to review · {uncertain} uncertain exons"),
    "run": ("Lauf", "Run"),
    "run_sum": ("Read-Prüfung, Werkzeuge, Befehle", "Read check, tools, commands"),
    "cite": ("Zitieren", "Cite"),
    "research": ("Forschungsbeitrag mit Methode und Prüfung", "Research article with method and validation"),
    "source": ("Quelltext", "Source code"),
    "files": ("Alle Einzelwerte stehen in brca12_results.json, brca12_report.md und brca12_report.xlsx neben dieser Seite; jeder ausgeführte Befehl in brca12_commands.sh.",
              "All single values are in brca12_results.json, brca12_report.md and brca12_report.xlsx next to this page; every executed command in brca12_commands.sh."),
}

CSS = """:root{--bg:#FAF9F5;--surface:#FFFFFF;--fg:#141413;--mut:#3F3E38;--soft-fg:#6E6C63;--line:#D8D5C7;--hair:rgba(216,213,199,.78);--acc:#D97757;--acc-t:#B34D30;--ok:#248A3D;--okdot:#34C759;--warn:#B34D30;--soft:#F2F0E8;--f1:#C2410C;--f2:#9A7400;--f3:#4D7C0F;--f4:#0E7490;--shadow:0 1px 2px rgba(20,20,19,.04),0 10px 28px rgba(20,20,19,.055);--r:14px;--ease:cubic-bezier(.22,.61,.36,1)}
@media (prefers-color-scheme:dark){:root{--bg:#1C1C1E;--surface:#2C2C2E;--fg:#FFFFFF;--mut:rgba(235,235,245,.72);--soft-fg:rgba(235,235,245,.55);--line:#38383A;--hair:rgba(255,255,255,.1);--acc-t:#E38B69;--ok:#30D158;--okdot:#30D158;--warn:#E38B69;--soft:#242426;--f1:#FB8A5C;--f2:#E3B341;--f3:#9BCB5A;--f4:#4FC3D9;--shadow:0 1px 2px rgba(0,0,0,.28),0 12px 30px rgba(0,0,0,.3)}}
*{box-sizing:border-box}
html{-webkit-text-size-adjust:100%}
body{margin:0;background:var(--bg);color:var(--fg);font:400 14px/1.6 'IBM Plex Mono',ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;-webkit-font-smoothing:antialiased}
.bar{position:sticky;top:0;z-index:5;background:color-mix(in srgb,var(--bg) 82%,transparent);-webkit-backdrop-filter:saturate(180%) blur(18px);backdrop-filter:saturate(180%) blur(18px);border-bottom:1px solid var(--hair)}
.bar-in{max-width:1000px;margin:0 auto;padding:12px 24px;display:flex;align-items:center;justify-content:space-between;gap:16px}
.brand{font-weight:500;letter-spacing:.02em;font-size:13px}.brand i{font-style:normal;color:var(--acc-t)}
.seg{display:inline-flex;background:var(--soft);border:1px solid var(--hair);border-radius:999px;padding:3px;gap:2px}
.seg button{font:inherit;font-size:12px;font-weight:500;background:none;border:0;border-radius:999px;color:var(--soft-fg);padding:5px 14px;cursor:pointer;transition:background 140ms var(--ease),color 140ms var(--ease)}
.seg button:hover{color:var(--fg)}
body[data-lang=de] nav.lang button[data-set=de],body[data-lang=en] nav.lang button[data-set=en],body[data-tab=r] nav.tabs button[data-tab=r],body[data-tab=d] nav.tabs button[data-tab=d]{background:var(--surface);color:var(--fg);box-shadow:0 1px 2px rgba(20,20,19,.12)}
main{max-width:1000px;margin:0 auto;padding:40px 24px 64px}
body[data-lang=de] main.l-en,body[data-lang=en] main.l-de{display:none}
.eyebrow{font-size:12px;font-weight:500;letter-spacing:.1em;text-transform:uppercase;color:var(--acc-t);margin-bottom:8px}
h1{font-size:30px;line-height:1.2;font-weight:500;margin:0;letter-spacing:-.01em}
h2{font-size:12px;font-weight:500;letter-spacing:.1em;text-transform:uppercase;color:var(--soft-fg);margin:40px 0 12px}
.sub{color:var(--soft-fg);margin-top:8px;font-size:13px}
nav.tabs{margin:28px 0 4px}
.cards{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin-top:28px}
.card{background:var(--surface);border:1px solid var(--hair);border-radius:var(--r);padding:18px 18px 16px;box-shadow:var(--shadow)}
.card b{display:block;font-size:26px;line-height:1.15;font-weight:500;letter-spacing:-.01em}
.card span{display:block;margin-top:6px;color:var(--soft-fg);font-size:12px;line-height:1.4}
.panel{background:var(--surface);border:1px solid var(--hair);border-radius:var(--r);box-shadow:var(--shadow);padding:6px 20px;overflow-x:auto}
table{width:100%;border-collapse:collapse}
th{font-weight:500;text-align:left;color:var(--soft-fg);font-size:11px;letter-spacing:.06em;text-transform:uppercase;border-bottom:1px solid var(--line);padding:12px 12px 10px 0}
td{border-bottom:1px solid var(--hair);padding:11px 12px 11px 0;vertical-align:top}
tr:last-child td{border-bottom:0}
td.g{font-weight:500;white-space:nowrap}
.call{font-weight:500;white-space:nowrap}
.full{display:block;color:var(--soft-fg);font-size:12px}
.d{white-space:nowrap;line-height:1.6}
table.tf th:nth-child(4),table.tf td:nth-child(4){width:96px}
.f{display:inline-block;width:8px;height:8px;border-radius:50%;border:1px solid var(--line);margin-right:4px}
.f.on{background:var(--okdot);border-color:var(--okdot)}
.note{color:var(--soft-fg);font-size:12px;line-height:1.55}
.warn{color:var(--warn)}
.pos{color:var(--warn);font-weight:500}
.cls td{background:none;color:var(--acc-t);font-size:11px;font-weight:500;letter-spacing:.1em;text-transform:uppercase;padding:20px 0 8px;border-bottom:1px solid var(--line)}
a{color:var(--acc-t);text-decoration:none}a:hover{text-decoration:underline}
details{background:var(--surface);border:1px solid var(--hair);border-radius:var(--r);margin:0 0 10px;padding:0 20px;box-shadow:0 1px 2px rgba(20,20,19,.04);transition:box-shadow 220ms var(--ease)}
details[open]{box-shadow:var(--shadow);padding-bottom:12px}
summary{cursor:pointer;list-style:none;padding:14px 0;display:flex;align-items:baseline;gap:12px}summary::-webkit-details-marker{display:none}
summary:after{content:'';margin-left:auto;width:7px;height:7px;border-right:1.5px solid var(--soft-fg);border-bottom:1.5px solid var(--soft-fg);transform:rotate(-45deg);transition:transform 220ms var(--ease);align-self:center}
details[open] summary:after{transform:rotate(45deg)}
summary b{font-weight:500;min-width:4.5em}
.kv{color:var(--mut);font-size:12px;margin:10px 0}
.kv b{color:var(--fg);font-weight:500}
.dt{font-size:12px;margin:6px 0 12px;background:var(--bg);border-radius:10px}
.dt th{padding:9px 12px 7px}.dt td{padding:7px 12px}.dt td:last-child{min-width:16em}
.reason{font-size:12px;color:var(--mut);border-left:2px solid var(--acc);padding:2px 0 2px 10px;margin:4px 0 14px}
.ok{color:var(--ok)}
.nm{white-space:nowrap;font-weight:500;letter-spacing:.01em}.nm .gn{color:var(--soft-fg);font-weight:400}
.fd{-webkit-print-color-adjust:exact;print-color-adjust:exact}
.fd.k1{color:var(--f1)}.fd.k2{color:var(--f2)}.fd.k3{color:var(--f3)}.fd.k4{color:var(--f4)}
.fd.off{color:var(--soft-fg);font-weight:400}
.c{color:var(--soft-fg)}
.legend{background:var(--surface);border:1px solid var(--hair);border-radius:var(--r);box-shadow:var(--shadow);padding:16px 20px 6px;margin:0 0 14px;overflow-x:auto}.legend .note{margin:2px 0 10px}h4.ah{font-size:12px;font-weight:500;letter-spacing:.1em;text-transform:uppercase;color:var(--acc-t);margin:18px 0 2px;padding-top:12px;border-top:1px solid var(--hair)}details h4.ah:first-of-type{border-top:0;padding-top:0}.legend h3{font-size:12px;font-weight:500;letter-spacing:.06em;text-transform:uppercase;color:var(--soft-fg);margin:0 0 6px}
.alert{margin:22px 0 0;padding:12px 16px;border-radius:var(--r);background:color-mix(in srgb,var(--acc) 12%,var(--surface));border:1px solid color-mix(in srgb,var(--acc) 40%,transparent);font-size:12px;line-height:1.55}.alert.soft{background:var(--surface);border-color:var(--hair)}
@media print{.bar,nav.tabs{display:none}summary:after{display:none}}
body[data-tab=r] section.tab-d,body[data-tab=d] section.tab-r{display:none}
.ruo{font-size:11px;font-weight:500;letter-spacing:.06em;text-transform:uppercase;color:var(--acc-t);border:1px solid color-mix(in srgb,var(--acc) 40%,transparent);border-radius:999px;display:inline-block;padding:3px 12px;margin:0 0 22px}
.head{display:grid;grid-template-columns:repeat(4,1fr);gap:0;margin-top:22px;background:var(--surface);border:1px solid var(--hair);border-radius:var(--r);box-shadow:var(--shadow);overflow:hidden}
.head div{padding:12px 16px;border-right:1px solid var(--hair);border-bottom:1px solid var(--hair);min-width:0}.head div:nth-child(4n){border-right:0}.head div:nth-last-child(-n+4){border-bottom:0}
.head span{display:block;font-size:11px;letter-spacing:.06em;text-transform:uppercase;color:var(--soft-fg)}.head b{display:block;font-weight:500;font-size:13px;margin-top:3px;overflow-wrap:anywhere}
table.clin td{padding-top:12px;padding-bottom:12px}table.clin td:nth-child(2),table.clin td:nth-child(3){width:40%}table.clin .call{font-size:15px}
.rw{font-size:11px;color:var(--soft-fg);margin-left:2px}
.gl{display:flex;align-items:center;gap:14px;padding:14px 20px}.gl code{flex:1;font:400 12px/1.6 'IBM Plex Mono',ui-monospace,monospace;overflow-wrap:anywhere;color:var(--mut)}
.copy{font:inherit;font-size:12px;font-weight:500;background:var(--soft);border:1px solid var(--hair);border-radius:999px;padding:5px 14px;color:var(--fg);cursor:pointer}
details.more{margin-top:14px}details.more table{margin-bottom:6px}
.prose{padding:8px 20px}.prose p{margin:10px 0;color:var(--mut);font-size:13px;line-height:1.65}
@page{size:A4;margin:16mm 14mm 18mm;@bottom-right{content:counter(page) ' / ' counter(pages);font:400 9px 'IBM Plex Mono',monospace;color:#6E6C63}@bottom-left{content:'brca12-wgs · RUO';font:400 9px 'IBM Plex Mono',monospace;color:#6E6C63}}
@media print{.copy{display:none}.head{grid-template-columns:repeat(4,1fr)}}
@media (max-width:720px){table.clin .call{white-space:normal;font-size:13px}table.clin td:nth-child(2),table.clin td:nth-child(3){width:auto}#tt{max-width:calc(100vw - 16px)}.head{grid-template-columns:repeat(2,1fr)}.head div:nth-child(4n){border-right:1px solid var(--hair)}.head div:nth-child(2n){border-right:0}.head div:nth-last-child(-n+4){border-bottom:1px solid var(--hair)}.head div:nth-last-child(-n+2){border-bottom:0}}
.tip{cursor:help;text-decoration:underline dotted color-mix(in srgb,var(--soft-fg) 70%,transparent);text-underline-offset:3px;text-decoration-thickness:1px;outline:none}.tip:focus-visible{outline:2px solid var(--acc);outline-offset:2px;border-radius:3px}
#tt{position:fixed;z-index:50;max-width:min(420px,calc(100vw - 24px));background:var(--fg);color:var(--bg);font:400 11.5px/1.55 'IBM Plex Mono',ui-monospace,monospace;letter-spacing:0;text-transform:none;white-space:pre-line;padding:10px 12px;border-radius:10px;box-shadow:0 8px 28px rgba(20,20,19,.22);pointer-events:none;opacity:0;transform:translateY(2px);transition:opacity 120ms var(--ease),transform 120ms var(--ease)}#tt.on{opacity:1;transform:none}
@media print{#tt{display:none}.tip{text-decoration:none}}
footer{margin-top:48px;border-top:1px solid var(--line);padding-top:18px;color:var(--soft-fg);font-size:12px}
footer p{margin:6px 0}
.cite{color:var(--fg)}
@media (max-width:720px){.cards{grid-template-columns:repeat(2,1fr)}.hide{display:none}h1{font-size:24px}.panel,details{padding-left:14px;padding-right:14px}}
@media print{:root{--bg:#fff;--surface:#fff;--shadow:none}main{padding:0;max-width:none}a{color:inherit}tr,summary,.kv{break-inside:avoid}h4{break-after:avoid}body{font-size:11px}.card,.panel,details{box-shadow:none}}
.report-figure{margin:6px 0 14px;padding:12px 14px;background:#fff;color:#1d1d1f;border:1px solid var(--hair);border-radius:10px;break-inside:avoid;-webkit-print-color-adjust:exact;print-color-adjust:exact}
.report-figure svg{max-width:100%;height:auto;display:block}.figure-note{font-size:11px;color:#636366;margin:4px 0 0}
tr.uncertain td{background:color-mix(in srgb,var(--acc) 9%,transparent);-webkit-print-color-adjust:exact;print-color-adjust:exact}
.verdict{margin:28px 0 0;padding:16px 20px;border-radius:var(--r);background:var(--surface);border:1px solid var(--hair);box-shadow:var(--shadow);border-left:4px solid var(--okdot)}
.verdict.v-unclear{border-left-color:var(--f2)}.verdict.v-attention{border-left-color:var(--acc)}
.verdict b{display:block;font-weight:500;font-size:15px;line-height:1.45}.verdict p{margin:8px 0 0;color:var(--mut);font-size:12px;line-height:1.6}
.dt{width:100%}.dt td{overflow-wrap:anywhere}
pre{background:var(--bg);padding:12px 14px;border-radius:10px;white-space:pre-wrap;overflow-wrap:anywhere;font:400 12px/1.6 'IBM Plex Mono',ui-monospace,monospace;margin:8px 0}
.head b em{font-style:normal}
.dt td:nth-child(n){min-width:3.5em}
details h3{font-size:12px;font-weight:500;letter-spacing:.06em;text-transform:uppercase;color:var(--soft-fg);margin:18px 0 6px}"""

SCRIPT = """
function aria(){document.querySelectorAll('nav.lang button').forEach(b=>b.setAttribute('aria-pressed',b.dataset.set===document.body.dataset.lang));document.querySelectorAll('nav.tabs button').forEach(b=>b.setAttribute('aria-pressed',b.dataset.tab===document.body.dataset.tab))}
function setL(l){document.body.dataset.lang=l;document.documentElement.lang=l;aria();try{localStorage.setItem('brca12-wgs-lang',l)}catch(e){}}
try{const s=localStorage.getItem('brca12-wgs-lang');if(s==='de'||s==='en')setL(s)}catch(e){}
document.querySelectorAll('nav.lang button').forEach(b=>b.addEventListener('click',()=>setL(b.dataset.set)));
document.querySelectorAll('nav.tabs button').forEach(b=>b.addEventListener('click',()=>{document.body.dataset.tab=b.dataset.tab;aria();scrollTo(0,0)}));
const tt=document.createElement('div');tt.id='tt';tt.setAttribute('role','tooltip');document.body.appendChild(tt);let cur=null;
function showT(el){cur=el;tt.textContent=el.dataset.tip;tt.classList.add('on');const r=el.getBoundingClientRect(),w=tt.offsetWidth,h=tt.offsetHeight;let x=Math.min(Math.max(8,r.left),innerWidth-w-8),y=r.bottom+8;if(y+h>innerHeight-8)y=Math.max(8,r.top-h-8);tt.style.left=x+'px';tt.style.top=y+'px'}
function hideT(){cur=null;tt.classList.remove('on')}
document.addEventListener('pointerover',e=>{if(e.pointerType!=='mouse')return;const el=e.target.closest('[data-tip]');if(el&&el!==cur)showT(el);else if(!el&&cur)hideT()});
document.addEventListener('focusin',e=>{const el=e.target.closest('[data-tip]');if(el)showT(el)});document.addEventListener('focusout',hideT);
document.addEventListener('click',e=>{const el=e.target.closest('[data-tip]');if(el&&el.tagName!=='A')showT(el);else if(!el)hideT()});
addEventListener('scroll',hideT,{passive:true});addEventListener('keydown',e=>{if(e.key==='Escape')hideT()});
document.querySelectorAll('button.copy').forEach(b=>b.addEventListener('click',()=>{const c=document.getElementById(b.dataset.copy);if(c&&navigator.clipboard)navigator.clipboard.writeText(c.textContent).then(()=>{const o=b.textContent;b.textContent='\\u2713';setTimeout(()=>b.textContent=o,1200)})}));
aria();
addEventListener('beforeprint',()=>document.querySelectorAll('details').forEach(d=>d.open=true));
"""


def _t(key, lang, **kwargs):
    value = TEXT[key][0 if lang == "de" else 1]
    return value.format(**kwargs) if kwargs else value


def _fonts():
    """@font-face rules with the embedded font; empty when the font files are missing (system font is used)."""
    rules = []
    for weight in (400, 500):
        try:
            with open(os.path.join(FONT_DIR, f"IBMPlexMono-Latin-{weight}.woff2"), "rb") as handle:
                data = base64.b64encode(handle.read()).decode("ascii")
        except OSError:
            return ""
        rules.append(f"@font-face{{font-family:'IBM Plex Mono';font-weight:{weight};font-style:normal;font-display:swap;"
                     f"src:url(data:font/woff2;base64,{data}) format('woff2')}}")
    return "\n".join(rules)


def _css_text(text):
    """Text that is safe inside a CSS string in a <style> block."""
    return "".join(c if c.isalnum() or c in " ._-" else "_" for c in str(text))


def _tip(label, tip, css="tip"):
    return f'<span class="{css}" tabindex="0" data-tip="{html.escape(str(tip), quote=True)}">{html.escape(str(label))}</span>'


def report_number(result):
    """FNV-1a (32 bit) over all results and the input file names; independent of time and paths."""
    keys = ("genes", "exonCoverage", "lowCoverageIntervals", "cnvs", "cnvSegmentsNotSupported", "exonDosageHints",
            "exonCopyNumber", "svs", "insertionSignals", "mosaicCandidates", "readQC")
    payload = json.dumps({k: result.get(k) for k in keys} | {"files": result["input"]["files"]},
                         sort_keys=True, ensure_ascii=False, default=str).encode("utf-8")
    h = 0x811C9DC5
    for byte in payload:
        h = ((h ^ byte) * 0x01000193) & 0xFFFFFFFF
    return f"BRCA-{h:08X}"


def _gene_depth(result):
    """{gene: (exon count, mean of exon mean depths, lowest exon mean depth, uncertain exons)}"""
    out = {}
    for gene in GENES:
        exons = [e for e in result["exonCoverage"] if e["gene"] == gene and e.get("meanDepth") is not None]
        if exons:
            depths = [e["meanDepth"] for e in exons]
            out[gene] = (len(exons), round(sum(depths) / len(depths), 1), min(depths),
                         sum(1 for e in exons if e["status"] != "ok"))
    return out


def quality(result):
    """(passed, [(reason key, arguments)]) for the report header; language-neutral."""
    minimum = result["thresholds"]["minExonDepth"]
    reasons = []
    uncertain = sum(1 for e in result["exonCoverage"] if e["status"] != "ok")
    if uncertain:
        reasons.append(("q_exons", {"n": uncertain}))
    for gene, (_, mean, _, _) in _gene_depth(result).items():
        if mean < minimum:
            reasons.append(("q_depth", {"gene": gene, "depth": mean}))
    flags = len((result.get("readQC") or {}).get("flags", []))
    if flags:
        reasons.append(("q_reads", {"n": flags}))
    return not reasons, reasons


def _is_plp(v):
    sig = v["clinvar"]["significance"] or ""
    return "athogenic" in sig and "onflicting" not in sig


def _is_vus(v):
    sig = v["clinvar"]["significance"] or ""
    return "ncertain" in sig and "onflicting" not in sig


def notable(result):
    """Entries to review as (tier, gene, kind key, item, source): tier 'attention' for findings that may
    indicate a predisposition (ClinVar P/LP, loss of function, CNV or exon copy number calls, SV with
    FILTER PASS), 'unclear' for the rest the pipeline reports as candidates."""
    out = []
    for gene in GENES:
        for v in result["genes"][gene]["smallVariants"]:
            if _is_plp(v):
                out.append(("attention", gene, "kind_plp", v, "small"))
            elif any(h["code"].startswith("PVS1") for h in v.get("acmgHints", [])):
                out.append(("attention", gene, "kind_lof", v, "small"))
            elif _is_vus(v):
                out.append(("unclear", gene, "kind_vus", v, "small"))
    for v in result.get("mosaicCandidates", []):
        out.append(("unclear", v["gene"], "kind_mosaic", v, "small"))
    for c in result["cnvs"]:
        out.append(("attention", c["gene"], "kind_cnv", c, "cnv"))
    for c in (result.get("exonCopyNumber") or {}).get("calls", []):
        out.append(("attention", c["gene"], "kind_exon", c, "cnv"))
    for s in result["svs"]:
        out.append(("attention" if s.get("filter") == "PASS" else "unclear", s["gene"], "kind_sv", s, "sv"))
    for i in result.get("insertionSignals", []):
        out.append(("unclear", i["gene"], "kind_ins", i, "ins"))
    return out


def _row(entry, lang):
    _, gene, kind, item, source = entry
    hint_key = "textDE" if lang == "de" else "textEN"

    def flags(x):
        return "; ".join(str(tr(f, lang)) for f in x.get("flags", []) + x.get("qualityFlags", []))

    def exons(x):
        return ", ".join(str(tr(e, lang)) for e in x.get("affectedExons", []))

    if source == "small":
        name = " ".join(x for x in (item.get("hgvsc"), item.get("hgvsp")) if x) or \
            f"{item['chrom']}:{item['pos']} {item['ref']}>{item['alt']}"
        notes = "; ".join([f"{h['code']}: {h[hint_key]}" for h in item.get("acmgHints", [])]
                          + ([flags(item)] if flags(item) else []))
        genotype = f"AF {item.get('alleleFraction')}" if kind == "kind_mosaic" else tr(item.get("genotype"), lang)
        return [gene, _t(kind, lang), name, genotype, item["clinvar"]["significance"], notes]
    if source == "cnv":
        return [gene, _t(kind, lang),
                f"{tr(item['type'], lang)} {item['chrom']}:{item['start']}-{item['end']} · CN {item['copyNumber']} · {exons(item)}",
                None, None, flags(item)]
    if source == "sv":
        pos = f"{item['chrom']}:{item['pos']}" + (f"-{item['end']}" if item.get("end") else "")
        return [gene, _t(kind, lang), f"{tr(item['type'], lang)} {pos} · {item.get('filter')}",
                tr(item.get("genotype"), lang), None, flags(item)]
    return [gene, _t(kind, lang), f"{tr(item['type'], lang)} {item['chrom']}:{item['pos']}-{item['end']}", None, None, flags(item)]


def notable_rows(result, lang):
    na = ui("na", lang)
    return [[na if c is None or c == "" else str(c) for c in _row(e, lang)] for e in notable(result)]


def verdict(result):
    """('none' | 'unclear' | 'attention', count) for the overall statement."""
    entries = notable(result)
    attention = [e for e in entries if e[0] == "attention"]
    if attention:
        return "attention", len(attention)
    if entries:
        return "unclear", len(entries)
    return "none", 0


# Cancer types with the strength of their link to BRCA1 and BRCA2 (NCCN Genetic/Familial High-Risk
# Assessment: Breast, Ovarian, Pancreatic, and Prostate; ClinGen gene-disease validity). Levels only,
# no percentages: the page states whether a finding exists, not a personal risk.
CANCERS = (
    ("c_breast_f", "l_high", "l_high"),
    ("c_breast_m", "l_up", "l_high"),
    ("c_ovarian", "l_high", "l_up"),
    ("c_prostate", "l_some", "l_high"),
    ("c_pancreas", "l_up", "l_up"),
    ("c_testis", "l_none", "l_none"),
)


def gene_levels(result):
    """Per gene 'none' | 'unclear' | 'attention', as the overall statement but for one gene."""
    levels = {gene: "none" for gene in GENES}
    for tier, gene, *_ in notable(result):
        if gene in levels and levels[gene] != "attention":
            levels[gene] = "attention" if tier == "attention" else "unclear"
    return levels


def cancer_rows(result, lang):
    levels = gene_levels(result)
    minimum = result["thresholds"]["minExonDepth"]
    depth = _gene_depth(result)
    limited = {g for g, (_, mean, _, uncertain) in depth.items() if uncertain or mean < minimum}
    rows = []
    for cancer, *links in CANCERS:
        linked = [g for g, link in zip(GENES, links) if link != "l_none"]
        cells = [_t(link, lang) for link in links]
        if not linked:
            outcome = _t("r_not_linked", lang)
        else:
            for level in ("attention", "unclear"):
                hit = [g for g in linked if levels[g] == level]
                if hit:
                    outcome = _t(f"r_{level}", lang, genes=" & ".join(hit))
                    break
            else:
                weak = [g for g in linked if g in limited]
                outcome = _t("r_none_limited", lang, genes=" & ".join(weak)) if weak else _t("r_none", lang)
        rows.append([_t(cancer, lang)] + cells + [_t("k" + cancer[1:], lang), outcome])
    return rows


def _table_html(headers, rows, mark=None, css=""):
    esc = html.escape
    out = [f"<table{f' class={chr(39)}{css}{chr(39)}' if css else ''}><thead><tr>"
           + "".join(f"<th>{esc(str(h))}</th>" for h in headers) + "</tr></thead><tbody>"]
    for row in rows:
        cls = " class='uncertain'" if mark and any(mark in str(c) for c in row) else ""
        out.append(f"<tr{cls}>" + "".join(f"<td>{esc(str(c))}</td>" for c in row) + "</tr>")
    out.append("</tbody></table>")
    return "".join(out)


def _gene_result(result, gene):
    """A copy of the result restricted to one gene, for the per-gene sections."""
    def keep(items):
        return [x for x in items if x.get("gene") == gene]
    g = dict(result)
    g["genes"] = {k: (v if k == gene else dict(v, smallVariants=[])) for k, v in result["genes"].items()}
    for key in ("exonCoverage", "lowCoverageIntervals", "cnvs", "exonDosageHints", "svs", "insertionSignals",
                "mosaicCandidates"):
        g[key] = keep(result.get(key) or [])
    ecn = result.get("exonCopyNumber") or {}
    g["exonCopyNumber"] = dict(ecn, calls=keep(ecn.get("calls", [])), hints=keep(ecn.get("hints", [])))
    g["readQC"] = None
    return g


def _main(result, lang, sections_for, disclaimer, limitations):
    esc = html.escape
    t = lambda key, **kw: _t(key, lang, **kw)  # noqa: E731
    pipeline = result.get("pipeline") or {}
    version = pipeline.get("version") or ui("unknown", lang)
    doi = pipeline.get("doi") or DOI
    minimum = result["thresholds"]["minExonDepth"]
    depth = _gene_depth(result)
    summary = result.get("summary") or {}
    passed, reasons = quality(result)
    reason_text = "; ".join(t(key, **args) for key, args in reasons)
    date = (result.get("generatedAt") or "")[:10]
    if lang == "de" and len(date) == 10:
        date = f"{date[8:10]}.{date[5:7]}.{date[0:4]}"
    sources = result.get("evidenceSources") or {}
    versions = result.get("toolVersions") or {}
    cells = [
        (_tip(t("sample"), t("sample_tip")), esc(result["sample"])),
        (_tip(t("report_no"), t("report_no_tip")), esc(report_number(result))),
        (_tip(t("input"), t("input_tip", mode=result["input"]["mode"], ref=os.path.basename(result["reference"]["path"]),
                            assembly=result["reference"].get("assembly"))), esc(", ".join(result["input"]["files"]))),
        (_tip(t("analysis"), t("analysis_tip", ts=result.get("generatedAt"))), esc(date)),
        (_tip(t("program"), t("program_tip", version=version, doi=doi, caller=result["snvCaller"],
                              session=result.get("nextflowSession"))), esc(f"brca12-wgs {version}")),
        (_tip(t("database"), t("database_tip", vep=versions.get("ensembl-vep") or ui("unknown", lang),
                               spec=sources.get("specification", "ClinGen ENIGMA BRCA1/BRCA2 VCEP v1.2"),
                               gnomad=", ".join(sources.get("gnomad") or []) or t("no"),
                               spliceai=t("yes") if sources.get("spliceai") else t("no"),
                               bayesdel=t("yes") if sources.get("bayesdel") else t("no"))),
         esc(f"VEP {versions.get('ensembl-vep') or '?'} · ClinVar {versions.get('clinvar') or ''}".rstrip())),
        (_tip(t("depth"), t("depth_tip", min=minimum, pct=result["thresholds"].get("minExonPct20x", 0))),
         esc(" · ".join(f"{gene} {v[1]:g}x" for gene, v in depth.items()) or ui("na", lang))),
        (_tip(t("quality"), t("quality_tip", min=minimum)),
         f"<em class='ok'>{esc(t('passed'))}</em>" if passed else
         f"<em class='warn'>{esc(t('limited'))}</em>: {esc(reason_text)}"),
    ]
    parts = [f"<main lang='{lang}' class='l-{lang}'>",
             f"<div class='ruo'>{esc(t('ruo'))}</div>",
             f"<div class='eyebrow'>{esc(t('eyebrow'))}</div><h1>{esc(result['sample'])}</h1>",
             "<div class='head'>" + "".join(f"<div><span>{k}</span><b>{v}</b></div>" for k, v in cells) + "</div>",
             f"<nav class='tabs seg'><button type='button' data-tab='r'>{esc(t('tab_r'))}</button>"
             f"<button type='button' data-tab='d'>{esc(t('tab_d'))}</button></nav>",
             "<section class='tab-r'>"]
    structural = (summary.get("cnvCandidateCount", len(result["cnvs"])) + summary.get("exonCopyNumberCallCount", 0)
                  + summary.get("svCount", len(result["svs"])) + summary.get("insertionSignalCount", 0))
    cards = [
        (summary.get("smallVariantCount", 0), _tip(t("card_small"), t("card_small_tip"))),
        (summary.get("clinvarPathogenicCount", 0),
         _tip(t("card_plp"), t("card_plp_tip", lof=summary.get("lossOfFunctionCandidateCount", 0)))),
        (structural, _tip(t("card_struct"), t("card_struct_tip", cnv=summary.get("cnvCandidateCount", len(result["cnvs"])),
                                               exon=summary.get("exonCopyNumberCallCount", 0),
                                               sv=summary.get("svCount", len(result["svs"])),
                                               ins=summary.get("insertionSignalCount", 0)))),
        (sum(v[3] for v in depth.values()),
         _tip(t("card_exons"), t("card_exons_tip", total=sum(v[0] for v in depth.values())))),
    ]
    parts.append("<div class='cards'>" + "".join(f"<div class='card'><b>{n}</b><span>{label}</span></div>"
                                                 for n, label in cards) + "</div>")
    level, count = verdict(result)
    mosaic = t("v_mosaic") if (result.get("mosaicCalling") or {}).get("enabled") else ""
    text = {"none": t("v_none"), "unclear": t("v_unclear", n=count), "attention": t("v_attention", n=count)}[level]
    parts.append(f"<h2>{esc(t('verdict_h'))}</h2><div class='verdict v-{level}'><b>{esc(text)}</b>"
                 f"<p>{esc(t('v_checked', mosaic=mosaic))}</p>"
                 + ("" if passed else f"<p class='warn'>{esc(t('v_limited', reasons=reason_text))}</p>") + "</div>")
    parts.append(f"<h2>{esc(t('cancer_h'))}</h2><div class='panel'>"
                 + _table_html(t("cancer_headers"), cancer_rows(result, lang), css="cancer") + "</div>"
                 f"<p class='note'>{esc(t('cancer_note'))}</p><p class='note'>{esc(t('cancer_refs'))}</p>")
    for gene, values in depth.items():
        if values[1] < minimum:
            parts.append(f"<div class='alert'>{esc(t('low_depth', gene=gene, depth=values[1], min=minimum))}</div>")
    for flag in (result.get("readQC") or {}).get("flags", []):
        parts.append(f"<div class='alert'>{esc(str(tr(flag, lang)))}</div>")
    parts.append(f"<h2>{esc(t('notable'))}</h2>")
    rows = notable_rows(result, lang)
    if rows:
        parts.append("<div class='panel'>" + _table_html(t("notable_headers"), rows) + "</div>"
                     f"<p class='note'>{esc(t('notable_note'))}</p>")
    else:
        parts.append(f"<div class='panel prose'><p>{esc(t('notable_none'))}</p></div>")
    if depth:
        parts.append(f"<h2>{esc(t('coverage'))}</h2><div class='panel'>"
                     + _table_html(t("cov_headers"), [[g, v[0], f"{v[1]:g}x", f"{v[2]:g}x", v[3]] for g, v in depth.items()])
                     + "</div>")
    clinvar_release = versions.get("clinvar") or ""
    if len(clinvar_release) == 8 and clinvar_release.isdigit():
        clinvar_release = f"{clinvar_release[:4]}-{clinvar_release[4:6]}-{clinvar_release[6:]}"
    parts.append(f"<h2>{esc(t('how_h'))}</h2><div class='panel prose'>"
                 + "".join(f"<p>{esc(x)}</p>" for x in (
                     t("how_sources", caller=result["snvCaller"], vep=versions.get("ensembl-vep") or ui("unknown", lang),
                       clinvar=clinvar_release or ui("unknown", lang),
                       gnomad=", ".join(sources.get("gnomad") or []) or "VEP cache"),
                     t("how_attention"), t("how_unclear"), t("how_rule", min=minimum),
                     t("how_trace", url=RESEARCH_URL)))
                 + "</div>")
    parts.append(f"<h2>{esc(t('method'))}</h2><div class='panel prose'><p>{esc(t('method_text', caller=result['snvCaller']))}</p>"
                 + "".join(f"<p>{esc(item)}</p>" for item in limitations)
                 + f"<p>{esc(disclaimer)}</p></div>")
    parts.append("</section><section class='tab-d'>")
    parts.append(f"<h2>{esc(t('tab_d'))}</h2><p class='note'>{esc(t('details_note'))}</p>")
    uncertain = tr("unsicher", lang)
    run_titles = {ui("qc_title", lang, assay=result.get("assay", "wgs"))}
    entries = notable(result)
    for gene in GENES:
        info = result["genes"][gene]
        g = _gene_result(result, gene)
        n_small = len(info["smallVariants"])
        d = depth.get(gene)
        strand = ui("strand_minus" if info["strand"] == "-" else "strand_plus", lang)
        parts.append(f"<details><summary><b>{gene}</b> <span class='note'>"
                     + esc(t("gene_sum", small=n_small, notable=sum(1 for e in entries if e[1] == gene),
                             uncertain=d[3] if d else 0)) + "</span></summary>")
        parts.append("<div class='kv'>" + esc(t("gene_kv", chrom=info["chrom"], strand=strand, transcript=info["transcript"],
                                                  depth=f"{d[1]:g}x" if d else ui("na", lang),
                                                  uncertain=d[3] if d else 0, small=n_small)) + "</div>")
        for svg in (figures.coverage_svg(result, gene, lang), figures.dosage_svg(result, gene, lang)):
            if svg:
                parts.append(f"<figure class='report-figure'>{svg}</figure>")
        for title, table, empty in sections_for(g, lang):
            if title in run_titles or (title.startswith(GENES[1 - GENES.index(gene)] + " ")):
                continue
            parts.append(f"<h3>{esc(title)}</h3>")
            if table["rows"]:
                parts.append(_table_html(table["headers"], table["rows"], uncertain, css="dt"))
            else:
                parts.append(f"<p class='note'>{esc(empty)}</p>")
        parts.append("</details>")
    parts.append(f"<details><summary><b>{esc(t('run'))}</b> <span class='note'>{esc(t('run_sum'))}</span></summary>")
    for title, table, empty in sections_for(result, lang):
        if title in run_titles or title.startswith(ui("weak_title", lang, depth=0).split("(")[0]):
            parts.append(f"<h3>{esc(title)}</h3>")
            parts.append(_table_html(table["headers"], table["rows"], css="dt") if table["rows"]
                         else f"<p class='note'>{esc(empty)}</p>")
    parts.append(f"<h3>{esc(ui('versions', lang))}</h3>"
                 + _table_html(["Tool", ui("version", lang)], [[k, v or ui("unknown", lang)] for k, v in versions.items()], css="dt"))
    parts.append(f"<h3>{esc(ui('commands', lang))}</h3><p class='note'>{esc(ui('commands_note', lang))} "
                 f"<code>brca12_commands.sh</code>.</p><pre id='cmd-{lang}'>{esc(result['invocation'])}</pre>"
                 f"<button type='button' class='copy' data-copy='cmd-{lang}'>{esc(t('copy'))}</button></details>")
    parts.append("</section>")
    parts.append(f"<footer><p class='cite'>{esc(t('cite'))}: Jungwirth P. brca12-wgs: BRCA1/BRCA2 germline analysis from "
                 f"whole-genome sequencing (GRCh38), version {esc(str(version))}. Zenodo. "
                 f"<a href='https://doi.org/{esc(doi)}'>https://doi.org/{esc(doi)}</a></p>"
                 f"<p>{esc(t('research'))}: <a href='{RESEARCH_URL}'>{RESEARCH_URL}</a> · {esc(t('source'))}: "
                 f"<a href='{SOURCE_URL}'>{SOURCE_URL}</a></p>"
                 f"<p>{esc(t('files'))}</p><p>{esc(disclaimer)}</p></footer></main>")
    return "\n".join(parts)


def html_page(result, lang, sections_for, disclaimer_for, limitations_for):
    """Whole page. `lang` is the language shown first; the other one is one click away."""
    esc = html.escape
    first = "de" if lang == "de" else "en"
    parts = [f"<!DOCTYPE html>\n<html lang='{first}'><head><meta charset='utf-8'>"
             "<meta name='viewport' content='width=device-width,initial-scale=1'>"
             f"<title>brca12-wgs · {esc(result['sample'])}</title>",
             f"<style>\n{_fonts()}\n{CSS}\n@page{{@bottom-left{{content:'{_css_text(result['sample'])} · "
             f"{report_number(result)} · brca12-wgs · RUO'}}}}\n</style></head>",
             f"<body data-lang='{first}' data-tab='r'>",
             "<header class='bar'><div class='bar-in'><div class='brand'>PJ Labs <i>·</i> brca12-wgs</div>"
             "<nav class='lang seg'><button type='button' data-set='en'>English</button>"
             "<button type='button' data-set='de'>Deutsch</button></nav></div></header>"]
    for code in (first, "en" if first == "de" else "de"):
        parts.append(_main(result, code, sections_for, disclaimer_for(code), limitations_for(result, code)))
    parts.append(f"<script>{SCRIPT}</script></body></html>")
    return "\n".join(parts)
