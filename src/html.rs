//! `PREFIX.html`: the result as one self-contained, interactive page for reading and printing, in
//! the form of the hla-wgs page (IBM Plex Mono embedded, German and English in one file with a
//! switch, tabs for the result and the details). It is built from the same values as PREFIX.json.
//! Default language from BRCA12_WGS_LANG, LC_ALL or LANG. Research use only.
use serde_json::Value;

pub const SOURCE: &str = "https://github.com/pjlabsdev/brca12-wgs";
const CSS: &str = include_str!("../assets/page.css");
// IBM Plex Mono (SIL Open Font License 1.1, assets/fonts/OFL.txt), embedded so the page needs no network
const FONT400: &[u8] = include_bytes!("../assets/fonts/IBMPlexMono-Latin-400.woff2");
const FONT500: &[u8] = include_bytes!("../assets/fonts/IBMPlexMono-Latin-500.woff2");

fn esc(s: &str) -> String {
    s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;").replace('"', "&quot;")
}

fn german() -> bool {
    let v = std::env::var("BRCA12_WGS_LANG").or_else(|_| std::env::var("LC_ALL")).or_else(|_| std::env::var("LANG")).unwrap_or_default();
    v.to_ascii_lowercase().starts_with("de")
}

fn b64(d: &[u8]) -> String {
    const A: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut o = String::with_capacity(d.len() / 3 * 4 + 4);
    for c in d.chunks(3) {
        let n = (c[0] as u32) << 16 | (*c.get(1).unwrap_or(&0) as u32) << 8 | *c.get(2).unwrap_or(&0) as u32;
        for i in 0..4 {
            o.push(if i <= c.len() { A[(n >> (18 - 6 * i) & 63) as usize] as char } else { '=' });
        }
    }
    o
}

/// calendar date (UTC) from seconds since 1970
fn date(secs: u64) -> String {
    let z = (secs / 86_400) as i64 + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m = (if mp < 10 { mp + 3 } else { mp - 9 }) as u32;
    format!("{:04}-{m:02}-{d:02}", yoe + era * 400 + i64::from(m <= 2))
}

/// FNV-1a over the results (genes, exons, variants, structural signatures) and the input names:
/// same result, same number
fn report_id(rep: &Value) -> String {
    let text = format!("{}{}{}{}", rep["genes"], rep["exons"], rep["variants"], rep["input"]["files"]);
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for b in text.bytes() {
        h ^= b as u64;
        h = h.wrapping_mul(0x0100_0000_01b3);
    }
    format!("B12-{:08X}", (h >> 32) as u32 ^ h as u32)
}

fn tip(label: &str, text: &str) -> String {
    format!("<span class=\"tip\" tabindex=\"0\" data-tip=\"{}\">{}</span>", esc(text), label)
}

fn s(v: &Value) -> String {
    match v {
        Value::Null => "–".into(),
        Value::String(x) if x.is_empty() => "–".into(),
        Value::String(x) => esc(x),
        x => esc(&x.to_string()),
    }
}

fn num(v: &Value, d: usize) -> String {
    v.as_f64().map(|x| format!("{x:.d$}")).unwrap_or_else(|| "–".into())
}

/// the page holds the result in German and in English; a switch at the top shows one of them
pub fn render(rep: &Value) -> String {
    let css = format!(
        "@font-face{{font-family:'IBM Plex Mono';font-weight:400;font-style:normal;font-display:swap;src:url(data:font/woff2;base64,{}) format('woff2')}}\n@font-face{{font-family:'IBM Plex Mono';font-weight:500;font-style:normal;font-display:swap;src:url(data:font/woff2;base64,{}) format('woff2')}}\n{}",
        b64(FONT400),
        b64(FONT500),
        CSS.lines().filter(|l| !l.starts_with("@font-face")).collect::<Vec<_>>().join("\n")
    );
    let first = if german() { "de" } else { "en" };
    let sample = rep["sample"].as_str().unwrap_or("sample");
    format!(
        "<!DOCTYPE html>\n<html lang=\"{first}\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">\n<title>brca12-wgs · {}</title>\n<style>{css}</style></head><body data-lang=\"{first}\" data-tab=\"r\">\n<header class=\"bar\"><div class=\"bar-in\"><div class=\"brand\">PJ Labs <i>·</i> brca12-wgs</div><nav class=\"lang seg\"><button type=\"button\" data-set=\"de\">Deutsch</button><button type=\"button\" data-set=\"en\">English</button></nav></div></header>\n<main lang=\"de\" class=\"l-de\">\n{}</main>\n<main lang=\"en\" class=\"l-en\">\n{}</main>\n<script>\n{SCRIPT}</script></body></html>\n",
        esc(sample),
        render_one(rep, true),
        render_one(rep, false)
    )
}

const SCRIPT: &str = "function aria(){document.querySelectorAll('nav.lang button').forEach(b=>b.setAttribute('aria-pressed',b.dataset.set===document.body.dataset.lang));document.querySelectorAll('nav.tabs button').forEach(b=>b.setAttribute('aria-pressed',b.dataset.tab===document.body.dataset.tab))}
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
aria();
addEventListener('beforeprint',()=>document.querySelectorAll('details').forEach(d=>d.open=true));
";

fn result_label(r: &str, de: bool) -> &'static str {
    match (r, de) {
        ("pathogenic", true) => "pathogen",
        ("pathogenic", false) => "pathogenic",
        ("candidate", true) => "zu bestätigen",
        ("candidate", false) => "to confirm",
        ("none", true) => "kein Befund",
        ("none", false) => "no finding",
        ("not_assessable", true) => "nicht beurteilbar",
        ("not_assessable", false) => "not assessable",
        (_, true) => "nicht bestimmt",
        (_, false) => "not determined",
    }
}

fn class_label(c: &str, de: bool) -> &'static str {
    match (c, de) {
        ("pathogenic", true) => "pathogen (ClinVar)",
        ("pathogenic", false) => "pathogenic (ClinVar)",
        ("candidate", true) => "Kandidat",
        ("candidate", false) => "candidate",
        ("benign", true) => "gutartig",
        ("benign", false) => "benign",
        ("known_common", true) => "häufiges Allel",
        ("known_common", false) => "common allele",
        ("signal", true) => "Signal",
        ("signal", false) => "signal",
        (_, true) => "übrige",
        (_, false) => "other",
    }
}

fn render_one(rep: &Value, de: bool) -> String {
    let t = |d: &'static str, e: &'static str| if de { d } else { e };
    let l = |v: &Value| -> String { s(&v[if de { "de" } else { "en" }]) };
    let sample = rep["sample"].as_str().unwrap_or("sample");
    let genes: Vec<&Value> = rep["genes"].as_array().map(|a| a.iter().collect()).unwrap_or_default();
    let exons: Vec<&Value> = rep["exons"].as_array().map(|a| a.iter().collect()).unwrap_or_default();
    let variants: Vec<&Value> = rep["variants"].as_array().map(|a| a.iter().collect()).unwrap_or_default();
    let ecn = &rep["exon_copy_number"];
    let mut h = String::with_capacity(120_000);
    h.push_str(&format!("<div class=\"ruo\">{}</div>\n", t(
        "Nur für Forschung (RUO). Kein medizinischer Befund und kein Medizinprodukt.",
        "Research use only (RUO). Not a medical report and not a medical device.")));
    h.push_str(&format!("<div class=\"eyebrow\">{}</div><h1>{}</h1>\n", t("BRCA1 und BRCA2 aus dem Genom", "BRCA1 and BRCA2 from the genome"), esc(sample)));

    // report head
    let uncertain: usize = exons.iter().filter(|e| e["status"] != "ok").count();
    let mut limits: Vec<&str> = vec![];
    if genes.iter().any(|g| g["result"] == "not_assessable") {
        limits.push(t("Gen nicht beurteilbar", "gene not assessable"));
    }
    if uncertain > 0 {
        limits.push(t("unsichere Exons", "uncertain exons"));
    }
    if ecn["enabled"] != true {
        limits.push(if ecn["gc_slope_per_10pct"].is_number() { t("Exon-Kopienzahl nicht beurteilbar", "exon copy number not assessable") }
                    else { t("Exon-Kopienzahl nicht bestimmt", "exon copy number not determined") });
    }
    let qc = if limits.is_empty() {
        format!("<span class=\"ok\">{}</span>", t("bestanden", "passed"))
    } else {
        format!("<span class=\"warn\">{}: {}</span>", t("eingeschränkt", "limited"), limits.join(", "))
    };
    let files = rep["input"]["files"].as_array().map(|a| {
        a.iter().filter_map(|x| x.as_str()).map(|p| esc(std::path::Path::new(p).file_name().and_then(|n| n.to_str()).unwrap_or(p))).collect::<Vec<_>>().join(", ")
    }).unwrap_or_else(|| "–".into());
    let cov: Vec<String> = genes.iter().map(|g| format!("{} {}x", s(&g["gene"]), num(&g["coverage"]["coding_mean_depth"], 1))).collect();
    let src = &rep["annotation_sources"];
    let kv: Vec<(&str, String, String)> = vec![
        (t("Probe", "Sample"), esc(sample), t("Name der Probe, aus dem Namen der Ausgabedateien.", "Name of the sample, from the name of the output files.").into()),
        (t("Berichtsnummer", "Report number"), report_id(rep), t("Prüfsumme (FNV-1a) über alle Ergebnisse und die Namen der Eingabedateien: gleiches Ergebnis, gleiche Nummer.", "Checksum (FNV-1a) over all results and the names of the input files: same result, same number.").into()),
        (t("Eingabe", "Input"), files, t("Eingabedateien ohne Ordner: BAM/CRAM oder ein FASTQ-Paar (dann ausgerichtet zu PREFIX.bam).", "Input files without folders: BAM/CRAM or a FASTQ pair (then aligned to PREFIX.bam).").into()),
        (t("Analyse", "Analysis"), date(rep["analysis_time"].as_u64().unwrap_or(0)), t("Datum des Laufs (UTC).", "Date of the run (UTC).").into()),
        (t("Programm", "Software"), format!("brca12-wgs {}", s(&rep["version"])), t("Ein Programm ohne externe Werkzeuge: Abdeckung, Variantenaufruf, Annotation, Exon-Kopienzahl und Strukturvarianten sind eingebaut.", "One program without external tools: coverage, variant calling, annotation, exon copy number and structural variants are built in.").into()),
        (t("Annotation", "Annotation"), "MANE · ClinVar · gnomAD".into(), format!("{} · {} · {}", s(&src["transcripts"]), s(&src["clinvar"]), s(&src["gnomad"]))),
        (t("Exon-Tiefe", "Exon depth"), cov.join(" · "), t("Mittlere Tiefe über die Exons (Mapping-Qualität ab 20), nach Exonlänge gewichtet.", "Mean depth over the exons (mapping quality 20 or more), weighted by exon length.").into()),
        (t("Qualität", "Quality"), qc, t("Bestanden, wenn beide Gene beurteilbar sind, kein Exon unsicher ist und die Exon-Kopienzahl bestimmt wurde. Sonst eingeschränkt, mit Grund.", "Passed when both genes are assessable, no exon is uncertain and the exon copy number was determined. Otherwise limited, with the reason.").into()),
    ];
    h.push_str("<div class=\"head\">");
    for (k, v, tp) in kv {
        h.push_str(&format!("<div><span>{}</span><b>{v}</b></div>", tip(k, &tp)));
    }
    h.push_str("</div>\n");
    for g in &genes {
        if g["result"] == "pathogenic" {
            h.push_str(&format!("<div class=\"alert\">{}</div>\n", l(&g["statement"])));
        }
    }
    for g in &genes {
        if g["result"] == "not_assessable" {
            h.push_str(&format!("<div class=\"alert soft\">{}</div>\n", l(&g["statement"])));
        }
    }
    if ecn["enabled"] != true && ecn["gc_slope_per_10pct"].is_number() {
        h.push_str(&format!("<div class=\"alert soft\">{}: {}</div>\n", t(
            "Exon-Kopienzahl nicht beurteilbar, mögliche Exon-Deletionen oder -Duplikationen werden nicht gemeldet",
            "Exon copy number not assessable, possible exon deletions or duplications are not reported"),
            if de { s(&ecn["reason_de"]) } else { s(&ecn["reason"]) }));
    } else if ecn["enabled"] != true {
        h.push_str(&format!("<div class=\"alert soft\">{} ({})</div>\n", t(
            "Exon-Kopienzahl nicht bestimmt: Die Kontrollfenster des Referenzpanels fehlen in der Eingabe, etwa bei einer BAM nur der Genregionen.",
            "Exon copy number not determined: the control windows of the reference panel are missing in the input, for example in a BAM of the gene regions only."),
            if de && ecn["reason"] == "too few control windows with depth" { "zu wenige Kontrollfenster mit Tiefe".to_string() } else { s(&ecn["reason"]) }));
    }

    h.push_str(&format!("<nav class=\"tabs seg\"><button type=\"button\" data-tab=\"r\">{}</button><button type=\"button\" data-tab=\"d\">{}</button></nav>\n<section class=\"tab-r\">\n",
        t("Ergebnis", "Result"), t("Details", "Details")));

    // cards
    let ok_exons = exons.len() - uncertain;
    let pass_in_genes: u64 = genes.iter().map(|g| g["small_variants"].as_u64().unwrap_or(0)).sum();
    h.push_str("<div class=\"cards\">");
    for g in &genes {
        let r = g["result"].as_str().unwrap_or("");
        h.push_str(&format!("<div class=\"card {r}\"><b>{}</b><span>{}</span></div>", result_label(r, de),
            tip(&format!("{} {}", s(&g["gene"]), t("Ergebnis", "result")), &if de { "Pathogen: ClinVar-pathogene Variante oder bekannte pathogene Gründermutation. Zu bestätigen: Funktionsverlust-Kandidat, Exon-Deletion oder -Duplikation, Strukturvariante in einem Exon. Kein Befund: nichts davon bei ausreichender Abdeckung." } else { "Pathogenic: ClinVar pathogenic variant or known pathogenic founder allele. To confirm: loss-of-function candidate, exon deletion or duplication, structural variant in an exon. No finding: none of these with enough coverage." })));
    }
    h.push_str(&format!("<div class=\"card\"><b>{ok_exons} / {}</b><span>{}</span></div>", exons.len(),
        tip(t("Exons sicher abgedeckt", "exons covered safely"), t("Exons mit mittlerer Tiefe ab 20x, mindestens 95 % der Basen ab 20x und höchstens 10 % Reads mit Mapping-Qualität unter 20.", "Exons with mean depth of 20x or more, at least 95 % of bases at 20x or more and at most 10 % of reads with mapping quality below 20."))));
    h.push_str(&format!("<div class=\"card\"><b>{pass_in_genes}</b><span>{}</span></div></div>\n",
        tip(t("Varianten in den Genen", "variants in the genes"), t("SNVs und Indels mit FILTER PASS in den Genregionen, gutartige und häufige eingeschlossen.", "SNVs and indels with FILTER PASS in the gene regions, benign and common ones included."))));

    // statement per gene
    h.push_str(&format!("<h2>{}</h2>\n<div class=\"panel\"><table class=\"st\"><thead><tr><th>{}</th><th>{}</th><th>{}</th></tr></thead><tbody>\n",
        t("Gesamtaussage je Gen", "Overall statement per gene"), t("Gen", "Gene"), t("Ergebnis", "Result"), t("Aussage", "Statement")));
    for g in &genes {
        let r = g["result"].as_str().unwrap_or("");
        h.push_str(&format!("<tr><td class=\"g\">{}<span class=\"full\">{}</span></td><td><span class=\"res {r}\">{}</span></td><td>{}</td></tr>\n",
            s(&g["gene"]), s(&g["transcript"]["refseq"]), result_label(r, de), l(&g["statement"])));
    }
    h.push_str("</tbody></table></div>\n");

    // findings
    let mut frows = String::new();
    for v in variants.iter().filter(|v| v["in_gene_region"] == true && (v["class"] == "pathogenic" || v["class"] == "candidate")) {
        let p = v["hgvs_p"].as_str().map(|p| format!(" p.({})", esc(p.split_once(":p.").map(|x| x.1).unwrap_or(p)))).unwrap_or_default();
        let cv = &v["clinvar"];
        let ev = if cv.is_null() { t("nicht in ClinVar", "not in ClinVar").to_string() } else {
            format!("ClinVar {} ({}★)", s(&cv["significance"]).replace('_', " "), s(&cv["stars"]))
        };
        frows.push_str(&format!("<tr><td class=\"g\">{}</td><td class=\"wrap\"><span class=\"call\">{}{p}</span><span class=\"full\">{}:{} {}&gt;{}</span></td><td>{}</td><td>{}<span class=\"full\">{} {} · {} {}</span></td><td><span class=\"res {}\">{}</span>{}</td></tr>\n",
            s(&v["gene"]), s(&v["hgvs_c"]), s(&v["chrom"]), s(&v["pos"]), s(&v["ref"]), s(&v["alt"]),
            v["consequence"].as_array().map(|a| a.iter().filter_map(|x| x.as_str()).map(|x| x.replace('_', " ")).collect::<Vec<_>>().join(", ")).unwrap_or_default(),
            ev, s(&v["zygosity"]), s(&v["genotype"]), t("Reads", "reads"), format!("{}/{}", s(&v["alt_reads"]), s(&v["depth"])),
            s(&v["class"]), class_label(v["class"].as_str().unwrap_or(""), de),
            if v["filter"] != "PASS" { format!("<span class=\"full warn\">FILTER {}</span>", s(&v["filter"])) } else { String::new() }));
    }
    for g in &genes {
        if let Some(calls) = g["exon_copy_number"].as_array() {
            for c in calls {
                let ex = c["affected_exons"].as_array().map(|a| a.iter().filter_map(|x| x.as_str()).collect::<Vec<_>>().join(", ")).unwrap_or_default();
                let del = c["type"] == "deletion";
                frows.push_str(&format!("<tr><td class=\"g\">{}</td><td class=\"wrap\"><span class=\"call\">{}</span><span class=\"full\">{}</span></td><td>{}</td><td>{} {} · {} {}</td><td><span class=\"res candidate\">{}</span></td></tr>\n",
                    s(&c["gene"]), if de { if del { "Exon-Deletion" } else { "Exon-Duplikation" } } else if del { "exon deletion" } else { "exon duplication" },
                    esc(&ex), t("Exon-Kopienzahl gegen 80 Referenzgenome", "exon copy number against 80 reference genomes"),
                    t("Kopienzahl", "copy number"), s(&c["copy_number"]), t("Qualität", "quality"), s(&c["quality"]), t("bestätigen (MLPA)", "confirm (MLPA)")));
            }
        }
        if let Some(svs) = g["structural_variants"].as_array() {
            for x in svs.iter().filter(|x| x["class"] == "pathogenic" || x["class"] == "candidate") {
                let ex = x["affected_exons"].as_array().map(|a| a.iter().filter_map(|y| y.as_str()).collect::<Vec<_>>().join(", ")).unwrap_or_default();
                frows.push_str(&format!("<tr><td class=\"g\">{}</td><td class=\"wrap\"><span class=\"call\">{}</span><span class=\"full\">{}:{}-{}{}</span></td><td>{}</td><td class=\"wrap\">{}</td><td><span class=\"res {}\">{}</span></td></tr>\n",
                    s(&x["gene"]), l(&x["type"]), s(&x["chrom"]), s(&x["start"]), s(&x["end"]), if ex.is_empty() { String::new() } else { format!(" · {}", esc(&ex)) },
                    t("Lesepaare, geteilte Lesungen, Neuausrichtung", "read pairs, split reads, realignment"), s(&x["evidence"]),
                    s(&x["class"]), class_label(x["class"].as_str().unwrap_or(""), de)));
            }
        }
    }
    h.push_str(&format!("<h2>{}</h2>\n", t("Befunde", "Findings")));
    if frows.is_empty() {
        h.push_str(&format!("<div class=\"panel prose\"><p>{}</p></div>\n", t(
            "Keine pathogene Variante und kein Befund, der bestätigt werden muss, in BRCA1 und BRCA2.",
            "No pathogenic variant and no finding that needs confirmation in BRCA1 and BRCA2.")));
    } else {
        h.push_str(&format!("<div class=\"panel\"><table><thead><tr><th>{}</th><th>{}</th><th>{}</th><th>{}</th><th>{}</th></tr></thead><tbody>\n{frows}</tbody></table></div>\n<p class=\"note\">{}</p>\n",
            t("Gen", "Gene"), t("Befund", "Finding"), t("Art", "Kind"), t("Belege", "Evidence"), t("Einstufung", "Class"),
            t("Jeder Befund braucht vor einer medizinischen Entscheidung die Bestätigung in einem akkreditierten Labor (Sanger-Sequenzierung, MLPA oder eine andere orthogonale Methode) und eine humangenetische Beratung.",
              "Every finding needs confirmation in an accredited laboratory (Sanger sequencing, MLPA or another orthogonal method) and genetic counselling before any medical decision.")));
    }

    // cancer types
    let cr = &rep["cancer_risks"];
    if let Some(cancers) = cr["cancers"].as_array() {
        h.push_str(&format!("<h2>{}</h2>\n<div class=\"panel\"><table class=\"clin\"><thead><tr><th>{}</th><th>BRCA1</th><th>BRCA2</th></tr></thead><tbody>\n",
            t("Krebsarten", "Cancer types"), t("Krebsart", "Cancer type")));
        for c in cancers {
            h.push_str(&format!("<tr><td class=\"g\">{}</td>", l(&c["name"])));
            for g in c["genes"].as_array().into_iter().flatten() {
                let srcs = g["sources"].as_array().map(|a| a.iter().filter_map(|x| x.as_str()).map(|k| {
                    let so = &cr["sources"][k];
                    format!("{}{}", so["citation"].as_str().unwrap_or(k), so["pmid"].as_str().map(|p| format!(" PMID {p}")).unwrap_or_default())
                }).collect::<Vec<_>>().join("; ")).unwrap_or_default();
                h.push_str(&format!("<td>{}<span class=\"full\">{}: {}</span><span class=\"note\">{}</span></td>",
                    tip(&l(&g["association_label"]), &srcs), t("Risiko für Träger", "risk for carriers"), l(&g["risk_for_carriers"]), l(&g["interpretation"])));
            }
            h.push_str("</tr>\n");
        }
        h.push_str(&format!("</tbody></table></div>\n<p class=\"note\">{} {}</p>\n", l(&cr["other_genes_note"]), l(&cr["disclaimer"])));
    }

    // exons
    h.push_str(&format!("<h2>{}</h2>\n", t("Exons", "Exons")));
    for g in &genes {
        let gene = g["gene"].as_str().unwrap_or("");
        let mine: Vec<&&Value> = exons.iter().filter(|e| e["gene"] == gene).collect();
        let unc = mine.iter().filter(|e| e["status"] != "ok").count();
        let mut rows = String::new();
        for e in &mine {
            let cn = &e["copy_number"];
            let reasons = e["reasons"][if de { "de" } else { "en" }].as_array().map(|a| a.iter().filter_map(|x| x.as_str()).map(esc).collect::<Vec<_>>().join("; ")).unwrap_or_default();
            rows.push_str(&format!("<tr><td class=\"mono\">{}</td><td class=\"num\">{}</td><td class=\"num\">{}</td><td class=\"num\">{}</td><td class=\"num\">{}</td><td>{}{}</td></tr>\n",
                s(&e["label"]), num(&e["mean_depth"], 1), num(&e["pct20x"], 1), num(&e["dosage_ratio"], 2),
                if cn.is_null() { "–".into() } else { s(&cn["copy_number"]) },
                if e["status"] == "ok" { format!("<span class=\"ok\">{}</span>", t("sicher", "ok")) } else { format!("<span class=\"warn\">{}</span>", t("unsicher", "uncertain")) },
                if reasons.is_empty() { String::new() } else { format!("<span class=\"full\">{reasons}</span>") }));
        }
        h.push_str(&format!("<details><summary><b>{gene}</b> <span class=\"note\">{} {}, {} {} · {} {}</span></summary><table><thead><tr><th>{}</th><th class=\"num\">{}</th><th class=\"num\">≥20x %</th><th class=\"num\">{}</th><th class=\"num\">{}</th><th>{}</th></tr></thead><tbody>\n{rows}</tbody></table></details>\n",
            mine.len(), t("Exons", "exons"), unc, t("unsicher", "uncertain"), num(&g["coverage"]["coding_mean_depth"], 1), t("x im Mittel", "x mean"),
            t("Exon", "Exon"), t("Tiefe", "Depth"), tip(t("Dosis", "Dosage"), t("Exontiefe relativ zum Median des Chromosoms (100-kb-Fenster).", "Exon depth relative to the median of the chromosome (100 kb windows).")),
            tip(t("Kopien", "Copies"), t("Kopienzahl des Exons im Hidden-Markov-Modell gegen 80 Referenzgenome; 2 ist normal.", "Copy number of the exon in the hidden Markov model against 80 reference genomes; 2 is normal.")),
            t("Status", "Status")));
    }

    // how the statements come about
    let rq = &rep["read_qc"];
    let steps: Vec<String> = vec![
        format!("{} {} · {} {} %, {} {} %, {} {}.", t("Reads aus", "Reads from"), files_plain(rep), t("davon ausgerichtet", "mapped"), num(&rq["mapped_pct"], 2), t("Duplikate", "duplicates"), num(&rq["duplicate_pct"], 2), t("mittlere Insertgröße", "mean insert size"), num(&rq["mean_insert_size"], 0)),
        t("Abdeckung je Exon mit Mapping-Qualität ab 20, die Überlappung zweier Lesepaar-Hälften einmal gezählt; Exons unter 20x im Mittel oder unter 95 % der Basen ab 20x gelten als unsicher.", "Coverage per exon with mapping quality 20 or more, the overlap of two mates counted once; exons below 20x mean or below 95 % of bases at 20x count as uncertain.").into(),
        t("Kleine Varianten (SNVs, Indels) mit dem eigenen Varianten-Caller: lokale Neuassemblierung, Paar-HMM, Genotypisierung je Allel wie GATK HaplotypeCaller; Filter QD und FS.", "Small variants (SNVs, indels) with the own variant caller: local reassembly, pair-HMM, per-allele genotyping as GATK HaplotypeCaller; filters QD and FS.").into(),
        format!("{} {} · {}.", t("Annotation: HGVS c./p. und Konsequenz auf den MANE-Select-Transkripten (gegen Ensembl VEP 116 geprüft, 99,97 % gleich), ClinVar und die gnomAD-Filterhäufigkeit FAF95. Quellen:", "Annotation: HGVS c./p. and consequence on the MANE Select transcripts (checked against Ensembl VEP 116, 99.97 % equal), ClinVar and the gnomAD filtering allele frequency FAF95. Sources:"), s(&src["clinvar"]), s(&src["gnomad"])),
        t("Bewertung: Als pathogen zählt eine Variante nur mit ClinVar-Einstufung pathogen oder wahrscheinlich pathogen (nicht widersprüchlich) oder als bekannte pathogene Gründermutation. Funktionsverlust (Frameshift, Stopp, Spleißstelle ±1/2, Startverlust), nicht gutartig und seltener als 0,1 % in gnomAD, gilt als Kandidat. Das Programm vergibt keine eigene ACMG-Einstufung.", "Evaluation: a variant counts as pathogenic only with a ClinVar classification pathogenic or likely pathogenic (not conflicting) or as a known pathogenic founder allele. Loss of function (frameshift, stop, splice site +/-1/2, start loss), not benign and rarer than 0.1 % in gnomAD, counts as a candidate. The program gives no own ACMG classification.").into(),
        t("Exon-Deletionen und -Duplikationen: Tiefe je Exon gegen 80 Referenzgenome (1000 Genomes, 30x) mit Hidden-Markov-Modell und Allelbalance; gilt für 30x-PCR-freie Illumina-Genome.", "Exon deletions and duplications: depth per exon against 80 reference genomes (1000 Genomes, 30x) with a hidden Markov model and allele balance; valid for 30x PCR-free Illumina genomes.").into(),
        t("Strukturvarianten: Deletionen, Inversionen und Tandemduplikationen aus Lesepaaren und geteilten Lesungen; Alu-Insertionen und zusätzliche Tandemkopien aus abgeschnittenen Leseenden; bekannte Allele (auch die pathogenen Gründermutationen BRCA2 c.156_157insAlu, BRCA1 Exon-13-Duplikation und 3.835-bp-Deletion) durch Neuausrichtung.", "Structural variants: deletions, inversions and tandem duplications from read pairs and split reads; Alu insertions and extra tandem copies from clipped read ends; known alleles (also the pathogenic founder alleles BRCA2 c.156_157insAlu, BRCA1 exon 13 duplication and 3,835 bp deletion) by realignment.").into(),
    ];
    h.push_str(&format!("<h2>{}</h2>\n<div class=\"panel\"><ol class=\"steps\">{}</ol></div>\n<p class=\"note\">{}</p>\n",
        t("Wie die Aussagen zustande kommen", "How the statements come about"),
        steps.iter().map(|x| format!("<li>{x}</li>")).collect::<String>(),
        t("Prüfung: GIAB HG001 (SNV 269/269, Indel 49/49), 39 HPRC-Genome (Exons ±20 bp: SNV 590/590 ohne Fehlruf), sechs bekannte Träger aus 1000 Genomes (alle erkannt), Spike-ins großer Ereignisse; Einzelheiten in docs/validation-2.md.",
          "Validation: GIAB HG001 (SNV 269/269, indel 49/49), 39 HPRC genomes (exons +/- 20 bp: SNV 590/590 without false call), six known carriers from 1000 Genomes (all found), spike-ins of large events; details in docs/validation-2.md.")));
    h.push_str("</section>\n<section class=\"tab-d\">\n");

    // details: variants per gene
    h.push_str(&format!("<h2>{}</h2>\n<p class=\"note\">{}</p>\n", t("Alle Varianten", "All variants"), t(
        "Alle aufgerufenen Allele der Genregionen ±200 bp, auch die gutartigen und die, die die Filter nicht bestanden haben. Ein Klick öffnet ein Gen; beim Drucken sind alle offen.",
        "All called alleles of the gene regions +/- 200 bp, also the benign ones and those that did not pass the filters. Click a gene to open it; all are open when printed.")));
    for g in &genes {
        let gene = g["gene"].as_str().unwrap_or("");
        let mine: Vec<&&Value> = variants.iter().filter(|v| v["gene"] == gene).collect();
        let mut rows = String::new();
        for v in &mine {
            let cv = &v["clinvar"];
            let gn = &v["gnomad"];
            rows.push_str(&format!("<tr><td class=\"mono\">{}<span class=\"full\">{}&gt;{}</span></td><td class=\"wrap\">{}<span class=\"full\">{}</span></td><td>{}</td><td>{}<span class=\"full\">{}</span></td><td>{}</td><td class=\"num\">{}</td><td>{}{}</td></tr>\n",
                s(&v["pos"]), short(&v["ref"]), short(&v["alt"]),
                v["hgvs_c"].as_str().map(|c| esc(c.split_once(":").map(|x| x.1).unwrap_or(c))).unwrap_or_else(|| "–".into()),
                v["hgvs_p"].as_str().map(|p| esc(p.split_once(":").map(|x| x.1).unwrap_or(p))).unwrap_or_default(),
                v["consequence"].as_array().and_then(|a| a.first()).and_then(|x| x.as_str()).map(|x| x.replace('_', " ")).unwrap_or_default(),
                s(&v["genotype"]), format!("{}/{}", s(&v["alt_reads"]), s(&v["depth"])),
                if cv.is_null() { "–".into() } else { format!("{} {}★", s(&cv["significance"]).replace('_', " "), s(&cv["stars"])) },
                gn["faf95"].as_f64().map(|f| format!("{f:.5}")).unwrap_or_else(|| if gn["present"] == true { t("vorhanden", "present").into() } else { "–".into() }),
                class_label(v["class"].as_str().unwrap_or(""), de),
                if v["filter"] != "PASS" { format!("<span class=\"full warn\">FILTER {}</span>", s(&v["filter"])) } else { String::new() }));
        }
        h.push_str(&format!("<details><summary><b>{gene}</b> <span class=\"note\">{} {}</span></summary><table><thead><tr><th>{}</th><th>HGVS</th><th>{}</th><th>{}</th><th>ClinVar</th><th class=\"num\">gnomAD FAF95</th><th>{}</th></tr></thead><tbody>\n{rows}</tbody></table></details>\n",
            mine.len(), t("Allele", "alleles"), t("Position", "Position"), t("Konsequenz", "Consequence"), t("Genotyp, Reads", "Genotype, reads"), t("Einstufung", "Class")));
    }

    // structural signatures
    let svs: Vec<&Value> = genes.iter().flat_map(|g| g["structural_variants"].as_array().into_iter().flatten()).collect();
    h.push_str(&format!("<h2>{}</h2>\n", t("Strukturvarianten und Insertionssignaturen", "Structural variants and insertion signatures")));
    if svs.is_empty() {
        h.push_str(&format!("<div class=\"panel prose\"><p>{}</p></div>\n", t("Keine Signatur.", "No signature.")));
    } else {
        h.push_str(&format!("<div class=\"panel\"><table><thead><tr><th>{}</th><th>{}</th><th>{}</th><th>{}</th><th>{}</th></tr></thead><tbody>\n", t("Gen", "Gene"), t("Art", "Kind"), t("Ort", "Location"), t("Belege", "Evidence"), t("Einstufung", "Class")));
        for x in svs {
            h.push_str(&format!("<tr><td class=\"g\">{}</td><td class=\"wrap\">{}</td><td class=\"mono\">{}:{}-{}<span class=\"full\">{}</span></td><td class=\"wrap note\">{}</td><td>{}</td></tr>\n",
                s(&x["gene"]), l(&x["type"]), s(&x["chrom"]), s(&x["start"]), s(&x["end"]),
                x["length"].as_i64().map(|n| format!("{n} bp")).unwrap_or_default(), s(&x["evidence"]), class_label(x["class"].as_str().unwrap_or(""), de)));
        }
        h.push_str("</tbody></table></div>\n");
    }

    // read QC and gaps
    h.push_str(&format!("<h2>{}</h2>\n<div class=\"panel\"><table><tbody>\n", t("Lese-Qualität", "Read quality")));
    for (k, key, d) in [
        (t("Reads in den Regionen", "Reads in the regions"), "reads_in_regions", 0), (t("ausgerichtet %", "mapped %"), "mapped_pct", 2),
        (t("richtig gepaart %", "properly paired %"), "properly_paired_pct", 2), (t("Duplikate %", "duplicates %"), "duplicate_pct", 2),
        (t("MAPQ 0 %", "MAPQ 0 %"), "mapq0_pct", 2), (t("Insertgröße Median", "insert size median"), "median_insert_size", 0),
        (t("Fehlerrate", "error rate"), "error_rate", 5), (t("Leselänge", "read length"), "average_read_length", 1),
        (t("Basenqualität", "base quality"), "average_base_quality", 1),
    ] {
        h.push_str(&format!("<tr><td>{k}</td><td class=\"num\">{}</td></tr>\n", num(&rq[key], d)));
    }
    h.push_str("</tbody></table></div>\n");
    for g in &genes {
        let iv = g["coverage"]["low_coverage_intervals"].as_array().cloned().unwrap_or_default();
        let rows: String = iv.iter().map(|i| format!("<tr><td class=\"mono\">{}:{}-{}</td><td class=\"num\">{}</td><td class=\"num\">{}</td></tr>\n", s(&i["chrom"]), s(&i["start"]), s(&i["end"]), s(&i["length"]), num(&i["min_depth"], 0))).collect();
        h.push_str(&format!("<details><summary><b>{}</b> <span class=\"note\">{} {} {} bp {}</span></summary><table><thead><tr><th>{}</th><th class=\"num\">bp</th><th class=\"num\">{}</th></tr></thead><tbody>\n{rows}</tbody></table></details>\n",
            s(&g["gene"]), iv.len(), t("Lücken,", "gaps,"), s(&g["coverage"]["low_coverage_bp"]), t("unter 10x", "below 10x"), t("Bereich", "Region"), t("min. Tiefe", "min. depth")));
    }

    // sources
    h.push_str(&format!("<h2>{}</h2>\n<div class=\"panel prose\">", t("Quellen", "Sources")));
    if let Some(m) = cr["sources"].as_object() {
        for (_, so) in m {
            h.push_str(&format!("<p class=\"src\">{}{}</p>", s(&so["citation"]), so["pmid"].as_str().map(|p| format!(" PMID <a href=\"https://pubmed.ncbi.nlm.nih.gov/{p}/\">{p}</a>")).unwrap_or_default()));
        }
    }
    h.push_str(&format!("<p class=\"src\">{}</p><p class=\"src\">{}</p></div>\n", s(&src["clinvar"]), s(&src["gnomad"])));
    h.push_str("</section>\n");
    h.push_str(&format!("<footer>\n<p>{}: <a href=\"{SOURCE}\">{SOURCE}</a></p>\n<p>{}</p>\n<p>{}</p>\n</footer>\n",
        t("Quelltext", "Source code"),
        t("Alle Einzelwerte stehen in den Dateien .tsv, .variants.tsv, .exons.tsv, .sv.tsv und .json neben dieser Seite.", "All values are in the files .tsv, .variants.tsv, .exons.tsv, .sv.tsv and .json next to this page."),
        t("Nur für Forschung (RUO). Kein Medizinprodukt. Vor jeder medizinischen Entscheidung sind die Bestätigung in einem akkreditierten Labor und eine humangenetische Beratung nötig.",
          "Research use only (RUO). Not a medical device. Any medical decision needs confirmation in an accredited laboratory and genetic counselling.")));
    h
}

fn short(v: &Value) -> String {
    let x = v.as_str().unwrap_or("");
    if x.len() > 12 { format!("{}…({} bp)", esc(&x[..8]), x.len()) } else { esc(x) }
}

fn files_plain(rep: &Value) -> String {
    rep["input"]["files"].as_array().map(|a| {
        a.iter().filter_map(|x| x.as_str()).map(|p| esc(std::path::Path::new(p).file_name().and_then(|n| n.to_str()).unwrap_or(p))).collect::<Vec<_>>().join(", ")
    }).unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn escapes_and_dates() {
        assert_eq!(esc("<a&\"b>"), "&lt;a&amp;&quot;b&gt;");
        assert_eq!(date(0), "1970-01-01");
        assert_eq!(date(1_795_426_512), "2026-11-23");
    }
}
