//! Result files, built like the hla-wgs output: `PREFIX.tsv` (one line per gene), `PREFIX.exons.tsv`
//! (one line per exon with its certainty), `PREFIX.json` (schema genome.brca12.v2). The HTML page is
//! made from PREFIX.json at the end. Columns that a version does not determine yet are `-` in the
//! TSV files and `null` in the JSON (listed in `analysed`). Research use only.
use crate::exons::{self, Settings};
use crate::readstats;
use crate::reads;
use serde_json::{Value, json};
use std::fs::File;
use std::io::{self, Write};

pub const SCHEMA: &str = "genome.brca12.v2";
pub const GENES: [(&str, &str, &str); 2] = [
    ("BRCA1", "NM_007294.4", "ENST00000357654"),
    ("BRCA2", "NM_000059.4", "ENST00000380152"),
];

type Bed = Vec<(String, u64, u64, String)>;

fn dash(v: &Value) -> String {
    match v {
        Value::Null => "-".into(),
        Value::String(s) if s.is_empty() => "-".into(),
        Value::String(s) => s.clone(),
        Value::Array(a) if a.is_empty() => "-".into(),
        Value::Array(a) => a.iter().map(dash).collect::<Vec<_>>().join(","),
        other => other.to_string(),
    }
}

/// gene summary as 1.1: exon depth and >= 20x share weighted by exon length (from the rounded exon values)
fn gene_summary(exons: &[Value], gene: &str) -> (usize, Option<f64>, Option<f64>, Vec<String>) {
    let rows: Vec<&Value> = exons.iter().filter(|e| e["gene"] == gene).collect();
    if rows.is_empty() {
        return (0, None, None, vec![]);
    }
    let total: f64 = rows.iter().map(|e| e["length"].as_f64().unwrap_or(0.0)).sum();
    let mean: f64 = rows.iter().map(|e| e["meanDepth"].as_f64().unwrap_or(0.0) * e["length"].as_f64().unwrap_or(0.0)).sum::<f64>() / total;
    let pct20: f64 = rows.iter().map(|e| e["pct20x"].as_f64().unwrap_or(0.0) * e["length"].as_f64().unwrap_or(0.0)).sum::<f64>() / total;
    let uncertain = rows.iter().filter(|e| e["status"] != "ok").map(|e| e["label"].as_str().unwrap_or("").to_string()).collect();
    (rows.len(), Some(exons::round(mean, 1)), Some(exons::round(pct20, 1)), uncertain)
}

/// One sentence per gene: what was found and how completely the gene was checked (English, German).
/// `result` is pathogenic / candidate / none once the variant analysis exists; until then None.
pub fn statement(gene: &str, result: Option<&str>, depth: Option<f64>, pct20: Option<f64>, uncertain: &[String]) -> Value {
    let cov_en = match (depth, pct20) {
        (Some(d), Some(p)) => format!("exons {d:.1}x, {p:.1} % of exon bases >= 20x"),
        _ => "exon coverage not determined".into(),
    };
    let cov_de = match (depth, pct20) {
        (Some(d), Some(p)) => format!("Exons {d:.1}x, {p:.1} % der Exonbasen >= 20x"),
        _ => "Exonabdeckung nicht bestimmt".into(),
    };
    let unc_en = if uncertain.is_empty() { "all exons checked".to_string() } else { format!("uncertain: {}", uncertain.join(", ")) };
    let unc_de = if uncertain.is_empty() { "alle Exons geprüft".to_string() } else { format!("unsicher: {}", uncertain.join(", ")) };
    let (en, de) = match result {
        Some("pathogenic") => (format!("{gene}: pathogenic variant found (confirmation in an accredited laboratory needed)"),
                               format!("{gene}: pathogene Variante gefunden (Bestätigung in einem akkreditierten Labor nötig)")),
        Some("candidate") => (format!("{gene}: finding that needs confirmation"), format!("{gene}: Befund, der bestätigt werden muss")),
        Some(_) => (format!("{gene}: no pathogenic variant found"), format!("{gene}: keine pathogene Variante gefunden")),
        None => (format!("{gene}: variants not analysed yet (coverage only)"), format!("{gene}: Varianten noch nicht ausgewertet (nur Abdeckung)")),
    };
    json!({"en": format!("{en}; {cov_en}; {unc_en}"), "de": format!("{de}; {cov_de}; {unc_de}")})
}

pub struct Run<'a> {
    pub input: &'a str,
    pub reference: Option<&'a str>,
    pub prefix: &'a str,
    pub exons: &'a Bed,
    pub regions: &'a Bed,
    pub calling_padding: u64,
    pub settings: &'a Settings,
}

pub fn run(r: &Run) -> io::Result<Value> {
    let contigs = reads::contigs(r.input, r.reference)?;
    let len = |c: &str| contigs.iter().find(|(n, _)| n == c).map(|&(_, l)| l);
    let table = exons::table(r.input, r.reference, r.exons, r.regions, &len, r.settings)?;
    let mut calling: Vec<(String, u64, u64)> =
        r.regions.iter().map(|(c, s, e, _)| (c.clone(), s.saturating_sub(r.calling_padding), e + r.calling_padding)).collect();
    calling.sort();
    let qc = readstats::region_stats(r.input, r.reference, &calling)?;
    let exon_rows: Vec<Value> = table["exonCoverage"].as_array().cloned().unwrap_or_default();
    let low: Vec<Value> = table["lowCoverageIntervals"].as_array().cloned().unwrap_or_default();
    let sample = std::path::Path::new(r.prefix).file_name().and_then(|n| n.to_str()).unwrap_or("sample").to_string();

    // PREFIX.exons.tsv
    let mut fe = File::create(format!("{}.exons.tsv", r.prefix))?;
    writeln!(fe, "gene\texon\tlabel\tchrom\tstart\tend\tlength\tmean_depth\tpct10x\tpct20x\tpct30x\tlow_mapq_fraction\tdosage_ratio\tstatus\treasons")?;
    let mut exons_json = vec![];
    for e in &exon_rows {
        let status = if e["status"] == "ok" { "ok" } else { "uncertain" };
        let reasons: Vec<String> = e["reasonsEN"].as_array().map(|a| a.iter().filter_map(|x| x.as_str().map(String::from)).collect()).unwrap_or_default();
        writeln!(
            fe,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            dash(&e["gene"]), dash(&e["exon"]), dash(&e["label"]), dash(&e["chrom"]), dash(&e["start"]), dash(&e["end"]),
            dash(&e["length"]), dash(&e["meanDepth"]), dash(&e["pct10x"]), dash(&e["pct20x"]), dash(&e["pct30x"]),
            dash(&e["lowMapqFraction"]), dash(&e["dosageRatio"]), status, if reasons.is_empty() { "-".into() } else { reasons.join("; ") }
        )?;
        exons_json.push(json!({
            "gene": e["gene"], "exon": e["exon"], "label": e["label"], "chrom": e["chrom"], "start": e["start"], "end": e["end"],
            "length": e["length"], "mean_depth": e["meanDepth"], "pct10x": e["pct10x"], "pct20x": e["pct20x"], "pct30x": e["pct30x"],
            "low_mapq_fraction": e["lowMapqFraction"], "dosage_ratio": e["dosageRatio"], "status": status,
            "reasons": {"en": e["reasonsEN"], "de": e["reasons"]},
        }));
    }

    // PREFIX.tsv: one line per gene
    let mut fg = File::create(format!("{}.tsv", r.prefix))?;
    writeln!(fg, "gene\ttranscript\tresult\tsmall_variants\tclinvar_plp\tloss_of_function\texon_copy_number\tcnv_candidates\tstructural_variants\tcoding_depth\tpct20x\tuncertain_exons\tlow_coverage_bp\tstatement")?;
    let mut genes_json = vec![];
    for (gene, nm, enst) in GENES {
        let (n, depth, pct20, uncertain) = gene_summary(&exon_rows, gene);
        let low_g: Vec<&Value> = low.iter().filter(|i| i["gene"] == gene).collect();
        let low_bp: u64 = low_g.iter().map(|i| i["length"].as_u64().unwrap_or(0)).sum();
        let st = statement(gene, None, depth, pct20, &uncertain);
        writeln!(
            fg,
            "{gene}\t{nm}\t-\t-\t-\t-\t-\t-\t-\t{}\t{}\t{}\t{low_bp}\t{}",
            depth.map(|d| format!("{d:.1}")).unwrap_or("-".into()),
            pct20.map(|d| format!("{d:.1}")).unwrap_or("-".into()),
            if uncertain.is_empty() { "-".into() } else { uncertain.join(",") },
            st["en"].as_str().unwrap_or("-"),
        )?;
        genes_json.push(json!({
            "gene": gene, "transcript": {"refseq": nm, "ensembl": enst}, "result": null, "statement": st,
            "small_variants": null, "clinvar_plp": null, "loss_of_function": null,
            "exon_copy_number": null, "cnv_candidates": null, "structural_variants": null,
            "coverage": {"exon_count": n, "coding_mean_depth": depth, "pct20x": pct20, "uncertain_exons": uncertain,
                         "low_coverage_bp": low_bp,
                         "low_coverage_intervals": low_g.iter().map(|i| json!({"chrom": i["chrom"], "start": i["start"], "end": i["end"], "length": i["length"], "min_depth": i["minDepth"]})).collect::<Vec<_>>()},
        }));
    }

    let st = r.settings;
    let rep = json!({
        "tool": "brca12-wgs", "version": env!("CARGO_PKG_VERSION"), "schema": SCHEMA,
        "sample": sample,
        "input": {"files": [r.input], "reference": r.reference},
        "analysis_time": std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0),
        "research_use_only": true,
        "analysed": ["coverage", "read_qc"],
        "not_yet_analysed": ["small_variants", "annotation", "exon_copy_number", "cnv", "structural_variants", "insertion_signals"],
        "thresholds": {"min_mapq": 20, "min_exon_depth": st.min_exon_depth, "min_exon_pct20x": st.min_exon_pct20,
                       "max_low_mapq_fraction": st.max_low_mapq, "low_coverage_depth": st.low_depth,
                       "calling_padding": r.calling_padding, "chromosome_window": st.window},
        "read_qc": {
            "reads_in_regions": qc["readsInRegions"], "mapped_pct": qc["mappedPct"], "properly_paired_pct": qc["properlyPairedPct"],
            "duplicate_pct": qc["duplicatePct"], "mapq0_pct": qc["mapq0Pct"], "median_insert_size": qc["medianInsertSize"],
            "mean_insert_size": qc["meanInsertSize"], "error_rate": qc["errorRate"], "average_read_length": qc["averageReadLength"],
            "average_base_quality": qc["averageBaseQuality"],
        },
        "chromosome_median": table["chromMedian"],
        "genes": genes_json,
        "exons": exons_json,
    });
    let mut fj = File::create(format!("{}.json", r.prefix))?;
    writeln!(fj, "{}", serde_json::to_string_pretty(&rep).map_err(io::Error::other)?)?;
    Ok(rep)
}
