//! Result files, built like the hla-wgs output: `PREFIX.tsv` (one line per gene), `PREFIX.exons.tsv`
//! (one line per exon with its certainty), `PREFIX.json` (schema genome.brca12.v2). The HTML page is
//! made from PREFIX.json at the end. Columns that a version does not determine yet are `-` in the
//! TSV files and `null` in the JSON (listed in `analysed`). Research use only.
use crate::annot::Tables;
use crate::caller;
use crate::exons::{self, Settings};
use crate::hgvs;
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
pub fn statement(gene: &str, result: Option<&str>, depth: Option<f64>, pct20: Option<f64>, uncertain: &[String], findings: &[(String, String)]) -> Value {
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
    let list_en = findings.iter().map(|f| f.0.clone()).collect::<Vec<_>>().join("; ");
    let list_de = findings.iter().map(|f| f.1.clone()).collect::<Vec<_>>().join("; ");
    let (en, de) = match result {
        Some("pathogenic") => (format!("{gene}: pathogenic variant found: {list_en} (confirmation in an accredited laboratory needed)"),
                               format!("{gene}: pathogene Variante gefunden: {list_de} (Bestätigung in einem akkreditierten Labor nötig)")),
        Some("candidate") => (format!("{gene}: finding that needs confirmation: {list_en}"), format!("{gene}: Befund, der bestätigt werden muss: {list_de}")),
        Some(_) => (format!("{gene}: no pathogenic variant found"), format!("{gene}: keine pathogene Variante gefunden")),
        None => (format!("{gene}: variants not analysed yet (coverage only)"), format!("{gene}: Varianten noch nicht ausgewertet (nur Abdeckung)")),
    };
    json!({"en": format!("{en}; {cov_en}; {unc_en}"), "de": format!("{de}; {cov_de}; {unc_de}")})
}

/// loss-of-function consequences (ENIGMA PVS1 territory)
const LOF: [&str; 5] = ["frameshift_variant", "stop_gained", "splice_donor_variant", "splice_acceptor_variant", "start_lost"];

/// one called allele with its annotation (a multi-allelic record gives two rows)
pub struct VariantRow {
    pub gene: String,
    pub in_gene: bool,
    pub chrom: String,
    pub pos: u64,
    pub r#ref: String,
    pub alt: String,
    pub genotype: String,
    pub zygosity: &'static str,
    pub filter: String,
    pub qual: f64,
    pub gq: u32,
    pub depth: u32,
    pub ref_reads: u32,
    pub alt_reads: u32,
    pub ann: hgvs::Annotation,
    pub clinvar: Option<crate::annot::ClinVar>,
    pub gnomad_present: bool,
    pub gnomad: Option<crate::annot::Gnomad>,
    pub class: &'static str,
}

impl VariantRow {
    pub fn lof(&self) -> bool {
        self.ann.consequence.split(',').any(|c| LOF.contains(&c))
    }
}

/// the tool's hint per variant (no ACMG classification): pathogenic = ClinVar pathogenic or likely
/// pathogenic (not conflicting); candidate = loss of function, not benign in ClinVar, gnomAD FAF95
/// below 0.1 %; benign = ClinVar benign/likely benign or FAF95 >= 0.1 % (ENIGMA BA1); otherwise other
fn classify(cv: Option<&crate::annot::ClinVar>, faf: Option<f64>, lof: bool) -> &'static str {
    if cv.is_some_and(|c| c.plp()) {
        "pathogenic"
    } else if cv.is_some_and(|c| c.benign()) || faf.is_some_and(|f| f >= 0.001) {
        "benign"
    } else if lof {
        "candidate"
    } else {
        "other"
    }
}

/// small variants of the gene regions (+/- padding) with HGVS, consequence, ClinVar and gnomAD
pub fn variants(input: &str, reference: &str, regions: &Bed, padding: u64, tables: &Tables) -> io::Result<Vec<VariantRow>> {
    let txs = hgvs::Tx::load(reference)?;
    let mut regs: Vec<(String, u64, u64)> = regions.iter().map(|(c, s, e, _)| (c.clone(), s.saturating_sub(padding), e + padding)).collect();
    regs.sort();
    let mut rows = vec![];
    for (c, s, e) in &regs {
        for call in caller::call_region(input, reference, c, *s, *e, &caller::DEFAULT)? {
            let pos = call.allele.pos + 1;
            let mut alts = vec![(call.allele.alt.clone(), 1u8)];
            if let Some(a2) = &call.alt2 {
                alts.push((a2.clone(), 2u8));
            }
            for (alt, idx) in alts {
                let ann = hgvs::annotate(&txs, &call.chrom, pos, &call.allele.r#ref, &alt).unwrap_or_default();
                let cv = tables.clinvar(&call.chrom, pos, &call.allele.r#ref, &alt).cloned();
                let (present, gn) = tables.gnomad(&call.chrom, pos, &call.allele.r#ref, &alt);
                let gn = gn.cloned();
                let in_gene = regions.iter().any(|(rc, rs, re, _)| rc == &call.chrom && pos > *rs && pos <= *re);
                let lof = ann.consequence.split(',').any(|c| LOF.contains(&c));
                let class = classify(cv.as_ref(), gn.as_ref().and_then(|g| g.faf95), lof);
                let zyg = if call.gt.0 == call.gt.1 { "homozygous" } else if call.alt2.is_some() { "heterozygous (two different alleles)" } else { "heterozygous" };
                rows.push(VariantRow {
                    gene: if ann.gene.is_empty() { "-".into() } else { ann.gene.clone() },
                    in_gene,
                    chrom: call.chrom.clone(),
                    pos,
                    r#ref: String::from_utf8_lossy(&call.allele.r#ref).into_owned(),
                    alt: String::from_utf8_lossy(&alt).into_owned(),
                    genotype: format!("{}/{}", call.gt.0, call.gt.1),
                    zygosity: zyg,
                    filter: call.filter.clone(),
                    qual: (call.qual * 100.0).round() / 100.0,
                    gq: call.gq,
                    depth: call.dp,
                    ref_reads: call.ad.0,
                    alt_reads: call.ad.1,
                    ann,
                    clinvar: cv,
                    gnomad_present: present,
                    gnomad: gn,
                    class,
                });
                let _ = idx;
            }
        }
    }
    rows.sort_by(|a, b| (a.chrom.as_str(), a.pos, a.alt.as_str()).cmp(&(b.chrom.as_str(), b.pos, b.alt.as_str())));
    Ok(rows)
}

const VARIANT_COLUMNS: &str = "gene\tchrom\tpos\tref\talt\tgenotype\tzygosity\tfilter\tqual\tgq\tdepth\tref_reads\talt_reads\tallele_fraction\thgvs_c\thgvs_p\tconsequence\texon\tintron\tclinvar_id\tclinvar_significance\tclinvar_review\tclinvar_stars\tgnomad_present\tgnomad_faf95\tgnomad_population\tgnomad_dataset\tclass\tin_gene_region";

fn af(v: &VariantRow) -> Option<f64> {
    let t = v.ref_reads + v.alt_reads;
    if t == 0 { None } else { Some(exons::round(v.alt_reads as f64 / t as f64, 3)) }
}

fn write_variants(path: &str, rows: &[VariantRow]) -> io::Result<()> {
    let mut f = File::create(path)?;
    writeln!(f, "{VARIANT_COLUMNS}")?;
    let d = |o: &Option<String>| o.clone().unwrap_or("-".into());
    for v in rows {
        let cv = v.clinvar.as_ref();
        let gn = v.gnomad.as_ref();
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            v.gene, v.chrom, v.pos, v.r#ref, v.alt, v.genotype, v.zygosity, v.filter, v.qual, v.gq, v.depth, v.ref_reads, v.alt_reads,
            af(v).map(|x| x.to_string()).unwrap_or("-".into()),
            d(&v.ann.hgvs_c), d(&v.ann.hgvs_p), if v.ann.consequence.is_empty() { "-".into() } else { v.ann.consequence.clone() },
            d(&v.ann.exon), d(&v.ann.intron),
            cv.map(|c| c.id.clone()).unwrap_or("-".into()), cv.map(|c| c.significance.clone()).unwrap_or("-".into()),
            cv.map(|c| c.review.clone()).unwrap_or("-".into()), cv.map(|c| c.stars().to_string()).unwrap_or("-".into()),
            if v.gnomad_present { "yes" } else { "no" },
            gn.and_then(|g| g.faf95).map(|x| x.to_string()).unwrap_or("-".into()),
            gn.map(|g| g.population.clone()).unwrap_or("-".into()), gn.map(|g| g.dataset.clone()).unwrap_or("-".into()),
            v.class, if v.in_gene { "yes" } else { "no" }
        )?;
    }
    Ok(())
}

fn variant_json(v: &VariantRow) -> Value {
    json!({
        "gene": v.gene, "chrom": v.chrom, "pos": v.pos, "ref": v.r#ref, "alt": v.alt, "genotype": v.genotype,
        "zygosity": v.zygosity, "filter": v.filter, "qual": v.qual, "gq": v.gq, "depth": v.depth,
        "ref_reads": v.ref_reads, "alt_reads": v.alt_reads, "allele_fraction": af(v),
        "hgvs_c": v.ann.hgvs_c, "hgvs_p": v.ann.hgvs_p, "consequence": v.ann.consequence.split(',').filter(|c| !c.is_empty()).collect::<Vec<_>>(),
        "exon": v.ann.exon, "intron": v.ann.intron,
        "clinvar": v.clinvar.as_ref().map(|c| json!({"id": c.id, "significance": c.significance, "review": c.review, "stars": c.stars()})),
        "gnomad": {"present": v.gnomad_present, "faf95": v.gnomad.as_ref().and_then(|g| g.faf95),
                   "population": v.gnomad.as_ref().map(|g| g.population.clone()), "dataset": v.gnomad.as_ref().map(|g| g.dataset.clone()),
                   "non_cancer": v.gnomad.as_ref().map(|g| g.non_cancer)},
        "class": v.class, "in_gene_region": v.in_gene,
    })
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
    // small variants with annotation (needs the reference FASTA)
    let reference = r.reference.ok_or_else(|| io::Error::other("--reference REF.fa is needed (variant calling and HGVS)"))?;
    let tables = Tables::load();
    let vrows = variants(r.input, reference, r.regions, r.calling_padding, &tables)?;
    write_variants(&format!("{}.variants.tsv", r.prefix), &vrows)?;

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
        // variant result of the gene: PASS calls inside the gene region
        let mine: Vec<&VariantRow> = vrows.iter().filter(|v| v.gene == gene && v.in_gene).collect();
        let pass: Vec<&&VariantRow> = mine.iter().filter(|v| v.filter == "PASS").collect();
        let plp: Vec<&&VariantRow> = pass.iter().filter(|v| v.class == "pathogenic").copied().collect();
        let lof: Vec<&&VariantRow> = pass.iter().filter(|v| v.lof() && v.class != "benign").copied().collect();
        let cand: Vec<&&VariantRow> = mine.iter().filter(|v| v.class == "candidate" || (v.class == "pathogenic" && v.filter != "PASS")).collect();
        let result = if !plp.is_empty() { "pathogenic" } else if !cand.is_empty() { "candidate" } else { "none" };
        let shown: Vec<&&VariantRow> = if result == "pathogenic" { plp.clone() } else { cand.clone() };
        let findings: Vec<(String, String)> = shown.iter().map(|v| {
            let name = format!("{}{}", v.ann.hgvs_c.clone().unwrap_or(format!("{}:{}{}>{}", v.chrom, v.pos, v.r#ref, v.alt)),
                               v.ann.hgvs_p.as_ref().map(|p| format!(" p.({})", p.split_once(":p.").map(|x| x.1).unwrap_or(p))).unwrap_or_default());
            let (zen, zde) = match v.zygosity { "homozygous" => ("homozygous", "homozygot"), _ => ("heterozygous", "heterozygot") };
            let fl = if v.filter == "PASS" { String::new() } else { format!(", filter {}", v.filter) };
            (format!("{name}, {zen}{fl}"), format!("{name}, {zde}{fl}"))
        }).collect();
        let st = statement(gene, Some(result), depth, pct20, &uncertain, &findings);
        writeln!(
            fg,
            "{gene}\t{nm}\t{result}\t{}\t{}\t{}\t-\t-\t-\t{}\t{}\t{}\t{low_bp}\t{}",
            pass.len(), plp.len(), lof.len(),
            depth.map(|d| format!("{d:.1}")).unwrap_or("-".into()),
            pct20.map(|d| format!("{d:.1}")).unwrap_or("-".into()),
            if uncertain.is_empty() { "-".into() } else { uncertain.join(",") },
            st["en"].as_str().unwrap_or("-"),
        )?;
        genes_json.push(json!({
            "gene": gene, "transcript": {"refseq": nm, "ensembl": enst}, "result": result, "statement": st,
            "small_variants": pass.len(), "clinvar_plp": plp.len(), "loss_of_function": lof.len(),
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
        "analysed": ["coverage", "read_qc", "small_variants", "annotation"],
        "not_yet_analysed": ["exon_copy_number", "cnv", "structural_variants", "insertion_signals"],
        "annotation_sources": {"transcripts": "MANE Select v1 (BRCA1 NM_007294.4 / ENST00000357654, BRCA2 NM_000059.4 / ENST00000380152)",
                               "clinvar": tables.clinvar_version, "gnomad": tables.gnomad_version},
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
        "variants": vrows.iter().map(variant_json).collect::<Vec<_>>(),
        "exons": exons_json,
    });
    let mut fj = File::create(format!("{}.json", r.prefix))?;
    writeln!(fj, "{}", serde_json::to_string_pretty(&rep).map_err(io::Error::other)?)?;
    Ok(rep)
}
