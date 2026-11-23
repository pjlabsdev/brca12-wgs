//! Exon copy number against a reference cohort, a one-to-one port of brca12-wgs 1.1
//! (bin/brca12lib/exon_cnv.py). Research use only.
//!
//! Per exon: mean depth over the exon +/- `pad` bases (as `samtools depth -a -Q 20`), divided by
//! the median depth of single-copy control windows of the same sample, log2, minus the exon's
//! centre in the panel of 80 1000 Genomes 30x genomes (assets/exon_cnv_panel_1kg30x.json). A hidden
//! Markov model over the exons of each gene (copy number 0-4, Gaussian emissions with the exon's
//! panel spread, allele balance at heterozygous SNVs) gives segments; quality = -10 log10 of the
//! posterior probability that the segment's exons have two copies; >= 10 a call, 5-9 a hint.
use crate::coverage::{self, Filter, SAMTOOLS_DEPTH};
use serde_json::{Value, json};
use std::collections::HashMap;
use std::io;

const PANEL: &str = include_str!("../assets/exon_cnv_panel_1kg30x.json");
const STATES: [usize; 5] = [0, 1, 2, 3, 4];
const HINT_QUALITY: i64 = 5;

fn shift(c: usize) -> f64 {
    [-3.0, -1.0, 0.0, 1.5f64.log2(), 1.0][c]
}

fn prior(c: usize) -> f64 {
    [1e-5, 5e-4, 1.0 - 1e-3, 5e-4 - 2e-5, 1e-5][c]
}

pub struct Exon {
    pub name: String,
    pub chrom: String,
    pub start: u64, // 0-based
    pub end: u64,
    pub gene: String,
    pub label: String,
    pub median: f64,
    pub sd: f64,
}

pub struct Panel {
    pub description: String,
    pub size: u64,
    pub pad: u64,
    pub p_change: f64,
    pub min_quality: i64,
    pub exons: Vec<Exon>,
    pub controls: Vec<(String, String, u64, u64)>, // name, chrom, start, end
}

pub fn panel() -> Panel {
    let v: Value = serde_json::from_str(PANEL).expect("bundled exon CNV panel");
    let exons = v["exons"]
        .as_object()
        .unwrap()
        .iter()
        .map(|(n, e)| Exon {
            name: n.clone(),
            chrom: e["chrom"].as_str().unwrap().into(),
            start: e["start"].as_u64().unwrap(),
            end: e["end"].as_u64().unwrap(),
            gene: e["gene"].as_str().unwrap().into(),
            label: e["label"].as_str().unwrap_or(n).into(),
            median: e["median"].as_f64().unwrap(),
            sd: e["sd"].as_f64().unwrap(),
        })
        .collect();
    let controls = v["controls"]
        .as_object()
        .unwrap()
        .iter()
        .map(|(n, c)| (n.clone(), c["chrom"].as_str().unwrap().into(), c["start"].as_u64().unwrap(), c["end"].as_u64().unwrap()))
        .collect();
    Panel {
        description: v["description"].as_str().unwrap_or("").into(),
        size: v["size"].as_u64().unwrap_or(0),
        pad: v["pad"].as_u64().unwrap(),
        p_change: v["pChange"].as_f64().unwrap(),
        min_quality: v["minQuality"].as_i64().unwrap(),
        exons,
        controls,
    }
}

/// mean depth per interval (exons +/- pad, control windows), samtools depth -a -Q 20
pub fn depths(path: &str, reference: Option<&str>, p: &Panel) -> io::Result<HashMap<String, Option<f64>>> {
    let flt = Filter { min_mapq: 20, ..SAMTOOLS_DEPTH };
    let mut out = HashMap::new();
    let mut rows: Vec<(String, String, u64, u64)> =
        p.exons.iter().map(|e| (e.name.clone(), e.chrom.clone(), e.start.saturating_sub(p.pad), e.end + p.pad)).collect();
    rows.extend(p.controls.iter().cloned());
    // one depth pass per contiguous block of a chromosome (exons of one gene lie close together)
    rows.sort_by(|a, b| (a.1.as_str(), a.2).cmp(&(b.1.as_str(), b.2)));
    let mut i = 0;
    while i < rows.len() {
        let mut j = i;
        let mut end = rows[i].3;
        while j + 1 < rows.len() && rows[j + 1].1 == rows[i].1 && rows[j + 1].2 <= end + 50_000 {
            j += 1;
            end = end.max(rows[j].3);
        }
        let (chrom, start) = (rows[i].1.clone(), rows[i].2);
        let d = coverage::depth(path, reference, &chrom, start, end, flt)?;
        for r in &rows[i..=j] {
            let s = (r.2 - start) as usize;
            let e = ((r.3 - start) as usize).min(d.len());
            let mean = if e > s { Some(d[s..e].iter().map(|&x| x as f64).sum::<f64>() / (e - s) as f64) } else { None };
            out.insert(r.0.clone(), mean);
        }
        i = j + 1;
    }
    Ok(out)
}

fn logpdf(x: f64, mu: f64, sd: f64) -> f64 {
    -0.5 * ((x - mu) / sd).powi(2) - (sd * (2.0 * std::f64::consts::PI).sqrt()).ln()
}

fn lse(v: &[f64]) -> f64 {
    let m = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    m + v.iter().map(|x| (x - m).exp()).sum::<f64>().ln()
}

fn ab_mix(c: usize) -> &'static [f64] {
    match c {
        0 | 1 => &[0.02, 0.98],
        2 => &[0.5],
        3 => &[1.0 / 3.0, 2.0 / 3.0],
        _ => &[0.25, 0.5, 0.75],
    }
}

fn allele_loglik(counts: &[(u32, u32)]) -> [f64; 5] {
    let mut out = [0.0; 5];
    for &(r, a) in counts {
        if r + a == 0 {
            continue;
        }
        for c in STATES {
            let fr = ab_mix(c);
            let v: Vec<f64> = fr.iter().map(|&f| a as f64 * f.ln() + r as f64 * (1.0 - f).ln() - (fr.len() as f64).ln()).collect();
            out[c] += lse(&v);
        }
    }
    out
}

/// forward-backward posteriors and the Viterbi path over the exons of one gene
fn hmm(ys: &[f64], sds: &[f64], p_change: f64, counts: Option<&[Vec<(u32, u32)>]>) -> (Vec<[f64; 5]>, Vec<usize>) {
    let n = ys.len();
    let mut emis: Vec<[f64; 5]> = (0..n)
        .map(|i| {
            let mut e = [0.0; 5];
            for c in STATES {
                e[c] = logpdf(ys[i], shift(c), if c != 0 { sds[i] } else { 1.0 });
            }
            e
        })
        .collect();
    if let Some(cs) = counts {
        for i in 0..n {
            if !cs[i].is_empty() {
                let ab = allele_loglik(&cs[i]);
                for c in STATES {
                    emis[i][c] += ab[c];
                }
            }
        }
    }
    let tr = |a: usize, b: usize| if a == b { (1.0 - p_change).ln() } else { (p_change * prior(b) / (1.0 - prior(a))).ln() };
    let start: [f64; 5] = std::array::from_fn(|c| prior(c).ln() + emis[0][c]);
    let mut fw = vec![start];
    for i in 1..n {
        let prev = fw[i - 1];
        fw.push(std::array::from_fn(|b| emis[i][b] + lse(&STATES.map(|a| prev[a] + tr(a, b)))));
    }
    let mut bw = vec![[0.0; 5]; n];
    for i in (0..n.saturating_sub(1)).rev() {
        let next = bw[i + 1];
        bw[i] = std::array::from_fn(|a| lse(&STATES.map(|b| tr(a, b) + emis[i + 1][b] + next[b])));
    }
    let total = lse(&fw[n - 1]);
    let post: Vec<[f64; 5]> = (0..n).map(|i| std::array::from_fn(|c| (fw[i][c] + bw[i][c] - total).exp())).collect();
    // Viterbi (ties: the first state, as Python's max)
    let argmax = |v: &dyn Fn(usize) -> f64| {
        let mut best = 0;
        for s in 1..5 {
            if v(s) > v(best) {
                best = s;
            }
        }
        best
    };
    let mut v = vec![start];
    let mut back: Vec<[usize; 5]> = vec![];
    for i in 1..n {
        let prev = v[i - 1];
        let mut row = [0.0; 5];
        let mut br = [0usize; 5];
        for b in STATES {
            let a = argmax(&|s| prev[s] + tr(s, b));
            row[b] = prev[a] + tr(a, b) + emis[i][b];
            br[b] = a;
        }
        v.push(row);
        back.push(br);
    }
    let last = v[n - 1];
    let mut path = vec![argmax(&|c| last[c])];
    for br in back.iter().rev() {
        path.push(br[*path.last().unwrap()]);
    }
    path.reverse();
    (post, path)
}

/// Python round(x, d)
fn round(x: f64, d: usize) -> f64 {
    crate::exons::round(x, d)
}

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n % 2 == 1 { v[n / 2] } else { (v[n / 2 - 1] + v[n / 2]) / 2.0 }
}

/// exon table and copy-number segments; `alleles`: per exon name (ref, alt) reads at heterozygous SNVs
pub fn call(p: &Panel, depths: &HashMap<String, Option<f64>>, alleles: Option<&HashMap<String, Vec<(u32, u32)>>>) -> Value {
    let mut controls: Vec<f64> = p.controls.iter().filter_map(|c| depths.get(&c.0).copied().flatten()).filter(|&d| d != 0.0).collect();
    let mut res = serde_json::Map::new();
    res.insert("enabled".into(), json!(true));
    res.insert("panel".into(), json!(p.description));
    res.insert("panel_size".into(), json!(p.size));
    res.insert("pad".into(), json!(p.pad));
    res.insert("min_quality".into(), json!(p.min_quality));
    res.insert("hint_quality".into(), json!(HINT_QUALITY));
    res.insert("allele_balance".into(), json!(alleles.is_some_and(|a| !a.is_empty())));
    if controls.len() < 3.max(p.controls.len() / 2) {
        res.insert("enabled".into(), json!(false));
        res.insert("reason".into(), json!("too few control windows with depth"));
        res.insert("exons".into(), json!([]));
        res.insert("calls".into(), json!([]));
        res.insert("hints".into(), json!([]));
        return Value::Object(res);
    }
    let norm = median(&mut controls);
    res.insert("normaliser_depth".into(), json!(round(norm, 2)));
    let alleles = alleles.filter(|a| !a.is_empty());
    let mut genes: Vec<String> = vec![];
    for e in &p.exons {
        if !genes.contains(&e.gene) {
            genes.push(e.gene.clone());
        }
    }
    let (mut exons_out, mut calls, mut hints) = (vec![], vec![], vec![]);
    for gene in genes {
        let mut rows: Vec<&Exon> = p.exons.iter().filter(|e| e.gene == gene && depths.get(&e.name).copied().flatten().is_some()).collect();
        rows.sort_by_key(|e| e.start);
        if rows.is_empty() {
            continue;
        }
        let ys: Vec<f64> = rows.iter().map(|e| (depths[&e.name].unwrap().max(0.05) / norm).log2() - e.median).collect();
        let sds: Vec<f64> = rows.iter().map(|e| e.sd).collect();
        let counts: Option<Vec<Vec<(u32, u32)>>> = alleles.map(|a| rows.iter().map(|e| a.get(&e.name).cloned().unwrap_or_default()).collect());
        let (post, path) = hmm(&ys, &sds, p.p_change, counts.as_deref());
        for (k, e) in rows.iter().enumerate() {
            exons_out.push(json!({
                "gene": gene, "exon": e.name, "label": e.label, "log2_deviation": round(ys[k], 3),
                "het_sites": alleles.and_then(|a| a.get(&e.name)).map(|v| v.len()).unwrap_or(0),
                "z": round(ys[k] / sds[k], 2), "copy_number": path[k], "p_two_copies": round(post[k][2], 4),
            }));
        }
        let mut i = 0;
        while i < rows.len() {
            if path[i] == 2 {
                i += 1;
                continue;
            }
            let mut j = i;
            while j + 1 < rows.len() && path[j + 1] == path[i] {
                j += 1;
            }
            let p2 = (i..=j).map(|k| post[k][2]).fold(f64::NEG_INFINITY, f64::max).max(1e-300);
            let quality = 99.min((-10.0 * p2.log10()) as i64);
            if quality >= HINT_QUALITY {
                let (first, last) = (rows[i], rows[j]);
                let seg = json!({
                    "gene": gene, "copy_number": path[i], "type": if path[i] < 2 { "deletion" } else { "duplication" },
                    "chrom": first.chrom, "start": first.start.min(last.start) + 1, "end": first.end.max(last.end),
                    "exons": rows[i..=j].iter().map(|e| e.name.clone()).collect::<Vec<_>>(),
                    "affected_exons": rows[i..=j].iter().map(|e| e.label.clone()).collect::<Vec<_>>(),
                    "quality": quality,
                    "mean_log2_deviation": round(ys[i..=j].iter().sum::<f64>() / (j - i + 1) as f64, 3),
                });
                if quality >= p.min_quality { calls.push(seg) } else { hints.push(seg) }
            }
            i = j + 1;
        }
    }
    res.insert("exons".into(), json!(exons_out));
    res.insert("calls".into(), json!(calls));
    res.insert("hints".into(), json!(hints));
    Value::Object(res)
}

/// (ref, alt) reads at heterozygous biallelic SNVs inside exon +/- pad; same site filter as 1.1
/// (>= 2 reads per allele, minor allele >= 0.15)
pub fn allele_counts(p: &Panel, snvs: &[(String, u64, u32, u32)]) -> HashMap<String, Vec<(u32, u32)>> {
    let mut out: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
    for (chrom, pos, r, a) in snvs {
        let (r, a) = (*r, *a);
        if r.min(a) < 2 || (r.min(a) as f64) / ((r + a) as f64) < 0.15 {
            continue;
        }
        for e in &p.exons {
            if &e.chrom == chrom && e.start.saturating_sub(p.pad) < *pos && *pos <= e.end + p.pad {
                out.entry(e.name.clone()).or_default().push((r, a));
            }
        }
    }
    out
}
