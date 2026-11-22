//! Exon table, chromosome baseline and low-coverage intervals, computed like brca12-wgs 1.1
//! (mosdepth for the exons and the chromosome windows, samtools depth -a for the per-base gene
//! depth) so the results can be compared one to one. Research use only.
use crate::coverage::{self, Filter, MOSDEPTH, SAMTOOLS_DEPTH};
use crate::reads::{self, Read};
use serde_json::{Value, json};
use std::io;

pub struct Settings {
    pub min_exon_depth: f64,  // 20
    pub min_exon_pct20: f64,  // 95
    pub max_low_mapq: f64,    // 0.10
    pub low_depth: f64,       // 10
    pub window: u64,          // 100,000 bp chromosome windows
}

pub const DEFAULT: Settings = Settings { min_exon_depth: 20.0, min_exon_pct20: 95.0, max_low_mapq: 0.10, low_depth: 10.0, window: 100_000 };

/// Exon label as in 1.1: BRCA1 exons 4-23 carry their legacy number (+1) in brackets.
pub fn exon_label(name: &str) -> (String, u32, String) {
    let (gene, num) = name.split_once("_ex").unwrap_or((name, "0"));
    let n: u32 = num.parse().unwrap_or(0);
    let label = if gene == "BRCA1" && n >= 4 { format!("Exon {n} (Legacy {})", n + 1) } else { format!("Exon {n}") };
    (gene.to_string(), n, label)
}

/// Python round(x, d) for the values used here (correctly rounded decimal of the binary value).
pub fn round(x: f64, d: usize) -> f64 {
    format!("{x:.d$}").parse().unwrap()
}

/// mosdepth --fast-mode --by WINDOW over one contig: mean depth per window from read spans
/// (alignment start to end, no CIGAR detail, no mate correction), printed with two decimals.
pub fn window_means(path: &str, reference: Option<&str>, chrom: &str, len: u64, window: u64, min_mapq: u8) -> io::Result<Vec<f64>> {
    let nwin = len.div_ceil(window) as usize;
    let mut sums = vec![0u64; nwin];
    reads::fetch(path, reference, &format!("{chrom}:1-{len}"), &mut |r: &Read| {
        if r.flags & 0x704 != 0 || r.mapq < min_mapq {
            return;
        }
        let (mut s, e) = (r.start, r.end.min(len));
        while s < e {
            let w = (s / window) as usize;
            let we = ((w as u64 + 1) * window).min(e);
            sums[w] += we - s;
            s = we;
        }
    })?;
    Ok(sums
        .iter()
        .enumerate()
        .map(|(i, &sum)| {
            let ws = i as u64 * window;
            let wl = (len.min(ws + window) - ws) as f64;
            round(sum as f64 / wl, 2)
        })
        .collect())
}

/// statistics.median of the window means above zero
pub fn median_positive(v: &[f64]) -> Option<f64> {
    let mut p: Vec<f64> = v.iter().copied().filter(|&x| x > 0.0).collect();
    if p.is_empty() {
        return None;
    }
    p.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = p.len();
    Some(if n % 2 == 1 { p[n / 2] } else { (p[n / 2 - 1] + p[n / 2]) / 2.0 })
}

type Bed = Vec<(String, u64, u64, String)>;

/// exonCoverage (1.1 schema), the chromosome medians and lowCoverageIntervals.
pub fn table(path: &str, reference: Option<&str>, exons: &Bed, regions: &Bed, contig_len: &dyn Fn(&str) -> Option<u64>, st: &Settings) -> io::Result<Value> {
    let q20 = Filter { min_mapq: 20, ..MOSDEPTH };
    let q0 = Filter { min_mapq: 0, ..MOSDEPTH };
    let mut chroms: Vec<String> = exons.iter().map(|b| b.0.clone()).collect();
    chroms.sort();
    chroms.dedup();
    let mut medians = serde_json::Map::new();
    let mut rows = vec![];
    for chrom in &chroms {
        let len = contig_len(chrom).ok_or_else(|| io::Error::other(format!("{chrom} is not in the alignment header")))?;
        let median = median_positive(&window_means(path, reference, chrom, len, st.window, 20)?).map(|m| round(m, 3));
        medians.insert(chrom.clone(), json!(median));
        let mine: Vec<_> = exons.iter().filter(|b| &b.0 == chrom).collect();
        let s = mine.iter().map(|b| b.1).min().unwrap();
        let e = mine.iter().map(|b| b.2).max().unwrap();
        let d20 = coverage::depth(path, reference, chrom, s, e, q20)?;
        let d0 = coverage::depth(path, reference, chrom, s, e, q0)?;
        for b in mine {
            let a = coverage::interval_stats(chrom, b.1, b.2, &b.3, &d20, s, &[10, 20, 30]);
            let z = coverage::interval_stats(chrom, b.1, b.2, &b.3, &d0, s, &[]);
            // 1.1 reads mosdepth's two-decimal means
            let mean = round(a.mean, 2);
            let all_mean = round(z.mean, 2);
            let length = (b.2 - b.1) as f64;
            let pct: Vec<f64> = a.at_least.iter().map(|&(_, n)| 100.0 * n as f64 / length).collect();
            let low_mapq = if all_mean > 0.0 { Some(1.0 - mean / all_mean) } else { None };
            let dosage = median.filter(|&m| m > 0.0).map(|m| mean / m);
            let mut reasons: Vec<String> = vec![];
            if mean < st.min_exon_depth {
                reasons.push(format!("mittlere Tiefe {mean:.1}x < {}x", st.min_exon_depth));
            }
            if pct[1] < st.min_exon_pct20 {
                reasons.push(format!("nur {:.1} % der Basen >= 20x", pct[1]));
            }
            if let Some(f) = low_mapq.filter(|&f| f > st.max_low_mapq) {
                reasons.push(format!("{:.0} % der Reads mit MAPQ < 20 (Pseudogen/Repeat)", 100.0 * f));
            }
            let (gene, number, label) = exon_label(&b.3);
            rows.push(json!({
                "gene": gene, "exon": number, "label": label, "chrom": chrom,
                "start": b.1 + 1, "end": b.2, "length": b.2 - b.1,
                "meanDepth": round(mean, 1),
                "pct10x": round(pct[0], 1), "pct20x": round(pct[1], 1), "pct30x": round(pct[2], 1),
                "lowMapqFraction": low_mapq.map(|f| round(f, 3)),
                "dosageRatio": dosage.map(|x| round(x, 2)),
                "status": if reasons.is_empty() { "ok" } else { "unsicher" },
                "reasons": reasons,
            }));
        }
    }
    rows.sort_by_key(|r| (r["gene"].as_str().unwrap_or("").to_string(), r["exon"].as_u64().unwrap_or(0)));
    // low-coverage intervals: samtools depth -a -Q 20 over the gene regions, depth < low_depth
    let mut low = vec![];
    let sd = Filter { min_mapq: 20, ..SAMTOOLS_DEPTH };
    let mut regs: Vec<_> = regions.iter().collect();
    regs.sort_by(|a, b| (a.0.as_str(), a.1).cmp(&(b.0.as_str(), b.1)));
    for (chrom, s, e, gene) in regs {
        let d = coverage::depth(path, reference, chrom, *s, *e, sd)?;
        let mut cur: Option<(u64, u64, u32)> = None;
        for (i, &v) in d.iter().enumerate() {
            let p = s + i as u64;
            if (v as f64) < st.low_depth {
                cur = Some(match cur {
                    Some((cs, _, m)) => (cs, p + 1, m.min(v)),
                    None => (p, p + 1, v),
                });
            } else if let Some((cs, ce, m)) = cur.take() {
                low.push(json!({"gene": gene, "chrom": chrom, "start": cs + 1, "end": ce, "minDepth": m as f64, "length": ce - cs}));
            }
        }
        if let Some((cs, ce, m)) = cur {
            low.push(json!({"gene": gene, "chrom": chrom, "start": cs + 1, "end": ce, "minDepth": m as f64, "length": ce - cs}));
        }
    }
    Ok(json!({"exonCoverage": rows, "chromMedian": medians, "lowCoverageIntervals": low}))
}
