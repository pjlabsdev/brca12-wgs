//! Structural-variant and insertion signatures inside BRCA1/BRCA2, a one-to-one port of brca12-wgs 1.1
//! (bin/brca12_insertion_signals.py, run with --gene-fasta and --read-pair-svs). Research use only.
//!
//! - Alu insertions: clusters of soft-clipped read ends from both sides of a site whose clipped
//!   sequence carries Alu consensus k-mers.
//! - Extra tandem repeat copies: left and right clip junctions that point to the same unit length.
//! - Known alleles from long-read assemblies (the 126 bp copy upstream of BRCA1 exon 1, the common
//!   intronic BRCA2 Alu, the 980 bp BRCA1 intron deletion, and the pathogenic founder alleles BRCA2
//!   c.156_157insAlu, BRCA1 exon 13 duplication and 3,835 bp deletion) by allele-aware realignment.
//! - Tandem duplications from outward-facing read pairs, deletions and inversions from discordant
//!   read pairs and split reads (SA tag).
//! Field names in the result are those of 1.1, so that both versions can be compared one to one.
use crate::reads::{self, Op, Read};
use indexmap::IndexMap;
use serde_json::{Value, json};
use std::collections::{HashMap, HashSet};
use std::io;

const ALU12: [&str; 10] = [
    "GGCCGGGCGCGG", "GGCTCACGCCTG", "CTGTAATCCCAG", "GCACTTTGGGAG", "GAGGCGGAGGTT", "GCGAGACTCCGT", "CAGCCTGGGCGA",
    "GGGAGGCTGAGG", "TACTCGGGAGGC", "GGGCGTGGTGGC",
];

fn revcomp(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            'N' => 'N',
            x => x,
        })
        .collect()
}

fn alu_kmers() -> Vec<String> {
    let mut v: Vec<String> = ALU12.iter().map(|s| s.to_string()).collect();
    v.extend(ALU12.iter().map(|s| revcomp(s)));
    v
}

/// chrom, marker position, reference base, base of the extra copy, unit start, unit end, unit length, label
type Copy = (&'static str, i64, char, char, i64, i64, i64, &'static str);
const KNOWN_COPIES: [Copy; 1] =
    [("chr17", 43128365, 'T', 'C', 43128247, 43128372, 126, "zusätzliche 126-bp-Kopie vor BRCA1 Exon 1 (häufiges Allel)")];

enum Ins {
    Seq(String),
    Dup(i64, i64),
}

struct Known {
    id: &'static str,
    chrom: &'static str,
    anchor: i64,
    ins: Ins,
    tandem: bool,
    label: &'static str,
    deleted: i64,
}

fn known_alleles() -> Vec<Known> {
    let alu_brca2 = format!(
        "{}{}{}{}",
        "AAAAGTTTTGGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATCGAGACCATCCCGGCTAAAAC",
        "GGTGAAACCCCGTCTCTACTAAAAATACAAAAAAATTAGCCGGGCGTAGTGGCGGGCGCCTGTAGTCCCAGCTACTTGGGAGGCTGAGGCAGGAGAATGGTGTGAAC",
        "CCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCCCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTC",
        "A".repeat(27)
    );
    let c156 = format!(
        "{}{}{}{}{}",
        "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATCGAGACCATCCCGGCTAAAACGG",
        "TGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCCGGGCGTAGTGGCGGGCGCCTGTAGTCCCAGCTACTTGGGAGGCTGAGGCAGGAGAATGGCGT",
        "GAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCCCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAA",
        "AAAAAAAAAAA",
        "AAGAATCTGAACAT"
    );
    vec![
        Known {
            id: "tandem126",
            chrom: "chr17",
            anchor: 43128264,
            ins: Ins::Seq(
                "GATGGGAGGAACGAACAAACTCCAGACATGCCGCCTTAAGAGCTGTAACACTCACCGCAAAGGTCTGCAGCTTCACTCCTGAAGCCAGCGAGACCACGAACCCACCAGGAGGAATGAACAACTCCA"
                    .into(),
            ),
            tandem: true,
            label: "zusätzliche 126-bp-Kopie vor BRCA1 Exon 1 (häufiges Allel)",
            deleted: 0,
        },
        Known {
            id: "alu_brca2",
            chrom: "chr13",
            anchor: 32348295,
            ins: Ins::Seq(alu_brca2),
            tandem: false,
            label: "bekannte intronische Alu-Insertion in BRCA2 (häufiges Allel)",
            deleted: 0,
        },
        Known {
            id: "del980_brca1",
            chrom: "chr17",
            anchor: 43071727,
            ins: Ins::Seq(String::new()),
            tandem: false,
            label: "bekannte intronische 980-bp-Deletion in BRCA1 (seltenes Allel)",
            deleted: 980,
        },
        Known {
            id: "alu_brca2_c156",
            chrom: "chr13",
            anchor: 32319165,
            ins: Ins::Seq(c156),
            tandem: false,
            label: "BRCA2 c.156_157insAlu (Gründermutation, Exon 3)",
            deleted: 0,
        },
        Known {
            id: "dup13_brca1",
            chrom: "chr17",
            anchor: 43084362,
            ins: Ins::Dup(43078282, 43084362),
            tandem: false,
            label: "BRCA1 Duplikation Exon 13 (alte Zählung; Gründermutation)",
            deleted: 0,
        },
        Known {
            id: "del3835_brca1",
            chrom: "chr17",
            anchor: 43080383,
            ins: Ins::Seq(String::new()),
            tandem: false,
            label: "BRCA1 Deletion Exon 13 (alte Zählung; 3.835 bp, Gründermutation)",
            deleted: 3835,
        },
    ]
}

const PATHOGENIC: [&str; 3] = ["alu_brca2_c156", "dup13_brca1", "del3835_brca1"];

/// SNPs on the haplotype of the 126 bp copy (chrom, position, reference, alternative)
const HAPLOTYPE_TAGS: [(&str, i64, char, char); 5] = [
    ("chr17", 43131837, 'G', 'A'),
    ("chr17", 43135863, 'A', 'G'),
    ("chr17", 43137747, 'A', 'G'),
    ("chr17", 43138596, 'G', 'T'),
    ("chr17", 43148782, 'C', 'A'),
];
const TAG_MARKER: i64 = 43128365;

/// one alignment record as `samtools view` prints it
#[derive(Clone)]
pub struct Rec {
    name: String,
    flag: u32,
    chrom: String,
    pos: i64, // 1-based
    mapq: u32,
    cigar: Vec<(i64, char)>,
    rnext_eq: bool,
    pnext: i64, // 1-based
    tlen: i64,
    seq: String,
    qual: Vec<u8>,
    sa: Option<String>,
}

fn op_char(o: Op) -> char {
    match o {
        Op::Match => 'M',
        Op::Ins => 'I',
        Op::Del => 'D',
        Op::Skip => 'N',
        Op::SoftClip => 'S',
        Op::HardClip => 'H',
        Op::Pad => 'P',
    }
}

fn to_rec(r: &Read, names: &[String]) -> Rec {
    let name_of = |t: i32| if t >= 0 { names.get(t as usize).cloned().unwrap_or_default() } else { "*".into() };
    Rec {
        name: String::from_utf8_lossy(&r.name).into_owned(),
        flag: r.flags as u32,
        chrom: name_of(r.tid),
        pos: r.start as i64 + 1,
        mapq: r.mapq as u32,
        cigar: r.cigar.iter().map(|&(o, l)| (l as i64, op_char(o))).collect(),
        rnext_eq: r.mate_tid >= 0 && r.mate_tid == r.tid,
        pnext: r.mate_start + 1,
        tlen: r.tlen,
        seq: String::from_utf8_lossy(&r.seq).to_ascii_uppercase(),
        qual: r.qual.clone(),
        sa: r.sa.clone(),
    }
}

fn parse_cigar(s: &str) -> Vec<(i64, char)> {
    let mut out = vec![];
    let mut n = 0i64;
    let mut any = false;
    for c in s.chars() {
        if let Some(d) = c.to_digit(10) {
            n = n * 10 + d as i64;
            any = true;
        } else if "MIDNSHP=X".contains(c) && any {
            out.push((n, c));
            n = 0;
            any = false;
        }
    }
    out
}

fn span_of(c: &[(i64, char)], ops: &str) -> i64 {
    c.iter().filter(|(_, o)| ops.contains(*o)).map(|(n, _)| n).sum()
}

type Region = (String, i64, i64, String);

fn bare(c: &str) -> &str {
    c.strip_prefix("chr").unwrap_or(c)
}

fn round2(x: f64) -> f64 {
    crate::exons::round(x, 2)
}

// ---------------------------------------------------------------------------------------------
// marker bases and haplotype tags

fn count_markers(recs: &[Rec], markers: &[(String, i64, char, char)], min_baseq: u8, min_mapq: u32) -> Vec<(i64, i64)> {
    let mut counts = vec![(0i64, 0i64); markers.len()];
    for r in recs {
        if r.flag & 0xF04 != 0 || r.cigar.is_empty() || r.mapq < min_mapq {
            continue;
        }
        let chrom = bare(&r.chrom);
        let hits: Vec<usize> = (0..markers.len()).filter(|&k| bare(&markers[k].0) == chrom && markers[k].1 >= r.pos).collect();
        if hits.is_empty() {
            continue;
        }
        let (mut rf, mut rd) = (r.pos, 0i64);
        for &(n, op) in &r.cigar {
            match op {
                'M' | '=' | 'X' => {
                    for &k in &hits {
                        let mp = markers[k].1;
                        if rf <= mp && mp < rf + n {
                            let i = (rd + mp - rf) as usize;
                            if r.qual.len() <= i || r.qual[i] >= min_baseq {
                                let base = r.seq.as_bytes().get(i).map(|b| *b as char).unwrap_or('N');
                                if base == markers[k].2 {
                                    counts[k].0 += 1;
                                } else if base == markers[k].3 {
                                    counts[k].1 += 1;
                                }
                            }
                        }
                    }
                    rf += n;
                    rd += n;
                }
                'D' | 'N' => rf += n,
                'I' | 'S' => rd += n,
                _ => {}
            }
        }
    }
    counts
}

/// (tags with the alternative base, tags tested, copies of the haplotype)
fn haplotype(counts: &[(i64, i64)]) -> (usize, usize, u8) {
    let alt: Vec<&(i64, i64)> = counts.iter().filter(|(r, a)| *a >= 3 && *a as f64 / (*r + *a) as f64 >= 0.2).collect();
    let hom = alt.iter().filter(|(r, a)| *a as f64 / (*r + *a) as f64 >= 0.8).count();
    let copies = if hom >= 3 { 2 } else if alt.len() >= 3 { 1 } else { 0 };
    (alt.len(), counts.len(), copies)
}

// ---------------------------------------------------------------------------------------------
// allele realignment

/// smallest edit distance of `pattern` to any substring of `text` (same result as Myers' algorithm)
fn infix_distance(pattern: &[u8], text: &[u8]) -> usize {
    let m = pattern.len();
    if m == 0 {
        return 0;
    }
    let mut col: Vec<usize> = (0..=m).collect();
    let mut best = m;
    for &t in text {
        let mut diag = col[0]; // free start: row 0 stays 0
        col[0] = 0;
        for i in 1..=m {
            let up = col[i];
            let v = (diag + if pattern[i - 1] == t { 0 } else { 1 }).min(col[i - 1] + 1).min(up + 1);
            diag = up;
            col[i] = v;
        }
        best = best.min(col[m]);
    }
    best
}

fn realign_known(recs: &[Rec], fasta: &HashMap<String, (i64, Vec<u8>)>, alleles: &[Known]) -> IndexMap<&'static str, (i64, i64)> {
    let (flank, before, after, min_diff, max_error, min_len) = (700i64, 400i64, 300i64, 2usize, 0.08, 50usize);
    struct Site<'a> {
        a: &'a Known,
        chrom: String,
        rf: Vec<u8>,
        alt: Vec<u8>,
    }
    let mut sites: Vec<Site> = vec![];
    for a in alleles {
        let Some(chrom) = fasta.keys().find(|c| bare(c) == bare(a.chrom)) else { continue };
        let (start, seq) = &fasta[chrom];
        let i = a.anchor - start;
        let deleted = a.deleted;
        if i - flank < 0 || i + deleted + flank > seq.len() as i64 {
            continue;
        }
        let (inserted, left): (Vec<u8>, i64) = match &a.ins {
            Ins::Seq(s) => (s.as_bytes().to_vec(), i - flank),
            Ins::Dup(ds, de) => {
                if ds - start - flank < 0 {
                    continue;
                }
                (seq[(ds - start) as usize..(de - start + 1) as usize].to_vec(), ds - start - flank)
            }
        };
        let rf = seq[left as usize..(i + deleted + flank) as usize].to_vec();
        let mut alt = seq[left as usize..(i + 1) as usize].to_vec();
        alt.extend(&inserted);
        alt.extend(&seq[(i + 1 + deleted) as usize..(i + deleted + flank) as usize]);
        sites.push(Site { a, chrom: bare(chrom).to_string(), rf, alt });
    }
    let mut counts: IndexMap<&'static str, (i64, i64)> = sites.iter().map(|s| (s.a.id, (0, 0))).collect();
    let mut seen: HashSet<(&'static str, String, u32)> = HashSet::new();
    for r in recs {
        if r.flag & 0xF00 != 0 || r.seq.len() < min_len {
            continue;
        }
        let chrom = bare(&r.chrom);
        for s in &sites {
            let a = s.a;
            let deleted = a.deleted;
            let dup_start = match a.ins {
                Ins::Dup(ds, _) => ds,
                _ => 0,
            };
            let near = a.anchor - before < r.pos && r.pos < a.anchor + after
                || deleted != 0 && a.anchor + deleted - before < r.pos && r.pos < a.anchor + deleted + after
                || dup_start != 0 && dup_start - before < r.pos && r.pos < dup_start + after;
            if s.chrom != chrom || !near {
                continue;
            }
            if !seen.insert((a.id, r.name.clone(), r.flag & 0xC0)) {
                continue;
            }
            let fwd = r.seq.as_bytes().to_vec();
            let mut variants = vec![fwd.clone()];
            if r.flag & 4 != 0 {
                variants.push(revcomp(&r.seq).into_bytes());
            }
            let (d_ref, d_alt) = variants
                .iter()
                .map(|x| (infix_distance(x, &s.rf), infix_distance(x, &s.alt)))
                .fold(None, |best: Option<(usize, usize)>, p| match best {
                    Some(b) if b.0.min(b.1) <= p.0.min(p.1) => Some(b),
                    _ => Some(p),
                })
                .unwrap();
            if d_ref.min(d_alt) as f64 > max_error * r.seq.len() as f64 {
                continue;
            }
            let c = counts.get_mut(a.id).unwrap();
            if d_alt + min_diff <= d_ref {
                c.1 += 1;
            } else if d_ref + min_diff <= d_alt {
                c.0 += 1;
            }
        }
    }
    counts
}

fn duplication_pairs(recs: &[Rec], alleles: &[Known]) -> IndexMap<&'static str, usize> {
    let (reach, min_mapq) = (700i64, 20u32);
    let sites: Vec<(&'static str, &str, i64, i64)> = alleles
        .iter()
        .filter_map(|a| match a.ins {
            Ins::Dup(s, e) => Some((a.id, bare(a.chrom), s, e)),
            _ => None,
        })
        .collect();
    let mut names: IndexMap<&'static str, HashSet<String>> = sites.iter().map(|s| (s.0, HashSet::new())).collect();
    for r in recs {
        let f = r.flag;
        if f & 0xF0C != 0 || r.mapq < min_mapq || !r.rnext_eq || f & 16 == 0 || f & 32 != 0 {
            continue;
        }
        let chrom = bare(&r.chrom);
        for &(id, c, start, end) in &sites {
            if c == chrom && start - 50 <= r.pos && r.pos <= start + reach && end - reach - 150 <= r.pnext && r.pnext <= end {
                names[id].insert(r.name.clone());
            }
        }
    }
    names.into_iter().map(|(k, v)| (k, v.len())).collect()
}

// ---------------------------------------------------------------------------------------------
// read pairs and split reads

fn in_regions(regions: &[Region], chrom: &str, pos: i64) -> bool {
    regions.iter().any(|(c, s, e, _)| c == chrom && *s < pos && pos <= *e)
}

fn gene_at(regions: &[Region], chrom: &str, pos: i64) -> String {
    regions.iter().find(|(c, s, e, _)| c == chrom && *s < pos && pos <= *e).map(|r| r.3.clone()).unwrap_or_default()
}

fn sa_entries(sa: &str) -> Vec<Vec<String>> {
    sa.trim_end_matches(';').split(';').map(|e| e.split(',').map(String::from).collect()).collect()
}

type Cluster = (String, i64, i64, Vec<String>, Vec<String>);

fn outward_pairs(recs: &[Rec], regions: &[Region]) -> Vec<Cluster> {
    let (min_mapq, window) = (20u32, 600i64);
    let mut pairs: Vec<(String, i64, i64, String)> = vec![];
    let mut clips: Vec<(String, i64, i64, String)> = vec![];
    let mut seen: HashSet<String> = HashSet::new();
    for r in recs {
        let f = r.flag;
        if f & 0xF0C != 0 || r.mapq < min_mapq || !in_regions(regions, &r.chrom, r.pos) {
            continue;
        }
        if f & 1 != 0 && r.rnext_eq && f & 16 != 0 && f & 32 == 0 && r.pnext > r.pos && !seen.contains(&r.name) {
            seen.insert(r.name.clone());
            pairs.push((r.chrom.clone(), r.pos, r.pnext + r.seq.len() as i64, r.name.clone()));
        }
        if let Some(sa) = &r.sa {
            let aligned = span_of(&r.cigar, "MD=X");
            for p in sa_entries(sa) {
                if p.len() >= 5 && p[0] == r.chrom && p[4].parse::<u32>().unwrap_or(0) >= min_mapq {
                    let other: i64 = p[1].parse().unwrap_or(0);
                    if r.cigar.last().is_some_and(|c| c.1 == 'S') && other < r.pos {
                        clips.push((r.chrom.clone(), other, r.pos + aligned, r.name.clone()));
                    } else if r.cigar.first().is_some_and(|c| c.1 == 'S') && other > r.pos {
                        clips.push((r.chrom.clone(), r.pos, other + aligned, r.name.clone()));
                    }
                }
            }
        }
    }
    pairs.sort();
    let mut clusters: Vec<Cluster> = vec![];
    for (chrom, start, end, name) in pairs {
        if let Some(c) = clusters.iter_mut().find(|c| c.0 == chrom && (c.1 - start).abs() <= window && (c.2 - end).abs() <= window) {
            c.1 = c.1.min(start);
            c.2 = c.2.max(end);
            if !c.3.contains(&name) {
                c.3.push(name);
            }
        } else {
            clusters.push((chrom, start, end, vec![name], vec![]));
        }
    }
    for (chrom, start, end, name) in clips {
        if let Some(c) = clusters.iter_mut().find(|c| c.0 == chrom && (c.1 - start).abs() <= window && (c.2 - end).abs() <= window) {
            if !c.4.contains(&name) {
                c.4.push(name);
            }
        }
    }
    clusters
}

fn tandem_duplications(recs: &[Rec], regions: &[Region]) -> Vec<Value> {
    let mut out = vec![];
    for (chrom, start, end, names, split) in outward_pairs(recs, regions) {
        let support: HashSet<&String> = names.iter().chain(split.iter()).collect();
        if support.len() < 4 || !(1000..=200000).contains(&(end - start)) {
            continue;
        }
        let mut ev = vec!["read pairs"];
        if !split.is_empty() {
            ev.push("split reads");
        }
        out.push(json!({"gene": gene_at(regions, &chrom, start), "chrom": chrom, "pos": start, "end": end,
            "duplicatedLength": end - start, "signature": ["tandem duplication", "read pairs"], "evidence": ev,
            "outwardPairs": names.len(), "splitReads": split.len(), "breakpointPrecision": "approximate (read-pair span)",
            "clippedReadsLeftJunction": 0, "clippedReadsRightJunction": 0}));
    }
    out
}

fn insert_stats(mut v: Vec<i64>) -> (f64, f64) {
    v.sort();
    if v.len() < 100 {
        return (450.0, 100.0);
    }
    let med = v[v.len() / 2];
    let step = (v.len() / 5000).max(1);
    let sample: Vec<i64> = v.iter().step_by(step).cloned().collect();
    let mut dev: Vec<i64> = sample.iter().map(|x| (x - med).abs()).collect();
    dev.sort();
    let mad = dev[sample.len() / 2];
    (med as f64, (1.4826 * mad as f64).max(20.0))
}

type Item = (String, i64, i64, String);

struct Signals {
    dpairs: Vec<Item>,
    ipairs: Vec<(String, i64, i64, String, &'static str)>,
    dsplit: Vec<Item>,
    isplit: Vec<Item>,
    med: f64,
    sd: f64,
}

fn sv_signals(recs: &[Rec], regions: &[Region]) -> Signals {
    let min_mapq = 20u32;
    let mut kept: Vec<(&Rec, i64)> = vec![];
    let mut tlens = vec![];
    for r in recs {
        if r.flag & 0xF0C != 0 || r.mapq < min_mapq || !in_regions(regions, &r.chrom, r.pos) {
            continue;
        }
        let span = span_of(&r.cigar, "MD=X");
        kept.push((r, span));
        if r.flag & 2 != 0 && r.rnext_eq && 0 < r.tlen && r.tlen < 2000 {
            tlens.push(r.tlen);
        }
    }
    let (med, sd) = insert_stats(tlens);
    let cut = med + 6.0 * sd;
    let mut s = Signals { dpairs: vec![], ipairs: vec![], dsplit: vec![], isplit: vec![], med, sd };
    for (r, span) in kept {
        let (chrom, pos, name, flag) = (&r.chrom, r.pos, &r.name, r.flag);
        if flag & 1 != 0 && r.rnext_eq && flag & 8 == 0 {
            let mpos = r.pnext;
            let (rev, mrev) = (flag & 16 != 0, flag & 32 != 0);
            if mpos > pos {
                if !rev && mrev && (mpos + r.seq.len() as i64 - pos) as f64 > cut {
                    s.dpairs.push((chrom.clone(), pos + span, mpos, name.clone()));
                } else if rev == mrev && mpos - pos > 300 {
                    s.ipairs.push((chrom.clone(), pos, mpos, name.clone(), if rev { "--" } else { "++" }));
                }
            }
        }
        let (Some(sa), false) = (&r.sa, r.cigar.is_empty()) else { continue };
        let strand = if flag & 16 != 0 { "-" } else { "+" };
        let cig = &r.cigar;
        let sh = |c: &(i64, char)| c.1 == 'S' || c.1 == 'H';
        for p in sa_entries(sa) {
            if p.len() < 6 || &p[0] != chrom || p[4].parse::<u32>().unwrap_or(0) < min_mapq {
                continue;
            }
            let opos: i64 = p[1].parse().unwrap_or(0);
            let ocig = parse_cigar(&p[3]);
            let ospan = span_of(&ocig, "MD=X");
            if ocig.is_empty() {
                continue;
            }
            if p[2] == strand {
                if sh(cig.last().unwrap()) && opos > pos + span + 30 && sh(&ocig[0]) {
                    s.dsplit.push((chrom.clone(), pos + span, opos, name.clone()));
                } else if sh(&cig[0]) && opos + ospan + 30 < pos && sh(ocig.last().unwrap()) {
                    s.dsplit.push((chrom.clone(), opos + ospan, pos, name.clone()));
                }
            } else {
                let junction = |start: i64, c: &[(i64, char)], length: i64| {
                    let right = sh(c.last().unwrap()) && !(sh(&c[0]) && c[0].0 > c.last().unwrap().0);
                    if right { start + length } else { start }
                };
                let (mut a, mut b) = (junction(pos, cig, span), junction(opos, &ocig, ospan));
                if a > b {
                    std::mem::swap(&mut a, &mut b);
                }
                if b - a > 300 {
                    s.isplit.push((chrom.clone(), a, b, name.clone()));
                }
            }
        }
    }
    s
}

/// clusters: (chrom, first item's start, first item's end, items) after sorting the items
fn clusters<T: Clone + Ord>(items: &[T], key: &dyn Fn(&T) -> (String, i64, i64), window: i64) -> Vec<(String, i64, i64, Vec<T>)> {
    let mut sorted = items.to_vec();
    sorted.sort();
    let mut out: Vec<(String, i64, i64, Vec<T>)> = vec![];
    for it in sorted {
        let (c, a, b) = key(&it);
        if let Some(h) = out.iter_mut().find(|h| h.0 == c && (h.1 - a).abs() <= window && (h.2 - b).abs() <= window) {
            h.3.push(it);
        } else {
            out.push((c, a, b, vec![it]));
        }
    }
    out
}

#[allow(clippy::too_many_arguments)]
fn sv_entry(regions: &[Region], chrom: &str, start: i64, end: i64, signature: Vec<&str>, pairs: usize, split: usize, precise: bool, kind_len: &str) -> Value {
    let mut ev = vec![];
    if pairs > 0 {
        ev.push("read pairs");
    }
    if split > 0 {
        ev.push("split reads");
    }
    let mut v = json!({"gene": gene_at(regions, chrom, start), "chrom": chrom, "pos": start, "end": end,
        "signature": signature, "evidence": ev, "discordantPairs": pairs, "splitReads": split,
        "breakpointPrecision": if precise { "exact (split reads)" } else { "approximate (read-pair span)" },
        "clippedReadsLeftJunction": 0, "clippedReadsRightJunction": 0});
    v[kind_len] = json!(end - start + 1);
    v
}

fn read_pair_svs(recs: &[Rec], regions: &[Region]) -> Vec<Value> {
    let (min_support, min_del, min_inv, max_size) = (4usize, 50i64, 300i64, 200000i64);
    let s = sv_signals(recs, regions);
    let w = (s.med + 6.0 * s.sd) as i64;
    let key = |x: &Item| (x.0.clone(), x.1, x.2);
    let mut out: Vec<Value> = vec![];
    let mut used: HashSet<String> = HashSet::new();
    for (chrom, _a, _b, its) in clusters(&s.dsplit, &key, 20) {
        // most common (start, end), first seen on ties
        let mut cnt: IndexMap<(i64, i64), usize> = IndexMap::new();
        for x in &its {
            *cnt.entry((x.1, x.2)).or_insert(0) += 1;
        }
        let mx = *cnt.values().max().unwrap();
        let (st, en) = *cnt.iter().find(|(_, v)| **v == mx).unwrap().0;
        let names: HashSet<&String> = its.iter().map(|x| &x.3).collect();
        let pnames: HashSet<&String> =
            s.dpairs.iter().filter(|x| x.0 == chrom && st - w <= x.1 && x.1 <= st + 20 && en - 20 <= x.2 && x.2 <= en + w).map(|x| &x.3).collect();
        if names.union(&pnames).count() >= min_support && min_del <= en - st && en - st <= max_size {
            out.push(sv_entry(regions, &chrom, st, en - 1, vec!["read-pair deletion"], pnames.len(), names.len(), true, "deletedLength"));
            used.extend(pnames.into_iter().cloned());
        }
    }
    let rest: Vec<Item> = s.dpairs.iter().filter(|x| !used.contains(&x.3)).cloned().collect();
    for (chrom, _a, _b, its) in clusters(&rest, &key, w) {
        let names: HashSet<&String> = its.iter().map(|x| &x.3).collect();
        let st = its.iter().map(|x| x.1).max().unwrap();
        let en = its.iter().map(|x| x.2).min().unwrap();
        let near = out.iter().any(|o| o["chrom"] == chrom.as_str() && (o["pos"].as_i64().unwrap() - st).abs() < w && (o["end"].as_i64().unwrap() - en).abs() < w);
        if names.len() >= min_support && min_del <= en - st && en - st <= max_size && !near {
            out.push(sv_entry(regions, &chrom, st, en - 1, vec!["read-pair deletion"], names.len(), 0, false, "deletedLength"));
        }
    }
    let wi = (s.med + 5.0 * s.sd) as i64;
    let mut invs: Vec<Value> = vec![];
    let ikey = |x: &(String, i64, i64, String, &'static str)| (x.0.clone(), x.1, x.2);
    for (chrom, a, b, its) in clusters(&s.ipairs, &ikey, wi) {
        let names: HashSet<&String> = its.iter().map(|x| &x.3).collect();
        let pp: Vec<_> = its.iter().filter(|x| x.4 == "++").collect();
        let mm: Vec<_> = its.iter().filter(|x| x.4 == "--").collect();
        let sp: Vec<&Item> = s.isplit.iter().filter(|x| x.0 == chrom && (x.1 - a).abs() <= wi && (x.2 - b).abs() <= wi).collect();
        let (mut aa, mut bb) = (a, b);
        if !pp.is_empty() && !mm.is_empty() {
            aa = (pp.iter().map(|x| x.1).max().unwrap() + 150 + mm.iter().map(|x| x.1).min().unwrap()).div_euclid(2);
            bb = (pp.iter().map(|x| x.2).max().unwrap() + 150 + mm.iter().map(|x| x.2).min().unwrap()).div_euclid(2);
        }
        if !sp.is_empty() {
            let mut xs: Vec<i64> = sp.iter().map(|x| x.1).collect();
            let mut ys: Vec<i64> = sp.iter().map(|x| x.2).collect();
            xs.sort();
            ys.sort();
            aa = xs[sp.len() / 2];
            bb = ys[sp.len() / 2];
        }
        if invs.iter().any(|o| o["chrom"] == chrom.as_str() && (o["pos"].as_i64().unwrap() - aa).abs() <= wi && (o["end"].as_i64().unwrap() - bb).abs() <= wi) {
            continue;
        }
        let spn: HashSet<&String> = sp.iter().map(|x| &x.3).collect();
        if names.union(&spn).count() >= min_support && min_inv <= bb - aa && bb - aa <= max_size {
            invs.push(sv_entry(regions, &chrom, aa, bb - 1, vec!["inversion", "read pairs"], names.len(), spn.len(), !sp.is_empty(), "invertedLength"));
        }
    }
    out.extend(invs);
    out
}

fn realigned_alleles(counts: &IndexMap<&'static str, (i64, i64)>, regions: &[Region], alleles: &[Known], pairs: &IndexMap<&'static str, usize>) -> Vec<Value> {
    let (min_reads, min_pairs) = (2i64, 3usize);
    let mut out = vec![];
    for a in alleles {
        let Some(&(n_ref, n_alt)) = counts.get(a.id) else { continue };
        if a.tandem {
            continue;
        }
        let n_pairs = pairs.get(a.id).copied().unwrap_or(0);
        let is_dup = matches!(a.ins, Ins::Dup(..));
        if n_alt < min_reads && !(is_dup && n_pairs >= min_pairs && n_alt >= 1) {
            continue;
        }
        let chrom = regions.iter().find(|r| bare(&r.0) == bare(a.chrom)).map(|r| r.0.clone()).unwrap_or(a.chrom.to_string());
        let gene = gene_at(regions, &chrom, a.anchor);
        let deleted = a.deleted;
        let mut kind: Vec<&str> = if a.id.starts_with("alu") { vec!["Alu", "known allele"] } else { vec!["known allele"] };
        if deleted != 0 {
            kind.push("deletion");
        }
        let duplicated = match a.ins {
            Ins::Dup(s, e) => e - s + 1,
            _ => 0,
        };
        if duplicated != 0 {
            kind.push("duplication");
        }
        let ins_len = match &a.ins {
            Ins::Seq(s) => s.len() as i64,
            Ins::Dup(..) => duplicated,
        };
        let (pos, end) = match a.ins {
            Ins::Dup(s, e) => (s, e),
            _ => (a.anchor + if deleted != 0 { 1 } else { 0 }, a.anchor + deleted.max(1)),
        };
        let mut ev = vec!["allele realignment"];
        if n_pairs >= min_pairs {
            ev.push("read pairs");
        }
        let zyg = if duplicated != 0 {
            "unknown"
        } else if n_ref as f64 <= 0.1 * (n_ref + n_alt) as f64 && n_alt >= 2 * min_reads {
            "homozygous"
        } else {
            "heterozygous"
        };
        let mut v = json!({"gene": gene, "chrom": chrom, "pos": pos, "end": end, "insertedLength": ins_len, "deletedLength": deleted,
            "targetSiteDuplication": 0, "clippedReadsLeftJunction": 0, "clippedReadsRightJunction": 0, "spanningReads": n_ref,
            "clippedFraction": 0.0, "label": a.label, "signature": kind, "clipsWithAluSequence": 0, "clipsWithPolyA": 0,
            "evidence": ev, "realignedReads": format!("{n_alt} allele / {n_ref} reference"),
            "clinicalSignificance": if PATHOGENIC.contains(&a.id) { "pathogenic" } else { "" }, "zygosity": zyg});
        if is_dup {
            v["discordantPairs"] = json!(n_pairs);
        }
        out.push(v);
    }
    out
}

// ---------------------------------------------------------------------------------------------
// clipped read ends

type Clips = IndexMap<(String, i64), Vec<String>>;

fn scan(recs: &[Rec], min_clip: i64) -> (Clips, Clips, HashMap<String, Vec<(i64, i64)>>) {
    let mut right: Clips = IndexMap::new();
    let mut left: Clips = IndexMap::new();
    let mut starts: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
    for r in recs {
        if r.flag & 0xF04 != 0 || r.mapq < 20 || r.cigar.is_empty() {
            continue;
        }
        let ref_len = span_of(&r.cigar, "MDN=X");
        let end = r.pos + ref_len - 1;
        starts.entry(r.chrom.clone()).or_default().push((r.pos, end));
        let first = r.cigar[0];
        if first.1 == 'S' && first.0 >= min_clip {
            left.entry((r.chrom.clone(), r.pos)).or_default().push(r.seq[..first.0 as usize].to_string());
        }
        let last = *r.cigar.last().unwrap();
        if last.1 == 'S' && last.0 >= min_clip {
            right.entry((r.chrom.clone(), end)).or_default().push(r.seq[r.seq.len() - last.0 as usize..].to_string());
        }
    }
    (right, left, starts)
}

fn describe(clips: &[String], kmers: &[String]) -> (usize, usize) {
    let alu = clips.iter().filter(|c| kmers.iter().any(|k| c.contains(k.as_str()))).count();
    let poly = clips.iter().filter(|c| c.contains(&"A".repeat(10)) || c.contains(&"T".repeat(10))).count();
    (alu, poly)
}

fn candidates(right: &Clips, left: &Clips, starts: &HashMap<String, Vec<(i64, i64)>>, regions: &[Region], kmers: &[String]) -> Vec<Value> {
    let (min_reads, max_tsd, min_fraction) = (3usize, 30i64, 0.15);
    let mut left_by_chrom: HashMap<String, Vec<i64>> = HashMap::new();
    for ((chrom, pos), clips) in left {
        if clips.len() >= min_reads {
            left_by_chrom.entry(chrom.clone()).or_default().push(*pos);
        }
    }
    let mut out = vec![];
    let mut used: HashSet<(String, i64)> = HashSet::new();
    let mut keys: Vec<&(String, i64)> = right.keys().collect();
    keys.sort();
    for k in keys {
        let (chrom, r_pos) = (&k.0, k.1);
        let r_clips = &right[k];
        if r_clips.len() < min_reads {
            continue;
        }
        let gene = gene_at(regions, chrom, r_pos);
        if gene.is_empty() {
            continue;
        }
        let partners: Vec<i64> = left_by_chrom.get(chrom).map(|v| v.iter().filter(|&&p| -1 <= r_pos - p && r_pos - p <= max_tsd).cloned().collect()).unwrap_or_default();
        if partners.is_empty() {
            continue;
        }
        let mut l_pos = partners[0];
        for &p in &partners {
            if left[&(chrom.clone(), p)].len() > left[&(chrom.clone(), l_pos)].len() {
                l_pos = p;
            }
        }
        if used.contains(&(chrom.clone(), l_pos)) {
            continue;
        }
        let l_clips = &left[&(chrom.clone(), l_pos)];
        let (lo, hi) = (l_pos.min(r_pos), l_pos.max(r_pos));
        let depth = starts.get(chrom).map(|v| v.iter().filter(|(s, e)| *s < lo && *e > hi).count()).unwrap_or(0) + r_clips.len() + l_clips.len();
        let fraction = if depth > 0 { (r_clips.len() + l_clips.len()) as f64 / depth as f64 } else { 0.0 };
        if fraction < min_fraction {
            continue;
        }
        used.insert((chrom.clone(), l_pos));
        let all: Vec<String> = r_clips.iter().chain(l_clips.iter()).cloned().collect();
        let (alu, poly) = describe(&all, kmers);
        let mut sig = vec![];
        if alu >= 2 {
            sig.push("Alu");
        }
        if poly >= 2 {
            sig.push("poly-A");
        }
        if !sig.contains(&"Alu") {
            continue;
        }
        out.push(json!({"gene": gene, "chrom": chrom, "pos": if l_pos <= r_pos { l_pos } else { r_pos }, "end": hi,
            "targetSiteDuplication": (r_pos - l_pos + 1).max(0), "clippedReadsLeftJunction": r_clips.len(),
            "clippedReadsRightJunction": l_clips.len(), "spanningReads": depth - r_clips.len() - l_clips.len(),
            "clippedFraction": round2(fraction), "signature": sig, "clipsWithAluSequence": alu, "clipsWithPolyA": poly}));
    }
    out
}

fn merge_nearby(clips: &Clips) -> Clips {
    let mut keys: Vec<&(String, i64)> = clips.keys().collect();
    keys.sort_by_key(|k| std::cmp::Reverse(clips[*k].len())); // stable: insertion order on ties
    let mut out: Clips = IndexMap::new();
    for key in keys {
        let home = out.keys().find(|k| k.0 == key.0 && (k.1 - key.1).abs() <= 10).cloned();
        match home {
            Some(h) => out[&h].extend(clips[key].iter().cloned()),
            None => {
                out.insert(key.clone(), clips[key].clone());
            }
        }
    }
    out
}

fn place(clip: &[u8], start: i64, rf: &[u8], lo: i64, hi: i64, upstream: bool) -> Option<(i64, f64)> {
    let k = clip.len() as i64;
    let mut best: Option<(i64, f64)> = None;
    let ps: Vec<i64> = if upstream { (lo..hi).rev().collect() } else { (lo..hi).collect() };
    for p in ps {
        let i = p - start;
        let (s, e) = if upstream { (i, i + k) } else { (i - k + 1, i + 1) };
        if s < 0 || e > rf.len() as i64 {
            continue;
        }
        let seg = &rf[s as usize..e as usize];
        let ident = clip.iter().zip(seg).filter(|(a, b)| a == b).count() as f64 / k as f64;
        if ident >= 0.9 && best.is_none_or(|b| ident > b.1) {
            best = Some((p, ident));
        }
    }
    best
}

fn lengths(clips: &Clips, fasta: &HashMap<String, (i64, Vec<u8>)>, right_side: bool) -> Vec<(String, i64, i64, usize)> {
    let (min_reads, min_clip, min_len, max_len, seed) = (2usize, 10usize, 50i64, 1000i64, 40usize);
    let mut out = vec![];
    for ((chrom, pos), cl) in clips {
        let cl: Vec<&String> = cl.iter().filter(|c| c.len() >= min_clip).collect();
        let Some((start, rf)) = fasta.get(chrom) else { continue };
        if cl.len() < min_reads {
            continue;
        }
        let mut ls: Vec<i64> = vec![];
        for c in cl {
            let b = c.as_bytes();
            if right_side {
                if let Some(h) = place(&b[..seed.min(b.len())], *start, rf, pos - max_len + 1, pos - min_len + 2, true) {
                    ls.push(pos - h.0 + 1);
                }
            } else if let Some(h) = place(&b[b.len().saturating_sub(seed)..], *start, rf, pos + min_len - 1, pos + max_len, false) {
                ls.push(h.0 - pos + 1);
            }
        }
        if ls.len() >= min_reads {
            let mut cnt: IndexMap<i64, usize> = IndexMap::new();
            for &l in &ls {
                *cnt.entry(l).or_insert(0) += 1;
            }
            let mx = *cnt.values().max().unwrap();
            let length = *cnt.iter().find(|(_, v)| **v == mx).unwrap().0;
            let support = ls.iter().filter(|&&x| ((x - length).abs() as f64) <= 0.1 * length as f64).count();
            if support >= min_reads {
                out.push((chrom.clone(), *pos, length, support));
            }
        }
    }
    out
}

fn duplications(right: &Clips, left: &Clips, starts: &HashMap<String, Vec<(i64, i64)>>, regions: &[Region], fasta: &HashMap<String, (i64, Vec<u8>)>) -> Vec<Value> {
    let mut rights = lengths(&merge_nearby(right), fasta, true);
    let lefts = lengths(&merge_nearby(left), fasta, false);
    rights.sort_by_key(|x| std::cmp::Reverse(x.3));
    let mut out = vec![];
    let mut used: HashSet<(String, i64)> = HashSet::new();
    for (chrom, b, length, n_r) in rights {
        let gene = gene_at(regions, &chrom, b);
        if gene.is_empty() {
            continue;
        }
        let partners: Vec<&(String, i64, i64, usize)> = lefts
            .iter()
            .filter(|x| x.0 == chrom && ((x.2 - length).abs() as f64) <= 0.1 * length as f64 && (b - x.1).abs() <= 2 * length && !used.contains(&(x.0.clone(), x.1)))
            .collect();
        if partners.is_empty() {
            continue;
        }
        let mut best = partners[0];
        for p in &partners {
            if p.3 > best.3 {
                best = p;
            }
        }
        let (a, n_l) = (best.1, best.3);
        used.insert((chrom.clone(), a));
        let (lo, hi) = (a.min(b - length + 1), (a + length - 1).max(b));
        let spanning = starts.get(&chrom).map(|v| v.iter().filter(|(s, e)| *s < b - 10 && *e > b + 10).count()).unwrap_or(0);
        let tot = spanning + n_r + n_l;
        out.push(json!({"gene": gene, "chrom": chrom, "pos": lo, "end": hi, "duplicatedLength": length,
            "targetSiteDuplication": 0, "clippedReadsLeftJunction": n_r, "clippedReadsRightJunction": n_l,
            "spanningReads": spanning, "clippedFraction": if tot > 0 { round2((n_r + n_l) as f64 / tot as f64) } else { 0.0 },
            "signature": ["tandem copy"], "clipsWithAluSequence": 0, "clipsWithPolyA": 0}));
    }
    out
}

fn has_sig(d: &Value, s: &str) -> bool {
    d["signature"].as_array().is_some_and(|a| a.iter().any(|x| x == s))
}

/// all signatures of the gene regions (the JSON list of 1.1's insertion_signals.json)
pub fn signals(path: &str, reference: &str, regions_bed: &[(String, u64, u64, String)]) -> io::Result<Vec<Value>> {
    let names: Vec<String> = reads::contigs(path, Some(reference))?.into_iter().map(|c| c.0).collect();
    let regions: Vec<Region> = regions_bed.iter().map(|(c, s, e, g)| (c.clone(), *s as i64, *e as i64, g.clone())).collect();
    let mut recs: Vec<Rec> = vec![];
    for (c, s, e, _) in &regions {
        reads::fetch(path, Some(reference), &format!("{c}:{}-{e}", s + 1), &mut |r: &Read| recs.push(to_rec(r, &names)))?;
    }
    // reference of the gene regions (as samtools faidx chrom:start-end)
    let mut fa = noodles::fasta::io::indexed_reader::Builder::default().build_from_path(reference)?;
    let mut fasta: HashMap<String, (i64, Vec<u8>)> = HashMap::new();
    for (c, s, e, _) in &regions {
        let region: noodles::core::Region = format!("{c}:{}-{e}", s + 1).parse().map_err(io::Error::other)?;
        let seq: Vec<u8> = fa.query(&region)?.sequence().as_ref().iter().map(|b| b.to_ascii_uppercase()).collect();
        fasta.insert(c.clone(), (s + 1, seq));
    }
    let kmers = alu_kmers();
    let alleles = known_alleles();
    let (right, left, starts) = scan(&recs, 10);
    let mut found = candidates(&right, &left, &starts, &regions, &kmers);
    let mut dups = duplications(&right, &left, &starts, &regions, &fasta);
    let same_copy = |d: &Value, m: &Copy| {
        bare(d["chrom"].as_str().unwrap_or("")) == bare(m.0)
            && (d["end"].as_i64().unwrap_or(0) - m.5).abs() <= 20
            && (d["duplicatedLength"].as_i64().unwrap_or(0) - m.6).abs() <= 10
    };
    let clipped: HashSet<i64> = KNOWN_COPIES.iter().filter(|m| dups.iter().any(|d| same_copy(d, m))).map(|m| m.1).collect();
    let tags: Vec<(String, i64, char, char)> = HAPLOTYPE_TAGS.iter().map(|t| (t.0.to_string(), t.1, t.2, t.3)).collect();
    let hap = haplotype(&count_markers(&recs, &tags, 13, 20));
    let realigned = realign_known(&recs, &fasta, &alleles);
    let by_marker: HashMap<i64, (i64, i64)> = realigned.get("tandem126").map(|v| (TAG_MARKER, *v)).into_iter().collect();
    let markers: Vec<(String, i64, char, char)> = KNOWN_COPIES.iter().map(|m| (m.0.to_string(), m.1, m.2, m.3)).collect();
    let mcounts = count_markers(&recs, &markers, 13, 0);
    for (k, m) in KNOWN_COPIES.iter().enumerate() {
        let (n_ref, n_alt) = mcounts[k];
        let depth = n_ref + n_alt;
        let marker = n_alt >= 2 && depth > 0 && n_alt as f64 / depth as f64 >= 0.05;
        let h = if m.1 == TAG_MARKER { Some(hap) } else { None };
        let mut evidence: Vec<&str> = vec![];
        if clipped.contains(&m.1) {
            evidence.push("clipped reads");
        }
        if marker {
            evidence.push("marker base");
        }
        if h.is_some_and(|x| x.2 > 0) {
            evidence.push("haplotype");
        }
        let n_re = by_marker.get(&m.1).map(|v| v.1).unwrap_or(0);
        if evidence.len() < 2 && n_re < 2 {
            continue;
        }
        if n_re >= 2 {
            evidence.insert(0, "allele realignment");
        }
        let chrom = regions.iter().find(|r| bare(&r.0) == bare(m.0)).map(|r| r.0.clone()).unwrap_or(m.0.to_string());
        let gene = gene_at(&regions, &chrom, m.1);
        let mut known = json!({"gene": gene, "chrom": chrom, "pos": m.4, "end": m.5, "duplicatedLength": m.6,
            "targetSiteDuplication": 0, "clippedReadsLeftJunction": 0, "clippedReadsRightJunction": 0,
            "spanningReads": depth, "clippedFraction": 0.0, "markerPosition": m.1, "markerReads": format!("{n_alt}/{depth}"),
            "label": m.7, "signature": ["tandem copy", "known allele"], "clipsWithAluSequence": 0, "clipsWithPolyA": 0,
            "evidence": evidence});
        if by_marker.contains_key(&m.1) {
            known["realignedReads"] = json!(format!("{n_re} allele"));
        }
        if let Some(x) = h {
            known["haplotypeTags"] = json!(format!("{}/{}", x.0, x.1));
            known["haplotypeCopies"] = json!(x.2);
        }
        if let Some(d) = dups.iter_mut().find(|d| same_copy(d, m)) {
            for key in ["markerPosition", "markerReads", "label", "signature", "evidence", "haplotypeTags", "haplotypeCopies", "realignedReads"] {
                if !known[key].is_null() {
                    d[key] = known[key].clone();
                }
            }
        } else {
            dups.push(known);
        }
    }
    found.extend(dups);
    let pairs = duplication_pairs(&recs, &alleles);
    for known in realigned_alleles(&realigned, &regions, &alleles, &pairs) {
        let kpos = known["pos"].as_i64().unwrap_or(0);
        let blocked = has_sig(&known, "deletion") || has_sig(&known, "duplication");
        let same = found.iter_mut().find(|d| {
            d["chrom"] == known["chrom"] && !has_sig(d, "tandem copy") && !blocked && (d["pos"].as_i64().unwrap_or(0) - kpos).abs() <= 30
        });
        match same {
            Some(d) => {
                for key in ["label", "evidence", "realignedReads", "zygosity", "insertedLength", "clinicalSignificance"] {
                    d[key] = known[key].clone();
                }
                d["evidence"] = json!(["allele realignment", "clipped reads"]);
                let mut sig: Vec<String> = d["signature"].as_array().unwrap().iter().filter_map(|x| x.as_str().map(String::from)).collect();
                sig.push("known allele".into());
                sig.sort();
                sig.dedup();
                d["signature"] = json!(sig);
            }
            None => found.push(known),
        }
    }
    for dup in tandem_duplications(&recs, &regions) {
        let covered = found.iter().any(|d| {
            d["chrom"] == dup["chrom"] && has_sig(d, "duplication")
                && (d["pos"].as_i64().unwrap_or(0) - dup["pos"].as_i64().unwrap_or(0)).abs() <= 600
                && (d["end"].as_i64().unwrap_or(0) - dup["end"].as_i64().unwrap_or(0)).abs() <= 600
        });
        if !covered {
            found.push(dup);
        }
    }
    for sv in read_pair_svs(&recs, &regions) {
        let covered = found.iter().any(|d| {
            d["chrom"] == sv["chrom"] && has_sig(d, "deletion")
                && (d["pos"].as_i64().unwrap_or(0) - sv["pos"].as_i64().unwrap_or(0)).abs() <= 600
                && (d["end"].as_i64().unwrap_or(0) - sv["end"].as_i64().unwrap_or(0)).abs() <= 600
        });
        if !covered {
            found.push(sv);
        }
    }
    Ok(found)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn infix_distance_matches_definition() {
        assert_eq!(infix_distance(b"ACGT", b"TTACGTTT"), 0);
        assert_eq!(infix_distance(b"ACGT", b"TTACTTT"), 1);
        assert_eq!(infix_distance(b"AAAA", b"CCCC"), 4);
        assert_eq!(infix_distance(b"ACGTACGT", b"ACGACGT"), 1);
    }
}
