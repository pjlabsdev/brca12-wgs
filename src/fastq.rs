//! Paired FASTQ input: the read pairs of BRCA1/BRCA2 are pulled out of a whole-genome FASTQ and
//! aligned inside the program, after hla-wgs. Research use only.
//!
//! 1. Targets: the gene regions +/- 2 kb (BRCA1 with the pseudogene BRCA1P1) and the control windows
//!    of the exon copy-number panel.
//! 2. One pass over the reference gives the genome-unique 31-mers and the decoys. Genome-unique 31-mers: every 31-mer of the targets that occurs exactly once in the whole
//!    reference (both strands). Reads of repeats from elsewhere (Alu, L1) do not share such a k-mer
//!    and are not pulled in, which replaces the mapping quality a whole-genome alignment would give.
//!    Decoys: loci elsewhere that share many low-copy (2-20 copies) 31-mers with the targets
//!    (pseudogenes, segmental duplications); they join the alignment panel, and reads that fit a decoy
//!    best are dropped, as a whole-genome aligner would place them there. Both are cached next to
//!    the reference (REF.brca12-kmers) after the first run.
//! 3. A pair is kept when one of its reads has at least two of these k-mers; kept pairs are aligned
//!    with the own aligner (src/aln.rs from hla-wgs: k-mer anchors, chaining, affine alignment; paired with mate rescue) to the targets, placed
//!    on GRCh38, duplicates marked as samtools markdup does for pairs, sorted and indexed. A read
//!    without two genome-unique 31-mers of its own, or whose mate does not align to the targets, gets
//!    mapping quality 0 (it may come from a copy elsewhere in the genome).
//! The result is OUTPREFIX.bam (+ .bai), which the analysis then reads like any BAM.
use crate::exon_cnv;
use ahash::{AHashMap, AHashSet};
use flate2::read::MultiGzDecoder;
use crate::aln;
use noodles::sam;
use rayon::prelude::*;
use std::io::{self, BufRead, BufReader, Write};

const K: usize = 31;
const FLANK: u64 = 2000;

pub struct Target {
    pub name: String,
    pub chrom: String,
    pub start: u64, // 0-based
    pub seq: Vec<u8>,
}

fn enc(c: u8) -> Option<u64> {
    match c {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// canonical 31-mers of a sequence (k-mers with N skipped)
fn kmers(s: &[u8], mut f: impl FnMut(u64)) {
    let mask = (1u64 << (2 * K)) - 1;
    let shift = 2 * (K - 1);
    let (mut fw, mut rv, mut valid) = (0u64, 0u64, 0usize);
    for &c in s {
        match enc(c) {
            Some(x) => {
                fw = ((fw << 2) | x) & mask;
                rv = (rv >> 2) | ((3 - x) << shift);
                valid += 1;
                if valid >= K {
                    f(fw.min(rv));
                }
            }
            None => valid = 0,
        }
    }
}

fn read_fasta_region(fa: &mut noodles::fasta::io::IndexedReader<noodles::fasta::io::BufReader<std::fs::File>>, chrom: &str, s: u64, e: u64) -> io::Result<Vec<u8>> {
    let region: noodles::core::Region = format!("{chrom}:{}-{e}", s + 1).parse().map_err(io::Error::other)?;
    Ok(fa.query(&region)?.sequence().as_ref().iter().map(|b| b.to_ascii_uppercase()).collect())
}

/// gene regions +/- 2 kb and the exon CNV control windows, merged per chromosome
pub fn targets(reference: &str, regions: &[(String, u64, u64, String)], lengths: &[(String, u64)]) -> io::Result<Vec<Target>> {
    let mut iv: Vec<(String, u64, u64)> = regions.iter().map(|(c, s, e, _)| (c.clone(), s.saturating_sub(FLANK), e + FLANK)).collect();
    for (_, c, s, e) in exon_cnv::panel().controls {
        iv.push((c, s.saturating_sub(500), e + 500));
    }
    iv.sort();
    let mut merged: Vec<(String, u64, u64)> = vec![];
    for (c, s, e) in iv {
        let len = lengths.iter().find(|x| x.0 == c).map(|x| x.1).unwrap_or(e);
        let e = e.min(len);
        match merged.last_mut() {
            Some(m) if m.0 == c && s <= m.2 => m.2 = m.2.max(e),
            _ => merged.push((c, s, e)),
        }
    }
    let mut fa = noodles::fasta::io::indexed_reader::Builder::default().build_from_path(reference)?;
    merged
        .into_iter()
        .map(|(c, s, e)| Ok(Target { name: format!("{c}:{}-{e}", s + 1), seq: read_fasta_region(&mut fa, &c, s, e)?, chrom: c, start: s }))
        .collect()
}

/// canonical 31-mers with their start offset
fn kmers_at(s: &[u8], mut f: impl FnMut(usize, u64)) {
    let mask = (1u64 << (2 * K)) - 1;
    let shift = 2 * (K - 1);
    let (mut fw, mut rv, mut valid) = (0u64, 0u64, 0usize);
    for (i, &c) in s.iter().enumerate() {
        match enc(c) {
            Some(x) => {
                fw = ((fw << 2) | x) & mask;
                rv = (rv >> 2) | ((3 - x) << shift);
                valid += 1;
                if valid >= K {
                    f(i + 1 - K, fw.min(rv));
                }
            }
            None => valid = 0,
        }
    }
}

/// k-mers seen in at most this many places count as paralog (not repeat) k-mers for the decoys
const DECOY_MAX_COPIES: u32 = 20;
/// a 1 kb window outside the targets with at least this many paralog k-mers becomes a decoy
const DECOY_MIN_KMERS: u32 = 30;

pub struct Scan {
    pub unique: AHashSet<u64>,
    /// loci elsewhere in the genome that share many low-copy 31-mers with the targets
    /// (pseudogenes, segmental duplications): chrom, start, end (0-based)
    pub decoys: Vec<(String, u64, u64)>,
}

/// One pass over the whole reference: the target 31-mers that occur once (the read filter), and
/// the paralogous loci of the targets (decoys). Cached next to the reference (REF.brca12-kmers).
pub fn scan_reference(reference: &str, targets: &[Target]) -> io::Result<Scan> {
    let mut target: AHashMap<u64, u32> = AHashMap::new();
    for t in targets {
        kmers(&t.seq, |k| {
            target.entry(k).or_insert(0);
        });
    }
    // cache: a checksum of the target k-mers identifies the target set
    let sig: u64 = target.keys().fold(0u64, |a, k| a.wrapping_add(k.wrapping_mul(0x9E37_79B9_7F4A_7C15))) ^ 0x2;
    let cache = format!("{reference}.brca12-kmers");
    if let Ok(bytes) = std::fs::read(&cache) {
        let u = |i: usize| u64::from_le_bytes(bytes[i..i + 8].try_into().unwrap());
        if bytes.len() >= 24 && u(0) == sig {
            let n = u(8) as usize;
            let mut unique = AHashSet::with_capacity(n);
            for i in 0..n {
                unique.insert(u(16 + 8 * i));
            }
            let mut at = 16 + 8 * n;
            let nd = u(at) as usize;
            at += 8;
            let mut decoys = vec![];
            for _ in 0..nd {
                let len = u(at) as usize;
                let name = String::from_utf8_lossy(&bytes[at + 8..at + 8 + len]).into_owned();
                at += 8 + len;
                decoys.push((name, u(at), u(at + 8)));
                at += 16;
            }
            eprintln!("brca12-wgs: unique k-mers and decoys from {cache}");
            return Ok(Scan { unique, decoys });
        }
    }
    eprintln!("brca12-wgs: scanning the whole reference for the target 31-mers (once; cached in {cache})");
    let mut inside: AHashMap<String, Vec<(u64, u64)>> = AHashMap::new();
    for t in targets {
        inside.entry(t.chrom.clone()).or_default().push((t.start, t.start + t.seq.len() as u64));
    }
    // per k-mer: occurrences, and up to DECOY_MAX_COPIES + 1 positions outside the targets
    let mut pos: AHashMap<u64, Vec<(u32, u64)>> = AHashMap::new();
    let mut names: Vec<String> = vec![];
    let file = BufReader::with_capacity(1 << 22, std::fs::File::open(reference)?);
    let mut contig: Vec<u8> = vec![];
    let mut flush = |name: &str, contig: &mut Vec<u8>, names: &mut Vec<String>, target: &mut AHashMap<u64, u32>, pos: &mut AHashMap<u64, Vec<(u32, u64)>>| {
        if contig.is_empty() {
            return;
        }
        let cid = names.len() as u32;
        names.push(name.to_string());
        let ins = inside.get(name).cloned().unwrap_or_default();
        let chunk = 4_000_000usize;
        let parts: Vec<(AHashMap<u64, u32>, AHashMap<u64, Vec<(u32, u64)>>)> = (0..contig.len().div_ceil(chunk))
            .into_par_iter()
            .map(|i| {
                let s = i * chunk;
                let e = ((i + 1) * chunk + K - 1).min(contig.len());
                let (mut cnt, mut ps): (AHashMap<u64, u32>, AHashMap<u64, Vec<(u32, u64)>>) = (AHashMap::new(), AHashMap::new());
                kmers_at(&contig[s..e], |o, k| {
                    if target.contains_key(&k) {
                        *cnt.entry(k).or_insert(0) += 1;
                        let p = (s + o) as u64;
                        if !ins.iter().any(|&(a, b)| p >= a && p < b) {
                            let v = ps.entry(k).or_default();
                            if v.len() <= DECOY_MAX_COPIES as usize {
                                v.push((cid, p));
                            }
                        }
                    }
                });
                (cnt, ps)
            })
            .collect();
        for (cnt, ps) in parts {
            for (k, n) in cnt {
                *target.get_mut(&k).unwrap() += n;
            }
            for (k, v) in ps {
                let e = pos.entry(k).or_default();
                for x in v {
                    if e.len() <= DECOY_MAX_COPIES as usize {
                        e.push(x);
                    }
                }
            }
        }
        contig.clear();
    };
    let mut name = String::new();
    for line in file.lines() {
        let line = line?;
        if let Some(h) = line.strip_prefix('>') {
            flush(&name, &mut contig, &mut names, &mut target, &mut pos);
            name = h.split_whitespace().next().unwrap_or("").to_string();
        } else {
            contig.extend(line.trim_end().as_bytes());
        }
    }
    flush(&name, &mut contig, &mut names, &mut target, &mut pos);
    // decoys: 1 kb windows with many low-copy (2..=20) target k-mers, merged and padded by 1 kb
    let mut win: AHashMap<(u32, u64), u32> = AHashMap::new();
    for (k, v) in &pos {
        let n = target[k];
        if (2..=DECOY_MAX_COPIES).contains(&n) {
            for &(c, p) in v {
                *win.entry((c, p / 1000)).or_insert(0) += 1;
            }
        }
    }
    let mut hot: Vec<(u32, u64)> = win.into_iter().filter(|&(_, n)| n >= DECOY_MIN_KMERS).map(|(w, _)| w).collect();
    hot.sort();
    let mut decoys: Vec<(String, u64, u64)> = vec![];
    for (c, w) in hot {
        let (s, e) = ((w * 1000).saturating_sub(1000), (w + 1) * 1000 + 1000);
        let name = &names[c as usize];
        match decoys.last_mut() {
            Some(d) if &d.0 == name && s <= d.2 => d.2 = d.2.max(e),
            _ => decoys.push((name.clone(), s, e)),
        }
    }
    let unique: AHashSet<u64> = target.into_iter().filter(|&(_, n)| n == 1).map(|(k, _)| k).collect();
    let mut out = Vec::with_capacity(32 + 8 * unique.len());
    out.extend(sig.to_le_bytes());
    out.extend((unique.len() as u64).to_le_bytes());
    for k in &unique {
        out.extend(k.to_le_bytes());
    }
    out.extend((decoys.len() as u64).to_le_bytes());
    for (n, s, e) in &decoys {
        out.extend((n.len() as u64).to_le_bytes());
        out.extend(n.as_bytes());
        out.extend(s.to_le_bytes());
        out.extend(e.to_le_bytes());
    }
    if std::fs::write(&cache, out).is_err() {
        eprintln!("brca12-wgs: could not write the k-mer cache {cache} (continuing)");
    }
    Ok(Scan { unique, decoys })
}

struct Fq {
    name: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

fn fq_reader(path: &str) -> io::Result<Box<dyn BufRead + Send>> {
    let f = std::fs::File::open(path).map_err(|e| io::Error::other(format!("{path}: {e}")))?;
    let mut head = [0u8; 2];
    let mut probe = BufReader::new(f);
    let buf = probe.fill_buf()?;
    head.copy_from_slice(&buf[..2.min(buf.len())]);
    Ok(if head == [0x1f, 0x8b] {
        Box::new(BufReader::with_capacity(1 << 22, MultiGzDecoder::new(probe)))
    } else {
        Box::new(probe)
    })
}

fn read_batch(r: &mut dyn BufRead, n: usize) -> io::Result<Vec<Fq>> {
    let mut out = Vec::with_capacity(n);
    let mut l = [vec![], vec![], vec![], vec![]];
    while out.len() < n {
        for x in l.iter_mut() {
            x.clear();
        }
        if r.read_until(b'\n', &mut l[0])? == 0 {
            break;
        }
        for x in l.iter_mut().skip(1) {
            r.read_until(b'\n', x)?;
        }
        let trim = |v: &[u8]| -> Vec<u8> { v.strip_suffix(b"\n").map(|v| v.strip_suffix(b"\r").unwrap_or(v)).unwrap_or(v).to_vec() };
        let mut name = trim(&l[0][1..]);
        if let Some(p) = name.iter().position(|&c| c == b' ' || c == b'\t') {
            name.truncate(p);
        }
        if name.ends_with(b"/1") || name.ends_with(b"/2") {
            name.truncate(name.len() - 2);
        }
        out.push(Fq { name, seq: trim(&l[1]).to_ascii_uppercase(), qual: trim(&l[3]).iter().map(|q| q.saturating_sub(33)).collect() });
    }
    Ok(out)
}

/// one alignment record before it is written
#[derive(Clone)]
struct Rec {
    name: Vec<u8>,
    flag: u16,
    tid: usize,
    pos: u64, // 0-based
    mapq: u8,
    cigar: Vec<(u32, char)>,
    mate_tid: Option<usize>,
    mate_pos: u64,
    tlen: i64,
    seq: Vec<u8>,
    qual: Vec<u8>,
    sa: Option<String>,
    nm: i32,
    end: u64,
}

fn parse_cigar(s: &str) -> Vec<(u32, char)> {
    let mut out = vec![];
    let mut n = 0u32;
    for c in s.chars() {
        if let Some(d) = c.to_digit(10) {
            n = n * 10 + d;
        } else {
            out.push((n, c));
            n = 0;
        }
    }
    out
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            x => x,
        })
        .collect()
}

struct Hit {
    tid: usize,
    pos: u64,
    end: u64,
    rev: bool,
    mapq: u8,
    cigar: Vec<(u32, char)>,
    supp: bool,
    nm: i32,
}

/// One mate's alignment as a record: soft clips, M/I/D CIGAR, position on the chromosome
fn to_hit(h: &aln::Hit, qlen: usize, targets: &[Target], tids: &AHashMap<String, usize>, mapq: u8, supp: bool) -> Hit {
    let t = &targets[h.seq];
    let lead = if h.fwd { h.q_st } else { qlen - h.q_en };
    let mut cigar: Vec<(u32, char)> = vec![];
    let mut push = |n: u32, c: char| {
        if n == 0 {
            return;
        }
        match cigar.last_mut() {
            Some(x) if x.1 == c => x.0 += n,
            _ => cigar.push((n, c)),
        }
    };
    push(lead as u32, 'S');
    for &(op, l) in &h.ops {
        push(l, match op { b'=' | b'X' => 'M', b'I' => 'I', _ => 'D' });
    }
    push((qlen - lead - (h.q_en - h.q_st)) as u32, 'S');
    Hit { tid: tids[&t.chrom], pos: t.start + h.r_st as u64, end: t.start + h.r_en as u64, rev: !h.fwd, mapq, cigar, supp, nm: h.nm as i32 }
}

/// mapping quality in the manner of minimap2 (sr): the best score against the best competing score on
/// the same part of the read, scaled down by the alignment's identity and, as minimap2 does through its
/// anchor count, for reads with few exact 21-mers (copies of a diverged repeat or paralog that is not
/// in the panel; minimap2 samples one minimizer per ~6 k-mers, hence the 1/6)
fn mapq_of(s1: i32, s2: i32, h: &aln::Hit) -> u8 {
    if s1 <= 0 {
        return 0;
    }
    if *OLD_MAPQ {
        if s2 <= 0 {
            return 60;
        }
        let x = 40.0 * (1.0 - s2 as f64 / s1 as f64) * (h.mlen as f64 / 10.0).min(1.0) * (s1 as f64).ln();
        return x.clamp(0.0, 60.0).round() as u8;
    }
    let cols: u32 = h.ops.iter().map(|&(_, l)| l).sum();
    let ident = h.mlen as f64 / cols.max(1) as f64;
    let kmers: u32 = h.ops.iter().filter(|&&(o, _)| o == b'=').map(|&(_, l)| l.saturating_sub(20)).sum();
    let cnt = kmers as f64 / 6.0;
    let pen = (0.1 * cnt).min(1.0).min(s1 as f64 / 100.0);
    let s2 = s2.max(40).min(s1);
    let x = ident * pen * 40.0 * (1.0 - s2 as f64 / s1 as f64) * (s1 as f64).ln();
    let q = x.clamp(0.0, 60.0) as u8;
    if q == 0 && s1 > s2 { 1 } else { q }
}

static OLD_MAPQ: std::sync::LazyLock<bool> = std::sync::LazyLock::new(|| std::env::var("BRCA12_MAPQ").is_ok_and(|v| v == "old"));
static OLD_RESCUE: std::sync::LazyLock<bool> = std::sync::LazyLock::new(|| std::env::var("BRCA12_RESCUE").is_ok_and(|v| v == "old"));
static NO_RESCUE: std::sync::LazyLock<bool> = std::sync::LazyLock::new(|| std::env::var("BRCA12_RESCUE").is_ok_and(|v| v == "0"));

/// alignments of a read pair with the own aligner (from hla-wgs): the best proper pair on one target
/// earns a bonus as in paired mapping (hla-wgs panel: 33); a mate without a hit next to its partner is
/// aligned in a window beside it (mate rescue, as in hla-wgs). The primary hit on a decoy means the
/// read belongs elsewhere in the genome: no record for that mate.
fn pair_hits(idx: &aln::Index, n_targets: usize, x: &Fq, y: &Fq, targets: &[Target], tids: &AHashMap<String, usize>) -> (Vec<Hit>, Vec<Hit>) {
    const BONUS: i32 = 33;
    const RESCUE_MIN: i32 = 40;
    let mut m = [idx.map(&x.seq), idx.map(&y.seq)];
    let qs = [&x.seq, &y.seq];
    // mate rescue next to the partner's best hit
    for mate in 0..2 {
        if *NO_RESCUE {
            break;
        }
        let other = 1 - mate;
        let Some(p) = m[other].first().cloned() else { continue };
        let near = |h: &aln::Hit| h.seq == p.seq && h.fwd != p.fwd && (h.r_st as i64 - p.r_st as i64).abs() <= 1500;
        if m[mate].iter().any(near) {
            continue;
        }
        let seq = &idx.seqs[p.seq];
        let (w0, w1) = (p.r_st.saturating_sub(1000), (p.r_en + 1000).min(seq.len()));
        let qo: Vec<u8> = if p.fwd { aln::revcomp(qs[mate]) } else { qs[mate].to_ascii_uppercase() };
        let sub = aln::Index::new(vec![String::new()], vec![seq[w0..w1].to_vec()], aln::Opts::sr_single());
        // a rescued mate must align as a whole and well: a partial piece of a mate that belongs to
        // another copy elsewhere would make the pair look aligned and keep the partner's MAPQ up
        let ql = qo.len();
        let good = |h: &aln::Hit| if *OLD_RESCUE { h.score >= RESCUE_MIN } else { (h.q_en - h.q_st) * 10 >= ql * 9 && h.score as usize >= ql };
        if let Some(mut h) = sub.map(&qo).into_iter().find(|h| h.fwd && good(h)) {
            // back to the read's own orientation
            h.seq = p.seq;
            h.r_st += w0;
            h.r_en += w0;
            h.fwd = !p.fwd;
            if p.fwd {
                let l = qo.len();
                let (a, b) = (l - h.q_en, l - h.q_st);
                h.q_st = a;
                h.q_en = b;
            }
            m[mate].push(h);
            m[mate].sort_by(|a, b| b.score.cmp(&a.score));
        }
    }
    // primary per mate: best single hits, or the best proper pair if it scores higher with the bonus
    let proper = |a: &aln::Hit, b: &aln::Hit| {
        a.seq == b.seq && a.fwd != b.fwd && {
            let (l, r) = if a.r_st <= b.r_st { (a, b) } else { (b, a) };
            l.fwd && !r.fwd && r.r_en.max(l.r_en) - l.r_st <= 2000
        }
    };
    let top = |k: usize| m[k].first().map_or(0, |h| h.score);
    let mut pairs: Vec<(i32, usize, usize)> = vec![];
    for (i, a) in m[0].iter().enumerate().take(20) {
        for (j, b) in m[1].iter().enumerate().take(20) {
            if proper(a, b) {
                pairs.push((a.score + b.score + BONUS, i, j));
            }
        }
    }
    pairs.sort_by(|a, b| b.0.cmp(&a.0));
    let paired = pairs.first().filter(|p| p.0 >= top(0) + top(1)).copied();
    let pick = match paired {
        Some((_, i, j)) => [Some(i), Some(j)],
        None => [(!m[0].is_empty()).then_some(0), (!m[1].is_empty()).then_some(0)],
    };
    // the pair's MAPQ as minimap2 (pe.c): the better mate's MAPQ, limited by how much the best pair
    // beats the next ones; a mate below it is raised to 0.2 own + 0.8 pair
    let single = |k: usize, pi: usize| {
        let p = &m[k][pi];
        let ov = |a: &aln::Hit, b: &aln::Hit| (a.q_en.min(b.q_en) as i64 - a.q_st.max(b.q_st) as i64).max(0) as usize;
        let s2 = m[k].iter().enumerate().filter(|(i, h)| *i != pi && ov(h, p) * 2 >= (p.q_en - p.q_st)).map(|(_, h)| h.score).max().unwrap_or(0);
        mapq_of(p.score, s2, p)
    };
    let mut qpe: [Option<u8>; 2] = [None, None];
    if let Some((ps, i, j)) = paired {
        if *OLD_MAPQ {
            let p2 = pairs.iter().skip(1).map(|x| x.0).next().unwrap_or(0);
            let v = |h: &aln::Hit| mapq_of(ps, p2, h);
            qpe = [Some(v(&m[0][i])), Some(v(&m[1][j]))];
        } else {
            let (q0, q1) = (single(0, i), single(1, j));
            let mut pe = q0.max(q1) as f64;
            if pairs.len() > 1 {
                let n_sub = pairs.iter().filter(|x| x.0 + 2 * 2 >= ps).count().max(1) as f64;
                let alt = 6.02 * (ps - pairs[1].0) as f64 / 2.0 - 4.343 * n_sub.ln();
                pe = pe.min(alt);
            }
            let lift = |q: u8| {
                let q = if (q as f64) < pe { (0.2 * q as f64 + 0.8 * pe + 0.499) as u8 } else { q };
                let floor = if pairs.len() == 1 { 2 } else if ps > pairs[1].0 { 1 } else { 0 };
                q.max(floor).min(60)
            };
            qpe = [Some(lift(q0)), Some(lift(q1))];
        }
    }
    let mut out: [Vec<Hit>; 2] = [vec![], vec![]];
    for k in 0..2 {
        let Some(pi) = pick[k] else { continue };
        let p = &m[k][pi];
        if p.seq >= n_targets {
            continue;
        }
        let qlen = qs[k].len();
        let ov = |a: &aln::Hit, b: &aln::Hit| (a.q_en.min(b.q_en) as i64 - a.q_st.max(b.q_st) as i64).max(0) as usize;
        let s2 = m[k].iter().enumerate().filter(|(i, h)| *i != pi && ov(h, p) * 2 >= (p.q_en - p.q_st)).map(|(_, h)| h.score).max().unwrap_or(0);
        let mut q = mapq_of(p.score, s2, p);
        if paired.is_some() {
            q = qpe[k].unwrap_or(q).max(q);
        }
        out[k].push(to_hit(p, qlen, targets, tids, q, false));
        // supplementary: hits on a part of the read the primary does not cover (split reads)
        let mut covered = vec![(p.q_st, p.q_en)];
        for (i, h) in m[k].iter().enumerate() {
            if i == pi || h.seq >= n_targets {
                continue;
            }
            let len = (h.q_en - h.q_st).max(1);
            let o = covered.iter().map(|&(a, b)| (h.q_en.min(b) as i64 - h.q_st.max(a) as i64).max(0) as usize).max().unwrap_or(0);
            if o * 2 >= len || len < 30 {
                continue;
            }
            // a piece that fits as well elsewhere (repeat copy, Alu) is no evidence of a split read
            let s2 = m[k].iter().enumerate().filter(|(j, g)| *j != i && ov(g, h) * 2 >= len).map(|(_, g)| g.score).max().unwrap_or(0);
            let sq = mapq_of(h.score, s2, h);
            if sq < 20 {
                continue;
            }
            covered.push((h.q_st, h.q_en));
            out[k].push(to_hit(h, qlen, targets, tids, q.min(sq), true));
        }
    }
    let [a, b] = out;
    (a, b)
}

/// SAM records of one pair (primary, supplementary, unmapped mate placed at its partner)
fn pair_records(name: &[u8], r1: &Fq, r2: &Fq, h1: &[Hit], h2: &[Hit], names: &[String]) -> Vec<Rec> {
    let p1 = h1.iter().find(|h| !h.supp);
    let p2 = h2.iter().find(|h| !h.supp);
    if p1.is_none() && p2.is_none() {
        return vec![];
    }
    let proper = match (p1, p2) {
        (Some(a), Some(b)) => {
            a.tid == b.tid && a.rev != b.rev && {
                let (l, r) = if a.pos <= b.pos { (a, b) } else { (b, a) };
                !l.rev && r.rev && r.end.max(l.end) - l.pos <= 2000
            }
        }
        _ => false,
    };
    let tlen = |a: &Hit, b: &Hit| -> i64 {
        if a.tid != b.tid {
            return 0;
        }
        let (s, e) = (a.pos.min(b.pos), a.end.max(b.end));
        let v = (e - s) as i64;
        if a.pos < b.pos || (a.pos == b.pos && !a.rev) { v } else { -v }
    };
    let sa_of = |hs: &[Hit], skip: usize| -> Option<String> {
        let other: Vec<String> = hs
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != skip)
            .map(|(_, h)| {
                let cig: String = h.cigar.iter().map(|(n, c)| format!("{n}{c}")).collect();
                format!("{},{},{},{},{},{};", names[h.tid], h.pos + 1, if h.rev { '-' } else { '+' }, cig, h.mapq, h.nm)
            })
            .collect();
        if other.is_empty() { None } else { Some(other.concat()) }
    };
    let mut out = vec![];
    for (mate, (r, hs, own, other)) in [(r1, h1, p1, p2), (r2, h2, p2, p1)].into_iter().enumerate() {
        let first: u16 = if mate == 0 { 0x40 } else { 0x80 };
        match own {
            None => {
                // unmapped read placed at its mapped partner
                let o = other.unwrap();
                let mut flag = 0x1 | 0x4 | first;
                if o.rev {
                    flag |= 0x20;
                }
                out.push(Rec {
                    name: name.to_vec(), flag, tid: o.tid, pos: o.pos, mapq: 0, cigar: vec![], mate_tid: Some(o.tid), mate_pos: o.pos,
                    tlen: 0, seq: r.seq.clone(), qual: r.qual.clone(), sa: None, nm: 0, end: o.pos + 1,
                });
            }
            Some(_) => {
                for (i, h) in hs.iter().enumerate() {
                    let mut flag = 0x1 | first;
                    if proper {
                        flag |= 0x2;
                    }
                    if h.rev {
                        flag |= 0x10;
                    }
                    if h.supp {
                        flag |= 0x800;
                    }
                    let (mate_tid, mate_pos, tl) = match other {
                        Some(o) => {
                            if o.rev {
                                flag |= 0x20;
                            }
                            (Some(o.tid), o.pos, if h.supp { 0 } else { tlen(h, o) })
                        }
                        None => {
                            flag |= 0x8;
                            (Some(h.tid), h.pos, 0)
                        }
                    };
                    let (seq, qual) = if h.rev { (revcomp(&r.seq), r.qual.iter().rev().cloned().collect()) } else { (r.seq.clone(), r.qual.clone()) };
                    out.push(Rec {
                        name: name.to_vec(), flag, tid: h.tid, pos: h.pos, mapq: h.mapq, cigar: h.cigar.clone(), mate_tid, mate_pos,
                        tlen: tl, seq, qual, sa: sa_of(hs, i), nm: h.nm, end: h.end,
                    });
                }
            }
        }
    }
    out
}

/// duplicates of pairs as samtools markdup: same chromosome, unclipped 5' positions and
/// orientations of both reads; the pair with the highest sum of base qualities >= 15 is kept
fn mark_duplicates(recs: &mut [Rec]) {
    let five = |r: &Rec| -> (usize, i64, bool) {
        let lead = |c: &[(u32, char)]| c.iter().take_while(|x| x.1 == 'S' || x.1 == 'H').map(|x| x.0 as i64).sum::<i64>();
        if r.flag & 0x10 != 0 {
            let rev: Vec<(u32, char)> = r.cigar.iter().rev().cloned().collect();
            (r.tid, r.end as i64 + lead(&rev), true)
        } else {
            (r.tid, r.pos as i64 - lead(&r.cigar), false)
        }
    };
    let mut by_name: AHashMap<Vec<u8>, Vec<usize>> = AHashMap::new();
    for (i, r) in recs.iter().enumerate() {
        if r.flag & 0x900 == 0 {
            by_name.entry(r.name.clone()).or_default().push(i);
        }
    }
    let mut groups: AHashMap<(usize, i64, bool, usize, i64, bool), Vec<(u64, Vec<u8>)>> = AHashMap::new();
    for (name, idx) in &by_name {
        if idx.len() != 2 || idx.iter().any(|&i| recs[i].flag & 0x4 != 0) {
            continue;
        }
        let (a, b) = (five(&recs[idx[0]]), five(&recs[idx[1]]));
        let key = if (a.0, a.1) <= (b.0, b.1) { (a.0, a.1, a.2, b.0, b.1, b.2) } else { (b.0, b.1, b.2, a.0, a.1, a.2) };
        let score: u64 = idx.iter().map(|&i| recs[i].qual.iter().filter(|&&q| q >= 15).map(|&q| q as u64).sum::<u64>()).sum();
        groups.entry(key).or_default().push((score, name.clone()));
    }
    let mut dups: AHashSet<Vec<u8>> = AHashSet::new();
    for (_, mut v) in groups {
        if v.len() > 1 {
            v.sort_by(|x, y| y.0.cmp(&x.0).then(x.1.cmp(&y.1)));
            dups.extend(v.into_iter().skip(1).map(|x| x.1));
        }
    }
    for r in recs.iter_mut() {
        if dups.contains(&r.name) {
            r.flag |= 0x400;
        }
    }
}

fn write_bam(path: &str, contigs: &[(String, u64)], recs: &[Rec]) -> io::Result<()> {
    use noodles::sam::alignment::record::Flags;
    use noodles::sam::alignment::record::cigar::op::{Kind, Op};
    use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence, data::field::Value};
    let mut text = format!("@HD\tVN:1.6\tSO:coordinate\n");
    for (n, l) in contigs {
        text.push_str(&format!("@SQ\tSN:{n}\tLN:{l}\n"));
    }
    text.push_str(&format!("@PG\tID:brca12-wgs\tPN:brca12-wgs\tVN:{}\n", env!("CARGO_PKG_VERSION")));
    let header: sam::Header = text.parse().map_err(io::Error::other)?;
    use noodles::sam::alignment::io::Write as _;
    let mut w = noodles::bam::io::Writer::new(std::fs::File::create(path)?);
    w.write_header(&header)?;
    for r in recs {
        let kind = |c: char| match c {
            'M' => Kind::Match,
            'I' => Kind::Insertion,
            'D' => Kind::Deletion,
            'N' => Kind::Skip,
            'S' => Kind::SoftClip,
            'H' => Kind::HardClip,
            'P' => Kind::Pad,
            '=' => Kind::SequenceMatch,
            _ => Kind::SequenceMismatch,
        };
        let cigar: Cigar = r.cigar.iter().map(|&(n, c)| Op::new(kind(c), n as usize)).collect();
        let mut data = noodles::sam::alignment::record_buf::Data::default();
        use noodles::sam::alignment::record::data::field::Tag;
        if r.flag & 0x4 == 0 {
            data.insert(Tag::EDIT_DISTANCE, Value::from(r.nm));
        }
        if let Some(sa) = &r.sa {
            data.insert(Tag::OTHER_ALIGNMENTS, Value::String(sa.as_bytes().into()));
        }
        let mut b = noodles::sam::alignment::RecordBuf::builder()
            .set_name(r.name.clone())
            .set_flags(Flags::from(r.flag))
            .set_reference_sequence_id(r.tid)
            .set_alignment_start(noodles::core::Position::try_from(r.pos as usize + 1).map_err(io::Error::other)?)
            .set_cigar(cigar)
            .set_template_length(r.tlen as i32)
            .set_sequence(Sequence::from(r.seq.clone()))
            .set_quality_scores(QualityScores::from(r.qual.clone()))
            .set_data(data);
        if let Some(q) = noodles::sam::alignment::record::MappingQuality::new(r.mapq) {
            b = b.set_mapping_quality(q);
        }
        if let Some(t) = r.mate_tid {
            b = b.set_mate_reference_sequence_id(t)
                .set_mate_alignment_start(noodles::core::Position::try_from(r.mate_pos as usize + 1).map_err(io::Error::other)?);
        }
        w.write_alignment_record(&header, &b.build())?;
    }
    w.try_finish()?;
    drop(w);
    match noodles::bam::fs::index(path)? {
        noodles::bam::Index::Bai(index) => noodles::bam::bai::fs::write(format!("{path}.bai"), &index)?,
        noodles::bam::Index::Csi(_) => return Err(io::Error::other("BAM index: expected BAI")),
    }
    Ok(())
}

/// extract, align and write OUTPREFIX.bam from paired FASTQ; returns the BAM path
pub fn to_bam(r1: &str, r2: &str, reference: &str, regions: &[(String, u64, u64, String)], prefix: &str, threads: usize) -> io::Result<String> {
    let t0 = std::time::Instant::now();
    let contigs: Vec<(String, u64)> = std::fs::read_to_string(format!("{reference}.fai"))
        .map_err(|e| io::Error::other(format!("{reference}.fai: {e}")))?
        .lines()
        .filter_map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            Some((f.first()?.to_string(), f.get(1)?.parse().ok()?))
        })
        .collect();
    let tids: AHashMap<String, usize> = contigs.iter().enumerate().map(|(i, c)| (c.0.clone(), i)).collect();
    let names: Vec<String> = contigs.iter().map(|c| c.0.clone()).collect();
    let tg = targets(reference, regions, &contigs)?;
    let scan = scan_reference(reference, &tg)?;
    let uniq = scan.unique;
    // decoys (paralogs elsewhere) join the alignment panel, so that their reads find their own locus
    let decoy_seqs: Vec<(String, Vec<u8>)> = {
        let mut fa = noodles::fasta::io::indexed_reader::Builder::default().build_from_path(reference)?;
        scan.decoys
            .iter()
            .map(|(c, s, e)| {
                let len = contigs.iter().find(|x| &x.0 == c).map(|x| x.1).unwrap_or(*e);
                Ok((format!("decoy|{c}:{}-{}", s + 1, (*e).min(len)), read_fasta_region(&mut fa, c, *s, (*e).min(len))?))
            })
            .collect::<io::Result<_>>()?
    };
    let decoy_bp: usize = decoy_seqs.iter().map(|d| d.1.len()).sum();
    eprintln!("brca12-wgs: {} targets, {} unique 31-mers, {} decoy loci ({} kb) ({:.0} s)", tg.len(), uniq.len(), decoy_seqs.len(), decoy_bp / 1000, t0.elapsed().as_secs_f64());
    // own aligner (from hla-wgs): targets first, then the decoys
    let n_targets = tg.len();
    let all_targets: Vec<Target> = tg
        .iter()
        .map(|t| Target { name: t.name.clone(), seq: vec![], chrom: t.chrom.clone(), start: t.start })
        .collect();
    let idx = aln::Index::new(
        tg.iter().map(|t| t.name.clone()).chain(decoy_seqs.iter().map(|d| d.0.clone())).collect(),
        tg.iter().map(|t| t.seq.clone()).chain(decoy_seqs.iter().map(|d| d.1.clone())).collect(),
        aln::Opts::sr(),
    );
    let _ = threads;
    // two reader threads (R1, R2) feed batches; k-mer screening and alignment run in parallel
    let (tx, rx) = std::sync::mpsc::sync_channel::<(Vec<Fq>, Vec<Fq>)>(4);
    let (p1, p2) = (r1.to_string(), r2.to_string());
    let reader = std::thread::spawn(move || -> io::Result<u64> {
        let (mut a, mut b) = (fq_reader(&p1)?, fq_reader(&p2)?);
        let mut n = 0u64;
        loop {
            let (b1, b2) = std::thread::scope(|s| {
                let h = s.spawn(|| read_batch(&mut *b, 200_000));
                let x = read_batch(&mut *a, 200_000);
                (x, h.join().unwrap())
            });
            let (b1, b2) = (b1?, b2?);
            if b1.is_empty() {
                break;
            }
            if b1.len() != b2.len() {
                return Err(io::Error::other("R1 and R2 have a different number of reads"));
            }
            n += b1.len() as u64;
            if tx.send((b1, b2)).is_err() {
                break;
            }
        }
        Ok(n)
    });
    let mut recs: Vec<Rec> = vec![];
    let mut kept = 0u64;
    for (b1, b2) in rx {
        let has = |s: &[u8]| {
            let mut n = 0;
            kmers(s, |k| {
                if uniq.contains(&k) {
                    n += 1;
                }
            });
            n >= 2
        };
        let got: Vec<Vec<Rec>> = b1
            .par_iter()
            .zip(b2.par_iter())
            .filter(|(x, y)| has(&x.seq) || has(&y.seq))
            .map(|(x, y)| {
                if x.name != y.name {
                    return vec![];
                }
                let (mut h1, mut h2) = pair_hits(&idx, n_targets, x, y, &all_targets, &tids);
                // a read without genome-unique 31-mers of its own could come from a copy outside the
                // targets (segmental duplication, repeat): mapping quality 0, as BWA gives it on the
                // whole genome
                // likewise a read whose mate does not align to the targets at all: the pair comes from
                // another locus (the targets have 2 kb flanks, so mates of gene reads lie inside them)
                let (u1, u2) = (h1.iter().all(|z| z.supp), h2.iter().all(|z| z.supp));
                for (h, r, mate_unaligned) in [(&mut h1, x, u2), (&mut h2, y, u1)] {
                    if !has(&r.seq) || mate_unaligned {
                        h.iter_mut().for_each(|z| z.mapq = 0);
                    }
                }
                pair_records(&x.name, x, y, &h1, &h2, &names)
            })
            .collect();
        for g in got {
            if !g.is_empty() {
                kept += 1;
                recs.extend(g);
            }
        }
    }
    let total = reader.join().map_err(|_| io::Error::other("FASTQ reader failed"))??;
    eprintln!("brca12-wgs: {total} read pairs read, {kept} aligned to BRCA1/BRCA2 and the control windows ({:.0} s)", t0.elapsed().as_secs_f64());
    if kept == 0 {
        return Err(io::Error::other("no read pair of BRCA1/BRCA2 in the FASTQ files (wrong files or reference?)"));
    }
    mark_duplicates(&mut recs);
    recs.sort_by(|a, b| (a.tid, a.pos, a.flag & 0x4, &a.name).cmp(&(b.tid, b.pos, b.flag & 0x4, &b.name)));
    let path = format!("{prefix}.bam");
    write_bam(&path, &contigs, &recs)?;
    std::io::stderr().flush().ok();
    Ok(path)
}
