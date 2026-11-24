//! Small-variant caller for the BRCA1/BRCA2 regions (SNVs and indels), built like the core of GATK
//! HaplotypeCaller without de-novo assembly. Research use only.
//!
//! 1. Reads: MAPQ >= 20; secondary, duplicate and QC-fail records skipped. Clipped as GATK and 1.1:
//!    bases past the fragment end (adapter), low-quality ends (< Q10), poly-G read-through tails.
//!    Where two mates of a fragment overlap and disagree, both bases get quality 0.
//! 2. Candidates from the alignments: mismatches with base quality >= 10, insertions and deletions
//!    (left-normalised), kept with >= 2 supporting fragments; candidates closer than 40 bp form a
//!    cluster (at most the 5 best supported), and every combination of non-overlapping candidates
//!    of a cluster is a haplotype.
//! 3. Every read over a cluster is aligned to every haplotype with a pair-HMM (global in the read,
//!    local in the haplotype): base qualities as error rates (below Q18 counted as Q6, capped at
//!    MAPQ, as GATK), gap open from GATK's CONSERVATIVE PCR indel model (Q40, lower in repeats of the
//!    read), gap continuation Q10; soft clips of 20 or more bases are uninformative.
//! 4. An allele supported only by reads of one strand while reads of the other strand cover it
//!    against it (Fisher > Q20) is an artifact (homopolymer phasing errors); its haplotypes go.
//! 5. Overlapping alleles form a site; per site the read likelihoods are marginalised over the
//!    haplotypes (best haplotype with the allele, without it), as GATK, and the diploid genotypes are
//!    scored with GATK's priors (heterozygosity 0.001 SNV, 1.25e-4 indel). QUAL = phred
//!    P(0/0 | reads); called when QUAL >= 30 and the best genotype is not 0/0; two different alleles
//!    at one site give one record with genotype 1/2.
//! 6. Filters as GATK hard filters: QD < 2, FS > 60 (SNV) / > 200 (indel); SOR only reported.
use crate::assemble;
use crate::reads::{self, Op, Read};
use std::collections::HashMap;
use std::io;

#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Allele {
    pub pos: u64, // 0-based position of the first REF base (VCF POS - 1)
    pub r#ref: Vec<u8>,
    pub alt: Vec<u8>,
}

#[derive(Clone, Debug)]
pub struct Call {
    pub chrom: String,
    pub allele: Allele,
    pub alt2: Option<Vec<u8>>, // second ALT of a multi-allelic record (genotype 1/2)
    pub filter: String,
    pub fs: f64,
    pub sor: f64,
    pub qual: f64,
    pub gt: (u8, u8),
    pub gq: u32,
    pub ad: (u32, u32),
    pub dp: u32,
    pub pl: [u32; 3],
}

pub struct Params {
    pub min_mapq: u8,
    pub min_base_q: u8,
    pub min_support: u32,
    pub call_qual: f64,
    pub het_snv: f64,
    pub het_indel: f64,
    pub flank: usize,
    pub low_base_q: u8,
    pub long_clip: u32,
    pub assemble: bool,          // local de-novo assembly of every cluster window
    pub max_haplotypes: usize,   // haplotypes per cluster at most (combinations first)          // soft clips of at least this length are uninformative          // base qualities below this count as Q6 in the pair-HMM (GATK: 18)
    pub strand_support_min: u32, // one-strand allele: at least this many supporting reads,
    pub strand_other_min: u32,   // at least this many reads of the other strand against it,
    pub strand_fs: f64,          // and Fisher strand test above this (phred)
    pub strand_fs_two_alleles: f64, // the same test for one allele of a 1/2 genotype (p < 0.05)
    pub stutter_tract: usize,    // indel sites in repeat tracts of at least this length (bp)
    pub stutter_cap: f64,        // get per-read evidence capped at this (log10)
}

pub const DEFAULT: Params =
    Params { min_mapq: 20, min_base_q: 10, min_support: 2, call_qual: 30.0, het_snv: 1e-3, het_indel: 1.25e-4, flank: 150,
             low_base_q: 18, long_clip: 20, assemble: true, max_haplotypes: 64, strand_support_min: 4, strand_other_min: 8, strand_fs: 20.0, strand_fs_two_alleles: 13.0, stutter_tract: 10, stutter_cap: 2.0 };

fn usable(r: &Read, p: &Params) -> bool {
    r.flags & 0x704 == 0 && r.mapq >= p.min_mapq && !r.seq.is_empty() && !r.cigar.is_empty()
}

/// left-normalise an indel (VCF style, with anchor base) against the reference `rf` (0-based at `off`)
fn normalise(mut a: Allele, rf: &[u8], off: u64) -> Allele {
    loop {
        // trim a common last base
        if a.r#ref.len() > 1 && a.alt.len() > 1 && a.r#ref.last() == a.alt.last() {
            a.r#ref.pop();
            a.alt.pop();
            continue;
        }
        // shift left when the last bases are equal and one allele would become empty
        if a.r#ref.last() == a.alt.last() && (a.r#ref.len() == 1 || a.alt.len() == 1) && a.pos > off {
            let prev = rf[(a.pos - 1 - off) as usize];
            a.pos -= 1;
            a.r#ref.pop();
            a.alt.pop();
            a.r#ref.insert(0, prev);
            a.alt.insert(0, prev);
            continue;
        }
        break;
    }
    // trim common leading bases beyond the anchor
    while a.r#ref.len() > 1 && a.alt.len() > 1 && a.r#ref[0] == a.alt[0] && a.r#ref[1] == a.alt[1] {
        a.r#ref.remove(0);
        a.alt.remove(0);
        a.pos += 1;
    }
    a
}

/// candidate alleles with the number of supporting fragments
fn candidates(reads: &[Read], rf: &[u8], off: u64, p: &Params) -> HashMap<Allele, (u32, Vec<Vec<u8>>)> {
    let mut out: HashMap<Allele, (u32, Vec<Vec<u8>>)> = HashMap::new();
    for r in reads {
        let mut events: Vec<Allele> = vec![];
        let (mut rp, mut qp) = (r.start, 0usize);
        for &(op, len) in &r.cigar {
            let len = len as usize;
            match op {
                Op::Match => {
                    for k in 0..len {
                        let ri = rp + k as u64;
                        if ri < off || (ri - off) as usize >= rf.len() || qp + k >= r.seq.len() {
                            continue;
                        }
                        let (rb, qb) = (rf[(ri - off) as usize].to_ascii_uppercase(), r.seq[qp + k].to_ascii_uppercase());
                        if rb != qb && qb != b'N' && rb != b'N' && r.qual.get(qp + k).copied().unwrap_or(0) >= p.min_base_q {
                            events.push(Allele { pos: ri, r#ref: vec![rb], alt: vec![qb] });
                        }
                    }
                    rp += len as u64;
                    qp += len;
                }
                Op::Ins => {
                    if rp > off && ((rp - off) as usize) <= rf.len() && qp > 0 {
                        let anchor = rf[(rp - 1 - off) as usize].to_ascii_uppercase();
                        let mut alt = vec![anchor];
                        alt.extend(r.seq[qp..(qp + len).min(r.seq.len())].iter().map(|b| b.to_ascii_uppercase()));
                        events.push(normalise(Allele { pos: rp - 1, r#ref: vec![anchor], alt }, rf, off));
                    }
                    qp += len;
                }
                Op::Del => {
                    if rp > off && ((rp - off) as usize + len) <= rf.len() {
                        let s = (rp - 1 - off) as usize;
                        let refb: Vec<u8> = rf[s..s + 1 + len].iter().map(|b| b.to_ascii_uppercase()).collect();
                        events.push(normalise(Allele { pos: rp - 1, r#ref: refb.clone(), alt: vec![refb[0]] }, rf, off));
                    }
                    rp += len as u64;
                }
                Op::Skip => rp += len as u64,
                Op::SoftClip => qp += len,
                _ => {}
            }
        }
        for e in events {
            let ent = out.entry(e).or_insert((0, vec![]));
            if !ent.1.iter().any(|n| n == &r.name) {
                ent.1.push(r.name.clone());
                ent.0 += 1;
            }
        }
    }
    out.retain(|_, v| v.0 >= p.min_support);
    out
}

/// Read clean-up before calling, as GATK HaplotypeCaller (and brca12-wgs 1.1) do: bases past the
/// fragment end (adapter read-through, from the insert size), low-quality ends (base quality < 10)
/// and soft-clipped poly-G read-through tails of two-colour chemistry (1.1: brca12_clip_polyg.py)
/// are removed from the read.
fn clip_read(r: &mut Read) {
    let n = r.seq.len();
    if n == 0 {
        return;
    }
    let (mut a, mut b) = (0usize, n);
    // query index of a reference position (first query base at or after it)
    let qpos_of = |r: &Read, refpos: u64| -> usize {
        let (mut rp, mut qp) = (r.start, 0usize);
        for &(op, len) in &r.cigar {
            let l = len as u64;
            match op {
                Op::Match => {
                    if refpos < rp + l {
                        return qp + refpos.saturating_sub(rp) as usize;
                    }
                    rp += l;
                    qp += len as usize;
                }
                Op::Ins | Op::SoftClip => qp += len as usize,
                Op::Del | Op::Skip => {
                    if refpos < rp + l {
                        return qp;
                    }
                    rp += l;
                }
                _ => {}
            }
        }
        qp
    };
    // adapter: forward read past mate start + |insert|, reverse read before its mate's start
    if r.flags & 0x1 != 0 && r.flags & 0x8 == 0 && r.mate_tid == r.tid && r.tlen != 0 && r.mate_start >= 0 {
        if r.flags & 0x10 == 0 {
            let frag_end = r.start + r.tlen.unsigned_abs();
            if frag_end < r.end {
                b = b.min(qpos_of(r, frag_end));
            }
        } else {
            let ms = r.mate_start as u64;
            if ms > r.start {
                a = a.max(qpos_of(r, ms));
            }
        }
    }
    // poly-G (C on reverse reads) soft-clipped tail at the sequencing 3' end
    let rev = r.flags & 0x10 != 0;
    if !rev {
        if let Some(&(Op::SoftClip, k)) = r.cigar.last() {
            let k = k as usize;
            let tail = &r.seq[n - k..(n - k + 10).min(n)];
            if k >= 5 && tail.iter().filter(|&&x| x == b'G' || x == b'g').count() as f64 >= 0.8 * tail.len() as f64 {
                b = b.min(n - k);
            }
        }
    } else if let Some(&(Op::SoftClip, k)) = r.cigar.first() {
        let k = k as usize;
        let head = &r.seq[k.saturating_sub(10)..k];
        if k >= 5 && head.iter().filter(|&&x| x == b'C' || x == b'c').count() as f64 >= 0.8 * head.len() as f64 {
            a = a.max(k);
        }
    }
    // low-quality ends
    while a < b && r.qual.get(a).is_some_and(|&q| q < 10) {
        a += 1;
    }
    while b > a && r.qual.get(b - 1).is_some_and(|&q| q < 10) {
        b -= 1;
    }
    if a == 0 && b == n {
        return;
    }
    if a >= b {
        r.seq.clear();
        r.qual.clear();
        r.cigar.clear();
        return;
    }
    // new CIGAR over query bases a..b; the start moves by the reference bases before a
    let mut cig: Vec<(Op, u32)> = vec![];
    let (mut qp, mut rp) = (0usize, r.start);
    let mut new_start = None;
    for &(op, len) in &r.cigar {
        let l = len as usize;
        let (cq, cr) = match op {
            Op::Match => (true, true),
            Op::Ins | Op::SoftClip => (true, false),
            Op::Del | Op::Skip => (false, true),
            _ => (false, false),
        };
        if cq {
            let s0 = qp.max(a);
            let e0 = (qp + l).min(b);
            if e0 > s0 {
                if new_start.is_none() {
                    new_start = Some(rp + if cr { (s0 - qp) as u64 } else { 0 });
                }
                cig.push((op, (e0 - s0) as u32));
            }
            qp += l;
            if cr {
                rp += l as u64;
            }
        } else if cr {
            // a deletion inside the kept part stays, at the edges it goes
            if qp > a && qp < b {
                cig.push((op, len));
            }
            rp += l as u64;
        }
    }
    while matches!(cig.last(), Some((Op::Del | Op::Skip | Op::Ins, _))) {
        cig.pop();
    }
    while matches!(cig.first(), Some((Op::Del | Op::Skip, _))) {
        cig.remove(0);
    }
    let span: u64 = cig.iter().filter(|(o, _)| matches!(o, Op::Match | Op::Del | Op::Skip)).map(|&(_, l)| l as u64).sum();
    r.start = new_start.unwrap_or(r.start);
    r.end = r.start + span;
    r.seq = r.seq[a..b].to_vec();
    r.qual = r.qual[a..b].to_vec();
    r.cigar = cig;
}

/// GATK-style handling of overlapping mates: where both mates of a fragment cover a reference
/// position and their bases disagree, the base quality of both is set to 0 (uninformative).
fn adjust_overlapping_mates(rs: &mut [Read]) {
    let mut by_name: HashMap<Vec<u8>, Vec<usize>> = HashMap::new();
    for (i, r) in rs.iter().enumerate() {
        if r.flags & 0x900 == 0 {
            by_name.entry(r.name.clone()).or_default().push(i);
        }
    }
    let ref_to_query = |r: &Read| -> HashMap<u64, usize> {
        let mut m = HashMap::new();
        let (mut rp, mut qp) = (r.start, 0usize);
        for &(op, len) in &r.cigar {
            let len = len as usize;
            match op {
                Op::Match => {
                    for k in 0..len {
                        m.insert(rp + k as u64, qp + k);
                    }
                    rp += len as u64;
                    qp += len;
                }
                Op::Ins | Op::SoftClip => qp += len,
                Op::Del | Op::Skip => rp += len as u64,
                _ => {}
            }
        }
        m
    };
    for idx in by_name.values().filter(|v| v.len() == 2) {
        let (a, b) = (idx[0], idx[1]);
        if rs[a].tid != rs[b].tid || rs[a].start >= rs[b].end || rs[b].start >= rs[a].end {
            continue;
        }
        let (ma, mb) = (ref_to_query(&rs[a]), ref_to_query(&rs[b]));
        let mut changes: Vec<(usize, usize, usize)> = vec![]; // (which read, query pos, new qual) encoded
        for (pos, &qa) in &ma {
            if let Some(&qb) = mb.get(pos) {
                let (ba, bb) = (rs[a].seq.get(qa).copied(), rs[b].seq.get(qb).copied());
                if ba.is_none() || bb.is_none() {
                    continue;
                }
                if ba != bb {
                    changes.push((a, qa, 0));
                    changes.push((b, qb, 0));
                }
            }
        }
        for (r, q, v) in changes {
            if let Some(x) = rs[r].qual.get_mut(q) {
                *x = v as u8;
            }
        }
    }
}

/// log10 P(read | haplotype) with a pair-HMM (global in the read, local in the haplotype), GATK style
/// GATK's PCR indel error model (CONSERVATIVE, the HaplotypeCaller default): gap-open quality per
/// read position from the tandem repeat around it in the read, Q = max(10, round(41 - e^(L/3pi))),
/// L = repeat units (unit 1-8 bases) capped at 20; Q40 outside repeats, about Q33 in long ones.
fn gap_open_quals(read: &[u8]) -> Vec<f64> {
    let n = read.len();
    let up = |b: u8| b.to_ascii_uppercase();
    (0..n)
        .map(|i| {
            // repeat unit ending at i (backward) with more than one copy, the shortest first
            let mut best: Option<(usize, usize)> = None; // (unit length, copies backward)
            for u in 1..=8usize {
                if i + 1 < u {
                    break;
                }
                let unit = &read[i + 1 - u..=i];
                let mut copies = 1;
                while i + 1 >= u * (copies + 1) && read[i + 1 - u * (copies + 1)..i + 1 - u * copies].iter().zip(unit).all(|(a, b)| up(*a) == up(*b)) {
                    copies += 1;
                }
                if copies > 1 {
                    best = Some((u, copies));
                    break;
                }
            }
            let (u, back) = best.unwrap_or((1, 1));
            let unit = &read[i + 1 - u..=i];
            let mut fwd = 0;
            while i + 1 + u * (fwd + 1) <= n && read[i + 1 + u * fwd..i + 1 + u * (fwd + 1)].iter().zip(unit).all(|(a, b)| up(*a) == up(*b)) {
                fwd += 1;
            }
            let l = ((back + fwd) as f64).min(20.0);
            (41.0 - (l / (3.0 * std::f64::consts::PI)).exp()).round().max(10.0)
        })
        .collect()
}

fn pair_hmm(read: &[u8], qual: &[u8], hap: &[u8]) -> f64 {
    let (n, m) = (read.len(), hap.len());
    if n == 0 || m == 0 {
        return 0.0;
    }
    let gq = gap_open_quals(read);
    let eps = 0.1; // gap continuation Q10
    let gm = 1.0 - eps;
    let init = 2f64.powi(1020) / m as f64;
    let mut pm = vec![0f64; m + 1];
    let mut pi = vec![0f64; m + 1];
    let mut pd = vec![init; m + 1];
    pd[0] = 0.0;
    let (mut cm, mut ci, mut cd) = (vec![0f64; m + 1], vec![0f64; m + 1], vec![0f64; m + 1]);
    let mut prev_delta = 10f64.powf(-gq[0] / 10.0);
    for i in 1..=n {
        let q = qual.get(i - 1).copied().unwrap_or(30) as f64;
        let e = 10f64.powf(-q / 10.0).min(0.75); // quality 0 (mates disagree, long clips): uninformative
        let rb = read[i - 1].to_ascii_uppercase();
        let delta = 10f64.powf(-gq[i - 1] / 10.0);
        let mm = 1.0 - 2.0 * prev_delta;
        cm[0] = 0.0;
        ci[0] = 0.0;
        cd[0] = 0.0;
        for j in 1..=m {
            let hb = hap[j - 1].to_ascii_uppercase();
            let prior = if rb == hb || rb == b'N' || hb == b'N' { 1.0 - e } else { e / 3.0 };
            cm[j] = prior * (pm[j - 1] * mm + pi[j - 1] * gm + pd[j - 1] * gm);
            ci[j] = pm[j] * delta + pi[j] * eps;
            cd[j] = cm[j - 1] * delta + cd[j - 1] * eps;
        }
        prev_delta = delta;
        std::mem::swap(&mut pm, &mut cm);
        std::mem::swap(&mut pi, &mut ci);
        std::mem::swap(&mut pd, &mut cd);
    }
    let sum: f64 = (1..=m).map(|j| pm[j] + pi[j]).sum();
    sum.log10() - 1020.0 * 2f64.log10()
}

fn log10_sum2(a: f64, b: f64) -> f64 {
    let m = a.max(b);
    m + (10f64.powf(a - m) + 10f64.powf(b - m)).log10()
}

/// ln of n! for the Fisher exact test
fn ln_fact(n: u32) -> f64 {
    (1..=n).map(|k| (k as f64).ln()).sum()
}

/// two-sided Fisher exact test p-value of the 2x2 table [[a, b], [c, d]]
fn fisher(a: u32, b: u32, c: u32, d: u32) -> f64 {
    let (r1, r2, c1, n) = (a + b, c + d, a + c, a + b + c + d);
    let p_of = |x: u32| -> f64 {
        let (aa, bb, cc, dd) = (x, r1 - x, c1 - x, r2 + x - c1);
        (ln_fact(r1) + ln_fact(r2) + ln_fact(c1) + ln_fact(n - c1) - ln_fact(n) - ln_fact(aa) - ln_fact(bb) - ln_fact(cc) - ln_fact(dd)).exp()
    };
    let lo = c1.saturating_sub(r2);
    let hi = r1.min(c1);
    let p0 = p_of(a);
    let mut p = 0.0;
    for x in lo..=hi {
        let px = p_of(x);
        if px <= p0 * (1.0 + 1e-7) {
            p += px;
        }
    }
    p.min(1.0)
}

/// GATK StrandOddsRatio from (ref fwd, ref rev, alt fwd, alt rev)
fn sor(rf: u32, rr: u32, af: u32, ar: u32) -> f64 {
    let (rf, rr, af, ar) = (rf as f64 + 1.0, rr as f64 + 1.0, af as f64 + 1.0, ar as f64 + 1.0);
    let r = (rf * ar) / (rr * af);
    let sym = r + 1.0 / r;
    let ref_ratio = rf.min(rr) / rf.max(rr);
    let alt_ratio = af.min(ar) / af.max(ar);
    sym.ln() + ref_ratio.ln() - alt_ratio.ln()
}

/// window of the reference [ws, we) with a set of non-overlapping alleles applied
fn haplotype(rf: &[u8], off: u64, ws: u64, we: u64, alleles: &[&Allele]) -> Vec<u8> {
    let mut h = vec![];
    let mut pos = ws;
    let mut sorted: Vec<&&Allele> = alleles.iter().collect();
    sorted.sort_by_key(|a| a.pos);
    for a in sorted {
        if a.pos < pos {
            continue;
        }
        h.extend(&rf[(pos - off) as usize..(a.pos - off) as usize]);
        h.extend(&a.alt);
        pos = a.pos + a.r#ref.len() as u64;
    }
    if pos < we {
        h.extend(&rf[(pos - off) as usize..(we - off) as usize]);
    }
    h
}

/// Length of the longest tandem repeat tract (unit 1 to 4 bases, at least 2 copies) in the reference
/// that touches [i, i + len] (0-based offsets into `rf`).
fn tract_len(rf: &[u8], i: usize, len: usize) -> usize {
    let up = |b: u8| b.to_ascii_uppercase();
    let hi = (i + len + 1).min(rf.len());
    let mut best = 0;
    for u in 1..=4usize {
        for s in i..hi {
            if s + u > rf.len() {
                break;
            }
            // the tract with the unit rf[s..s+u], extended right and left in phase
            let at = |x: usize| up(rf[s + (x + u * 1024 - s) % u]);
            let mut b = s + u;
            while b < rf.len() && up(rf[b]) == at(b) {
                b += 1;
            }
            let mut l = s;
            while l > 0 && up(rf[l - 1]) == at(l - 1) {
                l -= 1;
            }
            if b - l >= 2 * u {
                best = best.max(b - l);
            }
        }
    }
    best
}

fn overlaps(a: &Allele, b: &Allele) -> bool {
    a.pos < b.pos + b.r#ref.len() as u64 && b.pos < a.pos + a.r#ref.len() as u64
}

/// Call variants in [start, end) (0-based) of `chrom`: candidates are grouped into clusters (gaps
/// <= 10 bp); every read is aligned (pair-HMM) against the reference and the candidate haplotypes of
/// its cluster; the best diploid pair of haplotypes gives the genotypes.
static MISMAP: std::sync::LazyLock<f64> = std::sync::LazyLock::new(|| std::env::var("BRCA12_MISMAP").ok().and_then(|v| v.parse().ok()).unwrap_or(5.0));

pub fn call_region(path: &str, reference: &str, chrom: &str, start: u64, end: u64, p: &Params) -> io::Result<Vec<Call>> {
    let pad = 1000u64;
    let off = start.saturating_sub(pad);
    let mut fa = noodles::fasta::io::indexed_reader::Builder::default().build_from_path(reference)?;
    let region: noodles::core::Region = format!("{chrom}:{}-{}", off + 1, end + pad).parse().map_err(io::Error::other)?;
    let rf: Vec<u8> = fa.query(&region)?.sequence().as_ref().iter().map(|b| b.to_ascii_uppercase()).collect();
    let mut rs: Vec<Read> = vec![];
    reads::fetch(path, Some(reference), &format!("{chrom}:{}-{}", start + 1, end), &mut |r: &Read| {
        if usable(r, p) {
            let mut r = r.clone();
            clip_read(&mut r);
            if !r.seq.is_empty() && !r.cigar.is_empty() {
                rs.push(r);
            }
        }
    })?;
    adjust_overlapping_mates(&mut rs);
    let cands = candidates(&rs, &rf, off, p);
    let mut keys: Vec<Allele> = cands.keys().filter(|a| a.pos >= start && a.pos < end).cloned().collect();
    keys.sort();
    // clusters of candidates
    let mut clusters: Vec<Vec<Allele>> = vec![];
    for a in keys {
        match clusters.last_mut() {
            Some(c) if a.pos <= c.iter().map(|x| x.pos + x.r#ref.len() as u64).max().unwrap() + 40 => c.push(a),
            _ => clusters.push(vec![a]),
        }
    }
    let mut calls = vec![];
    for cl0 in clusters {
        // alignment candidates: the best supported 5 of a cluster form all combinations
        let mut cl = cl0.clone();
        if cl.len() > 5 {
            cl.sort_by_key(|a| std::cmp::Reverse(cands[a].0));
            cl.truncate(5);
            cl.sort();
        }
        let n_sub = cl.len();
        // local assembly over the cluster window: complex alleles the aligner placed differently
        let mut asm_haps: Vec<Vec<Allele>> = vec![];
        if p.assemble {
            let cs0 = cl0.iter().map(|a| a.pos).min().unwrap();
            let ce0 = cl0.iter().map(|a| a.pos + a.r#ref.len() as u64).max().unwrap();
            let ws0 = cs0.saturating_sub(p.flank as u64).max(off);
            let we0 = (ce0 + p.flank as u64).min(off + rf.len() as u64);
            let refw = &rf[(ws0 - off) as usize..(we0 - off) as usize];
            let reads: Vec<(&[u8], &[u8])> = rs
                .iter()
                .filter(|r| r.start < we0 && r.end > ws0)
                .map(|r| {
                    // long soft clips are not used
                    let a = match r.cigar.first() { Some(&(Op::SoftClip, k)) if k >= p.long_clip => k as usize, _ => 0 };
                    let b = match r.cigar.last() { Some(&(Op::SoftClip, k)) if k >= p.long_clip => r.seq.len() - k as usize, _ => r.seq.len() };
                    let b = b.max(a);
                    (&r.seq[a..b], &r.qual[a..b])
                })
                .collect();
            let dbg = std::env::var("BRCA12_DEBUG_POS").ok().and_then(|v| v.parse::<u64>().ok()).is_some_and(|d| d >= ws0 && d < we0);
            let assembled = assemble::assemble(refw, &reads, &assemble::DEFAULT);
            if dbg {
                eprintln!("assembly window {}-{}: {} haplotypes", ws0 + 1, we0, assembled.len());
                for h in &assembled {
                    eprintln!("  {}", String::from_utf8_lossy(h));
                }
            }
            for h in assembled.iter().skip(1) {
                let mut ev: Vec<Allele> = assemble::events(refw, h).into_iter().map(|(o, r, a)| Allele { pos: ws0 + o as u64, r#ref: r, alt: a }).collect();
                // left-align each event only where the haplotype stays the same (an indel must not
                // slide over a neighbouring event of the same haplotype)
                for i in 0..ev.len() {
                    let n = normalise(ev[i].clone(), &rf, off);
                    if n != ev[i] {
                        let mut trial = ev.clone();
                        trial[i] = n;
                        let ok = trial.iter().enumerate().all(|(x, a)| trial[x + 1..].iter().all(|b| !overlaps(a, b)))
                            && haplotype(&rf, off, ws0, we0, &trial.iter().collect::<Vec<_>>()) == *h;
                        if ok {
                            ev = trial;
                        }
                    }
                }
                if dbg {
                    eprintln!("events {:?}", ev.iter().map(|a| (a.pos + 1, String::from_utf8_lossy(&a.r#ref).into_owned(), String::from_utf8_lossy(&a.alt).into_owned())).collect::<Vec<_>>());
                }
                let ev: Vec<Allele> = ev.into_iter().filter(|a| a.pos + 20 >= cs0 && a.pos < ce0 + 20 && a.pos >= start && a.pos < end).collect();
                if !ev.is_empty() && !asm_haps.contains(&ev) {
                    asm_haps.push(ev);
                }
            }
        }
        for h in &asm_haps {
            for a in h {
                if !cl.contains(a) {
                    cl.push(a.clone());
                }
            }
        }
        let cs = cl.iter().map(|a| a.pos).min().unwrap();
        let ce = cl.iter().map(|a| a.pos + a.r#ref.len() as u64).max().unwrap();
        let ws = cs.saturating_sub(p.flank as u64).max(off);
        let we = (ce + p.flank as u64).min(off + rf.len() as u64);
        // haplotypes: every combination of the alignment candidates without overlapping alleles,
        // every assembled haplotype, and each assembled haplotype with one more alignment candidate.
        // An alignment candidate that overlaps an assembled event without being it is another
        // representation of that allele and is left out (representations must not be mixed).
        let no_overlap = |set: &[usize]| set.iter().enumerate().all(|(x, &i)| set[x + 1..].iter().all(|&j| !overlaps(&cl[i], &cl[j])));
        let conflict: Vec<bool> = (0..n_sub)
            .map(|i| asm_haps.iter().flatten().any(|e| e != &cl[i] && overlaps(e, &cl[i])))
            .collect();
        let free: Vec<usize> = (0..n_sub).filter(|&i| !conflict[i]).collect();
        let mut haps: Vec<Vec<usize>> = vec![];
        for mask in 0u32..(1u32 << free.len()) {
            let set: Vec<usize> = (0..free.len()).filter(|i| mask & (1 << i) != 0).map(|i| free[i]).collect();
            if no_overlap(&set) {
                haps.push(set);
            }
        }
        for h in &asm_haps {
            let mut set: Vec<usize> = h.iter().map(|a| cl.iter().position(|x| x == a).unwrap()).collect();
            set.sort();
            set.dedup();
            if !no_overlap(&set) {
                continue;
            }
            let mut add = vec![set.clone()];
            for &i in &free {
                if !set.contains(&i) {
                    let mut t = set.clone();
                    t.push(i);
                    t.sort();
                    if no_overlap(&t) {
                        add.push(t);
                    }
                }
            }
            for t in add {
                if haps.len() < p.max_haplotypes && !haps.contains(&t) {
                    haps.push(t);
                }
            }
        }
        let seqs: Vec<Vec<u8>> = haps.iter().map(|h| haplotype(&rf, off, ws, we, &h.iter().map(|&i| &cl[i]).collect::<Vec<_>>())).collect();
        let rds: Vec<&Read> = rs.iter().filter(|r| r.start < ce && r.end > cs).collect();
        if rds.is_empty() {
            continue;
        }
        // base qualities as GATK HaplotypeCaller: below Q18 lowered to Q6, capped at the read's MAPQ
        // (quality 0 of disagreeing mates stays uninformative)
        let quals: Vec<Vec<u8>> = rds
            .iter()
            .map(|r| {
                let mut q: Vec<u8> = r.qual.iter().map(|&q| if q == 0 { 0 } else if q < p.low_base_q { 6 } else { q.min(r.mapq) }).collect();
                // long soft clips (chimeric or junk read ends) are no evidence: uninformative; short
                // ones stay (a read ending inside a homopolymer is often clipped by a few bases)
                if let Some(&(Op::SoftClip, k)) = r.cigar.first().filter(|c| c.1 >= p.long_clip) {
                    q.iter_mut().take(k as usize).for_each(|x| *x = 0);
                }
                if let Some(&(Op::SoftClip, k)) = r.cigar.last().filter(|c| c.1 >= p.long_clip) {
                    let n = q.len();
                    q.iter_mut().skip(n.saturating_sub(k as usize)).for_each(|x| *x = 0);
                }
                q
            })
            .collect();
        let lik: Vec<Vec<f64>> = rds.iter().zip(&quals).map(|(r, q)| seqs.iter().map(|h| pair_hmm(&r.seq, q, h)).collect()).collect();
        // one-strand alleles: an allele that only reads of one strand support while the other strand
        // covers the site with reads against it is a sequencing artifact (e.g. phasing errors after a
        // long homopolymer); haplotypes containing it are not considered
        let mut artifact = vec![false; cl.len()];
        for i in 0..cl.len() {
            let (mut sf, mut sr, mut nf, mut nr) = (0u32, 0u32, 0u32, 0u32);
            for (k, l) in lik.iter().enumerate() {
                let with = haps.iter().enumerate().filter(|(_, h)| h.contains(&i)).map(|(x, _)| l[x]).fold(f64::NEG_INFINITY, f64::max);
                let without = haps.iter().enumerate().filter(|(_, h)| !h.contains(&i)).map(|(x, _)| l[x]).fold(f64::NEG_INFINITY, f64::max);
                let rev = rds[k].flags & 0x10 != 0;
                if with - without > 0.2 {
                    if rev { sr += 1 } else { sf += 1 }
                } else if without - with > 0.2 {
                    if rev { nr += 1 } else { nf += 1 }
                }
            }
            let one_sided = (sf == 0 && nf >= p.strand_other_min && sr >= p.strand_support_min)
                || (sr == 0 && nr >= p.strand_other_min && sf >= p.strand_support_min);
            if one_sided && -10.0 * fisher(nf, nr, sf, sr).max(1e-300).log10() > p.strand_fs {
                artifact[i] = true;
            }
        }
        let usable_h: Vec<usize> = (0..haps.len()).filter(|&x| !haps[x].iter().any(|&i| artifact[i])).collect();
        // sites: groups of overlapping alleles, genotyped one by one with the read likelihoods
        // marginalised over the haplotypes (as GATK): an allele draws its evidence only from the
        // reads that need it, not from other alleles on the same haplotype
        let alive: Vec<usize> = (0..cl.len()).filter(|&i| !artifact[i]).collect();
        let mut sites: Vec<Vec<usize>> = vec![];
        for &i in &alive {
            let hit: Vec<usize> = (0..sites.len()).filter(|&k| sites[k].iter().any(|&j| overlaps(&cl[i], &cl[j]))).collect();
            let mut merged = vec![i];
            for &k in hit.iter().rev() {
                merged.extend(sites.remove(k));
            }
            merged.sort();
            sites.push(merged);
        }
        sites.sort();
        for site in sites {
            // the site's alleles are the distinct combinations of its events on the haplotypes (an
            // assembled complex allele can hold several events of one site)
            let sig_of = |x: usize| -> Vec<usize> { site.iter().cloned().filter(|a| haps[x].contains(a)).collect() };
            let mut sigs: Vec<Vec<usize>> = vec![];
            for &x in &usable_h {
                let sg = sig_of(x);
                if !sg.is_empty() && !sigs.contains(&sg) {
                    sigs.push(sg);
                }
            }
            sigs.sort();
            let k = sigs.len();
            if k == 0 {
                continue;
            }
            // per read: log10 likelihood of each allele of the site (0 = none of them)
            let hap_idx: Vec<(usize, usize)> = usable_h
                .iter()
                .map(|&x| {
                    let sg = sig_of(x);
                    (x, if sg.is_empty() { 0 } else { sigs.iter().position(|y| y == &sg).unwrap() + 1 })
                })
                .collect();
            let la: Vec<Vec<f64>> = lik
                .iter()
                .map(|l| {
                    let mut v = vec![-1e9f64; k + 1];
                    for &(x, idx) in &hap_idx {
                        v[idx] = v[idx].max(l[x]);
                    }
                    // a read is at most this much (log10) more likely on one allele than on another
                    // (as GATK HaplotypeCaller's global read mismapping rate, there Q45): one misplaced
                    // read whose mismatches fit a large indel perfectly must not outweigh dozens of
                    // reference reads. 5.0 instead of 4.5: at 4.5 HG007 NovaSeq chr13:32347011
                    // T>TTTTATAAAA (2 of 7 reads, true 0/1) fell to QUAL 29.8, just below the call limit
                    if *MISMAP > 0.0 {
                        let mx = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                        for x in v.iter_mut() {
                            if *x > -1e8 {
                                *x = x.max(mx - *MISMAP);
                            }
                        }
                    }
                    v
                })
                .collect();
            if std::env::var("BRCA12_DEBUG_POS").ok().and_then(|v| v.parse::<u64>().ok()).is_some_and(|d| site.iter().any(|&a| cl[a].pos + 1 == d)) {
                eprintln!("site {:?}", site.iter().map(|&a| (cl[a].pos + 1, String::from_utf8_lossy(&cl[a].r#ref).into_owned(), String::from_utf8_lossy(&cl[a].alt).into_owned())).collect::<Vec<_>>());
                eprintln!("cluster {:?} haps {:?}", cl.iter().map(|a| a.pos + 1).collect::<Vec<_>>(), haps);
                for (r, v) in rds.iter().zip(&la) {
                    if std::env::var("BRCA12_DEBUG_ALL").is_ok() || v[1..].iter().any(|&x| x > v[0] + 0.2) {
                        eprintln!("{} flags={} start={} {:?}", String::from_utf8_lossy(&r.name), r.flags, r.start + 1, v.iter().map(|x| (x * 100.0).round() / 100.0).collect::<Vec<_>>());
                    }
                }
            }
            // stutter: in a long repeat tract a read can show a neighbouring tract length by
            // polymerase slippage, so one read argues at most `stutter_cap` (log10) against an
            // indel allele there; otherwise a few slipped reads turn a homozygous indel into 0/1
            let all_indel = sigs.iter().all(|sg| sg.iter().all(|&e| cl[e].r#ref.len() != cl[e].alt.len()));
            let tract = site.iter().map(|&a| tract_len(&rf, (cl[a].pos - off) as usize, cl[a].r#ref.len())).max().unwrap_or(0);
            let cap = std::env::var("BRCA12_STUTTER_CAP").ok().and_then(|v| v.parse::<f64>().ok()).unwrap_or(p.stutter_cap);
            // used only to choose among the non-reference genotypes (0/1 or 1/1), never for calling
            let la_c: Option<Vec<Vec<f64>>> = (all_indel && tract >= p.stutter_tract && cap > 0.0).then(|| {
                la.iter()
                    .map(|v| {
                        let mx = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                        v.iter().map(|&x| if x > -1e8 { x.max(mx - cap) } else { x }).collect()
                    })
                    .collect()
            });
            let het = |a: usize| if sigs[a - 1].iter().any(|&e| cl[e].r#ref.len() != cl[e].alt.len()) { p.het_indel } else { p.het_snv };
            let mut gts: Vec<(usize, usize, f64)> = vec![];
            for x in 0..=k {
                for y in x..=k {
                    let mut ll = 0.0;
                    for v in &la {
                        ll += log10_sum2(v[x], v[y]) - 2f64.log10();
                    }
                    let prior = match (x, y) {
                        (0, 0) => 0.0,
                        (0, b) => het(b).log10(),
                        (a, b) if a == b => (het(a) / 2.0).log10(),
                        (a, b) => het(a).log10() + het(b).log10(),
                    };
                    gts.push((x, y, ll + prior));
                }
            }
            // two different alleles: one that only reads of one strand prefer over the other allele and
            // REF, while reads of the other strand cover the site against it, is a sequencing artifact
            // (homopolymer phasing); the genotypes with it are dropped
            let first = gts.iter().cloned().max_by(|a, b| a.2.partial_cmp(&b.2).unwrap()).unwrap();
            if first.0 != 0 && first.0 != first.1 {
                for (a, b) in [(first.0, first.1), (first.1, first.0)] {
                    let (mut sf, mut sr, mut nf, mut nr) = (0u32, 0u32, 0u32, 0u32);
                    for (r, v) in rds.iter().zip(&la) {
                        let rev = r.flags & 0x10 != 0;
                        let other = v[b].max(v[0]);
                        if v[a] - other > 0.2 {
                            if rev { sr += 1 } else { sf += 1 }
                        } else if other - v[a] > 0.2 {
                            if rev { nr += 1 } else { nf += 1 }
                        }
                    }
                    let one_sided = (sf == 0 && nf >= p.strand_other_min && sr >= p.strand_support_min)
                        || (sr == 0 && nr >= p.strand_other_min && sf >= p.strand_support_min);
                    if std::env::var("BRCA12_DEBUG_POS").ok().and_then(|v| v.parse::<u64>().ok()).is_some_and(|d| site.iter().any(|&x| cl[x].pos + 1 == d)) {
                        eprintln!("strand allele {a}: support f{sf} r{sr} against f{nf} r{nr}");
                    }
                    if one_sided && -10.0 * fisher(nf, nr, sf, sr).max(1e-300).log10() > p.strand_fs_two_alleles {
                        gts.retain(|g| g.0 != a && g.1 != a);
                        break;
                    }
                }
            }
            let mx = gts.iter().map(|g| g.2).fold(f64::NEG_INFINITY, f64::max);
            let tot = mx + gts.iter().map(|g| 10f64.powf(g.2 - mx)).sum::<f64>().log10();
            let best = gts.iter().cloned().max_by(|a, b| a.2.partial_cmp(&b.2).unwrap()).unwrap();
            let (mut bx, mut by) = (best.0, best.1);
            if let (Some(lc), true) = (&la_c, by != 0) {
                let mut top = f64::NEG_INFINITY;
                for &(x, y, _) in gts.iter().filter(|g| g.1 != 0) {
                    let prior = match (x, y) {
                        (0, b) => het(b).log10(),
                        (a, b) if a == b => (het(a) / 2.0).log10(),
                        (a, b) => het(a).log10() + het(b).log10(),
                    };
                    let ll: f64 = lc.iter().map(|v| log10_sum2(v[x], v[y]) - 2f64.log10()).sum::<f64>() + prior;
                    if ll > top {
                        top = ll;
                        (bx, by) = (x, y);
                    }
                }
            }
            if std::env::var("BRCA12_DEBUG_POS").ok().and_then(|v| v.parse::<u64>().ok()).is_some_and(|d| site.iter().any(|&a| cl[a].pos + 1 == d)) {
                eprintln!("sigs {sigs:?} gts {gts:?} artifact {artifact:?}");
            }
            if by == 0 {
                continue;
            }
            let qual = (-10.0 * (gts[0].2 - tot)).min(10000.0);
            if qual < p.call_qual {
                continue;
            }
            let mut sorted: Vec<f64> = gts.iter().map(|g| g.2).collect();
            sorted.sort_by(|a, b| b.partial_cmp(a).unwrap());
            let gq = ((sorted[0] - sorted.get(1).copied().unwrap_or(f64::NEG_INFINITY)) * 10.0).round().clamp(0.0, 99.0) as u32;
            // records: one per event of the called allele; two different alleles at the site give one
            // record over their span with a common REF and two ALTs
            let span = |sg: &Vec<usize>| {
                (sg.iter().map(|&e| cl[e].pos).min().unwrap(), sg.iter().map(|&e| cl[e].pos + cl[e].r#ref.len() as u64).max().unwrap())
            };
            let records: Vec<(Allele, Option<Vec<u8>>, (u8, u8))> = if bx != 0 && bx != by {
                let (s1, s2) = (&sigs[bx - 1], &sigs[by - 1]);
                let ((p1, e1), (p2, e2)) = (span(s1), span(s2));
                let (ps, pe) = (p1.min(p2), e1.max(e2));
                let alt_of = |sg: &Vec<usize>| {
                    let mut evs: Vec<&Allele> = sg.iter().map(|&e| &cl[e]).collect();
                    evs.sort_by_key(|a| a.pos);
                    let mut v = vec![];
                    let mut pos = ps;
                    for a in evs {
                        v.extend(&rf[(pos - off) as usize..(a.pos - off) as usize]);
                        v.extend(&a.alt);
                        pos = a.pos + a.r#ref.len() as u64;
                    }
                    v.extend(&rf[(pos - off) as usize..(pe - off) as usize]);
                    v
                };
                let refs = rf[(ps - off) as usize..(pe - off) as usize].to_vec();
                vec![(Allele { pos: ps, r#ref: refs, alt: alt_of(s1) }, Some(alt_of(s2)), (1u8, 2u8))]
            } else {
                sigs[by - 1].iter().map(|&e| (cl[e].clone(), None, if bx == 0 { (0, 1) } else { (1, 1) })).collect()
            };
            // AD and strand: each read goes to the best of REF and the called alleles (margin 0.2 log10)
            let mut opts = vec![0usize, bx, by];
            opts.sort();
            opts.dedup();
            let (mut adr, mut ada, mut rf_, mut rr_, mut af_, mut ar_) = (0u32, 0u32, 0u32, 0u32, 0u32, 0u32);
            for (r, v) in rds.iter().zip(&la) {
                let mut o: Vec<(f64, usize)> = opts.iter().map(|&x| (v[x], x)).collect();
                o.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
                if o[0].0 - o[1].0 <= 0.2 {
                    continue;
                }
                let rev = r.flags & 0x10 != 0;
                if o[0].1 == 0 {
                    adr += 1;
                    if rev { rr_ += 1 } else { rf_ += 1 }
                } else {
                    ada += 1;
                    if rev { ar_ += 1 } else { af_ += 1 }
                }
            }
            let indel = records.iter().any(|(a, a2, _)| a.r#ref.len() != a.alt.len() || a2.is_some());
            let fs = -10.0 * fisher(rf_, rr_, af_, ar_).max(1e-300).log10();
            let so = sor(rf_, rr_, af_, ar_);
            let qd = if ada > 0 { qual / ada as f64 } else { 0.0 };
            let mut filters = vec![];
            if qd < 2.0 {
                filters.push("QD");
            }
            if (!indel && fs > 60.0) || (indel && fs > 200.0) {
                filters.push("FS");
            }
            // SOR is reported but not used as a filter on its own: it removed real SNVs in HG001
            for (allele, alts, gt) in records {
                calls.push(Call {
                    chrom: chrom.to_string(),
                    allele,
                    alt2: alts,
                    qual,
                    gt,
                    gq,
                    ad: (adr, ada),
                    dp: rds.len() as u32,
                    pl: [0, 0, 0],
                    filter: if filters.is_empty() { "PASS".into() } else { filters.join(";") },
                    fs,
                    sor: so,
                });
            }
        }
    }
    calls.sort_by(|a, b| a.allele.pos.cmp(&b.allele.pos));
    Ok(calls)
}

pub fn write_vcf(path: &str, sample: &str, contigs: &[(String, u64)], calls: &[Call]) -> io::Result<()> {
    use std::io::Write;
    let mut f = std::fs::File::create(path)?;
    writeln!(f, "##fileformat=VCFv4.2")?;
    writeln!(f, "##source=brca12-wgs {}", env!("CARGO_PKG_VERSION"))?;
    for (c, l) in contigs.iter().filter(|(c, _)| c == "chr13" || c == "chr17" || c == "13" || c == "17") {
        writeln!(f, "##contig=<ID={c},length={l}>")?;
    }
    writeln!(f, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")?;
    writeln!(f, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Reads supporting REF and ALT (pair-HMM likelihood ratio > 0.2 log10)\">")?;
    writeln!(f, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Reads overlapping the variant\">")?;
    writeln!(f, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">")?;
    writeln!(f, "##INFO=<ID=FS,Number=1,Type=Float,Description=\"Phred-scaled Fisher exact test for strand bias\">")?;
    writeln!(f, "##INFO=<ID=SOR,Number=1,Type=Float,Description=\"Strand odds ratio (GATK)\">")?;
    writeln!(f, "##FILTER=<ID=QD,Description=\"QUAL / alt reads < 2\">")?;
    writeln!(f, "##FILTER=<ID=FS,Description=\"FS > 60 (SNV) or > 200 (indel)\">")?;
    writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}")?;
    for c in calls {
        let alt = match &c.alt2 {
            Some(a2) => format!("{},{}", String::from_utf8_lossy(&c.allele.alt), String::from_utf8_lossy(a2)),
            None => String::from_utf8_lossy(&c.allele.alt).into_owned(),
        };
        writeln!(
            f,
            "{}\t{}\t.\t{}\t{}\t{:.2}\t{}\tFS={:.3};SOR={:.3}\tGT:AD:DP:GQ\t{}/{}:{},{}:{}:{}",
            c.chrom, c.allele.pos + 1, String::from_utf8_lossy(&c.allele.r#ref), alt, c.qual, c.filter, c.fs, c.sor,
            c.gt.0, c.gt.1, c.ad.0, c.ad.1, c.dp, c.gq
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod stutter_tests {
    use super::tract_len;

    #[test]
    fn tract_length_of_homopolymers_and_short_repeats() {
        // 0-based offsets: "GC" + 19 T + "GAG"
        let rf = b"GCTTTTTTTTTTTTTTTTTTTGAG";
        assert_eq!(tract_len(rf, 1, 3), 19);
        // dinucleotide repeat CA x 6
        let rf = b"GTCACACACACACAGT";
        assert_eq!(tract_len(rf, 1, 3), 12);
        // no repeat: the longest stretch is a 2-base run
        assert_eq!(tract_len(b"ACGTACGGTA", 5, 1), 2);
    }
}
