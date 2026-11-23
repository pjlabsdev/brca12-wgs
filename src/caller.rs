//! Small-variant caller for the BRCA1/BRCA2 regions (SNVs and indels), built like the core of GATK
//! HaplotypeCaller without de-novo assembly. Research use only.
//!
//! 1. Reads: MAPQ >= 20; secondary, duplicate and QC-fail records skipped.
//! 2. Candidates from the alignments: mismatches with base quality >= 10, insertions and deletions
//!    (left-normalised), kept with >= 2 supporting fragments.
//! 3. Each candidate: every read overlapping it is aligned with a pair-HMM (gap open from GATK's PCR indel model, gap
//!    continuation Q10, base qualities below Q18 as Q6 and capped at MAPQ as GATK, base qualities as error rates; global in the read, local in the haplotype)
//!    against the reference window and the window with the candidate applied.
//! 4. Genotype likelihoods for 0/0, 0/1, 1/1 from the read likelihoods; priors as GATK
//!    (heterozygosity 0.001 for SNVs, 1.25e-4 for indels); QUAL = phred P(0/0 | reads); called when
//!    QUAL >= 30 and the best genotype is not 0/0.
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
    pub long_clip: u32,          // soft clips of at least this length are uninformative          // base qualities below this count as Q6 in the pair-HMM (GATK: 18)
    pub strand_support_min: u32, // one-strand allele: at least this many supporting reads,
    pub strand_other_min: u32,   // at least this many reads of the other strand against it,
    pub strand_fs: f64,          // and Fisher strand test above this (phred)
}

pub const DEFAULT: Params =
    Params { min_mapq: 20, min_base_q: 10, min_support: 2, call_qual: 30.0, het_snv: 1e-3, het_indel: 1.25e-4, flank: 150,
             low_base_q: 18, long_clip: 20, strand_support_min: 4, strand_other_min: 8, strand_fs: 20.0 };

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

fn overlaps(a: &Allele, b: &Allele) -> bool {
    a.pos < b.pos + b.r#ref.len() as u64 && b.pos < a.pos + a.r#ref.len() as u64
}

/// Call variants in [start, end) (0-based) of `chrom`: candidates are grouped into clusters (gaps
/// <= 10 bp); every read is aligned (pair-HMM) against the reference and the candidate haplotypes of
/// its cluster; the best diploid pair of haplotypes gives the genotypes.
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
    for mut cl in clusters {
        // keep the best supported candidates of large clusters
        if cl.len() > 5 {
            cl.sort_by_key(|a| std::cmp::Reverse(cands[a].0));
            cl.truncate(5);
            cl.sort();
        }
        let cs = cl.iter().map(|a| a.pos).min().unwrap();
        let ce = cl.iter().map(|a| a.pos + a.r#ref.len() as u64).max().unwrap();
        let ws = cs.saturating_sub(p.flank as u64).max(off);
        let we = (ce + p.flank as u64).min(off + rf.len() as u64);
        // haplotypes: every subset of the cluster's candidates without overlapping alleles
        let mut haps: Vec<Vec<usize>> = vec![];
        for mask in 0u32..(1u32 << cl.len()) {
            let set: Vec<usize> = (0..cl.len()).filter(|i| mask & (1 << i) != 0).collect();
            if set.iter().enumerate().all(|(x, &i)| set[x + 1..].iter().all(|&j| !overlaps(&cl[i], &cl[j]))) {
                haps.push(set);
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
            let k = site.len();
            // per read: log10 likelihood of each allele of the site (0 = none of them)
            let la: Vec<Vec<f64>> = lik
                .iter()
                .map(|l| {
                    let mut v = vec![-1e9f64; k + 1];
                    for &x in &usable_h {
                        let idx = site.iter().position(|a| haps[x].contains(a)).map(|q| q + 1).unwrap_or(0);
                        v[idx] = v[idx].max(l[x]);
                    }
                    v
                })
                .collect();
            if std::env::var("BRCA12_DEBUG_POS").ok().and_then(|v| v.parse::<u64>().ok()).is_some_and(|d| site.iter().any(|&a| cl[a].pos + 1 == d)) {
                eprintln!("site {:?}", site.iter().map(|&a| (cl[a].pos + 1, String::from_utf8_lossy(&cl[a].r#ref).into_owned(), String::from_utf8_lossy(&cl[a].alt).into_owned())).collect::<Vec<_>>());
                eprintln!("cluster {:?} haps {:?}", cl.iter().map(|a| a.pos + 1).collect::<Vec<_>>(), haps);
                for (r, v) in rds.iter().zip(&la) {
                    if v[1..].iter().any(|&x| x > v[0] + 0.2) {
                        eprintln!("{} flags={} start={} {:?}", String::from_utf8_lossy(&r.name), r.flags, r.start + 1, v.iter().map(|x| (x * 100.0).round() / 100.0).collect::<Vec<_>>());
                    }
                }
            }
            let het = |a: usize| if cl[site[a - 1]].r#ref.len() != cl[site[a - 1]].alt.len() { p.het_indel } else { p.het_snv };
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
            let mx = gts.iter().map(|g| g.2).fold(f64::NEG_INFINITY, f64::max);
            let tot = mx + gts.iter().map(|g| 10f64.powf(g.2 - mx)).sum::<f64>().log10();
            let best = gts.iter().cloned().max_by(|a, b| a.2.partial_cmp(&b.2).unwrap()).unwrap();
            let (bx, by) = (best.0, best.1);
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
            let (allele, alts, gt) = if bx != 0 && bx != by {
                // two different alleles at one site: one record with a common REF and two ALTs
                let (a1, a2) = (&cl[site[bx - 1]], &cl[site[by - 1]]);
                let ps = a1.pos.min(a2.pos);
                let pe = (a1.pos + a1.r#ref.len() as u64).max(a2.pos + a2.r#ref.len() as u64);
                let refs = rf[(ps - off) as usize..(pe - off) as usize].to_vec();
                let alt_of = |a: &Allele| {
                    let mut v = rf[(ps - off) as usize..(a.pos - off) as usize].to_vec();
                    v.extend(&a.alt);
                    v.extend(&rf[(a.pos + a.r#ref.len() as u64 - off) as usize..(pe - off) as usize]);
                    v
                };
                (Allele { pos: ps, r#ref: refs, alt: alt_of(a1) }, Some(alt_of(a2)), (1u8, 2u8))
            } else {
                (cl[site[by - 1]].clone(), None, if bx == 0 { (0, 1) } else { (1, 1) })
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
            let indel = allele.r#ref.len() != allele.alt.len() || alts.is_some();
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
