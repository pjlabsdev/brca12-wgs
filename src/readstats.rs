//! Read statistics over the calling regions, computed like `samtools stats REF BAM regions...`
//! (samtools 1.23.1, stats.c) for the fields brca12-wgs 1.1 uses in its read QC. Research use only.
//!
//! Rules taken from stats.c: a record counts when it overlaps a region (pos < region end and alignment
//! end >= region begin); secondary records are skipped; per-read counts (total, mapped, properly
//! paired, MAPQ 0, length, quality, insert size) only for original records (not secondary, not
//! supplementary); duplicates and mismatches (NM) also for supplementary records; mapped bases from
//! the CIGAR only inside the region; each pair's insert size is counted for both mates and halved
//! (integer) per size; the insert-size mean uses the main 99 % of pairs.
use crate::reads::{self, Op, Read};
use serde_json::{Value, json};
use std::io;

const PAIRED: u16 = 0x1;
const PROPER: u16 = 0x2;
const UNMAP: u16 = 0x4;
const MUNMAP: u16 = 0x8;
const REVERSE: u16 = 0x10;
const MREVERSE: u16 = 0x20;
const READ1: u16 = 0x40;
const READ2: u16 = 0x80;
const SECONDARY: u16 = 0x100;
const DUP: u16 = 0x400;
const SUPPLEMENTARY: u16 = 0x800;
const NISIZE: usize = 8000;
const MAIN_BULK: f64 = 0.99;

#[derive(Default)]
struct Acc {
    reads_1st: u64,
    reads_2nd: u64,
    reads_other: u64,
    mapped_paired: u64,
    mapped_single: u64,
    properly_paired: u64,
    duplicated: u64,
    mq0: u64,
    total_len: u64,
    sum_qual: f64,
    mismatches: u64,
    bases_mapped_cigar: u64,
    inward: Vec<u64>,
    outward: Vec<u64>,
    other: Vec<u64>,
}

impl Acc {
    fn new() -> Self {
        Acc { inward: vec![0; NISIZE + 1], outward: vec![0; NISIZE + 1], other: vec![0; NISIZE + 1], ..Default::default() }
    }

    /// one record that overlaps the region [beg, end] (1-based, inclusive)
    fn add(&mut self, r: &Read, beg: u64, end: u64) {
        let f = r.flags;
        if f & SECONDARY != 0 {
            return;
        }
        let seq_len = r.seq.len() as u64;
        if seq_len == 0 {
            return;
        }
        if f & DUP != 0 {
            self.duplicated += 1;
        }
        let original = f & SUPPLEMENTARY == 0;
        let paired_and_mapped = f & PAIRED != 0 && f & UNMAP == 0 && f & MUNMAP == 0;
        if original {
            self.total_len += seq_len;
            let first = f & PAIRED == 0 || f & READ1 != 0;
            let last = f & PAIRED != 0 && f & READ2 != 0;
            // order: unpaired counts as first; READ1|READ2 together is neither (other)
            let order = if f & PAIRED == 0 { 1 } else { (if f & READ1 != 0 { 1 } else { 0 }) + (if f & READ2 != 0 { 2 } else { 0 }) };
            match order {
                1 => self.reads_1st += 1,
                2 => self.reads_2nd += 1,
                _ => self.reads_other += 1,
            }
            if (order == 1 && first) || (order == 2 && last) {
                self.sum_qual += r.qual.iter().map(|&q| q as f64).sum::<f64>();
            }
            if f & UNMAP == 0 {
                if r.mapq == 0 {
                    self.mq0 += 1;
                }
                if paired_and_mapped {
                    self.mapped_paired += 1;
                    if f & PROPER != 0 {
                        self.properly_paired += 1;
                    }
                } else {
                    self.mapped_single += 1;
                }
            }
        }
        if f & UNMAP != 0 {
            return;
        }
        if paired_and_mapped && original {
            let mut isize = r.tlen.unsigned_abs() as usize;
            if isize > NISIZE {
                isize = NISIZE;
            }
            if isize > 0 || r.tid == r.mate_tid {
                let pos_fst = r.mate_start - r.start as i64;
                let is_fst: i64 = if f & READ1 != 0 { 1 } else { -1 };
                let is_fwd: i64 = if f & REVERSE != 0 { -1 } else { 1 };
                let is_mfwd: i64 = if f & MREVERSE != 0 { -1 } else { 1 };
                if is_fwd * is_mfwd > 0 {
                    self.other[isize] += 1;
                } else if is_fst * pos_fst > 0 {
                    if is_fst * is_fwd > 0 { self.inward[isize] += 1 } else { self.outward[isize] += 1 }
                } else if is_fst * pos_fst < 0 {
                    if is_fst * is_fwd > 0 { self.outward[isize] += 1 } else { self.inward[isize] += 1 }
                } else {
                    self.inward[isize] += 1;
                }
            }
        }
        if let Some(nm) = r.nm {
            self.mismatches += nm.max(0) as u64;
        }
        // mapped bases from the CIGAR, only inside the region (stats.c with regions)
        let mut iref = r.start + 1; // 1-based
        for &(op, len) in &r.cigar {
            let len = len as i64;
            if len == 0 {
                continue;
            }
            match op {
                Op::Match => {
                    let mut n = len;
                    if (iref as i64) < beg as i64 {
                        n -= beg as i64 - iref as i64;
                    } else if iref as i64 + len - 1 > end as i64 {
                        n -= iref as i64 + len - 1 - end as i64;
                    }
                    if n > 0 {
                        self.bases_mapped_cigar += n as u64;
                    }
                    iref += len as u64;
                }
                Op::Ins => {
                    iref += len as u64;
                    if iref >= beg && iref <= end {
                        self.bases_mapped_cigar += len as u64;
                    }
                }
                _ => {}
            }
        }
    }
}

/// samtools-stats-like summary over regions [(chrom, beg0, end0)] (0-based half-open BED intervals)
pub fn region_stats(path: &str, reference: Option<&str>, regions: &[(String, u64, u64)]) -> io::Result<Value> {
    let mut acc = Acc::new();
    for (chrom, s, e) in regions {
        let (beg, end) = (s + 1, *e); // 1-based inclusive, as samtools keeps them
        // reference around the region, for NM when the records carry none (CRAM: htslib computes it)
        let pad = 5000u64;
        let ref_start = beg.saturating_sub(pad).max(1);
        let refseq: Option<Vec<u8>> = match reference {
            Some(fa) => {
                let mut rd = noodles::fasta::io::indexed_reader::Builder::default().build_from_path(fa)?;
                let region: noodles::core::Region = format!("{chrom}:{ref_start}-{}", end + pad).parse().map_err(io::Error::other)?;
                Some(rd.query(&region)?.sequence().as_ref().to_vec())
            }
            None => None,
        };
        reads::fetch(path, reference, &format!("{chrom}:{beg}-{end}"), &mut |r: &Read| {
            // is_in_regions: pos < end and alignment end >= begin; unmapped records span one base
            let endpos = if r.end > r.start { r.end } else { r.start + 1 };
            if r.start < end && endpos >= beg {
                if r.nm.is_none() && r.flags & UNMAP == 0 {
                    if let Some(rs) = &refseq {
                        let mut rr = r.clone();
                        rr.nm = edit_distance(r, rs, ref_start - 1);
                        acc.add(&rr, beg, end);
                        return;
                    }
                }
                acc.add(r, beg, end);
            }
        })?;
    }
    Ok(summary(&acc))
}

/// NM as htslib computes it (bam_fillmd): mismatches in M/=/X (N counts as a mismatch) plus the
/// inserted and deleted bases. `ref0` is the 0-based position of refseq[0].
pub fn edit_distance(r: &Read, refseq: &[u8], ref0: u64) -> Option<i64> {
    let code = |b: u8| match b.to_ascii_uppercase() { b'A' => 1, b'C' => 2, b'G' => 4, b'T' => 8, _ => 15 };
    let (mut rp, mut qp, mut nm) = (r.start, 0usize, 0i64);
    for &(op, len) in &r.cigar {
        let len = len as u64;
        match op {
            Op::Match => {
                for k in 0..len {
                    let ri = (rp + k).checked_sub(ref0)? as usize;
                    let (c1, c2) = (code(*refseq.get(ri)?), code(*r.seq.get(qp + k as usize)?));
                    if !(c1 == c2 && c1 != 15) {
                        nm += 1;
                    }
                }
                rp += len;
                qp += len as usize;
            }
            Op::Ins => {
                nm += len as i64;
                qp += len as usize;
            }
            Op::Del => {
                nm += len as i64;
                rp += len;
            }
            Op::Skip => rp += len,
            Op::SoftClip => qp += len as usize,
            _ => {}
        }
    }
    Some(nm)
}

fn summary(a: &Acc) -> Value {
    // insert sizes: halve per size (each pair was counted twice), then the main 99 % bulk
    let half = |v: &Vec<u64>| -> Vec<u64> { v.iter().map(|&x| (x as f64 * 0.5) as u64).collect() };
    let (inw, out, oth) = (half(&a.inward), half(&a.outward), half(&a.other));
    let total: Vec<u64> = (0..inw.len()).map(|i| inw[i] + out[i] + oth[i]).collect();
    let mut nisize: u64 = total.iter().sum();
    let mut ibulk = 0usize;
    let (mut bulk, mut avg) = (0f64, 0f64);
    for (i, &num) in total.iter().enumerate() {
        if num > 0 {
            ibulk = i + 1;
        }
        bulk += num as f64;
        avg += i as f64 * num as f64;
        if bulk / nisize as f64 > MAIN_BULK {
            ibulk = i + 1;
            nisize = bulk as u64;
            break;
        }
    }
    let avg_isize = avg / if nisize > 0 { nisize as f64 } else { 1.0 };
    // median from the printed IS lines (0 .. ibulk-1), as the 1.1 read QC does
    let printed: u64 = total[..ibulk].iter().sum();
    let mut median = None;
    let mut seen = 0u64;
    for (i, &n) in total[..ibulk].iter().enumerate() {
        seen += n;
        if printed > 0 && seen as f64 >= printed as f64 / 2.0 {
            median = Some(i);
            break;
        }
    }
    let raw = a.reads_1st + a.reads_2nd + a.reads_other;
    let mapped = a.mapped_paired + a.mapped_single;
    let pct = |x: u64, y: u64| if y > 0 { format!("{:.2}", 100.0 * x as f64 / y as f64).parse::<f64>().ok() } else { None };
    // samtools prints the error rate with %e (float); 1.1 rounds the printed value to 5 decimals
    let err = if a.bases_mapped_cigar > 0 { (a.mismatches as f32 / a.bases_mapped_cigar as f32) as f64 } else { 0.0 };
    let err_printed: f64 = format!("{err:e}").parse().unwrap_or(err);
    let err_printed: f64 = format!("{:.6e}", err_printed).parse().unwrap_or(err);
    let avg_len = if raw > 0 { (a.total_len as f32 / raw as f32) as f64 } else { 0.0 };
    let avg_qual = if a.total_len > 0 { a.sum_qual / a.total_len as f64 } else { 0.0 };
    json!({
        "readsInRegions": raw,
        "mappedPct": pct(mapped, raw),
        "properlyPairedPct": pct(a.properly_paired, mapped),
        "duplicatePct": pct(a.duplicated, mapped),
        "mapq0Pct": pct(a.mq0, mapped),
        "medianInsertSize": median,
        "meanInsertSize": format!("{avg_isize:.1}").parse::<f64>().ok(),
        "errorRate": format!("{:.5}", err_printed).parse::<f64>().ok(),
        "averageReadLength": format!("{avg_len:.0}").parse::<f64>().ok(),
        "averageBaseQuality": format!("{avg_qual:.1}").parse::<f64>().ok(),
    })
}
