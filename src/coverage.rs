//! Per-base depth over a window, with the read filters of mosdepth (default) and samtools depth.
//! mosdepth: MAPQ >= min, records with flags 0x704 (unmapped, secondary, QC fail, duplicate)
//! skipped, aligned bases (M/=/X) counted, and the overlap of the two mates of a properly paired
//! fragment counted once. samtools depth: the same flags skipped, no mate-overlap correction. Research use only.
use crate::reads::{self, Op, Read};
use std::collections::HashMap;
use std::io;

#[derive(Clone, Copy, Debug)]
pub struct Filter {
    pub min_mapq: u8,
    pub exclude: u16,
    pub mate_overlap_once: bool,
}

pub const MOSDEPTH: Filter = Filter { min_mapq: 0, exclude: 0x704, mate_overlap_once: true };
pub const SAMTOOLS_DEPTH: Filter = Filter { min_mapq: 0, exclude: 0x704, mate_overlap_once: false };

/// aligned reference blocks [s, e) of a read (M/=/X only; deletions and skips are not covered)
pub fn blocks(r: &Read) -> Vec<(u64, u64)> {
    let mut out = vec![];
    let mut pos = r.start;
    for &(op, len) in &r.cigar {
        let len = len as u64;
        match op {
            Op::Match => {
                if let Some(last) = out.last_mut() {
                    let (_, e): &mut (u64, u64) = last;
                    if *e == pos {
                        *e = pos + len;
                        pos += len;
                        continue;
                    }
                }
                out.push((pos, pos + len));
                pos += len;
            }
            Op::Del | Op::Skip => pos += len,
            _ => {}
        }
    }
    out
}

/// Depth for positions [start, end) (0-based) of contig `chrom`.
pub fn depth(path: &str, reference: Option<&str>, chrom: &str, start: u64, end: u64, flt: Filter) -> io::Result<Vec<u32>> {
    let mut d = vec![0u32; (end - start) as usize];
    let mut add = |s: u64, e: u64, d: &mut Vec<u32>| {
        let (s, e) = (s.max(start), e.min(end));
        for p in s..e {
            d[(p - start) as usize] += 1;
        }
    };
    // fragments whose two mates may overlap: blocks kept until both mates are seen
    let mut open: HashMap<Vec<u8>, Vec<(u64, u64)>> = HashMap::new();
    let region = format!("{chrom}:{}-{}", start + 1, end);
    reads::fetch(path, reference, &region, &mut |r: &Read| {
        if r.flags & flt.exclude != 0 {
            return;
        }
        if r.mapq < flt.min_mapq {
            return;
        }
        let b = blocks(r);
        // mosdepth's mate-overlap rule, reproduced so the numbers match 1.1: only properly paired
        // reads (flag 0x2) take part. A read whose mate starts inside it (start <= mate start <
        // end) is remembered under its name; when the next record with that name comes, the bases
        // both cover are counted once. Mates with the same start are corrected too.
        let paired_same_contig = flt.mate_overlap_once
            && r.flags & reads::PROPER_PAIR != 0
            && r.mate_tid == r.tid
            && r.mate_start >= 0;
        if paired_same_contig {
            if let Some(first) = open.remove(&r.name) {
                for (s, e) in b {
                    // count only the bases of this read that the remembered one does not cover
                    let mut cur = s;
                    for &(fs, fe) in &first {
                        if fe <= cur || fs >= e {
                            continue;
                        }
                        if fs > cur {
                            add(cur, fs.min(e), &mut d);
                        }
                        cur = cur.max(fe);
                        if cur >= e {
                            break;
                        }
                    }
                    if cur < e {
                        add(cur, e, &mut d);
                    }
                }
                return;
            }
        }
        for &(s, e) in &b {
            add(s, e, &mut d);
        }
        if paired_same_contig && r.start <= r.mate_start as u64 && (r.mate_start as u64) < r.end {
            open.insert(r.name.clone(), b);
        }
    })?;
    // remembered mates are already counted; nothing to add for mates outside the window
    drop(open);
    Ok(d)
}

#[derive(Clone, Debug, serde::Serialize)]
pub struct IntervalStats {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub name: String,
    pub mean: f64,
    pub at_least: Vec<(u32, u64)>, // threshold, bases
}

pub fn interval_stats(chrom: &str, start: u64, end: u64, name: &str, d: &[u32], d_start: u64, thresholds: &[u32]) -> IntervalStats {
    let slice = &d[(start - d_start) as usize..(end - d_start) as usize];
    let sum: u64 = slice.iter().map(|&x| x as u64).sum();
    IntervalStats {
        chrom: chrom.to_string(),
        start,
        end,
        name: name.to_string(),
        mean: if slice.is_empty() { 0.0 } else { sum as f64 / slice.len() as f64 },
        at_least: thresholds.iter().map(|&t| (t, slice.iter().filter(|&&x| x >= t).count() as u64)).collect(),
    }
}
