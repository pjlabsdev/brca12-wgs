//! Paired FASTQ input: the read pairs of BRCA1/BRCA2 are pulled out of a whole-genome FASTQ and
//! aligned inside the program, after hla-wgs. Research use only.
//!
//! 1. Targets: the gene regions +/- 2 kb (BRCA1 with the pseudogene BRCA1P1) and the control windows
//!    of the exon copy-number panel.
//! 2. One pass over the reference gives the genome-unique 31-mers and the decoys. Genome-unique 31-mers: every 31-mer of the targets that occurs exactly once in the whole
//!    reference (both strands). Reads of repeats from elsewhere (Alu, L1) do not share such a k-mer
//!    and are not pulled in, which replaces the mapping quality a whole-genome alignment would give.
//!    Decoys: loci elsewhere that share many low-copy (2-20 copies) 31-mers with the targets
//!    (pseudogenes, segmental duplications); they join the minimap2 panel, and reads that fit a decoy
//!    best are dropped, as a whole-genome aligner would place them there. Both are cached next to
//!    the reference (REF.brca12-kmers) after the first run.
//! 3. A pair is kept when one of its reads has at least two of these k-mers; kept pairs are aligned
//!    with minimap2 (short-read preset, paired; compiled into the program) to the targets, placed
//!    on GRCh38, duplicates marked as samtools markdup does for pairs, sorted and indexed. A read
//!    without two genome-unique 31-mers of its own, or whose mate does not align to the targets, gets
//!    mapping quality 0 (it may come from a copy elsewhere in the genome).
//! The result is OUTPREFIX.bam (+ .bai), which the analysis then reads like any BAM.
use crate::exon_cnv;
use ahash::{AHashMap, AHashSet};
use minimap2::Aligner;
use noodles::sam;
use rayon::prelude::*;
use crate::util::{Fq, PairedFq};
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

fn hits(ms: &[minimap2::Mapping], targets: &[Target], by_name: &AHashMap<String, usize>, tids: &AHashMap<String, usize>) -> Vec<Hit> {
    // a read whose primary alignment is on a decoy belongs elsewhere in the genome
    if ms.iter().any(|m| m.is_primary && m.target_name.as_ref().is_some_and(|n| n.starts_with("decoy|"))) {
        return vec![];
    }
    ms.iter()
        .filter(|m| m.is_primary || m.is_supplementary)
        .filter_map(|m| {
            let t = &targets[*by_name.get(m.target_name.as_ref()?.as_str())?];
            let a = m.alignment.as_ref()?;
            let cigar = parse_cigar(a.cigar_str.as_deref()?);
            let pos = t.start + m.target_start as u64;
            Some(Hit {
                tid: tids[&t.chrom],
                pos,
                end: t.start + m.target_end as u64,
                rev: m.strand == minimap2::Strand::Reverse,
                mapq: m.mapq.min(60) as u8,
                cigar,
                supp: m.is_supplementary,
                nm: a.nm,
            })
        })
        .collect()
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
    // decoys (paralogs elsewhere) join the minimap2 panel, so that their reads find their own locus
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
    // minimap2 index from a FASTA of the targets (its in-memory index works for one sequence only)
    let panel = format!("{prefix}.targets.fa");
    {
        let mut f = std::io::BufWriter::new(std::fs::File::create(&panel)?);
        for (n, seq) in tg.iter().map(|t| (&t.name, &t.seq)).chain(decoy_seqs.iter().map(|d| (&d.0, &d.1))) {
            writeln!(f, ">{n}")?;
            for c in seq.chunks(80) {
                f.write_all(c)?;
                f.write_all(b"\n")?;
            }
        }
    }
    let al = Aligner::builder()
        .sr()
        .with_cigar()
        .with_index_threads(threads)
        .with_index(&panel, None)
        .map_err(io::Error::other)?;
    let _ = std::fs::remove_file(&panel);
    // target order in the index = order in the FASTA
    let by_name: AHashMap<String, usize> = tg.iter().enumerate().map(|(i, t)| (t.name.clone(), i)).collect();
    // the paired reader of hla-wgs feeds batches (it survives I/O errors: pause, reopen, continue);
    // k-mer screening and alignment run in parallel
    let (tx, rx) = std::sync::mpsc::sync_channel::<(Vec<Fq>, Vec<Fq>)>(4);
    let (p1, p2) = (r1.to_string(), r2.to_string());
    let reader = std::thread::spawn(move || -> io::Result<u64> {
        let mut src = PairedFq::new(&p1, &p2);
        let mut n = 0u64;
        loop {
            let (b1, b2) = src.batch(200_000);
            if b1.is_empty() {
                break;
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
                let Ok((m1, m2)) = al.map_pair(&x.seq, &y.seq, false, false, None, None, Some(x.name.as_bytes())) else { return vec![] };
                let (mut h1, mut h2) = (hits(&m1, &tg, &by_name, &tids), hits(&m2, &tg, &by_name, &tids));
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
                pair_records(x.name.as_bytes(), x, y, &h1, &h2, &names)
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
