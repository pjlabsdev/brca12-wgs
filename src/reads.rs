//! Alignment input: indexed BAM or CRAM, read per region through the index (no whole-file pass).
//! Every record is handed on as a `Read` with the fields the analyses need. Research use only.
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record::cigar::op::Kind;
use std::io;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Op {
    Match,     // M, =, X: consumes reference and read
    Ins,       // I
    Del,       // D
    Skip,      // N
    SoftClip,  // S
    HardClip,  // H
    Pad,       // P
}

#[derive(Clone, Debug)]
pub struct Read {
    pub name: Vec<u8>,
    pub flags: u16,
    pub tid: i32,
    pub start: u64, // 0-based alignment start
    pub end: u64,   // 0-based exclusive alignment end
    pub mapq: u8,
    pub cigar: Vec<(Op, u32)>,
    pub mate_tid: i32,
    pub mate_start: i64, // 0-based, -1 if none
    pub tlen: i64,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    pub nm: Option<i64>, // NM tag (edit distance)
    pub sa: Option<String>, // SA tag (supplementary alignments), as in SAM
}

pub const UNMAPPED: u16 = 0x4;
pub const SECONDARY: u16 = 0x100;
pub const QCFAIL: u16 = 0x200;
pub const DUPLICATE: u16 = 0x400;
pub const SUPPLEMENTARY: u16 = 0x800;
pub const PAIRED: u16 = 0x1;
pub const PROPER_PAIR: u16 = 0x2;
pub const MATE_UNMAPPED: u16 = 0x8;
pub const REVERSE: u16 = 0x10;
pub const MATE_REVERSE: u16 = 0x20;
pub const FIRST: u16 = 0x40;

fn err<E: std::fmt::Display>(e: E) -> io::Error {
    io::Error::other(e.to_string())
}

fn convert(r: &RecordBuf, tid_of: &dyn Fn(Option<usize>) -> i32) -> Option<Read> {
    let start = r.alignment_start().map(|p| p.get() as u64 - 1);
    let cigar: Vec<(Op, u32)> = r
        .cigar()
        .as_ref()
        .iter()
        .map(|o| {
            let k = match o.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => Op::Match,
                Kind::Insertion => Op::Ins,
                Kind::Deletion => Op::Del,
                Kind::Skip => Op::Skip,
                Kind::SoftClip => Op::SoftClip,
                Kind::HardClip => Op::HardClip,
                Kind::Pad => Op::Pad,
            };
            (k, o.len() as u32)
        })
        .collect();
    let start = start?;
    let span: u64 = cigar.iter().filter(|(k, _)| matches!(k, Op::Match | Op::Del | Op::Skip)).map(|&(_, l)| l as u64).sum();
    Some(Read {
        name: r.name().map(|n| n.to_vec()).unwrap_or_default(),
        flags: u16::from(r.flags()),
        tid: tid_of(r.reference_sequence_id()),
        start,
        end: start + span,
        mapq: r.mapping_quality().map(|m| m.get()).unwrap_or(255),
        cigar,
        mate_tid: tid_of(r.mate_reference_sequence_id()),
        mate_start: r.mate_alignment_start().map(|p| p.get() as i64 - 1).unwrap_or(-1),
        tlen: r.template_length() as i64,
        seq: r.sequence().as_ref().to_vec(),
        qual: r.quality_scores().as_ref().to_vec(),
        nm: r.data().get(&noodles::sam::alignment::record::data::field::Tag::EDIT_DISTANCE).and_then(|v| v.as_int()),
        sa: r.data().get(&noodles::sam::alignment::record::data::field::Tag::OTHER_ALIGNMENTS).and_then(|v| match v {
            noodles::sam::alignment::record_buf::data::field::Value::String(s) => Some(String::from_utf8_lossy(s.as_ref()).into_owned()),
            _ => None,
        }),
    })
}

/// Contig names and lengths from the header of an alignment file.
pub fn contigs(path: &str, reference: Option<&str>) -> io::Result<Vec<(String, u64)>> {
    let mut out = vec![];
    with_reader(path, reference, &mut |hdr, _| {
        out = hdr.reference_sequences().iter().map(|(n, m)| (n.to_string(), m.length().get() as u64)).collect();
        Ok(false)
    })?;
    Ok(out)
}

/// Raw header text (for @RG and @PG).
pub fn header_text(path: &str, reference: Option<&str>) -> io::Result<String> {
    let mut text = String::new();
    with_reader(path, reference, &mut |hdr, _| {
        let mut w = noodles::sam::io::Writer::new(Vec::new());
        w.write_header(hdr)?;
        text = String::from_utf8_lossy(&w.into_inner()).into_owned();
        Ok(false)
    })?;
    Ok(text)
}

enum Rd {
    Bam(noodles::bam::io::Reader<noodles::bgzf::io::Reader<std::fs::File>>, noodles::bam::Index),
    Cram(noodles::cram::io::Reader<std::fs::File>, noodles::cram::crai::Index),
}

fn open(path: &str, reference: Option<&str>) -> io::Result<(Rd, noodles::sam::Header)> {
    if path.to_ascii_lowercase().ends_with(".cram") {
        use noodles::{cram, fasta};
        let repo = match reference {
            Some(src) => fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(
                fasta::io::indexed_reader::Builder::default().build_from_path(src).map_err(err)?,
            )),
            None => fasta::Repository::default(),
        };
        let mut rd = cram::io::reader::Builder::default().set_reference_sequence_repository(repo).build_from_path(path)?;
        let header = rd.read_header()?;
        let index = cram::fs::read_associated_index(path)?;
        Ok((Rd::Cram(rd, index), header))
    } else {
        use noodles::bam;
        let mut rd = bam::io::reader::Builder::default().build_from_path(path)?;
        let header = rd.read_header()?;
        let index = bam::fs::read_associated_index(path)
            .or_else(|_| {
                let alt = path.strip_suffix(".bam").map(|b| format!("{b}.bai")).unwrap_or_default();
                bam::bai::fs::read(&alt).map(bam::Index::Bai)
            })
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("index .bam.bai/.bai/.csi needed next to {path}: {e}")))?;
        Ok((Rd::Bam(rd, index), header))
    }
}

fn with_reader(path: &str, reference: Option<&str>, f: &mut dyn FnMut(&noodles::sam::Header, &mut Rd) -> io::Result<bool>) -> io::Result<()> {
    let (mut rd, header) = open(path, reference)?;
    f(&header, &mut rd)?;
    Ok(())
}

/// Hand every record overlapping `region` ("chr17:43044295-43170327", 1-based inclusive) to `f`.
pub fn fetch(path: &str, reference: Option<&str>, region: &str, f: &mut dyn FnMut(&Read)) -> io::Result<u64> {
    let (rd, header) = open(path, reference)?;
    let tid_of = |id: Option<usize>| id.map(|i| i as i32).unwrap_or(-1);
    let region: noodles::core::Region = region.parse().map_err(err)?;
    let mut n = 0u64;
    match rd {
        Rd::Bam(mut rd, index) => {
            for r in rd.query(&header, &index, &region)?.records() {
                let r = r?;
                let buf = RecordBuf::try_from_alignment_record(&header, &r)?;
                if let Some(read) = convert(&buf, &tid_of) {
                    n += 1;
                    f(&read);
                }
            }
        }
        Rd::Cram(mut rd, index) => {
            for r in rd.query(&header, &index, &region)?.records() {
                let r = r?;
                let buf = RecordBuf::try_from_alignment_record(&header, &r)?;
                if let Some(read) = convert(&buf, &tid_of) {
                    n += 1;
                    f(&read);
                }
            }
        }
    }
    Ok(n)
}
