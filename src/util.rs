//! Input helpers taken over from hla-wgs (src/util.rs): FASTQ reading that survives I/O errors
//! (pause, reopen, skip what was delivered), plain error exit. The FASTQ record also keeps the base
//! qualities here, which the variant caller needs.
use flate2::read::MultiGzDecoder;
use std::io::{BufRead, BufReader};

pub fn open_reader(p: &str) -> Box<dyn BufRead> {
    let f = std::fs::File::open(p).unwrap_or_else(|e| panic!("{p}: {e}"));
    if p.ends_with(".gz") {
        Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(f)))
    } else {
        Box::new(BufReader::with_capacity(1 << 20, f))
    }
}

/// end the run with a plain message (no panic trace) and exit code 1
pub fn die(msg: &str) -> ! {
    eprintln!("error: {msg}");
    std::process::exit(1)
}

/// errors that a pause cannot cure: the content itself is wrong
pub fn permanent(e: &std::io::Error) -> bool {
    matches!(e.kind(), std::io::ErrorKind::InvalidData | std::io::ErrorKind::InvalidInput)
}

fn bad(msg: String) -> std::io::Error {
    std::io::Error::new(std::io::ErrorKind::InvalidData, msg)
}

pub struct Fq {
    pub name: String,
    pub seq: Vec<u8>,
    /// base qualities (phred, without the +33 offset)
    pub qual: Vec<u8>,
}

/// how often an I/O error on the input is retried (env BRCA12_IO_RETRIES, default 6)
pub fn io_retries() -> u32 {
    std::env::var("BRCA12_IO_RETRIES").ok().and_then(|v| v.parse().ok()).unwrap_or(6)
}

/// Test hook: env BRCA12_TEST_IO_FAULT=N makes the first read pass fail once after N records
/// (BAM/CRAM) or pairs (FASTQ), so that the retry path can be checked without pulling a disk.
pub fn test_fault(n: u64) -> bool {
    use std::sync::atomic::{AtomicBool, Ordering};
    static FIRED: AtomicBool = AtomicBool::new(false);
    match std::env::var("BRCA12_TEST_IO_FAULT").ok().and_then(|v| v.parse::<u64>().ok()) {
        Some(at) if n >= at && !FIRED.swap(true, Ordering::SeqCst) => true,
        _ => false,
    }
}

/// pause before retry k (1-based): 5, 10, 20, 40, 80, 160 s, at most 300 s
pub fn backoff(k: u32) -> u64 {
    if std::env::var("BRCA12_TEST_IO_FAULT").is_ok() {
        return 1;
    }
    (5u64 << (k.saturating_sub(1)).min(6)).min(300)
}

/// next n FASTQ records (name = header up to first space, /1 /2 stripped); I/O errors are returned
pub fn read_fq_res(r: &mut dyn BufRead, n: usize) -> std::io::Result<Vec<Fq>> {
    let mut out = Vec::with_capacity(n);
    let mut l = [Vec::new(), Vec::new(), Vec::new(), Vec::new()];
    for _ in 0..n {
        for (i, x) in l.iter_mut().enumerate() {
            x.clear();
            if r.read_until(b'\n', x)? == 0 {
                if i == 0 {
                    return Ok(out);
                }
                return Err(bad(format!("FASTQ ends in the middle of a record (record {}); the file is truncated", out.len() + 1)));
            }
            while x.last().map_or(false, |c| *c == b'\n' || *c == b'\r') {
                x.pop();
            }
            if i == 0 && x.first() != Some(&b'@') {
                return Err(bad(format!("not a FASTQ file: record {} does not start with '@'", out.len() + 1)));
            }
        }
        if l[2].first() != Some(&b'+') {
            return Err(bad(format!("not a FASTQ record (record {}): the third line must start with '+'", out.len() + 1)));
        }
        let h = &l[0][1..];
        let e = h.iter().position(|&c| c == b' ' || c == b'\t').unwrap_or(h.len());
        let mut nm = &h[..e];
        if nm.len() > 2 && nm[nm.len() - 2] == b'/' {
            nm = &nm[..nm.len() - 2];
        }
        let seq: Vec<u8> = l[1].iter().map(|b| b.to_ascii_uppercase()).collect();
        let qual: Vec<u8> = l[3].iter().map(|q| q.saturating_sub(33)).collect();
        if qual.len() != seq.len() {
            return Err(bad(format!("FASTQ record {}: sequence and quality have different lengths", out.len() + 1)));
        }
        out.push(Fq { name: String::from_utf8_lossy(nm).into_owned(), seq, qual });
    }
    Ok(out)
}

/// open_reader that reports a missing or unreadable file instead of aborting
pub fn try_open_reader(p: &str) -> std::io::Result<Box<dyn BufRead>> {
    let f = std::fs::File::open(p)?;
    Ok(if p.ends_with(".gz") {
        Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(f)))
    } else {
        Box::new(BufReader::with_capacity(1 << 20, f))
    })
}

/// Paired FASTQ reader that survives I/O errors: after a pause both files are reopened and the
/// records already delivered are skipped, so the caller sees the same sequence of batches as
/// in an uninterrupted read.
pub struct PairedFq {
    paths: [String; 2],
    rd: Option<[Box<dyn BufRead>; 2]>,
    done: u64,
}

impl PairedFq {
    pub fn new(r1: &str, r2: &str) -> Self {
        PairedFq { paths: [r1.to_string(), r2.to_string()], rd: None, done: 0 }
    }

    fn try_batch(&mut self, n: usize) -> std::io::Result<(Vec<Fq>, Vec<Fq>)> {
        if self.rd.is_none() {
            let mut a = try_open_reader(&self.paths[0])?;
            let mut b = try_open_reader(&self.paths[1])?;
            let mut left = self.done;
            while left > 0 {
                let k = left.min(50_000) as usize;
                let (x, y) = (read_fq_res(&mut *a, k)?, read_fq_res(&mut *b, k)?);
                if x.len() != k || y.len() != k {
                    return Err(std::io::Error::other("input shorter than before the I/O error"));
                }
                left -= k as u64;
            }
            self.rd = Some([a, b]);
        }
        let [a, b] = self.rd.as_mut().unwrap();
        let b1 = read_fq_res(&mut **a, n)?;
        let b2 = read_fq_res(&mut **b, n)?;
        if b1.len() != b2.len() {
            return Err(bad(format!("R1 and R2 hold different numbers of reads (after {} pairs)", self.done + b1.len().min(b2.len()) as u64)));
        }
        self.done += b1.len() as u64;
        if test_fault(self.done) {
            self.done -= b1.len() as u64;
            return Err(std::io::Error::other("injected test fault (BRCA12_TEST_IO_FAULT)"));
        }
        Ok((b1, b2))
    }

    /// next batch of up to n pairs; empty at the end of the input
    pub fn batch(&mut self, n: usize) -> (Vec<Fq>, Vec<Fq>) {
        let tries = io_retries();
        let mut attempt = 0;
        loop {
            match self.try_batch(n) {
                Ok(x) => return x,
                Err(e) => {
                    if permanent(&e) || e.kind() == std::io::ErrorKind::UnexpectedEof {
                        die(&format!("{} / {}: {e}", self.paths[0], self.paths[1]));
                    }
                    self.rd = None;
                    attempt += 1;
                    if attempt > tries {
                        die(&format!("{}: {e} (gave up after {tries} retries)", self.paths[0]));
                    }
                    let wait = backoff(attempt);
                    eprintln!("[input] {}: I/O error after {} pairs: {e}; retry {attempt}/{tries} in {wait}s", self.paths[0], self.done);
                    std::thread::sleep(std::time::Duration::from_secs(wait));
                }
            }
        }
    }
}

pub fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&c| match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'a' => b't',
            b'c' => b'g',
            b'g' => b'c',
            b't' => b'a',
            x => x,
        })
        .collect()
}
