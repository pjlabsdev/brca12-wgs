//! HGVS c. and p. names and the variant consequence on the MANE Select transcripts of BRCA1
//! (NM_007294.4) and BRCA2 (NM_000059.4), from the bundled transcript models and the reference
//! FASTA. Deletions, duplications and insertions are shifted 3' in transcript direction (HGVS 3'
//! rule), protein changes are named with three-letter amino acids as Ensembl VEP does (Ter for the
//! stop codon). Research use only.
use std::io;

const MODELS: &str = include_str!("../assets/mane_brca12.tsv");
const PAD: u64 = 6000; // reference kept around the transcript (3' shift, upstream/downstream)

pub struct Tx {
    pub gene: String,
    pub refseq: String,
    pub ensembl: String,
    pub protein: &'static str,
    pub chrom: String,
    pub plus: bool,
    exons: Vec<(u64, u64)>, // genomic, 1-based inclusive, in transcript order
    cds: (u64, u64),        // genomic min..max including the stop codon
    seq: Vec<u8>,           // plus-strand reference from seq_off
    seq_off: u64,           // genomic 1-based position of seq[0]
    cdna: Vec<u8>,          // transcript sequence (5' to 3')
    cs: i64,                // cDNA index (1-based) of the A of ATG
    ce: i64,                // cDNA index of the last base of the stop codon
}

#[derive(Clone, Debug, Default)]
pub struct Annotation {
    pub gene: String,
    pub hgvs_c: Option<String>,
    pub hgvs_p: Option<String>,
    pub consequence: String,
    pub exon: Option<String>,   // "11/23"
    pub intron: Option<String>, // "10/22"
}

fn comp(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        x => x,
    }
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| comp(b)).collect()
}

fn aa3(a: u8) -> &'static str {
    match a {
        b'A' => "Ala", b'R' => "Arg", b'N' => "Asn", b'D' => "Asp", b'C' => "Cys", b'Q' => "Gln", b'E' => "Glu",
        b'G' => "Gly", b'H' => "His", b'I' => "Ile", b'L' => "Leu", b'K' => "Lys", b'M' => "Met", b'F' => "Phe",
        b'P' => "Pro", b'S' => "Ser", b'T' => "Thr", b'W' => "Trp", b'Y' => "Tyr", b'V' => "Val", b'*' => "Ter",
        _ => "Xaa",
    }
}

fn aas(s: &[u8]) -> String {
    s.iter().map(|&a| aa3(a)).collect()
}

fn codon(c: &[u8]) -> u8 {
    let idx = |b: u8| match b {
        b'T' => Some(0),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    };
    const T: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    match (idx(c[0]), idx(c[1]), idx(c[2])) {
        (Some(a), Some(b), Some(d)) => T[a * 16 + b * 4 + d],
        _ => b'X',
    }
}

/// translation from the first base until the first stop codon (included) or the end
fn translate(s: &[u8]) -> Vec<u8> {
    let mut p = vec![];
    for c in s.chunks(3) {
        if c.len() < 3 {
            break;
        }
        let a = codon(c);
        p.push(a);
        if a == b'*' {
            break;
        }
    }
    p
}

/// coding consequence terms (Sequence Ontology, as VEP)
const CODING: [&str; 11] = ["stop_gained", "frameshift_variant", "stop_lost", "start_lost", "inframe_insertion", "inframe_deletion",
    "missense_variant", "protein_altering_variant", "synonymous_variant", "stop_retained_variant", "start_retained_variant"];

/// most severe first (Ensembl VEP order)
const SEVERITY: [&str; 22] = ["splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant", "stop_lost",
    "start_lost", "inframe_insertion", "inframe_deletion", "missense_variant", "protein_altering_variant",
    "splice_donor_5th_base_variant", "splice_region_variant", "splice_donor_region_variant", "splice_polypyrimidine_tract_variant",
    "start_retained_variant", "stop_retained_variant", "synonymous_variant", "coding_sequence_variant", "5_prime_UTR_variant",
    "3_prime_UTR_variant", "intron_variant", "upstream_gene_variant"];

/// position on the transcript: cDNA index of an exonic base, or the nearest exon base and the
/// distance into the intron; None outside the transcript
#[derive(Clone, Copy, Debug, PartialEq)]
struct TPos {
    t: i64,
    off: i64,
}

impl Tx {
    pub fn load(reference: &str) -> io::Result<Vec<Tx>> {
        let mut fa = noodles::fasta::io::indexed_reader::Builder::default().build_from_path(reference)?;
        let mut out = vec![];
        for (gene, protein, chroms) in [("BRCA1", "NP_009225.1", ["chr17", "17"]), ("BRCA2", "NP_000050.3", ["chr13", "13"])] {
            let rows: Vec<Vec<&str>> =
                MODELS.lines().filter(|l| !l.starts_with('#')).map(|l| l.split('\t').collect::<Vec<_>>()).filter(|f| f[0] == gene).collect();
            let plus = rows[0][3] == "+";
            let mut exons: Vec<(u64, u64)> = rows.iter().filter(|f| f[4] == "exon").map(|f| (f[5].parse().unwrap(), f[6].parse().unwrap())).collect();
            exons.sort();
            if !plus {
                exons.reverse();
            }
            let cds_rows: Vec<(u64, u64)> = rows
                .iter()
                .filter(|f| f[4] == "CDS" || f[4] == "stop_codon")
                .map(|f| (f[5].parse().unwrap(), f[6].parse().unwrap()))
                .collect();
            let cds = (cds_rows.iter().map(|x| x.0).min().unwrap(), cds_rows.iter().map(|x| x.1).max().unwrap());
            let lo = exons.iter().map(|e| e.0).min().unwrap().saturating_sub(PAD).max(1);
            let hi = exons.iter().map(|e| e.1).max().unwrap() + PAD;
            let mut seq = None;
            let mut chrom = String::new();
            for c in chroms {
                if let Ok(region) = format!("{c}:{lo}-{hi}").parse::<noodles::core::Region>() {
                    if let Ok(r) = fa.query(&region) {
                        seq = Some(r.sequence().as_ref().iter().map(|b| b.to_ascii_uppercase()).collect::<Vec<u8>>());
                        chrom = c.to_string();
                        break;
                    }
                }
            }
            let seq = seq.ok_or_else(|| io::Error::other(format!("{gene}: chromosome not in the reference")))?;
            let mut tx = Tx {
                gene: gene.into(),
                refseq: rows[0][1].into(),
                ensembl: rows[0][2].into(),
                protein,
                chrom,
                plus,
                exons,
                cds,
                seq,
                seq_off: lo,
                cdna: vec![],
                cs: 0,
                ce: 0,
            };
            let mut cdna = vec![];
            for &(s, e) in &tx.exons {
                let part = tx.gseq(s, e);
                cdna.extend(if plus { part } else { revcomp(&part) });
            }
            tx.cdna = cdna;
            let (a, b) = if plus { (cds.0, cds.1) } else { (cds.1, cds.0) };
            tx.cs = tx.tpos(a).unwrap().t;
            tx.ce = tx.tpos(b).unwrap().t;
            let p = translate(&tx.cdna[(tx.cs - 1) as usize..tx.ce as usize]);
            if p.first() != Some(&b'M') || p.last() != Some(&b'*') || p.len() as i64 * 3 != tx.ce - tx.cs + 1 {
                return Err(io::Error::other(format!("{gene}: transcript model does not translate (wrong reference?)")));
            }
            out.push(tx);
        }
        Ok(out)
    }

    /// plus-strand reference bases s..=e (1-based)
    fn gseq(&self, s: u64, e: u64) -> Vec<u8> {
        if e < s {
            return vec![];
        }
        self.seq[(s - self.seq_off) as usize..=(e - self.seq_off) as usize].to_vec()
    }

    fn base(&self, g: u64) -> u8 {
        self.seq.get((g.wrapping_sub(self.seq_off)) as usize).copied().unwrap_or(b'N')
    }

    pub fn span(&self) -> (u64, u64) {
        (self.exons.iter().map(|e| e.0).min().unwrap(), self.exons.iter().map(|e| e.1).max().unwrap())
    }

    fn tpos(&self, g: u64) -> Option<TPos> {
        let mut t = 0i64;
        for (k, &(s, e)) in self.exons.iter().enumerate() {
            if g >= s && g <= e {
                let inside = if self.plus { g - s } else { e - g } as i64;
                return Some(TPos { t: t + inside + 1, off: 0 });
            }
            let len = (e - s + 1) as i64;
            if let Some(&(ns, ne)) = self.exons.get(k + 1) {
                // intron between exon k and k+1 (transcript order)
                let (in_intron, up_d, down_d) =
                    if self.plus { (g > e && g < ns, g as i64 - e as i64, ns as i64 - g as i64) } else { (g < s && g > ne, s as i64 - g as i64, g as i64 - ne as i64) };
                if in_intron {
                    return Some(if up_d <= down_d { TPos { t: t + len, off: up_d } } else { TPos { t: t + len + 1, off: -down_d } });
                }
            }
            t += len;
        }
        None
    }

    /// distance of an exonic cDNA position to the nearest exon border that faces an intron
    fn exon_border_distance(&self, t: i64) -> Option<i64> {
        let mut start = 1i64;
        let n = self.exons.len();
        for (k, &(s, e)) in self.exons.iter().enumerate() {
            let end = start + (e - s) as i64;
            if t >= start && t <= end {
                let mut d: Option<i64> = None;
                if k > 0 { d = Some(t - start); }
                if k + 1 < n { d = Some(d.map_or(end - t, |x| x.min(end - t))); }
                return d;
            }
            start = end + 1;
        }
        None
    }

    fn cnum(&self, p: TPos) -> String {
        let base = if p.t < self.cs {
            format!("-{}", self.cs - p.t)
        } else if p.t > self.ce {
            format!("*{}", p.t - self.ce)
        } else {
            format!("{}", p.t - self.cs + 1)
        };
        match p.off {
            0 => base,
            o if o > 0 => format!("{base}+{o}"),
            o => format!("{base}{o}"),
        }
    }

    /// exon or intron number (k of n) of a genomic position
    fn feature(&self, g: u64) -> (Option<String>, Option<String>) {
        let n = self.exons.len();
        for (k, &(s, e)) in self.exons.iter().enumerate() {
            if g >= s && g <= e {
                return (Some(format!("{}/{n}", k + 1)), None);
            }
            if let Some(&(ns, ne)) = self.exons.get(k + 1) {
                if (self.plus && g > e && g < ns) || (!self.plus && g < s && g > ne) {
                    return (None, Some(format!("{}/{}", k + 1, n - 1)));
                }
            }
        }
        (None, None)
    }

    /// annotation of one VCF allele (1-based POS, REF, ALT)
    pub fn annotate(&self, pos: u64, r: &[u8], a: &[u8]) -> Option<Annotation> {
        let (span_s, span_e) = self.span();
        let mut r = r.to_ascii_uppercase();
        let mut a = a.to_ascii_uppercase();
        let mut pos = pos;
        // minimal representation: common suffix, then common prefix
        while r.len() > 0 && a.len() > 0 && r.last() == a.last() {
            r.pop();
            a.pop();
        }
        while r.len() > 0 && a.len() > 0 && r[0] == a[0] {
            r.remove(0);
            a.remove(0);
            pos += 1;
        }
        if r.is_empty() && a.is_empty() {
            return None;
        }
        // deleted genomic bases s..=e (e = s - 1 for an insertion before s)
        let mut s = pos;
        let mut e = if r.is_empty() { s - 1 } else { pos + r.len() as u64 - 1 };
        if e + 1 + 2 * PAD / 3 < span_s || s > span_e + 2 * PAD / 3 {
            return None;
        }
        let mut ins = a.clone();
        // consequence and protein change come from the left-aligned (unshifted) allele, as in VEP
        let (s0, e0) = (s, e);
        // HGVS 3' rule in transcript direction
        if r.is_empty() || a.is_empty() {
            if self.plus {
                if a.is_empty() {
                    while self.base(e + 1) == self.base(s) && e + 1 < self.seq_off + self.seq.len() as u64 {
                        s += 1;
                        e += 1;
                    }
                } else {
                    while self.base(s) == ins[0] && s < self.seq_off + self.seq.len() as u64 {
                        ins.rotate_left(1);
                        s += 1;
                        e += 1;
                    }
                }
            } else if a.is_empty() {
                while s > self.seq_off && self.base(s - 1) == self.base(e) {
                    s -= 1;
                    e -= 1;
                }
            } else {
                while s > self.seq_off && self.base(s - 1) == *ins.last().unwrap() {
                    ins.rotate_right(1);
                    s -= 1;
                    e -= 1;
                }
            }
        }
        let gene = self.gene.clone();
        // feature and consequence ranges
        let outside = |g: u64| g < span_s || g > span_e;
        let (lo_g, hi_g) = if e >= s { (s, e) } else { (e, s) }; // insertion: flanks e and s
        if outside(lo_g) && outside(hi_g) && (hi_g < span_s || lo_g > span_e) {
            let up = if self.plus { hi_g < span_s } else { lo_g > span_e };
            return Some(Annotation { gene, consequence: if up { "upstream_gene_variant" } else { "downstream_gene_variant" }.into(), ..Default::default() });
        }
        // c. name
        let tp = |g: u64| self.tpos(g);
        let kind;
        let mut name = String::new();
        let tr = |x: &[u8]| -> String { String::from_utf8_lossy(&if self.plus { x.to_vec() } else { revcomp(x) }).into_owned() };
        let range = |g1: u64, g2: u64| -> Option<String> {
            let (a1, a2) = if self.plus { (g1, g2) } else { (g2, g1) };
            let (p1, p2) = (tp(a1)?, tp(a2)?);
            Some(if g1 == g2 { self.cnum(p1) } else { format!("{}_{}", self.cnum(p1), self.cnum(p2)) })
        };
        if a.is_empty() {
            kind = "del";
            if let Some(x) = range(s, e) {
                name = format!("{x}del");
            }
        } else if r.is_empty() {
            let l = ins.len() as u64;
            // duplication: the inserted bases equal the bases just 5' of the insertion (transcript)
            let dup = if self.plus { s >= l + self.seq_off && self.gseq(s - l, s - 1) == ins } else { self.gseq(s, s + l - 1) == ins };
            if dup {
                kind = "dup";
                let (g1, g2) = if self.plus { (s - l, s - 1) } else { (s, s + l - 1) };
                if let Some(x) = range(g1, g2) {
                    name = format!("{x}dup");
                }
            } else {
                kind = "ins";
                if let (Some(p1), Some(p2)) = if self.plus { (tp(s - 1), tp(s)) } else { (tp(s), tp(s - 1)) } {
                    name = format!("{}_{}ins{}", self.cnum(p1), self.cnum(p2), tr(&ins));
                }
            }
        } else if r.len() == 1 && a.len() == 1 {
            kind = "snv";
            if let Some(p) = tp(s) {
                name = format!("{}{}>{}", self.cnum(p), tr(&r), tr(&a));
            }
        } else {
            kind = "delins";
            if let Some(x) = range(s, e) {
                // an inversion: ALT is the reverse complement of REF (HGVS inv)
                name = if r.len() > 1 && revcomp(&r) == a { format!("{x}inv") } else { format!("{x}delins{}", tr(&a)) };
            }
        }
        let (exon, intron) = {
            let (e1, i1) = self.feature(lo_g.max(span_s).min(span_e));
            (e1, i1)
        };
        let kind0 = if r.is_empty() { "ins" } else if a.is_empty() { "del" } else if r.len() == 1 && a.len() == 1 { "snv" } else { "delins" };
        // positions on the transcript touched by the unshifted allele (an insertion: its two flanks)
        let affected: Vec<Option<TPos>> = if kind0 == "ins" { vec![tp(s0 - 1), tp(s0)] } else { (s0..=e0).map(tp).collect() };
        let in_tx: Vec<TPos> = affected.iter().flatten().cloned().collect();
        let cds_ok = |p: &TPos| p.off == 0 && p.t >= self.cs && p.t <= self.ce;
        // coding: both ends in the coding exons (a deletion may span whole exons); an insertion
        // needs a coding base on its 5' side (transcript direction), as Ensembl VEP
        let coding_at = |k: &str, a: u64, b: u64| -> bool {
            if k == "ins" {
                tp(if self.plus { a - 1 } else { a }).is_some_and(|p| cds_ok(&p))
            } else {
                tp(a).is_some_and(|p| cds_ok(&p)) && tp(b).is_some_and(|p| cds_ok(&p))
            }
        };
        let kind_s = if kind == "dup" { "ins" } else { kind };
        let mut hgvs_p = None;
        let mut terms: Vec<&'static str> = vec![];
        // p. from the 3'-shifted allele, the consequence from the unshifted one (as Ensembl VEP)
        let shifted = if coding_at(kind_s, s, e) { Some(self.protein_change(kind_s, s, e, &ins)) } else { None };
        // the p. name needs the unshifted allele in the coding exons or in an intron next to them
        // (the protein is the same); not in the UTR, where a shift into the start codon would
        // move the start (VEP gives no p. for a few intronic cases, see docs/output-2.md)
        let utr0 = (if kind0 == "ins" { vec![tp(s0 - 1), tp(s0)] } else { vec![tp(s0), tp(e0)] })
            .iter()
            .any(|p| p.is_none_or(|p| p.off == 0 && (p.t < self.cs || p.t > self.ce)));
        if let Some((pn, _)) = &shifted {
            if coding_at(kind0, s0, e0) || !utr0 {
                hgvs_p = pn.clone();
            }
        }
        if coding_at(kind0, s0, e0) {
            terms.extend(match &shifted {
                Some((_, t)) => t.clone(),
                None => self.protein_change(kind0, s0, e0, &a).1,
            });
        }
        let offs: Vec<i64> = if kind0 == "ins" {
            // an insertion counts at an intron position only when both flanks are in the intron
            if in_tx.len() == 2 && in_tx[0].off != 0 && in_tx[1].off != 0 { vec![in_tx[0].off, in_tx[1].off] } else { vec![] }
        } else {
            in_tx.iter().filter(|p| p.off != 0).map(|p| p.off).collect()
        };
        let core = |o: i64| o.abs() <= 2;
        let donor = if kind0 == "ins" { offs.len() == 2 && offs.iter().all(|&o| o > 0 && core(o)) } else { offs.iter().any(|&o| o > 0 && core(o)) };
        let acceptor = if kind0 == "ins" { offs.len() == 2 && offs.iter().all(|&o| o < 0 && core(o)) } else { offs.iter().any(|&o| o < 0 && core(o)) };
        if donor { terms.push("splice_donor_variant"); }
        if acceptor { terms.push("splice_acceptor_variant"); }

        // exonic bases next to an intron (1-3 bases) and intron bases 3-8 form the splice region
        let near_border = in_tx.iter().filter(|p| p.off == 0).any(|p| self.exon_border_distance(p.t).is_some_and(|d| d < 3));
        let ins_at_border = kind0 == "ins" && in_tx.len() == 2 && (in_tx[0].off == 0) != (in_tx[1].off == 0);
        // an insertion lies in a region only with both flanks in it
        let has = |f: &dyn Fn(i64) -> bool| if kind0 == "ins" { !offs.is_empty() && offs.iter().all(|&o| f(o)) } else { offs.iter().any(|&o| f(o)) };
        let fifth = kind0 != "ins" && has(&|o| o == 5);
        let donor_region = kind0 != "ins" && has(&|o| (3..=6).contains(&o));
        if fifth { terms.push("splice_donor_5th_base_variant"); }
        if !(donor || acceptor || fifth || donor_region) && (near_border || ins_at_border || has(&|o| o.abs() >= 3 && o.abs() <= 8)) {
            terms.push("splice_region_variant");
        }
        if donor_region && !fifth { terms.push("splice_donor_region_variant"); }
        if has(&|o| (-17..=-3).contains(&o)) { terms.push("splice_polypyrimidine_tract_variant"); }
        let ex: Vec<&TPos> = in_tx.iter().filter(|p| p.off == 0).collect();
        if !terms.iter().any(|t| CODING.contains(t)) && ex.iter().any(|p| cds_ok(p)) {
            terms.push("coding_sequence_variant");
        }
        if ex.iter().any(|p| p.t < self.cs) { terms.push("5_prime_UTR_variant"); }
        if ex.iter().any(|p| p.t > self.ce) { terms.push("3_prime_UTR_variant"); }
        if offs.iter().any(|&o| o.abs() > 2) { terms.push("intron_variant"); }
        if terms.is_empty() {
            terms.push("intron_variant");
        }
        terms.sort_by_key(|t| SEVERITY.iter().position(|x| x == t).unwrap_or(99));
        terms.dedup();
        let consequence = terms.join(",");
        Some(Annotation {
            gene,
            hgvs_c: if name.is_empty() { None } else { Some(format!("{}:c.{name}", self.refseq)) },
            hgvs_p: hgvs_p.map(|p| format!("{}:p.{p}", self.protein)),
            consequence,
            exon,
            intron,
        })
    }

    /// p. name and coding consequence terms of a change inside the coding exons (unshifted allele)
    fn protein_change(&self, kind: &str, s: u64, e: u64, alt: &[u8]) -> (Option<String>, Vec<&'static str>) {
        let _ = CODING;
        // the change on the cDNA: 0-based half-open [t0, t1) replaced by alt_t
        let alt_t = if self.plus { alt.to_vec() } else { revcomp(alt) };
        let (t0, t1) = if kind == "ins" {
            // insertion between genomic s-1 and s: bases of the transcript 5' of it
            let (f1, f2) = if self.plus { (self.tpos(s - 1).unwrap(), self.tpos(s).unwrap()) } else { (self.tpos(s).unwrap(), self.tpos(s - 1).unwrap()) };
            let t = if f1.off == 0 { f1.t } else { f2.t - 1 } as usize;
            (t, t)
        } else {
            let (a, b) = if self.plus { (s, e) } else { (e, s) };
            (self.tpos(a).unwrap().t as usize - 1, self.tpos(b).unwrap().t as usize)
        };
        let cs0 = (self.cs - 1) as usize;
        let refp = translate(&self.cdna[cs0..self.ce as usize]);
        let mut alt_cdna = self.cdna[..t0].to_vec();
        alt_cdna.extend(&alt_t);
        alt_cdna.extend(&self.cdna[t1..]);
        let altp = translate(&alt_cdna[cs0..]);
        let diff = alt_t.len() as i64 - (t1 - t0) as i64;
        // codons touched in the reference: [c0, c1) (an insertion at a codon border touches none)
        let c0 = (t0 - cs0) / 3;
        let c1 = if t1 > t0 { (t1 - cs0).div_ceil(3) } else if (t0 - cs0) % 3 == 0 { c0 } else { c0 + 1 };
        if c0 == 0 && (c1 > 0 || kind == "ins" && t0 > cs0) && altp.first() != Some(&b'M') {
            return (Some("Met1?".into()), vec!["start_lost"]);
        }
        if c0 == 0 && t1 > cs0 && diff != 0 && altp.first() == Some(&b'M') && t0 < cs0 + 3 {
            return (Some("Met1?".into()), vec!["start_lost"]);
        }
        let i = refp.iter().zip(altp.iter()).take_while(|(x, y)| x == y).count();
        if i >= refp.len() {
            let k = c0.min(refp.len() - 1);
            let term = if refp[k] == b'*' { "stop_retained_variant" } else { "synonymous_variant" };
            return (Some(format!("{}{}=", aa3(refp[k]), k + 1)), vec![term]);
        }
        if refp[i] == b'*' {
            // the stop codon is lost: extension to the next stop
            let a = altp.get(i).copied().unwrap_or(b'X');
            let ter = altp[i..].iter().skip(1).position(|&x| x == b'*').map(|j| (j + 1).to_string()).unwrap_or("?".into());
            let mut t = vec!["stop_lost"];
            if diff.rem_euclid(3) != 0 { t.insert(0, "frameshift_variant"); }
            return (Some(format!("Ter{}{}extTer{}", i + 1, aa3(a), ter)), t);
        }
        if diff.rem_euclid(3) != 0 {
            // stop_gained when the changed codons themselves read a stop (VEP)
            let w0 = c0 * 3 + cs0;
            let wl = ((c1 - c0) as i64 * 3 + diff).max(0) as usize;
            let wl = wl / 3 * 3; // complete codons only
            let win_stop = alt_cdna.get(w0..(w0 + wl).min(alt_cdna.len())).is_some_and(|w| w.chunks(3).filter(|c| c.len() == 3).any(|c| codon(c) == b'*'));
            let mut t = vec!["frameshift_variant"];
            if win_stop { t.push("stop_gained"); }
            if altp.get(i) == Some(&b'*') {
                return (Some(format!("{}{}Ter", aa3(refp[i]), i + 1)), t);
            }
            let ter = altp[i..].iter().position(|&x| x == b'*').map(|j| (j + 1).to_string()).unwrap_or("?".into());
            let a = altp.get(i).copied().unwrap_or(b'X');
            return (Some(format!("{}{}{}fsTer{}", aa3(refp[i]), i + 1, aa3(a), ter)), t);
        }
        // in frame: the touched codons in reference and alternative, common ends trimmed (VEP)
        let alt_end = (c1 * 3 + cs0) as i64 + diff;
        let rw: Vec<u8> = refp[c0..c1.min(refp.len())].to_vec();
        let aw: Vec<u8> = if alt_end as usize >= c0 * 3 + cs0 {
            alt_cdna[c0 * 3 + cs0..alt_end as usize].chunks(3).filter(|c| c.len() == 3).map(codon).collect()
        } else {
            vec![]
        };
        if aw.contains(&b'*') {
            let (mut lo, mut rh, mut ah) = (0usize, rw.len(), aw.len());
            while lo < rh && lo < ah && rw[lo] == aw[lo] { lo += 1; }
            while rh > lo && ah > lo && rw[rh - 1] == aw[ah - 1] { rh -= 1; ah -= 1; }
            let rr = &rw[lo..rh];
            let mut xx = aw[lo..ah].to_vec();
            if let Some(k) = xx.iter().position(|&x| x == b'*') { xx.truncate(k + 1); }
            let at = c0 + lo; // 0-based index of the first reference amino acid in rr
            let name = if rr.is_empty() {
                format!("{}{}_{}{}ins{}", aa3(refp[at - 1]), at, aa3(refp[at]), at + 1, aas(&xx))
            } else if rr.len() == 1 && xx.len() == 1 && diff == 0 {
                format!("{}{}{}", aa3(rr[0]), at + 1, aa3(xx[0]))
            } else if rr.len() == 1 {
                format!("{}{}delins{}", aa3(rr[0]), at + 1, aas(&xx))
            } else {
                format!("{}{}_{}{}delins{}", aa3(rr[0]), at + 1, aa3(*rr.last().unwrap()), at + rr.len(), aas(&xx))
            };
            let _ = (rh, ah);
            return (Some(name), vec!["stop_gained"]);
        }
        // otherwise the whole proteins, which shifts deletions and insertions 3' (HGVS)
        let mut j = 0;
        while j < refp.len() - i && j < altp.len() - i && refp[refp.len() - 1 - j] == altp[altp.len() - 1 - j] {
            j += 1;
        }
        let rr = &refp[i..refp.len() - j];
        let xx = &altp[i..altp.len() - j];
        let (name, term) = if rr.len() == 1 && xx.len() == 1 {
            (format!("{}{}{}", aa3(rr[0]), i + 1, aa3(xx[0])), "missense_variant")
        } else if xx.is_empty() {
            let name = if rr.len() == 1 {
                format!("{}{}del", aa3(rr[0]), i + 1)
            } else {
                format!("{}{}_{}{}del", aa3(rr[0]), i + 1, aa3(*rr.last().unwrap()), i + rr.len())
            };
            (name, "inframe_deletion")
        } else if rr.is_empty() {
            let l = xx.len();
            if i >= l && &refp[i - l..i] == xx {
                let name = if l == 1 { format!("{}{}dup", aa3(xx[0]), i) } else { format!("{}{}_{}{}dup", aa3(refp[i - l]), i - l + 1, aa3(refp[i - 1]), i) };
                (name, "inframe_insertion")
            } else {
                (format!("{}{}_{}{}ins{}", aa3(refp[i - 1]), i, aa3(refp[i]), i + 1, aas(xx)), "inframe_insertion")
            }
        } else {
            let name = if rr.len() == 1 {
                format!("{}{}delins{}", aa3(rr[0]), i + 1, aas(xx))
            } else {
                format!("{}{}_{}{}delins{}", aa3(rr[0]), i + 1, aa3(*rr.last().unwrap()), i + rr.len(), aas(xx))
            };
            (name, if diff < 0 { "inframe_deletion" } else if diff == 0 { "missense_variant" } else { "protein_altering_variant" })
        };
        (Some(name), vec![term])
    }
}

/// the annotation on the transcript whose gene region holds the variant (None outside both)
pub fn annotate(txs: &[Tx], chrom: &str, pos: u64, r: &[u8], a: &[u8]) -> Option<Annotation> {
    let c = chrom.trim_start_matches("chr");
    let on: Vec<&Tx> = txs.iter().filter(|t| t.chrom.trim_start_matches("chr") == c).collect();
    if let Some(x) = on.iter().filter_map(|t| t.annotate(pos, r, a)).next() {
        return Some(x);
    }
    // farther than the kept reference from the transcript: up- or downstream of the nearest one
    let t = on.iter().min_by_key(|t| {
        let (s, e) = t.span();
        if pos < s { s - pos } else { pos.saturating_sub(e) }
    })?;
    let (s, e) = t.span();
    let up = if t.plus { pos < s } else { pos > e };
    Some(Annotation {
        gene: t.gene.clone(),
        consequence: if up { "upstream_gene_variant" } else { "downstream_gene_variant" }.into(),
        ..Default::default()
    })
}
