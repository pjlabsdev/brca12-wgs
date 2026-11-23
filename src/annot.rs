//! ClinVar and gnomAD for the BRCA1/BRCA2 regions, carried in the binary (built by
//! tools/build_annotation_tables.py): exact allele match after trimming REF/ALT to the minimal
//! representation. Research use only.
use flate2::read::GzDecoder;
use std::collections::HashMap;
use std::io::Read as _;

const CLINVAR: &[u8] = include_bytes!("../assets/clinvar_brca12.tsv.gz");
const GNOMAD: &[u8] = include_bytes!("../assets/gnomad_brca12.tsv.gz");

type Key = (String, u64, Vec<u8>, Vec<u8>);

#[derive(Clone, Debug)]
pub struct ClinVar {
    pub id: String,
    pub significance: String, // CLNSIG, e.g. Pathogenic, Likely_pathogenic, Conflicting_classifications_of_pathogenicity
    pub review: String,       // CLNREVSTAT
}

impl ClinVar {
    /// pathogenic or likely pathogenic (not conflicting)
    pub fn plp(&self) -> bool {
        let s = self.significance.to_ascii_lowercase();
        (s.starts_with("pathogenic") || s.starts_with("likely_pathogenic")) && !s.contains("conflicting")
    }
    pub fn benign(&self) -> bool {
        let s = self.significance.to_ascii_lowercase();
        (s.starts_with("benign") || s.starts_with("likely_benign")) && !s.contains("conflicting")
    }
    /// review stars as ClinVar shows them
    pub fn stars(&self) -> u8 {
        match self.review.as_str() {
            "practice_guideline" => 4,
            "reviewed_by_expert_panel" => 3,
            "criteria_provided,_multiple_submitters,_no_conflicts" => 2,
            "criteria_provided,_single_submitter" | "criteria_provided,_conflicting_classifications" => 1,
            _ => 0,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Gnomad {
    pub dataset: String,
    pub present: bool,
    pub faf95: Option<f64>,
    pub population: String,
    pub non_cancer: bool,
}

pub struct Tables {
    pub clinvar_version: String,
    pub gnomad_version: String,
    clinvar: HashMap<Key, ClinVar>,
    gnomad: HashMap<Key, Vec<Gnomad>>,
}

fn text(gz: &[u8]) -> String {
    let mut s = String::new();
    GzDecoder::new(gz).read_to_string(&mut s).expect("bundled table");
    s
}

/// minimal representation of an allele (common suffix, then common prefix removed)
pub fn key(chrom: &str, pos: u64, r: &[u8], a: &[u8]) -> Key {
    let mut r = r.to_ascii_uppercase();
    let mut a = a.to_ascii_uppercase();
    let mut pos = pos;
    while r.len() > 1 && a.len() > 1 && r.last() == a.last() {
        r.pop();
        a.pop();
    }
    while r.len() > 1 && a.len() > 1 && r[0] == a[0] {
        r.remove(0);
        a.remove(0);
        pos += 1;
    }
    let c = if chrom.starts_with("chr") { chrom.to_string() } else { format!("chr{chrom}") };
    (c, pos, r, a)
}

impl Tables {
    pub fn load() -> Tables {
        let cv = text(CLINVAR);
        let gn = text(GNOMAD);
        let clinvar_version = cv.lines().next().unwrap_or("").trim_start_matches("# ").to_string();
        let gnomad_version = gn.lines().next().unwrap_or("").trim_start_matches("# ").to_string();
        let mut clinvar = HashMap::new();
        for l in cv.lines().filter(|l| !l.starts_with('#')) {
            let f: Vec<&str> = l.split('\t').collect();
            if f.len() < 7 {
                continue;
            }
            let k = key(f[0], f[1].parse().unwrap_or(0), f[2].as_bytes(), f[3].as_bytes());
            clinvar.insert(k, ClinVar { id: f[4].into(), significance: f[5].into(), review: f[6].into() });
        }
        let mut gnomad: HashMap<Key, Vec<Gnomad>> = HashMap::new();
        for l in gn.lines().filter(|l| !l.starts_with('#')) {
            let f: Vec<&str> = l.split('\t').collect();
            if f.len() < 9 {
                continue;
            }
            let k = key(f[0], f[1].parse().unwrap_or(0), f[2].as_bytes(), f[3].as_bytes());
            gnomad.entry(k).or_default().push(Gnomad {
                dataset: f[4].into(),
                present: f[5] == "1",
                faf95: f[6].parse().ok(),
                population: f[7].into(),
                non_cancer: f[8] == "1",
            });
        }
        Tables { clinvar_version, gnomad_version, clinvar, gnomad }
    }

    pub fn clinvar(&self, chrom: &str, pos: u64, r: &[u8], a: &[u8]) -> Option<&ClinVar> {
        self.clinvar.get(&key(chrom, pos, r, a))
    }

    /// highest FAF95 over the datasets (as 1.1), and whether the allele is present in any of them
    pub fn gnomad(&self, chrom: &str, pos: u64, r: &[u8], a: &[u8]) -> (bool, Option<&Gnomad>) {
        match self.gnomad.get(&key(chrom, pos, r, a)) {
            None => (false, None),
            Some(v) => {
                let present = v.iter().any(|g| g.present);
                let best = v.iter().filter(|g| g.faf95.is_some()).max_by(|x, y| x.faf95.partial_cmp(&y.faf95).unwrap());
                (present, best)
            }
        }
    }
}
