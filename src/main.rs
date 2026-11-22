//! brca12-wgs 2.0 (in development): BRCA1/BRCA2 germline analysis from whole-genome sequencing in
//! one program, without external tools. Research use only (RUO).
mod coverage;
mod reads;

use std::io::{self, BufRead};

const USAGE: &str = "brca12-wgs 2.0.0-dev (research use only)

usage:
  brca12-wgs coverage <BAM|CRAM> [--reference REF.fa] [--exons EXONS.bed] [--min-mapq 20]
      exon depth like mosdepth (MAPQ filter, flags 0x704 skipped, mate overlap once):
      name, mean depth and bases >= 10/20/30x per exon, tab-separated on stdout
  brca12-wgs depth <BAM|CRAM> <chrom:start-end> [--reference REF.fa] [--min-mapq 20] [--no-mate-overlap]
      per-base depth (1-based positions), like mosdepth per-base or (--no-mate-overlap) samtools depth -a
";

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter().position(|a| a == key).and_then(|i| args.get(i + 1)).cloned()
}

fn read_bed(path: &str) -> io::Result<Vec<(String, u64, u64, String)>> {
    let mut out = vec![];
    for line in io::BufReader::new(std::fs::File::open(path)?).lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') || line.starts_with("track") {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 {
            continue;
        }
        let name = f.get(3).copied().unwrap_or("").to_string();
        out.push((f[0].to_string(), f[1].parse().unwrap_or(0), f[2].parse().unwrap_or(0), name));
    }
    Ok(out)
}

fn cmd_coverage(args: &[String]) -> io::Result<()> {
    let input = args.first().filter(|a| !a.starts_with("--")).cloned().ok_or_else(|| io::Error::other("input BAM/CRAM missing"))?;
    let reference = arg(args, "--reference");
    let exons = arg(args, "--exons").unwrap_or_else(|| "assets/brca12_exons.bed".into());
    let min_mapq: u8 = arg(args, "--min-mapq").and_then(|v| v.parse().ok()).unwrap_or(20);
    let bed = read_bed(&exons)?;
    let flt = coverage::Filter { min_mapq, ..coverage::MOSDEPTH };
    // one window per contig over all its exons: each read is fetched once
    let mut chroms: Vec<String> = bed.iter().map(|b| b.0.clone()).collect();
    chroms.dedup();
    for chrom in chroms {
        let rows: Vec<_> = bed.iter().filter(|b| b.0 == chrom).collect();
        let s = rows.iter().map(|b| b.1).min().unwrap();
        let e = rows.iter().map(|b| b.2).max().unwrap();
        let d = coverage::depth(&input, reference.as_deref(), &chrom, s, e, flt)?;
        for b in rows {
            let st = coverage::interval_stats(&chrom, b.1, b.2, &b.3, &d, s, &[10, 20, 30]);
            let th: Vec<String> = st.at_least.iter().map(|(_, n)| n.to_string()).collect();
            println!("{}\t{}\t{}\t{}\t{:.2}\t{}", st.chrom, st.start, st.end, st.name, st.mean, th.join("\t"));
        }
    }
    Ok(())
}

fn cmd_depth(args: &[String]) -> io::Result<()> {
    let input = args.first().cloned().ok_or_else(|| io::Error::other("input BAM/CRAM missing"))?;
    let region = args.get(1).cloned().ok_or_else(|| io::Error::other("region chrom:start-end missing"))?;
    let (chrom, range) = region.rsplit_once(':').ok_or_else(|| io::Error::other("region must be chrom:start-end"))?;
    let (s, e) = range.split_once('-').ok_or_else(|| io::Error::other("region must be chrom:start-end"))?;
    let (s, e): (u64, u64) = (s.replace(',', "").parse().map_err(io::Error::other)?, e.replace(',', "").parse().map_err(io::Error::other)?);
    let min_mapq: u8 = arg(args, "--min-mapq").and_then(|v| v.parse().ok()).unwrap_or(20);
    let mut flt = coverage::Filter { min_mapq, ..coverage::MOSDEPTH };
    if args.iter().any(|a| a == "--no-mate-overlap") {
        flt.mate_overlap_once = false;
    }
    let d = coverage::depth(&input, arg(args, "--reference").as_deref(), chrom, s - 1, e, flt)?;
    for (i, v) in d.iter().enumerate() {
        println!("{chrom}\t{}\t{v}", s + i as u64);
    }
    Ok(())
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let res = match args.first().map(|s| s.as_str()) {
        Some("coverage") => cmd_coverage(&args[1..]),
        Some("depth") => cmd_depth(&args[1..]),
        Some("--version") | Some("-V") => {
            println!("brca12-wgs {}", env!("CARGO_PKG_VERSION"));
            Ok(())
        }
        _ => {
            eprint!("{USAGE}");
            std::process::exit(2);
        }
    };
    if let Err(e) = res {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}
