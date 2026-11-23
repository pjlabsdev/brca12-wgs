//! brca12-wgs 2.0 (in development): BRCA1/BRCA2 germline analysis from whole-genome sequencing in
//! one program, without external tools. Research use only (RUO).
mod annot;
mod assemble;
mod caller;
mod coverage;
mod exon_cnv;
mod exons;
mod hgvs;
mod output;
mod reads;
mod readstats;
mod sv;

use std::io::{self, BufRead};

const USAGE: &str = "brca12-wgs 2.0.0-dev (research use only)

usage:
  brca12-wgs INPUT.bam|INPUT.cram OUTPREFIX [--reference REF.fa] [--exons BED] [--regions BED]
      the analysis: writes OUTPREFIX.tsv (one line per gene), OUTPREFIX.variants.tsv (one line per
      small variant with HGVS, ClinVar and gnomAD), OUTPREFIX.exons.tsv (one line per exon with its
      certainty) and OUTPREFIX.json (schema genome.brca12.v2); see docs/output-2.md
  brca12-wgs coverage <BAM|CRAM> [--reference REF.fa] [--exons EXONS.bed] [--min-mapq 20]
      exon depth like mosdepth (MAPQ filter, flags 0x704 skipped, mate overlap once):
      name, mean depth and bases >= 10/20/30x per exon, tab-separated on stdout
  brca12-wgs exons <BAM|CRAM> [--reference REF.fa] [--exons EXONS.bed] [--regions REGIONS.bed]
      exon table (depth, >= 10/20/30x, share of reads with MAPQ < 20, dosage against the chromosome
      median, status), chromosome medians and low-coverage intervals as JSON (1.1 schema)
  brca12-wgs readqc <BAM|CRAM> [--reference REF.fa] [--regions REGIONS.bed] [--padding 200]
      read statistics over the calling regions (gene regions +/- padding) like samtools stats
  brca12-wgs call <BAM|CRAM> OUT.vcf --reference REF.fa [--regions BED] [--padding 200] [--sample NAME]
      SNVs and indels in the gene regions (+/- padding) with the own caller, as VCF
  brca12-wgs annotate <VCF> --reference REF.fa
      HGVS c./p. and consequence on the MANE Select transcripts for every allele of a VCF
  brca12-wgs depth <BAM|CRAM> <chrom:start-end> [--reference REF.fa] [--min-mapq 20] [--no-mate-overlap]
      per-base depth (1-based positions), like mosdepth per-base or (--no-mate-overlap) samtools depth -a
";

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter().position(|a| a == key).and_then(|i| args.get(i + 1)).cloned()
}

/// exon and gene-region BEDs carried in the binary (used when --exons / --regions are not given)
const EXONS_BED: &str = include_str!("../assets/brca12_exons.bed");
const REGIONS_BED: &str = include_str!("../assets/brca12_regions.bed");

fn read_bed(path: &str) -> io::Result<Vec<(String, u64, u64, String)>> {
    let text = match path {
        "assets/brca12_exons.bed" => EXONS_BED.to_string(),
        "assets/brca12_regions.bed" => REGIONS_BED.to_string(),
        p => std::fs::read_to_string(p).map_err(|e| io::Error::other(format!("{p}: {e}")))?,
    };
    let mut out = vec![];
    for line in text.lines() {
        let line = line.to_string();
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

fn cmd_exons(args: &[String]) -> io::Result<()> {
    let input = args.first().filter(|a| !a.starts_with("--")).cloned().ok_or_else(|| io::Error::other("input BAM/CRAM missing"))?;
    let reference = arg(args, "--reference");
    let exons = read_bed(&arg(args, "--exons").unwrap_or_else(|| "assets/brca12_exons.bed".into()))?;
    let regions = read_bed(&arg(args, "--regions").unwrap_or_else(|| "assets/brca12_regions.bed".into()))?;
    let contigs = reads::contigs(&input, reference.as_deref())?;
    let len = |c: &str| contigs.iter().find(|(n, _)| n == c).map(|&(_, l)| l);
    let v = exons::table(&input, reference.as_deref(), &exons, &regions, &len, &exons::DEFAULT)?;
    println!("{}", serde_json::to_string_pretty(&v).map_err(io::Error::other)?);
    Ok(())
}

fn cmd_readqc(args: &[String]) -> io::Result<()> {
    let input = args.first().filter(|a| !a.starts_with("--")).cloned().ok_or_else(|| io::Error::other("input BAM/CRAM missing"))?;
    let reference = arg(args, "--reference");
    let pad: u64 = arg(args, "--padding").and_then(|v| v.parse().ok()).unwrap_or(200);
    let mut regs: Vec<(String, u64, u64)> = read_bed(&arg(args, "--regions").unwrap_or_else(|| "assets/brca12_regions.bed".into()))?
        .into_iter()
        .map(|(c, s, e, _)| (c, s.saturating_sub(pad), e + pad))
        .collect();
    regs.sort();
    let v = readstats::region_stats(&input, reference.as_deref(), &regs)?;
    println!("{}", serde_json::to_string_pretty(&v).map_err(io::Error::other)?);
    Ok(())
}

fn cmd_run(args: &[String]) -> io::Result<()> {
    let pos: Vec<&String> = {
        let mut v = vec![];
        let mut i = 0;
        while i < args.len() {
            if args[i].starts_with("--") { i += 2; continue; }
            v.push(&args[i]);
            i += 1;
        }
        v
    };
    let (input, prefix) = match pos.as_slice() {
        [i, p, ..] => (i.as_str(), p.as_str()),
        _ => return Err(io::Error::other("usage: brca12-wgs INPUT.bam|INPUT.cram OUTPREFIX [--reference REF.fa]")),
    };
    let reference = arg(args, "--reference");
    let exons_bed = read_bed(&arg(args, "--exons").unwrap_or_else(|| "assets/brca12_exons.bed".into()))?;
    let regions = read_bed(&arg(args, "--regions").unwrap_or_else(|| "assets/brca12_regions.bed".into()))?;
    if let Some(dir) = std::path::Path::new(prefix).parent() {
        if !dir.as_os_str().is_empty() { std::fs::create_dir_all(dir)?; }
    }
    let t0 = std::time::Instant::now();
    output::run(&output::Run { input, reference: reference.as_deref(), prefix, exons: &exons_bed, regions: &regions,
                               calling_padding: 200, settings: &exons::DEFAULT })?;
    eprintln!("brca12-wgs: {prefix}.tsv, {prefix}.variants.tsv, {prefix}.exons.tsv, {prefix}.json ({:.1} s)", t0.elapsed().as_secs_f64());
    Ok(())
}

fn cmd_call(args: &[String]) -> io::Result<()> {
    let input = args.first().cloned().ok_or_else(|| io::Error::other("input BAM/CRAM missing"))?;
    let out = args.get(1).cloned().ok_or_else(|| io::Error::other("output VCF missing"))?;
    let reference = arg(args, "--reference").ok_or_else(|| io::Error::other("--reference REF.fa is needed"))?;
    let pad: u64 = arg(args, "--padding").and_then(|v| v.parse().ok()).unwrap_or(200);
    let regions = read_bed(&arg(args, "--regions").unwrap_or_else(|| "assets/brca12_regions.bed".into()))?;
    let sample = arg(args, "--sample").unwrap_or_else(|| {
        std::path::Path::new(&input).file_stem().and_then(|s| s.to_str()).unwrap_or("sample").to_string()
    });
    let mut regs: Vec<_> = regions.iter().map(|(c, s, e, _)| (c.clone(), s.saturating_sub(pad), e + pad)).collect();
    regs.sort();
    let t0 = std::time::Instant::now();
    let mut calls = vec![];
    for (c, s, e) in &regs {
        calls.extend(caller::call_region(&input, &reference, c, *s, *e, &caller::DEFAULT)?);
    }
    let contigs = reads::contigs(&input, Some(&reference))?;
    caller::write_vcf(&out, &sample, &contigs, &calls)?;
    eprintln!("brca12-wgs call: {} variants -> {out} ({:.1} s)", calls.len(), t0.elapsed().as_secs_f64());
    Ok(())
}

/// HGVS c./p. and consequence for every allele of a VCF (validation against Ensembl VEP)
fn cmd_annotate(args: &[String]) -> io::Result<()> {
    let input = args.first().cloned().ok_or_else(|| io::Error::other("input VCF missing"))?;
    let reference = arg(args, "--reference").ok_or_else(|| io::Error::other("--reference REF.fa is needed"))?;
    let txs = hgvs::Tx::load(&reference)?;
    println!("id\tchrom\tpos\tref\talt\tgene\texon\tintron\tconsequence\thgvs_c\thgvs_p");
    for line in io::BufReader::new(std::fs::File::open(&input)?).lines() {
        let line = line?;
        if line.starts_with('#') { continue; }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 5 { continue; }
        let pos: u64 = f[1].parse().map_err(io::Error::other)?;
        for alt in f[4].split(',') {
            let a = hgvs::annotate(&txs, f[0], pos, f[3].as_bytes(), alt.as_bytes()).unwrap_or_default();
            let d = |o: &Option<String>| o.clone().unwrap_or("-".into());
            println!("{}\t{}\t{}\t{}\t{alt}\t{}\t{}\t{}\t{}\t{}\t{}", f[2], f[0], f[1], f[3], if a.gene.is_empty() { "-" } else { &a.gene },
                     d(&a.exon), d(&a.intron), if a.consequence.is_empty() { "-" } else { &a.consequence }, d(&a.hgvs_c), d(&a.hgvs_p));
        }
    }
    Ok(())
}

/// exon copy number alone (depths and the HMM result as JSON); --depths-only prints the depths
fn cmd_exon_cnv(args: &[String]) -> io::Result<()> {
    let input = args.first().filter(|a| !a.starts_with("--")).cloned().ok_or_else(|| io::Error::other("input BAM/CRAM missing"))?;
    let p = exon_cnv::panel();
    let d = exon_cnv::depths(&input, arg(args, "--reference").as_deref(), &p)?;
    if args.iter().any(|a| a == "--depths-only") {
        let mut keys: Vec<&String> = d.keys().collect();
        keys.sort();
        for k in keys {
            println!("{k}\t{}", d[k].map(|x| format!("{x:.6}")).unwrap_or("-".into()));
        }
        return Ok(());
    }
    println!("{}", serde_json::to_string_pretty(&exon_cnv::call(&p, &d, None)).map_err(io::Error::other)?);
    Ok(())
}

/// structural-variant and insertion signatures of the gene regions as JSON (1.1 field names)
fn cmd_sv(args: &[String]) -> io::Result<()> {
    let input = args.first().filter(|a| !a.starts_with("--")).cloned().ok_or_else(|| io::Error::other("input BAM/CRAM missing"))?;
    let reference = arg(args, "--reference").ok_or_else(|| io::Error::other("--reference REF.fa is needed"))?;
    let regions = read_bed(&arg(args, "--regions").unwrap_or_else(|| "assets/brca12_regions.bed".into()))?;
    let v = sv::signals(&input, &reference, &regions)?;
    println!("{}", serde_json::to_string_pretty(&v).map_err(io::Error::other)?);
    Ok(())
}

/// debug: bases of every read at one position (flags, strand, MAPQ, base, base quality)
fn cmd_pileup(args: &[String]) -> io::Result<()> {
    let input = args.first().cloned().ok_or_else(|| io::Error::other("input BAM/CRAM missing"))?;
    let pos = args.get(1).cloned().ok_or_else(|| io::Error::other("chrom:pos missing"))?;
    let (chrom, p) = pos.rsplit_once(':').ok_or_else(|| io::Error::other("chrom:pos"))?;
    let p: u64 = p.parse::<u64>().map_err(io::Error::other)? - 1;
    let mut counts: std::collections::BTreeMap<String, u32> = Default::default();
    reads::fetch(&input, arg(args, "--reference").as_deref(), &format!("{chrom}:{}-{}", p + 1, p + 1), &mut |r: &reads::Read| {
        let (mut rp, mut qp) = (r.start, 0usize);
        let mut base = None;
        for &(op, n) in &r.cigar {
            let n = n as u64;
            match op {
                reads::Op::Match => { if p >= rp && p < rp + n { let i = qp + (p - rp) as usize; base = Some(((r.seq[i] as char).to_string(), r.qual[i])); break; } rp += n; qp += n as usize; }
                reads::Op::Ins | reads::Op::SoftClip => qp += n as usize,
                reads::Op::Del | reads::Op::Skip => { if p >= rp && p < rp + n { base = Some(("-".into(), 0)); break; } rp += n; }
                _ => {}
            }
        }
        if let Some((b, q)) = base {
            let strand = if r.flags & 0x10 != 0 { "rev" } else { "fwd" };
            println!("{}	{}	{strand}	mapq={}	{b}	bq={q}	start={}", String::from_utf8_lossy(&r.name), r.flags, r.mapq, r.start + 1);
            if r.flags & 0x704 == 0 && r.mapq >= 20 { *counts.entry(format!("{b} {strand}")).or_default() += 1; }
        }
    })?;
    eprintln!("{counts:?}");
    Ok(())
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let res = match args.first().map(|s| s.as_str()) {
        Some("coverage") => cmd_coverage(&args[1..]),
        Some("depth") => cmd_depth(&args[1..]),
        Some("exons") => cmd_exons(&args[1..]),
        Some("readqc") => cmd_readqc(&args[1..]),
        Some("call") => cmd_call(&args[1..]),
        Some("pileup") => cmd_pileup(&args[1..]),
        Some("annotate") => cmd_annotate(&args[1..]),
        Some("exon-cnv") => cmd_exon_cnv(&args[1..]),
        Some("sv") => cmd_sv(&args[1..]),
        Some("--version") | Some("-V") => {
            println!("brca12-wgs {}", env!("CARGO_PKG_VERSION"));
            Ok(())
        }
        Some(a) if !a.starts_with('-') && args.len() >= 2 => cmd_run(&args),
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
