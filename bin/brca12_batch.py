#!/usr/bin/env python3
"""Run brca12-wgs for several samples from a sample sheet and write a cohort summary.

    python3 bin/brca12_batch.py --samplesheet samples.csv --outdir batch_results --run-name RUN01 \\
        --jobs 2 -- -profile conda --reference GRCh38.fa --vep_cache ~/.vep --clinvar_vcf clinvar.vcf.gz

samples.csv (comma or tab separated, header required):
    sample,reads_r1,reads_r2,bam
    S01,S01_R1.fastq.gz,S01_R2.fastq.gz,
    S02,,,S02.cram

Every sample runs as its own Nextflow run in <outdir>/<sample>/ (work directory <outdir>/work/<sample>),
so finished samples are kept and a failed sample can be resumed alone. Arguments after "--" are
passed to every run. Afterwards bin/brca12_cohort.py writes cohort_samples.tsv, cohort_variants.tsv,
cohort_contamination.tsv and cohort.xlsx into <outdir>/cohort/.
"""
import argparse
import csv
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
NAME = re.compile(r"^[A-Za-z0-9._-]+$")


def read_samplesheet(path):
    with open(path, newline="") as handle:
        text = handle.read()
    dialect = "excel-tab" if "\t" in text.splitlines()[0] else "excel"
    rows = list(csv.DictReader(text.splitlines(), dialect=dialect))
    samples, seen = [], set()
    for i, row in enumerate(rows, start=2):
        row = {k.strip().lower(): (v or "").strip() for k, v in row.items() if k}
        name = row.get("sample", "")
        if not NAME.match(name):
            raise SystemExit(f"[brca12] line {i}: invalid or missing sample name {name!r}")
        if name in seen:
            raise SystemExit(f"[brca12] line {i}: sample {name} appears twice")
        seen.add(name)
        bam = row.get("bam") or row.get("cram") or ""
        r1, r2 = row.get("reads_r1", ""), row.get("reads_r2", "")
        if bool(bam) == bool(r1 or r2) or (not bam and not (r1 and r2)):
            raise SystemExit(f"[brca12] line {i}: give either reads_r1 and reads_r2, or bam")
        samples.append({"sample": name, "reads_r1": r1, "reads_r2": r2, "bam": bam})
    if not samples:
        raise SystemExit("[brca12] sample sheet is empty")
    return samples


def command(sample, outdir, nextflow, extra):
    cmd = [nextflow, "run", os.path.join(ROOT, "main.nf"), "--sample", sample["sample"],
           "--outdir", os.path.join(outdir, sample["sample"]), "-work-dir", os.path.join(outdir, "work", sample["sample"])]
    if sample["bam"]:
        cmd += ["--bam", sample["bam"]]
    else:
        cmd += ["--reads_r1", sample["reads_r1"], "--reads_r2", sample["reads_r2"]]
    return cmd + list(extra)


def run(sample, outdir, nextflow, extra, resume):
    cmd = command(sample, outdir, nextflow, extra) + (["-resume"] if resume else [])
    log = os.path.join(outdir, f"{sample['sample']}.nextflow.log")
    with open(log, "w") as handle:
        code = subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT,
                              cwd=os.path.join(outdir, sample["sample"])).returncode
    return sample["sample"], code, log


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    extra = []
    if "--" in argv:
        split = argv.index("--")
        argv, extra = argv[:split], argv[split + 1:]
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--outdir", default="batch_results")
    parser.add_argument("--run-name", default="BRCA12")
    parser.add_argument("--jobs", type=int, default=1, help="samples processed in parallel")
    parser.add_argument("--nextflow", default="nextflow")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--dry-run", action="store_true", help="only print the Nextflow commands")
    args = parser.parse_args(argv)
    samples = read_samplesheet(args.samplesheet)
    outdir = os.path.abspath(args.outdir)
    if args.dry_run:
        for sample in samples:
            print(" ".join(command(sample, outdir, args.nextflow, extra)))
        return 0
    for sample in samples:
        os.makedirs(os.path.join(outdir, sample["sample"]), exist_ok=True)
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        outcomes = list(pool.map(lambda s: run(s, outdir, args.nextflow, extra, args.resume), samples))
    failed = [(name, log) for name, code, log in outcomes if code != 0]
    for name, log in failed:
        print(f"[brca12] {name} failed, see {log}", file=sys.stderr)
    done = [os.path.join(outdir, name) for name, code, _ in outcomes if code == 0]
    if done:
        cohort = [sys.executable, os.path.join(HERE, "brca12_cohort.py"), "--results", *done,
                  "--run-name", args.run_name, "--outdir", os.path.join(outdir, "cohort")]
        subprocess.run(cohort, check=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
