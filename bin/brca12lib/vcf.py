"""Minimaler VCF-Leser ohne Fremdabhängigkeiten (bgzip/gzip/plain)."""
import gzip
import re

_CSQ_FORMAT = re.compile(r'Format: ([^"]+)"')


def _open(path):
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_info(text):
    info = {}
    if text in ("", "."):
        return info
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = True
    return info


class VCF:
    def __init__(self, path):
        self.path = path
        self.csq_fields = []
        self.samples = []
        self.records = []
        with _open(path) as handle:
            for line in handle:
                if line.startswith("##"):
                    if line.startswith("##INFO=<ID=CSQ"):
                        match = _CSQ_FORMAT.search(line)
                        if match:
                            self.csq_fields = match.group(1).strip().split("|")
                    continue
                if line.startswith("#CHROM"):
                    self.samples = line.rstrip("\n").split("\t")[9:]
                    continue
                if line.strip():
                    self.records.append(self._record(line.rstrip("\n").split("\t")))

    def _record(self, fields):
        record = {
            "chrom": fields[0], "pos": int(fields[1]), "id": fields[2], "ref": fields[3],
            "alts": fields[4].split(","), "qual": fields[5], "filter": fields[6],
            "info": parse_info(fields[7]), "sample": {},
        }
        if len(fields) > 9:
            keys = fields[8].split(":")
            values = fields[9].split(":")
            record["sample"] = dict(zip(keys, values))
        if self.csq_fields and "CSQ" in record["info"]:
            record["csq"] = [dict(zip(self.csq_fields, entry.split("|")))
                             for entry in str(record["info"]["CSQ"]).split(",")]
        else:
            record["csq"] = []
        return record


def to_int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def genotype_label(gt):
    if not gt or gt in (".", "./.", ".|."):
        return "nicht bestimmt"
    alleles = re.split(r"[/|]", gt)
    if all(a == "0" for a in alleles):
        return "homozygot Referenz"
    if len(set(alleles)) == 1:
        return "homozygot"
    if "0" in alleles:
        return "heterozygot"
    return "compound/multiallelisch"
