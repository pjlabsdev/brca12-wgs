"""Single-copy control windows (chr13/chr17, away from BRCA1/BRCA2) used as depth normaliser."""
CONTROLS = {
    "c17_TP53": {"chrom": "chr17", "start": 7668000, "end": 7672000},
    "c17_q12": {"chrom": "chr17", "start": 39700000, "end": 39704000},
    "c17_q21b": {"chrom": "chr17", "start": 50000000, "end": 50004000},
    "c17_q23": {"chrom": "chr17", "start": 61000000, "end": 61004000},
    "c13_FLT3": {"chrom": "chr13", "start": 28000000, "end": 28004000},
    "c13_RB1": {"chrom": "chr13", "start": 48300000, "end": 48304000},
    "c13_q31": {"chrom": "chr13", "start": 90000000, "end": 90004000},
    "c13_q33": {"chrom": "chr13", "start": 100000000, "end": 100004000},
}
