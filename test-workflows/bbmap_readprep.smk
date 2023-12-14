#!/usr/bin/env python3

from pathlib import Path

test_data = Path("test-data", "star")
adaptors = Path("test-data", "bbmap_readprep", "bbmap_39.01_adaptors.fa")
test_outdir = Path("test-output", "bbmap_readprep")

# bbmap_snakefile = "../modules/bbmap_readprep/Snakefile"
bbmap_snakefile = github(
    "tomharrop/smk-modules",
    path="modules/bbmap_readprep/Snakefile",
    tag="0.0.24",
)

# outer keys are sample names, inner keys must be r1 and r2
sample_dict = {
    "GSM461177": {
        "r1": Path(test_data, "GSM461177_subsampled.r1.fq.gz"),
        "r2": Path(test_data, "GSM461177_subsampled.r2.fq.gz"),
    },
    "GSM461180": {
        "r1": Path(test_data, "GSM461180_subsampled.r1.fq.gz"),
        "r2": Path(test_data, "GSM461180_subsampled.r2.fq.gz"),
    },
}


module bbmap_readprep:
    snakefile:
        bbmap_snakefile
    config:
        {"outdir": test_outdir, "samples": sample_dict, "adaptors": adaptors}


use rule * from bbmap_readprep as bbmap_readprep_*
