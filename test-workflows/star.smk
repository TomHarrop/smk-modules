#!/usr/bin/env python3

from pathlib import Path

test_data = Path("test-data", "star")
test_outdir = Path("test-output", "star")

reference = Path(test_data, "NT_033779.5.fna.gz")
annotation = Path(test_data, "NT_033779.5.gff.gz")

star_snakefile = "../modules/star/Snakefile"

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


module star:
    snakefile:
        star_snakefile
    config:
        {
            "reference": reference,
            "annotation": annotation,
            "outdir": test_outdir,
            "run_tmpdir": Path(test_outdir, "tmp"),
            "samples": sample_dict,
        }


use rule * from star as star_*
