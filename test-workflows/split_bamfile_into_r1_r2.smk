#!/usr/bin/env python3

from pathlib import Path
import tempfile

# my_snakefile = "../modules/split_bamfile_into_r1_r2/Snakefile"
my_snakefile = github(
    "tomharrop/smk-modules",
    path="modules/split_bamfile_into_r1_r2/Snakefile",
    tag="0.0.38",
)

my_outdir = Path("test-output", "split_bamfile_into_r1_r2")


module split_bamfile:
    snakefile:
        my_snakefile
    config:
        {
            "outdir": my_outdir,
            "run_tmpdir": Path(my_outdir, "tmp"),
            "bamfile": Path("test-data", "braker3", "RNAseq.bam"),
        }


use rule * from split_bamfile as split_bamfile_*
