#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import tempfile

sample_data = Path("test-data", "hybpiper", "samples.csv")
target_file = Path("test-data", "hybpiper", "combined_targetfiles.fixed.fasta")
read_directory = Path("test-data", "captus", "reads")
output_directory = Path(
    "test-output",
    "captus",
)

# just run on whatever read files we have
all_samples = sorted(
    set(x.stem.split(".")[0] for x in read_directory.glob("*r1.fastq.gz"))
)

# captus_snakefile = "../modules/captus/Snakefile"
captus_snakefile = github(
    "tomharrop/smk-modules",
    path="modules/captus/Snakefile",
    tag="0.0.53",
)


module captus:
    snakefile:
        captus_snakefile
    config:
        {
            "sample_list": all_samples,
            "outdir": output_directory,
            "read_directory": read_directory,
            "run_tmpdir": Path(output_directory, "tmp"),
            "target_file": target_file,
        }


use rule * from captus as captus_*
