#! /usr/bin/env python3

from pathlib import Path

read_directory = Path("test-data", "tcsnps", "reads")
output_directory = Path(
    "test-output",
    "captus",
)

# Eucalyptus grandis Chr1
ref_genome = Path("test-data", "tcsnps", "reference", "NC_052612.1.fa")

# just run on whatever read files we have
all_samples = sorted(
    set(x.stem.split(".")[0] for x in read_directory.glob("*r1.fastq.gz"))
)

tcsnps_snakefile = "../modules/tcsnps/Snakefile"
# tcsnps_snakefile = github(
#     "tomharrop/smk-modules",
#     path="modules/tcsnps/Snakefile",
#     tag="0.5.1",
# )


module tcsnps:
    snakefile:
        tcsnps_snakefile
    config:
        {
            "sample_list": all_samples,
            "read_directory": read_directory,
            "ref_genome": ref_genome,
            "outdir": Path("test-output", "tcsnps"),
        }


use rule * from tcsnps as tcsnps_*
