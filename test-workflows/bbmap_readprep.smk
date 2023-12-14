#!/usr/bin/env python3

from pathlib import Path

test_data = Path("test-data", "star")
adaptors = Path("test-data", "bbmap_readprep", "bbmap_39.01_adaptors.fa")
test_outdir = Path("test-output", "bbmap_readprep")

# bbmap_snakefile = "../modules/bbmap_readprep/Snakefile"
bbmap_snakefile = github(
    "tomharrop/smk-modules",
    path="modules/bbmap_readprep/Snakefile",
    tag="0.0.25",
)

rule target:
    input:
        expand(
            Path(test_outdir, "{sample}.{readset}.fastq.gz"),
            sample=["GSM461177_subsampled","GSM461180_subsampled"],
            readset=["r1", "r2"] # can also ask for unpaired

        )

module bbmap_readprep:
    snakefile:
        bbmap_snakefile
    config:
        {
            "outdir": test_outdir,
            "reads_directory": test_data,
            "adaptors": adaptors,
        }


use rule * from bbmap_readprep as bbmap_readprep_*
