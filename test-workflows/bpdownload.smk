#!/usr/bin/env python3

from pathlib import Path

accession = "PRJNA975329"

outdir = Path(
    "test-output",
    "bpdownload",
)


module bpdownload:
    snakefile:
        # github(
        #     "tomharrop/smk-modules",
        #     path="modules/braker3/Snakefile",
        #     tag="0.0.14",
        # )
        "../modules/bpdownload/Snakefile"
    config:
        {"accession": accession, "outdir": outdir}


use rule * from bpdownload
