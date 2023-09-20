#!/usr/bin/env python3

from pathlib import Path
import tempfile


# from http://bioinf.uni-greifswald.de/augustus/datasets/RNAseq.bam
rnaseq = Path("test-data", "braker3", "RNAseq.bam")
genome = Path("test-data", "braker3", "genome.fa.gz")

outdir = Path(
    "test-output",
    "funannotate",
)


module funannotate:
    snakefile:
        # github(
        #     "tomharrop/smk-modules",
        #     path="modules/braker3/Snakefile",
        #     tag="0.0.7",
        # )
        "../modules/funannotate/Snakefile"
    config:
        {
            "outdir": outdir,
            "rnaseq": rnaseq,
            "run_tmpdir": tempfile.mkdtemp(),
            "query_genome": genome,
            "species_name": "testspecies"
        }


use rule * from funannotate
