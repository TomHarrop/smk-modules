#!/usr/bin/env python3

from pathlib import Path
import tempfile

query_genome = Path("test-data", "repeatmasker", "genome.fa.gz")
outdir = Path(
    "test-output",
    "repeatmasker",
)

rm_output = Path(outdir, "masked.fa.gz")


# configure the run like this, or in a yaml file
if "repeatmasker" not in config.keys():
    config["repeatmasker"] = {}

repeatmasker_config = config["repeatmasker"]
repeatmasker_config["outdir"] = outdir
repeatmasker_config["query_genome"] = query_genome
repeatmasker_config["rm_output"] = rm_output
config["repeatmasker"] = repeatmasker_config


module repeatmasker:
    snakefile:
        github(
            "tomharrop/smk-modules",
            path="modules/repeatmasker/Snakefile",
            tag="0.0.6"
        )
    config:
        config["repeatmasker"]


use rule * from repeatmasker
