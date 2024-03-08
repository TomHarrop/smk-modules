#!/usr/bin/env python3

from pathlib import Path
import tempfile

query_genome = Path("test-data", "repeatmasker", "genome.fa.gz")
outdir = Path(
    "test-output",
    "repeatmasker",
)


# configure the run like this, or in a yaml file
if "repeatmasker" not in config.keys():
    config["repeatmasker"] = {}

repeatmasker_config = config["repeatmasker"]
repeatmasker_config["outdir"] = outdir
repeatmasker_config["query_genome"] = query_genome
config["repeatmasker"] = repeatmasker_config


rm_snakefile = github(
    "tomharrop/smk-modules",
    path="modules/repeatmasker/Snakefile",
    tag="0.0.46",
)
# rm_snakefile = "../modules/repeatmasker/Snakefile"


module repeatmasker:
    snakefile:
        rm_snakefile
    config:
        config["repeatmasker"]


use rule * from repeatmasker
