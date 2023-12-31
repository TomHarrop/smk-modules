#!/usr/bin/env python3

from pathlib import Path

proteins = Path("test-data", "braker3", "proteins.fa.gz")
query_genome = Path("test-data", "braker3", "genome.fa.gz")
# from http://bioinf.uni-greifswald.de/augustus/datasets/RNAseq.bam
rnaseq = Path("test-data", "braker3", "RNAseq.bam")

species_name = "test_species"

outdir = Path(
    "test-output",
    "braker3",
)

# configure the run like this, or in a yaml file
if "braker3" not in config.keys():
    config["braker3"] = {}

braker3_config = config["braker3"]

braker3_config["outdir"] = outdir
braker3_config["proteins"] = proteins
braker3_config["query_genome"] = query_genome
braker3_config["rnaseq"] = rnaseq
braker3_config["species_name"] = species_name

config["braker3"] = braker3_config


rule target:
    input:
        Path(outdir, "braker.gff3"),


module braker3:
    snakefile:
        github(
            "tomharrop/smk-modules",
            path="modules/braker3/Snakefile",
            tag="0.0.31",
        )
    config:
        config["braker3"]


use rule * from braker3
