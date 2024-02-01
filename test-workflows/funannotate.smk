#!/usr/bin/env python3

from pathlib import Path
import tempfile

########
# NOTE #
########

# The funannotate container doesn't work with the usual --containall --cleanenv
# --writable-tmpfs arguments.

###########
# GLOBALS #
###########

bbmap = "docker://quay.io/biocontainers/bbmap:39.01--h92535d8_1"

# from http://bioinf.uni-greifswald.de/augustus/datasets/RNAseq.bam
rnaseq = Path("test-data", "braker3", "RNAseq.bam")
genome = Path("test-data", "braker3", "genome.fa.gz")
db_path = Path("test-data", "funannotate", "db")
dmnd_db = Path("test-data", "funannotate", "eggnog", "eggnog_proteins.dmnd")
eggnog_db = Path("test-data", "funannotate", "eggnog", "eggnog.db")

# use output from braker as evidence (just for testing)
protein_evidence = Path("test-data", "funannotate", "braker.aa")
transcript_evidence = Path("test-data", "funannotate", "braker.codingseq")

# this is the path to the included adaptors file in bbmap
bbmap_adaptors = Path(
    "/usr", "local", "opt", "bbmap-39.01-1", "resources", "adapters.fa"
)

outdir = Path(
    "test-output",
    "funannotate",
)
logdir = Path(outdir, "logs")
# avoid rerunning steps
run_tmpdir = Path(outdir, "tmp")

################################################################################
# Example configuration for the funannotate module
################################################################################
# Set interproscan_container to False to disable interproscan. If you don't
# provide gm_key, the module will try to use one from the container. This will
# fail if that key has expired.
#
# Optional inputs:
#   - protein_evidence
#   - transctipt_evidence
#   - rnaseq
#   - min_training_models (default 200, set lower for test data)
#
# Configuring busco:
#   - run funannotate species to find a list of species for busco_seed_species
#   - run funannotate database --show-buscos to find a list of species you can
#     use for busco_db
#   - you can add busco lineages to the db folder to make them available

fa_config = {
    "db_path": db_path,
    "dmnd_db": dmnd_db,
    "eggnog_db": eggnog_db,
    "gm_key": Path("test-data", "funannotate", "gm_key_64"),
    # "interproscan_container": False,
    "interproscan_container": "interproscan_5.65-97.0_cv3.sif",
    "min_training_models": 20,
    "outdir": outdir,
    "protein_evidence": protein_evidence,
    "query_genome": genome,
    "rnaseq_r1": Path(outdir, "reads", "reads.trimmed.r1.fq.gz"),
    "rnaseq_r2": Path(outdir, "reads", "reads.trimmed.r2.fq.gz"),
    "run_tmpdir": run_tmpdir,
    "species_name": "testspecies",
    "transcript_evidence": transcript_evidence,
    "busco_seed_species": "arabidopsis",
    "busco_db": "embryophyta",
}

################################################################################

#########
# RULES #
#########

# fa_snakefile = "../modules/funannotate/Snakefile"
fa_snakefile = github(
    "tomharrop/smk-modules",
    path="modules/funannotate/Snakefile",
    tag="0.0.39",
)


module funannotate:
    snakefile:
        fa_snakefile
    config:
        fa_config


use rule * from funannotate as funannotate_*


module split_bamfile:
    snakefile:
        github(
            "tomharrop/smk-modules",
            path="modules/split_bamfile_into_r1_r2/Snakefile",
            tag="0.0.38",
        )
    config:
        {
            "outdir": Path(outdir, "reads"),
            "run_tmpdir": Path(outdir, "tmp"),
            "bamfile": rnaseq,
        }

use rule * from split_bamfile as split_bamfile_*