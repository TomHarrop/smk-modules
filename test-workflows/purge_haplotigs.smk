#!/usr/bin/env python3

from pathlib import Path
import tempfile

###########
# GLOBALS #
###########

purge_snakefile = "../modules/purge_haplotigs/Snakefile"
# purge_snakefile = github(
#     "tomharrop/smk-modules",
#     path="modules/purge_haplotigs/Snakefile",
#     tag="0.0.40",
# )


# from https://bitbucket.org/mroachawri/purge_haplotigs/src/master/test/
bamfile = Path("test-data", "purge_haplotigs", "aligned.bam")
contigs = Path("test-data", "purge_haplotigs", "contigs.fa")

outdir = Path(
    "test-output",
    "purge_haplotigs",
)
# avoid rerunning steps
run_tmpdir = Path(outdir, "tmp")
logdir = Path(outdir, "logs")

#########
# RULES #
#########


# If you run the module without the cutoffs, it just generates the histogram.
# Add the cutoffs to the config to run the following rules.
module purge_haplotigs:
    snakefile:
        purge_snakefile
    config:
        {
            "bamfile": bamfile,
            "contigs": contigs,
            "outdir": outdir,
            "run_tmpdir": run_tmpdir,
            "cutoffs": {
                "low": 10,
                "mid": 30,
                "high": 45,
            },
        }


use rule * from purge_haplotigs as purge_haplotigs_*
