#!/usr/bin/env python3

from pathlib import Path

paralog_sequences = Path("test-data", "paragone", "paralog_input")
# the external outgroups can be a target file with a known outgroup ID
external_outgroups = Path("test-data", "paragone", "external_outgroups.fasta")
internal_outgroup = "80974"  # taxon id?


# paragone_snakefile = github(
#     "tomharrop/smk-modules",
#     path="modules/paragone/Snakefile",
#     tag="0.0.22",
# )
paragone_snakefile = "../modules/paragone/Snakefile"


rule target:
    input:
        expand(
            "test-output/paragone/{run}/intermediate_files.tar.gz",
            run=["internal", "external"],
        ),


module paragone_external:
    snakefile:
        paragone_snakefile
    config:
        {
            "external_outgroups": external_outgroups,
            "paralog_sequences": paralog_sequences,
            "outdir": Path("test-output", "paragone", "external"),
            "pool": 3,
        }


use rule * from paragone_external as pgext_*


module paragone_internal:
    snakefile:
        paragone_snakefile
    config:
        {
            "internal_outgroup": internal_outgroup,
            "paralog_sequences": paralog_sequences,
            "outdir": Path("test-output", "paragone", "internal"),
            "pool": 3,
        }


use rule * from paragone_internal as pgint_*
