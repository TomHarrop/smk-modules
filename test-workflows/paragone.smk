#!/usr/bin/env python3

from pathlib import Path

paralog_sequences = Path("test-data", "paragone", "paralog_input")
# the external outgroups can be a target file with a known outgroup ID
external_outgroups = Path("test-data", "paragone", "external_outgroups.fasta")
internal_outgroup = "80974"  # taxon id?


module paragone:
    snakefile:
        # github(
        #     "tomharrop/smk-modules",
        #     path="modules/hybpiper/Snakefile",
        #     tag="0.0.17",
        # )
        "../modules/paragone/Snakefile"
    config:
        {
            "external_outgroups": external_outgroups,
            "internal_outgroup": internal_outgroup,
            "paralog_sequences": paralog_sequences,
            "outdir": Path("test-output", "paragone"),
            "pool": 3
        }


use rule * from paragone
