#!/usr/bin/env python3

from pathlib import Path

paralog_sequences = Path("test-data", "paragone", "test-orient2")
# the external outgroups can be a target file with a known outgroup ID
external_outgroups = Path("test-data", "paragone", "external_outgroups.fasta")


# paragone_snakefile = github(
#     "tomharrop/smk-modules",
#     path="modules/paragone/Snakefile",
#     tag="0.2.18",
# )
paragone_snakefile = "../modules/paragone/Snakefile"


module paragone:
    snakefile:
        paragone_snakefile
    config:
        {
            "external_outgroups": "external_outgroups.fasta",
            "paralog_sequences": "paralog_input",
            "pool": 3,
        }
    prefix:
        Path("test-output", "paragone", "orient")


use rule * from paragone as paragone_*


rule set_up_paragone_inputs:
    input:
        external_outgroups=external_outgroups,
        paralog_sequences=paralog_sequences,
    output:
        paralog_sequences=temp(
            directory(
                Path("test-output", "paragone", "orient", "paralog_input")
            )
        ),
        external_outgroups=temp(
            Path(
                "test-output",
                "paragone",
                "orient",
                "external_outgroups.fasta",
            )
        ),
    shell:
        "ln -s "
        "$(readlink -f {input.paralog_sequences}) "
        "$(readlink -f {output.paralog_sequences} ) ; "
        "ln -s "
        "$(readlink -f {input.external_outgroups} ) "
        "$(readlink -f {output.external_outgroups} ) ; "


rule target:
    input:
        rules.paragone_target.input,
    default_target: True
