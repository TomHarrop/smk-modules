#!/usr/bin/env python3

from pathlib import Path

paralog_sequences = Path("test-data", "paragone", "paralog_input")
# the external outgroups can be a target file with a known outgroup ID
external_outgroups = Path("test-data", "paragone", "external_outgroups.fasta")
internal_outgroup = "80974"  # taxon id?


paragone_snakefile = github(
    "tomharrop/smk-modules",
    path="modules/paragone/Snakefile",
    tag="0.2.10",
)
# paragone_snakefile = "../modules/paragone/Snakefile"


module paragone_external:
    snakefile:
        paragone_snakefile
    config:
        {
            "external_outgroups": "external_outgroups.fasta",
            "paralog_sequences": "paralog_input",
            "pool": 3,
        }
    prefix:
        Path("test-output", "paragone", "external")


use rule * from paragone_external as pgext_*


rule set_up_paragone_inputs:
    input:
        external_outgroups=external_outgroups,
        paralog_sequences=paralog_sequences,
    output:
        paralog_sequences=temp(
            directory(
                Path("test-output", "paragone", "external", "paralog_input")
            )
        ),
        external_outgroups=temp(
            Path(
                "test-output",
                "paragone",
                "external",
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


rule target:
    input:
        rules.pgint_target.input,
        rules.pgext_target.input,
    default_target: True
