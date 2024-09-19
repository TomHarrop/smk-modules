#!/usr/bin/env python3

from pathlib import Path


test_alignments = Path("test-data", "iqtree")
output_directory = Path(
    "test-output",
    "trimal",
)

# For now, don't run trimal on the already-trimalled alignment. Trimal removes
# all sequences from some alignments on the second pass.
datasets = [
    x.name
    for x in test_alignments.glob("*")
    if (x.is_dir() and not "trimal" in x.name)
]

# trimal_snakefile = "../modules/trimal/Snakefile"
trimal_snakefile = github(
    "tomharrop/smk-modules",
    path="modules/trimal/Snakefile",
    tag="0.6.3",
)


module trimal:
    snakefile:
        trimal_snakefile
    config:
        {
            "alignment_directory": Path(test_alignments, "{dataset}"),
            "outdir": Path(output_directory, "{dataset}"),
        }


use rule * from trimal as trimal_*


rule target:
    default_target: True
    input:
        expand([x for x in rules.trimal_target.input], dataset=datasets),
