from pathlib import Path
import pandas as pd
import tempfile

test_alignments = Path("test-data", "iqtree")
output_directory = Path(
    "test-output",
    "iqtree",
)

datasets = [x.name for x in test_alignments.glob("*") if x.is_dir()]

iqtree_snakefile = "../modules/iqtree/Snakefile"
# iqtree_snakefile = github(
#     "tomharrop/smk-modules",
#     path="modules/iqtree/Snakefile",
#     tag="0.1.01",
# )


rule target:
    input:
        expand(
            Path(output_directory, "{dataset}", "tree.treefile"),
            dataset=datasets,
        ),


module iqtree:
    snakefile:
        iqtree_snakefile
    config:
        {
            "alignment_directory": Path(test_alignments, "{dataset}"),
            "outdir": Path(output_directory, "{dataset}"),
            "run_tmpdir": Path(output_directory, "tmp", "{dataset}"),
        }


use rule * from iqtree as iqtree_*
