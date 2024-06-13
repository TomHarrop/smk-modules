from pathlib import Path
import pandas as pd
import tempfile

test_alignments = Path("test-data", "iqtree")
output_directory = Path(
    "test-output",
    "iqtree",
)

# iqtree_snakefile = "../modules/iqtree/Snakefile"
iqtree_snakefile = github(
    "tomharrop/smk-modules",
    path="modules/iqtree/Snakefile",
    tag="0.0.56",
)


module iqtree_captus:
    snakefile:
        iqtree_snakefile
    config:
        {
            "alignment_directory": Path(test_alignments, "captus"),
            "outdir": Path(output_directory, "captus"),
            "run_tmpdir": Path(output_directory, "tmp", "captus"),
        }


use rule * from iqtree_captus as iqtree_captus_*


module iqtree_paragone:
    snakefile:
        iqtree_snakefile
    config:
        {
            "alignment_directory": Path(test_alignments, "paragone"),
            "outdir": Path(output_directory, "paragone"),
            "run_tmpdir": Path(output_directory, "tmp", "paragone"),
        }


use rule * from iqtree_paragone as iqtree_paragone_*


rule target:
    input:
        rules.iqtree_paragone_target.input,
        rules.iqtree_captus_target.input,
    default_target: True
