#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import tempfile

sample_data = Path("test-data", "hybpiper", "samples.csv")
target_file = Path("test-data", "hybpiper", "combined_targetfiles.fixed.fasta")
read_directory = Path("test-data", "hybpiper", "reads")
output_directory = Path(
    "test-output",
    "hybpiper",
)

samples = pd.read_csv(sample_data, index_col="name")
all_samples = sorted(set(samples.index))


hybpiper_snakefile = github(
    "tomharrop/smk-modules",
    path="modules/hybpiper/Snakefile",
    tag="0.2.04",
)
# hybpiper_snakefile = "../modules/hybpiper/Snakefile"


module hybpiper:
    snakefile:
        hybpiper_snakefile
    config:
        {
            "namelist": Path("test-data", "hybpiper", "namelist.txt"),
            "outdir": output_directory,
            "read_directory": read_directory,
            "target_file": target_file,
        }


use rule * from hybpiper as hybpiper_*


checkpoint generate_namelist:
    output:
        Path("test-data", "hybpiper", "namelist.txt"),
    run:
        with open(output[0], "wt") as f:
            f.write("\n".join(all_samples))


rule target:
    default_target: True
    input:
        rules.hybpiper_target.input,
