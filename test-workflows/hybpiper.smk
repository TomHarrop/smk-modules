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

# configure the run like this, or in a yaml file
if "hybpiper" not in config.keys():
    config["hybpiper"] = {}

hybpiper_config = config["hybpiper"]
hybpiper_config["run_tmpdir"] = tempfile.mkdtemp()
hybpiper_config["sample_list"] = all_samples
hybpiper_config["read_directory"] = read_directory
hybpiper_config["target_file"] = target_file
hybpiper_config["outdir"] = output_directory
config["hybpiper"] = hybpiper_config


module hybpiper:
    snakefile:
        github(
            "tomharrop/smk-modules",
            path="modules/hybpiper/Snakefile",
            tag="0.0.16",
        )
    config:
        config["hybpiper"]


use rule * from hybpiper
