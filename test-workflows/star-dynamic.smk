#!/usr/bin/env python3

from pathlib import Path
import csv
import pandas as pd
import re

test_data = Path("test-data", "star", "reads")
test_outdir = Path("test-output", "star_dynamic")

reference = Path(test_data, "NT_033779.5.fna.gz")
annotation = Path(test_data, "NT_033779.5.gff.gz")

tempdir = Path(test_outdir, "tmp")
csv_file = Path(tempdir, "samples.csv")
csv_file.parent.mkdir(parents=True, exist_ok=True)

# generate a CSV of sample name, r1, r2
pattern = re.compile(
    r"^(?P<sample>[a-zA-Z0-9_]+?)\.(?P<read_type>r[12])\.(?P<extension>fq|fastq)(?P<gz_extension>\.gz)?$"
)

files = test_data.glob("*")

# Initialize the dictionary to store the paths
file_dict = {}

for file_path in files:
    file_name = file_path.name
    match = pattern.match(file_name)
    if match:
        groups = match.groupdict()
        sample_name = groups["sample"]
        read_type = groups["read_type"]
        extension = groups["extension"]
        gz_extension = groups["gz_extension"] or ""
        # Create sub-dictionary if sample_name is not already in the main dictionary
        if sample_name not in file_dict:
            file_dict[sample_name] = {}
        # Add the Path object to the sub-dictionary
        read_key = f"r{read_type[-1]}"
        file_dict[sample_name][read_key] = file_path


with open(csv_file, "w", newline="") as f:
    fieldnames = ["sample_name", "r1", "r2"]
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    # Write the header
    writer.writeheader()
    # Write each row
    for sample_name, read_paths in file_dict.items():
        row = {
            "sample_name": sample_name,
            "r1": read_paths.get("r1", ""),
            "r2": read_paths.get("r2", ""),
        }
        writer.writerow(row)


module star_dynamic:
    snakefile:
        "../modules/star/wait_for_input.smk"
    config:
        {"sample_data": csv_file}
