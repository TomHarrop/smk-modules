#!/usr/bin/env python3

from pathlib import Path
import os
import sys


def create_nexus_partition(directory, outfile):
    partition_content = []
    file_list = sorted(os.listdir(directory))

    for filename in file_list:
        if not filename.startswith("."):
            alignment_name = os.path.splitext(filename)[0]
            filepath = os.path.join(directory, filename)
            partition_content.append(f"charset {alignment_name} = {filepath}: *;")

    nexus_content = "#NEXUS\n\nbegin sets;\n"

    for partition_line in partition_content:
        nexus_content += f"    {partition_line}\n"

    nexus_content += "end;"

    with open(outfile, "wt") as f:
        print(nexus_content, file=f)


if __name__ == "__main__":
    directory_path = Path(snakemake.input["alignment_directory"])
    outfile = snakemake.output[0]
    create_nexus_partition(directory_path, outfile)
