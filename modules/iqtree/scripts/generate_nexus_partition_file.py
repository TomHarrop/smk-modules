#!/usr/bin/env python3

from Bio import SeqIO
from pathlib import Path
import os
import sys


def count_fasta_records(fasta_file):
    return len(SeqIO.index(fasta_file, "fasta"))


def create_nexus_partition(directory, outfile, reportfile):
    partition_content = []
    file_list = sorted(os.listdir(directory))

    with open(reportfile, "a") as f:
        for filename in file_list:
            if not filename.startswith("."):
                alignment_name = os.path.splitext(filename)[0]
                filepath = os.path.join(directory, filename)
                if count_fasta_records(filepath) < 3:
                    f.write(f"{filepath}\n")
                else:
                    partition_content.append(
                        f"charset {alignment_name} = {filepath}: *;"
                    )

    nexus_content = "#NEXUS\n\nbegin sets;\n"

    for partition_line in partition_content:
        nexus_content += f"    {partition_line}\n"

    nexus_content += "end;"

    with open(outfile, "wt") as f:
        print(nexus_content, file=f)


if __name__ == "__main__":
    directory_path = Path(snakemake.params["alignment_directory"])
    outfile = snakemake.output["nexus"]
    reportfile = snakemake.output["report"]
    create_nexus_partition(directory_path, outfile, reportfile)
    with open(snakemake.output["pathfile"], "wt") as f:
        f.write(directory_path.as_posix())
