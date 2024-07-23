#!/usr/bin/env python3

from Bio import SeqIO
from pathlib import Path
import os
import tarfile


def count_fasta_records(fasta_file):
    return len(SeqIO.index(fasta_file, "fasta"))


def create_nexus_partition(directory, nexus_file, output_tarfile, excluded_samples):
    partition_content = []
    file_list = sorted(os.listdir(directory))

    with tarfile.open(output_tarfile, "w") as tar:
        with open(excluded_samples, "a") as f_excluded:
            for filename in file_list:
                if not filename.startswith("."):
                    alignment_name = os.path.splitext(filename)[0]
                    filepath = os.path.join(directory, filename)
                    if count_fasta_records(filepath) < 3:
                        f_excluded.write(f"{filename}\n")
                    else:
                        partition_content.append(
                            f"charset {alignment_name} = input_alignments/{filename}: *;"
                        )
                        tar.add(filepath, arcname=filename)

    nexus_content = "#NEXUS\n\nbegin sets;\n"

    for partition_line in partition_content:
        nexus_content += f"    {partition_line}\n"

    nexus_content += "end;"

    with open(nexus_file, "wt") as f:
        print(nexus_content, file=f)


if __name__ == "__main__":
    directory_path = Path(snakemake.input["alignment_directory"])
    output_tarfile = snakemake.output["kept_alignments"]
    nexus_file = snakemake.output["nexus"]
    excluded_samples = snakemake.output["excluded_samples"]
    create_nexus_partition(directory_path, nexus_file, output_tarfile, excluded_samples)
