#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import tempfile

sample_data = Path("test-data", "hybpiper", "samples.csv")
target_file = Path("test-data", "hybpiper", "combined_targetfiles.fixed.fasta")
read_directory = Path("test-data", "captus", "reads")
output_directory = Path(
    "test-output",
    "captus",
)

# just run on whatever read files we have
all_samples = sorted(
    set(x.stem.split(".")[0] for x in read_directory.glob("*r1.fastq.gz"))
)

# these are the stats files we want
statsfiles = [
    "captus-assembly_align.alignments",
    "captus-assembly_align.paralogs",
    "captus-assembly_align.samples",
]

captus_snakefile = "../modules/captus/Snakefile"
# captus_snakefile = github(
#     "tomharrop/smk-modules",
#     path="modules/captus/Snakefile",
#     tag="0.1.06",
# )


module captus:
    snakefile:
        captus_snakefile
    config:
        {
            "namelist": Path("namelist.txt"),
            "read_directory": Path("inputs"),
            "target_file": "targetfile.fasta",
        }
    prefix:
        Path("test-output", "captus")


use rule * from captus as captus_*


rule set_up_captus_inputs:
    input:
        expand(
            Path(
                "test-output",
                "captus",
                "inputs",
                "{sample}.r{r}.fastq.gz",
            ),
            sample=all_samples,
            r=["1", "2"],
        ),
        target_file=target_file,
        namelist=Path("test-data", "captus", "namelist.txt"),
    output:
        target_file=temp(
            Path(
                "test-output",
                "captus",
                "targetfile.fasta",
            )
        ),
        namelist=temp(Path("test-output", "captus", "namelist.txt")),
    shell:
        "ln -s "
        "$(readlink -f {input.target_file} ) "
        "$(readlink -f {output.target_file} ) ; "
        "ln -s "
        "$(readlink -f {input.namelist} ) "
        "$(readlink -f {output.namelist} ) ; "


rule set_up_read_files:
    input:
        read_directory=Path(read_directory, "{sample}.r{r}.fastq.gz"),
    output:
        temp(Path("test-output", "captus", "inputs", "{sample}.r{r}.fastq.gz")),
    shell:
        "ln -s $( readlink -f {input} ) $( readlink -f {output} )"


rule generate_namelist:
    output:
        Path("test-data", "captus", "namelist.txt"),
    run:
        with open(output[0], "wt") as f:
            f.write("\n".join(all_samples))


# rule target:
#     input:
#         rules.captus_target.input,
#         # Path(
#         #     "test-output",
#         #     "captus",
#         #     "04_alignments",
#         #     "04_alignments",
#         #     "02_untrimmed",
#         #     "06_informed",
#         #     "03_coding_MIT",
#         #     "02_NT",
#         # ),
#     default_target: True
