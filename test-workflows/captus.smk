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
            "sample_list": all_samples,
            "read_directory": Path("inputs"),
            "target_file": "targetfile.fasta",
        }
    prefix:
        Path("test-output", "captus")


use rule * from captus as captus_*


rule set_up_captus_inputs:
    input:
        read_directory=read_directory,
        target_file=target_file,
    output:
        reads=expand(
            Path(
                "test-output",
                "captus",
                "inputs",
                "{sample}.r{r}.fastq.gz",
            ),
            sample=all_samples,
            r=["1", "2"],
        ),
        target_file=Path(
            "test-output",
            "captus",
            "targetfile.fasta",
        ),
    params:
        read_directory=lambda wildcards, output: Path(output.reads[0]).parent,
    shell:
        "rm -r {params.read_directory} ; "
        "ln -s "
        "$(readlink -f {input.read_directory}) "
        "$(readlink -f {params.read_directory} ) ; "
        "ln -s "
        "$(readlink -f {input.target_file} ) "
        "$(readlink -f {output.target_file} ) ; "


rule target:
    input:
        rules.captus_target.input,
        # Path(
        #     "test-output",
        #     "captus",
        #     "04_alignments",
        #     "04_alignments",
        #     "02_untrimmed",
        #     "06_informed",
        #     "03_coding_MIT",
        #     "02_NT",
        # ),
    default_target: True
