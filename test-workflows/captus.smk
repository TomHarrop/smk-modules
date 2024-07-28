#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import tempfile


# demonstrate how to use a param to find the input directory for e.g. iqtree
def get_alignment_dir(wildcards):
    captus_align_outdir = rules.captus_align.output["outdir"]
    captus_align_output = [x for x in Path(captus_align_outdir).glob("**")]
    trimpath = Path(
        "03_trimmed",
        "06_informed",
        "03_coding_MIT",
        f"{wildcards.alignment_type}",
    )
    untrimpath = Path(
        "02_untrimmed",
        "06_informed",
        "03_coding_MIT",
        f"{wildcards.alignment_type}",
    )
    for path in captus_align_output:
        if path.match(Path("**", trimpath).as_posix()) and path.is_dir():
            return path
    for path in captus_align_output:
        if path.match(Path("**", untrimpath).as_posix()) and path.is_dir():
            return path
    return ""


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

# captus_snakefile = "../modules/captus/Snakefile"
captus_snakefile = github(
    "tomharrop/smk-modules",
    path="modules/captus/Snakefile",
    tag="0.3.02",
)


captus_alignments = ["01_AA", "02_NT", "03_genes"]


rule target:
    input:
        expand(
            Path(
                "test-output", "captus", "99_post-captus", "{alignment_type}"
            ),
            alignment_type=captus_alignments,
        ),


rule post_captus:
    input:
        Path("test-output", "captus", "04_alignments"),
    output:
        directory(
            Path(
                "test-output", "captus", "99_post-captus", "{alignment_type}"
            )
        ),
    params:
        alignment_dir=get_alignment_dir,
    shell:
        "cp -r {params.alignment_dir} {output}"


# NOTE. this does NOT work well with `prefix:` because of the way captus parses
# directory arguments. better to use the outdir where possible.
# Another NOTE. It's very hard to make this work without a checkpoint. If you
# pass the samples as a list, you have to define the module once for each list.
# If you pass the samples in a namelist, you can't start the workflow until
# you've read the list.
module captus:
    snakefile:
        captus_snakefile
    config:
        {
            "namelist": Path("test-data", "captus", "namelist.txt"),
            "read_directory": Path("test-output", "captus", "inputs"),
            "target_file": target_file,
            "outdir": Path("test-output", "captus"),
        }


use rule * from captus as captus_*


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
