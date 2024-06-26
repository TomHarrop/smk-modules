#!/usr/bin/env python3

from pathlib import Path

# quick to download, but only one sample
accession = "PRJNA268659"

# a ton of small samples
# accession = "PRJNA307264"

# an actual dataset of interest
# accession = "PRJNA549921"

outdir = Path("test-output", "bpdownload", accession)


def get_bpdownload_output(wildcards):
    """
    The bpdownload rule called "target" is a checkpoint. Running this after the
    module completes allows you to collect the list of downloaded samples.
    """
    checkpoints.bpdownload_target.get(**wildcards)
    all_samples = glob_wildcards(Path(outdir, "{sample}.r1.fastq.gz")).sample
    return expand(Path(outdir, "{sample}.done"), sample=all_samples)


# the following two rules allow you to use the output for subsequent steps
# without knowing what they are in advance.
rule target:
    input:
        get_bpdownload_output,


rule use_bpdownload_output:
    input:
        r1=Path(outdir, "{sample}.r1.fastq.gz"),
        r2=Path(outdir, "{sample}.r2.fastq.gz"),
    output:
        touch(Path(outdir, "{sample}.done")),


# Call the module.
bpdownload_snakefile = github(
    "tomharrop/smk-modules",
    path="modules/bpdownload/Snakefile",
    tag="0.1.04",
)
# bpdownload_snakefile = "../modules/bpdownload/Snakefile"


module bpdownload:
    snakefile:
        bpdownload_snakefile
    config:
        {
            "accession": accession,
            "outdir": outdir,
            "run_tmpdir": Path(outdir, "tmp"),
        }


use rule * from bpdownload as bpdownload_*
