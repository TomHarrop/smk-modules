#!/usr/bin/env python3

from pathlib import Path
import tempfile

########
# NOTE #
########

# The funannotate container doesn't work with the usual --containall --cleanenv
# --writable-tmpfs arguments.

###########
# GLOBALS #
###########

bbmap = "docker://quay.io/biocontainers/bbmap:39.01--h92535d8_1"

# from http://bioinf.uni-greifswald.de/augustus/datasets/RNAseq.bam
rnaseq = Path("test-data", "braker3", "RNAseq.bam")
genome = Path("test-data", "braker3", "genome.fa.gz")
db_path = Path("test-data", "funannotate", "db")
dmnd_db = Path("test-data", "funannotate", "eggnog", "eggnog_proteins.dmnd")
eggnog_db = Path("test-data", "funannotate", "eggnog", "eggnog.db")

# this is the path to the included adaptors file in bbmap
bbmap_adaptors = Path(
    "/usr", "local", "opt", "bbmap-39.01-1", "resources", "adapters.fa"
)

outdir = Path(
    "test-output",
    "funannotate",
)
logdir = Path(outdir, "logs")
# avoid rerunning steps
run_tmpdir = Path(outdir, "tmp")

fa_config = {
    "db_path": db_path,
    "dmnd_db": dmnd_db,
    "eggnog_db": eggnog_db,
    "outdir": outdir,
    "query_genome": genome,
    "rnaseq_r1": Path(outdir, "reads", "reads.trimmed.r1.fq.gz"),
    "rnaseq_r2": Path(outdir, "reads", "reads.trimmed.r2.fq.gz"),
    "run_tmpdir": run_tmpdir,
    "species_name": "testspecies",
}

#########
# RULES #
#########


module funannotate:
    snakefile:
        # github(
        #     "tomharrop/smk-modules",
        #     path="modules/braker3/Snakefile",
        #     tag="0.0.7",
        # )
        "../modules/funannotate/Snakefile"
    config:
        fa_config


use rule * from funannotate as funannotate_*


rule target:
    input:
        r1=Path(outdir, "reads", "reads.trimmed.r1.fq.gz"),
        r2=Path(outdir, "reads", "reads.trimmed.r2.fq.gz"),


# Demonstrate how to split the bamfile into R1 and R2 before running
# funannotate.
rule split:
    input:
        Path(run_tmpdir, "reads.trimmed.paired.fq"),
    output:
        r1=Path(outdir, "reads", "reads.trimmed.r1.fq.gz"),
        r2=Path(outdir, "reads", "reads.trimmed.r2.fq.gz"),
    log:
        Path(logdir, "split.log"),
    threads: 1
    resources:
        time=lambda wildcards, attempt: 10 * attempt,
        mem_mb=lambda wildcards, attempt: 2e3 * attempt,
    container:
        bbmap
    shell:
        "reformat.sh "
        "-Xmx{resources.mem_mb}m "
        "-Xms100m "
        "in=stdin.fastq "
        "int=t "
        "addcolon=t "
        "zl=9 "
        "out={output.r1} "
        "out2={output.r2} "
        "< {input} "
        "2> {log}"


rule trim:
    input:
        reads=Path(run_tmpdir, "reads.{type}.fq"),
    output:
        pipe=pipe(Path(run_tmpdir, "reads.trimmed.{type}.fq")),
    params:
        interleaved=lambda wildcards: "t"
        if wildcards.type == "paired"
        else "f",
        adaptors=bbmap_adaptors,
    log:
        log=Path(logdir, "trim.{type}.log"),
        stats=Path(logdir, "trim.{type}.stats.txt"),
    threads: 2
    resources:
        time=lambda wildcards, attempt: 10 * attempt,
        mem_mb=lambda wildcards, attempt: 2e3 * attempt,
    container:
        bbmap
    shell:
        "bbduk.sh "
        "-Xmx{resources.mem_mb}m "
        "-Xms100m "
        "threads={threads} "
        "in={input.reads} "
        "int={params.interleaved} "
        "out=stdout.fastq "
        "outs=/dev/null "
        "ref={params.adaptors} "
        "ktrim=r k=23 mink=11 hdist=1 tpe tbo "
        "forcetrimmod=5 "
        "stats={log.stats} "
        ">> {output.pipe} "
        "2> {log.log} "


rule repair:
    input:
        Path(run_tmpdir, "reads.fq"),
    output:
        out=pipe(Path(run_tmpdir, "reads.paired.fq")),
        outs=temp(Path(run_tmpdir, "reads.unpaired.fq")),
    log:
        Path(logdir, "repair.log"),
    threads: 1
    resources:
        time=lambda wildcards, attempt: 10 * attempt,
        mem_mb=lambda wildcards, attempt: 4e3 * attempt,
    container:
        bbmap
    shell:
        "repair.sh "
        "-Xmx{resources.mem_mb}m "
        "-Xms100m "
        "in=stdin.fastq "
        "int=t "
        "out=stdout.fastq "
        "outs={output.outs} "
        "repair=t "
        "tossbrokenreads=t "
        "tossjunk=t "
        "< {input} "
        ">> {output.out} "
        "2> {log}"


rule bam2fastq:
    input:
        rnaseq,
    output:
        pipe=pipe(Path(run_tmpdir, "reads.fq")),
    log:
        Path(logdir, "bam2fastq.log"),
    threads: 1
    resources:
        time=lambda wildcards, attempt: 10 * attempt,
        mem_mb=lambda wildcards, attempt: 4e3 * attempt,
    container:
        bbmap
    shell:
        "reformat.sh "
        "-Xmx{resources.mem_mb}m "
        "-Xms100m "
        "in={input} "
        "out=stdout.fq "
        "trimrname=t "
        "tossjunk=t "
        "primaryonly=t "
        ">> {output.pipe} "
        "2> {log}"
