#!/usr/bin/env python3


captus = "docker://quay.io/biocontainers/captus:1.0.1--pyhdfd78af_2"

# WHAT HAPPENS IF LOCI ARE DUPLICATED IN THE TARGET FILE?
# NC_045138.2 has the most hits so could use that. First try NC_045147.1, which has a hit for 5858, which is the first locus in the target file
target_file_with_duplicated_loci = Path(
    "test-data", "captus", "target_file_with_duplicated_loci.fasta"
)
single_chr_fasta_file = Path("test-data", "captus", "NC_045147.1.fasta")
outdir = Path("test-output", "captus_with_duplicate_loci")
logdir = Path(outdir, "logs")


rule captus_extract:
    input:
        external_fasta=single_chr_fasta_file,
        targets=target_file_with_duplicated_loci,
    output:
        outdir=directory(
            Path(
                outdir,
                "03_extractions",
            )
        ),
    log:
        Path(logdir, "extract.log"),
    benchmark:
        Path(logdir, "benchmark.extract.log")
    threads: lambda wildcards, attempt: 32
    resources:
        time=lambda wildcards, attempt: "5-00",
        mem_mb=lambda wildcards, attempt: 32e3,
    shadow:
        "minimal"
    container:
        captus
    shell:
        "captus_assembly extract "
        "--captus_assemblies_dir 02_assemblies "
        "--fastas {input.external_fasta} "
        "--out {output.outdir}/. "
        "--nuc_refs {input.targets} "
        '--ram "$(( {resources.mem_mb}/1000 ))" '
        "--threads {threads} "
        "--concurrent 3 "
        "&> {log} "
