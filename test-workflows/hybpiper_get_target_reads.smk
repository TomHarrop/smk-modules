#!/usr/bin/env python3

from pathlib import Path
import re

test_data = Path("test-data", "hybpiper_get_target_reads")
test_outdir = Path("test-output", "hybpiper_get_target_reads")
Path("data", "read_files")


run_tmpdir = Path(test_outdir, "tmp")


my_snakefile = github(
    "tomharrop/smk-modules",
    path="modules/hybpiper_get_target_reads/Snakefile",
    tag="0.0.44",
)
# my_snakefile = "../modules/hybpiper_get_target_reads/Snakefile"


all_samples = sorted(
    set(x.name.split(".")[0] for x in test_data.glob("*.tar.gz"))
)

###########
# MODULES #
###########


rule target:
    input:
        expand(
            Path(
                test_outdir,
                "sample_files",
                "{sample}",
                "{sample}.done",
            ),
            sample=all_samples,
        ),


module hybpiper_get_target_reads:
    snakefile:
        my_snakefile
    config:
        {
            "archive_directory": test_data,
            "outdir": test_outdir,
            "raw_read_dir": Path(run_tmpdir, "read_files"),
            "run_tmpdir": run_tmpdir,
        }


use rule * from hybpiper_get_target_reads as hybpiper_get_target_reads_*


# stick to the lower-case naming convention
def get_input(wildcards):
    sample_name = wildcards.sample
    readfile_name = re.sub("subsample1m_", "", sample_name)
    return (
        Path(
            test_data, "read_files", f"{readfile_name}.r{wildcards.r}.fastq.gz"
        ),
    )


rule collect_raw_reads:
    input:
        get_input,
    output:
        Path(run_tmpdir, "read_files", "{sample}.r{r}.fastq.gz"),
    shell:
        "ln -s "
        "$(readlink -f {input} ) "
        "$(readlink -f {output} ) "
