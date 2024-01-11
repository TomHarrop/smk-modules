#!/usr/bin/env python3

# docker pull quay.io/biocontainers/hisat2:2.2.1--hdbdd923_6

from pathlib import Path
import csv

test_data = Path("test-data", "star")
test_outdir = Path("test-output", "hisat2")

run_tmpdir = Path(test_outdir, "tmp")

# The module waits for the sample_csv file to exist so you can use
# checkpoints.  See the generate_csv checkpoint below.
sample_csv = Path(run_tmpdir, "samples.csv")

# hisat2 requires GTF (not gff)
reference = Path(test_data, "NT_033779.5.fna.gz")
annotation = Path(test_data, "NT_033779.5.gtf.gz")

# hisat2_snakefile = "../modules/hisat2/Snakefile"
hisat2_snakefile = github(
    "tomharrop/smk-modules",
    path="modules/hisat2/Snakefile",
    tag="0.0.32",
)


# outer keys are sample names, inner keys must be r1 and r2
sample_dict = {
    "GSM461177": {
        "r1": Path(test_data, "GSM461177_subsampled.r1.fq.gz"),
        "r2": Path(test_data, "GSM461177_subsampled.r2.fq.gz"),
    },
    "GSM461180": {
        "r1": Path(test_data, "GSM461180_subsampled.r1.fq.gz"),
        "r2": Path(test_data, "GSM461180_subsampled.r2.fq.gz"),
    },
}

###########
# MODULES #
###########

# these have to come before local rules or they don't get added to rules


module hisat2_withannot:
    snakefile:
        hisat2_snakefile
    config:
        {
            "reference": reference,
            "annotation": annotation,
            "outdir": Path(test_outdir, "withannot"),
            "run_tmpdir": Path(run_tmpdir, "withannot"),
            "sample_csv": sample_csv,
        }


use rule * from hisat2_withannot as hisat2_withannot_*


module hisat2_noannot:
    snakefile:
        hisat2_snakefile
    config:
        {
            "reference": reference,
            "outdir": Path(test_outdir, "noannot"),
            "run_tmpdir": Path(run_tmpdir, "noannot"),
            "sample_csv": sample_csv,
        }


use rule * from hisat2_noannot as hisat2_noannot_*


#########
# RULES #
#########


rule target:
    input:
        rules.hisat2_withannot_target.input,
        rules.hisat2_noannot_target.input,
    default_target: True


checkpoint generate_csv:
    output:
        sample_csv,
    run:
        with open(sample_csv, "w", newline="") as f:
            fieldnames = ["sample_name", "r1", "r2"]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            # Write the header
            writer.writeheader()
            # Write each row
            for sample_name, read_paths in sample_dict.items():
                row = {
                    "sample_name": sample_name,
                    "r1": read_paths.get("r1", ""),
                    "r2": read_paths.get("r2", ""),
                }
                writer.writerow(row)
