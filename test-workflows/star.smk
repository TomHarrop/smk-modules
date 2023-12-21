#!/usr/bin/env python3

from pathlib import Path
import csv

test_data = Path("test-data", "star")
test_outdir = Path("test-output", "star")

run_tmpdir = Path(test_outdir, "tmp")

# The star module waits for the sample_csv file to exist so you can use
# checkpoints.  See the generate_csv checkpoint below.
sample_csv = Path(run_tmpdir, "samples.csv")

reference = Path(test_data, "NT_033779.5.fna.gz")
annotation = Path(test_data, "NT_033779.5.gff.gz")

# star_snakefile = "../modules/star/Snakefile"
star_snakefile = github(
    "tomharrop/smk-modules",
    path="modules/star/Snakefile",
    tag="0.0.27",
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


module star:
    snakefile:
        star_snakefile
    config:
        {
            "reference": reference,
            "annotation": annotation,
            "outdir": test_outdir,
            "run_tmpdir": run_tmpdir,
            "sample_csv": sample_csv,
        }


use rule * from star as star_*


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
