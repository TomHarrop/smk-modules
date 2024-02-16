#!/usr/bin/env python3

from Bio import SeqIO
from pathlib import Path
import logging

def fastq_streamer(read_file, handles):
    """
    Iterate through the reads. Look up each read id in the handles dict, and
    write it to the appropriate file.
    """
    fastq_iterator = SeqIO.parse(read_file, "fastq")
    for seq_rec in fastq_iterator:
        try:
            write_targets = read_to_target[seq_rec.id]
            for write_target in write_targets:
                SeqIO.write(seq_rec, handles[write_target], "fastq")
        except KeyError as e:
            logging.debug(f"Read {seq_rec.id} skipped")
    for handle in handles.values():
        handle.close()


# inputs
read_lists = snakemake.input["read_lists"]
r1 = snakemake.input["r1"]
r2 = snakemake.input["r2"]

# params from snakemake
my_outdir = Path(snakemake.output["outdir"])
my_sample = snakemake.wildcards["sample"]

# Which target was each read used for. NOTE! An individual read can have more
# than one target. Well done hybpiper.
read_to_target = {}
all_targets = []

for read_list in read_lists:
    target_name = Path(read_list).stem.rstrip("_reads")
    all_targets.append(target_name)
    with open(read_list, "rt") as f:
        for line in f.readlines():
            read = line.rstrip("\n")
            if read in read_to_target:
                read_to_target[read].append(target_name)
            else:
                read_to_target[read] = [target_name]

# Generate handles for each target.
all_targets = sorted(set(all_targets))
r1_outfiles = {x: Path(my_outdir, x, f"{my_sample}.r1.fastq") for x in all_targets}
r2_outfiles = {x: Path(my_outdir, x, f"{my_sample}.r2.fastq") for x in all_targets}

for outfile in sorted(set(list(r1_outfiles.values()) + list(r2_outfiles.values()))):
    outfile.parent.mkdir(parents=True, exist_ok=True)

r1_handles = {x: open(r1_outfiles[x], "wt") for x in r1_outfiles}
r2_handles = {x: open(r2_outfiles[x], "wt") for x in r2_outfiles}

read_file = r1
handles = r1_handles

# Stream the reads
fastq_streamer(r1, r1_handles)
fastq_streamer(r2, r2_handles)
