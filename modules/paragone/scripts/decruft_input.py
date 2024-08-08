#!/usr/bin/env python3

from Bio import SeqIO
from multiprocessing import Pool
from pathlib import Path
from snakemake.logging import logger
import logging
import os
import gzip


def list_fasta_files(directory):
    return [f for f in directory.glob("*.fasta") if not f.name.startswith(".")]


def check_empty(files):
    empty_files = []
    non_empty_files = []
    for file in files:
        if os.path.getsize(file) == 0:
            empty_files.append(file)
        else:
            non_empty_files.append(file)
    return empty_files, non_empty_files


def calculate_gap_fraction(seq):
    return (seq.count("-") + seq.count("N")) / len(seq)


def filter_and_write_sequences(file, max_gap_fract, output_directory):
    output_fasta = Path(output_directory, file.name)
    sequences = list(SeqIO.parse(file, "fasta"))
    if len(sequences) < 3:
        logger.warning(
            f"Only {len(sequences)} record(s) in {file}, it will be excluded."
        )
        return file, {}, len(sequences)

    gap_stats = {}
    filtered_sequences = []
    kept_sequences = []

    for sequence in sequences:
        my_gap_frac = calculate_gap_fraction(sequence.seq)
        gap_stats[sequence.id] = my_gap_frac
        if my_gap_frac > max_gap_fract:
            filtered_sequences.append(sequence.id)
        else:
            kept_sequences.append(sequence)

    for seq_id in filtered_sequences:
        logger.warning(
            f"Excluding {seq_id} from {file} with gap_fraction {gap_stats[seq_id]}"
        )

    if len(kept_sequences) < 3:
        logger.warning(
            f"{len(kept_sequences)} record(s) left in {file} after filtering for gappy sequences, it will be excluded."
        )
        return file, gap_stats, len(sequences)

    else:
        with open(output_fasta, "w") as output_handle:
            for sequence in kept_sequences:
                SeqIO.write(sequence, output_handle, "fasta")

    return file, gap_stats, len(filtered_sequences)


def process_fasta_file(args):
    file, max_gap_fract, output_directory = args
    return filter_and_write_sequences(file, max_gap_fract, output_directory)


def main(
    input_directory,
    max_gap_fract,
    output_directory,
    seq_stats,
    removed_files,
    n_threads,
):
    input_directory = Path(input_directory)
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    logger.info(f"Looking for fasta files in {input_directory}")
    fasta_files = list_fasta_files(input_directory)
    empty_files, non_empty_files = check_empty(fasta_files)
    logger.info(f"Found {len(fasta_files)} files")

    args = [(file, max_gap_fract, output_directory) for file in non_empty_files]

    logger.info(f"Processing fasta files with {n_threads} processes")
    with Pool(processes=n_threads) as pool:
        results = pool.map(process_fasta_file, args)
    logger.info(f"Processed {len(results)} files")

    total_filtered_sequences = 0
    files_without_enough_records = []

    with gzip.open(seq_stats, "wt") as csv:
        csv.write("filename,sequence,gap_fraction\n")
        for file, gap_stats, filtered_sequences in results:
            if gap_stats:
                total_filtered_sequences += filtered_sequences
                logger.info(f"Logging stats for {file} to {seq_stats}")
                for k, v in gap_stats.items():
                    csv.write(f"{file},{k},{v}\n")
            else:
                files_without_enough_records.append(file.as_posix())

            if len(gap_stats) == filtered_sequences:
                files_without_enough_records.append(file.as_posix())

    logger.info(f"Recording {len(empty_files)} empty files in {removed_files}")
    with open(removed_files, "w") as ef:
        ef.write("Empty files:\n")
        ef.write("\n".join(empty_files))
        ef.write("\nFiles with less than three records:\n")
        ef.write("\n".join(files_without_enough_records))
        ef.write("\n")

    logger.info(
        f"\nExcluded {len(empty_files)} empty file(s).\n"
        f"Excluded {len(files_without_enough_records)} file(s) with < 3 records.\n"
        f"Removed a total of {total_filtered_sequences} sequences with gap fraction above {max_gap_fract} from the remaining files."
    )


if __name__ == "__main__":
    file_handler = logging.FileHandler(snakemake.log[0])
    logger.logfile_handler = file_handler
    logger.logger.addHandler(logger.logfile_handler)

    input_directory = snakemake.input["input_directory"]
    max_gap_fract = snakemake.params["max_gap_fract"]

    output_directory = snakemake.output["output_directory"]
    seq_stats = snakemake.output["seq_stats"]
    removed_files = snakemake.output["removed_files"]

    n_threads = snakemake.threads

    main(
        input_directory,
        max_gap_fract,
        output_directory,
        seq_stats,
        removed_files,
        n_threads,
    )
