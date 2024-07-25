#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool
from pathlib import Path
from snakemake.logging import logger
import logging
import tarfile
import tempfile


def best_orientation(seq1, seq2):
    aligner = PairwiseAligner()
    orientations = [
        ("forward", seq2),
        ("reverse", seq2[::-1]),
        ("forward complement", seq2.complement()),
        ("reverse complement", seq2.reverse_complement()),
    ]
    best_score = float("-inf")
    best_orientation = None
    for name, oriented_seq in orientations:
        score = aligner.score(seq1, oriented_seq)
        if score > best_score:
            best_score = score
            best_orientation = oriented_seq
            best_name = name
    return best_orientation, best_name


def process_fasta(input_fasta):
    logger.info(f"Processing FASTA file {input_fasta}")
    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as tmp_output_fasta:
        logger.info(f"Writing output for {input_fasta} to {tmp_output_fasta.name}")
        best_orientations = reorient_sequences(input_fasta, tmp_output_fasta.name)
        return tmp_output_fasta.name, input_fasta.name, best_orientations


def reorient_sequences(input_fasta, output_fasta):
    input_handle = open(input_fasta, "r")
    output_handle = open(output_fasta, "w")
    sequences = SeqIO.parse(input_handle, "fasta")
    # record the orientations
    best_orientations = {}
    # Get the first sequence as ref and write it to the output. By definition
    # this is the forward orientation
    reference_seq = next(sequences)
    SeqIO.write(reference_seq, output_handle, "fasta")
    best_orientations[reference_seq.id] = "forward"
    for seq_record in sequences:
        best_oriented_seq, best_name = best_orientation(
            reference_seq.seq, seq_record.seq
        )
        best_orientations[seq_record.id] = best_name
        reoriented_record = SeqRecord(
            best_oriented_seq, id=seq_record.id, description=seq_record.description
        )
        SeqIO.write(reoriented_record, output_handle, "fasta")
    input_handle.close()
    output_handle.close()
    return best_orientations


def main():
    file_handler = logging.FileHandler(snakemake.log[0])
    logger.logfile_handler = file_handler
    logger.logger.addHandler(logger.logfile_handler)

    input_directory = snakemake.input[0]
    output_tarfile = snakemake.output["tarfile"]
    output_csv = snakemake.output["csv"]
    n_threads = snakemake.threads

    logger.info(f"Looking for fasta files in {input_directory}")
    fasta_files = [
        file
        for file in Path(input_directory).glob("*.fasta")
        if not file.name.startswith(".")
    ]

    logger.info(f"Processing {len(fasta_files)} fasta files with {n_threads} processes")
    with Pool(processes=n_threads) as pool:
        results = pool.map(process_fasta, fasta_files)

    logger.info(
        f"Finished processing {len(results)} files"
    )
    with tarfile.open(output_tarfile, "w") as tar, open(output_csv, "wt") as csv:
        csv.write("filename,sequence_id,orientation\n")
        for tmp_output_fasta, original_name, best_orientations in results:
            arcname = original_name
            logger.info(f"Collecting results for {original_name}")
            logger.info(f"Adding {tmp_output_fasta} to {output_tarfile}")
            tar.add(tmp_output_fasta, arcname=original_name)
            logger.info(f"Recording orientations in {output_csv}")
            for k, v in best_orientations.items():
                csv.write(f"{original_name},{k},{v}\n")

    logger.info("\nDone")

if __name__ == "__main__":
    main()
