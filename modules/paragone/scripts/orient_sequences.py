#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from snakemake.logging import logger
import logging
import os
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


def reorient_sequences(input_fasta, output_fasta):
    input_handle = open(input_fasta, "r")
    output_handle = open(output_fasta, "w")
    sequences = SeqIO.parse(input_handle, "fasta")
    # Get the first sequence as ref and write it to the output
    reference_seq = next(sequences)
    SeqIO.write(reference_seq, output_handle, "fasta")
    for seq_record in sequences:
        best_oriented_seq, best_name = best_orientation(
            reference_seq.seq, seq_record.seq
        )
        logger.info(f"Found orientation {best_name} for sequence {seq_record.id}")
        reoriented_record = SeqRecord(
            best_oriented_seq, id=seq_record.id, description=seq_record.description
        )
        SeqIO.write(reoriented_record, output_handle, "fasta")
    input_handle.close()
    output_handle.close()


def main():
    file_handler = logging.FileHandler(snakemake.log[0])
    logger.logfile_handler = file_handler
    logger.logger.addHandler(logger.logfile_handler)
    
    input_directory = snakemake.input[0]
    output_tarfile = snakemake.output[0]

    logger.info(f"Looking for fasta files in {input_directory}")

    fasta_files = [
        file
        for file in Path(input_directory).glob("*.fasta")
        if not file.name.startswith(".")
    ]

    with tarfile.open(output_tarfile, "w") as tar:
        for input_fasta in fasta_files:
            logger.info(f"Opening fasta file {input_fasta}")
            with tempfile.NamedTemporaryFile(
                delete=False, suffix=".fasta"
            ) as tmp_output_fasta:
                reorient_sequences(input_fasta, tmp_output_fasta.name)
                tar.add(tmp_output_fasta.name, arcname=input_fasta.name)


if __name__ == "__main__":
    main()
