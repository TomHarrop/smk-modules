#!/usr/bin/env python3


from pathlib import Path
from snakemake.logging import logger
import json
import logging


def find_cluster_file(directory):
    return [
        Path(x).resolve() for x in Path(directory).rglob("*_captus_cluster_refs.fasta")
    ]


def load_json(file_path):
    with open(file_path, "r") as file:
        data = json.load(file)
    return data


def replace_keys(data, keys_to_replace, new_value):
    if "CLR" in data and isinstance(data["CLR"], dict):
        for key in keys_to_replace:
            if key in data["CLR"]:
                data["CLR"][key] = new_value


def main():
    file_handler = logging.FileHandler(snakemake.log[0])
    logger.logfile_handler = file_handler
    logger.logger.addHandler(logger.logfile_handler)

    tmp_json_path = snakemake.input["tmp_json_path"]
    extraction_dir = snakemake.input["extraction_dir"]
    fixed_json_path = snakemake.output["fixed_json"]

    logger.info(f"Searching in {extraction_dir}")
    cluster_files = find_cluster_file(extraction_dir)
    logger.info(f"Found files: {cluster_files}")

    keys_to_replace = ["NT_path", "NT_msg"]
    full_cluster_path = Path(cluster_files[0]).resolve().as_posix()
    logger.info(f"Using full_cluster_path {full_cluster_path}")

    data = load_json(tmp_json_path)
    logger.info(f"\nOriginal json:\n{data}\n")
    replace_keys(data, keys_to_replace, full_cluster_path)
    logger.info(f"Fixed json:\n{data}\n")

    logger.info(f"Writing output to {fixed_json_path}")
    with open(fixed_json_path, "w") as f:
        json.dump(data, f, indent=4)

    logger.info(f"Done")


if __name__ == "__main__":
    main()
