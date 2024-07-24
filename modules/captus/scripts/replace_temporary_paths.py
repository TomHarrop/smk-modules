#!/usr/bin/env python3

import json
from pathlib import Path


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
    tmp_json_path = snakemake.input[0]
    fixed_json_path = snakemake.output[0]
    cluster_path = snakemake.params["correct_path"][0]

    keys_to_replace = ["NT_path", "NT_msg"]
    full_cluster_path = Path(cluster_path).resolve().as_posix()

    data = load_json(tmp_json_path)
    replace_keys(data, keys_to_replace, full_cluster_path)

    with open(fixed_json_path, "w") as f:
        json.dump(data, f, indent=4)


if __name__ == "__main__":
    main()
