#!/usr/bin/env python

import csv
import sys
import xml.etree.ElementTree as ET

run_xml = ET.parse(snakemake.input["xml"])
root = run_xml.getroot()

run_to_url = {}
for EXPERIMENT_PACKAGE in root.findall("EXPERIMENT_PACKAGE"):
    for RUN_SET in EXPERIMENT_PACKAGE.findall("RUN_SET"):
        for RUN in RUN_SET.findall("RUN"):
            # the XML also has a public AWS url to the full SRA file
            # (not SRAlite)
            for SRAFiles in RUN.findall("SRAFiles"):
                for SRAFile in SRAFiles.findall("SRAFile"):
                    if SRAFile.get("url"):
                        my_filename = SRAFile.get("filename")
                        my_url = SRAFile.get("url")
                        run_to_url[my_filename] = my_url

full_files_only = {k: v for k, v in run_to_url.items() if "lite" not in k}

with open(snakemake.output["csv"], "w") as f:
    for k, v in full_files_only.items():
        f.write(f"{k},{v}\n")
