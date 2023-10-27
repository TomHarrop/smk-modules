#!/usr/bin/env python3

import uuid

config_dict = {"guid": uuid.getnode(), "run_tmpdir": snakemake.params["run_tmpdir"]}

with open(snakemake.input[0], "r") as template:
    template_content = template.read()
    formatted_content = template_content.format(**config_dict)

with open(snakemake.output[0], "w") as output:
    output.write(formatted_content)
