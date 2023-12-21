#!/usr/bin/env python3

# This module gets a read directory, and a csv file to wait for. Once the flag file is ready, it gets the paths from the csv to generate a sample dict, and passes tthat to the other STAR moduel


module star:
    snakefile:
        "Snakefile"
    config:
        {
            "reference": reference,
            "annotation": annotation,
            "outdir": test_outdir,
            
        }


checkpoint get_samples:
    input:
        config["sample_data"],
    output:
        touch(".samples.ready"),
