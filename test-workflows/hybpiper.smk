#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

sample_data = Path(
    'test-data',
    'hybpiper',
    'samples.csv'
    )
target_file = Path(
    'test-data',
    'hybpiper',
    'combined_targetfiles.fixed.fasta'
    )
read_directory = Path(
    'test-data',
    'hybpiper',
    'sample_reads'
    )
output_directory = Path(
    'test-output',
    'hybpiper',
    )


samples = pd.read_csv(sample_data, index_col='name')
all_samples = sorted(set(samples.index))

# set up the hybpiper run
if 'hybpiper' not in config.keys():
    config['hybpiper'] = {}

# configure the run like this, or in a yaml file
hybpiper_config = config['hybpiper']

hybpiper_config['sample_list'] = all_samples
hybpiper_config['read_directory'] = read_directory
hybpiper_config['target_file'] = target_file
hybpiper_config['outdir'] = output_directory

config['hybpiper'] = hybpiper_config


module hybpiper:
    snakefile:
        'modules/hybpiper/Snakefile'
        # github(
        #     "tomharrop/smk-modules",
        #     path="modules/hybpiper/Snakefile",
        #     commit="043ebde"
        # )
    config:
        config['hybpiper']

use rule * from hybpiper
