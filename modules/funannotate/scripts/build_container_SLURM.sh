#!/bin/bash

#SBATCH --job-name=ipr_build
#SBATCH --time=2-00
#SBATCH --ntasks=1
#SBATCH --mem=8g
#SBATCH --output=ipr_build.slurm.out
#SBATCH --error=ipr_build.slurm.err
#SBATCH --mail-type=ALL
#SBATCH --export=ALL
#SBATCH --partition=io

# Dependencies
module load apptainer/1.1.5-suid

git clone https://github.com/TomHarrop/container-interproscan.git

VERSION="$(cat container-interproscan/VERSION)"

apptainer build \
    "interproscan_${VERSION}.sif"
    container-interproscan/Singularity
