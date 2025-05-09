#!/bin/bash

# Define some basic environmental variables before launching the suite

# Load the analysis3 conda environment
module use /g/data/xp65/public/modules
module load conda/analysis3-25.02

# Root directory for this repo
export ROOT=/scratch/k10/${USER}/hk25-AusNode-ConvTrop	# Change this to match where you clone this repo
export MODULES=${ROOT}/filter_mode

# Append to our python path
export PYTHONPATH=${MODULES}:${PYTHONPATH}
