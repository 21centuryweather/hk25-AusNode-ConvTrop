from pathlib import Path
import os

user_id = os.environ['USER']
ROOT_DIR = os.environ['ROOT']

# Location of input ISCCP OLR directory
OLR_DIR=Path('/scratch/gb02/mr4682/data/txuptp')

# Location of output directory
FILTERED_OLR_DIR=Path('/scratch/gb02/mr4682/data/txuptp/filtered')
