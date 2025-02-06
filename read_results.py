import os
import glob
import argparse
import pandas as pd
import logging

# Set up argument parser 
parser = argparse.ArgumentParser(prog='results_processing', description='Get directory for data processing')
parser.add_argument('results_dir',help='directory of affinity result files')
parser.add_argument('output_name',help='name of output csv file')
args = parser.parse_args()

# Parse results directory for file names
files = glob.glob(f'{args.results_dir}/*.pdbqt',recursive=True)
filenames = []
affinities = []

# Obtain docking score from pdbqt file
for file in files:
    output = open(file, 'r')
    lines = output.readlines()
    line = lines[1]
    filenames.append(file)
    affinities.append(float(line.split()[3]))

# Compile scores and files in csv format
df = pd.DataFrame({"file":filenames, "affinity":affinities})
df_sorted = df.sort_values(by="affinity")
df_sorted.to_csv(args.output_name,index=False)