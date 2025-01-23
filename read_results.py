import os
import glob
import argparse
import pandas as pd
import logging

parser = argparse.ArgumentParser(prog='results_processing', description='Get directory for data processing')
parser.add_argument('results_dir',help='directory of result files')
args = parser.parse_args()

files = glob.glob(f'{args.results_dir}/*.pdbqt',recursive=True)
filenames = []
affinities = []

for file in files:
    output = open(file, 'r')
    lines = output.readlines()
    line = lines[1]
    filenames.append(file)
    affinities.append(float(line.split()[3]))

df = pd.DataFrame({"file":filenames, "affinities":affinities})
# df = pd.Series(data).to_frame()
df.to_csv('Docking_affinities.csv')
print(df) 