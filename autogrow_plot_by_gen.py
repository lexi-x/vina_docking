# import matplotlib.pyplot as plt
# import argparse
# import pandas as pd
# import os

# # Use argparse to allow multiple input files
# parser = argparse.ArgumentParser(description='Plot average binding affinity across generations from multiple docking runs.')
# parser.add_argument('files', nargs='+', help='CSV files to process')
# args = parser.parse_args()

# # Plot setup
# plt.figure(figsize=(10, 6))

# for file in args.files:
#     df = pd.read_csv(file)

#     # Ensure the required columns are present
#     if df.shape[1] < 4:
#         print(f"Skipping {file}: Not enough columns.")
#         continue

#     # Extract Generation and Affinity Score columns
#     generations = df.iloc[:, 3]
#     affinities = df.iloc[:, 2]

#     # Group by generation and compute mean affinity
#     gen_avg = df.groupby(df.columns[3])[df.columns[2]].mean()

#     # Sort by generation for line plot consistency
#     gen_avg = gen_avg.sort_index()

#     # Plot
#     label = os.path.splitext(os.path.basename(file))[0]  # Use filename as label
#     plt.plot(gen_avg.index, gen_avg.values, marker='o', label=label)

# # Final plot formatting
# plt.xlabel("Generation")
# plt.ylabel("Average Best 50 Binding Affinity")
# plt.title("Average Binding Affinity of 50 Best Compounds Per Generation")
# plt.legend(title="Docking Runs", loc='best')
# plt.grid(True)
# plt.tight_layout()

# plt.savefig("average_affinity_plot.png", dpi=300)
# print("Plot saved as average_affinity_plot.png")

import matplotlib.pyplot as plt
import argparse
import pandas as pd
import os

# Use argparse to allow directory input
parser = argparse.ArgumentParser(description='Plot average binding affinity across generations from docking runs in a directory.')
parser.add_argument('directory', type=str, help='Directory containing CSV files to process')
args = parser.parse_args()

# Ensure the provided directory exists
if not os.path.isdir(args.directory):
    print(f"The directory {args.directory} does not exist.")
    exit()

# List all CSV files in the directory
csv_files = [f for f in os.listdir(args.directory) if f.endswith('.csv')]

# Check if there are any CSV files in the directory
if not csv_files:
    print(f"No CSV files found in the directory {args.directory}.")
    exit()

# Plot setup
plt.figure(figsize=(10, 6))

for file in csv_files:
    file_path = os.path.join(args.directory, file)
    df = pd.read_csv(file_path)

    # Ensure the required columns are present
    if df.shape[1] < 4:
        print(f"Skipping {file}: Not enough columns.")
        continue

    # Extract Generation and Affinity Score columns
    generations = df.iloc[:, 3]
    affinities = df.iloc[:, 2]

    # Group by generation and compute mean affinity
    gen_avg = df.groupby(df.columns[3])[df.columns[2]].mean()

    # Sort by generation for line plot consistency
    gen_avg = gen_avg.sort_index()

    # Plot
    label = os.path.splitext(file)[0]  # Use filename as label
    plt.plot(gen_avg.index, gen_avg.values, marker='o', label=label)

# Final plot formatting
plt.xlabel("Generation")
plt.ylabel("Average Best 50 Binding Affinity")
plt.title("Average Binding Affinity of 50 Best Compounds Per Generation")
plt.legend(title="Docking Runs", loc='best')
plt.grid(True)
plt.tight_layout()

# Save the plot as a PNG image
plt.savefig("average_affinity_plot.png", dpi=300)
print("Plot saved as average_affinity_plot.png")