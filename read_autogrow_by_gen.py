import os
import glob
import pandas as pd
import re
import argparse

# Define a function to extract SMILES string, name, and affinity score
def extract_data_from_line(line):
    # Split the line into parts based on tab characters
    parts = line.strip().split('\t')
    
    # SMILES string is the first element
    smiles = parts[0]
    
    # The name is the second element (which seems to have some kind of identifier)
    name = parts[1]
    
    # The affinity score is the last element (assuming it's always the last one)
    affinity_score = float(parts[4])
    
    return smiles, name, affinity_score

# Function to get the generation number from the file path
def get_generation_from_path(path):
    match = re.search(r'generation_(\d+)', path)
    return int(match.group(1)) if match else None

# Function to process a single file
def process_file(file_path, vset):
    extracted_data = []
    visited = vset
    generation = get_generation_from_path(file_path)
    # Open the file for reading
    with open(file_path, 'r') as file:
        # Read each line in the file
        for line in file:
            # Extract data from the line
            smiles, name, affinity_score = extract_data_from_line(line)
            # Store the extracted data
            if smiles not in vset:
                extracted_data.append((smiles, name, affinity_score, generation))
            vset.add(smiles)
    return extracted_data, vset

# Function to parse a directory and process all files ending with '_ranked.smi'
def process_directory(directory_path):
    # Use glob to find all files ending with '_ranked.smi' in the directory and subdirectories
    file_pattern = os.path.join(directory_path, '**', '*_ranked.smi')
    files = glob.glob(file_pattern, recursive=True)
    visited = set()
    
    # Prepare a list to hold all extracted data
    all_data = []

    # Loop through all found files
    for file_path in files:
        # print(f"Processing file: {file_path}")
        extracted_data, visited = process_file(file_path, visited)
        
        # Add the extracted data to the list
        all_data.extend(extracted_data)

    # If data was extracted, create and return the DataFrame
    df = pd.DataFrame(all_data, columns=['SMILES', 'Name', 'Affinity Score', 'Generation'])
    
    return df


# Main function to execute the extraction on a given directory
def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Process docking run files and generate a sorted CSV of the best 50 compounds per generation.')
    parser.add_argument('input_directory', type=str, help='Input directory containing docking run files')
    parser.add_argument('output_csv', type=str, help='Output CSV file to save the results')
    
    args = parser.parse_args()
    
    # Process the directory and get the DataFrame with the extracted data
    df = process_directory(args.input_directory)
    
    # Sorting the dataframe to select the top 50 compounds for each generation
    df_sorted = df.groupby(df.columns[3], group_keys=False).apply(
        lambda x: x.sort_values(by=x.columns[2], ascending=True).head(50)
    )

    if df_sorted is not None:
        # Store the DataFrame in a CSV file
        df_sorted.to_csv(args.output_csv, index=False)
        print(f"Data has been saved to {args.output_csv}")
    else:
        print("No data to save.")

if __name__ == '__main__':
    main()
