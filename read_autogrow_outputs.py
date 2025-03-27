import os
import glob
import pandas as pd

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

# Function to process a single file
def process_file(file_path):
    extracted_data = []
    
    # Open the file for reading
    with open(file_path, 'r') as file:
        # Read each line in the file
        for line in file:
            # Extract data from the line
            smiles, name, affinity_score = extract_data_from_line(line)
            # Store the extracted data
            if affinity_score < -7.5:
                extracted_data.append((smiles, name, affinity_score))
    return extracted_data

# Function to parse a directory and process all files ending with '_ranked.smi'
def process_directory(directory_path):
    # Use glob to find all files ending with '_ranked.smi' in the directory and subdirectories
    file_pattern = os.path.join(directory_path, '**', '*_ranked.smi')
    files = glob.glob(file_pattern, recursive=True)
    
    # Prepare a list to hold all extracted data
    all_data = []

    # Loop through all found files
    for file_path in files:
        # print(f"Processing file: {file_path}")
        extracted_data = process_file(file_path)
        
        # Add the extracted data to the list
        all_data.extend(extracted_data)

    # If data was extracted, create and return the DataFrame
    df = pd.DataFrame(all_data, columns=['SMILES', 'Name', 'Affinity Score'])
    
    return df


# Main function to execute the extraction on a given directory
def main():
    # Specify the directory path to process 
    directory_path = '/home/lexi/autogrow4/autogrow_results' 
    output_csv = 'autogrow_best_compounds.csv'  # Output CSV file

    # Process the directory and get the DataFrame with the extracted data
    df = process_directory(directory_path)
    df_sorted = df.sort_values(by='Affinity Score', ascending=True)

    if df is not None:
        # Store the DataFrame in a CSV file
        df_sorted.to_csv(output_csv, index=False)
        print(f"Data has been saved to {output_csv}")
    else:
        print("No data to save.")

if __name__ == '__main__':
    main()