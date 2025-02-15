import pandas as pd
import re

# Load the CSV file that needs to be filtered
data_file = "./data/barcodes_mpxv_large.csv"
data_df = pd.read_csv(data_file)

# Load the CSV file with the column numbers to be filtered
filter_file = "./data/mpxv_homoplasic_sites.csv"
filter_df = pd.read_csv(filter_file, header=None)

# Extract the numbers from the filter file
filter_numbers = filter_df.iloc[:, 0].astype(str).tolist()

# Function to extract numbers from column headers
def extract_number(column_name):
    match = re.search(r'(\d+)', column_name)
    return match.group(1) if match else None

# Identify columns to drop based on the extracted numbers
columns_to_drop = [col for col in data_df.columns if extract_number(col) in filter_numbers]

# Drop the identified columns
filtered_df = data_df.drop(columns=columns_to_drop)

# Save the filtered dataframe to a new CSV file
filtered_df.to_csv("./data/barcodes_mpxv_large_no_homoplasic_sites.csv", index=False)

print("Columns filtered successfully")
