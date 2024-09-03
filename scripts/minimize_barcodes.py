import pandas as pd
from itertools import combinations

def get_unique_combinations(df, columns):
    """
    Get the unique combinations of values for the given columns.
    """
    return df[columns].drop_duplicates()

def find_greedy_minimum_sites(df):
    """
    Find a minimal set of nucleotide sites that can distinguish between all MPX lineages using a greedy approach.
    """
    # Extract the columns and lineages
    columns = df.columns
    all_lineages = df.index.tolist()
    
    # Initialize the set of selected sites
    selected_sites = []
    
    # Iterate until all lineages are distinguishable
    while True:
        best_site = None
        best_unique_count = 0
        
        for site in columns:
            if site in selected_sites:
                continue
            
            # Check if adding this site helps
            current_sites = selected_sites + [site]
            unique_combinations = get_unique_combinations(df, current_sites)
            
            if len(unique_combinations) > best_unique_count:
                best_unique_count = len(unique_combinations)
                best_site = site
        
        if best_site is None:
            break  # Exit if no improvement
        
        selected_sites.append(best_site)
        
        # Check if all lineages are distinguishable
        if len(get_unique_combinations(df, selected_sites)) == len(all_lineages):
            break
    
    return selected_sites

def modify_site_names(site_names):
    """
    Remove the first and last letters from each site name.
    """
    return [site[1:-1] for site in site_names]

# File paths
file_path = './data/barcodes_mpxv_merged_no_homoplasic_sites_v3.csv'
output_min_sites_file = './minimal_sites_v3.csv'
output_min_barcodes_file = './minimal_barcodes_v3.csv'

# Load the data
df = pd.read_csv(file_path, index_col=0)

# Find minimal nucleotide sites
min_sites = find_greedy_minimum_sites(df)

# Save the original minimal sites for the barcodes
min_barcodes_df = df[min_sites]
min_barcodes_df.to_csv(output_min_barcodes_file)

# Modify the site names for the minimal sites output
modified_sites = modify_site_names(min_sites)

# Save the modified minimal sites to a CSV
output_sites_df = pd.DataFrame({'Minimal Sites': modified_sites})
output_sites_df.to_csv(output_min_sites_file, index=False)

print(f"Modified minimal set of nucleotide sites saved to: {output_min_sites_file}")
print(f"Minimal barcodes saved to: {output_min_barcodes_file}")
