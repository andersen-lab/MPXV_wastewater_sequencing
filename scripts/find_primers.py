import primer3
import pandas as pd
import csv
from Bio import SeqIO

def design_primers(reference_seq, positions, product_size_range=(225, 285), annealing_temp=60.0):
    primers = []

    for position in positions:
        # Define the target region around the given position
        start = max(0, position - product_size_range[1])
        end = min(len(reference_seq), position + product_size_range[1])
        
        # Adjust the target region if necessary to meet the minimum product size
        if end - start < product_size_range[0]:
            if start > 0:
                end = min(len(reference_seq), start + product_size_range[0])
            else:
                end = min(len(reference_seq), product_size_range[0])
        
        target_region = reference_seq[start:end]
        target_relative_position = position - start

        primer_design = primer3.design_primers(
            {
                'SEQUENCE_ID': 'target',
                'SEQUENCE_TEMPLATE': target_region,
                'SEQUENCE_INCLUDED_REGION': [0, len(target_region)],
                'SEQUENCE_TARGET': [(target_relative_position, 1)]
            },
            {
                'PRIMER_OPT_SIZE': 20,
                'PRIMER_PICK_INTERNAL_OLIGO': 1,
                'PRIMER_INTERNAL_MAX_SELF_END': 8,
                'PRIMER_MIN_SIZE': 18,
                'PRIMER_MAX_SIZE': 35,
                'PRIMER_OPT_TM': annealing_temp,
                'PRIMER_MIN_TM': annealing_temp - 2.0,
                'PRIMER_MAX_TM': annealing_temp + 2.0,
                'PRIMER_MIN_GC': 20.0,
                'PRIMER_MAX_GC': 80.0,
                'PRIMER_PRODUCT_SIZE_RANGE': [product_size_range],
                'PRIMER_NUM_RETURN': 1
            }
        )

        if 'PRIMER_LEFT_0_SEQUENCE' in primer_design and 'PRIMER_RIGHT_0_SEQUENCE' in primer_design:
            forward_start = start + primer_design['PRIMER_LEFT_0'][0]
            forward_end = forward_start + primer_design['PRIMER_LEFT_0'][1]
            
            amplicon_size = primer_design['PRIMER_PAIR_0_PRODUCT_SIZE']
            reverse_end = start +  primer_design['PRIMER_RIGHT_0'][0]
            reverse_start = reverse_end - primer_design['PRIMER_RIGHT_0'][1]
            
            primers.append({
                'Position': position,
                'Forward_Primer': primer_design['PRIMER_LEFT_0_SEQUENCE'],
                'Reverse_Primer': primer_design['PRIMER_RIGHT_0_SEQUENCE'],
                'Forward_Start': forward_start,
                'Forward_End': forward_end,
                'Reverse_Start': reverse_end,
                'Reverse_End': reverse_start,
                'Amplicon_Size': amplicon_size,
                'Forward_Primer_Tm': primer_design['PRIMER_LEFT_0_TM'],
                'Reverse_Primer_Tm': primer_design['PRIMER_RIGHT_0_TM']
            })
        else:
            print(f"Warning: Primer design failed for position {position}. Skipping.")

    return primers


def read_fasta(fasta_file):
    with open(fasta_file, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            return str(record.seq)

def read_positions_from_csv(csv_file):
    positions = []
    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # Skip header if there's one
        for row in reader:
            positions.append(int(row[0]))  # Assuming positions are in the first column
    return positions

def write_primers_to_csv(primers, output_file):
    df = pd.DataFrame(primers)
    df.to_csv(output_file, index=False)

def write_bed_file(primers, bed_file, chromosome='chr1'):
    """
    Write the primer start and end positions to a BED file.
    """
    with open(bed_file, 'w') as f:
        for primer in primers:
            forward_start = primer['Forward_Start']
            forward_end = primer['Forward_End']
            reverse_start = primer['Reverse_Start']
            reverse_end = primer['Reverse_End']
            
            # Write forward primer information
            f.write(f"{chromosome}\t{forward_start}\t{forward_end}\tforward_primer_{primer['Position']}\n")
            
            # Write reverse primer information
            f.write(f"{chromosome}\t{reverse_start}\t{reverse_end}\treverse_primer_{primer['Position']}\n")


# Example usage
fasta_file = "./data/mpxv_reference.fasta"
positions_csv = "./minimal_sites_v3_added_Iab_seq.csv"
output_csv = "./primers_output_v2.csv"
output_bed = "./primers_output_v2.bed"

reference_sequence = read_fasta(fasta_file)
positions = read_positions_from_csv(positions_csv)
primers = design_primers(reference_sequence, positions)

write_primers_to_csv(primers, output_csv)
write_bed_file(primers, output_bed)

print(f"Primers have been written to {output_csv} and {output_bed}")

