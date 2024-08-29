import pandas as pd
from Bio import AlignIO
from Bio.Seq import Seq
from collections import defaultdict

def parse_primers_csv(csv_file):
    """
    Parse the CSV file to extract primer positions and sequences.
    """
    primers = pd.read_csv(csv_file)
    return primers

def clean_sequence(seq):
    """
    Remove gaps (-) and N's from a sequence.
    """
    return seq.replace('-', '').replace('N', '')

def get_reverse_complement(seq):
    """
    Return the reverse complement of a DNA sequence.
    """
    return str(Seq(seq).reverse_complement())

def check_primer_variations(msa_file, primers):
    """
    Check the multiple sequence alignment for variations in the primer binding sites.
    """
    alignment = AlignIO.read(msa_file, "fasta")
    num_sequences = len(alignment)  # Total number of sequences in the alignment
    variations = defaultdict(lambda: defaultdict(int))

    for _, primer in primers.iterrows():
        forward_start = primer['Forward_Start']
        forward_end = forward_start + len(primer['Forward_Primer'])
        reverse_start = (primer['Reverse_End'] + 1)
        reverse_end = (primer['Reverse_End'] + len(primer['Reverse_Primer']) + 1)

        # Extract and clean the primer binding site sequences from the alignment
        forward_binding_site = alignment[:, forward_start:forward_end]
        reverse_binding_site = alignment[:, reverse_start:reverse_end]

        forward_reference_seq = clean_sequence(str(forward_binding_site[0].seq))
        reverse_reference_seq = clean_sequence(str(reverse_binding_site[0].seq))

        for i, record in enumerate(alignment):
            forward_seq = clean_sequence(str(forward_binding_site[i].seq))
            reverse_seq = clean_sequence(str(reverse_binding_site[i].seq))

            # Ignore sequences shorter than the reference sequence
            if len(forward_seq) >= len(forward_reference_seq) and forward_seq != forward_reference_seq:
                variations[primer['Position']]['Forward', forward_seq, primer['Forward_Primer']] += 1

            if len(reverse_seq) >= len(reverse_reference_seq) and reverse_seq != reverse_reference_seq:
                rev_complement_seq = get_reverse_complement(reverse_seq)
                variations[primer['Position']]['Reverse', rev_complement_seq, primer['Reverse_Primer']] += 1

    return variations, num_sequences

def write_variations_summary_to_csv(variations, num_sequences, output_csv):
    """
    Write a summary of observed variations to a CSV file with occurrence percentages.
    """
    summarized_variations = []

    # Write the summarized variations with percentages
    for position, variant_dict in variations.items():
        for (strand, variant_seq, primer_seq), count in variant_dict.items():
            percentage = (count / num_sequences) * 100
            summarized_variations.append({
                'Primer': position,
                'Strand': strand,
                'Reference_Primer_Sequence': primer_seq,
                'Variant_Sequence': variant_seq,
                'Occurrence_Count': count,
                'Percentage': percentage
            })

    df = pd.DataFrame(summarized_variations)

    # Sort by Primer and then by Strand
    df = df.sort_values(by=['Primer', 'Strand'])

    df.to_csv(output_csv, index=False)

if __name__ == "__main__":
    primers_csv = "./primers_output.csv"
    msa_file = "./aligned 2.fasta"
    output_csv = "./primer_variations.csv"

    # Parse the CSV file to get primer binding sites
    primers = parse_primers_csv(primers_csv)
    
    # Check the MSA for variations in the primer binding sites
    variations, num_sequences = check_primer_variations(msa_file, primers)
    
    # Write a summary of the variations to a CSV file
    write_variations_summary_to_csv(variations, num_sequences, output_csv)

    print(f"Summary of variations has been written to {output_csv}")


