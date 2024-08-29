# Establishing a MXPV wastewater sequencing workflow

The goal of this project to identify and track MPXV clades and lineages in wastewater. Viral load will be determined using qPCR or ddPCR, while limited amplicon sequencing is used to identify MPXV clades and linages.

## qPCR/ddPCR

Pan-MPXV primers and probes used for viral quantification (from [Li et al](https://doi.org/10.1016/j.jviromet.2010.07.012)):

| 	Oligos	                                    | sequences			 						                   |
|:----------------------------------------------|---------------------------------------------:|
| G2R_G_forward                                 |  5′-GGAAAATGTAAAGACAACGAATACAG 				       |
| G2R_G_reverse                                	|  5′-GCTATCACATAATCTGGAAGCGTA 				         |
| G2R_G_probe                                	  |  5′FAM-AAGCCGTAATCTATGTTGTCTATCGTGTCC-3′BHQ1 |


## Amplicon sequencing

To limit the number of amplicons, we first determine the minimum set of nucleotide positions needed to discriminate between all clades/lineages
```sh
python ./scripts/minimize_barcodes.py
```
This results in a list of 36 nucleotide positions (minimal_positions.csv) that cover all lineages (minimal_barcodes.csv).

Next, we search for forward and reverse primers around these nucleotide positions.
```sh
python ./scripts/find_primers.py
```
These primer sets are non-overlapping. Two amplicons span 2 nucleotide positions, reducing the total number of amplicons to 34 (primers_output.csv).

We then cross-check in the alignment if the primer-binding locations contain mutations that could lead to amplicon drop out.
```sh
python ./scripts/primer_variation_check.py
```
Most primer-binding sites show no or very limited diversity. Only, one mismatch in primer 25644 was resolved by introducing an ambiguous nucleotide (primer_variations.csv).

## To do
- add additional amplicons in variable regions
