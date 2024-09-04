# MPXV wastewater qPCR & sequencing

The goal of this project to identify and track MPXV clades and lineages in wastewater. Viral load will be determined using qPCR or ddPCR, while limited amplicon sequencing is used to identify MPXV clades and lineages.

## qPCR/ddPCR

Pan-MPXV primers and probes used for viral quantification (from [Li et al.](https://doi.org/10.1016/j.jviromet.2010.07.012)):

| Oligo names	                                  |  Sequences			 						                 |
|:----------------------------------------------|---------------------------------------------:|
| G2R_G_forward                                 |  5′-GGAAAATGTAAAGACAACGAATACAG 				       |
| G2R_G_reverse                                	|  5′-GCTATCACATAATCTGGAAGCGTA 				         |
| G2R_G_probe                                	  |  5′FAM-AAGCCGTAATCTATGTTGTCTATCGTGTCC-3′BHQ1 |


## Amplicon sequencing

To limit the number of amplicons, we first determine the minimum set of nucleotide positions needed to discriminate between all clades/lineages as annotated in the [NextStrain tree](https://nextstrain.org/mpox/all-clades). For this we excluded all homoplasic sites. We also excluded sites within 1000 nts of the 5' and 3' end as relatively few genomes have coverage in those regions.
```sh
python ./scripts/minimize_barcodes.py
```
This results in a list of 44 nucleotide positions (minimal_positions.csv) that cover all currently assigned lineages (minimal_barcodes.csv) and some variation in Clade Ia and Ib (no lineages have been defined for these clades yet). 

Next, we create amplicons by searching for forward and reverse primers around these 44 nucleotide positions.
```sh
python ./scripts/find_primers.py
```
Several amplicons span multiple nucleotide positions, reducing the total number of amplicons to 38 (primers_output.csv). Amplicons are non-overlapping.

We then do a quick cross-check in the alignment if the primer-binding locations contain mutations that could lead to amplicon drop-out.
```sh
python ./scripts/primer_variation_check.py
```
Most primer-binding sites show no or very limited diversity (primer_variations.csv).

## To do
- add additional amplicons in variable regions
