# MPXV wastewater monitoring

The goal of this project to identify and track mpox clades and lineages in wastewater. A detailed guide to mpox wastewater monitoring can be found [here](https://docs.google.com/document/d/1nAviDZ1zYuch_lw5V3apHTP8yuYBijM-YpmJa1pnH38/edit?usp=sharing). Briefly, we will determine viral load using qPCR or ddPCR, followed by limited amplicon sequencing to identify mpox clades and lineages.

## qPCR/ddPCR

Details for mpox qPCR can be found [here](https://docs.google.com/document/d/1gJjlMpCT4WedxmWK1IeYNG631m9tKNoL1fQC2yteVvw/edit?usp=sharing)

## Amplicon sequencing

We have developed a limited amplicon panel to track mpox in wastewater, which is currently undergoing validation. Full details of the Illumina and ONT sequencing protocols, including primer sequences and PCR conditions, can be found [here](https://docs.google.com/document/d/16wW9Z8x56SuXh_-_rdKi04XwALk-fom-WNRF0cfFlOA/edit?usp=sharing) and [here](https://docs.google.com/document/d/18B24xHGjn-GOMYlefXGzfzYsH5kJM3vz_ndUm467Yq8/edit?usp=sharing). Both protocols use the same MPXV primer sets, but use a different barcoding system.

### Amplicon selection

Limited amplicon sequencing targets higly informative nucleotide sites in the virus genome, and thereby ensure cost-effective and sensitive sequencing method. To limit the number of amplicons in the reactions, we first have determined the minimum set of nucleotide positions needed to discriminate between all clades/lineages as annotated in the [NextStrain tree](https://nextstrain.org/mpox/all-clades). For this we excluded all homoplasic sites and excluded sites within 1000 nucleotides of the 5' and 3' end as relatively few genomes have coverage in those regions.
```sh
python ./scripts/minimize_barcodes.py
```
This results in a list of 45 nucleotide positions (minimal_positions.csv) that cover all currently assigned lineages (minimal_barcodes.csv) in addition to some variation in Clade Ia and Ib (no lineages have been assigned for these clades yet). 

Next, we create amplicons by searching for forward and reverse primers around these 45 nucleotide positions.
```sh
python ./scripts/find_primers.py
```
Several amplicons span multiple nucleotide positions, reducing the total number of amplicons to 38 (mpxv_output.csv). Amplicons are non-overlapping so the assay can be run in a single tube.

We did a quick cross-check in the alignment if the primer-binding locations contain mutations that could lead to amplicon drop-out.
```sh
python ./scripts/primer_variation_check.py
```
Most primer-binding sites showed no or very limited diversity, and the onces that did were manually corrected (mpxv_variations.csv).

## To do
- Add crAssPhage and Synthetic primers.
