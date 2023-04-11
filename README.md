# Estimating transmission rates of different *Mtb* lineages in Tanzania using phylodynamics
This repository contains the code associated with the phylodynamic analyses performed in Zwyer, Rutaihwa, Windels et al. (2023) Back-to-Africa introductions of *Mycobacterium tuberculosis* as the main cause of tuberculosis in Dar es Salaam, Tanzania. *PLOS Pathogens* 19, 4 (2023). https://doi.org/10.1371/journal.ppat.1010893

These analyses aimed at estimating the transmission rate of different *Mycobacterium tuberculosis* lineages (introductions), by fitting a birth-death model to genomic sequencing data collected from TB patients in Dar es Salaam, Tanzania (2013-2019).

The raw sequencing data used in this study are available under project accession number [PRJEB49562](https://www.ebi.ac.uk/ena/browser/view/PRJEB49562).

- The folder `analyses/` contains the XML files used for the birth-death analyses in BEAST2.
- The folder `scripts/` contains R scripts used for pre-processing of the alignments and post-processing of the BEAST2 output.
