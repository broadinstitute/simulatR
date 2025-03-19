# simulatR

This R package simulates the transmission and evolution of viral pathogens. It is used to validate and benchmark `JUNIPER`, a transmission inference tool by Specht et al. (2025), available at doi.org/10.1101/2025.03.02.25323192.

## Installation

To install `simulatR`, first install the `devtools` package via `install.packages("devtools")`. Then, run `devtools::install_github("broadinstitute/simulatR")`.

## Synthetic Outbreak Generation

To simulate an outbreak, load the `simulatR` library, and then run `epi_sim(...)`, where `...` represents optional arguments passed to `epi_sim()`. For details of these arguments, run `?epi_sim` to view documentation. 

The `epi_sim()` function generates a new subdirctory of the current working directory with name specified by the `outdir` argument to `epi_sim()`, which defaults to `my_epidemic`. This directory contains an aligned FASTA file `aligned.fasta	of all sampled cases, a file `ref.FASTA` containing the genome assigned to the index case, a file `transmission.csv` recording who infects whom (for both sampled and unsampled cases) and when (in time units post the infection time of the index case), and, if the argument `write_adjacency` of `epi_sim()` is set to `TRUE`, an adjacency matrix of transmissions between sampled cases. Additionally, this directory contains a subdirectory named `vcf`, which holds one VCF file for each sampled case.
