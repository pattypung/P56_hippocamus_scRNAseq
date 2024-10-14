# P56_hippocamus_scRNAseq
Analysis of public adult mouse hippocampus single-cell RNA sequence data to extract the excitatory neuron pseudo-bulk RNA data. 

The code was written under R version 4.4.1 and Rstudio version 2023.12.1+402. 

To focus on the target involved in Emx1 ineage hippocampus neuron in our bulk RNA data, we extract the excitatory neuron gene list from the public database https://knowledge.brain-map.org/data/Z0GBA7V12N4J4NNSUHA/summary (Anoushka Joglekar et al., 2023)
The processed 10X v3 singel-cell counts data were downloaded from NeMo database (Identifier: nemo:dat-qwqfftg).
We integrated two replicated experiment, P56_M1_HIP and P56_M2_HIP. The pseudobulk counts were than extract from the excitatory marker positive cluster and were merged with our bulk RNA data as described in the result. 

Reference:

Single-cell long-read mRNA isoform regulation is pervasive across mammalian brain regions, cell types, and development
Anoushka Joglekar, Wen Hu, Bei Zhang, Oleksandr Narykov, Mark Diekhans, Jennifer Balacco, Lishomwa C Ndhlovu, Teresa A Milner, Olivier Fedrigo, Erich D Jarvis, Gloria Sheynkman, Dmitry Korkin, M. Elizabeth Ross, Hagen U. Tilgner
bioRxiv 2023.04.02.535281; doi: https://doi.org/10.1101/2023.04.02.535281
