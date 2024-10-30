# P56_hippocamus_scRNAseq
Analysis of publicly available adult mouse hippocampus single-cell RNA sequence data to extract the excitatory neuron specific genes.

The code was written under R version 4.4.1 and Rstudio version 2023.12.1+402.

To focus on the excitatory neuron specific genes in our bulk RNA data, we extracted the excitatory neuron gene list from the public database https://knowledge.brain-map.org/data/Z0GBA7V12N4J4NNSUHA/summary (Anoushka Joglekar et al., 2024 ) The processed 10X v3 single-cell counts data were downloaded from NeMo database (Identifier: nemo:dat-qwqfftg). We integrated two replicated experiment, P56_M1_HIP and P56_M2_HIP. The pseudobulk counts were extracted from the clusters triple positive for _Camk2a_; _Grin1_ and _Slc17a7_ and merged with our bulk RNA data as described in the result.

Reference:

Joglekar, A., Hu, W., Zhang, B., Narykov, O., Diekhans, M., Marrocco, J., Balacco, J., Ndhlovu, L. C., Milner, T. A., Fedrigo, O., Jarvis, E. D., Sheynkman, G., Korkin, D., Ross, M. E., & Tilgner, H. U. (2024). Single-cell long-read sequencing-based mapping reveals specialized splicing patterns in developing and adult mouse and human brain. Nature neuroscience, 27(6), 1051–1063. https://doi.org/10.1038/s41593-024-01616-4
