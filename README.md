# MDS analysis

**Brief**: This R script is used for clustering significant pathways identified by using multidimensional scalling (MDS) and Jaccard distance

**Usage**: When conducting genome-wide RNA expression analysis, GWASs, or systematic genomics integrative analyses, we could identify a number of significant genes associated with traits of interest, including smoking behavior, schizophrenia, major depressive disorder, coronary artery disease and so on. By using these identified genes, many researchers would like to perform a pathway-based enrichment analysis to get insight into biological mechanisms underlying diseases. Generally, there would be numerous pathways enriched by these genes, which may lead to a certain extent confusion for subsequent explanation. Thus, we could use this R script to narrow down the number of significant pathways by using the Jaccard distance and obtain several of imporant cluters for these significant pathways.

**Example**: We give an example for using this R script to cluster the significant pathways enriched by Sherlock-identified genes for coronary artery disease. 

