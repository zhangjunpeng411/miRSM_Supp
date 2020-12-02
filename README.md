# The supplementary material of miRSM package
In the supplementary material, we provide a case study to infer and analyze miRNA sponge modules by using miRSM package.

## Background
MicroRNAs (miRNAs) play key roles in many biological processes including cancers. Thus, uncovering miRNA functions and regulatory mechanisms is important for gene diagnosis and therapy. Previous studies have shown that a pool of coding and non-coding RNAs that shares common miRNA biding sites competes with each other, thus alter miRNA activity. The corresponding regulatory mechanism is named competing endogenous RNA (ceRNA) hypothesis. These RNAs are called ceRNAs or miRNA sponges or miRNA decoys, and include long non-coding RNAs (lncRNAs), pseudogenes, circular RNAs (circRNAs) and messenger RNAs (mRNAs), etc. To study the module-level properties of miRNA sponges, it is necessary to identify miRNA sponge modules. The miRNA sponge modules will help to reveal the biological mechanism in cancer.

To speed up the research of miRNA sponge modules, we develop an R/Bioconductor package ‘miRSM’ to infer miRNA sponge modules. Unlike the existing R/Bioconductor packages (‘miRspongeR’ and ‘SPONGE’), ‘miRSM’ focuses on identifying miRNA sponge modules by integrating expression data and miRNA-target binding information instead of miRNA sponge interaction networks. 

## Description of each file
BRCA_SampleData.rda: Matched miRNA, lncRNA and mRNA expression data of BRCA, putative miRNA-target interactions, BRCA related genes, and experimentally validated lncRNA related miRNA sponge interactions.

miRSM_scripts.R: Scripts for a case study to infer and analyze miRNA sponge modules by using miRSM package.

## The usage
Paste all files into a single folder (set the folder as the directory of R environment). Users can simply run the scripts as follows. The usage of functions in miRSM package can be seen at https://bioconductor.org/packages/miRSM/.

```{r echo=FALSE, results='hide', message=FALSE}
source("miRSM_scripts.R")
```

