## Load data including miRNA, lncRNA, mRNA expression data, and miRNA-target binding information
load("BRCA_SampleData.rda")

## Load miRSM package, and reimport NMF package to avoid conflicts with DelayedArray package
suppressMessages(library(miRSM))
suppressMessages(library(NMF))

## Identification of lncRNA-mRNA co-expression modules using different methods
# The number of gene co-expression modules
num.modules <- 72
set.seed(12345)

# WGCNA method in WGCNA package
modulegenes_WGCNA <- module_WGCNA(ceRExp, mRExp)

# GFA method in GFA package
modulegenes_GFA <- module_GFA(ceRExp, mRExp)

# Method in igraph package, and we select greedy
modulegenes_igraph_greedy <- module_igraph(ceRExp, mRExp, cor.method = "pearson", 
                                           pos.p.value.cutoff = 1e-10, 
				           cluster.method = "greedy")

# Method in ProNet package, and we select MCL
modulegenes_ProNet_MCL <- module_ProNet(ceRExp, mRExp, cor.method = "pearson", 
                                        pos.p.value.cutoff = 1e-10, 
				        cluster.method = "MCL")

# Method in NMF package, and default method is brunet
modulegenes_NMF_brunet <- module_NMF(ceRExp, mRExp, num.modules = num.modules)

# Method in stats package, and we select kmeans. 
modulegenes_clust_kmeans <- module_clust(ceRExp, mRExp, cluster.method = "kmeans", num.modules = num.modules)

# Method in fabia package, and we select fabia.
modulegenes_biclust_fabia <- module_biclust(ceRExp, mRExp, BCmethod = "fabia", num.modules = num.modules)

## Identification of lncRNA related miRNA sponge modules using Sensitivity Canonical Correlation (SCC), 
## Sensitivity Distance Correlation (SDC) and Sensitivity RV Coefficient (SRVC) methods.

## WGCNA 
miRSM_SCC_WGCNA <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_WGCNA, 
                       num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                       method = "SCC",  MC.cutoff = 0.8, SMC.cutoff = 0.1)

miRSM_SDC_WGCNA <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_WGCNA, 
                       num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                       method = "SDC",  MC.cutoff = 0.8, SMC.cutoff = 0.1)

miRSM_SRVC_WGCNA <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_WGCNA, 
                        num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                        method = "SRVC",  MC.cutoff = 0.8, SMC.cutoff = 0.1,
                        RV_method = "RV")

## GFA 
miRSM_SCC_GFA <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_GFA, 
                       num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                       method = "SCC",  MC.cutoff = 0.8, SMC.cutoff = 0.1)

miRSM_SDC_GFA <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_GFA, 
                       num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                       method = "SDC",  MC.cutoff = 0.8, SMC.cutoff = 0.1)

miRSM_SRVC_GFA <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_GFA, 
                        num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                        method = "SRVC",  MC.cutoff = 0.8, SMC.cutoff = 0.1,
                        RV_method = "RV")

## greedy 
miRSM_SCC_greedy <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_igraph_greedy, 
                       num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                       method = "SCC",  MC.cutoff = 0.8, SMC.cutoff = 0.1)

miRSM_SDC_greedy <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_igraph_greedy, 
                       num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                       method = "SDC",  MC.cutoff = 0.8, SMC.cutoff = 0.1)

miRSM_SRVC_greedy <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_igraph_greedy, 
                        num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                        method = "SRVC",  MC.cutoff = 0.8, SMC.cutoff = 0.1,
                        RV_method = "RV")

## MCL 
miRSM_SCC_MCL <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_ProNet_MCL, 
                       num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                       method = "SCC",  MC.cutoff = 0.8, SMC.cutoff = 0.1)

miRSM_SDC_MCL <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_ProNet_MCL, 
                       num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                       method = "SDC",  MC.cutoff = 0.8, SMC.cutoff = 0.1)

miRSM_SRVC_MCL <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_ProNet_MCL, 
                        num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                        method = "SRVC",  MC.cutoff = 0.8, SMC.cutoff = 0.1,
                        RV_method = "RV")

## NMF 
miRSM_SCC_NMF <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_NMF_brunet, 
                       num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                       method = "SCC",  MC.cutoff = 0.8, SMC.cutoff = 0.1)

miRSM_SDC_NMF <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_NMF_brunet, 
                       num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                       method = "SDC",  MC.cutoff = 0.8, SMC.cutoff = 0.1)

miRSM_SRVC_NMF <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_NMF_brunet, 
                        num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                        method = "SRVC",  MC.cutoff = 0.8, SMC.cutoff = 0.1,
                        RV_method = "RV")

## kmeans 
miRSM_SCC_kmeans <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_clust_kmeans, 
                       num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                       method = "SCC",  MC.cutoff = 0.8, SMC.cutoff = 0.1)

miRSM_SDC_kmeans <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_clust_kmeans, 
                       num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                       method = "SDC",  MC.cutoff = 0.8, SMC.cutoff = 0.1)

miRSM_SRVC_kmeans <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_clust_kmeans, 
                        num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                        method = "SRVC",  MC.cutoff = 0.8, SMC.cutoff = 0.1,
                        RV_method = "RV")

## fabia 
miRSM_SCC_fabia <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_biclust_fabia, 
                       num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                       method = "SCC",  MC.cutoff = 0.8, SMC.cutoff = 0.1)

miRSM_SDC_fabia <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_biclust_fabia, 
                       num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                       method = "SDC",  MC.cutoff = 0.8, SMC.cutoff = 0.1)

miRSM_SRVC_fabia <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_biclust_fabia, 
                        num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                        method = "SRVC",  MC.cutoff = 0.8, SMC.cutoff = 0.1,
                        RV_method = "RV")

## Modular analysis of the identified miRNA sponge modules. 
# Functional analysis of miRNA sponge modules
miRSM_SCC_fabia_FEA <- module_FA(miRSM_SCC_fabia[[2]], Analysis.type = "FEA")
miRSM_SCC_fabia_DEA <- module_FA(miRSM_SCC_fabia[[2]], Analysis.type = "DEA")
miRSM_SDC_fabia_FEA <- module_FA(miRSM_SDC_fabia[[2]], Analysis.type = "FEA")
miRSM_SDC_fabia_DEA <- module_FA(miRSM_SDC_fabia[[2]], Analysis.type = "DEA")
miRSM_SRVC_fabia_FEA <- module_FA(miRSM_SRVC_fabia[[2]], Analysis.type = "FEA")
miRSM_SRVC_fabia_DEA <- module_FA(miRSM_SRVC_fabia[[2]], Analysis.type = "DEA")

# BRCA enrichment analysis of miRNA sponge modules
miRSM_SCC_fabia_pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, miRSM_SCC_fabia[[2]])
miRSM_SDC_fabia_pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, miRSM_SDC_fabia[[2]])
miRSM_SRVC_fabia_pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, miRSM_SRVC_fabia[[2]])

# Validation of miRNA sponge interactions in miRNA sponge modules
miRSM_SCC_fabia_Validate <- module_Validate(miRSM_SCC_fabia[[2]], Groundtruth)
miRSM_SDC_fabia_Validate <- module_Validate(miRSM_SDC_fabia[[2]], Groundtruth)
miRSM_SRVC_fabia_Validate <- module_Validate(miRSM_SRVC_fabia[[2]], Groundtruth)

# Co-expression analysis of miRNA sponge modules
miRSM_SCC_fabia_Coexpress <-  module_Coexpress(ceRExp, mRExp, miRSM_SCC_fabia[[2]], resample = 1000, method = "mean", test.method = "t.test")
miRSM_SDC_fabia_Coexpress <-  module_Coexpress(ceRExp, mRExp, miRSM_SDC_fabia[[2]], resample = 1000, method = "mean", test.method = "t.test")
miRSM_SRVC_fabia_Coexpress <-  module_Coexpress(ceRExp, mRExp, miRSM_SRVC_fabia[[2]], resample = 1000, method = "mean", test.method = "t.test")

# miRNA distribution analysis of sharing miRNAs
miRSM_SCC_fabia_share_miRs <-  share_miRs(miRExp, ceRExp, mRExp, miRTarget, miRSM_SCC_fabia[[2]])
miRSM_SCC_fabia_miRdistribute <- module_miRdistribute(miRSM_SCC_fabia_share_miRs)
miRSM_SDC_fabia_share_miRs <-  share_miRs(miRExp, ceRExp, mRExp, miRTarget, miRSM_SDC_fabia[[2]])
miRSM_SDC_fabia_miRdistribute <- module_miRdistribute(miRSM_SDC_fabia_share_miRs)
miRSM_SRVC_fabia_share_miRs <-  share_miRs(miRExp, ceRExp, mRExp, miRTarget, miRSM_SRVC_fabia[[2]])
miRSM_SRVC_fabia_miRdistribute <- module_miRdistribute(miRSM_SRVC_fabia_share_miRs)

# Predict miRNA-target interactions
miRSM_SCC_fabia_miRtarget <- module_miRtarget(miRSM_SCC_fabia_share_miRs, miRSM_SCC_fabia[[2]])
miRSM_SDC_fabia_miRtarget <- module_miRtarget(miRSM_SDC_fabia_share_miRs, miRSM_SCC_fabia[[2]])
miRSM_SRVC_fabia_miRtarget <- module_miRtarget(miRSM_SRVC_fabia_share_miRs, miRSM_SCC_fabia[[2]])

# Identify miRNA sponge interactions
miRSM_SCC_fabia_miRsponge <- module_miRsponge(ceRExp, mRExp, miRSM_SCC_fabia[[2]])
miRSM_SDC_fabia_miRsponge <- module_miRsponge(ceRExp, mRExp, miRSM_SDC_fabia[[2]])
miRSM_SRVC_fabia_miRsponge <- module_miRsponge(ceRExp, mRExp, miRSM_SRVC_fabia[[2]])

save.image("TCGA_BRCA_miRSM.RData")
