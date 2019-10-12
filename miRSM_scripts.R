## Load data including miRNA, lncRNA, mRNA expression data, and miRNA-target binding information
load("BRCA_SampleData.rda")

## Load miRSM package, and reimport NMF package to avoid conflicts with DelayedArray package
suppressMessages(library(miRSM))
suppressMessages(library(miRspongeR))
suppressMessages(library(NMF))

## Identification of lncRNA-mRNA co-expression modules using different methods
# The number of gene co-expression modules
num.modules <- 72
set.seed(12345)

# WGCNA method in WGCNA package
modulegenes_WGCNA <- module_WGCNA(ceRExp, mRExp)

# GFA method in GFA package
modulegenes_GFA <- module_GFA(ceRExp, mRExp)

# Methods in igraph package, and default method is greedy
modulegenes_igraph_greedy <- module_igraph(ceRExp, mRExp, cor.method = "pearson", 
                                           pos.p.value.cutoff = 1e-10, 
				           cluster.method = "greedy")

# Methods in ProNet package, and default method is MCL
modulegenes_ProNet_MCL <- module_ProNet(ceRExp, mRExp, cor.method = "pearson", 
                                        pos.p.value.cutoff = 1e-10, 
				        cluster.method = "MCL")

# Methods in NMF package, and default method is brunet
modulegenes_NMF_brunet <- module_NMF(ceRExp, mRExp, NMF.algorithm = "brunet", num.modules = num.modules)

# Methods in stats, flashClust, dbscan, subspace, mclust, SOMbrero and ppclust packages, and default method is kmeans. 
# We select kmeans, hclust, dbscan, clique, gmm, som and fcm.
modulegenes_clust_kmeans <- module_clust(ceRExp, mRExp, cluster.method = "kmeans", num.modules = num.modules)
modulegenes_clust_hclust <- module_clust(ceRExp, mRExp, cluster.method = "hclust", num.modules = num.modules)
modulegenes_clust_dbscan <- module_clust(ceRExp, mRExp, cluster.method = "dbscan")
modulegenes_clust_clique <- module_clust(ceRExp, mRExp, cluster.method = "clique", num.modules = num.modules)
modulegenes_clust_gmm <- module_clust(ceRExp, mRExp, cluster.method = "gmm", num.modules = num.modules)
modulegenes_clust_som <- module_clust(ceRExp, mRExp, cluster.method = "som")
modulegenes_clust_fcm <- module_clust(ceRExp, mRExp, cluster.method = "fcm", num.modules = num.modules)

# Methods in biclust, runibic, iBBiG, fabia, BicARE, isa2, s4vd, BiBitR and rqubic packages, and default method is fabia. 
# We select BCPlaid, BCUnibic, iBBiG, fabia, FLOC, isa, BCs4vd, bibit and quBicluster.
modulegenes_biclust_BCSpectral <- module_biclust(ceRExp, mRExp, BCmethod = "BCSpectral")
modulegenes_biclust_BCUnibic <- module_biclust(ceRExp, mRExp, BCmethod = "BCUnibic")
modulegenes_biclust_iBBiG <- module_biclust(ceRExp, mRExp, BCmethod = "iBBiG", num.modules = num.modules)
modulegenes_biclust_fabia <- module_biclust(ceRExp, mRExp, BCmethod = "fabia", num.modules = num.modules)
modulegenes_biclust_FLOC <- module_biclust(ceRExp, mRExp, BCmethod = "FLOC", num.modules = num.modules)
modulegenes_biclust_isa <- module_biclust(ceRExp, mRExp, BCmethod = "isa")
modulegenes_biclust_BCs4vd <- module_biclust(ceRExp, mRExp, BCmethod = "BCs4vd", num.modules = num.modules)
modulegenes_biclust_bibit <- module_biclust(ceRExp, mRExp, BCmethod = "bibit", num.modules = num.modules)
modulegenes_biclust_quBicluster <- module_biclust(ceRExp, mRExp, BCmethod = "quBicluster", num.modules = num.modules)

## Identification of lncRNA related miRNA sponge modules by integrating the canonical correlation (CC) and sensitivity canonical correlation (SCC) methods.
## As a result, 17 out 21 (bi)cluster methods of identifying lncRNA-mRNA co-expression modules can generate miRNA sponge modules.
miRSM_WGCNA <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_WGCNA, 
                        num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                        method = "CCplusSCC",  CC.cutoff = 0.8, SCC.cutoff = 0.1)
miRSM_GFA <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_GFA, 
                        num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                        method = "CCplusSCC",  CC.cutoff = 0.8, SCC.cutoff = 0.1)
miRSM_igraph_greedy <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_igraph_greedy, 
                        num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                        method = "CCplusSCC",  CC.cutoff = 0.8, SCC.cutoff = 0.1)
miRSM_ProNet_MCL <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_ProNet_MCL, 
                        num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                        method = "CCplusSCC",  CC.cutoff = 0.8, SCC.cutoff = 0.1)
miRSM_NMF_brunet <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_NMF_brunet, 
                        num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                        method = "CCplusSCC",  CC.cutoff = 0.8, SCC.cutoff = 0.1)
miRSM_clust_kmeans <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_clust_kmeans, 
                         num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                         method = "CCplusSCC",  CC.cutoff = 0.8, SCC.cutoff = 0.1)
miRSM_clust_hclust <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_clust_hclust, 
                         num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                         method = "CCplusSCC",  CC.cutoff = 0.8, SCC.cutoff = 0.1)
miRSM_clust_dbscan <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_clust_dbscan, 
                         num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                         method = "CCplusSCC",  CC.cutoff = 0.8, SCC.cutoff = 0.1)
miRSM_clust_clique <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_clust_clique, 
                         num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                         method = "CCplusSCC",  CC.cutoff = 0.8, SCC.cutoff = 0.1)
miRSM_clust_gmm <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_clust_gmm, 
                         num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                         method = "CCplusSCC",  CC.cutoff = 0.8, SCC.cutoff = 0.1)
miRSM_clust_som <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_clust_som, 
                         num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                         method = "CCplusSCC",  CC.cutoff = 0.8, SCC.cutoff = 0.1)
miRSM_clust_fcm <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_clust_fcm, 
                         num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                         method = "CCplusSCC",  CC.cutoff = 0.8, SCC.cutoff = 0.1)
miRSM_biclust_BCSpectral <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_biclust_BCSpectral, 
                         num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                         method = "CCplusSCC",  CC.cutoff = 0.8, SCC.cutoff = 0.1)
miRSM_biclust_BCUnibic <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_biclust_BCUnibic, 
                        num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                        method = "CCplusSCC",  CC.cutoff = 0.8, SCC.cutoff = 0.1)
miRSM_biclust_iBBiG <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_biclust_iBBiG, 
                        num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                        method = "CCplusSCC",  CC.cutoff = 0.8, SCC.cutoff = 0.1)
miRSM_biclust_fabia <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_biclust_fabia, 
                        num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                        method = "CCplusSCC",  CC.cutoff = 0.8, SCC.cutoff = 0.1)
miRSM_biclust_FLOC <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_biclust_FLOC, 
                        num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                        method = "CCplusSCC",  CC.cutoff = 0.8, SCC.cutoff = 0.1)
miRSM_biclust_isa <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_biclust_isa, 
                        num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                        method = "CCplusSCC",  CC.cutoff = 0.8, SCC.cutoff = 0.1)
miRSM_biclust_BCs4vd <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_biclust_BCs4vd, 
                        num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                        method = "CCplusSCC",  CC.cutoff = 0.8, SCC.cutoff = 0.1)
miRSM_biclust_bibit <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_biclust_bibit, 
                        num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                        method = "CCplusSCC",  CC.cutoff = 0.8, SCC.cutoff = 0.1)
miRSM_biclust_quBicluster <- miRSM(miRExp, ceRExp, mRExp, miRTarget, modulegenes_biclust_quBicluster, 
                        num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                        method = "CCplusSCC",  CC.cutoff = 0.8, SCC.cutoff = 0.1)

## Modular analysis of the identified miRNA sponge modules. 
# Functional analysis of miRNA sponge modules
miRSM_WGCNA_FEA <- module_FA(miRSM_WGCNA[[2]], Analysis.type = "FEA")
miRSM_WGCNA_DEA <- module_FA(miRSM_WGCNA[[2]], Analysis.type = "DEA")
miRSM_GFA_FEA <- module_FA(miRSM_GFA[[2]], Analysis.type = "FEA")
miRSM_GFA_DEA <- module_FA(miRSM_GFA[[2]], Analysis.type = "DEA")
miRSM_igraph_greedy_FEA <- module_FA(miRSM_igraph_greedy[[2]], Analysis.type = "FEA")
miRSM_igraph_greedy_DEA <- module_FA(miRSM_igraph_greedy[[2]], Analysis.type = "DEA")
miRSM_ProNet_MCL_FEA <- module_FA(miRSM_ProNet_MCL[[2]], Analysis.type = "FEA")
miRSM_ProNet_MCL_DEA <- module_FA(miRSM_ProNet_MCL[[2]], Analysis.type = "DEA")
miRSM_NMF_brunet_FEA <- module_FA(miRSM_NMF_brunet[[2]], Analysis.type = "FEA")
miRSM_NMF_brunet_DEA <- module_FA(miRSM_NMF_brunet[[2]], Analysis.type = "DEA")
miRSM_clust_kmeans_FEA <- module_FA(miRSM_clust_kmeans[[2]], Analysis.type = "FEA")
miRSM_clust_kmeans_DEA <- module_FA(miRSM_clust_kmeans[[2]], Analysis.type = "DEA")
miRSM_clust_hclust_FEA <- module_FA(miRSM_clust_hclust[[2]], Analysis.type = "FEA")
miRSM_clust_hclust_DEA <- module_FA(miRSM_clust_hclust[[2]], Analysis.type = "DEA")
miRSM_clust_gmm_FEA <- module_FA(miRSM_clust_gmm[[2]], Analysis.type = "FEA")
miRSM_clust_gmm_DEA <- module_FA(miRSM_clust_gmm[[2]], Analysis.type = "DEA")
miRSM_clust_som_FEA <- module_FA(miRSM_clust_som[[2]], Analysis.type = "FEA")
miRSM_clust_som_DEA <- module_FA(miRSM_clust_som[[2]], Analysis.type = "DEA")
miRSM_clust_fcm_FEA <- module_FA(miRSM_clust_fcm[[2]], Analysis.type = "FEA")
miRSM_clust_fcm_DEA <- module_FA(miRSM_clust_fcm[[2]], Analysis.type = "DEA")
miRSM_biclust_BCSpectral_FEA <- module_FA(miRSM_biclust_BCSpectral[[2]], Analysis.type = "FEA")
miRSM_biclust_BCSpectral_DEA <- module_FA(miRSM_biclust_BCSpectral[[2]], Analysis.type = "DEA")
miRSM_biclust_iBBiG_FEA <- module_FA(miRSM_biclust_iBBiG[[2]], Analysis.type = "FEA")
miRSM_biclust_iBBiG_DEA <- module_FA(miRSM_biclust_iBBiG[[2]], Analysis.type = "DEA")
miRSM_biclust_fabia_FEA <- module_FA(miRSM_biclust_fabia[[2]], Analysis.type = "FEA")
miRSM_biclust_fabia_DEA <- module_FA(miRSM_biclust_fabia[[2]], Analysis.type = "DEA")
miRSM_biclust_FLOC_FEA <- module_FA(miRSM_biclust_FLOC[[2]], Analysis.type = "FEA")
miRSM_biclust_FLOC_DEA <- module_FA(miRSM_biclust_FLOC[[2]], Analysis.type = "DEA")
miRSM_biclust_isa_FEA <- module_FA(miRSM_biclust_isa[[2]], Analysis.type = "FEA")
miRSM_biclust_isa_DEA <- module_FA(miRSM_biclust_isa[[2]], Analysis.type = "DEA")
miRSM_biclust_bibit_FEA <- module_FA(miRSM_biclust_bibit[[2]], Analysis.type = "FEA")
miRSM_biclust_bibit_DEA <- module_FA(miRSM_biclust_bibit[[2]], Analysis.type = "DEA")
miRSM_biclust_quBicluster_FEA <- module_FA(miRSM_biclust_quBicluster[[2]], Analysis.type = "FEA")
miRSM_biclust_quBicluster_DEA <- module_FA(miRSM_biclust_quBicluster[[2]], Analysis.type = "DEA")

# BRCA enrichment analysis of miRNA sponge modules
miRSM_WGCNA_pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, miRSM_WGCNA[[2]])
miRSM_GFA_pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, miRSM_GFA[[2]])
miRSM_igraph_greedy_pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, miRSM_igraph_greedy[[2]])
miRSM_ProNet_MCL_pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, miRSM_ProNet_MCL[[2]])
miRSM_NMF_brunet_pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, miRSM_NMF_brunet[[2]])
miRSM_clust_kmeans_pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, miRSM_clust_kmeans[[2]])
miRSM_clust_hclust_pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, miRSM_clust_hclust[[2]])
miRSM_clust_gmm_pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, miRSM_clust_gmm[[2]])
miRSM_clust_som_pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, miRSM_clust_som[[2]])
miRSM_clust_fcm_pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, miRSM_clust_fcm[[2]])
miRSM_biclust_BCSpectral_pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, miRSM_biclust_BCSpectral[[2]])
miRSM_biclust_iBBiG_pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, miRSM_biclust_iBBiG[[2]])
miRSM_biclust_fabia_pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, miRSM_biclust_fabia[[2]])
miRSM_biclust_FLOC_pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, miRSM_biclust_FLOC[[2]])
miRSM_biclust_isa_pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, miRSM_biclust_isa[[2]])
miRSM_biclust_bibit_pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, miRSM_biclust_bibit[[2]])
miRSM_biclust_quBicluster_pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, miRSM_biclust_quBicluster[[2]])

# Validation of miRNA sponge interactions in miRNA sponge modules
miRSM_WGCNA_Validate <- module_Validate(miRSM_WGCNA[[2]], Groundtruth)
miRSM_GFA_Validate <- module_Validate(miRSM_GFA[[2]], Groundtruth)
miRSM_igraph_greedy_Validate <- module_Validate(miRSM_igraph_greedy[[2]], Groundtruth)
miRSM_ProNet_MCL_Validate <- module_Validate(miRSM_ProNet_MCL[[2]], Groundtruth)
miRSM_NMF_brunet_Validate <- module_Validate(miRSM_NMF_brunet[[2]], Groundtruth)
miRSM_clust_kmeans_Validate <- module_Validate(miRSM_clust_kmeans[[2]], Groundtruth)
miRSM_clust_hclust_Validate <- module_Validate(miRSM_clust_hclust[[2]], Groundtruth)
miRSM_clust_gmm_Validate <- module_Validate(miRSM_clust_gmm[[2]], Groundtruth)
miRSM_clust_som_Validate <- module_Validate(miRSM_clust_som[[2]], Groundtruth)
miRSM_clust_fcm_Validate <- module_Validate(miRSM_clust_fcm[[2]], Groundtruth)
miRSM_biclust_BCSpectral_Validate <- module_Validate(miRSM_biclust_BCSpectral[[2]], Groundtruth)
miRSM_biclust_iBBiG_Validate <- module_Validate(miRSM_biclust_iBBiG[[2]], Groundtruth)
miRSM_biclust_fabia_Validate <- module_Validate(miRSM_biclust_fabia[[2]], Groundtruth)
miRSM_biclust_FLOC_Validate <- module_Validate(miRSM_biclust_FLOC[[2]], Groundtruth)
miRSM_biclust_isa_Validate <- module_Validate(miRSM_biclust_isa[[2]], Groundtruth)
miRSM_biclust_bibit_Validate <- module_Validate(miRSM_biclust_bibit[[2]], Groundtruth)
miRSM_biclust_quBicluster_Validate <- module_Validate(miRSM_biclust_quBicluster[[2]], Groundtruth)

# Co-expression analysis of miRNA sponge modules
miRSM_WGCNA_Coexpress <-  module_Coexpress(ceRExp, mRExp, miRSM_WGCNA[[2]], resample = 1000, method = "mean")
miRSM_GFA_Coexpress <-  module_Coexpress(ceRExp, mRExp, miRSM_GFA[[2]], resample = 1000, method = "mean")
miRSM_igraph_greedy_Coexpress <-  module_Coexpress(ceRExp, mRExp, miRSM_igraph_greedy[[2]], resample = 1000, method = "mean")
miRSM_ProNet_MCL_Coexpress <-  module_Coexpress(ceRExp, mRExp, miRSM_ProNet_MCL[[2]], resample = 1000, method = "mean")
miRSM_NMF_brunet_Coexpress <-  module_Coexpress(ceRExp, mRExp, miRSM_NMF_brunet[[2]], resample = 1000, method = "mean")
miRSM_clust_kmeans_Coexpress <-  module_Coexpress(ceRExp, mRExp, miRSM_clust_kmeans[[2]], resample = 1000, method = "mean")
miRSM_clust_hclust_Coexpress <-  module_Coexpress(ceRExp, mRExp, miRSM_clust_hclust[[2]], resample = 1000, method = "mean")
miRSM_clust_gmm_Coexpress <-  module_Coexpress(ceRExp, mRExp, miRSM_clust_gmm[[2]], resample = 1000, method = "mean")
miRSM_clust_som_Coexpress <-  module_Coexpress(ceRExp, mRExp, miRSM_clust_som[[2]], resample = 1000, method = "mean")
miRSM_clust_fcm_Coexpress <-  module_Coexpress(ceRExp, mRExp, miRSM_clust_fcm[[2]], resample = 1000, method = "mean")
miRSM_biclust_BCSpectral_Coexpress <-  module_Coexpress(ceRExp, mRExp, miRSM_biclust_BCSpectral[[2]], resample = 1000, method = "mean")
miRSM_biclust_iBBiG_Coexpress <-  module_Coexpress(ceRExp, mRExp, miRSM_biclust_iBBiG[[2]], resample = 1000, method = "mean")
miRSM_biclust_fabia_Coexpress <-  module_Coexpress(ceRExp, mRExp, miRSM_biclust_fabia[[2]], resample = 1000, method = "mean")
miRSM_biclust_FLOC_Coexpress <-  module_Coexpress(ceRExp, mRExp, miRSM_biclust_FLOC[[2]], resample = 1000, method = "mean")
miRSM_biclust_isa_Coexpress <-  module_Coexpress(ceRExp, mRExp, miRSM_biclust_isa[[2]], resample = 1000, method = "mean")
miRSM_biclust_bibit_Coexpress <-  module_Coexpress(ceRExp, mRExp, miRSM_biclust_bibit[[2]], resample = 1000, method = "mean")
miRSM_biclust_quBicluster_Coexpress <-  module_Coexpress(ceRExp, mRExp, miRSM_biclust_quBicluster[[2]], resample = 1000, method = "mean")

# miRNA distribution analysis of sharing miRNAs
miRSM_WGCNA_share_miRs <-  share_miRs(miRExp, ceRExp, mRExp, miRTarget, miRSM_WGCNA[[2]])
miRSM_WGCNA_miRdistribute <- module_miRdistribute(miRSM_WGCNA_share_miRs)
miRSM_GFA_share_miRs <-  share_miRs(miRExp, ceRExp, mRExp, miRTarget, miRSM_GFA[[2]])
miRSM_GFA_miRdistribute <- module_miRdistribute(miRSM_GFA_share_miRs)
miRSM_igraph_greedy_share_miRs <-  share_miRs(miRExp, ceRExp, mRExp, miRTarget, miRSM_igraph_greedy[[2]])
miRSM_igraph_greedy_miRdistribute <- module_miRdistribute(miRSM_igraph_greedy_share_miRs)
miRSM_ProNet_MCL_share_miRs <-  share_miRs(miRExp, ceRExp, mRExp, miRTarget, miRSM_ProNet_MCL[[2]])
miRSM_ProNet_MCL_miRdistribute <- module_miRdistribute(miRSM_ProNet_MCL_share_miRs)
miRSM_NMF_brunet_share_miRs <-  share_miRs(miRExp, ceRExp, mRExp, miRTarget, miRSM_NMF_brunet[[2]])
miRSM_NMF_brunet_miRdistribute <- module_miRdistribute(miRSM_NMF_brunet_share_miRs)
miRSM_clust_kmeans_share_miRs <-  share_miRs(miRExp, ceRExp, mRExp, miRTarget, miRSM_clust_kmeans[[2]])
miRSM_clust_kmeans_miRdistribute <- module_miRdistribute(miRSM_clust_kmeans_share_miRs)
miRSM_clust_hclust_share_miRs <-  share_miRs(miRExp, ceRExp, mRExp, miRTarget, miRSM_clust_hclust[[2]])
miRSM_clust_hclust_miRdistribute <- module_miRdistribute(miRSM_clust_hclust_share_miRs)
miRSM_clust_gmm_share_miRs <-  share_miRs(miRExp, ceRExp, mRExp, miRTarget, miRSM_clust_gmm[[2]])
miRSM_clust_gmm_miRdistribute <- module_miRdistribute(miRSM_clust_gmm_share_miRs)
miRSM_clust_som_share_miRs <-  share_miRs(miRExp, ceRExp, mRExp, miRTarget, miRSM_clust_som[[2]])
miRSM_clust_som_miRdistribute <- module_miRdistribute(miRSM_clust_som_share_miRs)
miRSM_clust_fcm_share_miRs <-  share_miRs(miRExp, ceRExp, mRExp, miRTarget, miRSM_clust_fcm[[2]])
miRSM_clust_fcm_miRdistribute <- module_miRdistribute(miRSM_clust_fcm_share_miRs)
miRSM_biclust_BCSpectral_share_miRs <-  share_miRs(miRExp, ceRExp, mRExp, miRTarget, miRSM_biclust_BCSpectral[[2]])
miRSM_biclust_BCSpectral_miRdistribute <- module_miRdistribute(miRSM_biclust_BCSpectral_share_miRs)
miRSM_biclust_iBBiG_share_miRs <-  share_miRs(miRExp, ceRExp, mRExp, miRTarget, miRSM_biclust_iBBiG[[2]])
miRSM_biclust_iBBiG_miRdistribute <- module_miRdistribute(miRSM_biclust_iBBiG_share_miRs)
miRSM_biclust_fabia_share_miRs <-  share_miRs(miRExp, ceRExp, mRExp, miRTarget, miRSM_biclust_fabia[[2]])
miRSM_biclust_fabia_miRdistribute <- module_miRdistribute(miRSM_biclust_fabia_share_miRs)
miRSM_biclust_FLOC_share_miRs <-  share_miRs(miRExp, ceRExp, mRExp, miRTarget, miRSM_biclust_FLOC[[2]])
miRSM_biclust_FLOC_miRdistribute <- module_miRdistribute(miRSM_biclust_FLOC_share_miRs)
miRSM_biclust_isa_share_miRs <-  share_miRs(miRExp, ceRExp, mRExp, miRTarget, miRSM_biclust_isa[[2]])
miRSM_biclust_isa_miRdistribute <- module_miRdistribute(miRSM_biclust_isa_share_miRs)
miRSM_biclust_bibit_share_miRs <-  share_miRs(miRExp, ceRExp, mRExp, miRTarget, miRSM_biclust_bibit[[2]])
miRSM_biclust_bibit_miRdistribute <- module_miRdistribute(miRSM_biclust_bibit_share_miRs)
miRSM_biclust_quBicluster_share_miRs <-  share_miRs(miRExp, ceRExp, mRExp, miRTarget, miRSM_biclust_quBicluster[[2]])
miRSM_biclust_quBicluster_miRdistribute <- module_miRdistribute(miRSM_biclust_quBicluster_share_miRs)

# Predict miRNA-target interactions
miRSM_WGCNA_miRtarget <- module_miRtarget(miRSM_WGCNA_share_miRs, miRSM_WGCNA[[2]])
miRSM_GFA_miRtarget <- module_miRtarget(miRSM_GFA_share_miRs, miRSM_GFA[[2]])
miRSM_igraph_greedy_miRtarget <- module_miRtarget(miRSM_igraph_greedy_share_miRs, miRSM_igraph_greedy[[2]])
miRSM_ProNet_MCL_miRtarget <- module_miRtarget(miRSM_ProNet_MCL_share_miRs, miRSM_ProNet_MCL[[2]])
miRSM_NMF_brunet_miRtarget <- module_miRtarget(miRSM_NMF_brunet_share_miRs, miRSM_NMF_brunet[[2]])
miRSM_clust_kmeans_miRtarget <- module_miRtarget(miRSM_clust_kmeans_share_miRs, miRSM_clust_kmeans[[2]])
miRSM_clust_hclust_miRtarget <- module_miRtarget(miRSM_clust_hclust_share_miRs, miRSM_clust_hclust[[2]])
miRSM_clust_gmm_miRtarget <- module_miRtarget(miRSM_clust_gmm_share_miRs, miRSM_clust_gmm[[2]])
miRSM_clust_som_miRtarget <- module_miRtarget(miRSM_clust_som_share_miRs, miRSM_clust_som[[2]])
miRSM_clust_fcm_miRtarget <- module_miRtarget(miRSM_clust_fcm_share_miRs, miRSM_clust_fcm[[2]])
miRSM_biclust_BCSpectral_miRtarget <- module_miRtarget(miRSM_biclust_BCSpectral_share_miRs, miRSM_biclust_BCSpectral[[2]])
miRSM_biclust_iBBiG_miRtarget <- module_miRtarget(miRSM_biclust_iBBiG_share_miRs, miRSM_biclust_iBBiG[[2]])
miRSM_biclust_fabia_miRtarget <- module_miRtarget(miRSM_biclust_fabia_share_miRs, miRSM_biclust_fabia[[2]])
miRSM_biclust_FLOC_miRtarget <- module_miRtarget(miRSM_biclust_FLOC_share_miRs, miRSM_biclust_FLOC[[2]])
miRSM_biclust_isa_miRtarget <- module_miRtarget(miRSM_biclust_isa_share_miRs, miRSM_biclust_isa[[2]])
miRSM_biclust_bibit_miRtarget <- module_miRtarget(miRSM_biclust_bibit_share_miRs, miRSM_biclust_bibit[[2]])
miRSM_biclust_quBicluster_miRtarget <- module_miRtarget(miRSM_biclust_quBicluster_share_miRs, miRSM_biclust_quBicluster[[2]])

# Identify miRNA sponge interactions
miRSM_WGCNA_miRsponge <- module_miRsponge(ceRExp, mRExp, miRSM_WGCNA[[2]])
miRSM_GFA_miRsponge <- module_miRsponge(ceRExp, mRExp, miRSM_GFA[[2]])
miRSM_igraph_greedy_miRsponge <- module_miRsponge(ceRExp, mRExp, miRSM_igraph_greedy[[2]])
miRSM_ProNet_MCL_miRsponge <- module_miRsponge(ceRExp, mRExp, miRSM_ProNet_MCL[[2]])
miRSM_NMF_brunet_miRsponge <- module_miRsponge(ceRExp, mRExp, miRSM_NMF_brunet[[2]])
miRSM_clust_kmeans_miRsponge <- module_miRsponge(ceRExp, mRExp, miRSM_clust_kmeans[[2]])
miRSM_clust_hclust_miRsponge <- module_miRsponge(ceRExp, mRExp, miRSM_clust_hclust[[2]])
miRSM_clust_gmm_miRsponge <- module_miRsponge(ceRExp, mRExp, miRSM_clust_gmm[[2]])
miRSM_clust_som_miRsponge <- module_miRsponge(ceRExp, mRExp, miRSM_clust_som[[2]])
miRSM_clust_fcm_miRsponge <- module_miRsponge(ceRExp, mRExp, miRSM_clust_fcm[[2]])
miRSM_biclust_BCSpectral_miRsponge <- module_miRsponge(ceRExp, mRExp, miRSM_biclust_BCSpectral[[2]])
miRSM_biclust_iBBiG_miRsponge <- module_miRsponge(ceRExp, mRExp, miRSM_biclust_iBBiG[[2]])
miRSM_biclust_fabia_miRsponge <- module_miRsponge(ceRExp, mRExp, miRSM_biclust_fabia[[2]])
miRSM_biclust_FLOC_miRsponge <- module_miRsponge(ceRExp, mRExp, miRSM_biclust_FLOC[[2]])
miRSM_biclust_isa_miRsponge <- module_miRsponge(ceRExp, mRExp, miRSM_biclust_isa[[2]])
miRSM_biclust_bibit_miRsponge <- module_miRsponge(ceRExp, mRExp, miRSM_biclust_bibit[[2]])
miRSM_biclust_quBicluster_miRsponge <- module_miRsponge(ceRExp, mRExp, miRSM_biclust_quBicluster[[2]])

save.image("TCGA_BRCA_miRSM.RData")
