##########################################################################################################################
######################################## Utility functions for the case study ############################################
##########################################################################################################################
## Function of SVM k-fold classification
SVMClassifier <- function(ExpData, tumor.index, normal.index, nfolds = 10){
    
    set.seed(123)
    ExpData1 <- ExpData[tumor.index, ]
    ExpData2 <- ExpData[normal.index, ]

    #Randomly shuffle the data
    ExpData1 <- ExpData1[sample(nrow(ExpData1)), ]
    ExpData2 <- ExpData2[sample(nrow(ExpData2)), ]
    ExpData <- rbind(ExpData1, ExpData2)

    #Create nfolds equally size folds
    folds1 <- cut(seq(1, nrow(ExpData1)), breaks = nfolds, labels = FALSE)
    folds2 <- cut(seq(1, nrow(ExpData2)), breaks = nfolds, labels = FALSE)

    #Perform nfolds fold cross validation
    testIndexes <- lapply(seq(nfolds), function(i) c(which(folds1 == i, arr.ind = TRUE), which(folds2 == i, arr.ind = TRUE) + length(tumor.index)))
    testData <- lapply(seq(nfolds), function(i) ExpData[testIndexes[[i]], ])
    trainData <- lapply(seq(nfolds), function(i) ExpData[-testIndexes[[i]], ])

    #learning from training
    svmmodel <- lapply(seq(nfolds), function(i) svm(Class~., data = trainData[[i]], probability=TRUE))

    #predicting the test data
    svmmodel.predict <- lapply(seq(nfolds), function(i) predict(svmmodel[[i]], testData[[i]], decision.values=TRUE))
    svmmodel.probs <- lapply(seq(nfolds), function(i) attr(svmmodel.predict[[i]], "decision.values"))
    svmmodel.labels <- lapply(seq(nfolds), function(i) testData[[i]]$Class)

    #ROC analysis for test data
    svmmodel.prediction <- lapply(seq(nfolds), function(i) prediction(svmmodel.probs[[i]], svmmodel.labels[[i]]))
    # svmmodel.performance <- lapply(seq(nfolds), function(i) performance(svmmodel.prediction[[i]], "tpr", "fpr"))
    svmmodel.auc <- lapply(seq(nfolds), function(i) performance(svmmodel.prediction[[i]],"auc")@y.values[[1]])    
    svmAUC <- mean(unlist(svmmodel.auc))    
    
return(svmAUC)
}

## Identifying miRNA sponge interactions using Sensitivity Correlation (SC) method
# miRExp: miRNA expression data, rows are samples and columns are miRNAs.
# ceRExp: ceRNA expression data, rows are samples and columns are ceRNAs (lncRNAs, pseudogenes, etc.).
# mRExp: miRNA expression data, rows are samples and columns are mRNAs.
# miRTarget: Putative miRNA-target interactions.
# minSharedmiR: Minimum number of sharing miRNAs for rach ceRNA-mRNA pair.
# pvaluecutoff: Cutoff value for p-values.
# poscorcutoff: Cutoff value of positive correlation.
# senscorcutoff: Cutoff value of sensitivity correlation.
SC <- function(miRExp, ceRExp, mRExp, miRTarget, minSharedmiR = 1, 
                 pvaluecutoff = 0.05, poscorcutoff = 0, senscorcutoff = 0.1,
		 num.cores = 6){   
         
        miRTarget <- as.matrix(miRTarget)	
        miRceR <- miRTarget[intersect(which(miRTarget[, 1] %in% colnames(miRExp)),
	                      which(miRTarget[, 2] %in% colnames(ceRExp))), ]
        miRmR <- miRTarget[intersect(which(miRTarget[, 1] %in% colnames(miRExp)),
	                      which(miRTarget[, 2] %in% colnames(mRExp))), ]

	ceRSym <- unique(miRceR[, 2])
        mRSym <- unique(miRmR[, 2])
	miRSym <- unique(c(miRceR[, 1], miRmR[, 1]))

	m <- length(ceRSym)
	n <- length(mRSym)

        ceRExp_query <- ceRExp[, which(colnames(ceRExp) %in% ceRSym)]
        mRExp_query <- mRExp[, which(colnames(mRExp) %in% mRSym)]
	
        Cor.Pvalue <- corAndPvalue(ceRExp_query, mRExp_query)
        
	index <- which(Cor.Pvalue$cor > poscorcutoff & Cor.Pvalue$p < pvaluecutoff, arr.ind = TRUE)        
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)

        Res <- foreach(i = seq_len(nrow(index)), .packages = "corpcor") %dopar% {
	    
	        Interin1 <- miRceR[which(miRceR[, 2] %in% ceRSym[index[i, 1]]), 1]
                Interin2 <- miRmR[which(miRmR[, 2] %in% mRSym[index[i, 2]]), 1]

		M1 <- length(Interin1)
                M2 <- length(Interin2)
                M3 <- length(intersect(Interin1, Interin2))
                M4 <- length(miRSym)
                M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)

                if (M3 >= minSharedmiR & M5 < pvaluecutoff) {
                    
                    C1 <- ceRSym[index[i, 1]]
                    C2 <- mRSym[index[i, 2]]

                    ceRExpIdx <- which(colnames(ceRExp) %in% ceRSym[index[i, 1]])
                    mRExpIdx <- which(colnames(mRExp) %in% mRSym[index[i, 2]])
                    miRExpIdx <- which(colnames(miRExp) %in% intersect(Interin1, Interin2))

		    M6 <- Cor.Pvalue$cor[index[i, 1], index[i, 2]]
		    M7 <- Cor.Pvalue$p[index[i, 1], index[i, 2]]
                    M8 <- M6 - corpcor::pcor.shrink(cbind(ceRExp[, ceRExpIdx], mRExp[, mRExpIdx],
                                                    miRExp[, miRExpIdx]), verbose = FALSE)[1, 2]
                } else {

                C1 <- NA; C2 <- NA; M6 <- NA; M7 <- NA; M8 <- NA 
	       
                }
	       
	        tmp <- c(C1, C2, M3, M5, M6, M7, M8)    
                return(tmp)
	}
        
	# shut down the workers
        stopCluster(cl)
        stopImplicitCluster()
        
	Res <- do.call(rbind, Res)
        Res <- Res[which(as.numeric(Res[, 7]) > senscorcutoff), ]

        colnames(Res) <- c("sponge_1", "sponge_2", "#Shared miRNAs", 
                           "Sig. p.value of sharing miRNAs", 
			   "Correlation", 
                           "Sig. p.value of correlation", 			   
			   "Sensitivity correlation")
			   
       rownames(Res) <- seq_len(nrow(Res))

return(Res)

}

## Identifying miRNA sponge interactions using Sparse Partial correlation ON Gene Expression (SPONGE) method
# miRExp: miRNA expression data, rows are samples and columns are miRNAs.
# ceRExp: ceRNA expression data, rows are samples and columns are ceRNAs (lncRNAs, pseudogenes, etc.).
# mRExp: miRNA expression data, rows are samples and columns are mRNAs.
# miRTarget: Putative miRNA-target interactions.
# coefficient.threshold: Threshold to cross for a regression coefficient to be called significant.
# min.cor: Consider only ceRNA-mRNA pairs with a minimum correlation specified here.
# minSharedmiR: Minimum number of sharing miRNAs for rach ceRNA-mRNA pair.
# pvaluecutoff: Cutoff value for p-values.
# number_of_datasets: The number of datesets defining the precision of the p-value when building null model for p-value computation.
SPONGE <- function(miRExp, ceRExp, mRExp, miRTarget, coefficient.threshold = -0.05, 
                   min.cor = 0.1, minSharedmiR = 1, pvaluecutoff = 0.05, 
		   null_model, num.cores = 6){   
        
	miRTarget <- as.matrix(miRTarget)	
        miRTarget <- miRTarget[intersect(which(miRTarget[, 1] %in% colnames(miRExp)),
	                      which(miRTarget[, 2] %in% colnames(cbind(ceRExp, mRExp)))), ]
        
	miRlist <- unique(miRTarget[, 1])
	Tarlist <- unique(miRTarget[, 2])
        miRTarget_graph <- make_graph(c(t(miRTarget)), directed = FALSE)
	miRTarget_adjacency <- as_adjacency_matrix(miRTarget_graph, sparse = FALSE)
        miRTarget_adjacency_matrix <- miRTarget_adjacency[which(rownames(miRTarget_adjacency) %in% Tarlist), which(colnames(miRTarget_adjacency) %in% miRlist)]

	# get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)
        
	miRTarget_candidates <- sponge_gene_miRNA_interaction_filter(gene_expr = cbind(ceRExp, mRExp),
                                                                     mir_expr = miRExp,
                                                                     mir_predicted_targets = miRTarget_adjacency_matrix,
								     coefficient.threshold = coefficient.threshold)
	
	ceRlist <- names(miRTarget_candidates)[which(names(miRTarget_candidates) %in% colnames(ceRExp))]
	mRlist <- names(miRTarget_candidates)[which(names(miRTarget_candidates) %in% colnames(mRExp))]
	gene.combinations <- expand.grid(ceRlist, mRlist)

        ceRNA_interactions <- sponge(gene_expr = cbind(ceRExp, mRExp),
                                     mir_expr = miRExp,
                                     mir_interactions = miRTarget_candidates,
				     gene.combinations = gene.combinations,
				     min.cor = min.cor)

	# shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

	ceRNA_interactions_filter <- ceRNA_interactions[which(ceRNA_interactions$df >= minSharedmiR), ]
	
	ceRNA_interactions_sign <- sponge_compute_p_values(sponge_result = ceRNA_interactions_filter, 
                                                   null_model = null_model)

	Res <- ceRNA_interactions_sign[which(ceRNA_interactions_sign$p.adj < pvaluecutoff),]

return(Res)

}

## Extracting modules with at least num.ModuleceRs (e.g. 2) ceRNAs and num.ModulemRs (e.g. 2) mRNAs
CandModgenes <- function(ceRExp, mRExp, Modulegenes, num.ModuleceRs = 2, num.ModulemRs = 2){
  
    ceR_Num <- lapply(seq_along(Modulegenes), function(i) length(which(Modulegenes[[i]] %in%
        colnames(ceRExp))))
    mR_Num <- lapply(seq_along(Modulegenes), function(i) length(which(Modulegenes[[i]] %in%
        colnames(mRExp))))

    index <- which(ceR_Num >= num.ModuleceRs & mR_Num >= num.ModulemRs)
    CandidateModulegenes <- lapply(index, function(i) Modulegenes[[i]])
    
return(CandidateModulegenes)
}

##########################################################################################################################
######################################## Using miRSM package for the case study ##########################################
##########################################################################################################################
## Load data including miRNA, lncRNA, mRNA expression data, and miRNA-target binding information
load("BRCA_SampleData.rda")

## Load miRSM package, reimport NMF package to avoid conflicts with DelayedArray package, and some other packages
suppressMessages(library(miRSM))
suppressMessages(library(NMF)))
suppressMessages(library(SPONGE))
suppressMessages(library(igraph))
suppressMessages(library(WGCNA))
suppressMessages(library(corpcor))
suppressMessages(library(miRspongeR))
suppressMessages(library(e1071))
suppressMessages(library(ROCR))
suppressMessages(library(SummarizedExperiment))

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

## Evaluate the significance of each module by using null model
SPONGE_null_model <- sponge_build_null_model(number_of_datasets = 1e+06, number_of_samples = nrow(miRExp), cov_matrices = precomputed_cov_matrices,  ks = seq(0.8, 0.9, 0.1), m_max = 1)
cutoff_modules <- data.frame(geneA = paste("module_1", seq(1)), 
                                 geneB = paste("module_2", seq(1)), 
				 df = replicate(1, 1), 
				 cor = 0.8, 
				 pcor = 0.7, 
				 mscor = 0.1)
cutoff_modules_p.values <- sponge_compute_p_values(sponge_result = cutoff_modules, null_model = SPONGE_null_model)

SPONGE_null_model_1 <- sponge_build_null_model(number_of_datasets = 1e+06, number_of_samples = 100, cov_matrices = precomputed_cov_matrices,  ks = seq(0.8, 0.9, 0.1), m_max = 1)
cutoff_modules <- data.frame(geneA = paste("module_1", seq(1)), 
                                 geneB = paste("module_2", seq(1)), 
				 df = replicate(1, 1), 
				 cor = 0.8, 
				 pcor = 0.7, 
				 mscor = 0.1)
cutoff_modules_p.values_1 <- sponge_compute_p_values(sponge_result = cutoff_modules, null_model = SPONGE_null_model_1)

cutoff_modules <- data.frame(geneA = paste("module_1", seq(1)), 
                                 geneB = paste("module_2", seq(1)), 
				 df = replicate(1, 1), 
				 cor = 0.85, 
				 pcor = 0.7, 
				 mscor = 0.15)
cutoff_modules_p.values_2 <- sponge_compute_p_values(sponge_result = cutoff_modules, null_model = SPONGE_null_model_1)
# WGCNA
miRSM_SCC_WGCNA_modules <- data.frame(geneA = paste("module_1", seq(1)), 
                                 geneB = paste("module_2", seq(1)), 
				 df = replicate(1, 1), 
				 cor = miRSM_SCC_WGCNA[[1]][6], 
				 pcor = miRSM_SCC_WGCNA[[1]][9], 
				 mscor = miRSM_SCC_WGCNA[[1]][10])
miRSM_SCC_WGCNA_modules_p.values <- sponge_compute_p_values(sponge_result = miRSM_SCC_WGCNA_modules, null_model = SPONGE_null_model)

miRSM_SDC_WGCNA_modules <- data.frame(geneA = paste("module_1", seq(1)), 
                                 geneB = paste("module_2", seq(1)), 
				 df = replicate(1, 1), 
				 cor = miRSM_SDC_WGCNA[[1]][6], 
				 pcor = miRSM_SDC_WGCNA[[1]][7], 
				 mscor = miRSM_SDC_WGCNA[[1]][8])
miRSM_SDC_WGCNA_modules_p.values <- sponge_compute_p_values(sponge_result = miRSM_SDC_WGCNA_modules, null_model = SPONGE_null_model)

# GFA
miRSM_SCC_GFA_modules <- data.frame(geneA = paste("module_1", seq(nrow(miRSM_SCC_GFA[[1]]))), 
                                 geneB = paste("module_2", seq(nrow(miRSM_SCC_GFA[[1]]))), 
				 df = replicate(nrow(miRSM_SCC_GFA[[1]]), 1), 
				 cor = miRSM_SCC_GFA[[1]][, 6], 
				 pcor = miRSM_SCC_GFA[[1]][, 9], 
				 mscor = miRSM_SCC_GFA[[1]][, 10])
miRSM_SCC_GFA_modules_p.values <- sponge_compute_p_values(sponge_result = miRSM_SCC_GFA_modules, null_model = SPONGE_null_model)

miRSM_SDC_GFA_modules <- data.frame(geneA = paste("module_1", seq(nrow(miRSM_SDC_GFA[[1]]))), 
                                 geneB = paste("module_2", seq(nrow(miRSM_SDC_GFA[[1]]))), 
				 df = replicate(nrow(miRSM_SDC_GFA[[1]]), 1), 
				 cor = miRSM_SDC_GFA[[1]][, 6], 
				 pcor = miRSM_SDC_GFA[[1]][, 7], 
				 mscor = miRSM_SDC_GFA[[1]][, 8])
miRSM_SDC_GFA_modules_p.values <- sponge_compute_p_values(sponge_result = miRSM_SDC_GFA_modules, null_model = SPONGE_null_model)

miRSM_SRVC_GFA_modules <- data.frame(geneA = paste("module_1", seq(nrow(miRSM_SRVC_GFA[[1]]))), 
                                 geneB = paste("module_2", seq(nrow(miRSM_SRVC_GFA[[1]]))), 
				 df = replicate(nrow(miRSM_SRVC_GFA[[1]]), 1), 
				 cor = miRSM_SRVC_GFA[[1]][, 6], 
				 pcor = miRSM_SRVC_GFA[[1]][, 9], 
				 mscor = miRSM_SRVC_GFA[[1]][, 10])
miRSM_SRVC_GFA_modules_p.values <- sponge_compute_p_values(sponge_result = miRSM_SRVC_GFA_modules, null_model = SPONGE_null_model)

# greedy
miRSM_SCC_greedy_modules <- data.frame(geneA = paste("module_1", seq(nrow(miRSM_SCC_greedy[[1]]))), 
                                 geneB = paste("module_2", seq(nrow(miRSM_SCC_greedy[[1]]))), 
				 df = replicate(nrow(miRSM_SCC_greedy[[1]]), 1), 
				 cor = miRSM_SCC_greedy[[1]][, 6], 
				 pcor = miRSM_SCC_greedy[[1]][, 9], 
				 mscor = miRSM_SCC_greedy[[1]][, 10])
miRSM_SCC_greedy_modules_p.values <- sponge_compute_p_values(sponge_result = miRSM_SCC_greedy_modules, null_model = SPONGE_null_model)

miRSM_SDC_greedy_modules <- data.frame(geneA = paste("module_1", seq(nrow(miRSM_SDC_greedy[[1]]))), 
                                 geneB = paste("module_2", seq(nrow(miRSM_SDC_greedy[[1]]))), 
				 df = replicate(nrow(miRSM_SDC_greedy[[1]]), 1), 
				 cor = miRSM_SDC_greedy[[1]][, 6], 
				 pcor = miRSM_SDC_greedy[[1]][, 7], 
				 mscor = miRSM_SDC_greedy[[1]][, 8])
miRSM_SDC_greedy_modules_p.values <- sponge_compute_p_values(sponge_result = miRSM_SDC_greedy_modules, null_model = SPONGE_null_model)

miRSM_SRVC_greedy_modules <- data.frame(geneA = paste("module_1", seq(1)), 
                                 geneB = paste("module_2", seq(1)), 
				 df = replicate(1, 1), 
				 cor = miRSM_SRVC_greedy[[1]][6], 
				 pcor = miRSM_SRVC_greedy[[1]][9], 
				 mscor = miRSM_SRVC_greedy[[1]][10])
miRSM_SRVC_greedy_modules_p.values <- sponge_compute_p_values(sponge_result = miRSM_SRVC_greedy_modules, null_model = SPONGE_null_model)


# MCL
miRSM_SCC_MCL_modules <- data.frame(geneA = paste("module_1", seq(nrow(miRSM_SCC_MCL[[1]]))), 
                                 geneB = paste("module_2", seq(nrow(miRSM_SCC_MCL[[1]]))), 
				 df = replicate(nrow(miRSM_SCC_MCL[[1]]), 1), 
				 cor = miRSM_SCC_MCL[[1]][, 6], 
				 pcor = miRSM_SCC_MCL[[1]][, 9], 
				 mscor = miRSM_SCC_MCL[[1]][, 10])
miRSM_SCC_MCL_modules_p.values <- sponge_compute_p_values(sponge_result = miRSM_SCC_MCL_modules, null_model = SPONGE_null_model)

miRSM_SDC_MCL_modules <- data.frame(geneA = paste("module_1", seq(nrow(miRSM_SDC_MCL[[1]]))), 
                                 geneB = paste("module_2", seq(nrow(miRSM_SDC_MCL[[1]]))), 
				 df = replicate(nrow(miRSM_SDC_MCL[[1]]), 1), 
				 cor = miRSM_SDC_MCL[[1]][, 6], 
				 pcor = miRSM_SDC_MCL[[1]][, 7], 
				 mscor = miRSM_SDC_MCL[[1]][, 8])
miRSM_SDC_MCL_modules_p.values <- sponge_compute_p_values(sponge_result = miRSM_SDC_MCL_modules, null_model = SPONGE_null_model)

# NMF
miRSM_SCC_NMF_modules <- data.frame(geneA = paste("module_1", seq(nrow(miRSM_SCC_NMF[[1]]))), 
                                 geneB = paste("module_2", seq(nrow(miRSM_SCC_NMF[[1]]))), 
				 df = replicate(nrow(miRSM_SCC_NMF[[1]]), 1), 
				 cor = miRSM_SCC_NMF[[1]][, 6], 
				 pcor = miRSM_SCC_NMF[[1]][, 9], 
				 mscor = miRSM_SCC_NMF[[1]][, 10])
miRSM_SCC_NMF_modules_p.values <- sponge_compute_p_values(sponge_result = miRSM_SCC_NMF_modules, null_model = SPONGE_null_model)

miRSM_SDC_NMF_modules <- data.frame(geneA = paste("module_1", seq(nrow(miRSM_SDC_NMF[[1]]))), 
                                 geneB = paste("module_2", seq(nrow(miRSM_SDC_NMF[[1]]))), 
				 df = replicate(nrow(miRSM_SDC_NMF[[1]]), 1), 
				 cor = miRSM_SDC_NMF[[1]][, 6], 
				 pcor = miRSM_SDC_NMF[[1]][, 7], 
				 mscor = miRSM_SDC_NMF[[1]][, 8])
miRSM_SDC_NMF_modules_p.values <- sponge_compute_p_values(sponge_result = miRSM_SDC_NMF_modules, null_model = SPONGE_null_model)

miRSM_SRVC_NMF_modules <- data.frame(geneA = paste("module_1", seq(nrow(miRSM_SRVC_NMF[[1]]))), 
                                 geneB = paste("module_2", seq(nrow(miRSM_SRVC_NMF[[1]]))), 
				 df = replicate(nrow(miRSM_SRVC_NMF[[1]]), 1), 
				 cor = miRSM_SRVC_NMF[[1]][, 6], 
				 pcor = miRSM_SRVC_NMF[[1]][, 9], 
				 mscor = miRSM_SRVC_NMF[[1]][, 10])
miRSM_SRVC_NMF_modules_p.values <- sponge_compute_p_values(sponge_result = miRSM_SRVC_NMF_modules, null_model = SPONGE_null_model)

# kmeans
miRSM_SCC_kmeans_modules <- data.frame(geneA = paste("module_1", seq(nrow(miRSM_SCC_kmeans[[1]]))), 
                                 geneB = paste("module_2", seq(nrow(miRSM_SCC_kmeans[[1]]))), 
				 df = replicate(nrow(miRSM_SCC_kmeans[[1]]), 1), 
				 cor = miRSM_SCC_kmeans[[1]][, 6], 
				 pcor = miRSM_SCC_kmeans[[1]][, 9], 
				 mscor = miRSM_SCC_kmeans[[1]][, 10])
miRSM_SCC_kmeans_modules_p.values <- sponge_compute_p_values(sponge_result = miRSM_SCC_kmeans_modules, null_model = SPONGE_null_model)

miRSM_SDC_kmeans_modules <- data.frame(geneA = paste("module_1", seq(nrow(miRSM_SDC_kmeans[[1]]))), 
                                 geneB = paste("module_2", seq(nrow(miRSM_SDC_kmeans[[1]]))), 
				 df = replicate(nrow(miRSM_SDC_kmeans[[1]]), 1), 
				 cor = miRSM_SDC_kmeans[[1]][, 6], 
				 pcor = miRSM_SDC_kmeans[[1]][, 7], 
				 mscor = miRSM_SDC_kmeans[[1]][, 8])
miRSM_SDC_kmeans_modules_p.values <- sponge_compute_p_values(sponge_result = miRSM_SDC_kmeans_modules, null_model = SPONGE_null_model)

miRSM_SRVC_kmeans_modules <- data.frame(geneA = paste("module_1", seq(nrow(miRSM_SRVC_kmeans[[1]]))), 
                                 geneB = paste("module_2", seq(nrow(miRSM_SRVC_kmeans[[1]]))), 
				 df = replicate(nrow(miRSM_SRVC_kmeans[[1]]), 1), 
				 cor = miRSM_SRVC_kmeans[[1]][, 6], 
				 pcor = miRSM_SRVC_kmeans[[1]][, 9], 
				 mscor = miRSM_SRVC_kmeans[[1]][, 10])
miRSM_SRVC_kmeans_modules_p.values <- sponge_compute_p_values(sponge_result = miRSM_SRVC_kmeans_modules, null_model = SPONGE_null_model)

# fabia
miRSM_SCC_fabia_modules <- data.frame(geneA = paste("module_1", seq(nrow(miRSM_SCC_fabia[[1]]))), 
                                 geneB = paste("module_2", seq(nrow(miRSM_SCC_fabia[[1]]))), 
				 df = replicate(nrow(miRSM_SCC_fabia[[1]]), 1), 
				 cor = miRSM_SCC_fabia[[1]][, 6], 
				 pcor = miRSM_SCC_fabia[[1]][, 9], 
				 mscor = miRSM_SCC_fabia[[1]][, 10])
miRSM_SCC_fabia_modules_p.values <- sponge_compute_p_values(sponge_result = miRSM_SCC_fabia_modules, null_model = SPONGE_null_model)

miRSM_SDC_fabia_modules <- data.frame(geneA = paste("module_1", seq(nrow(miRSM_SDC_fabia[[1]]))), 
                                 geneB = paste("module_2", seq(nrow(miRSM_SDC_fabia[[1]]))), 
				 df = replicate(nrow(miRSM_SDC_fabia[[1]]), 1), 
				 cor = miRSM_SDC_fabia[[1]][, 6], 
				 pcor = miRSM_SDC_fabia[[1]][, 7], 
				 mscor = miRSM_SDC_fabia[[1]][, 8])
miRSM_SDC_fabia_modules_p.values <- sponge_compute_p_values(sponge_result = miRSM_SDC_fabia_modules, null_model = SPONGE_null_model)

miRSM_SRVC_fabia_modules <- data.frame(geneA = paste("module_1", seq(nrow(miRSM_SRVC_fabia[[1]]))), 
                                 geneB = paste("module_2", seq(nrow(miRSM_SRVC_fabia[[1]]))), 
				 df = replicate(nrow(miRSM_SRVC_fabia[[1]]), 1), 
				 cor = miRSM_SRVC_fabia[[1]][, 6], 
				 pcor = miRSM_SRVC_fabia[[1]][, 9], 
				 mscor = miRSM_SRVC_fabia[[1]][, 10])
miRSM_SRVC_fabia_modules_p.values <- sponge_compute_p_values(sponge_result = miRSM_SRVC_fabia_modules, null_model = SPONGE_null_model)

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

# Classification performation using SVM classifier
Class <- as.matrix(c(rep(1, 72),rep(0, 72)))
colnames(Class) <- "Class"
miRSM_SDC_fabia_Modulelist <- miRSM_SDC_fabia[[2]]
miRSM_SDC_fabia_RNA1Exp <- lapply(seq_along(miRSM_SDC_fabia_Modulelist), function(i) ceRExp[, which(colnames(ceRExp) %in% miRSM_SDC_fabia_Modulelist[[i]])])
miRSM_SDC_fabia_RNA2Exp <- lapply(seq_along(miRSM_SDC_fabia_Modulelist), function(i) mRExp[, which(colnames(mRExp) %in% miRSM_SDC_fabia_Modulelist[[i]])])
miRSM_SDC_fabia_ExpData <- lapply(seq_along(miRSM_SDC_fabia_Modulelist), function(i) as.data.frame(cbind(assay(miRSM_SDC_fabia_RNA1Exp[[i]]), assay(miRSM_SDC_fabia_RNA2Exp[[i]]), Class)))
miRSM_SDC_fabia_SVMClassify_res <- lapply(seq_along(miRSM_SDC_fabia_Modulelist), function(i) SVMClassifier(miRSM_SDC_fabia_ExpData[[i]], 1:72, 73:144, nfolds = 10))

miRSM_SCC_fabia_Modulelist <- miRSM_SCC_fabia[[2]]
miRSM_SCC_fabia_RNA1Exp <- lapply(seq_along(miRSM_SCC_fabia_Modulelist), function(i) ceRExp[, which(colnames(ceRExp) %in% miRSM_SCC_fabia_Modulelist[[i]])])
miRSM_SCC_fabia_RNA2Exp <- lapply(seq_along(miRSM_SCC_fabia_Modulelist), function(i) mRExp[, which(colnames(mRExp) %in% miRSM_SCC_fabia_Modulelist[[i]])])
miRSM_SCC_fabia_ExpData <- lapply(seq_along(miRSM_SCC_fabia_Modulelist), function(i) as.data.frame(cbind(assay(miRSM_SCC_fabia_RNA1Exp[[i]]), assay(miRSM_SCC_fabia_RNA2Exp[[i]]), Class)))
miRSM_SCC_fabia_SVMClassify_res <- lapply(seq_along(miRSM_SCC_fabia_Modulelist), function(i) SVMClassifier(miRSM_SCC_fabia_ExpData[[i]], 1:72, 73:144, nfolds = 10))

miRSM_SRVC_fabia_Modulelist <- miRSM_SRVC_fabia[[2]]
miRSM_SRVC_fabia_RNA1Exp <- lapply(seq_along(miRSM_SRVC_fabia_Modulelist), function(i) ceRExp[, which(colnames(ceRExp) %in% miRSM_SRVC_fabia_Modulelist[[i]])])
miRSM_SRVC_fabia_RNA2Exp <- lapply(seq_along(miRSM_SRVC_fabia_Modulelist), function(i) mRExp[, which(colnames(mRExp) %in% miRSM_SRVC_fabia_Modulelist[[i]])])
miRSM_SRVC_fabia_ExpData <- lapply(seq_along(miRSM_SRVC_fabia_Modulelist), function(i) as.data.frame(cbind(assay(miRSM_SRVC_fabia_RNA1Exp[[i]]), assay(miRSM_SRVC_fabia_RNA2Exp[[i]]), Class)))
miRSM_SRVC_fabia_SVMClassify_res <- lapply(seq_along(miRSM_SRVC_fabia_Modulelist), function(i) SVMClassifier(miRSM_SRVC_fabia_ExpData[[i]], 1:72, 73:144, nfolds = 10))

##########################################################################################################################
########################### Comparison with graph-based clustering methods in the case study #############################
##########################################################################################################################
## Identifying miRNA sponge modules by using SC and SPONGE methods
SC_ceRNA_network <- SC(assay(miRExp), assay(ceRExp), assay(mRExp), assay(miRTarget), 
		       pvaluecutoff = 0.05, senscorcutoff = 0.1)

SPONGE_ceRNA_network <- SPONGE(assay(miRExp), assay(ceRExp), assay(mRExp), assay(miRTarget), 
			       pvaluecutoff = 0.05, null_model = SPONGE_null_model)

SPONGE_ceRNA_module_MCL <- netModule(SPONGE_ceRNA_network[, 1:2], method = "MCL", modulesize = 4)
SPONGE_ceRNA_module_MCODE <- netModule(SPONGE_ceRNA_network[, 1:2], method = "MCODE", modulesize = 4)
SPONGE_ceRNA_module_FN <- netModule(SPONGE_ceRNA_network[, 1:2], method = "FN", modulesize = 4)
SPONGE_ceRNA_module_LINKCOMM <- netModule(SPONGE_ceRNA_network[, 1:2], method = "LINKCOMM", modulesize = 4)

SPONGE_ceRNA_module_MCL <- CandModgenes(assay(ceRExp), assay(mRExp), SPONGE_ceRNA_module_MCL, num.ModuleceRs = 2, num.ModulemRs = 2)
SPONGE_ceRNA_module_MCODE <- CandModgenes(assay(ceRExp), assay(mRExp), SPONGE_ceRNA_module_MCODE, num.ModuleceRs = 2, num.ModulemRs = 2)
SPONGE_ceRNA_module_FN <- CandModgenes(assay(ceRExp), assay(mRExp), SPONGE_ceRNA_module_FN, num.ModuleceRs = 2, num.ModulemRs = 2)
SPONGE_ceRNA_module_LINKCOMM <- CandModgenes(assay(ceRExp), assay(mRExp), SPONGE_ceRNA_module_LINKCOMM, num.ModuleceRs = 2, num.ModulemRs = 2)

SC_ceRNA_module_MCL <- netModule(SC_ceRNA_network[, 1:2], method = "MCL", modulesize = 4)
SC_ceRNA_module_MCODE <- netModule(SC_ceRNA_network[, 1:2], method = "MCODE", modulesize = 4)
SC_ceRNA_module_FN <- netModule(SC_ceRNA_network[, 1:2], method = "FN", modulesize = 4)
SC_ceRNA_module_LINKCOMM <- netModule(SC_ceRNA_network[, 1:2], method = "LINKCOMM", modulesize = 4)

SC_ceRNA_module_MCL <- CandModgenes(assay(ceRExp), assay(mRExp), SC_ceRNA_module_MCL, num.ModuleceRs = 2, num.ModulemRs = 2)
SC_ceRNA_module_MCODE <- CandModgenes(assay(ceRExp), assay(mRExp), SC_ceRNA_module_MCODE, num.ModuleceRs = 2, num.ModulemRs = 2)
SC_ceRNA_module_FN <- CandModgenes(assay(ceRExp), assay(mRExp), SC_ceRNA_module_FN, num.ModuleceRs = 2, num.ModulemRs = 2)
SC_ceRNA_module_LINKCOMM <- CandModgenes(assay(ceRExp), assay(mRExp), SC_ceRNA_module_LINKCOMM, num.ModuleceRs = 2, num.ModulemRs = 2)

## Modular analysis of the identified miRNA sponge modules by graph-based clustering methods
# Functional analysis of miRNA sponge modules
SPONGE_ceRNA_module_MCL_FEA <- module_FA(SPONGE_ceRNA_module_MCL, Analysis.type = "FEA")
SPONGE_ceRNA_module_MCL_DEA <- module_FA(SPONGE_ceRNA_module_MCL, Analysis.type = "DEA")
SPONGE_ceRNA_module_MCODE_FEA <- module_FA(SPONGE_ceRNA_module_MCODE, Analysis.type = "FEA")
SPONGE_ceRNA_module_MCODE_DEA <- module_FA(SPONGE_ceRNA_module_MCODE, Analysis.type = "DEA")
SPONGE_ceRNA_module_FN_FEA <- module_FA(SPONGE_ceRNA_module_FN, Analysis.type = "FEA")
SPONGE_ceRNA_module_FN_DEA <- module_FA(SPONGE_ceRNA_module_FN, Analysis.type = "DEA")
SPONGE_ceRNA_module_LINKCOMM_FEA <- module_FA(SPONGE_ceRNA_module_LINKCOMM, Analysis.type = "FEA")
SPONGE_ceRNA_module_LINKCOMM_DEA <- module_FA(SPONGE_ceRNA_module_LINKCOMM, Analysis.type = "DEA")

SC_ceRNA_module_MCL_FEA <- module_FA(SC_ceRNA_module_MCL, Analysis.type = "FEA")
SC_ceRNA_module_MCL_DEA <- module_FA(SC_ceRNA_module_MCL, Analysis.type = "DEA")
SC_ceRNA_module_MCODE_FEA <- module_FA(SC_ceRNA_module_MCODE, Analysis.type = "FEA")
SC_ceRNA_module_MCODE_DEA <- module_FA(SC_ceRNA_module_MCODE, Analysis.type = "DEA")
SC_ceRNA_module_FN_FEA <- module_FA(SC_ceRNA_module_FN, Analysis.type = "FEA")
SC_ceRNA_module_FN_DEA <- module_FA(SC_ceRNA_module_FN, Analysis.type = "DEA")
SC_ceRNA_module_LINKCOMM_FEA <- module_FA(SC_ceRNA_module_LINKCOMM, Analysis.type = "FEA")
SC_ceRNA_module_LINKCOMM_DEA <- module_FA(SC_ceRNA_module_LINKCOMM, Analysis.type = "DEA")

# BRCA enrichment analysis of miRNA sponge modules
SPONGE_ceRNA_module_MCL_pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, SPONGE_ceRNA_module_MCL)
SPONGE_ceRNA_module_MCODE_pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, SPONGE_ceRNA_module_MCODE)
SPONGE_ceRNA_module_FN_pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, SPONGE_ceRNA_module_FN)
SPONGE_ceRNA_module_LINKCOMM_pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, SPONGE_ceRNA_module_LINKCOMM)

SC_ceRNA_module_MCL_pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, SC_ceRNA_module_MCL)
SC_ceRNA_module_MCODE_pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, SC_ceRNA_module_MCODE)
SC_ceRNA_module_FN_pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, SC_ceRNA_module_FN)
SC_ceRNA_module_LINKCOMM_pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, SC_ceRNA_module_LINKCOMM)

# Validation of miRNA sponge interactions in miRNA sponge modules
SPONGE_ceRNA_module_MCL_validate <- module_Validate(SPONGE_ceRNA_module_MCL, Groundtruth)
SPONGE_ceRNA_module_MCODE_validate <- module_Validate(SPONGE_ceRNA_module_MCODE, Groundtruth)
SPONGE_ceRNA_module_FN_validate <- module_Validate(SPONGE_ceRNA_module_FN, Groundtruth)
SPONGE_ceRNA_module_LINKCOMM_validate <- module_Validate(SPONGE_ceRNA_module_LINKCOMM, Groundtruth)

SC_ceRNA_module_MCL_validate <- module_Validate(SC_ceRNA_module_MCL, Groundtruth)
SC_ceRNA_module_MCODE_validate <- module_Validate(SC_ceRNA_module_MCODE, Groundtruth)
SC_ceRNA_module_FN_validate <- module_Validate(SC_ceRNA_module_FN, Groundtruth)
SC_ceRNA_module_LINKCOMM_validate <- module_Validate(SC_ceRNA_module_LINKCOMM, Groundtruth)

# Co-expression analysis of miRNA sponge modules
SPONGE_ceRNA_module_MCL_Coexpress <- module_Coexpress(ceRExp, mRExp, SPONGE_ceRNA_module_MCL, resample = 1000, method = "mean", test.method = "t.test")
SPONGE_ceRNA_module_MCODE_Coexpress <- module_Coexpress(ceRExp, mRExp, SPONGE_ceRNA_module_MCODE, resample = 1000, method = "mean", test.method = "t.test")
SPONGE_ceRNA_module_FN_Coexpress <- module_Coexpress(ceRExp, mRExp, SPONGE_ceRNA_module_FN, resample = 1000, method = "mean", test.method = "t.test")
SPONGE_ceRNA_module_LINKCOMM_Coexpress <- module_Coexpress(ceRExp, mRExp, SPONGE_ceRNA_module_LINKCOMM, resample = 1000, method = "mean", test.method = "t.test")

SC_ceRNA_module_MCL_Coexpress <- module_Coexpress(ceRExp, mRExp, SC_ceRNA_module_MCL, resample = 1000, method = "mean", test.method = "t.test")
SC_ceRNA_module_MCODE_Coexpress <- module_Coexpress(ceRExp, mRExp, SC_ceRNA_module_MCODE, resample = 1000, method = "mean", test.method = "t.test")
SC_ceRNA_module_FN_Coexpress <- module_Coexpress(ceRExp, mRExp, SC_ceRNA_module_FN, resample = 1000, method = "mean", test.method = "t.test")
SC_ceRNA_module_LINKCOMM_Coexpress <- module_Coexpress(ceRExp, mRExp, SC_ceRNA_module_LINKCOMM, resample = 1000, method = "mean", test.method = "t.test")

# SVM classification
Modulelist_MCL <- SPONGE_ceRNA_module_MCL
SPONGE_ceRNA_module_MCL_RNA1Exp <- lapply(seq_along(Modulelist_MCL), function(i) ceRExp[, which(colnames(ceRExp) %in% Modulelist_MCL[[i]])])
SPONGE_ceRNA_module_MCL_RNA2Exp <- lapply(seq_along(Modulelist_MCL), function(i) mRExp[, which(colnames(mRExp) %in% Modulelist_MCL[[i]])])
SPONGE_ceRNA_module_MCL_ExpData <- lapply(seq_along(Modulelist_MCL), function(i) as.data.frame(cbind(assay(SPONGE_ceRNA_module_MCL_RNA1Exp[[i]]), assay(SPONGE_ceRNA_module_MCL_RNA2Exp[[i]]), Class)))
SPONGE_ceRNA_module_MCL_SVMClassify_res <- lapply(seq_along(Modulelist_MCL), function(i) SVMClassifier(SPONGE_ceRNA_module_MCL_ExpData[[i]], 1:72, 73:144, nfolds = 10))

Modulelist_MCODE <- SPONGE_ceRNA_module_MCODE
SPONGE_ceRNA_module_MCODE_RNA1Exp <- lapply(seq_along(Modulelist_MCODE), function(i) ceRExp[, which(colnames(ceRExp) %in% Modulelist_MCODE[[i]])])
SPONGE_ceRNA_module_MCODE_RNA2Exp <- lapply(seq_along(Modulelist_MCODE), function(i) mRExp[, which(colnames(mRExp) %in% Modulelist_MCODE[[i]])])
SPONGE_ceRNA_module_MCODE_ExpData <- lapply(seq_along(Modulelist_MCODE), function(i) as.data.frame(cbind(assay(SPONGE_ceRNA_module_MCODE_RNA1Exp[[i]]), assay(SPONGE_ceRNA_module_MCODE_RNA2Exp[[i]]), Class)))
SPONGE_ceRNA_module_MCODE_SVMClassify_res <- lapply(seq_along(Modulelist_MCODE), function(i) SVMClassifier(SPONGE_ceRNA_module_MCODE_ExpData[[i]], 1:72, 73:144, nfolds = 10))

Modulelist_FN <- SPONGE_ceRNA_module_FN
SPONGE_ceRNA_module_FN_RNA1Exp <- lapply(seq_along(Modulelist_FN), function(i) ceRExp[, which(colnames(ceRExp) %in% Modulelist_FN[[i]])])
SPONGE_ceRNA_module_FN_RNA2Exp <- lapply(seq_along(Modulelist_FN), function(i) mRExp[, which(colnames(mRExp) %in% Modulelist_FN[[i]])])
SPONGE_ceRNA_module_FN_ExpData <- lapply(seq_along(Modulelist_FN), function(i) as.data.frame(cbind(assay(SPONGE_ceRNA_module_FN_RNA1Exp[[i]]), assay(SPONGE_ceRNA_module_FN_RNA2Exp[[i]]), Class)))
SPONGE_ceRNA_module_FN_SVMClassify_res <- lapply(seq_along(Modulelist_FN), function(i) SVMClassifier(SPONGE_ceRNA_module_FN_ExpData[[i]], 1:72, 73:144, nfolds = 10))

Modulelist_LINKCOMM <- SPONGE_ceRNA_module_LINKCOMM
SPONGE_ceRNA_module_LINKCOMM_RNA1Exp <- lapply(seq_along(Modulelist_LINKCOMM), function(i) ceRExp[, which(colnames(ceRExp) %in% Modulelist_LINKCOMM[[i]])])
SPONGE_ceRNA_module_LINKCOMM_RNA2Exp <- lapply(seq_along(Modulelist_LINKCOMM), function(i) mRExp[, which(colnames(mRExp) %in% Modulelist_LINKCOMM[[i]])])
SPONGE_ceRNA_module_LINKCOMM_ExpData <- lapply(seq_along(Modulelist_LINKCOMM), function(i) as.data.frame(cbind(assay(SPONGE_ceRNA_module_LINKCOMM_RNA1Exp[[i]]), assay(SPONGE_ceRNA_module_LINKCOMM_RNA2Exp[[i]]), Class)))
SPONGE_ceRNA_module_LINKCOMM_SVMClassify_res <- lapply(seq_along(Modulelist_LINKCOMM), function(i) SVMClassifier(SPONGE_ceRNA_module_LINKCOMM_ExpData[[i]], 1:72, 73:144, nfolds = 10))

Modulelist_MCL <- SC_ceRNA_module_MCL
SC_ceRNA_module_MCL_RNA1Exp <- lapply(seq_along(Modulelist_MCL), function(i) ceRExp[, which(colnames(ceRExp) %in% Modulelist_MCL[[i]])])
SC_ceRNA_module_MCL_RNA2Exp <- lapply(seq_along(Modulelist_MCL), function(i) mRExp[, which(colnames(mRExp) %in% Modulelist_MCL[[i]])])
SC_ceRNA_module_MCL_ExpData <- lapply(seq_along(Modulelist_MCL), function(i) as.data.frame(cbind(assay(SC_ceRNA_module_MCL_RNA1Exp[[i]]), assay(SC_ceRNA_module_MCL_RNA2Exp[[i]]), Class)))
SC_ceRNA_module_MCL_SVMClassify_res <- lapply(seq_along(Modulelist_MCL), function(i) SVMClassifier(SC_ceRNA_module_MCL_ExpData[[i]], 1:72, 73:144, nfolds = 10))

Modulelist_MCODE <- SC_ceRNA_module_MCODE
SC_ceRNA_module_MCODE_RNA1Exp <- lapply(seq_along(Modulelist_MCODE), function(i) ceRExp[, which(colnames(ceRExp) %in% Modulelist_MCODE[[i]])])
SC_ceRNA_module_MCODE_RNA2Exp <- lapply(seq_along(Modulelist_MCODE), function(i) mRExp[, which(colnames(mRExp) %in% Modulelist_MCODE[[i]])])
SC_ceRNA_module_MCODE_ExpData <- lapply(seq_along(Modulelist_MCODE), function(i) as.data.frame(cbind(assay(SC_ceRNA_module_MCODE_RNA1Exp[[i]]), assay(SC_ceRNA_module_MCODE_RNA2Exp[[i]]), Class)))
SC_ceRNA_module_MCODE_SVMClassify_res <- lapply(seq_along(Modulelist_MCODE), function(i) SVMClassifier(SC_ceRNA_module_MCODE_ExpData[[i]], 1:72, 73:144, nfolds = 10))

Modulelist_FN <- SC_ceRNA_module_FN
SC_ceRNA_module_FN_RNA1Exp <- lapply(seq_along(Modulelist_FN), function(i) ceRExp[, which(colnames(ceRExp) %in% Modulelist_FN[[i]])])
SC_ceRNA_module_FN_RNA2Exp <- lapply(seq_along(Modulelist_FN), function(i) mRExp[, which(colnames(mRExp) %in% Modulelist_FN[[i]])])
SC_ceRNA_module_FN_ExpData <- lapply(seq_along(Modulelist_FN), function(i) as.data.frame(cbind(assay(SC_ceRNA_module_FN_RNA1Exp[[i]]), assay(SC_ceRNA_module_FN_RNA2Exp[[i]]), Class)))
SC_ceRNA_module_FN_SVMClassify_res <- lapply(seq_along(Modulelist_FN), function(i) SVMClassifier(SC_ceRNA_module_FN_ExpData[[i]], 1:72, 73:144, nfolds = 10))

Modulelist_LINKCOMM <- SC_ceRNA_module_LINKCOMM
SC_ceRNA_module_LINKCOMM_RNA1Exp <- lapply(seq_along(Modulelist_LINKCOMM), function(i) ceRExp[, which(colnames(ceRExp) %in% Modulelist_LINKCOMM[[i]])])
SC_ceRNA_module_LINKCOMM_RNA2Exp <- lapply(seq_along(Modulelist_LINKCOMM), function(i) mRExp[, which(colnames(mRExp) %in% Modulelist_LINKCOMM[[i]])])
SC_ceRNA_module_LINKCOMM_ExpData <- lapply(seq_along(Modulelist_LINKCOMM), function(i) as.data.frame(cbind(assay(SC_ceRNA_module_LINKCOMM_RNA1Exp[[i]]), assay(SC_ceRNA_module_LINKCOMM_RNA2Exp[[i]]), Class)))
SC_ceRNA_module_LINKCOMM_SVMClassify_res <- lapply(seq_along(Modulelist_LINKCOMM), function(i) SVMClassifier(SC_ceRNA_module_LINKCOMM_ExpData[[i]], 1:72, 73:144, nfolds = 10))

save.image("TCGA_BRCA_miRSM.RData")
