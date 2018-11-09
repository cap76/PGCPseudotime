library(monocle)
set.seed(12345)

#Read in the data
HSMM_gene_annotation <- read.delim("/Users/christopher_penfold/Desktop/BranchingGPs/demos/PGC_transcriptomics_featureData.txt",row.names = 1)
HSMM_sample_sheet <- read.delim("/Users/christopher_penfold/Desktop/BranchingGPs/demos/PGC_transcriptomics_data_PhenoData.txt", header=TRUE)
HSMM_expr_matrix <- read.table("/Users/christopher_penfold/Desktop/BranchingGPs/demos/PGC_transcriptomics_data.txt",row.names = 1, header=TRUE)

#HSMM_expre_matrix <- HSMM_expre_matrix + matrix( rnorm(3573*418,mean=0,sd=.001),3573,418) 

colnames(HSMM_gene_annotation) <- c("gene_short_name")
rownames(HSMM_sample_sheet) <- colnames(HSMM_expr_matrix)

pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)

#Converting TPM/FPKM values into mRNA counts
HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0.0,
                       expressionFamily = tobit(Lower = 0.0))
#rpc_matrix <- relative2abs(HSMM, method = "num_genes")
#HSMM <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
#                       phenoData = pd,
#                       featureData = fd,
#                       lowerDetectionLimit = 0.5,
#                       expressionFamily = negbinomial.size())

#Estimate size factors and dispersions
HSMM <- estimateSizeFactors(HSMM)
#HSMM <- estimateDispersions(HSMM)

#HSMM <- detectGenes(HSMM, min_expr = 0.1)
#fData(HSMM)$use_for_ordering <- fData(HSMM)$num_cells_expressed > 0.05 * ncol(HSMM)
#plot_pc_variance_explained(HSMM, return_all = F)

#HSMM <- reduceDimension(HSMM,
#                            max_components = 2,
#                            norm_method = 'log',
#                            num_dim = 3,
#                            reduction_method = 'tSNE',
#                            verbose = T)


cth <- newCellTypeHierarchy()
SOX17_id <- row.names(subset(fData(HSMM), gene_short_name == "SOX17"))
WT1_id <- row.names(subset(fData(HSMM), gene_short_name == "WT1"))
DPPA3_id <- row.names(subset(fData(HSMM), gene_short_name == "DPPA3"))
cth <- addCellType(cth, "PGC", classify_func = function(x) { x[SOX17_id,] >= .3})
cth <- addCellType(cth, "Soma", classify_func = function(x) { x[WT1_id,] >= .2})
cth <- addCellType(cth, "ES", classify_func = function(x) { x[SOX17_id,] < .4 & x[DPPA3_id,] >= 2})

HSMM <- classifyCells(HSMM, cth, 0.1)

print(head(pData(HSMM)))

#marker_diff <- markerDiffTable(HSMM,
#                               cth,
#                               residualModelFormulaStr = "~Type + num_genes_expressed",
#                               cores = 1)

#candidate_clustering_genes <- row.names(subset(marker_diff, qval < 0.01))
#marker_spec <- calculateMarkerSpecificity(HSMM[candidate_clustering_genes,], cth)
marker_diff <- markerDiffTable(HSMM,cth,cores = 1)

candidate_clustering_genes <-row.names(subset(marker_diff, qval < 0.01))
marker_spec <-calculateMarkerSpecificity(HSMM[candidate_clustering_genes,], cth)
head(selectTopMarkers(marker_spec, 3))
semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 20)$gene_id)
HSMM <- setOrderingFilter(HSMM, semisup_clustering_genes)
#plot_ordering_genes(HSMM)
plot_pc_variance_explained(HSMM, return_all = F)

HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 3,
                        norm_method = 'log',
                        reduction_method = 'tSNE',
                        method = 'DDRTree',
                        verbose = T)
HSMM <- clusterCells(HSMM, num_clusters = 5)

plot_cell_clusters(HSMM, 1, 2, color = "CellType")
plot_cell_clusters(HSMM, 1, 2, color_by = 'as.factor(Cluster)')
plot_cell_clusters(HSMM, 1, 2, color_by = 'as.factor(Time)')
plot_cell_clusters(HSMM, 1, 2, color_by = 'as.factor(Type)')


#HSMM <- clusterCells(HSMM,
#                     num_clusters = 2,
#                     frequency_thresh = 0.1,
#                     cell_type_hierarchy = cth)
#plot_cell_clusters(HSMM, 1, 2, color = "CellType",
#                   markers = c("SOX17", "DPPA3"))



#clustering_DEG_genes <- differentialGeneTest(HSMM,fullModelFormulaStr = '~Cluster',cores = 1)
##Pseudotime (approach 1)
##diff_test_res <- differentialGeneTest(HSMM)
##ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
##HSMM <- setOrderingFilter(HSMM, ordering_genes)

##HSMM <- reduceDimension(HSMM, max_components = 2,method = 'DDRTree')
##HSMM <- orderCells(HSMM)
##plot_cell_trajectory(HSMM, color_by = "Time")

GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$Type)[,"0"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
##HSMM <- orderCells(HSMM, root_state = GM_state(HSMM))
##plot_cell_trajectory(HSMM, color_by = "Pseudotime")

##plot_cell_trajectory(HSMM, color_by = "State") +
##  facet_wrap(~State, nrow = 1)

##my_genes <- row.names(subset(fData(HSMM),
##                             gene_short_name %in% c("SOX17", "PRDM1", "DPPA3")))
##cds_subset <- HSMM[my_genes,]
##plot_genes_in_pseudotime(cds_subset, color_by = "Time")

#Pseudotime 2
clustering_DEG_genes <- differentialGeneTest(HSMM,fullModelFormulaStr = '~Cluster',cores = 1)
HSMM_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
HSMM <- setOrderingFilter(HSMM,ordering_genes = HSMM_ordering_genes)
HSMM <-reduceDimension(HSMM, method = 'ICA')
HSMM <- orderCells(HSMM)
HSMM <- orderCells(HSMM, root_state = GM_state(HSMM), num_paths=2)
plot_cell_trajectory(HSMM, color_by = "Time")



#HSMM_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:20]
#HSMM <- setOrderingFilter(HSMM, ordering_genes = HSMM_ordering_genes)
#HSMM <- reduceDimension(HSMM, method = 'DDRTree')
#HSMM <- orderCells(HSMM)
#HSMM <-orderCells(HSMM, root_state = GM_state(HSMM))
#plot_cell_trajectory(HSMM, color_by = "Time")

my_genes <- row.names(subset(fData(HSMM),
                             gene_short_name %in% c("SOX17", "DPPA3", "WT1")))
HSMM_subset <- HSMM[my_genes,]
plot_genes_in_pseudotime(HSMM_subset, color_by = "Time")

meh<-print((pData(HSMM)))
write.table(as.data.frame(meh), "~/Desktop/Monocle_allprocesses2.txt", sep="\t")

#to_be_tested <- row.names(subset(fData(HSMM),gene_short_name %in% HSMM_gene_annotation$gene_short_name))
#cds_subset <- HSMM[to_be_tested,]
#diff_test_res <- differentialGeneTest(cds_subset, fullModelFormulaStr = "~Type")
#diff_test_res <- differentialGeneTest(cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)")


#Analysing branches in singe-cell trajectories
#plot_cell_trajectory(HSMM, color_by = "Time")
BEAM_res <- BEAM(HSMM, branch_point = 2, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
#HSMM_gene_annotation$gene_short_name
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
PGC_genes <- row.names(subset(fData(HSMM),
                              gene_short_name %in% c("SOX17","PRDM1","WT1","DPPA3","SOX2")))
plot_genes_branched_pseudotime(HSMM[PGC_genes,],
                               branch_point = 2,
                               color_by = "Time",
                               ncol = 1)

#branch_ratios <- calILRs(cds, trend_formula = "~sm.ns(Pseudotime)*Branch",
#        branch_point = 1)
#branch_times <- detectBifurcationPoint(branch_ratios)
