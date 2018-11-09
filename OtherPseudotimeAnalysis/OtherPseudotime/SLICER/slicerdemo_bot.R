#Slicer run
library(SLICER)
require("lle")

#Use the same data as other methods
HSMM_gene_annotation <- read.delim("/Users/christopher_penfold/Desktop/BranchingGPs/demos/Botrytis_featureData.txt",row.names = 1)
HSMM_sample_sheet    <- read.delim("/Users/christopher_penfold/Desktop/BranchingGPs/demos/Botrytis_PhenoData.txt", header=TRUE)
HSMM_expr_matrix     <- read.table("/Users/christopher_penfold/Desktop/BranchingGPs/demos/Botrytis_data.txt",row.names = 1, header=TRUE)
colnames(HSMM_gene_annotation) <- c("gene_short_name")
rownames(HSMM_sample_sheet) <- colnames(HSMM_expr_matrix)

slicer_genes <- select_genes(t(HSMM_expr_matrix))
k <- select_k(t(HSMM_expr_matrix), kmin = 30, kmax=60)
slicer_traj_lle <- lle(t(HSMM_expr_matrix), m = 2, k)$Y
#plot(slicer_traj_lle, xlab = "LLE Comp 1", ylab = "LLE Comp 2",
#     main = "Locally linear embedding of cells from SLICER", 
#      pch=16)

slicer_traj_graph <- conn_knn_graph(slicer_traj_lle, 10)
#plot(slicer_traj_graph, main = "Fully connected kNN graph from SLICER")
ends <- find_extreme_cells(slicer_traj_graph, slicer_traj_lle)
start <- ends[2]

pseudotime_order_slicer <- cell_order(slicer_traj_graph, start)
df <- data.frame(colnames(HSMM_expr_matrix), pseudotime_order_slicer)
write.table(df, "~/Desktop/slicer_pseudotime.txt", sep="\t")

start <- ends[1]
pseudotime_order_slicer <- cell_order(slicer_traj_graph, start)
df <- data.frame(colnames(HSMM_expr_matrix), pseudotime_order_slicer)
write.table(df, "~/Desktop/slicer_pseudotimer.txt", sep="\t")

#Individual runs for PGCs
HSMM_gene_annotation <- read.delim("/Users/christopher_penfold/Desktop/BranchingGPs/demos/PGC_transcriptomics_featureData.txt",row.names = 1)
HSMM_sample_sheet    <- read.delim("/Users/christopher_penfold/Desktop/BranchingGPs/demos/PGC_transcriptomics_data_PhenoData.txt", header=TRUE)
HSMM_expr_matrix     <- read.table("/Users/christopher_penfold/Desktop/BranchingGPs/demos/PGC_transcriptomics_data.txt",row.names = 1, header=TRUE)
colnames(HSMM_gene_annotation) <- c("gene_short_name")
rownames(HSMM_sample_sheet) <- colnames(HSMM_expr_matrix)
HSMM_expr_matrix[,which(HSMM_sample_sheet$Type != 2)]
HSMM_sample_sheet[which(HSMM_sample_sheet$Type != 2),]

#slicer_genes <- select_genes(t(HSMM_expr_matrix))
k <- select_k(t(HSMM_expr_matrix), kmin = 30, kmax=60)
slicer_traj_lle <- lle(t(HSMM_expr_matrix), m = 2, k)$Y
#plot(slicer_traj_lle, xlab = "LLE Comp 1", ylab = "LLE Comp 2",
#     main = "Locally linear embedding of cells from SLICER", 
#     pch=16)

slicer_traj_graph <- conn_knn_graph(slicer_traj_lle, 10)
#plot(slicer_traj_graph, main = "Fully connected kNN graph from SLICER")
ends <- find_extreme_cells(slicer_traj_graph, slicer_traj_lle)
start <- ends[2]
pseudotime_order_slicer <- cell_order(slicer_traj_graph, start)
df <- data.frame(colnames(HSMM_expr_matrix), pseudotime_order_slicer)
write.table(df, "~/Desktop/slicer_pseudotime_PGC.txt", sep="\t")
start <- ends[1]

pseudotime_order_slicer <- cell_order(slicer_traj_graph, start)
df <- data.frame(colnames(HSMM_expr_matrix), pseudotime_order_slicer)
write.table(df, "~/Desktop/slicer_pseudotime_PGCr.txt", sep="\t")


#Individual runs for soma
HSMM_gene_annotation <- read.delim("/Users/christopher_penfold/Desktop/BranchingGPs/demos/PGC_transcriptomics_featureData.txt",row.names = 1)
HSMM_sample_sheet    <- read.delim("/Users/christopher_penfold/Desktop/BranchingGPs/demos/PGC_transcriptomics_data_PhenoData.txt", header=TRUE)
HSMM_expr_matrix     <- read.table("/Users/christopher_penfold/Desktop/BranchingGPs/demos/PGC_transcriptomics_data.txt",row.names = 1, header=TRUE)
colnames(HSMM_gene_annotation) <- c("gene_short_name")
rownames(HSMM_sample_sheet) <- colnames(HSMM_expr_matrix)
HSMM_expr_matrix[,which(HSMM_sample_sheet$Type != 1)]
HSMM_sample_sheet[which(HSMM_sample_sheet$Type != 1),]

#slicer_genes <- select_genes(t(HSMM_expr_matrix))
k <- select_k(t(HSMM_expr_matrix), kmin = 30, kmax=60)
slicer_traj_lle <- lle(t(HSMM_expr_matrix), m = 2, k)$Y
#plot(slicer_traj_lle, xlab = "LLE Comp 1", ylab = "LLE Comp 2",
#     main = "Locally linear embedding of cells from SLICER", 
#     pch=16)


slicer_traj_graph <- conn_knn_graph(slicer_traj_lle, 10)
plot(slicer_traj_graph, main = "Fully connected kNN graph from SLICER")
ends <- find_extreme_cells(slicer_traj_graph, slicer_traj_lle)
start <- ends[2]
pseudotime_order_slicer <- cell_order(slicer_traj_graph, start)
df <- data.frame(colnames(HSMM_expr_matrix), pseudotime_order_slicer)
write.table(df, "~/Desktop/slicer_pseudotime_soma.txt", sep="\t")

start <- ends[1]
pseudotime_order_slicer <- cell_order(slicer_traj_graph, start)
df <- data.frame(colnames(HSMM_expr_matrix), pseudotime_order_slicer)
write.table(df, "~/Desktop/slicer_pseudotime_somar.txt", sep="\t")