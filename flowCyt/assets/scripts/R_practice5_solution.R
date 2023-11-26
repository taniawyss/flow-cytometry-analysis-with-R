#############################
# Let’s practice – 5        #
# Differential testing      #
# Solutions                 #
#############################


# clear the environment
rm(list = ls())


# load libraries
library(CATALYST)


# 1) Load the sce object with UMAP from previous exercise
load("course_datasets/FR_FCM_Z3WR/sce_UMAP.RData")

# 2) Apply FlowSOM clustering + ConsensusClusterPlus metaclustering

set.seed(1234)
sce <- cluster(sce, 
               features = "type")

names(cluster_codes(sce))
# [1] "som100" "meta2"  "meta3"  "meta4"  "meta5"  "meta6"  "meta7"  "meta8"  "meta9"  "meta10"
# [11] "meta11" "meta12" "meta13" "meta14" "meta15" "meta16" "meta17" "meta18" "meta19" "meta20"

# Plot which shows the change in area under the Consensus Cumulative Distribution function (CDF) per k
sce@metadata$delta_area
# Based on this plot, after k=8, there is not much change anymore.
# We could select a k between 5 and 8.

# 3) Plot UMAP with clusters, expression heatmap and ridge plots

# UMAP
plotDR(sce,
       dr = "UMAP",
       color_by = "meta8" )

# Heatmap 

plotExprHeatmap(sce,row_clust = F, col_clust = F,
                features = "type",
                by="cluster_id",
                k="meta8")


# Ridge plots

plotClusterExprs(sce,
                 k="meta10")

# 4) Rename / merge clusters
# create a merging table
merging_table <- data.frame(old_cluster = 1:8,
                            new_cluster = c("CD4 T cells","Other","Monocytes","Other","CD8 T cells","CD8 T cells","Other","B cells"))

# write to file
write.csv(merging_table,file = "course_datasets/FR_FCM_Z3WR/merging_table.csv",quote = F,row.names = F)


# annotate clusters
sce <- mergeClusters(sce,
                     k="meta8",
                     table = merging_table,
                     id = "Major_cell_populations")

# 5) UMAP with major cell populations
plotDR(sce,
       dr = "UMAP",
       color_by = "Major_cell_populations")

plotExprHeatmap(sce,row_clust = F, col_clust = F,
                features = "type",
                by="cluster_id",
                k="Major_cell_populations")


# save
save(sce,file = "course_datasets/FR_FCM_Z3WR/sce_annotated.RData" )

