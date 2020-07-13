library(ggplot2)
library(readxl)
library(data.table)
library(tidyverse)

g <- glimpse 
s <- summary
h <- head

################################## Importing and Inspecting Data ################################## 

# Importing data
ag <- as.data.table(read_excel("counts.xlsx", skip = 1))

colnames(ag) <- c("gene", colnames(ag[, -1]))


# agMat: count matrix
agMat <- ag[, -1]


# Library and Coverage Table 
Cov_Table <- data.table(coverage = colMeans(agMat > 0),
                        cell_name = colnames(agMat),
                        LibSize = colSums(agMat)) %>%
        separate(cell_name, c("age",
                              "strain",
                              "cell_type",
                              "id"), 
                 sep = "_")

Cov_Table <- Cov_Table[cell_type != "population", 
                       cell_name := paste(age, strain, cell_type, id, sep = "_")][
                               !is.na(cell_name)
                       ]

# Plotting distribution of library size per cell
Dist_LibSize <- ggplot(Cov_Table,
                       aes(x = LibSize,
                           fill = "blue")) + 
        geom_density(alpha = 0.5) +
        geom_vline(xintercept = 20000, color = "blue", size = 1) + 
        theme_bw() + 
        theme(legend.position = "none") +
        ggtitle("Distribution of Library Size") + 
        xlab("Total Counts per Cell (Log)") + 
        ylab("Density") +
        scale_x_log10()

# Plotting coverage
coverage <- ggplot(Cov_Table,
                   aes(x = coverage,
                       fill = "red")) +
        geom_density(alpha = 0.5) +
        theme_bw() + 
        theme(legend.position = "none") +
        ggtitle("Cell Coverage") + 
        xlab("Coverage") +
        ylab("Density")


ag <- select(ag, c("gene", Cov_Table$cell_name))

erccRow <- ag[str_detect(ag$gene, "Ercc")]


################################## Creating an SCE object ################################## 

library(SingleCellExperiment)


# sce object 
sce <- SingleCellExperiment(assays = list(counts = as.matrix(ag[, -1])),
                            rowData = data.frame(gene_names = ag$gene),
                            colData = data.frame(cell_names = colnames(ag[, -1])))


rownames(sce) <- rowData(sce)$gene_names

################################## Quality Control ################################## 

# find non-zero counts
nonZero <- (counts(sce) > 0)

# find rows/genes with at least one non-zero count 
keep <- (rowSums(nonZero) > 0)

table(keep)

# keep only the genes with non-zero count the SCE object 
sce_2 <- sce[keep, ]

# spike-ins ERCC
sce_2 <- splitAltExps(sce_2, ifelse(str_detect(rowData(sce_2)$gene_names, "Ercc"), "ERCC", "gene"))


            

################### Library size QC: filter out cells whose library size is less than threshold

# threshold was determined by library distribution 
threshold <- 20000

# tabulate the cells that were kept on the previous step
keep <- colSums(assay(sce)) > threshold 
table(keep)
sce <- sce[ , keep]


################### Batch QC: filter out outlier batch 


library(scater)

# calculate QCs (set a positive control gene to "ERCC")
GCmatrix <- perCellQCMetrics(sce, 
                             subsets = list(ERCC = str_detect(rowData(sce)$gene_names, 
                                                              "Ercc")))

colData(sce) <- GCmatrix  

# extract cell data into a data frame
cDataFrame <- as.data.frame(colData(sce)) %>%
        mutate(cell_name = rownames(colData(sce))) %>%
        separate(cell_name, 
                 c("Age", "Strain", "Lineage", "id"))

# plot cell data
GCByGroup <- ggplot(cDataFrame, aes(x = total, 
                                      y = subsets_ERCC_sum, 
                                      color = Lineage,
                                      shape = Age)) + 
        geom_point(alpha = 0.5, size = 2) + 
        theme_bw() +
        xlab("Total Counts per Cell") +
        ylab("ERCC Counts per Cell (Log)") +
        ggtitle("ERCC Counts over Total Counts by Age and Lineage") +
        scale_y_log10()






################### Gene QC: filter out genes 

# filter genes with at least a count of 1 in at least 2 cells
filter_genes <- apply(counts(sce), 
                      1, 
                      function(x){length(x[x >= 1]) >= 2})

# tabulate the results of filter_genes
sce <- sce[filter_genes, ]


#################################### Normalization ####################################


library(scran)

# Normalization w scran
sce <- computeSumFactors(sce)

# Distribution of Size factor
SizeFactorDist <- ggplot(data.frame(sf = sizeFactors(sce)), 
                         aes(x = sf)) + 
        geom_histogram(binwidth = 0.1,
                       fill = "pink",
                       color = "black") + 
        theme_bw() + 
        xlab("Size Factors") + 
        ylab("Frequency") +
        ggtitle("Distribution of Size Factors")


# view the assays for sce
assays(sce)

# normalize sce
norm_sce <- logNormCounts(sce)

# view the assays for normalized_sce
assays(norm_sce)


# examine the counts matrix of sce and the logcounts matrix of normalized_sce
counts(sce)[1:3, 1:3]
logcounts(norm_sce)[1:3, 1:3]



#################################### Dimensionality Reduction ####################################


CellType <- data.frame(name = colnames(norm_sce)) %>% 
        separate(name, c("age", "strain", "lineage", "id")) %>%
        mutate(age = as.factor(age),
               lineage = as.factor(lineage))

colData(norm_sce) <- cbind(colData(norm_sce), CellType[, c(1, 3)])
        

# row = cells, column = genes
TrLog <- t(logcounts(norm_sce))


# PCA
TrPCA <- prcomp(TrLog)


# tSNE 
library(Rtsne)

set.seed(21)
TrTSNE <- Rtsne(TrLog, PCA = F, perplexity = 50, max_iter = 2000)

plot(TrTSNE$itercosts, 
     type = "l", 
     main = "K-L Divergence Change",
     xlab = "Gradient Descent (Every 50 Steps)",
     ylab = "Total K-L Divergence Cost")

plot(TrTSNE$costs, 
     type = "l",
     main = "K-L Divergence Change",
     xlab = "Record Index",
     ylab = "Total K-L Divergence Cost per Record")

# Adding reduced dims to my sce object
reducedDims(norm_sce) <- list(PCA = TrPCA$x,
                              TSNE = TrTSNE$Y)


# data cleaning for plotting PCA and tSNE results 
DimRedTab <- cbind(as.data.table(reducedDim(norm_sce, "PCA")[, 1:2]),
                   as.data.table(reducedDim(norm_sce, "TSNE")[, 1:2]))

names(DimRedTab) <- c("PC1", "PC2", "tSNE_x", "tSNE_y")
DimRedTab <- DimRedTab[, c("Age", 
                           "Lineage") := 
                               .(colData(norm_sce)$age, 
                                 colData(norm_sce)$lineage)]

DimRedPlot_fn <- function(data, 
                          xcol, 
                          ycol, 
                          Lineage, 
                          Age,
                          title,
                          xtitle,
                          ytitle) {
        ggplot(data,
               aes(x = xcol, 
                   y = ycol, 
                   color = Lineage,
                   shape = Age)) + 
                geom_point(alpha = 0.5, size = 2) + 
                theme_bw() +
                ggtitle(title) + 
                xlab(xtitle) + 
                ylab(ytitle)
}



Total_PCA_plot <- DimRedPlot_fn(DimRedTab,
                                DimRedTab$PC1,
                                DimRedTab$PC2,
                                DimRedTab$Lineage,
                                DimRedTab$Age,
                                "PCA of Total Cell Population",
                                "PC1",
                                "PC2")


Total_tSNE_plot <- DimRedPlot_fn(DimRedTab,
                                 DimRedTab$tSNE_x,
                                 DimRedTab$tSNE_y,
                                 DimRedTab$Lineage,
                                 DimRedTab$Age,
                                 "t-SNE of Total Cell Population", 
                                 "tSNE_x",
                                 "tSNE_y")


#################################### Subsetting High Variance Genes ####################################

# Compute variance by gene
vars <- rowVars(logcounts(norm_sce))

#rename vars
names(vars) <- rownames(norm_sce)

vars <- sort(vars, decreasing = TRUE)

# subsetting sce
sub_sce <- norm_sce[names(vars[1:100])]


# row = cells, column = genes
s_TrLog <- t(logcounts(sub_sce))


# PCA
s_TrPCA <- prcomp(TrLog)


# tSNE 
set.seed(21)
s_TrTSNE <- Rtsne(s_TrLog, PCA = T, perplexity = 50, max_iter = 2000)

# Adding reduced dims to my sce object
reducedDims(sub_sce) <- list(PCA = s_TrPCA$x,
                             TSNE = s_TrTSNE$Y)


# data cleaning for plotting PCA and tSNE results 
sub_DimRedTab <- cbind(as.data.table(reducedDim(sub_sce, "PCA")[, 1:2]),
                   as.data.table(reducedDim(sub_sce, "TSNE")[, 1:2]))

names(sub_DimRedTab) <- c("PC1", "PC2", "tSNE_x", "tSNE_y")
sub_DimRedTab <- sub_DimRedTab[, c("Age", 
                           "Lineage") := 
                               .(colData(sub_sce)$age, 
                                 colData(sub_sce)$lineage)]

Total_SubPCA_plot <- DimRedPlot_fn(sub_DimRedTab,
                                   sub_DimRedTab$PC1,
                                   sub_DimRedTab$PC2,
                                   sub_DimRedTab$Lineage,
                                   sub_DimRedTab$Age,
                                   "PCA of Total Cell Population",
                                   "PC1",
                                   "PC2")



Total_SubtSNE_plot <- DimRedPlot_fn(sub_DimRedTab,
                                    sub_DimRedTab$tSNE_x,
                                    sub_DimRedTab$tSNE_y,
                                    sub_DimRedTab$Lineage,
                                    sub_DimRedTab$Age,
                                    "t-SNE of Total Cell Population", 
                                    "tSNE_x",
                                    "tSNE_y")


#################################### Clustering ####################################

# https://satijalab.org/seurat/essential_commands.html <-- seurat standard workflow 

library(Seurat)

colData(sce) <- cbind(colData(sce), CellType[, c("age", "lineage")])

#create an seurat object
seuset <- CreateSeuratObject(assay(sce),
                             meta.data = as.data.frame(colData(sce)))

seuset_sub <- CreateSeuratObject(assay(sub_sce),
                                 meta.data = as.data.frame(colData(sub_sce)))

# NormalizeData: counts are divided by total counts per cell, multiplied by the scale factor,
# and natural-log transformed using log1p 
seuset <- NormalizeData(seuset,
                        normalization.method = "LogNormalize")
seuset_sub <- NormalizeData(seuset_sub,
                            normalization.method = "LogNormalize")
        
seuset <- FindVariableFeatures(seuset)
seuset_sub <- FindVariableFeatures(seuset_sub)

seuset <- ScaleData(seuset)
seuset_sub <- ScaleData(seuset_sub)

# PCA, tSNE, and Clustering
seuset <- RunPCA(seuset)
seuset_sub <- RunPCA(seuset_sub)

seuset <- FindNeighbors(seuset)
seuset_sub <- FindNeighbors(seuset_sub)

seuset <- FindClusters(seuset)
seuset_sub <- FindClusters(seuset_sub)


seuset <- RunTSNE(seuset)
seuset_sub <- RunTSNE(seuset_sub)

# PCA and tSNE plots 
DimPlot(object = seuset, reduction = "pca")
DimPlot(object = seuset, reduction = "tsne")
DimPlot(object = seuset_sub, reduction = "pca")
DimPlot(object = seuset_sub, reduction = "tsne")


# A single gene in  dimensional reduction plot 
# CD48: MPP marker
# CD150 (Slamf1): LT-HSC marker 
FeaturePlot(seuset, features = "Cd48")
FeaturePlot(seuset, features = "Slamf1")
FeaturePlot(seuset, features = "Flt3")
FeaturePlot(seuset, features = "Cd34")

FeaturePlot(seuset_sub, features = "tSNE_1")


# Comparing specific genes in single cells 

FeatureScatter(object = seuset, feature1 = "Eef1a1", feature2 = "Cd48")
FeatureScatter(object = seuset, feature1 = "Eef1a1", feature2 = "Slamf1")
FeatureScatter(object = seuset, feature1 = "Eef1a1", feature2 = "Slamf1")
FeatureScatter(object = seuset, feature1 = "PC_1", feature2 = "Slamf1")


# Comparing cells across genes 
CellScatter(object = seuset, cell1 = "old_DBA_STHSC_30", cell2 = "old_DBA_MPP_14")

# The number of significant vs insignificant variables (genes)
# Black = insignificant, Red = significant 
VariableFeaturePlot(seuset)
VariableFeaturePlot(seuset_sub)

# Gene expression plot by cluster
VlnPlot(seuset, features = c("Cd48", "Slamf1", "Flt3", "Cd34"))
VlnPlot(seuset, features = "PC_1")
VlnPlot(seuset, features = "tSNE_1")
VlnPlot(seuset_sub, features = "PC_1")
VlnPlot(seuset_sub, features = "tSNE_1")


# Gene expression plot by lineage
VlnPlot(seuset, features = c("Cd48", "Slamf1", "Flt3", "Cd34"), split.by = "lineage")
VlnPlot(seuset, features = c("Cd48", "Slamf1", "Flt3", "Cd34"), split.by = "age")
VlnPlot(seuset, features = "tSNE_1", split.by = "lineage")
VlnPlot(seuset, features = "tSNE_1", split.by = "age")
VlnPlot(seuset, features = "PC_1", split.by = "lineage")
VlnPlot(seuset, features = "PC_1", split.by = "age")

# Gene expression plot by age
DotPlot(object = seuset, 
        features = c("Cd48", "Slamf1", "Flt3", "Cd34"), 
        split.by = "age")


RidgePlot(seuset, features = c("Cd48", "Slamf1", "Flt3", "Cd34"))




# Heatmaps 
DoHeatmap(seuset)
heatmap_mostchanged100 <- DoHeatmap(seuset_sub) + 
        ggtitle("100 Most Changed Genes") + 
        xlab("Cluster") + 
        ylab("Gene") + 
        theme(axis.text.y = element_text(size = 5),
              axis.text.x = element_blank(),
              axis.title = element_text(size = 15))

DimHeatmap(seuset, reduction = "pca", cells = 504)
DimHeatmap(seuset_sub, reduction = "pca", cells = 504)


# Distribution of gene across the low dimensions 
FeaturePlot(object = seuset, 
               features = c("Flt3", "Cd48"),
               blend = TRUE)


clusters <- seuset@active.ident
age <- seuset@meta.data[, "age"]
lineage <- seuset@meta.data[, "lineage"]

cluster_age_ratio <- as.data.table(table(clusters, age))[
        , Proportion := N/sum(N) * 100, by = "age"][
                , Proportion := round(Proportion, digits = 2)][]


cluster_lineage_ratio <- as.data.table(table(clusters, lineage))[
        , Proportion := N/sum(N) * 100, by = "lineage"][
                , Proportion := round(Proportion, digits = 2)][]

age_lineage_ratio <- as.data.table(table(age, lineage))[
        , Proportion := N/sum(N) * 100, by = "age"][
                , Proportion := round(Proportion, digits = 2)][
                        , .SD, by = "age"][]

cluster_celltype <- as.data.table(seuset@meta.data[, c("age", 
                                                       "lineage", 
                                                       "seurat_clusters")]) %>%
        unite(cell_type, age, lineage) 

cluster_celltype_dt <- round(100* prop.table(table(cluster_celltype$cell_type, 
                                                   cluster_celltype$seurat_clusters), 
                                             2), 
                             digits = 2) %>%
        as.data.table() %>%
        mutate(V1 = factor(V1, levels = c("young_LTHSC",
                                          "old_LTHSC",
                                          "young_STHSC",
                                          "old_STHSC",
                                          "young_MPP",
                                          "old_MPP")),
               V2 = factor(V2))

names(cluster_celltype_dt) <- c("Cell_Type", 
                                "Clusters",
                                "Proportion")
        
                            
ggplot(cluster_celltype_dt, 
       aes(x = Clusters,
           y = Cell_Type,
           fill = Proportion)) + 
        geom_tile() +
        theme_bw() + 
        ylab("Cell Type") +
        ggtitle("Population Proportion (%) of Each Cell Type per cluster")

library(pheatmap)

clusters_sub <- seuset_sub@active.ident
colData(sub_sce) <- cbind(colData(sub_sce), clusters_sub)
cmat <- logcounts(sub_sce)
meta <- as.data.frame(colData(sub_sce)[, c("age", "lineage", "clusters_sub")]) 
names(meta) <- c("Age", "Cell_Type", "Cluster")



pheatmap(cmat, 
         annotation = meta[, 1:2],
         show_rownames = F, 
         show_colnames = F,
         main = "Gene Expression Profile across Cells")