data1 <- Read10X(data.dir = "../data/SCC0120_1_Sample_2/outs/filtered_feature_bc_matrix/")

sample1 <- CreateSeuratObject(counts = data1$`Gene Expression`)
sample1[['HTO']] <- CreateAssayObject(counts = data1$`Antibody Capture`[c(1:8),])
sample1[['CITE']] <- CreateAssayObject(counts = data1$`Antibody Capture`[c(9:10),])
sample1[['PROT']] <- CreateAssayObject(counts = data1$`Antibody Capture`)
sample1[["percent.mt"]] <- PercentageFeatureSet(sample1, pattern = "^MT-")
sample1 <- subset(sample1, subset = nFeature_RNA > 1000 & percent.mt < 10 & percent.mt > 0.5 & nCount_RNA < 45000)
sample1@meta.data$sample <- "SCC0120_1_Sample_2"
sample1 <- Seurat::NormalizeData(sample1, assay = "HTO", normalization.method = "CLR")
sample1 <- Seurat::NormalizeData(sample1, assay = "CITE", normalization.method = "CLR")
DefaultAssay(sample1) <- "HTO"
rownames(sample1@assays$HTO)
RidgePlot(object = sample1, rownames(sample1@assays$HTO)[1:8], group.by = , ncol = 4)
sample1 <- Seurat::HTODemux(sample1, assay = "HTO", positive.quantile = 0.99)
RidgePlot(object = sample1, rownames(sample1@assays$HTO)[1:8], group.by = "hash.ID", ncol = 4)
table(sample1$HTO_classification.global)
sample1@meta.data
table(sample1$hash.ID)
Seurat::Idents(sample1) <- "hash.ID"
#sample1 <- subset(x = sample1, idents = c("Doublet"), invert = TRUE)
sample1$hashtag<-sample1$hash.ID
sample1@meta.data$unique<-sample1@meta.data$hashtag
sample1@meta.data$unique <- gsub("Hashtag1-TotalA", "hs_1", sample1@meta.data$unique)
sample1@meta.data$unique <- gsub("Hashtag4-TotalA", "hs_2", sample1@meta.data$unique)
sample1@meta.data$unique <- gsub("Hashtag7-TotalA", "hs_3", sample1@meta.data$unique)
sample1@meta.data$unique <- gsub("Hashtag2-TotalA", "is_1", sample1@meta.data$unique)
sample1@meta.data$unique <- gsub("Hashtag5-TotalA", "is_2", sample1@meta.data$unique)
sample1@meta.data$unique <- gsub("Hashtag8-TotalA", "is_3", sample1@meta.data$unique)
sample1@meta.data$unique <- gsub("Hashtag3-TotalA", "pbmc_1", sample1@meta.data$unique)
sample1@meta.data$unique <- gsub("Hashtag6-TotalA", "pbmc_2", sample1@meta.data$unique)
sample1@meta.data$group<-sample1@meta.data$hashtag
sample1@meta.data$group <- gsub("Hashtag1-TotalA", "hs", sample1@meta.data$group)
sample1@meta.data$group <- gsub("Hashtag4-TotalA", "hs", sample1@meta.data$group)
sample1@meta.data$group <- gsub("Hashtag7-TotalA", "hs", sample1@meta.data$group)
sample1@meta.data$group <- gsub("Hashtag2-TotalA", "is", sample1@meta.data$group)
sample1@meta.data$group <- gsub("Hashtag5-TotalA", "is", sample1@meta.data$group)
sample1@meta.data$group <- gsub("Hashtag8-TotalA", "is", sample1@meta.data$group)
sample1@meta.data$group <- gsub("Hashtag3-TotalA", "pbmc", sample1@meta.data$group)
sample1@meta.data$group <- gsub("Hashtag6-TotalA", "pbmc", sample1@meta.data$group)
