setwd("E:/LabWork/Project/CRC_NGS_ATAC/CRC_Epi_scATAC/Results/07.Independent_Validation")

source("../../code/00.Requirements.R")
proj_Epi <- loadArchRProject(project.dir.epi, force = TRUE)


# 1. Primary analysis of CRC 10X atlas -----
dir.create("NG_10X")

## 1.1 load data ----
rc <- as.data.frame(fread("data/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt"),
  header = TRUE
)
sample.info.10x <- as.data.frame(fread("data/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt"))

rownames(rc) <- rc$Index
rc <- rc[, -1]

identical(colnames(rc.SM), sample.info.10x$Index)

rownames(sample.info.10x) <- sample.info.10x$Index

rc.SM_1 <- as(as.matrix(rc[, 1:30000]), "dgCMatrix")
rc.SM_2 <- as(as.matrix(rc[, 30001:63689]), "dgCMatrix")

rc.SM <- cbind(rc.SM_1, rc.SM_2)
dim(rc.SM)
rm(rc, rc.SM_1, rc.SM_2)
saveRDS(rc.SM, "data/rc.SM.RDS")

# patient info
patient.info <- read.table("patient.info.tsv", sep = "\t", header = TRUE)
rownames(patient.info) <- patient.info$Patient

sample.info.10x$Gender <- patient.info[sample.info.10x$Patient, ]$Gender
sample.info.10x$Stage <- patient.info[sample.info.10x$Patient, ]$Stage
sample.info.10x$Region <- patient.info[sample.info.10x$Patient, ]$Anatomic.region
sample.info.10x$Side <- patient.info[sample.info.10x$Patient, ]$Left.Right.sided
sample.info.10x$MSI <- patient.info[sample.info.10x$Patient, ]$MSI

## 1.2. Suerat Pipeline ----
CRC.obj <- CreateSeuratObject(rc.SM, project = "CRC", meta.data = sample.info.10x)
rm(rc.SM)

CRC.obj <- NormalizeData(CRC.obj) %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  ScaleData() %>%
  RunPCA()

ElbowPlot(CRC.obj)
pc.use <- 1:10

CRC.obj <- FindNeighbors(CRC.obj)

CRC.obj <- RunTSNE(CRC.obj, dims = pc.use)

pdf("NG_10X/tSNE.Cell.type.all.pdf", 6, 5)
DimPlot(CRC.obj,
  group.by = "Cell_type",
  cols = brewer.pal(6, "Set2")
) +
  theme_bw() + theme(panel.grid = element_blank()) +
  coord_fixed()
dev.off()

pdf("NG_10X/tSNE.Cell.subtype.all.pdf", 7, 6)
DimPlot(CRC.obj, group.by = "Cell_subtype") +
  theme_bw() + theme(panel.grid = element_blank()) +
  coord_fixed()
dev.off()

pdf("NG_10X/tSNE.sample.info.10x.all.pdf", 6, 5)
DimPlot(CRC.obj, group.by = "Patient") +
  theme_bw() + theme(panel.grid = element_blank()) +
  coord_fixed()
DimPlot(CRC.obj, group.by = "Class") +
  theme_bw() + theme(panel.grid = element_blank()) +
  coord_fixed()
DimPlot(CRC.obj, group.by = "Sample") +
  theme_bw() + theme(panel.grid = element_blank()) +
  coord_fixed()
DimPlot(CRC.obj, group.by = "Stage") +
  theme_bw() + theme(panel.grid = element_blank()) +
  coord_fixed()
DimPlot(CRC.obj, group.by = "Gender") +
  theme_bw() + theme(panel.grid = element_blank()) +
  coord_fixed()
DimPlot(CRC.obj, group.by = "Region") +
  theme_bw() + theme(panel.grid = element_blank()) +
  coord_fixed()
DimPlot(CRC.obj, group.by = "Side") +
  theme_bw() + theme(panel.grid = element_blank()) +
  coord_fixed()
DimPlot(CRC.obj, group.by = "MSI") +
  theme_bw() + theme(panel.grid = element_blank()) +
  coord_fixed()
dev.off()

# 2. subset epithelial cells ----
## 2.1. primary analysis ----
CRC.obj <- readRDS("E:/LabWork/Project/CRC_NGS_ATAC/new_analysis/10X/epithelial/CRC.obj.RDS")
Idents(CRC.obj) <- CRC.obj$Cell_subtype
table(CRC.obj$Cell_subtype, CRC.obj$Cell_type)

epi.obj <- subset(CRC.obj, subset = Cell_type == "Epithelial cells")
epi.obj <- subset(epi.obj, subset = nFeature_RNA > 500)
rm(CRC.obj)
quantile(epi.obj$nFeature_RNA)

epi.obj <- FindVariableFeatures(epi.obj, selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:15)

pdf("NG_10X/UMAP.epi.Cell.subtype.pdf", 7, 6)
DimPlot(epi.obj, group.by = "Cell_subtype") +
    theme_classic() +
    scale_color_manual(values = ArchR::paletteDiscrete(epi.obj$Cell_subtype)) +
    coord_fixed()
dev.off()

pdf("NG_10X/UMAP.epi.Class.pdf", 7, 6)
DimPlot(epi.obj, group.by = "Class") +
    theme_classic() +
    scale_color_manual(values = ArchR::paletteDiscrete(epi.obj$Class)) +
    coord_fixed()
dev.off()

pdf("NG_10X/UMAP.epi.Sample.pdf", 7, 6)
DimPlot(epi.obj, group.by = "Sample") +
    theme_classic() +
    scale_color_manual(values = ArchR::paletteDiscrete(epi.obj$Sample)) +
    coord_fixed()
dev.off()

pdf("NG_10X/UMAP.epi.Patient.pdf", 7, 6)
DimPlot(epi.obj, group.by = "Patient") +
    theme_classic() +
    scale_color_manual(values = ArchR::paletteDiscrete(epi.obj$Patient)) +
    coord_fixed()
dev.off()

## 2.1. Malignant vs Normal ----
table(epi.obj$Cell_subtype, epi.obj$Class)
marker.10X <- FindMarkers(
  epi.obj,
  group.by = "Class",
  ident.1 = "Tumor",
  ident.2 = "Normal",
  min.pct = 0.05,
  logfc.threshold = 0.01,
)

# 3. Diff gene ----
## 3.1. TCGA RNA data ----
marker.TCGA <- read.table("TCGA/TCGA.GEPIA.Limma.txt", sep = "\t", header = TRUE)
colnames(marker.TCGA) <- c("gene", "gene_id", "mean_tumor", "mean_normal", "log2fc", "p_adj")
marker.TCGA <- marker.TCGA[!duplicated(marker.TCGA$gene), ]
rownames(marker.TCGA) <- marker.TCGA$gene

marker.TCGA$DEG <- "none"
marker.TCGA$DEG[marker.TCGA$log2fc > 0.25 & marker.TCGA$p_adj < 0.05] <- "up"
marker.TCGA$DEG[marker.TCGA$log2fc < -0.25 & marker.TCGA$p_adj < 0.05] <- "down"
table(marker.TCGA$DEG)

pdf("Volcano.TCGA.pdf", 3.5, 3)
ggplot(marker.TCGA, aes(x = log2fc, y = -log10(p_adj), color = DEG)) +
    geom_point(size = 0.1) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed") +
    ylim(c(0, 200)) +
    scale_color_manual(values = c("red", "gray", "blue")) +
    theme_classic() +
    xlab("Log2 fold change") +
    ylab("-log10 adjusted p-value")
dev.off()

marker.TCGA.list <- list(
    up = marker.TCGA$gene[marker.TCGA$DEG == "up"],
    down = marker.TCGA$gene[marker.TCGA$DEG == "down"]
)

## 3.2. 10X RNA-seq data ----
marker.10X <- readRDS("10X/marker.10X.RDS")
marker.10X$DEG <- "none"
marker.10X$DEG[marker.10X$avg_log2FC > 0.25 & marker.10X$p_val_adj < 0.05] <- "up"
marker.10X$DEG[marker.10X$avg_log2FC < -0.25 & marker.10X$p_val_adj < 0.05] <- "down"
table(marker.10X$DEG)

pdf("Volcano.10X.pdf", 3.5, 3)
ggplot(marker.10X, aes(x = avg_log2FC, y = -log10(p_val_adj), color = DEG)) +
    geom_point(size = 0.1) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed") +
    # ylim(c(0, 200)) +
    xlim(c(-4, 4)) +
    scale_color_manual(values = c("red", "gray", "blue")) +
    theme_classic() +
    xlab("Log2 fold change") +
    ylab("-log10 adjusted p-value")
dev.off()

marker.10X.list <- list(
    up = marker.10X$Gene[marker.10X$DEG == "up"],
    down = marker.10X$Gene[marker.10X$DEG == "down"]
)

## 3.3. ATAC data ----
DEG.ATAC <- read.table("DEG.ATAC.txt", sep = "\t", header = TRUE)
DEG.ATAC.list <- list(
    up = c(DEG.ATAC$C_high, DEG.ATAC$NM_low) %>% setdiff(""),
    down = c(DEG.ATAC$NM_high, DEG.ATAC$C_low) %>% setdiff("")
)
sapply(DEG.ATAC.list, length)

table(proj_Epi$Epi_type)
marker.ATAC <- getMarkerFeatures(
    ArchRProj = proj_Epi,
    useMatrix = "GeneScoreMatrix",
    groupBy = "Epi_type",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = c("Malignant"),
    bgdGroups = "Normal"
)

marker.ATAC.df <- do.call(cbind, marker.ATAC@assays@data)
rownames(marker.ATAC.df) <- marker.ATAC@elementMetadata$name
colnames(marker.ATAC.df) <- names(marker.ATAC@assays@data)

## 3.4. compare DEGs ----
genes.ATAC <- ArchR::getFeatures(proj_Epi) %>% as.character()
genes.TCGA <- rownames(marker.TCGA)
genes.10X <- rownames(marker.10X)
gene.common <- intersect(genes.ATAC, genes.TCGA) %>%
    intersect(genes.10X)

v <- Vennerable::Venn(list(
    "ATAC" = DEG.ATAC.list$up %>% intersect(gene.common),
    "TCGA" = marker.TCGA.list$up %>% intersect(gene.common),
    "10X" = marker.10X.list$up %>% intersect(gene.common)
))
intersect(marker.TCGA.list$up, DEG.ATAC.list$up) %>% sort()

pdf("Venn.ATAC.Up.pdf", 10, 10)
plot(v, doWeights = TRUE, doEuler = TRUE, show = list(Faces = FALSE))
dev.off()

v <- Vennerable::Venn(list(
    "ATAC" = DEG.ATAC.list$down %>% intersect(gene.common),
    "TCGA" = marker.TCGA.list$down %>% intersect(gene.common),
    "10X" = marker.10X.list$down %>% intersect(gene.common)
))
pdf("Venn.ATAC.Down.pdf", 10, 10)
plot(v, doWeights = TRUE, doEuler = TRUE, show = list(Faces = FALSE))
dev.off()
rm(p, plot.data, foo, v)

rm(proj_Epi)
save.image("Independent_Validation.RNA.RData")
