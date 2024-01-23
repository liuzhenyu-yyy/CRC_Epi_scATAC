#####################################################
#
# Re-analysis of Nature Genetics CRC 10X atlas
# Data source: NCBI Gene Expression Omnibus GSE132465
# Citation: https://www.nature.com/articles/s41588-020-0636-z
#
#####################################################

setwd("E:/LabWork/Project/CRC_NGS_ATAC/CRC_Epi_scATAC/Results/07.Independent_Validation")
load("Independent_Validation.RNA.RData")
source("../../code/00.Requirements.R")
library(Seurat)
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

## 2.2. Malignant vs Normal ----
table(epi.obj$Cell_subtype, epi.obj$Class)
marker.10X <- FindMarkers(
  epi.obj,
  group.by = "Class",
  ident.1 = "Tumor",
  ident.2 = "Normal",
  min.pct = 0.05,
  logfc.threshold = 0.01,
)
saveRDS(marker.10X, "NG_10X/marker.10X.RDS")
marker.10X$Gene <- rownames(marker.10X)

marker.10X$DEG <- "none"
marker.10X$DEG[marker.10X$avg_log2FC > 0.25 & marker.10X$p_val_adj < 0.05] <- "up"
marker.10X$DEG[marker.10X$avg_log2FC < -0.25 & marker.10X$p_val_adj < 0.05] <- "down"
table(marker.10X$DEG)

pdf("NG_10X/Volcano.10X.pdf", 3.5, 3)
ggplot(marker.10X, aes(x = avg_log2FC, y = -log10(p_val_adj), color = DEG)) +
    geom_point(size = 0.1) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed") +
    # ylim(c(0, 200)) +
    xlim(c(-4, 4)) +
    scale_color_manual(values = c("blue", "gray", "red")) +
    theme_classic() +
    xlab("Log2 fold change") +
    ylab("-log10 adjusted p-value")
dev.off()

marker.10X.list <- list(
    up = marker.10X$Gene[marker.10X$DEG == "up"],
    down = marker.10X$Gene[marker.10X$DEG == "down"]
)

table(epi.obj$Class, epi.obj$Patient)

## 2.3. marker my patient ----
epi.obj$iCMS <- "none"
epi.obj$iCMS[grep("CMS1", epi.obj$Cell_subtype)] <- "iCMS3"
epi.obj$iCMS[grep("CMS2", epi.obj$Cell_subtype)] <- "iCMS2"
epi.obj$iCMS[grep("CMS3", epi.obj$Cell_subtype)] <- "iCMS3"
epi.obj$iCMS[grep("CMS", epi.obj$Cell_subtype, invert = TRUE)] <- "Normal"
table(epi.obj$iCMS)

epi.obj$Patient_Tumor <- epi.obj$Patient
epi.obj$Patient_Tumor[epi.obj$Class == "Normal"] <- "Normal"
table(epi.obj$Patient_Tumor)

patients <- unique(epi.obj$Patient_Tumor) %>% setdiff("Normal")

marker.10x.patient.up <- list()
marker.10x.patient.down <- list()

for (one in patients) {
    message(Sys.time(), ": Find markers for ", one, "...")
    temp <- FindMarkers(
        epi.obj,
        group.by = "Patient_Tumor",
        ident.1 = one,
        ident.2 = "Normal"
    )
    marker.10x.patient.up[[one]] <- rownames(temp)[temp$avg_log2FC > 0.25 & temp$p_val_adj < 0.05]
    marker.10x.patient.down[[one]] <- rownames(temp)[temp$avg_log2FC < (-0.25) & temp$p_val_adj < 0.05]
}

length(marker.10x.patient.up %>% unlist() %>% unique())
table(marker.10X$DEG)

temp <- table(epi.obj$Patient_Tumor, epi.obj$iCMS)
temp2 <- colnames(temp)[apply(temp, 1, function(x) {
    which.max(x)
})]
names(temp2) <- rownames(temp)
temp2 <- sort(temp2)
temp2 <- temp2[c(24, 1:23)]

# patient heatmap
gene.selected <- c(
    unlist(marker.10x.patient.up[names(temp2[2:24])]),
    unlist(marker.10x.patient.down[names(temp2[2:24])])
) %>% unique()
length(unlist(marker.10x.patient.up) %>% unique())
length(unlist(marker.10x.patient.down) %>% unique())
# gene.selected <- gene.selected[order(marker.10X[gene.selected, "avg_log2FC"], decreasing = TRUE)]

plot.data <- AverageExpression(epi.obj,
    features = gene.selected,
    group.by = "Patient_Tumor"
)
plot.data <- plot.data$RNA %>% as.data.frame()
identical(rownames(plot.data), gene.selected)

anno.row <- data.frame(
    row.names = gene.selected,
    "Gene" = ifelse(gene.selected %in% unlist(marker.10x.patient.up), "up", "down")
)
table(anno.row$Gene)

anno.col <- data.frame(
    row.names = names(temp2),
    "iCMS" = temp2
)

plot.data <- plot.data[rownames(anno.row), rownames(anno.col)]
plot.data.nor <- apply(plot.data, 1, function(x) {
    x <- (x - mean(x)) / sd(x)
}) %>% t()
plot.data.nor[plot.data.nor > 2] <- 2
plot.data.nor[plot.data.nor < -2] <- -2
png("NG_10X/Heatmap.DEG.patient.png", 600, 600)
pheatmap::pheatmap(
    plot.data.nor,
    annotation_col = anno.col,
    annotation_row = anno.row,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    show_colnames = TRUE,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    annotation_colors = list(
        iCMS = c("Normal" = "#208a42", "iCMS2" = "#283891", "iCMS3" = "#62b7e6"),
        Gene = c("up" = "#cd2525", "down" = "#1774cd")
    )
)
dev.off()

## 2.4. known iCMS signature ----
markers.iCMS <- read.table("E:/LabWork/Project/CRC_NGS_ATAC/iCMS markers.txt",
    stringsAsFactors = FALSE,
    sep = "\t", header = TRUE
)
markers.iCMS <- base::as.list(markers.iCMS) %>% sapply(function(x) setdiff(x, ""))
lapply(markers.iCMS, length)
(308 + 74) / 5160 # 7.4%
(279 + 54) / 1468 # 22.7%

markers.iCMS <- lapply(markers.iCMS, function(x) intersect(x, rownames(epi.obj)))
lapply(markers.iCMS, length)

gene.selected <- unlist(markers.iCMS)

plot.data <- AverageExpression(epi.obj,
    features = gene.selected,
    group.by = "Patient_Tumor"
)
plot.data <- plot.data$RNA %>% as.data.frame()

anno.row <- data.frame(
    row.names = gene.selected,
    "Gene" = names(gene.selected) %>%
        gsub("[0-9]$", "", .) %>%
        gsub("[0-9]$", "", .) %>%
        gsub("[0-9]$", "", .)
)
table(anno.row$Gene)

anno.col <- data.frame(
    row.names = names(temp2),
    "iCMS" = temp2
)

plot.data <- plot.data[rownames(anno.row), rownames(anno.col)]
plot.data.nor <- apply(plot.data, 1, function(x) {
    x <- (x - mean(x)) / sd(x)
}) %>% t()
plot.data.nor[plot.data.nor > 2] <- 2
plot.data.nor[plot.data.nor < -2] <- -2
png("NG_10X/Heatmap.iCMS.geme.patient.png", 600, 600)
pheatmap::pheatmap(
    plot.data.nor,
    annotation_col = anno.col,
    annotation_row = anno.row,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    show_colnames = TRUE,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    annotation_colors = list(
        iCMS = c("Normal" = "#208a42", "iCMS2" = "#283891", "iCMS3" = "#62b7e6"),
        Gene = c("iCMS2_Up" = "#c80505", "iCMS2_Down" = "#064e91",  "iCMS3_Up" = "#DC7C7C", "iCMS3_Down" = "#66aced")
    ),
    gaps_row = sum(length(markers.iCMS$iCMS2_Up), length(markers.iCMS$iCMS2_Down)),
    gaps_col = c(1, 15)
)
dev.off()

length(markers.iCMS$iCMS2_Up) / length(unlist(marker.10x.patient.up) %>% unique())
rm(plot.data, plot.data.nor, anno.row, anno.col, patients, one, temp, temp2)

## 2.5 expression of selected TFs ----
# iCMS TFs
TF.selected <- c(
    "NFE2L2", "MAFB", "MAFK", "FOXA3", "FOXA2", "FOXM1", "SOX2", "SOX4", "ELF1", "EHF",
    "ETS1", "JUN", "LEF1", "TCF3", "NR4A1", "HNF1A", "CDX2", "PPARA", "NR2C2", "HNF4A"
)

epi.obj$iCMS
p <- Seurat::DotPlot(epi.obj, features = TF.selected, group.by = "iCMS", scale = FALSE)
plot.data <- p$data

plot.data.1 <- plot.data[plot.data$id == "iCMS2", ] %>%
    .[TF.selected, ] %>%
    mutate(id = "iCMS2")
plot.data.2 <- plot.data[plot.data$id == "iCMS3", ] %>%
    .[TF.selected, ] %>%
    mutate(id = "iCMS3")
plot.data <- rbind(plot.data.1, plot.data.2)
plot.data$features.plot <- factor(c(TF.selected, TF.selected), levels = TF.selected)
plot.data[plot.data$avg.exp > 4, ]$avg.exp <- 4

pdf("NG_10X/Dot.motif.sig.iCMS.selected.pdf", 2.5, 4)
ggplot(plot.data, aes(x = id, y = features.plot, fill = avg.exp, size = pct.exp)) +
    geom_point(pch = 21) +
    theme_bw() +
    xlab("iCMS") +
    scale_fill_gradientn(colors = ArchR::paletteContinuous()) +
    scale_size_continuous(limits = c(0, 90)) +
    ylab("TF")
dev.off()

# CIMP TFs
TF.selected <- c(
    "TCF7L2", "TCF7", "TCF3", "LEF1", "CREB5", "CDX4", "CDX2", "ATF2", "TEAD2", "TEAD1"
)

epi.obj$Epi_type <- ifelse(epi.obj$Patient_Tumor == "Normal", "Normal", "Malignant")
table(epi.obj$Epi_type)

p <- Seurat::DotPlot(epi.obj, features = TF.selected, group.by = "Epi_type", scale = FALSE)
plot.data <- p$data

plot.data <- plot.data[plot.data$id == "Malignant", ] %>%
    .[TF.selected, ]

plot.data$features.plot <- factor(c(TF.selected), levels = TF.selected)
plot.data[plot.data$avg.exp > 4, ]$avg.exp <- 4

pdf("NG_10X/Dot.motif.sig.CIMP.selected.pdf", 5, 1)
ggplot(plot.data, aes(x = features.plot, y = id, fill = avg.exp, size = pct.exp)) +
    geom_point(pch = 21) +
    theme_bw() +
    xlab("iCMS") +
    scale_fill_gradientn(colors = ArchR::paletteContinuous()) +
    scale_size_continuous(limits = c(0, 90)) +
    ylab("TF")
dev.off()

# 3. Diff gene ----
## 3.1. TCGA RNA data ----
marker.TCGA <- read.table("data/TCGA/TCGA.GEPIA.Limma.txt", sep = "\t", header = TRUE)
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
    scale_color_manual(values = c("blue", "gray", "red")) +
    theme_classic() +
    xlab("Log2 fold change") +
    ylab("-log10 adjusted p-value")
dev.off()

marker.TCGA.list <- list(
    up = marker.TCGA$gene[marker.TCGA$DEG == "up"],
    down = marker.TCGA$gene[marker.TCGA$DEG == "down"]
)

## 3.2. ATAC data ----
marker.ATAC <- getMarkerFeatures(
    ArchRProj = proj_Epi,
    useMatrix = "GeneScoreMatrix",
    groupBy = "Epi_type",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = c("Malignant"),
    bgdGroups = "Normal"
)

marker.ATAC.df <- getMarkers(marker.ATAC,
    cutOff = "FDR <= 0.05 & abs(Log2FC) >= 0.25"
)
marker.ATAC.df <- marker.ATAC.df$Malignant %>% as.data.frame()
dim(marker.ATAC.df)

marker.ATAC.df$DEG <- "none"
marker.ATAC.df[marker.ATAC.df$Log2FC > 0, "DEG"] <- "up"
marker.ATAC.df[marker.ATAC.df$Log2FC < 0, "DEG"] <- "down"
table(marker.ATAC.df$DEG)

DEG.ATAC.list <- list(
    up = marker.ATAC.df$name[marker.ATAC.df$DEG == "up"],
    down = marker.ATAC.df$name[marker.ATAC.df$DEG == "down"]
)

## 3.3. compare DEGs ----
genes.ATAC <- ArchR::getFeatures(proj_Epi) %>% as.character()
genes.TCGA <- rownames(marker.TCGA)
genes.10X <- rownames(epi.obj)
gene.common <- intersect(genes.ATAC, genes.TCGA) %>%
    intersect(genes.10X)

v <- Vennerable::Venn(list(
    "scATAC" = DEG.ATAC.list$up %>% intersect(gene.common),
    "TCGA" = marker.TCGA.list$up %>% intersect(gene.common),
    "NG_10X_scRNA" = marker.10X.list$up %>% intersect(gene.common)
))
intersect(marker.TCGA.list$up, DEG.ATAC.list$up) %>% sort()

pdf("Venn.ATAC.merge.Up.pdf", 10, 10)
plot(v, doWeights = TRUE, doEuler = TRUE, show = list(Faces = FALSE))
dev.off()

v <- Vennerable::Venn(list(
    "scATAC" = DEG.ATAC.list$up %>% intersect(gene.common),
    "TCGA" = marker.TCGA.list$up %>% intersect(gene.common),
    "NG_10X_scRNA" = unlist(marker.10x.patient.up) %>%
        unique() %>%
        intersect(gene.common)
))
pdf("Venn.ATAC.patient.Up.pdf", 10, 10)
plot(v, doWeights = TRUE, doEuler = TRUE, show = list(Faces = FALSE))
dev.off()

v <- Vennerable::Venn(list(
    "scATAC" = DEG.ATAC.list$down %>% intersect(gene.common),
    "TCGA" = marker.TCGA.list$down %>% intersect(gene.common),
    "NG_10X_scRNA" = marker.10X.list$down %>% intersect(gene.common)
))
pdf("Venn.ATAC.merge.Down.pdf", 10, 10)
plot(v, doWeights = TRUE, doEuler = TRUE, show = list(Faces = FALSE))
dev.off()

v <- Vennerable::Venn(list(
    "scATAC" = DEG.ATAC.list$down %>% intersect(gene.common),
    "TCGA" = marker.TCGA.list$down %>% intersect(gene.common),
    "NG_10X_scRNA" = unlist(marker.10x.patient.down) %>%
        unique() %>%
        intersect(gene.common)
))
pdf("Venn.ATAC.patient.Down.pdf", 10, 10)
plot(v, doWeights = TRUE, doEuler = TRUE, show = list(Faces = FALSE))
dev.off()

rm(p, plot.data, foo, v)

rm(proj_Epi)
save.image("Independent_Validation.RNA.RData")
