setwd("E:/LabWork/Project/CRC_NGS_ATAC/CRC_Epi_scATAC/Results/07.Independent_Validation")

source("../../code/00.Requirements.R")
proj_Epi <- loadArchRProject(project.dir.epi, force = TRUE)

# motif match matrix
motif.match <- getMatches(ArchRProj = proj_Epi, name = "Motif")
temp <- rowRanges(motif.match)
motif.match <- motif.match@assays@data$matches

rownames(motif.match) <- paste(seqnames(temp), start(temp), end(temp), sep = "_")
colnames(motif.match) <- gsub("_.+?$", "", colnames(motif.match))
motif.match[1:5, 1:5]
dim(motif.match)
temp <- colSums(motif.match)
temp[c("HNF4A", "PPARA", "TCF3", "LEF1", "FOXA3", "MAFK")]

# 1. re-analysis of NG CRC samples ----
dir.create("NG_ATAC")
## 1.1. create ArchR project ----
# nohup /usr/bin/Rscript Create_Arrow.R hg38 GSM6058842_CRC-1-8810-D_20200917_fragments.tsv.gz > CRC-1.log 2>&1 &

proj_NG <- ArchRProject(
    ArrowFiles = dir("arrows", full.names = TRUE, pattern = ".arrow$"),
    outputDirectory = "./Project_Dir_NG",
    copyArrows = TRUE
)
colnames(proj_NG@cellColData)
quantile(proj_NG$TSSEnrichment)
quantile(proj_NG$nFrags)

proj_NG <- saveArchRProject(proj_NG, load = TRUE)
proj_NG <- loadArchRProject("./Project_Dir_NG")

## 1.2. quality control ----
# QC
cell.selected <- proj_NG@cellColData %>%
    as.data.frame() %>%
    filter(
        nFrags > 3000,
        TSSEnrichment > 4
    ) %>%
    rownames()
length(cell.selected)
proj_NG <- proj_NG[cell.selected, ]

# add patient metadata
table(proj_NG$Sample)
proj_NG$Sample_2 <- gsub("GSM.*_|-D", "", proj_NG$Sample)
table(proj_NG$Sample_2)
sample.info.NG <- proj_NG@cellColData %>% as.data.frame()

patient.info.NG <- read.table("patient.info.NG.txt", sep = "\t", header = TRUE)
rownames(patient.info.NG) <- patient.info.NG$Sample
table(sample.info.NG$Sample_2 %in% rownames(patient.info.NG))

sample.info.NG$Patient <- patient.info.NG[sample.info.NG$Sample_2, "Donor"]
sample.info.NG$FAP <- patient.info.NG[sample.info.NG$Sample_2, "FAP"]
sample.info.NG$GrossPathology <- patient.info.NG[sample.info.NG$Sample_2, "GrossPathology"]

proj_NG$Patient <- sample.info.NG$Patient
proj_NG$FAP <- sample.info.NG$FAP
proj_NG$GrossPathology <- sample.info.NG$GrossPathology

## 1.3 add peak matrx, reduction and clustering ----
# add peak matrix
proj_NG <- addPeakSet(
    ArchRProj = proj_NG,
    peakSet = proj_Epi@peakSet,
    force = TRUE
)

proj_NG <- addPeakMatrix(
    ArchRProj = proj_NG,
    binarize = FALSE,
    verbose = TRUE,
)
proj_NG <- saveArchRProject(proj_NG, load = TRUE)
proj_NG <- loadArchRProject("./Project_Dir_NG")

ArchR::getAvailableMatrices(proj_NG)

# reduction
proj_NG <- addIterativeLSI(
    ArchRProj = proj_NG,
    useMatrix = "PeakMatrix",
    name = "IterativeLSI"
)

proj_NG <- addUMAP(
    ArchRProj = proj_NG,
    force = TRUE,
    reducedDims = "IterativeLSI",
    name = "UMAP",
    nNeighbors = 35,
    minDist = 0.5,
    metric = "cosine"
)

# clustering
proj_NG <- addClusters(
    input = proj_NG,
    reducedDims = "IterativeLSI",
    dimsToUse = 1:30,
    # maxClusters=50,
    resolution = 0.8,
    force = TRUE
)
table(proj_NG$Clusters)

p <- plotEmbedding(
    ArchRProj = proj_NG, colorBy = "cellColData",
    name = "Sample_2", embedding = "UMAP",
    size = 0.2, plotAs = "points",
    labelMeans = FALSE
)
pdf("NG_ATAC/UMAP.Sample2.pdf", 7, 7)
plot(p)
dev.off()

p <- plotEmbedding(
    ArchRProj = proj_NG, colorBy = "cellColData",
    name = "Patient", embedding = "UMAP",
    size = 0.2, plotAs = "points",
    labelMeans = FALSE
)
pdf("NG_ATAC/UMAP.Patient.pdf", 7, 7)
plot(p)
dev.off()

p <- plotEmbedding(
    ArchRProj = proj_NG, colorBy = "cellColData",
    name = "Clusters", embedding = "UMAP",
    size = 0.2, plotAs = "points",
    labelMeans = TRUE
)
pdf("NG_ATAC/UMAP.Clusters.pdf", 7, 7)
plot(p)
dev.off()

p <- plotEmbedding(
    ArchRProj = proj_NG, colorBy = "cellColData",
    name = "Patient", embedding = "UMAP",
    size = 0.2, plotAs = "points",
    labelMeans = FALSE
)
pdf("NG_ATAC/UMAP.Patient.pdf", 7, 7)
plot(p)
dev.off()

table(proj_NG$GrossPathology)
p <- plotEmbedding(
    ArchRProj = proj_NG, colorBy = "cellColData",
    name = "GrossPathology", embedding = "UMAP",
    size = 0.2, plotAs = "points",
    pal = c(Normal = "#1f8942", Adenocarcinoma = "#272d6a", Polyp = "#d41f25"),
    labelMeans = FALSE
)
pdf("NG_ATAC/UMAP.GrossPathology.pdf", 7, 7)
plot(p)
dev.off()

p <- plotEmbedding(
    ArchRProj = proj_NG, colorBy = "cellColData",
    name = "FAP", embedding = "UMAP",
    size = 0.2, plotAs = "points",
    labelMeans = FALSE
)
pdf("NG_ATAC/UMAP.FAP.pdf", 7, 7)
plot(p)
dev.off()

## 1.4. cell typing ----
# marker gene 
markers <- c(
    "PTPRC", "CD3D", "CD3G", "CD79A", "MS4A1", "CD14", "CD68",
    "EPCAM", "PECAM1", "CDH5", "VIM", "COL2A1"
)
proj_NG <- addImputeWeights(proj_NG)

p <- plotEmbedding(
    ArchRProj = proj_NG, colorBy = "GeneScoreMatrix",
    name = markers, embedding = "UMAP",
    imputeWeights = getImputeWeights(proj_NG),
    size = 0.2, plotAs = "points"
)
pdf("NG_ATAC/UMAP.Marker.Gene.pdf", 12, 12)
patchwork::wrap_plots(p, ncol = 4)
dev.off()

# gene module ----
markers <- c(
    "NXPE4", "TNP2", "HAGLR", "HOXA2",
    "ATP2A3", "MIR3185", "CDH7",
    "PLA2G2E", "TFAP2A", "CLDN1",
    "LOC101928540", "MIR4463"
)
proj_NG <- addImputeWeights(proj_NG)

p <- plotEmbedding(
    ArchRProj = proj_NG, colorBy = "GeneScoreMatrix",
    name = markers, embedding = "UMAP",
    imputeWeights = getImputeWeights(proj_NG),
    size = 0.2, plotAs = "points"
)
pdf("NG_ATAC/UMAP.Marker.AD.Gene.pdf", 12, 12)
patchwork::wrap_plots(p, ncol = 4)
dev.off()

p <- plotEmbedding(
    ArchRProj = proj_NG, colorBy = "GeneScoreMatrix",
    name = markers, embedding = "UMAP",
    imputeWeights = getImputeWeights(proj_NG),
    size = 0.2, plotAs = "points"
)
pdf("NG_ATAC/UMAP.Marker.Gene.pdf", 12, 12)
patchwork::wrap_plots(p, ncol = 4)
dev.off()

plot.data <- table(proj_NG$Clusters, proj_NG$GrossPathology) %>%
    as.data.frame()

pdf("NG_ATAC/Bar.Cluster.GrossPathology.pdf", 10, 3)
ggplot(plot.data, aes(x = Var1, y = Freq, fill = Var2)) +
    geom_bar(stat = "identity", position = "fill") +
    # geom_text(aes(label = Freq), position = position_fill(vjust = 0.5)) +
    theme_classic() +
    xlab("Clusters") +
    ylab("Proportion") +
    scale_fill_manual(values = c(Normal = "#1f8942", Adenocarcinoma = "#272d6a", Polyp = "#d41f25"))
dev.off()

plot.data <- table(proj_NG$Clusters, proj_NG$Patient) %>%
    as.data.frame()
pdf("NG_ATAC/Bar.Cluster.Paatient.pdf", 10, 3)
ggplot(plot.data, aes(x = Var1, y = Freq, fill = Var2)) +
    geom_bar(stat = "identity", position = "fill") +
    theme_classic() +
    xlab("Clusters") +
    ylab("Proportion") +
    scale_fill_manual(values = ArchR::paletteDiscrete(proj_NG$Patient))
dev.off()

sample.info.NG  <- proj_NG@cellColData %>% as.data.frame()
sample.info.NG$Cell_type_1 <- "none"
sample.info.NG[sample.info.NG$Clusters %in% paste0("C", c(6)), "Cell_type_1"] <- "T"
sample.info.NG[sample.info.NG$Clusters %in% paste0("C", c(7, 8)), "Cell_type_1"] <- "B"
sample.info.NG[sample.info.NG$Clusters %in% paste0("C", c(4, 5)), "Cell_type_1"] <- "Myeloid"
sample.info.NG[sample.info.NG$Clusters %in% paste0("C", c(1:3)), "Cell_type_1"] <- "Stromal"
sample.info.NG[sample.info.NG$Clusters %in% paste0("C", c(9:12, 14:16)), "Cell_type_1"] <- "Malignant"
sample.info.NG[sample.info.NG$Clusters %in% paste0("C", c(13, 17:19, 21:23)), "Cell_type_1"] <- "Normal"
sample.info.NG[sample.info.NG$Clusters %in% paste0("C", c(25, 24, 20)), "Cell_type_1"] <- "Adenoma"
# sample.info.NG[sample.info.NG$Clusters %in% paste0("C", c(17, 18, 19, 21, 22)), "Cell_type_1"] <- "Intermediate"
table(sample.info.NG$Cell_type_1)

proj_NG$Cell_type_1 <- sample.info.NG$Cell_type_1

p <- plotEmbedding(
    ArchRProj = proj_NG, colorBy = "cellColData",
    name = "Cell_type_1", embedding = "UMAP",
    size = 0.2, plotAs = "points",
    labelMeans = FALSE
)
pdf("NG_ATAC/UMAP.Cell_type_1.pdf", 7, 7)
plot(p)
dev.off()

# 2. diff peaks ----
## 2.1.identify from NG data ----
marker.peak.NG <- getMarkerFeatures(
    ArchRProj = proj_NG,
    useMatrix = "PeakMatrix",
    groupBy = "Cell_type_1",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = c("Adenoma", "Malignant"),
    bgdGroups = "Normal"
)

for (one in c("Adenoma", "Malignant")) {
    message(paste("plot markers Volcano for ", one, " ...", sep = ""))
    pv <- plotMarkers(
        seMarker = marker.peak.NG,
        name = c(one),
        cutOff = "FDR <= 0.05 & abs(Log2FC) >= 0.25",
        plotAs = "Volcano"
    )
    pdf(paste("NG_ATAC/Volcano.markers.", one, ".pdf", sep = ""), 5, 4)
    plot(pv)
    dev.off()
    rm(pv, one)
}

# Diff gene ----

## 1. TCGA RNA data ----
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

## 3. 10X RNA-seq data ----
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

## 2. ATAC data ----
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

## 3. compare DEGs ----
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
