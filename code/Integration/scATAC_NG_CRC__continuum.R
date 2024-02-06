#######################################################
#
# Re-analysis of Nature Genetics CRC 10X scATAC continuum
# Data source: NCBI Gene Expression Omnibus GSE201349
# Citation: https://www.nature.com/articles/s41588-022-01088-x
#
#####################################################

setwd("E:/LabWork/Project/CRC_NGS_ATAC/CRC_Epi_scATAC/Results/07.Independent_Validation")
load("Independent_Validation.ATAC.RData")

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
proj_NG <- saveArchRProject(proj_NG, load = TRUE)
proj_NG <- loadArchRProject("./Project_Dir_NG")

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
    maxClusters = 50,
    resolution = 1.5,
    force = TRUE
)
table(proj_NG$Clusters)

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
    "EPCAM", "LOC101928540", "MIR4463", "PCSK9", "F2RL3",
    "PECAM1", "CDH5", "VIM", "COL2A1", "LGR5", "ASCL2"
)
proj_NG <- addImputeWeights(proj_NG)

p <- plotEmbedding(
    ArchRProj = proj_NG, colorBy = "GeneScoreMatrix",
    pal = ArchRPalettes$blueYellow,
    name = markers, embedding = "UMAP",
    imputeWeights = getImputeWeights(proj_NG),
    size = 0.2, plotAs = "points"
)
pdf("NG_ATAC/UMAP.Marker.Gene.pdf", 25, 15)
patchwork::wrap_plots(p, ncol = 6)
dev.off()

# gene module
DEG.ATAC <- read.table("DEG.ATAC.txt", sep = "\t", header = TRUE) %>%
    as.list() %>%
    sapply(function(x) x %>% setdiff(""))
sapply(DEG.ATAC, length)

colnames(proj_NG@cellColData)
proj_NG <- addModuleScore(proj_NG,
    useMatrix = "GeneScoreMatrix",
    name = "Module",
    features = DEG.ATAC[5:6]
)

p <- plotEmbedding(
    ArchRProj = proj_NG, colorBy = "cellColData",
    pal = ArchRPalettes$blueYellow,
    name = c("Module.AD_high", "Module.AD_low"), embedding = "UMAP",
    imputeWeights = getImputeWeights(proj_NG),
    size = 0.2, plotAs = "points"
)
pdf("NG_ATAC/UMAP.AD.module.pdf", 10, 5)
patchwork::wrap_plots(p, ncol = 2)
dev.off()

# composition
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

sample.info.NG <- proj_NG@cellColData %>% as.data.frame()
sample.info.NG$Cell_type_1 <- "none"
sample.info.NG[sample.info.NG$Clusters %in% paste0("C", c(34)), "Cell_type_1"] <- "T"
sample.info.NG[sample.info.NG$Clusters %in% paste0("C", c(35:37)), "Cell_type_1"] <- "B"
sample.info.NG[sample.info.NG$Clusters %in% paste0("C", c(32, 33)), "Cell_type_1"] <- "Myeloid"
sample.info.NG[sample.info.NG$Clusters %in% paste0("C", c(29:31)), "Cell_type_1"] <- "Stromal"
sample.info.NG[sample.info.NG$Clusters %in% paste0("C", c(1:5, 22:24)), "Cell_type_1"] <- "Malignant"
sample.info.NG[sample.info.NG$Clusters %in% paste0("C", c(12, 18:19)), "Cell_type_1"] <- "Normal"
sample.info.NG[sample.info.NG$Clusters %in% paste0("C", c(11, 13, 20, 21, 17)), "Cell_type_1"] <- "Adenoma_Early"
sample.info.NG[sample.info.NG$Clusters %in% paste0("C", c(14:16, 7:10)), "Cell_type_1"] <- "Adenoma_Late"
sample.info.NG[sample.info.NG$Clusters %in% paste0("C", c(6, 25:28)), "Cell_type_1"] <- "Intermediate"
sample.info.NG[sample.info.NG$Clusters %in% paste0("C", c(19)), "Cell_type_1"] <- "Stem-like"
table(sample.info.NG$Cell_type_1)
table(sample.info.NG$Cell_type_1, sample.info.NG$Clusters)
proj_NG$Cell_type_1 <- sample.info.NG$Cell_type_1

p <- plotEmbedding(
    ArchRProj = proj_NG, colorBy = "cellColData",
    name = "Cell_type_1", embedding = "UMAP",
    size = 0.2, plotAs = "points",
    labelMeans = FALSE
)
pdf("NG_ATAC/UMAP.Cell_type_12.pdf", 7, 7)
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
    useGroups = c("Adenoma_Early", "Adenoma_Late", "Malignant"),
    bgdGroups = "Stem-like"
)

for (one in c("Adenoma_Early", "Adenoma_Late", "Malignant")) {
    message(paste("plot markers Volcano for ", one, " ...", sep = ""))
    pv <- plotMarkers(
        seMarker = marker.peak.NG,
        name = c(one),
        cutOff = "FDR <= 0.05 & abs(Log2FC) >= 0.1",
        plotAs = "Volcano"
    )
    pdf(paste("NG_ATAC/Volcano.markers.stem.", one, ".pdf", sep = ""), 5, 4)
    plot(pv)
    dev.off()
    rm(pv, one)
}

marker.NG.up.list <- getMarkers(marker.peak.NG,
    cutOff = "FDR <= 0.05 & Log2FC >= 0.1"
) %>% lapply(function(x) {
    x <- x %>%
        as.data.frame() %>%
        mutate(
            peak_id = paste(seqnames, start, end, sep = "_")
        )
    return(x$peak_id)
})
marker.NG.down.list <- getMarkers(marker.peak.NG,
    cutOff = "FDR <= 0.05 & Log2FC <= (-0.1)"
) %>% lapply(function(x) {
    x <- x %>%
        as.data.frame() %>%
        mutate(
            peak_id = paste(seqnames, start, end, sep = "_")
        )
    return(x$peak_id)
})

sapply(marker.NG.up.list, length)
sapply(marker.NG.down.list, length)

## 2.2 compare with out data ----
marker.up.list <- list(
    "Adenoma" = read.table(
        file.path("data/T_AD_C_Diff_Peak", "T_vs_AD_ADhigh_peak.txt"),
        header = TRUE, row.names = 1
    ) %>%
        mutate(peak_id = paste(seqnames, start, end, sep = "_")) %>%
        pull(peak_id),
    "Malignant" = read.table(
        file.path("data/T_AD_C_Diff_Peak", "T_vs_C_Chigh_peak.txt"),
        header = TRUE, row.names = 1
    ) %>%
        mutate(peak_id = paste(seqnames, start, end, sep = "_")) %>%
        pull(peak_id)
)
marker.down.list <- list(
    "Adenoma" = read.table(
        file.path("data/T_AD_C_Diff_Peak", "T_vs_AD_Thigh_peak.txt"),
        header = TRUE, row.names = 1
    ) %>%
        mutate(peak_id = paste(seqnames, start, end, sep = "_")) %>%
        pull(peak_id),
    "Malignant" = read.table(
        file.path("data/T_AD_C_Diff_Peak", "T_vs_C_Thigh_peak.txt"),
        header = TRUE, row.names = 1
    ) %>%
        mutate(peak_id = paste(seqnames, start, end, sep = "_")) %>%
        pull(peak_id)
)
sapply(marker.up.list, length)
sapply(marker.down.list, length)

v <- Vennerable::Venn(list(
    "NG_continuum" = marker.NG.up.list$Malignant,
    "scATAC" = marker.up.list$Malignant
))
pdf("NG_ATAC/Venn.Malignant.stem.Up.pdf", 10, 10)
plot(v, doWeights = TRUE, doEuler = TRUE, show = list(Faces = FALSE))
dev.off()

v <- Vennerable::Venn(list(
    "NG_continuum" = marker.NG.down.list$Malignant,
    "scATAC" = marker.down.list$Malignant
))
pdf("NG_ATAC/Venn.Malignant.Down.stem.pdf", 10, 10)
plot(v, doWeights = TRUE, doEuler = TRUE, show = list(Faces = FALSE))
dev.off()

v <- Vennerable::Venn(list(
    "NG_continuum" = marker.NG.up.list$Adenoma_Late,
    "scATAC" = marker.up.list$Adenoma
))
pdf("NG_ATAC/Venn.Adenoma.stem.Up.pdf", 10, 10)
plot(v, doWeights = TRUE, doEuler = TRUE, show = list(Faces = FALSE))
dev.off()

v <- Vennerable::Venn(list(
    "NG_continuum" = marker.NG.down.list$Adenoma_Late,
    "scATAC" = marker.down.list$Adenoma
))
pdf("NG_ATAC/Venn.Adenoma.stem.Down.pdf", 10, 10)
plot(v, doWeights = TRUE, doEuler = TRUE, show = list(Faces = FALSE))
dev.off()

rm(v, p, plot.data)
rm(marker.peak.NG)

# 3. TF modules ----
library(WGCNA)

## 3.1. Motif matrix and module eigengene ----
proj_NG <- saveArchRProject(proj_NG, load = TRUE)
proj_NG <- loadArchRProject("./Project_Dir_NG")
proj_NG@peakAnnotation
proj_Epi@peakAnnotation$Motif$Positions
proj_Epi@peakAnnotation$Motif$Matches

proj_NG <- addMotifAnnotations(
    ArchRProj = proj_NG,
    motifSet = "cisbp", name = "Motif"
)
proj_NG <- addBgdPeaks(proj_NG)

proj_NG <- addDeviationsMatrix(
    ArchRProj = proj_NG,
    peakAnnotation = "Motif",
    force = TRUE
)

# group by clusters
getAvailableMatrices(proj_NG)
seMotif.cluster <- getGroupSE(
    ArchRProj = proj_NG,
    useMatrix = "MotifMatrix",
    groupBy = "Clusters",
    divideN = TRUE,
    scaleTo = NULL
)

dim(seMotif.cluster)
rowData(seMotif.cluster)
seMotif.cluster <- seMotif.cluster[rowData(seMotif.cluster)$seqnames == "deviations", ]
rownames(seMotif.cluster@assays@data$MotifMatrix) <- rowData(seMotif.cluster)$name

MotifMat.cluster <- seMotif.cluster@assays@data$MotifMatrix
rownames(MotifMat.cluster)

MotifMat.cluster["ENSG00000250542_156", ]
MotifMat.cluster <- MotifMat.cluster[rownames(MotifMat.cluster) != "ENSG00000250542_156", ]
MotifMat.cluster <- t(MotifMat.cluster)
MotifMat.cluster <- as.data.frame(MotifMat.cluster)
rm(seMotif.cluster)

net <- readRDS("../05.Epi_TF_Clustering/WGCNA.net.rds")
ME.NG <- moduleEigengenes(MotifMat.cluster,
    colors = net$colors
)
ME.NG <- ME.NG$eigengenes

## 3.2 identify iCMS subtypes ----
markers.iCMS <- read.table("E:/LabWork/Project/CRC_NGS_ATAC/iCMS markers.txt",
    stringsAsFactors = FALSE,
    sep = "\t", header = TRUE
)
markers.iCMS <- base::as.list(markers.iCMS) %>% sapply(function(x) setdiff(x, ""))
markers.iCMS <- lapply(markers.iCMS, function(x) intersect(x, getFeatures(proj_Epi, useMatrix = "GeneScoreMatrix")))
lapply(markers.iCMS, length)

proj_NG <- addModuleScore(proj_NG,
    useMatrix = "GeneScoreMatrix",
    name = "Module",
    features = markers.iCMS
)
proj_NG$Module.iCMS2 <- proj_NG$Module.iCMS2_Up - proj_NG$Module.iCMS2_Down
proj_NG$Module.iCMS3 <- proj_NG$Module.iCMS3_Up - proj_NG$Module.iCMS3_Down

cluster.malignant <- paste0("C", c(1:5, 22:24))

sample.info.NG <- proj_NG@cellColData %>% as.data.frame()
colnames(sample.info.NG)

temp <- sample.info.NG %>%
    filter(Clusters %in% cluster.malignant)
temp <- table(temp$Clusters, temp$Patient)
colnames(temp)[apply(temp, 1, which.max)]

pdf("NG_ATAC/Violin.Module.iCMS.pdf", 6, 3)
sample.info.NG %>%
    subset(Clusters %in% cluster.malignant) %>%
    select(Clusters, Module.iCMS2, Module.iCMS3) %>%
    reshape2::melt(id.vars = "Clusters") %>%
    mutate(
        Clusters = factor(Clusters, levels = c("C2", "C3", "C5","C1", "C4", "C22", "C23", "C24")),
        iCMS = ifelse(Clusters %in% c("C1", "C2", "C3", "C4", "C5"), "iCMS2", "iCMS3")
    ) %>%
    ggplot(aes(x = Clusters, y = value, fill = iCMS)) +
        geom_violin() +
        scale_fill_manual(values = c("iCMS2" = "#283891", "iCMS3" = "#62b7e6")) +
        facet_wrap(~variable, scales = "free_y") +
        ylab("iCMS module signature") +
    ylim(-1, 1) +
    theme_classic()
dev.off()

cluster.info.NG <- data.frame(
    row.names = cluster.malignant,
    Clusters = cluster.malignant,
    iCMS = ifelse(cluster.malignant %in% c("C1", "C2", "C3", "C4", "C5"), "iCMS2", "iCMS3")
)

cluster.info.NG <- cbind(cluster.info.NG, ME.NG[cluster.malignant, ])
cluster.info.NG$ME5 <- 0 - cluster.info.NG$ME5

pdf("NG_ATAC/Box.TF.module.all.iCMS.pdf", 10, 8)
cluster.info.NG %>%
    mutate(
        iCMS = factor(iCMS, levels = c("iCMS2", "iCMS3")),
        Clusters = factor(Clusters, levels = cluster.malignant)
    ) %>%
    reshape2::melt(id.vars = c("Clusters", "iCMS")) %>%
    ggplot(aes(x = iCMS, y = value, fill = iCMS)) +
    geom_boxplot(outlier.color = NA) +
    geom_jitter(width = 0.2) +
    scale_fill_manual(values = c("iCMS2" = "#283891", "iCMS3" = "#62b7e6")) +
    ggpubr::stat_compare_means() +
    facet_wrap(~variable, ncol = 5) +
    ylab("Module eigengene") +
    ylim(-0.5, 0.5) +
    theme_classic()
dev.off()

pdf("NG_ATAC/Box.TF.module.iCMS.pdf", 8, 3)
cluster.info.NG %>%
    mutate(
        iCMS = factor(iCMS, levels = c("iCMS2", "iCMS3")),
        Clusters = factor(Clusters, levels = cluster.malignant)
    ) %>%
    select(c("Clusters", "iCMS", "ME5", "ME8", "ME10", "ME11")) %>%
    reshape2::melt(id.vars = c("Clusters", "iCMS")) %>%
    ggplot(aes(x = iCMS, y = value, fill = iCMS)) +
    geom_boxplot(outlier.color = NA) +
    geom_jitter(width = 0.2) +
    scale_fill_manual(values = c("iCMS2" = "#283891", "iCMS3" = "#62b7e6")) +
    ggpubr::stat_compare_means() +
    facet_wrap(~variable, ncol = 4) +
    ylab("Module eigengene") +
    ylim(-0.5, 0.5) +
    theme_classic()
dev.off()

rm(p, plot.data, net)
gc()
save.image("Independent_Validation.ATAC.RData")
