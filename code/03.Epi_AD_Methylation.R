setwd("E:/LabWork/Project/CRC_NGS_ATAC/CRC_Epi_scATAC/Results/03.Epi_AD_Methylation")
# load("03.Epi_AD_Methylation.RData")
source("../../code/00.Requirements.R")
library(ggbeeswarm)

# 1. load data ----
## 1.1 load saved project ----
proj_Epi <- loadArchRProject(project.dir.epi, force = TRUE)
table(proj_Epi$Epi_type)

## 1.2 load diff peaks ----
diff_peaks <- list(
    "AD_vs_NA" = list(
        "Up" = read.table(
            file.path("T_AD_C_Diff_Peak", "T_vs_AD_ADhigh_peak.txt"),
            header = TRUE, row.names = 1
        ) %>%
            mutate(peak_id = paste(seqnames, start, end, sep = "_")) %>%
            GRanges(),
        "Down" = read.table(
            file.path("T_AD_C_Diff_Peak", "T_vs_AD_Thigh_peak.txt"),
            header = TRUE, row.names = 1
        ) %>%
            mutate(peak_id = paste(seqnames, start, end, sep = "_")) %>%
            GRanges()
    ),
    "Ca_vs_NA" = list(
        "Up" = read.table(
            file.path("T_AD_C_Diff_Peak", "T_vs_C_Chigh_peak.txt"),
            header = TRUE, row.names = 1
        ) %>%
            mutate(peak_id = paste(seqnames, start, end, sep = "_")) %>%
            GRanges(),
        "Down" = read.table(
            file.path("T_AD_C_Diff_Peak", "T_vs_C_Thigh_peak.txt"),
            header = TRUE, row.names = 1
        ) %>%
            mutate(peak_id = paste(seqnames, start, end, sep = "_")) %>%
            GRanges()
    ),
    "Ca_vs_AD" = list(
        "Up" = read.table(
            file.path("T_AD_C_Diff_Peak", "C_vs_AD_Chigh_peak.txt"),
            header = TRUE, row.names = 1
        ) %>%
            mutate(peak_id = paste(seqnames, start, end, sep = "_")) %>%
            GRanges(),
        "Down" = read.table(
            file.path("T_AD_C_Diff_Peak", "C_vs_AD_ADhigh_peak.txt"),
            header = TRUE, row.names = 1
        ) %>%
            mutate(peak_id = paste(seqnames, start, end, sep = "_")) %>%
            GRanges()
    )
)

table(duplicated(diff_peaks$AD_vs_NA$Up$idx))
table(duplicated(diff_peaks$AD_vs_NA$Up$peak_id))

# 2. AD diff peak keep in cancer ----
## 2.1. overlap with cancer peaks ----
v <- Venn(list(
    "AD_vs_NA" = diff_peaks$AD_vs_NA$Up$peak_id,
    "Ca_vs_NA" = diff_peaks$Ca_vs_NA$Up$peak_id
))
pdf("Venn.Up.peak.vs.NA.pdf", 6, 6)
plot(v, doWeights = TRUE, show = list(Faces = FALSE))
dev.off()

v <- Venn(list(
    "AD_vs_NA" = diff_peaks$AD_vs_NA$Down$peak_id,
    "Ca_vs_NA" = diff_peaks$Ca_vs_NA$Down$peak_id
))
pdf("Venn.Down.peak.vs.NA.pdf", 6, 6)
plot(v, doWeights = TRUE, show = list(Faces = FALSE))
dev.off()
rm(v)

## 2.2. heatmap of AD peaks ----
mycolor <- list(
    "Epi_type" = c(
        "Normal" = "#208a42", "Adenoma" = "#d51f26",
        "Malignant" = "#272e6a"
    ),
    "Epi_Group" = c(
        "Normal" = "#208a42", "Adenoma" = "#d51f26",
        "Group_1" = "#62b7e6", "Group_2" = "#283891"
    )
)

sePeaks <- getGroupSE(
    ArchRProj = proj_Epi,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters",
    divideN = TRUE,
    scaleTo = NULL
)

rownames(sePeaks) <- paste(rowData(sePeaks)$seqnames,
    rowData(sePeaks)$start,
    rowData(sePeaks)$end,
    sep = "_"
)
rownames(sePeaks@assays@data$PeakMatrix) <- rownames(sePeaks)
peakmat.cluster <- sePeaks@assays@data$PeakMatrix
peakmat.cluster <- peakmat.cluster[c(diff_peaks$AD_vs_NA$Up$peak_id, diff_peaks$AD_vs_NA$Down$peak_id), ]

cluster.info <- colData(sePeaks) %>% as.data.frame()
temp <- readRDS("../02.Epi_Molecular_Subtype/cluster.info.rds")

cluster.info$Epi_type <- "Malignant"
cluster.info[c("C3", "C4"), ]$Epi_type <- "Normal"
cluster.info[c("C28", "C9"), ]$Epi_type <- "Adenoma"
cluster.info$Epi_type <- factor(cluster.info$Epi_type,
    levels = c("Normal", "Adenoma", "Malignant")
)

cluster.info$Epi_Group <- cluster.info$Epi_type
cluster.info[rownames(temp), ]$Epi_Group <- temp$Epi_Group
cluster.info$Epi_Group <- factor(cluster.info$Epi_Group,
    levels = c("Normal", "Adenoma", "Group_1", "Group_2")
)
rm(temp, sePeaks)

anno.col <- cluster.info %>% select(c("Epi_type", "Epi_Group"))

# up peak
plot.data <- peakmat.cluster[diff_peaks$AD_vs_NA$Up$peak_id, ]
plot.data <- plot.data[, order(cluster.info$Epi_Group)]

pdf("Heatmap.AD_vs_NA.Up.pdf", 6, 6)
pheatmap(plot.data,
    scale = "row",
    show_rownames = FALSE,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    annotation_col  = anno.col,
    annotation_colors = mycolor
)
dev.off()

# down peak
plot.data <- peakmat.cluster[diff_peaks$AD_vs_NA$Down$peak_id, ]
plot.data <- plot.data[, order(cluster.info$Epi_Group)]

pdf("Heatmap.AD_vs_NA.Down.pdf", 6, 6)
pheatmap(plot.data,
    scale = "row",
    show_rownames = FALSE,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    annotation_col = anno.col,
    annotation_colors = mycolor
)
dev.off()
rm(plot.data, anno.col)

## 2.3. Line chart of AD peaks ----
# cluter resolution
identical(rownames(cluster.info), colnames(peakmat.cluster))

cluster.info$Mean_AD_Up <- colMeans(
    peakmat.cluster[diff_peaks$AD_vs_NA$Up$peak_id, ]
)
cluster.info$Mean_AD_Down <- colMeans(
    peakmat.cluster[diff_peaks$AD_vs_NA$Down$peak_id, ]
)

plot.data <- data.frame(
    "Epi_type" = levels(cluster.info$Epi_type),
    "Mean" = aggregate(cluster.info$Mean_AD_Up, list(cluster.info$Epi_type), mean)$x,
    "SD" = aggregate(cluster.info$Mean_AD_Up, list(cluster.info$Epi_type), sd)$x
)
plot.data$Epi_type <- factor(plot.data$Epi_type,
    levels = c("Normal", "Adenoma", "Malignant")
)
pdf("Line.Cluster.AD_vs_NA.Up.pdf", 5, 3)
ggplot(plot.data, aes(x = Epi_type, y = Mean)) +
    geom_line(aes(group = 1), size = 1) +
    geom_ribbon(aes(x = 1:3, ymin = Mean - SD, ymax = Mean + SD),
        alpha = 0.2
    ) +
    geom_point(aes(color = Epi_type), size = 2, show.legend = FALSE) +
    theme_classic() +
    scale_color_manual(values = mycolor$Epi_type) +
    scale_x_discrete(expand = c(0.02, 0.05)) +
    ylab("Mean ATAC signal")
dev.off()

plot.data <- data.frame(
    "Epi_type" = levels(cluster.info$Epi_type),
    "Mean" = aggregate(cluster.info$Mean_AD_Down, list(cluster.info$Epi_type), mean)$x,
    "SD" = aggregate(cluster.info$Mean_AD_Down, list(cluster.info$Epi_type), sd)$x
)
plot.data$Epi_type <- factor(plot.data$Epi_type,
    levels = c("Normal", "Adenoma", "Malignant")
)
pdf("Line.Cluster.AD_vs_NA.Down.pdf", 5, 3)
ggplot(plot.data, aes(x = Epi_type, y = Mean)) +
    geom_line(aes(group = 1), size = 1) +
    geom_ribbon(aes(x = 1:3, ymin = Mean - SD, ymax = Mean + SD),
        alpha = 0.2
    ) +
    geom_point(aes(color = Epi_type), size = 2, show.legend = FALSE) +
    theme_classic() +
    scale_color_manual(values = mycolor$Epi_type) +
    scale_x_discrete(expand = c(0.02, 0.05)) +
    ylab("Mean ATAC signal")
dev.off()
rm(plot.data)

# # single cell resolution
# PeakMat <- getMatrixFromProject(
#     ArchRProj = proj_Epi,
#     useMatrix = "PeakMatrix"
# )
# sample.info <- colData(PeakMat) %>% as.data.frame()

# peakset <- proj_Epi@peakSet
# peakset$peak_id <- paste(seqnames(peakset),
#     start(peakset), end(peakset),
#     sep = "_"
# )
# rownames(PeakMat) <- peakset$peak_id
# rownames(PeakMat@assays@data$PeakMatrix) <- peakset$peak_id
# PeakMat <- PeakMat@assays@data$PeakMatrix
# PeakMat[1:5, 1:5]
# PeakMat <- PeakMat[c(diff_peaks$AD_vs_NA$Up$peak_id, diff_peaks$AD_vs_NA$Down$peak_id), ]

# PeakMat <- imputeMatrix(PeakMat,
#     imputeWeights = getImputeWeights(proj_Epi)
# )

# identical(colnames(PeakMat), rownames(sample.info))

# sample.info$Mean_AD_Up <- colMeans(
#     PeakMat[diff_peaks$AD_vs_NA$Up$peak_id, ]
# )
# sample.info$Mean_AD_Down <- colMeans(
#     PeakMat[diff_peaks$AD_vs_NA$Down$peak_id, ]
# )
# sample.info$Epi_type <- factor(sample.info$Epi_type,
#     levels = c("Normal", "Adenoma", "Malignant")
# )
# # rm(PeakMat, sample.info, peakset, plot.data)

# pdf("Line.cell.AD_vs_NA.Up.pdf", 5, 3)
# plot.data <- data.frame(
#     "Epi_type" = levels(sample.info$Epi_type),
#     "Mean" = aggregate(sample.info$Mean_AD_Up, list(sample.info$Epi_type), mean)$x,
#     "SD" = aggregate(sample.info$Mean_AD_Up, list(sample.info$Epi_type), sd)$x
# )
# plot.data$Epi_type <- factor(plot.data$Epi_type,
#     levels = c("Normal", "Adenoma", "Malignant")
# )
# ggplot(plot.data, aes(x = Epi_type, y = Mean)) +
#     geom_point(aes(color = Epi_type)) +
#     geom_line(aes(group = 1)) +
#     geom_ribbon(aes(x = 1:3, ymin = Mean - SD, ymax = Mean + SD), alpha = 0.2) +
#     theme_classic()
# dev.off()

# pdf("Line.cell.AD_vs_NA.Down.pdf", 5, 3)
# plot.data <- data.frame(
#     "Epi_type" = levels(sample.info$Epi_type),
#     "Mean" = aggregate(sample.info$Mean_AD_Down, list(sample.info$Epi_type), mean)$x,
#     "SD" = aggregate(sample.info$Mean_AD_Down, list(sample.info$Epi_type), sd)$x
# )
# plot.data$Epi_type <- factor(plot.data$Epi_type,
#     levels = c("Normal", "Adenoma", "Malignant")
# )
# ggplot(plot.data, aes(x = Epi_type, y = Mean)) +
#     geom_point(aes(color = Epi_type)) +
#     geom_line(aes(group = 1)) +
#     geom_ribbon(aes(x = 1:3, ymin = Mean - SD, ymax = Mean + SD), alpha = 0.2) +
#     theme_classic()
# dev.off()


# 3. correlation of diff peak and methylation ----

save.image("03.Epi_AD_Methylation.RData")
