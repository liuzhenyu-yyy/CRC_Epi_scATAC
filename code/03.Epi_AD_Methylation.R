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
peakset <- proj_Epi@peakSet
peakset$peak_id <- paste(seqnames(peakset),
    start(peakset), end(peakset),
    sep = "_"
)
names(peakset) <- peakset$peak_id
table(peakset$peakType)
peakset <- peakset[names(peakset) %in%
    c(diff_peaks$AD_vs_NA$Up$peak_id, diff_peaks$AD_vs_NA$Down$peak_id)]

mycolor <- list(
    "Epi_type" = c(
        "Normal" = "#208a42", "Adenoma" = "#d51f26",
        "Malignant" = "#272e6a"
    ),
    "Epi_Group" = c(
        "Normal" = "#208a42", "Adenoma" = "#d51f26",
        "Group_1" = "#62b7e6", "Group_2" = "#283891"
    ),
    "peakType" = c(
        "Distal" = "#59b795", "Exonic" = "#f8784f",
        "Intronic" = "#7a8dbf", "Promoter" = "#df73b7"
    )
)

# cell type level
sePeaks.epi.type <- getGroupSE(
    ArchRProj = proj_Epi,
    useMatrix = "PeakMatrix",
    groupBy = "Epi_type",
    divideN = TRUE,
    scaleTo = NULL
)
rownames(sePeaks.epi.type) <- paste(rowData(sePeaks.epi.type)$seqnames,
    rowData(sePeaks.epi.type)$start,
    rowData(sePeaks.epi.type)$end,
    sep = "_"
)
rownames(sePeaks.epi.type@assays@data$PeakMatrix) <- rownames(sePeaks.epi.type)
peakmat.type <- sePeaks.epi.type@assays@data$PeakMatrix
peakmat.type <- peakmat.type[c(diff_peaks$AD_vs_NA$Up$peak_id, diff_peaks$AD_vs_NA$Down$peak_id), ]
peakmat.type <- peakmat.type[!duplicated(rownames(peakmat.type)), ]

anno.row <- data.frame(
    row.names = rownames(peakmat.type),
    peak = rep("Down", length(rownames(peakmat.type))),
    peakType = peakset[rownames(peakmat.type)]$peakType
)

anno.row[diff_peaks$AD_vs_NA$Up$peak_id, "peak"] <- "Up"
table(anno.row$peak)

anno.col <- data.frame(
    row.names = c("Normal", "Adenoma", "Malignant"),
    Epi_type = c("Normal", "Adenoma", "Malignant")
)
mycolor$peak <- c("Up" = "#cd2525", "Down" = "#1774cd")

pdf("Heatmap.AD_vs_NA.all.type1.pdf", 5.5, 6)
pheatmap(peakmat.type,
    scale = "row",
    show_rownames = FALSE,
    clustering_method = "ward.D2",
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    annotation_col = anno.col,
    annotation_row = anno.row,
    annotation_colors = mycolor
)
dev.off()
rm(anno.col, anno.row, sePeaks.epi.type)

# cluster level
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
peakmat.cluster <- peakmat.cluster[!duplicated(rownames(peakmat.cluster)), ]

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

## 2.3. Line chart of aggregated peaks ----
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

# 3. correlation of diff peak and methylation ----
load("Methylation_data/methylation.rda")
identical(colnames(beta.mat), rownames(sample.info.array))
sample.info.array$Location <- gsub("CRC", "Malignant", sample.info.array$Location)
sample.info.array$Location <- factor(sample.info.array$Location,
    levels = c("Normal", "Adenoma", "Malignant")
)

## 3.1. align probs to peaks ----
diff_peak.AD <- c(diff_peaks$AD_vs_NA$Up, diff_peaks$AD_vs_NA$Down) %>%
    .[!duplicated(.$peak_id)]
names(diff_peak.AD) <- diff_peak.AD$peak_id

sum(countOverlaps(diff_peaks$AD_vs_NA$Down, probe.bed.hg38) > 0) # 1258
sum(countOverlaps(diff_peaks$AD_vs_NA$Up, probe.bed.hg38) > 0) # 3479

diff_peak.AD.probe <- diff_peak.AD[countOverlaps(diff_peak.AD, probe.bed.hg38) > 0] # 4736
quantile(countOverlaps(diff_peak.AD.probe, probe.bed.hg38))
rm(diff_peak.AD)

diff_peak.AD.probe$nProb <- countOverlaps(diff_peak.AD.probe, probe.bed.hg38)
diff_peak.AD.probe$Prboes <- "none"
hist(diff_peak.AD.probe$nProb)

betamat.peak <- data.frame(row.names = rownames(sample.info.array))
length(names(diff_peak.AD.probe))

# which(names(diff_peak.AD) == one)
# one <- "chr5_150640354_150640854"
for (one in names(diff_peak.AD.probe)) {
    peak <- diff_peak.AD.probe[one]
    probe.selected <- names(probe.bed.hg38)[countOverlaps(probe.bed.hg38, peak) > 0]
    # diff_peak.AD.probe[one, ]$nProb <- length(probe.selected)
    diff_peak.AD.probe[one, ]$Prboes <- paste(probe.selected, collapse = ",")
    if (length(probe.selected) == 1) {
        betamat.peak[, one] <- beta.mat[probe.selected, ] %>% as.numeric()
    } else if (length(probe.selected) > 1) {
        betamat.peak[, one] <- apply(beta.mat[probe.selected, ], 2, mean, na.rm = TRUE) %>% as.numeric()
    }
}
rm(one, peak, probe.selected, probe.bed.hg38)
rm(beta.mat)

betamat.peak <- t(betamat.peak) %>% as.data.frame()
dim(betamat.peak) # 4736 peaks with probe

## 3.2. peak level heatmap ----
anno.row <- data.frame(
    row.names = rownames(betamat.peak),
    peak = rep("Down", length(rownames(betamat.peak))),
    peakType = peakset[rownames(betamat.peak)]$peakType
)
anno.row[rownames(anno.row) %in% diff_peaks$AD_vs_NA$Up$peak_id, "peak"] <- "Up"
table(anno.row$peak)

anno.col <- sample.info.array %>%
    select(Location)

p <- pheatmap(betamat.peak,
    scale = "none",
    show_rownames = FALSE,
    clustering_method = "ward.D2",
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    annotation_col = anno.col,
    annotation_row = anno.row,
    annotation_colors = list(
        "Location" = mycolor$Epi_type,
        "peak" = mycolor$peak
    )
)

plot.data <- betamat.peak[p$tree_row$order, p$tree_col$order]
plot.data <- plot.data[order(rownames(plot.data) %in% diff_peaks$AD_vs_NA$Up$peak_id), ]
plot.data <- plot.data[, order(sample.info.array[colnames(plot.data), ]$Location)]
pdf("Heatmap.beta.AD_vs_NA.all.type1.pdf", 5.5, 6)
pheatmap(plot.data,
    scale = "none",
    color = colorRampPalette(c("#3a2f99", "#d2cc02"))(100),
    show_rownames = FALSE,
    show_colnames = FALSE,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    annotation_col = anno.col,
    annotation_row = anno.row,
    annotation_colors = list(
        "Location" = mycolor$Epi_type,
        "peak" = mycolor$peak,
        "peakType" = mycolor$peakType
    )
)
dev.off()

rm(anno.col, anno.row, plot.data, p)

## 3.3. beesworm of aggregated peaks ----
sample.info.array$Mean_AD_Up <- colMeans(
    betamat.peak[rownames(betamat.peak) %in% diff_peaks$AD_vs_NA$Up$peak_id, ],
    na.rm = TRUE
)
sample.info.array$Mean_AD_Down <- colMeans(
    betamat.peak[rownames(betamat.peak) %in% diff_peaks$AD_vs_NA$Down$peak_id, ],
    na.rm = TRUE
)

plot.data <- data.frame(
    "Location" = levels(sample.info.array$Location),
    "Mean" = aggregate(sample.info.array$Mean_AD_Up, list(sample.info.array$Location), mean)$x,
    "SD" = aggregate(sample.info.array$Mean_AD_Up, list(sample.info.array$Location), sd)$x
)
plot.data$Location <- factor(plot.data$Location,
    levels = c("Normal", "Adenoma", "Malignant")
)
pdf("Line.beta.Cluster.AD_vs_NA.Up.pdf", 5, 3)
ggplot(plot.data, aes(x = Location, y = Mean)) +
    geom_line(aes(group = 1), size = 1) +
    geom_ribbon(aes(x = 1:3, ymin = Mean - SD, ymax = Mean + SD),
        alpha = 0.2
    ) +
    geom_point(aes(color = Location), size = 2, show.legend = FALSE) +
    theme_classic() +
    scale_color_manual(values = mycolor$Epi_type) +
    scale_x_discrete(expand = c(0.02, 0.05)) +
    ylab("Mean beta value")
dev.off()

plot.data <- data.frame(
    "Location" = levels(sample.info.array$Location),
    "Mean" = aggregate(sample.info.array$Mean_AD_Down, list(sample.info.array$Location), mean)$x,
    "SD" = aggregate(sample.info.array$Mean_AD_Down, list(sample.info.array$Location), sd)$x
)
plot.data$Location <- factor(plot.data$Location,
    levels = c("Normal", "Adenoma", "Malignant")
)
pdf("Line.beta.Cluster.AD_vs_NA.Down.pdf", 5, 3)
ggplot(plot.data, aes(x = Location, y = Mean)) +
    geom_line(aes(group = 1), size = 1) +
    geom_ribbon(aes(x = 1:3, ymin = Mean - SD, ymax = Mean + SD),
        alpha = 0.2
    ) +
    geom_point(aes(color = Location), size = 2, show.legend = FALSE) +
    theme_classic() +
    scale_color_manual(values = mycolor$Epi_type) +
    scale_x_discrete(expand = c(0.02, 0.05)) +
    ylab("Mean beta value")
dev.off()
rm(plot.data)

pdf("Beasworm.AD_vs_NA.Up.pdf", 4, 3.5)
ggplot(sample.info.array, aes(x = Location, y = Mean_AD_Up)) +
    geom_beeswarm(aes(color = Location), cex = 2.5, show.legend = FALSE) +
    geom_violin(aes(fill = Location), alpha = 0.2, show.legend = FALSE) +
    scale_color_manual(values = mycolor$Epi_type) +
    ggpubr::stat_compare_means(
        comparisons = list(
            c("Normal", "Adenoma"),
            c("Normal", "Malignant"),
            c("Adenoma", "Malignant")
        )
    ) +
    ylab("Mean beta value") +
    scale_fill_manual(values = mycolor$Epi_type) +
    theme_classic()
dev.off()
pdf("Beasworm.AD_vs_NA.Down.pdf", 4, 3.5)
ggplot(sample.info.array, aes(x = Location, y = Mean_AD_Down)) +
    geom_beeswarm(aes(color = Location), cex = 2.5, show.legend = FALSE) +
    geom_violin(aes(fill = Location), alpha = 0.2, show.legend = FALSE) +
    scale_color_manual(values = mycolor$Epi_type) +
    scale_fill_manual(values = mycolor$Epi_type) +
    ggpubr::stat_compare_means(
        comparisons = list(
            c("Normal", "Adenoma"),
            c("Normal", "Malignant"),
            c("Adenoma", "Malignant")
        )
    ) +
    ylab("Mean beta value") +
    theme_classic()
dev.off()

## 3.4. single gene visulization ----
temp <- peakset[names(diff_peak.AD.probe)] %>%
    subset(peak_id %in% diff_peaks$AD_vs_NA$Down$peak_id) %>%
    subset(peakType %in% c("Promoter", "Exonic", "Intronic")) %>%
    .$nearestGene %>%
    table() %>%
    sort(decreasing = TRUE)
temp %>% head(50)
temp["MAPK11"]

gene.selected <- c(
    "ATP11A", "AXIN2", "KIF26B", "KRT80", "MIR135B", "NKD2", # up
    "NR5A2", "HOXA2", "HOXA5", "EDIL3", "EPB41L3", "CHST2" # down
)
gene.selected <- data.frame(
    row.names = gene.selected,
    gene = gene.selected
)

peakset.selected <- peakset[names(diff_peak.AD.probe)] %>%
    subset(peakType %in% c("Promoter", "Exonic", "Intronic")) %>%
    subset(nearestGene %in% gene.selected$gene)

betamat.gene.selected <- data.frame(row.names = rownames(sample.info.array))
gene.selected$nPeak <- table(peakset.selected$nearestGene)[gene.selected$gene] %>% as.numeric()
gene.selected$nProbe <- 0
hist(gene.selected$nPeak)

# one <- "ATP11A"
for (one in gene.selected$gene) {
    peak <- peakset.selected[peakset.selected$nearestGene == one]
    gene.selected[one, ]$nProbe <- sum(diff_peak.AD.probe[names(peak)]$nProb)

    betamat.gene.selected[, one] <- apply(
        betamat.peak[names(peak), ],
        2, mean,
        na.rm = TRUE
    ) %>% as.numeric()
}
rm(peak, peakset.selected, one, temp)
betamat.gene.selected <- t(betamat.gene.selected) %>% as.data.frame()

p.list <- list()
identical(colnames(betamat.gene.selected), rownames(sample.info.array))

for (one in gene.selected$gene) {
    plot.data <- sample.info.array %>%
        select(c("ID", "Location"))
    plot.data$beta <- betamat.gene.selected[one, ] %>% as.numeric()

    p <- ggplot(plot.data, aes(x = Location, y = beta)) +
        geom_beeswarm(aes(color = Location), cex = 2.5, show.legend = FALSE) +
        geom_violin(aes(fill = Location), alpha = 0.2, show.legend = FALSE) +
        scale_color_manual(values = mycolor$Epi_type) +
        scale_fill_manual(values = mycolor$Epi_type) +
        ggpubr::stat_compare_means(
            comparisons = list(
                c("Normal", "Adenoma"),
                c("Normal", "Malignant"),
                c("Adenoma", "Malignant")
            ),
            size = 3, label.y.npc = 0.5,
        ) +
        ylab("Mean beta value") +
            ggtitle(paste0(one, " (", gene.selected[one, ]$nPrboe, " probes)")) +
            theme_classic() +
            theme(axis.title.x = element_blank())
    p.list[[one]] <- p
}

pdf("Beasworm.AD_vs_NA.selected.gene.pdf", 15, 5)
patchwork::wrap_plots(p.list, ncol = 6)
dev.off()
rm(p.list, plot.data, p, one)

save.image("03.Epi_AD_Methylation.RData")
