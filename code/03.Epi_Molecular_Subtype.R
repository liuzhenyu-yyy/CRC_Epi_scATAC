setwd("E:/LabWork/Project/CRC_NGS_ATAC/CRC_Epi_scATAC/Results/03.Epi_Molecular_Subtype")
source("../../code/00.Requirements.R")
load("Epi_Molecular_Subtype.RData")

# 1. prepare ArchR project ----
## 1.1 load saved project ----
proj_Epi <- loadArchRProject(project.dir.epi, force = TRUE)

proj_Epi$Patient <- gsub("-nofacs", "", proj_Epi$Sample)
table(proj_Epi$Patient)

proj_Epi$Epi_type <- "Malignant"
proj_Epi$Epi_type[proj_Epi$Clusters %in% c("C3", "C4")] <- "Normal"
proj_Epi$Epi_type[proj_Epi$Clusters %in% c("C28", "C9")] <- "Adenoma"

table(proj_Epi$Epi_type, proj_Epi$new_location)

p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    name = "Epi_type", embedding = "UMAP",
    size = 0.2, plotAs = "points",
    labelMeans = FALSE
)

pdf("UAMP.Epi.Epi_type.pdf", 7, 7)
plot(p)
dev.off()

p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    name = "Clusters", embedding = "UMAP",
    size = 0.2, plotAs = "points"
)

pdf("UAMP.Epi.Clusters.pdf", 7, 7)
plot(p)
dev.off()

sample.info.epi <- proj_Epi@cellColData %>% as.data.frame()

AMI.cluster <- c(
    Normal = aricode::AMI(
        sample.info.epi %>% filter(Epi_type == "Normal") %>% pull(Clusters),
        sample.info.epi %>% filter(Epi_type == "Normal") %>%  pull(Sample)
    ),
    Adenoma = aricode::AMI(
        sample.info.epi %>% filter(Epi_type == "Adenoma") %>% pull(Clusters),
        sample.info.epi %>%
            filter(Epi_type == "Adenoma") %>%
            pull(Sample)
    ),
    Malignant = aricode::AMI(
        sample.info.epi %>% filter(Epi_type == "Malignant") %>% pull(Clusters),
        sample.info.epi %>% filter(Epi_type == "Malignant") %>% pull(Sample)
    )
)

## 1.2. add clinical information ----
patient.info <- read.table("E:/LabWork/Project/CRC_NGS_ATAC/patient.info.txt",
    header = TRUE, sep = "\t"
)
patient.info[is.na(patient.info)] <- "NA"
colnames(patient.info)
rownames(patient.info) <- patient.info$Patient

table(sample.info.epi$Patient %in% rownames(patient.info))

temp <- patient.info[sample.info.epi$Patient, ]
rownames(temp) <- rownames(sample.info.epi)
temp$Patient <- NULL

for (one in colnames(temp)) {
    proj_Epi <- addCellColData(proj_Epi,
        data = temp[, one],
        name = one,
        cells = rownames(temp),
        force = TRUE
    )
}
proj_Epi$MSI_Status[proj_Epi$Epi_type == "Normal"] <- "MSS"
sample.info.epi <- proj_Epi@cellColData %>% as.data.frame()
rm(temp, one)

# color palette
RColorBrewer::display.brewer.all()
mycolor <- list(
    "Gender" = c("Female" = "#cab2d6", "Male" = "#6a3d9a"),
    "MSI" = c("MSS" = "#ffff99", "MSI-H" = "#b15928", "NA" = "gray70"),
    "Side" = c("Right" = "#b2df8a", "Left" = "#33a02c"),
    "Group" = c("Group_1" = "#62b7e6", "Group_2" = "#283891")
)

p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    name = "Gender", embedding = "UMAP",
    pal = mycolor$Gender,
    size = 0.2, plotAs = "points",
    labelMeans = FALSE
)
pdf("UAMP.Epi.Gender.pdf", 7, 7)
plot(p)
dev.off()

p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    name = "MSI_Status", embedding = "UMAP",
    pal = mycolor$MSI,
    size = 0.2, plotAs = "points",
    labelMeans = FALSE
)
pdf("UAMP.Epi.MSI_Status1.pdf", 7, 7)
plot(p)
dev.off()

p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    name = "Side", embedding = "UMAP",
    pal = mycolor$Side,
    size = 0.2, plotAs = "points",
    labelMeans = FALSE
)
pdf("UAMP.Epi.Side.pdf", 7, 7)
plot(p)
dev.off()

## 1.3. pseudo-bulk repicates ----
proj_Epi@projectMetadata$GroupCoverages$Clusters$Param
proj_Epi@projectMetadata$GroupCoverages$Clusters$coverageMetadata[, 3] # all cell??

proj_Epi <- addGroupCoverages(
    ArchRProj = proj_Epi,
    groupBy = "Clusters",
    useLabels = FALSE,
    force = TRUE
)
table(proj_Epi$Clusters)

getGroupBW(
    ArchRProj = proj_Epi,
    groupBy = "Clusters",
    normMethod = "ReadsInTSS",
    tileSize = 30,
    maxCells = 2000,
    threads = 4
)
getAvailableMatrices(proj_Epi)

sePeaks <- getGroupSE(
    ArchRProj = proj_Epi,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters",
    divideN = TRUE,
    scaleTo = NULL
)
saveRDS(sePeaks, "sePeaks.cluster.rds")

write.csv(as.data.frame(colData(sePeaks)), "cluster.info.all.csv")
sample.selected <- setdiff(rownames(sePeaks@colData), c("C3", "C4", "C28", "C9"))

## 1.4. annotate clusters ----
cluster.info <- colData(sePeaks) %>%
    as.data.frame() %>%
    .[sample.selected, ]

# annotate clusters
for (one in c("Patient", "Gender", "MSI_Status", "Side")) {
    temp <- aggregate(sample.info.epi[, one],
        by = list(sample.info.epi$Clusters),
        FUN = function(x) {
            sort(table(x), decreasing = TRUE)[1] %>% names()
        }
    )
    rownames(temp) <- temp$Group.1
    cluster.info[, paste0(one, "_Major")] <- temp[sample.selected, ]$x

    temp <- aggregate(sample.info.epi[, one],
        by = list(sample.info.epi$Clusters),
        FUN = function(x) {
            sort(table(x), decreasing = TRUE)[1]
        }
    )
    rownames(temp) <- temp$Group.1
    cluster.info[, paste0(one, "_Major_Count")] <- temp[sample.selected, ]$x
    cluster.info[, paste0(one, "_Major_Pct")] <- cluster.info[, paste0(one, "_Major_Count")] / cluster.info$nCells
}
rm(temp, one)

# 2. clustering in peak space ----
library(NMF)

## 2.1. NMF of malignant clusters ----
# feature selection by SD
SDs <- apply(sePeaks@assays@data$PeakMatrix, 1, sd)
quantile(SDs)
feature.selcted <- sort(SDs, decreasing = TRUE) %>%
    head(10000) %>%
    names()

NMF.test <- nmf(
    sePeaks@assays@data$PeakMatrix[feature.selcted, sample.selected],
    rank = 2:10,
    nrun = 30
)
pdf("NMF.rank.test.pdf", 8, 6)
plot(NMF.test)
dev.off()

# run NMF with 2 clusters
NMF.res <- nmf(
    sePeaks@assays@data$PeakMatrix[feature.selcted, sample.selected],
    rank = 2,
    nrun = 200
)
rm(sePeaks)

# ref <- c(1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2)
# names(ref) <- sample.selected
# table(ref, predict(NMF.res))
# rm(ref)

group.res <- predict(NMF.res, what = "consensus") %>%
    as.numeric()
group.res <- group.res %>%
    paste("Group_", ., sep = "")
names(group.res) <- sample.selected

cluster.info$Epi_Group <- group.res[rownames(cluster.info)]
sample.info.epi$Epi_Group <- group.res[sample.info.epi$Clusters]

proj_Epi$Epi_Group <- proj_Epi$Epi_type
proj_Epi$Epi_Group[proj_Epi$Epi_Group == "Malignant"] <-
    group.res[proj_Epi$Clusters[proj_Epi$Epi_Group == "Malignant"]]
table(proj_Epi$Epi_Group)

mycolor[["Epi_Group"]] <- c(
    "Normal" = "#208a42", "Adenoma" = "#d51f26",
    "Group_1" = "#62b7e6", "Group_2" = "#283891"
)
p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    name = "Epi_Group", embedding = "UMAP",
    pal = mycolor$Epi_Group,
    size = 0.2, plotAs = "points",
    labelMeans = FALSE
)
pdf("UAMP.Epi.Epi_Group.pdf", 7, 7)
plot(p)
dev.off()

saveRDS(cluster.info, "cluster.info.rds")

## 2.2. visualize NMF results ----
# conseusus matrix
con.mat <- NMF.res@consensus
sil <- silhouette(NMF.res, what = "consensus")

anno.col <- cluster.info %>%
    select(c("Gender_Major", "MSI_Status_Major", "Side_Major", "Epi_Group"))
colnames(anno.col) <- names(mycolor)[c(1:4)]
anno.col$silhouette <- sil[, 3][rownames(anno.col)]

pdf("NMF.consensus.clusters.pdf", 7, 6)
pheatmap(con.mat,
    annotation_col = anno.col[c(5, 4, 3, 2, 1)],
    annotation_colors = mycolor[colnames(anno.col)],
    border_color = NA,
    cutree_rows = 2,
    cutree_cols = 2
)
dev.off()

# coefficient matrix
coefmap(NMF.res)
coef.mat <- coef(NMF.res)
rownames(coef.mat) <- c("Basis1", "Basis2")
write.csv(coef.mat, "NMF.coefficient.clusters.csv")
pdf("NMF.coefficient.clusters.pdf", 6, 3)
pheatmap(coef.mat[, order(coef.mat[1, ] - coef.mat[2, ])],
    annotation_col = anno.col[c(4)],
    annotation_colors = mycolor[colnames(anno.col)],
    cluster.rows = FALSE,
    cluster_cols = FALSE,
    border_color = NA,
    clustering_method = "complete"
)
dev.off()

rm(con.mat, coef.mat, sil, p, anno.col)

# 3. compare with known iCMS features ----
## 3.1. with clinical information ----
group.count <- table(cluster.info$Epi_Group)

# MSI
plot.data <- table(cluster.info %>%
    select(c("Epi_Group", "MSI_Status_Major"))) %>%
    as.data.frame()
plot.data$Epi_Group <- factor(plot.data$Epi_Group,
    levels = c("Group_2", "Group_1")
)
p1 <- ggplot(plot.data) +
    geom_bar(aes(x = Epi_Group, y = Freq, fill = MSI_Status_Major),
        stat = "identity", position = position_fill(reverse = TRUE)
    ) +
    geom_text(aes(x = Epi_Group, y = Freq, label = Freq),
        position = position_fill(vjust = 0.5)
    ) +
    scale_fill_manual(values = mycolor$MSI) +
        coord_flip() +
        theme_classic() +
        theme(axis.title = element_blank(), axis.text.x = element_blank())

cluster.info %>%
    # filter(MSI_Status_Major != "NA") %>%
    select(c("Epi_Group", "MSI_Status_Major")) %>%
    table() %>%
    chisq.test()

## Side
plot.data <- table(cluster.info %>%
    select(c("Epi_Group", "Side_Major"))) %>%
    as.data.frame()
plot.data$Epi_Group <- factor(plot.data$Epi_Group,
    levels = c("Group_2", "Group_1")
)
p2 <- ggplot(plot.data) +
    geom_bar(aes(x = Epi_Group, y = Freq, fill = Side_Major),
        stat = "identity", position = position_fill(reverse = TRUE)
    ) +
    geom_text(aes(x = Epi_Group, y = Freq, label = Freq),
        position = position_fill(vjust = 0.5)
    ) +
    scale_fill_manual(values = mycolor$Side) +
        coord_flip() +
        theme_classic() +
        theme(axis.title = element_blank(), axis.text.x = element_blank())

cluster.info %>%
    select(c("Epi_Group", "Side_Major")) %>%
    table() %>%
    chisq.test()

## Gender
plot.data <- table(cluster.info %>%
    select(c("Epi_Group", "Gender_Major"))) %>%
    as.data.frame()
plot.data$Epi_Group <- factor(plot.data$Epi_Group,
    levels = c("Group_2", "Group_1")
)
p3 <- ggplot(plot.data) +
    geom_bar(aes(x = Epi_Group, y = Freq, fill = Gender_Major),
        stat = "identity", position = position_fill(reverse = TRUE)
    ) +
    geom_text(aes(x = Epi_Group, y = Freq, label = Freq),
        position = position_fill(vjust = 0.5)
    ) +
    scale_fill_manual(values = mycolor$Gender) +
        coord_flip() +
        theme_classic() +
        theme(axis.title.y = element_blank()) +
        ylab("Fraction")
cluster.info %>%
    select(c("Epi_Group", "Gender_Major")) %>%
    table() %>%
    chisq.test()

pdf("Bar.Clinical.Group.pdf", 6, 4)
wrap_plots(p1, p2, p3, ncol = 1)
dev.off()
rm(plot.data, p1, p2, p3, group.count)

## 3.2. CNV profile ----
load("../01.All_scCNV/CRC_CNV.rda")
rm(sample.info, anno.row)

cell.selected <- sample.info.epi %>%
    filter(Clusters %in% sample.selected) %>%
    rownames()

anno.row <- sample.info.epi[cell.selected, ] %>%
    select(c("Side", "MSI_Status", "Epi_Group"))
colnames(anno.row)[2] <- "MSI"
CNV.FC <- CNV.FC[cell.selected, ]

png("Heatmap.CNV.cell.wardd2.png", 1000, 750)
p <- pheatmap::pheatmap(CNV.FC,
    color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
    border_color = FALSE,
    cluster_rows = TRUE, cluster_cols = FALSE,
    clustering_method = "ward.D2",
    annotation_col = anno.col,
    annotation_row = anno.row,
    annotation_colors = c(anno.color, mycolor)[c(colnames(anno.row), "seqnames")],
    annotation_legend = FALSE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    gaps_col = which(!duplicated(anno.col$seqnames)) - 1
)
dev.off()
CNV.cluster <- cutree(p$tree_row, k = 2)
identical(names(CNV.cluster), rownames(anno.row))
table(CNV.cluster, anno.row$Epi_Group) %>% chisq.test()

rm(CNV.FC, anno.col, anno.row, anno.color)

## 3.3. iCMS signatures ----
markers.iCMS <- read.table("E:/LabWork/Project/CRC_NGS_ATAC/iCMS markers.txt",
    stringsAsFactors = FALSE,
    sep = "\t", header = TRUE
)
markers.iCMS <- base::as.list(markers.iCMS)

genes <- getFeatures(proj_Epi, useMatrix = "GeneScoreMatrix")

markers.iCMS <- lapply(markers.iCMS, function(x) intersect(x, genes))
lapply(markers.iCMS, length)

proj_Epi <- addModuleScore(proj_Epi,
    useMatrix = "GeneScoreMatrix",
    name = "Module",
    features = markers.iCMS
)

proj_Epi$Module.iCMS2 <- proj_Epi$Module.iCMS2_Up - proj_Epi$Module.iCMS2_Down
proj_Epi$Module.iCMS3 <- proj_Epi$Module.iCMS3_Up - proj_Epi$Module.iCMS3_Down

sample.info.epi <- proj_Epi@cellColData %>% as.data.frame()
colnames(sample.info.epi)

sePeaks <- readRDS("sePeaks.cluster.rds")
plot.data <- sePeaks@colData %>%
    as.data.frame() %>%
    mutate("Epi_Group" = "none")
plot.data[c("C3", "C4"), ]$Epi_Group <- "Normal"
plot.data[rownames(cluster.info), ]$Epi_Group <- cluster.info$Epi_Group
plot.data <- plot.data[plot.data$Epi_Group != "none", ]

pdf("Dot.Module.iCMS.clusters.pdf", 5, 3.7)
ggplot(plot.data, aes(x = Module.iCMS2, y = Module.iCMS3)) +
    geom_point(aes(color = Epi_Group), size = 2) +
    stat_ellipse(aes(group = Epi_Group, color = Epi_Group), level = 0.9) +
    scale_color_manual(values = mycolor$Epi_Group) +
    ggpubr::stat_cor() +
    theme_classic()
dev.off()

plot.data <- sample.info.epi %>%
    filter(Epi_type != "Adenoma") %>%
    select(matches("Module|Group"))

pdf("Dot.Module.iCMS.cells2.pdf", 5, 3.7)
ggplot(plot.data, aes(x = Module.iCMS2, y = Module.iCMS3)) +
    geom_point(aes(color = Epi_Group), size = 0.6) +
    stat_ellipse(aes(group = Epi_Group), level = 0.9) +
    scale_color_manual(values = mycolor$Epi_Group) +
    ggpubr::stat_cor(label.x.npc = 0.5) +
    xlim(c(-30, 40)) +
    ylim(c(-50, 60)) +
    theme_classic()
dev.off()

plot.data <- melt(plot.data[, c("Epi_Group", "Module.iCMS3", "Module.iCMS2")], id.vars = c("Epi_Group"))
plot.data$Epi_Group <- factor(plot.data$Epi_Group, levels = c("Normal", "Group_2", "Group_1"))
plot.data$variable <- factor(plot.data$variable, levels = c("Module.iCMS2", "Module.iCMS3"))

pdf("Violin.Module.iCMS.pdf", 4.5, 3.5)
ggplot(plot.data, aes(x = Epi_Group, y = value)) +
    geom_violin(aes(fill = Epi_Group),
        show.legend = FALSE, size = 0.5
    ) +
    scale_fill_manual(values = mycolor$Epi_Group) +
    ggpubr::stat_compare_means(comparisons = list(
        c("Group_1", "Group_2"), c("Normal", "Group_2"),
        c("Group_1", "Normal")
    )) +
    facet_wrap(~variable, scales = "free_y") +
    ylab("Module Score") +
    theme_classic()
dev.off()

p1 <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    pal = ArchRPalettes$blueYellow,
    name = "Module.iCMS2", embedding = "UMAP",
    size = 0.2, plotAs = "points"
)
p2 <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    pal = ArchRPalettes$blueYellow,
    name = "Module.iCMS3", embedding = "UMAP",
    size = 0.2, plotAs = "points"
)
pdf("UAMP.Epi.Module.iCMS.pdf", 8, 4)
wrap_plots(p1, p2, ncol = 2)
dev.off()

proj_Epi$Module.iCMS <- proj_Epi$Module.iCMS3 - proj_Epi$Module.iCMS2

p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    pal = ArchRPalettes$blueYellow,
    name = "Module.iCMS", embedding = "UMAP",
    size = 0.2, plotAs = "points"
)
pdf("UAMP.Epi.Module.iCMS.diff.pdf", 5, 5)
plot(p)
dev.off()

rm(genes, plot.data, p, p1, p2)

gene.selected <- c(
    "EREG", "TIMP3", "MYC", "EIF6", "CTSA",
    # "KRT23", "AREG", "CPNE1", "MYC", "PTPRO","FN1", "CXCL14",
    "TNFRSF10B", "PFKP", "CDKN2A", "ATL3", "CLU"
    # , "DNAJC10", "TSTA3", "MYOF", "DUSP4", "TM4SF4"
)
proj_Epi <- addImputeWeights(proj_Epi, reducedDims = "IterativeLSI_merge")
p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "GeneScoreMatrix",
    pal = ArchRPalettes$blueYellow,
    name = gene.selected, embedding = "UMAP",
    imputeWeights = getImputeWeights(proj_Epi),
    size = 0.2, plotAs = "points"
)

pdf("UMAP.markers.iCMS.pdf", 20, 8)
patchwork::wrap_plots(plotlist = p, nrow = 2, ncol = 5, byrow = TRUE)
dev.off()

lapply(markers.iCMS, length)

# 4. differential peak of iCMS ----
## 4.1. identify diff peaks ----
# call diff peaks
table(proj_Epi$Epi_Group)

marker.peak.vs.Normal <- getMarkerFeatures(
    ArchRProj = proj_Epi,
    useMatrix = "PeakMatrix",
    groupBy = "Epi_Group",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = c("Group_1", "Group_2"),
    bgdGroups = "Normal"
)
# volcano plot
for (one in c("Group_1", "Group_2")) {
    message(paste("plot markers Volcano for ", one, " ...", sep = ""))
    pv <- plotMarkers(
        seMarker = marker.peak.vs.Normal,
        name = c(one),
        cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1",
        plotAs = "Volcano"
    )
    pdf(paste("Volcano.markers.", one, ".pdf", sep = ""), 5, 4)
    plot(pv)
    dev.off()
    rm(pv, one)
}

# heatmap
sePeaks <- getGroupSE(
    ArchRProj = proj_Epi,
    useMatrix = "PeakMatrix",
    groupBy = "Epi_Group",
    divideN = TRUE,
    scaleTo = NULL
)

rownames(marker.peak.vs.Normal) <- paste0("f", rownames(marker.peak.vs.Normal))
colnames(marker.peak.vs.Normal)

temp <- marker.peak.vs.Normal@assays@data
peak.G1.up <- rownames(marker.peak.vs.Normal)[temp$Log2FC[, 1] >= 1 & temp$FDR[, 1] <= 0.01]
peak.G1.down <- rownames(marker.peak.vs.Normal)[temp$Log2FC[, 1] <= (-1) & temp$FDR[, 1] <= 0.01]
peak.G2.up <- rownames(marker.peak.vs.Normal)[temp$Log2FC[, 2] >= 1 & temp$FDR[, 2] <= 0.01]
peak.G2.down <- rownames(marker.peak.vs.Normal)[temp$Log2FC[, 2] <= (-1) & temp$FDR[, 2] <= 0.01]
rm(temp)

pdf("Upset.markers.group_vs_normal.pdf", 6, 4.5)
upset(
    fromList(
        list(
            "G1.up" = peak.G1.up,
            "G2.up" = peak.G2.up,
            "G1.down" = peak.G1.down,
            "G2.down" = peak.G2.down
        )
    ),
    nset = 4,
    order.by = "freq"
)
dev.off()

anno.row <- data.frame(
    row.names = unique(c(peak.G1.up, peak.G1.down, peak.G2.up, peak.G2.down))
)
anno.row$Group <-  "none"
anno.row$Direction <- "none"

anno.row[unique(c(peak.G1.up, peak.G2.up)), ]$Direction <- "Up"
anno.row[unique(c(peak.G1.down, peak.G2.down)), ]$Direction <- "Down"

anno.row[unique(c(peak.G1.up, peak.G1.down)), ]$Group <- "Group_1"
anno.row[unique(c(peak.G2.up, peak.G2.down)), ]$Group <- "Group_2"

anno.row[intersect(peak.G1.up, peak.G2.up), ]$Group <- "Both"
anno.row[intersect(peak.G1.down, peak.G2.down), ]$Group <- "Both"

anno.row[intersect(peak.G2.down, peak.G1.up), ]$Group <- "none"
anno.row[intersect(peak.G1.down, peak.G2.up), ]$Group <- "none"
anno.row <- anno.row[anno.row$Group != "none", ]
table(anno.row$Group, anno.row$Direction)

anno.row$Group <- factor(anno.row$Group, levels = c("Group_1", "Group_2", "Both"))
anno.row$Direction <- factor(anno.row$Direction, levels = c("Up", "Down"))

anno.row <- anno.row[order(
    anno.row$Direction,
    anno.row$Group,
    sePeaks@assays@data$PeakMatrix[rownames(anno.row), "Normal"]
), ]

plot.data <- sePeaks@assays@data$PeakMatrix[
    rownames(anno.row)[anno.row$Direction == "Up"],
    c("Normal", "Group_1", "Group_2")
]
png("Heatmap.marker.tmuor.cluster.Up.png", 600, 800)
pheatmap(plot.data,
    scale = "row",
    annotation_row = anno.row,
    annotation_colors = list(
        Group = c(Group_1 = "#62b7e6", Group_2 = "#283891", Both = "#86d786"),
        Direction = c(Up = "#f57474", Down = "#90cdf0")),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = FALSE
)
dev.off()

plot.data <- sePeaks@assays@data$PeakMatrix[
    rownames(anno.row)[anno.row$Direction == "Down"],
    c("Normal", "Group_1", "Group_2")
]

png("Heatmap.marker.tmuor.cluster.Down.png", 600, 800)
pheatmap(plot.data,
    scale = "row",
    annotation_row = anno.row,
    annotation_colors = list(
        Group = c(Group_1 = "#62b7e6", Group_2 = "#283891", Both = "#86d786"),
        Direction = c(Up = "#f57474", Down = "#90cdf0")),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = FALSE
)
dev.off()

rm(plot.data, peak.G1.up, peak.G1.down, peak.G2.up, peak.G2.down)

## 4.2. GO enrichment of diff peaks ----
# marker peak to nearest gene
marker.peak.list <- list(
    "Group_1" = rownames(anno.row)[anno.row$Group == "Group_1" & anno.row$Direction == "Up"],
    "Group_2" = rownames(anno.row)[anno.row$Group == "Group_2" & anno.row$Direction == "Up"],
    "Both" = rownames(anno.row)[anno.row$Group == "Both" & anno.row$Direction == "Up"]
)
lapply(marker.peak.list, length)
rm(anno.row)

peakset <- proj_Epi@peakSet
names(peakset) <- paste(
    seqnames(peakset), start(peakset), end(peakset),
    sep = "_"
)

temp <- rowData(sePeaks) %>% as.data.frame()
temp$Name <- paste(temp$seqnames, temp$start, temp$end, sep = "_")
table(temp$Name %in% names(peakset))

marker.peak.list <- lapply(
    marker.peak.list,
    function(x) temp[x, ]$Name
)

table(peakset$peakType)
marker.peak2gene.list <- lapply(
    marker.peak.list,
    function(x) {
        peakset[x] %>%
           #.[.$peakType != "Distal"] %>%
            .$nearestGene %>%
            unique()
    }
)
common <- intersect(marker.peak2gene.list$Group_1, marker.peak2gene.list$Group_2)
marker.peak2gene.list$Both <- c(marker.peak2gene.list$Both, common) %>% unique()
marker.peak2gene.list$Group_1 <- setdiff(marker.peak2gene.list$Group_1, common)
marker.peak2gene.list$Group_2 <- setdiff(marker.peak2gene.list$Group_2, common)

lapply(marker.peak2gene.list, length)
rm(temp, sePeaks, peakset, common)

# GO analysis
library(clusterProfiler)
GO.marker.list <- lapply(
    marker.peak2gene.list,
    function(x) {
        res <- clusterProfiler::enrichGO(x,
            OrgDb = org.Hs.eg.db,
            keyType = "SYMBOL",
            ont = "BP",
            pAdjustMethod = "BH"
        )
        res <- clusterProfiler::simplify(res)
        return(res@result)
    }
)
object.size(marker.peak.list) / 1e6

write.csv(GO.marker.list$Group_1, "GO.marker.real.Group_1.csv")
write.csv(GO.marker.list$Group_2, "GO.marker.real.Group_2.csv")
write.csv(GO.marker.list$Both, "GO.marker.real.Both.csv")
plot.data <- rbind(
    GO.marker.list$Group_1[c("GO:0033674", "GO:0043410", "GO:0032956", "GO:0034329", "GO:0050727"), ],
    GO.marker.list$Group_2[c("GO:0001763", "GO:0198738", "GO:0016055", "GO:0019827", "GO:0051591"), ],
    GO.marker.list$Both[c("GO:0007409", "GO:0048568", "GO:0090132", "GO:0045785", "GO:2001236"), ]
)
plot.data$Group <- factor(
    rep(c("Group_1", "Group_2", "Both"), each = 5),
    levels = c("Group_1", "Group_2", "Both")
)

plot.data <- plot.data[order(plot.data$Group,
    plot.data$p.adjust,
    decreasing = TRUE
), ]
plot.data$Description <- factor(plot.data$Description, levels = plot.data$Description)

pdf("GO.marker.group.pdf", 7, 5)
ggplot(plot.data) +
    geom_bar(aes(x = (0 - log10(p.adjust)), y = Description, fill = Group),
        stat = "identity", position = "dodge", width = 0.4
    ) +
    scale_fill_manual(values = c("#62b7e6", "#283891", "#86d786")) +
    facet_wrap(~Group, ncol = 1, scales = "free_y") +
        theme_classic() +
        xlab("minus log10 adjusted p-value") +
        ylab("GO term")
dev.off()
rm(plot.data)

## 4.3. track analysis of important genes ----
# export bed files
dir.create("bed")
# marker.peak.tumor.gr <- getMarkers(marker.peak.vs.Normal,
#     cutOff = "FDR <= 0.01 & Log2FC >= 1",
#     returnGR = TRUE
# )
# lapply(marker.peak.tumor.gr, length)
# write.table(as.data.frame(marker.peak.tumor.gr$Group_1)[, 1:3],
#     "marker.peak.tumor.All.Group_1.bed",
#     col.names = FALSE,
#     sep = "\t", quote = FALSE, row.names = FALSE
# )
# write.table(as.data.frame(marker.peak.tumor.gr$Group_2)[, 1:3],
#     "marker.peak.tumor.All.Group_2.bed",
#     col.names = FALSE,
#     sep = "\t", quote = FALSE, row.names = FALSE
# )

write.table(
    marker.peak.list$Group_1 %>%
        gsub("_", "\t", .) %>%
        as.data.frame() %>%
        mutate(
            "Name" = marker.peak.list$Group_1,
            "length" = 500,
            "strand" = "."
        ),
    "bed/marker.peak.tumor.Group_1.specific.bed",
    col.names = FALSE,
    sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
    marker.peak.list$Group_2 %>%
        gsub("_", "\t", .) %>%
        as.data.frame() %>%
        mutate(
            "Name" = marker.peak.list$Group_2,
            "length" = 500,
            "strand" = "."
        ),
    "bed/marker.peak.tumor.Group_2.specific.bed",
    col.names = FALSE,
    sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
    marker.peak.list$Both %>%
        gsub("_", "\t", .) %>%
        gsub("_", "\t", .) %>%
        as.data.frame() %>%
        mutate(
            "Name" = marker.peak.list$Both,
            "length" = 500,
            "strand" = "."
        ),
    "bed/marker.peak.tumor.common.bed",
    col.names = FALSE,
    sep = "\t", quote = FALSE, row.names = FALSE
)
temp <- unlist(marker.peak.list) %>%
    gsub("_", "\t", .) %>%
    gsub("_", "\t", .) %>%
    as.data.frame() %>%
    mutate("a" = unlist(marker.peak.list), "b" = 0, "c" = "-", "d" = 1, "e" = 1)
c(Group_1 = "#62b7e6", Group_2 = "#283891", Both = "#86d786")
temp$itemRgb <- rep(c(Group_1 = "98,183,230", Group_2 = "40,56,145", Both = "134,215,134"),
    times = c(length(marker.peak.list$Group_1), length(marker.peak.list$Group_2), length(marker.peak.list$Both))
)
table(temp$itemRgb)
write.table(temp, "bed/marker.peak.tumor.all.bed",
    col.names = FALSE,
    sep = "\t", quote = FALSE, row.names = FALSE
)

# gene.selected <- c(
#     "CXCL14", "EREG", "TIMP3", "KRT23", "FN1", "AREG", "CPNE1", "MYC", "PTPRO", "EIF6", "CTSA",
#     "TNFRSF10B", "DNAJC10", "TSTA3", "PFKP", "CDKN2A", "ATL3", "CLU", "MYOF", "DUSP4", "TM4SF4"
# )
# track.subtype <- plotBrowserTrack(
#     ArchRProj = proj_Epi, groupBy = "Epi_Group",
#     geneSymbol = gene.selected,
#     features = marker.peak.tumor.gr,
#     upstream = 50000, downstream = 50000,
#     minCells = 10,
#     pal = mycolor$Epi_Group,
#     normMethod = "ReadsInTSS"
# )

# pdf("Track.Marker.iCMS.pdf", 30, 16)
# patchwork::wrap_plots(plotlist = track.subtype, nrow = 4, ncol = 6, byrow = TRUE)
# dev.off()
# change to cluster view in IGV

rm(gene, marker.peak.tumor.gr, p, track.subtype)

## 4.4. GREAT enrichment of diff peaks ----
GREAT.marker.list <- list(
    "Group_1" = read.table("GREAT/Group1_greatExportAll.tsv", header = FALSE, sep = "\t"),
    "Group_2" = read.table("GREAT/Group2_greatExportAll.tsv", , header = FALSE, sep = "\t"),
    "Common" = read.table("GREAT/Common_greatExportAll.tsv", , header = FALSE, sep = "\t")
)
sapply(GREAT.marker.list, dim)
colnames(GREAT.marker.list$Group_1) <- c(
    "Ontology", "ID", "Desc",
    "BinomRank", "BinomP", "BinomBonfP", "BinomFdrQ", "RegionFoldEnrich", "ExpRegions",
    "ObsRegions", "GenomeFrac", "SetCov", "HyperRank", "HyperP", "HyperBonfP",
    "HyperFdrQ", "GeneFoldEnrich", "ExpGenes", "ObsGenes", "TotalGenes", "GeneSetCov",
    "TermCov", "Regions", "Genes"
)
colnames(GREAT.marker.list$Group_2) <- colnames(GREAT.marker.list$Group_1)
colnames(GREAT.marker.list$Common) <- colnames(GREAT.marker.list$Group_1)
rownames(GREAT.marker.list$Group_1) <- GREAT.marker.list$Group_1$ID
rownames(GREAT.marker.list$Group_2) <- GREAT.marker.list$Group_2$ID
rownames(GREAT.marker.list$Common) <- GREAT.marker.list$Common$ID

plot(Vennerable::Venn(list(
    "Group_1" = GREAT.marker.list$Group_1$ID,
    "Group_2" = GREAT.marker.list$Group_2$ID,
    "Common" = GREAT.marker.list$Common$ID
)))
View(GREAT.marker.list$Group_1)

plot.data <- rbind(
    GREAT.marker.list$Group_1[c("GO:0043069", "GO:0002274", "GO:1903037", "GO:0045321", "GO:0050863", "GO:0043408"), ],
    GREAT.marker.list$Group_2[c("GO:0034969", "GO:0018216", "GO:1903465", "GO:1902153", "GO:1902808", "GO:0048378"), ],
    GREAT.marker.list$Common[c("GO:0030856", "GO:0002009", "GO:0043065", "GO:0006366", "GO:0035295", "GO:0035019"), ]
)
plot.data$Group <- factor(
    rep(c("Group_1", "Group_2", "Common"), each = 6),
    levels = c("Group_2", "Group_1", "Common")
)

plot.data <- plot.data[order(plot.data$Group,
    plot.data$BinomBonfP,
    decreasing = TRUE
), ]
plot.data$Desc <- factor(plot.data$Desc, levels = plot.data$Desc)
# View(plot.data)
pdf("GREAT.marker.group.pdf", 8, 6)
ggplot(plot.data) +
    geom_bar(aes(x = (0 - log10(BinomBonfP)), y = Desc, fill = Group),
        stat = "identity", position = "dodge", width = 0.4
    ) +
    scale_fill_manual(values = c("#283891", "#62b7e6",  "#86d786")) +
    facet_wrap(~Group, ncol = 1, scales = "free") +
    theme_classic() +
    xlab("minus log10 adjusted p-value") +
    ylab("GO term")
dev.off()
rm(plot.data, GREAT.marker.list, GO.marker.list)

# 5. TF enrichment in iCMS marker peaks ----
dir.create("TF_motif")
## 5.1. run HOMER ----
# findMotifsGenome.pl bed/marker.peak.tumor.Group_1.specific.bed hg38 homer/Group_1 -size 200

homer.res <- list(
    "Group_1" = homer.parser("homer/Group_1/knownResults.txt") %>%
        mutate(Group = "Group_1"),
    "Group_2" = homer.parser("homer/Group_2/knownResults.txt") %>%
        mutate(Group = "Group_2"),
    "Common" = homer.parser("homer/Common/knownResults.txt") %>%
        mutate(Group = "Common")
)
sapply(homer.res, dim) # 419 TFs
identical(homer.res$Group_1$TF, homer.res$Group_2$TF)

## 5.2. identify significant TFs ----
pdf("TF_motif/Dot.motif.Group1.pdf", 5, 4)
ggplot(homer.res$Group_1, aes(x = Log2_Enrichment, y = log.p.value)) +
    geom_point(aes(color = Diff), size = 1) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    geom_hline(yintercept = 50, linetype = "dashed") +
    scale_color_manual(values = c("none" = "grey", "up" = "red")) +
    ggrepel::geom_text_repel(aes(label = TF), size = 2, max.overlaps = 30) +
    theme_classic()
dev.off()
pdf("TF_motif/Dot.motif.Group2.pdf", 5, 4)
ggplot(homer.res$Group_2, aes(x = Log2_Enrichment, y = log.p.value)) +
    geom_point(aes(color = Diff), size = 1) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    geom_hline(yintercept = 50, linetype = "dashed") +
    scale_color_manual(values = c("none" = "grey", "up" = "red")) +
    ggrepel::geom_text_repel(aes(label = TF), size = 2, max.overlaps = 30) +
    theme_classic()
dev.off()
pdf("TF_motif/Dot.motif.common.pdf", 5, 4)
ggplot(homer.res$Common, aes(x = Log2_Enrichment, y = log.p.value)) +
    geom_point(aes(color = Diff), size = 1) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    geom_hline(yintercept = 50, linetype = "dashed") +
    scale_color_manual(values = c("none" = "grey", "up" = "red")) +
    ggrepel::geom_text_repel(aes(label = TF), size = 2, max.overlaps = 30) +
    theme_classic()
dev.off()

TF.sig <- lapply(homer.res, function(x) {
    x <- x[x$Diff == "up", ]
    return(x$TF)
}) %>%
    do.call(c, .) %>%
    unname() %>%
    unique()
TF.sig <- grep(":|-", TF.sig, invert = TRUE, value = TRUE)
TF.sig <- c("AP-1", TF.sig)

FC.mat <- data.frame(
    row.names = homer.res$Group_1$TF,
    "Group_1" = homer.res$Group_1$Log2_Enrichment,
    "Common" = homer.res$Common$Log2_Enrichment,
    "Group_2" = homer.res$Group_2$Log2_Enrichment
)

pdf("TF_motif/Heatmap.motif.sig.pdf", 10, 4)
pheatmap(t(FC.mat[TF.sig, ]),
    color = colorRampPalette(ArchRPalettes$comet)(100),
    scale = "column",
    cutree_cols = 4,
    cluster_rows = FALSE, clustering_method = "ward.D2"
)
dev.off()

# sort by Order
p <- pheatmap(t(FC.mat[TF.sig, ]),
    scale = "column",
    cluster_rows = FALSE, clustering_method = "ward.D2"
)
TF.sig <- TF.sig[p$tree_col$order]

plot.data <- lapply(homer.res, function(x) {
    x <- x[x$TF %in% TF.sig, ]
    return(x)
}) %>% do.call(rbind, .)

plot.data$Group <- factor(plot.data$Group, levels = c("Group_1", "Common", "Group_2"))
plot.data$TF <- factor(plot.data$TF, levels = TF.sig)

plot.data[plot.data$log.p.value > 2000, ]$log.p.value <- 2000
plot.data[plot.data$Log2_Enrichment > 2, ]$Log2_Enrichment <- 2

pdf("TF_motif/Dot.motif.sig.pdf", 5, 7)
ggplot(plot.data) +
    geom_point(aes(
        x = Group, y = TF,
        fill = log.p.value, size = Log2_Enrichment
    ), pch = 21) +
    scale_fill_viridis_c() +
    theme_bw()
dev.off()

## 5.3. combine with Motif deviation & gene score ----
table(proj_Epi$Epi_Group)
# marker of motif deviation
marker.motif.vs.Normal <- getMarkerFeatures(
    ArchRProj = proj_Epi,
    useMatrix = "MotifMatrix",
    groupBy = "Epi_Group",
    useGroups = c("Group_1", "Group_2"),
    bgdGroups = "Normal"
)

marker.motif <- getMarkers(marker.motif.vs.Normal,
    cutOff = "FDR <= 0.01 & MeanDiff > 0.05"
) %>% lapply(., as.data.frame)
marker.motif$Group_1

TF.CISBP <- getFeatures(proj_Epi, useMatrix = "MotifMatrix")

# align CISBP and HOMER
TF.sig.align <- sapply(TF.sig, function(x) {
    res <- paste("z:", x, sep = "")
    res <- grep(res, TF.CISBP, value = TRUE)
    return(res)
})  %>% unlist()
names(TF.sig.align) <- TF.sig.align %>%
    gsub("z:", "", .) %>%
    gsub("_.+?$", "", .)
sapply(TF.sig, function(x) {
    res <- paste("z:", x, sep = "")
    res <- grep(res, TF.CISBP, value = TRUE)
    return(res)
})["NRF2"]

# plot gene score
p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "GeneScoreMatrix",
    name = names(TF.sig.align), embedding = "UMAP",
    imputeWeights = getImputeWeights(proj_Epi),
    size = 0.2, plotAs = "points"
)
pdf("TF_motif/UMAP.TF.sig.GeneScore.pdf", 5, 5)
for (one in p) {
    plot(one)
}
dev.off()

# plot deviation
p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "MotifMatrix",
    name = TF.sig.align, embedding = "UMAP",
    imputeWeights = getImputeWeights(proj_Epi),
    size = 0.2, plotAs = "points"
)
pdf("TF_motif/UMAP.TF.sig.MotifZ.pdf", 5, 5)
for (one in p) {
    plot(one)
}
dev.off()

# correlation of motif deviation and gene score
marker.gene.vs.Normal <- getMarkerFeatures(
    ArchRProj = proj_Epi,
    useMatrix = "GeneScoreMatrix",
    groupBy = "Epi_Group",
    useGroups = c("Group_1", "Group_2"),
    bgdGroups = "Normal"
)
TF.all <- getFeatures(proj_Epi, useMatrix = "MotifMatrix") %>%
    grep("deviations:", ., value = TRUE) %>%
        gsub("deviations:", "", .)

names(TF.all) <- TF.all  %>%
    gsub("_.+?$", "", .)

TF.all <- TF.all[names(TF.all) %in% marker.gene.vs.Normal@elementMetadata$name]

MeanDiff <- marker.motif.vs.Normal@assays@data$MeanDiff
colnames(MeanDiff) <- c("Group_1", "Group_2")
rownames(MeanDiff) <- marker.motif.vs.Normal@elementMetadata$name
MeanDiff <- MeanDiff[TF.all, ]

FDR <- marker.motif.vs.Normal@assays@data$FDR
colnames(FDR) <- c("Group_1", "Group_2")
rownames(FDR) <- marker.motif.vs.Normal@elementMetadata$name
FDR <- FDR[TF.all, ]

Log2FC <- marker.gene.vs.Normal@assays@data$Log2FC
colnames(Log2FC) <- c("Group_1", "Group_2")
rownames(Log2FC) <- marker.gene.vs.Normal@elementMetadata$name
Log2FC <- Log2FC[names(TF.all), ]

plot(MeanDiff$Group_1, Log2FC$Group_1)
plot.data.1 <- data.frame(
    row.names = names(TF.all),
    TF = names(TF.all),
    motif.diff = MeanDiff$Group_1,
    FDR = FDR$Group_1,
    gene.diff = Log2FC$Group_1,
    iCMS = "iCMS3"
)

plot.data.2 <- data.frame(
    row.names = names(TF.all),
    TF = names(TF.all),
    motif.diff = MeanDiff$Group_2,
    FDR = FDR$Group_2,
    gene.diff = Log2FC$Group_2,
    iCMS = "iCMS2"
)
plot.data.1["FOXA3", ]
plot.data <- rbind(plot.data.1, plot.data.2)
TF.selected <- c(
    "NRF2", "MAFB",
    "MAFK", "FOXA3", "FOXA2", "FOXM1", "SOX2", "SOX4",
    "ELF1", "EHF", "ETS1",
    "AP-1", "LEF1", "TCF3",
    "NUR77", "HNF1", "CDX2", "PPARA", "TR4", "HNF4A"
)
TF.selected <- c(TF.selected, consensus.TF$iCMS2, consensus.TF$iCMS3)
plot.data$Label <- plot.data$TF
plot.data$Label[!plot.data$TF %in% TF.selected] <- NA

plot.data$m_log10_FDR <- -log10(plot.data$FDR)
plot.data$m_log10_FDR[plot.data$m_log10_FDR > 50] <- 50

pdf("TF_motif/Scatter.TF.GS.deviation.pdf", 8, 3.7)
ggplot(plot.data, aes(x = gene.diff, y = motif.diff)) +
    geom_point(aes(color = m_log10_FDR)) +
    ggrepel::geom_text_repel(aes(label = Label), size = 2.5, max.overlaps = 50) +
    facet_wrap(~iCMS, ncol = 2, scales = "free") +
    scale_color_viridis_c() +
    theme_bw()
dev.off()

pdf("TF_motif/Scatter.TF.GS.deviation.all.pdf", 8, 3.7)
ggplot(plot.data, aes(x = gene.diff, y = motif.diff)) +
    geom_point(aes(color = m_log10_FDR)) +
    ggrepel::geom_text_repel(aes(label = TF), size = 2.5, max.overlaps = 15) +
    facet_wrap(~iCMS, ncol = 2, scales = "free") +
    scale_color_viridis_c() +
    theme_bw()
dev.off()

rm(p, one, TF.CISBP)

## 5.4. select for visulization ----
# dot plot for TFs
TF.selected <- c(
    "ELF1", "EHF", "ETS1",
    "FOXM1", "MAFK", "FOXA3",
    "AP-1", "P53", "LEF1", "TCF3",
    "CDX2", "PPARA", "TR4", "HNF4A", "CTCF"
)

pdf("TF_motif/Heatmap.motif.sig.selected.pdf", 10, 4)
pheatmap(t(FC.mat[TF.selected, ]),
    cluster_cols = FALSE,
    scale = "column",
    cluster_rows = FALSE
)
dev.off()

TF.selected <- c(
    "NRF2", "MAFB",
    "MAFK", "FOXA3", "FOXA2", "FOXM1", "SOX2", "SOX4",
    "ELF1", "EHF", "ETS1",
    "AP-1", "LEF1", "TCF3",
    "NUR77", "HNF1", "CDX2", "PPARA", "TR4", "HNF4A"
)

TF.selected <- setdiff(TF.selected, "CTCF")
plot.data <- lapply(homer.res, function(x) {
    x <- x[x$TF %in% TF.selected, ]
    return(x)
}) %>% do.call(rbind, .)

plot.data$Group <- factor(plot.data$Group, levels = c("Group_2", "Common", "Group_1"))
plot.data$TF <- factor(plot.data$TF, levels = (TF.selected))

plot.data[plot.data$log.p.value > 500, ]$log.p.value <- 500
#plot.data[plot.data$Log2_Enrichment > 3, ]$Log2_Enrichment <- 3

pdf("TF_motif/Dot.motif.sig.selected.pdf", 4, 4)
ggplot(plot.data) +
    geom_point(aes(
        x = Group, y = TF,
        fill = log.p.value, size = Log2_Enrichment
    ), pch = 21) +
    scale_fill_viridis_c() +
    theme_bw()
dev.off()

rm(FC.mat, plot.data)

# umap
TF.selected <- TF.sig.align[c("MAFK", "FOXA3", "LEF1", "TCF3", "PPARA", "HNF4A")]
p1 <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "GeneScoreMatrix",
    pal = ArchRPalettes$blueYellow,
    name = names(TF.selected), embedding = "UMAP",
    imputeWeights = getImputeWeights(proj_Epi),
    size = 0.2, plotAs = "points"
)

p2 <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "MotifMatrix",
    name = TF.selected, embedding = "UMAP",
    imputeWeights = getImputeWeights(proj_Epi),
    size = 0.2, plotAs = "points"
)
pdf("TF_motif/UMAP.TF.sig.selected.pdf", 30, 14)
wrap_plots(c(p1, p2), ncol = 6)
dev.off()
rm(p1, p2)

# footprint
# temp <- readRDS(file.path(project.dir.epi, "Save-ArchR-Project_raw.rds"))
# temp@peakAnnotation$Motif$Positions <- paste0(project.dir.epi, "/Annotations/Motif-Positions-In-Peaks.rds")
# temp@peakAnnotation$Motif$Matches <- paste0(project.dir.epi, "/Annotations/Motif-Matches-In-Peaks.rds")
# proj_Epi@peakAnnotation <- temp@peakAnnotation
# rm(temp)

motifPositions <- getPositions(proj_Epi)

proj_Epi <- addGroupCoverages(
    ArchRProj = proj_Epi,
    groupBy = "Epi_Group"
)
source("E:/LabWork/code/archr_footprint_no_ribbon.R")

seFoot <- getFootprints(
    ArchRProj = proj_Epi,
    positions = motifPositions[gsub("z:", "", TF.selected)],
    groupBy = "Epi_Group",
    useGroups = c("Normal", "Group_1", "Group_2"),
)

plotFootprints(
    seFoot = seFoot,
    ArchRProj = proj_Epi,
    normMethod = "Subtract",
    addDOC = FALSE,
    plotName = "Plot-Footprints-Subtract.sw10.no_ribbon.pdf",
    pal = mycolor$Epi_Group,
    smoothWindow = 10
)

rm(motifPositions, seFoot)

## 5.5. compare significant TF with Absea antiboty ----
# Up peaks enriched
TF.sig.list <- lapply(homer.res, function(x) {
    x <- x[x$Diff == "up", ]
    return(x$TF)
})
TF.sig.list$Absea <- read.table("E:/LabWork/Project/CRC_NGS_ATAC/Absea.TF.txt")$V1

pdf("TF_motif/Upset.Absea.TF.sig.pdf", 6, 4.5)
upset(fromList(TF.sig.list), nsets = 4, order.by = "freq")
dev.off()

TF.sig.list$Absea %>%
    intersect(TF.sig.list$Group_1) %>%
    intersect(TF.sig.list$Group_2) %>%
    intersect(TF.sig.list$Common)

# up regulated TFs
Up.gene <- Gepia.DEG %>%
    filter(adjp < .01 & `Log2(Fold Change)` > 1) %>%
    pull(`Gene Symbol`)
v <- Vennerable::Venn(list(
    "Absea" = TF.Absea,
    "Up" = Up.gene
))
pdf("TF_motif/Venn.Absea.Up.pdf", 6, 4.5)
plot(v, doWeights = FALSE, show = list(Faces = FALSE))
dev.off()
intersect(TF.Absea, Up.gene)

# down regulated TFs
TF.Absea <- read.table("E:/LabWork/Project/CRC_NGS_ATAC/Absea.TF.txt")$V1
Gepia.DEG <- data.table::fread("E:/LabWork/Project/CRC_NGS_ATAC/GEPIA.COAD.DEG.txt", header = TRUE) %>%
    as.data.frame()
colnames(Gepia.DEG)

Down.gene <- Gepia.DEG %>%
    filter(adjp < .01 & `Log2(Fold Change)` < -1) %>%
    pull(`Gene Symbol`)
v <- Vennerable::Venn(list(
    "Absea" = TF.Absea,
    "Down" = Down.gene
))
pdf("TF_motif/Venn.Absea.Down.pdf", 6, 4.5)
plot(v, doWeights = FALSE, show = list(Faces = FALSE))
dev.off()
intersect(TF.Absea, Down.gene)

# 6. Cluster-level analysis ----
dir.create("Cluster_level")
cluster.info$iCMS <- ifelse(cluster.info$Epi_Group == "Group_1", "iCMS3", "iCMS2")
mycolor$iCMS <- c("iCMS2" = "#283891", "iCMS3" = "#62b7e6")

## 6.1 diff peaks in each cluster & run Homers ----
table(proj_Epi$Clusters)
dir.create("diff_peak_cluster")
dir.create("diff_peak_cluster/bed")
dir.create("diff_peak_cluster/down")

for (one in rownames(cluster.info)) {
    marker.peak.one <- getMarkerFeatures(
        ArchRProj = proj_Epi,
        useMatrix = "PeakMatrix",
        groupBy = "Clusters",
        testMethod = "wilcoxon",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        useGroups = one,
        bgdGroups = c("C3", "C4")
    )

    # up peak data frame
    marker.peak.one.df <- getMarkers(marker.peak.one,
        cutOff = "FDR <= 0.01 & Log2FC >= 1"
    )[[1]] %>% as.data.frame()
    write.table(
        marker.peak.one.df,
        paste0("diff_peak_cluster/marker.peak.tumor.", one, ".tsv"),
        col.names = TRUE,
        sep = "\t", quote = FALSE, row.names = FALSE
    )

    # up peak bed file
    marker.peak.one.df <- marker.peak.one.df %>%
        select(seqnames, start, end) %>%
        mutate(
            id = paste(seqnames, start, end, sep = "_"),
            length = 500,
            strand = "."
        )
    write.table(
        marker.peak.one.df,
        paste0("diff_peak_cluster/bed/marker.peak.tumor.", one, ".bed"),
        col.names = FALSE,
        sep = "\t", quote = FALSE, row.names = FALSE
    )

    # down peak data frame
    marker.peak.one.df <- getMarkers(marker.peak.one,
        cutOff = "FDR <= 0.01 & Log2FC <= -1"
    )[[1]] %>% as.data.frame()
    write.table(
        marker.peak.one.df,
        paste0("diff_peak_cluster/down/marker.peak.tumor.", one, ".tsv"),
        col.names = TRUE,
        sep = "\t", quote = FALSE, row.names = FALSE
    )
}
rm(one, marker.peak.one.df, marker.peak.one)

# find motifs in each patient
# sh Run.Homer.Motif.sh

peaks.clusters.up <- sapply(rownames(cluster.info), function(one) {
    temp <- read.table(paste0("diff_peak_cluster/marker.peak.tumor.", one, ".tsv"),
        header = TRUE, sep = "\t", stringsAsFactors = FALSE
    )
    return(paste(temp$seqnames, temp$start, temp$end, sep = "_"))
})
peaks.clusters.down <- sapply(rownames(cluster.info), function(one) {
    temp <- read.table(paste0("diff_peak_cluster/down/marker.peak.tumor.", one, ".tsv"),
        header = TRUE, sep = "\t", stringsAsFactors = FALSE
    )
    return(paste(temp$seqnames, temp$start, temp$end, sep = "_"))
})
identical(names(peaks.clusters.up), rownames(cluster.info))
identical(names(peaks.clusters.down), rownames(cluster.info))

cluster.info$nPeak_Up <- sapply(peaks.clusters.up, length)
cluster.info$nPeak_Down <- sapply(peaks.clusters.down, length)
cluster.info$Cluster <- rownames(cluster.info)

pdf("Cluster_level/Barplot.diff_peak.pdf", 7, 3)
cluster.info %>%
    arrange(-nPeak_Up) %>%
    mutate(Cluster = factor(Cluster, levels = Cluster)) %>%
    mutate(label = as.character(round(nPeak_Up / nPeak_Down, 2))) %>%
    ggplot() +
    geom_bar(aes(x = Cluster, y = nPeak_Up, fill = "Gained"),
        stat = "identity"
    ) +
    geom_bar(aes(x = Cluster, y = -nPeak_Down, fill = "Lost"),
        stat = "identity"
    ) +
    # geom_text(aes(x = Cluster, y = nPeak_Up + 3000, label = label), size = 2) +
    facet_grid(cols = vars(iCMS), scales = "free", space = "free") +
    scale_fill_manual(values = c("Gained" = "#cd2525", "Lost" = "#1774cd")) +
    ylab("Number of peaks") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()

## 6.2 compare cluster peaks ----
Jaccard_Sim <- function(x, y) {
    x <- unique(x)
    y <- unique(y)
    return(length(intersect(x, y)) / length(union(x, y)))
}

cluster.info2 <- readRDS("../04.Epi_CIMP/cluster.info.rds")
identical(rownames(cluster.info), rownames(cluster.info2))
diag(plot.data) <- NA
anno.col <- cluster.info %>%
    select(c("iCMS")) %>%
    mutate("CIMP_Group" = cluster.info2$CIMP_Group)
rm(cluster.info2)

colnames(anno.col) <- names(mycolor)[c(1:4)]
mycolor$iCMS <- c("iCMS2" = "#283891", "iCMS3" = "#62b7e6")
mycolor$CIMP_Group <- c("CIMP_High" = "#f76960", "CIMP_Low" = "#fbe625", "CIMP_Negative" = "#89dc0e")

plot.data <- matrix(0, nrow = length(peaks.clusters.up), ncol = length(peaks.clusters.up))
rownames(plot.data) <- colnames(plot.data) <- rownames(cluster.info)

for (i in seq_len(length(peaks.clusters.up))) {
    for (j in seq_len(length(peaks.clusters.up))) {
        plot.data[i, j] <- Jaccard_Sim(peaks.clusters.up[[i]], peaks.clusters.up[[j]])
    }
}
diag(plot.data) <- NA

rownames(plot.data) <- colnames(plot.data) <- cluster.rename[colnames(plot.data), ]$manual
rownames(anno.col) <- cluster.rename[rownames(anno.col), ]$manual

pdf("Cluster_level/Heatmap.Jaccard.peak.up.pdf", 7, 5)
pheatmap::pheatmap(plot.data,
    annotation_col = anno.col,
    annotation_row = anno.col,
    cluster_rows = TRUE, cluster_cols = TRUE,
    annotation_colors = mycolor[colnames(anno.col)],
    cluster_method = "ward.D",
    border_color = NA,
)
dev.off()

plot.data <- matrix(0, nrow = length(peaks.clusters.up), ncol = length(peaks.clusters.up))
rownames(plot.data) <- colnames(plot.data) <- rownames(cluster.info)

for (i in seq_len(length(peaks.clusters.down))) {
    for (j in seq_len(length(peaks.clusters.down))) {
        plot.data[i, j] <- Jaccard_Sim(peaks.clusters.down[[i]], peaks.clusters.down[[j]])
    }
}
diag(plot.data) <- NA

rownames(plot.data) <- colnames(plot.data) <- cluster.rename[colnames(plot.data), ]$manual

pdf("Cluster_level/Heatmap.Jaccard.peak.down.pdf", 7, 5)
pheatmap::pheatmap(plot.data,
    annotation_col = anno.col,
    annotation_row = anno.col,
    cluster_rows = TRUE, cluster_cols = TRUE,
    annotation_colors = mycolor[colnames(anno.col)],
    cluster_method = "ward.D",
    border_color = NA,
)
dev.off()

iCMS2.up <- union(marker.peak.list$Group_1, marker.peak.list$Both) %>%
    intersect(do.call(c, peaks.clusters.up[cluster.info$iCMS == "iCMS2"]))
temp <- table(do.call(c, peaks.clusters.up[cluster.info$iCMS == "iCMS2"]))[iCMS2.up]
rm(iCMS2.up, temp, peaks.clusters.up, peaks.clusters.down, plot.data, anno.col)

## 6.3 cluster level heatmap ----
sePeaks <- readRDS("sePeaks.cluster.rds")
cluster.info$iCMS <- ifelse(cluster.info$Epi_Group == "Group_1", "iCMS3", "iCMS2")

PeakMatrix <- sePeaks@assays@data$PeakMatrix

rownames(PeakMatrix) <- sePeaks@elementMetadata %>%
    as.data.frame() %>%
    mutate(name = paste(seqnames, start, end, sep = "_")) %>%
    pull(name)

plot.data <- PeakMatrix[
    c(marker.peak.list$Group_2, marker.peak.list$Both, marker.peak.list$Group_1),
    c("C3", "C4", rownames(cluster.info))
]
plot.data <- apply(plot.data, 1, function(x) {
    x <- (x - mean(x)) / sd(x)
    return(x)
}) %>% t()
plot.data[plot.data > 1.5] <- 1.5
plot.data[plot.data < (-1.5)] <- -1.5

names(marker.peak.list)
anno.row <- data.frame(
    row.names = rownames(plot.data),
    Group = factor(rep(c("iCMS2", "Both", "iCMS3"),
        times = c(length(marker.peak.list$Group_2), length(marker.peak.list$Both), length(marker.peak.list$Group_1))
    ), levels = c("iCMS2", "Both", "iCMS3")),
    Peak = "Up"
)
# anno.row <- anno.row[order(anno.row$Group, rowMeans(plot.data[, 18:27]) - rowMeans(plot.data[, 3:17])), ]
anno.row <- anno.row[order(anno.row$Group, rnorm(nrow(plot.data))), ]

ann.col <- data.frame(
    row.names = colnames(plot.data),
    iCMS = factor(c("Normal", "Normal", cluster.info$iCMS), levels = c("Normal", "iCMS2", "iCMS3")),
    temp = "temp"
)

# ann.col <- ann.col[order(ann.col$iCMS, rnorm(ncol(plot.data))), ]
ann.col <- ann.col[order(ann.col$iCMS, colSums(plot.data)), ]
plot.data <- plot.data[row.names(anno.row), rownames(ann.col)]

png("Cluster_level/Heatmap.marker.tmuor.cluster.Up.png", 400, 600)
pheatmap(plot.data,
    # color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")[2:8]))(100),
    scale = "none",
    annotation_row = anno.row %>% select(Group),
    annotation_col = ann.col %>% select(iCMS),
    annotation_colors = list(
        Group = c(iCMS3 = "#62b7e6", iCMS2 = "#283891", Both = "#86d786"),
        iCMS = c("Normal" = "#208a42", "iCMS2" = "#283891", "iCMS3" = "#62b7e6")
    ),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    gaps_col = c(2, 17),
    gaps_row = c(length(marker.peak.list$Group_2), length(marker.peak.list$Group_2) + length(marker.peak.list$Both))
)
dev.off()
rm(plot.data, anno.row, ann.col, PeakMatrix)

## 6.4 cluster-wise TF enrichment ----
homer.res.cluster <- list()

sapply(rownames(cluster.info), function(one) {
    homer.res <- homer.parser(paste0("homer/Cluster/", one, "/knownResults.txt"))
    homer.res$Cluster <- one
    homer.res.cluster[[one]] <<- homer.res
    return(1)
})
sapply(homer.res.cluster, dim)

homer.res.cluster <- do.call(rbind, homer.res.cluster)
table(homer.res.cluster$Cluster)

# iCMS TFs
TF.selected <- c("HNF4A", "PPARA", "HNF1", "TR4", "CDX2", "SOX4", "SOX2", "MAFB", "MAFK", "FOXA2", "FOXA3")

plot.data <- homer.res.cluster %>%
    filter(TF %in% TF.selected) %>%
    mutate(
        iCMS = cluster.info[.$Cluster, "iCMS"],
        CIMP_Group = cluster.info[.$Cluster, "CIMP_Group"],
        MSI_Status = cluster.info[.$Cluster, "MSI_Status_Major"]
    )

plot.data$TF <- factor(plot.data$TF, levels = TF.selected)
plot.data$Group <- ifelse(plot.data$TF %in% c("HNF4A", "PPARA", "CDX2", "HNF1", "TR4"), "iCMS2 TF", "iCMS3 TF")
if (max(plot.data$log.p.value) > 500) {
    plot.data[plot.data$log.p.value > 500, ]$log.p.value <- 500
}

cluster.rename <- read.table("../01.All_scCNV/cluster_rename.txt", header = TRUE)
rownames(cluster.rename) <- cluster.rename$Cluster
plot.data$Cluster_rename <- cluster.rename[plot.data$Cluster, "manual"]
plot.data$Cluster_rename <- factor(plot.data$Cluster_rename,
    levels = paste("C", 5:29, sep = "") %>% rev()
)
pdf("Cluster_level/Dot.TF.iCMS.cluster.rename.pdf", 5, 6)
ggplot(plot.data) +
    geom_point(aes(
        x = TF, y = Cluster_rename,
        fill = log.p.value, size = Log2_Enrichment
    ), pch = 21) +
    scale_fill_viridis_c() +
    facet_grid(
        rows = vars(iCMS),
        cols = vars(Group),
        scales = "free", space = "free"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

rm(p, plot.data, i, gene.selected)

# 7. TF and intra-subtype heterogeneity ----

## 7.1. get data: diff peak & motif match ----
# significant peaks of each clusters
peaks.clusters <- list()
sapply(rownames(cluster.info), function(one) {
    temp <- read.table(paste0("diff_peak_cluster/marker.peak.tumor.", one, ".tsv"), header = TRUE) %>%
        mutate(
            id = paste(seqnames, start, end, sep = "_")
        )
    peaks.clusters[[one]] <<- temp$id
    return(1)
})
sapply(peaks.clusters, length)

# motif match matrix
motif.match <- getMatches(ArchRProj = proj_Epi, name = "Motif")
temp <- rowRanges(motif.match)
motif.match <- motif.match@assays@data$matches

rownames(motif.match) <- paste(seqnames(temp), start(temp), end(temp), sep = "_")
colnames(motif.match) <- gsub("_.+?$", "", colnames(motif.match))
motif.match[1:5, 1:5]
rm(temp)

## 7.2. plot diff peak with motif in each cluster ----
peaks.common <- list()
sapply(peaks.common, length)
# Group2 TFs: HNF4A, PPARA
cluster.selected <- cluster.info %>%
    filter(Epi_Group == "Group_2") %>%
    rownames()

# HNF4A
peak.selected <- motif.match[, "HNF4A"] %>%
    .[.] %>%
    names() %>%
    intersect(., unlist(peaks.clusters[cluster.selected]) %>% unique())
length(peak.selected) # 20059 HNF4A peaks

plot.data <- matrix(0, nrow = length(peak.selected), ncol = length(cluster.selected))
rownames(plot.data) <- peak.selected
colnames(plot.data) <- cluster.selected
for (one in cluster.selected) {
    plot.data[peak.selected %in% peaks.clusters[[one]], one] <- 1
}

pdf("Cluster_level/Heatmap.HNF4A.clusterPeaks1.pdf", 10, 6)
p <- pheatmap(
    t(plot.data)[c(
        "C18", "C2", "C16", "C13", "C22", "C17", "C20",
        "C12", "C21", "C19", "C7", "C23", "C14", "C15", "C8"
    ), ],
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_row = cluster.info %>% select(c("Epi_Group")),
    annotation_color = mycolor,
    color = colorRampPalette(c("gray80", "#cd2525"))(2),
    clustering_method = "ward.D2",
    cutree_cols = 17
)
dev.off()

temp <- cutree(p$tree_col, 17)
rowSums(plot.data[names(temp), ]) %>%
    aggregate(., by = list(temp), FUN = mean) # 4 and 15
peaks.common[["HNF4A"]] <- names(temp)[temp %in% c(4, 15)]

anno.col <- data.frame(
    row.names = rownames(plot.data),
    "con" = ifelse(rownames(plot.data) %in% peaks.common[["HNF4A"]], "Consensus", "Specific")
)

pdf("Cluster_level/Heatmap.HNF4A.clusterPeaks.pdf", 10, 6)
p <- pheatmap(
    t(plot.data)[c(
        "C18", "C2", "C16", "C13", "C22", "C17", "C20",
        "C12", "C21", "C19", "C7", "C23", "C14", "C15", "C8"
    ), ],
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_row = cluster.info %>% select(c("Epi_Group")),
    annotation_col = anno.col,
    annotation_color = c(mycolor, list(con = c("Consensus" = "#86d786", "Specific" = "#f6be43"))),
    color = colorRampPalette(c("gray80", "#cd2525"))(2),
    clustering_method = "ward.D2",
    cutree_cols = 17
)
dev.off()

# PPARA
peak.selected <- motif.match[, "PPARA"] %>%
    .[.] %>%
    names() %>%
    intersect(., unlist(peaks.clusters[cluster.selected]) %>% unique())
length(peak.selected) # 9720 PPARA peaks

plot.data <- matrix(0, nrow = length(peak.selected), ncol = length(cluster.selected))
rownames(plot.data) <- peak.selected
colnames(plot.data) <- cluster.selected
for (one in cluster.selected) {
    plot.data[peak.selected %in% peaks.clusters[[one]], one] <- 1
}

pdf("Cluster_level/Heatmap.PPARA.clusterPeaks.pdf", 10, 6)
p <- pheatmap(
    t(plot.data)[c(
        "C2", "C16", "C13", "C12", "C14", "C22", "C18", "C23",
        "C15", "C19", "C17", "C20", "C21", "C7", "C8"
    ), ],
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_row = cluster.info %>% select(c("Epi_Group")),
    annotation_color = mycolor,
    color = colorRampPalette(c("gray80", "#cd2525"))(2),
    clustering_method = "ward.D2",
    cutree_cols = 15
)
dev.off()

temp <- cutree(p$tree_col, 15)
rowSums(plot.data[names(temp), ]) %>%
    aggregate(., by = list(temp), FUN = mean) # 1 and 4
peaks.common[["PPARA"]] <- names(temp)[temp %in% c(1, 4)]

anno.col <- data.frame(
    row.names = rownames(plot.data),
    "con" = ifelse(rownames(plot.data) %in% peaks.common[["PPARA"]], "Consensus", "Specific")
)

pdf("Cluster_level/Heatmap.PPARA.clusterPeaks.pdf", 10, 6)
p <- pheatmap(
    t(plot.data)[c(
        "C2", "C16", "C13", "C12", "C14", "C22", "C18", "C23",
        "C15", "C19", "C17", "C20", "C21", "C7", "C8"
    ), ],
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_row = cluster.info %>% select(c("Epi_Group")),
    annotation_col = anno.col,
    annotation_color = c(mycolor, list(con = c("Consensus" = "#86d786", "Specific" = "#f6be43"))),
    color = colorRampPalette(c("gray80", "#cd2525"))(2),
    clustering_method = "ward.D2",
    cutree_cols = 15
)
dev.off()

# Group1 TFs: SOX2, FOXA3
# SOX2
cluster.selected <- cluster.info %>%
    filter(Epi_Group == "Group_1") %>%
    rownames()

peak.selected <- motif.match[, "SOX2"] %>%
    .[.] %>%
    names() %>%
    intersect(., unlist(peaks.clusters[cluster.selected]) %>% unique())
length(peak.selected) # 6048 SOX2 peaks

plot.data <- matrix(0, nrow = length(peak.selected), ncol = length(cluster.selected))
rownames(plot.data) <- peak.selected
colnames(plot.data) <- cluster.selected
for (one in cluster.selected) {
    plot.data[peak.selected %in% peaks.clusters[[one]], one] <- 1
}

pdf("Cluster_level/Heatmap.SOX2.clusterPeaks.pdf", 10, 4.5)
p <- pheatmap(
    t(plot.data)[c(
        "C5", "C29", "C11", "C10", "C24", "C6", "C26", "C25",
        "C27", "C1"
    ), ],
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_row = cluster.info %>% select(c("Epi_Group")),
    annotation_color = mycolor,
    color = colorRampPalette(c("gray80", "#cd2525"))(2),
    clustering_method = "ward.D2",
    cutree_cols = 11
)
dev.off()

temp <- cutree(p$tree_col, 11)
rowSums(plot.data[names(temp), ]) %>%
    aggregate(., by = list(temp), FUN = mean) # 1
peaks.common[["SOX2"]] <- names(temp)[temp %in% c(1)]

anno.col <- data.frame(
    row.names = rownames(plot.data),
    "con" = ifelse(rownames(plot.data) %in% peaks.common[["SOX2"]], "Consensus", "Specific")
)

pdf("Cluster_level/Heatmap.SOX2.clusterPeaks.pdf", 10, 4.5)
p <- pheatmap(
    t(plot.data)[c(
        "C5", "C29", "C11", "C10", "C24", "C6", "C26", "C25",
        "C27", "C1"
    ), ],
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_row = cluster.info %>% select(c("Epi_Group")),
    annotation_col = anno.col,
    annotation_color = c(mycolor, list(con = c("Consensus" = "#86d786", "Specific" = "#f6be43"))),
    color = colorRampPalette(c("gray80", "#cd2525"))(2),
    clustering_method = "ward.D2",
    cutree_cols = 11
)
dev.off()

# FOXA3
peak.selected <- motif.match[, "FOXA3"] %>%
    .[.] %>%
    names() %>%
    intersect(., unlist(peaks.clusters[cluster.selected]) %>% unique())
length(peak.selected) # 18275 FOXA3 peaks

plot.data <- matrix(0, nrow = length(peak.selected), ncol = length(cluster.selected))
rownames(plot.data) <- peak.selected
colnames(plot.data) <- cluster.selected
for (one in cluster.selected) {
    plot.data[peak.selected %in% peaks.clusters[[one]], one] <- 1
}

pdf("Cluster_level/Heatmap.FOXA3.clusterPeaks.pdf", 10, 4.5)
p <- pheatmap(
    t(plot.data)[c(
        "C5", "C29", "C11", "C10", "C24", "C6", "C26", "C25",
        "C27", "C1"
    ), ],
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_row = cluster.info %>% select(c("Epi_Group")),
    annotation_color = mycolor,
    color = colorRampPalette(c("gray80", "#cd2525"))(2),
    clustering_method = "ward.D2",
    cutree_cols = 11
)
dev.off()

temp <- cutree(p$tree_col, 11)
rowSums(plot.data[names(temp), ]) %>%
    aggregate(., by = list(temp), FUN = mean) # 8
peaks.common[["FOXA3"]] <- names(temp)[temp %in% c(8)]

anno.col <- data.frame(
    row.names = rownames(plot.data),
    "con" = ifelse(rownames(plot.data) %in% peaks.common[["FOXA3"]], "Consensus", "Specific")
)

pdf("Cluster_level/Heatmap.FOXA3.clusterPeaks.pdf", 10, 4.5)
p <- pheatmap(
    t(plot.data)[c(
        "C5", "C24", "C6", "C25", "C29", "C1", "C10", "C11", "C26",
        "C27"
    ), ],
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_row = cluster.info %>% select(c("Epi_Group")),
    annotation_col = anno.col,
    annotation_color = c(mycolor, list(con = c("Consensus" = "#86d786", "Specific" = "#f6be43"))),
    color = colorRampPalette(c("gray80", "#cd2525"))(2),
    clustering_method = "ward.D2",
    cutree_cols = 11
)
dev.off()

rm(p, peak.selected, cluster.selected, plot.data, one, temp, anno.col)

## 7.3. identify iCMS-consensus TFs ----
homer.res.cluster <- list()

sapply(rownames(cluster.info), function(one) {
    homer.res <- homer.parser(paste0("homer/Cluster/", one, "/knownResults.txt"), log2.enrichment = 0.8)
    homer.res$Cluster <- one
    homer.res.cluster[[one]] <<- homer.res
    return(1)
})
sapply(homer.res.cluster, dim)
# sapply(homer.res.cluster, function(x) {
#     return(x[x$TF == "PPARA", c("p.value", "Log2_Enrichment")])
# })

TF.sig.cluster <- sapply(homer.res.cluster, function(x) {
    x <- x[x$Diff == "up", ]
    return(x$TF)
})

cluster.selected <- cluster.info %>%
    filter(iCMS == "iCMS2") %>%
    rownames()
plot.data <- table(unlist(TF.sig.cluster[cluster.selected])) %>%
    as.data.frame() %>%
    mutate(iCMS = "iCMS2")
p1 <- ggplot(plot.data) +
    geom_bar(aes(x = Freq, fill = iCMS),
        stat = "count", show.legend = FALSE
    ) +
    facet_grid(cols = vars(iCMS), scales = "free", space = "free") +
        scale_fill_manual(values = c("iCMS2" = "#283891", "iCMS3" = "#62b7e6")) +
        theme_bw() +
        geom_vline(xintercept = 11.5, linetype = 2) +
        xlab("Number of clusters") +
        ylab("Number of TFs")
consensus.TF <- list()
consensus.TF$iCMS2 <- as.character(plot.data$Var1)[plot.data$Freq > 11.5]

cluster.selected <- cluster.info %>%
    filter(iCMS == "iCMS3") %>%
    rownames()
plot.data <- table(unlist(TF.sig.cluster[cluster.selected])) %>%
    as.data.frame() %>%
    mutate(iCMS = "iCMS3")
p2 <- ggplot(plot.data) +
    geom_bar(aes(x = Freq, fill = iCMS),
        stat = "count", show.legend = FALSE
    ) +
    facet_grid(cols = vars(iCMS), scales = "free", space = "free") +
    scale_fill_manual(values = c("iCMS2" = "#283891", "iCMS3" = "#62b7e6")) +
    theme_bw() +
    geom_vline(xintercept = 7.5, linetype = 2) +
        xlab("Number of clusters") +
        ylab("Number of TFs")

consensus.TF$iCMS3 <- as.character(plot.data$Var1)[plot.data$Freq > 7.5]

pdf("Cluster_level/Barplot.TF.iCMS.cluster.pdf", 7, 3)
patchwork::wrap_plots(p1, p2, ncol = 2)
dev.off()

pdf("Cluster_level/Venn.consensus.TF.iCMS.pdf", 5, 4)
plot(Vennerable::Venn(list(
    "iCMS2" = consensus.TF$iCMS2,
    "iCMS3" = consensus.TF$iCMS3
)), doWeights = FALSE, show = list(Faces = FALSE))
dev.off()
setdiff(consensus.TF$iCMS2, consensus.TF$iCMS3)
# "BORIS" "CDX2" "CDX4" "HNF1" "HNF1B" "HNF4A" "HOXB13" "JUND" "LEF1" "NUR77" "SP1" "TR4"
setdiff(consensus.TF$iCMS3, consensus.TF$iCMS2)
# "EKLF" "ELF3" "ELF5" "ETS1-DISTAL" "FOX:EBOX" "FOXA1" "FOXM1""GRHL2" "KLF3" "MAFB" "NF1:FOXA1"
intersect(consensus.TF$iCMS2, consensus.TF$iCMS3)

rm(p, plot.data, p1, p2, cluster.selected, one, temp)

## 7.4. target of consensus TFs ----
consensus.TF.selected <-
    sapply(consensus.TF, function(x) {
        x[x %in% colnames(motif.match)]
    })
sapply(consensus.TF.selected, length)

# iCMS2 consensu TFs
cluster.selected <- cluster.info %>%
    filter(iCMS == "iCMS2") %>%
    rownames()
cTF.peak.stat <- data.frame(
    row.names = consensus.TF.selected$iCMS2,
    TF = consensus.TF.selected$iCMS2,
    "n_Peak" = 0,
    "mean_n_cluster" = 0,
    "n_cPeak" = 0,
    "frac_cPeak" = 0,
    "iCMS" = "iCMS2"
)
for (one in consensus.TF.selected$iCMS2) {
    peak.selected <- motif.match[, one] %>%
        .[.] %>%
        names() %>%
        intersect(., unlist(peaks.clusters[cluster.selected]) %>% unique())

    cTF.peak.stat[one, "n_Peak"] <- length(peak.selected)
    temp <- peaks.clusters[cluster.selected] %>%
        unlist() %>%
        .[. %in% peak.selected] %>%
        table() %>%
        sort(decreasing = TRUE)
    cTF.peak.stat[one, "mean_n_cluster"] <- mean(temp)
    cTF.peak.stat[one, "n_cPeak"] <- sum(temp > length(cluster.selected) * 0.4)
    cTF.peak.stat[one, "frac_cPeak"] <- cTF.peak.stat[one, "n_cPeak"] / length(peak.selected)
}
cTF.peak.stat.iCM2 <- cTF.peak.stat

# iCMS3 consensu TFs
cluster.selected <- cluster.info %>%
    filter(iCMS == "iCMS3") %>%
    rownames()
cTF.peak.stat <- data.frame(
    row.names = consensus.TF.selected$iCMS3,
    TF = consensus.TF.selected$iCMS3,
    "n_Peak" = 0,
    "mean_n_cluster" = 0,
    "n_cPeak" = 0,
    "frac_cPeak" = 0,
    "iCMS" = "iCMS3"
)

for (one in consensus.TF.selected$iCMS3){
    peak.selected <- motif.match[, one] %>%
        .[.] %>%
        names() %>%
        intersect(., unlist(peaks.clusters[cluster.selected]) %>% unique())

    cTF.peak.stat[one, "n_Peak"] <- length(peak.selected)
    temp <- peaks.clusters[cluster.selected] %>%
        unlist() %>%
        .[. %in% peak.selected] %>%
        table() %>%
        sort(decreasing = TRUE)
    cTF.peak.stat[one, "mean_n_cluster"] <- mean(temp)
    cTF.peak.stat[one, "n_cPeak"] <- sum(temp > length(cluster.selected) * 0.4)
    cTF.peak.stat[one, "frac_cPeak"] <- cTF.peak.stat[one, "n_cPeak"] / length(peak.selected)
}
cTF.peak.stat.iCM3 <- cTF.peak.stat
cTF.peak.stat <- rbind(cTF.peak.stat.iCM2, cTF.peak.stat.iCM3)
saveRDS(cTF.peak.stat, "cTF.peak.stat.rds")

pdf("Cluster_level/Scatter.frac_cTF.pdf", 9, 4)
ggplot(cTF.peak.stat, aes(x = n_Peak, y = frac_cPeak)) +
    geom_point(aes(color = mean_n_cluster), size = 2.5) +
    ggrepel::geom_text_repel(aes(label = TF), size = 3) +
    facet_grid(cols = vars(iCMS), scales = "free") +
    xlab("Number of gained targeting peaks") +
    ylab("Fraction of consensus targeting peaks") +
    scale_color_viridis_c() +
    theme_bw()
dev.off()

# 8. Test 6-rank NMF ----
dir.create("NMF_r6")
## 8.1 run NMF with 6 clusters ----
sePeaks <- readRDS("sePeaks.cluster.rds")
NMF.res.6 <- nmf(
    sePeaks@assays@data$PeakMatrix[feature.selcted, sample.selected],
    rank = 6,
    nrun = 200
)
rm(sePeaks)

group.res.6 <- predict(NMF.res.6, what = "consensus") %>%
    as.numeric()
group.res.6 <- group.res.6 %>%
    paste("Group_", ., sep = "")
names(group.res.6) <- sample.selected

## 8.2. visualize NMF results ----
# conseusus matrix
con.mat <- NMF.res.6@consensus
sil <- silhouette(NMF.res.6, what = "consensus")

cluster.info2 <- read.table("../04.Epi_CIMP/cluster.info.csv", header = TRUE, sep = ",", row.names = 1)
identical(rownames(cluster.info2), rownames(cluster.info))
cluster.info$CIMP_Group <- cluster.info2[rownames(cluster.info), "CIMP_Group"]
rm(cluster.info2)

anno.col <- cluster.info %>%
    select(c("Gender_Major", "MSI_Status_Major", "Side_Major", "iCMS", "CIMP_Group"))
colnames(anno.col) <- c(names(mycolor)[c(1:3, 6)], "CIMP")

anno.col$silhouette <- sil[, 3][rownames(anno.col)]
anno.col$NMF_r6 <- group.res.6[rownames(anno.col)]
mycolor$CIMP <- c("CIMP_Negative" = "#89dc0e", "CIMP_Low" = "#fbe625", "CIMP_High" = "#f76960")
mycolor$NMF_r6 <- ArchR::paletteDiscrete(anno.col$NMF_r6)

cluster.rename <- read.table("../01.All_scCNV/cluster_rename.txt", header = TRUE)
rownames(cluster.rename) <- cluster.rename$Cluster

rownames(anno.col) <- cluster.rename[rownames(anno.col), "manual"]
colnames(con.mat) <- cluster.rename[colnames(con.mat), "manual"]
rownames(con.mat) <- cluster.rename[rownames(con.mat), "manual"]

pdf("NMF_r6/NMF.consensus.clusters.pdf", 7, 6)
pheatmap::pheatmap(con.mat,
    annotation_col = anno.col[, c(6, 2, 3, 5, 4, 7)],
    annotation_colors = mycolor[colnames(anno.col)],
    border_color = NA,
    cutree_rows = 6,
    cutree_cols = 6
)
dev.off()

# coefficient matrix
coefmap(NMF.res.6)
coef.mat <- coef(NMF.res.6)
rownames(coef.mat) <- paste("Basis", 1:6, sep = "")
write.csv(coef.mat, "NMF_r6/NMF.coefficient.clusters.csv")
colnames(coef.mat) <- cluster.rename[colnames(coef.mat), "manual"]

pdf("NMF_r6/NMF.coefficient.clusters.pdf", 6, 3)
pheatmap(coef.mat,
    annotation_col = anno.col[c(4, 5, 7)],
    annotation_colors = mycolor[colnames(anno.col)],
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    border_color = NA,
    clustering_method = "ward.D2",
)
dev.off()

rm(con.mat, coef.mat, sil, p, anno.col)


gc()
write.csv(cluster.info, "cluster.info.csv")
proj_Epi <- saveArchRProject(ArchRProj = proj_Epi, load = TRUE)
save.image("Epi_Molecular_Subtype.RData")