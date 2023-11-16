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

group.res <- predict(NMF.res) %>%
    as.numeric()
group.res <- (3 - group.res) %>%
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
        stat = "identity", position =  position_fill(reverse = TRUE)
    ) +
    geom_text(aes(x = Epi_Group, y = Freq, label = Freq),
        position = position_fill(vjust = 0.5)
    ) +
    scale_fill_manual(values = mycolor$MSI) +
        coord_flip() +
        theme_classic() +
        theme(axis.title = element_blank(), axis.text.x = element_blank())

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
        theme(axis.title.y = element_blank(), axis.title.y = element_blank()) +
        ylab("Fraction")

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
pheatmap::pheatmap(CNV.FC,
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
plot.data <- sample.info.epi %>%
    filter(Epi_type == "Malignant") %>%
    select(matches("Module|Group"))

get_density <- function(x, y, nbins = 100, ...) {
    dens <- MASS::kde2d(x, y, n = nbins, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
}

plot.data$Density <- get_density(plot.data$Module.iCMS2, plot.data$Module.iCMS3)

pdf("Dot.Module.iCMS.pdf", 5, 3.7)
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

pdf("Box.Module.iCMS.pdf", 4, 3.5)
ggplot(plot.data, aes(x = Epi_Group, y = value)) +
    geom_violin(aes(fill = Epi_Group),
        show.legend = FALSE, size = 0.2
    ) +
    scale_fill_manual(values = mycolor$Epi_Group) +
        ggpubr::stat_compare_means(comparisons = list(c("Group_1", "Group_2"))) +
        facet_grid(cols = vars(variable), scales = "free_y") +
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
# 4. differential peak analysis ----
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

# 5. TF enrichment in marker peaks ----
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

TF.selected <- setdiff(TF.selected, "CTCF")
plot.data <- lapply(homer.res, function(x) {
    x <- x[x$TF %in% TF.selected, ]
    return(x)
}) %>% do.call(rbind, .)

plot.data$Group <- factor(plot.data$Group, levels = c("Group_2", "Common", "Group_1"))
plot.data$TF <- factor(plot.data$TF, levels = (TF.selected))

plot.data[plot.data$log.p.value > 1000, ]$log.p.value <- 1000
#plot.data[plot.data$Log2_Enrichment > 3, ]$Log2_Enrichment <- 3

pdf("TF_motif/Dot.motif.sig.selected1.pdf", 4, 4)
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

write.csv(cluster.info, "cluster.info.csv")
proj_Epi <- saveArchRProject(ArchRProj = proj_Epi, load = TRUE)
save.image("Epi_Molecular_Subtype.RData")
