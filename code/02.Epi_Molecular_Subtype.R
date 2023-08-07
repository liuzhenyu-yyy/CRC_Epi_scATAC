setwd("E:/LabWork/Project/CRC_NGS_ATAC/CRC_Epi_scATAC/Results/02.Epi_Molecular_Subtype")
source("../../code/00.Requirements.R")

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

getGroupBW(
    ArchRProj = proj_Epi,
    groupBy = "Clusters",
    normMethod = "ReadsInTSS",
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

dim(sePeaks)
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
    size = 0.2, plotAs = "points"
)
pdf("UAMP.Epi.Epi_Group.pdf", 7, 7)
plot(p)
dev.off()

## 2.3. visualize NMF results ----
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
        theme(axis.title = element_blank())

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
        theme(axis.title = element_blank())

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

pdf("Dot.Module.iCMS.pdf", 9, 4)
p1 <- ggplot(plot.data) +
    geom_point(aes(x = Module.iCMS2, y = Module.iCMS3, color = Density),
        show.legend = FALSE) +
    scale_color_viridis_c() +
    theme_classic()
p2 <- ggplot(plot.data) +
    geom_point(aes(x = Module.iCMS2, y = Module.iCMS3, color = Epi_Group)) +
    scale_color_manual(values = mycolor$Epi_Group) +
    theme_classic()
wrap_plots(p1, p2, ncol = 2)
dev.off()

pdf("Box.Module.iCMS.pdf", 5, 4)
p1 <- ggplot(plot.data, aes(x = Epi_Group, y = Module.iCMS2)) +
    geom_boxplot(aes(fill = Epi_Group),
        show.legend = FALSE
    ) +
    scale_fill_manual(values = mycolor$Epi_Group) +
    ggpubr::stat_compare_means(comparisons = list(c("Group_1", "Group_2"))) +
    theme_classic()
p2 <- ggplot(plot.data, aes(x = Epi_Group, y = Module.iCMS3)) +
    geom_boxplot(aes(fill = Epi_Group),
        show.legend = FALSE
    ) +
    scale_fill_manual(values = mycolor$Epi_Group) +
    ggpubr::stat_compare_means(comparisons = list(c("Group_1", "Group_2"))) +
    theme_classic()
wrap_plots(p1, p2, ncol = 2)
dev.off()

proj_Epi$Module.iCMS <- proj_Epi$Module.iCMS3 - proj_Epi$Module.iCMS2
p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    name = "Module.iCMS", embedding = "UMAP",
    size = 0.2, plotAs = "points"
)
pdf("UAMP.Epi.Module.iCMS.pdf", 7, 7)
plot(p)
dev.off()

rm(genes, plot.data, p, p1, p2)

proj_Epi <- saveArchRProject(ArchRProj = proj_Epi, load = TRUE)
save.image("Epi_Molecular_Subtype.RData")
