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
    size = 0.2, plotAs = "points"
)

pdf("UAMP.Epi.Epi_type.pdf", 7, 7)
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
sample.info.epi <- proj_Epi@cellColData %>% as.data.frame()
rm(temp)

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
    size = 0.2, plotAs = "points"
)
pdf("UAMP.Epi.Gender.pdf", 7, 7)
plot(p)
dev.off()

p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    name = "MSI_Status", embedding = "UMAP",
    pal = mycolor$MSI,
    size = 0.2, plotAs = "points"
)
pdf("UAMP.Epi.MSI_Status.pdf", 7, 7)
plot(p)
dev.off()

p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    name = "Side", embedding = "UMAP",
    pal = mycolor$Side,
    size = 0.2, plotAs = "points"
)
pdf("UAMP.Epi.Side.pdf", 7, 7)
plot(p)
dev.off()

## 1.3 pseudo-bulk repicates ----
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

# 2. clustering in peak space ----
## 2.1. feature selection by DIP test ----
dip.pvalue <- apply(sePeaks@assays@data$PeakMatrix, 1, function(x) {
    return(diptest::dip.test(x)$p.value)
})
quantile(dip.pvalue)
table(dip.pvalue < 0.05)
feature.selcted <- dip.pvalue[dip.pvalue < 0.05] %>% names()

## 2.2. hclust of malignant clusters ----
input.data <- sePeaks@assays@data$PeakMatrix[feature.selcted, sample.selected]
hc.tumor <- hclust(
    d = dist(t(input.data)),
    method = "ward.D2"
)
rm(input.data)

plot(hc.tumor)
rect.hclust(
    tree = hc.tumor,
    k = 2, which = 1:2, border = c("#62b7e6", "#283891"),
    cluster = cutree(tree = hc.tumor, k = 2)
)
group.res <- cutree(tree = hc.tumor, k = 2)

proj_Epi$Epi_Group <- proj_Epi$Epi_type
proj_Epi$Epi_Group[proj_Epi$Epi_Group == "Malignant"] <-
    paste("Group_",
        (group.res[proj_Epi$Clusters[proj_Epi$Epi_Group == "Malignant"]]),
        sep = ""
    )
rm(group.res)
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

## 2.3. plot hclust ----
cluster.info <- as.data.frame(cutree(hc.tumor, 2))
identical(rownames(cluster.info), sample.selected)
colnames(cluster.info) <- "Group"
cluster.info$Group <- paste("Group_", cluster.info$Group, sep = "")
cluster.info$Count <- table(sample.info.epi$Clusters)[rownames(cluster.info)]

# annotate clusters
for (one in c("Patient", "Gender", "MSI_Status", "Side")){
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
    cluster.info[, paste0(one, "_Major_Pct")] <- cluster.info[, paste0(one, "_Major_Count")] / cluster.info$Count
}
rm(temp, one)

dend <- as.dendrogram(hc.tumor)
dend.color <- data.frame(
    row.names = rownames(cluster.info),
    "Gender" = mycolor$Gender[cluster.info$Gender_Major],
    "Side" = mycolor$Side[cluster.info$Side_Major],
    "MSI" = mycolor$MSI[cluster.info$MSI_Status_Major],
    "Group" = mycolor$Group[cluster.info$Group]
)

temp <- mycolor[c(4, 2, 3, 1)] %>% unlist()
names(temp) <- gsub("^.+?\\.", "", names(temp))
temp <- temp[c(1, 3, 6, 8, 5, 2, 4, 7, 9)]

pdf("Dendrogram.Tumor.clusters.pdf", 6, 4)
plot(dend, main = "Hierarchical clustering of malignant clusters")
dendextend::colored_bars(
    colors = dend.color, dend = dend,
    y_scale = 10, y_shift = -8
)
legend("topright",
    legend = names(temp), col = temp,
    pch = 15, bty = "n",
    pt.cex = 1.2, cex = 0.8,
    ncol = 2
)
dev.off()
rm(dend, dend.color, temp, p)

save.image("Epi_Molecular_Subtype.RData")

# 4. cluster in CNV space ----
load("../01.All_scCNV/CRC_CNV.rda")

CNV.FC <- CNV.FC[rownames(sample.info.epi), ]

CNV.FC <- aggregate(CNV.FC,
    by = list(sample.info.epi$Clusters),
    FUN = mean
)
rownames(CNV.FC) <- CNV.FC$Group.1
CNV.FC$Group.1 <- NULL

pdf("Heatmap.CNV.merge2.pdf", 8, 6)
pheatmap::pheatmap(CNV.FC[sample.selected, ],
    color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
    cluster_rows = TRUE, cluster_cols = FALSE,
    clustering_method = "ward.D2",
    annotation_col = anno.col,
    annotation_row = clusters.res,
    annotation_colors = anno.color,
    annotation_legend = FALSE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    cutree_rows = 2,
    gaps_col = which(!duplicated(anno.col$seqnames)) - 1
)
dev.off()

# hclust of malignant clusters
hc.tumor.CNV <- hclust(
    d = dist(CNV.FC[sample.selected, ]),
    method = "ward.D2"
)
pdf("Dendrogram.Tumor.CNV.clusters.pdf", 6, 5)
plot(hc.tumor.CNV, main = "Hierarchical clustering of malignant clusters")
dev.off()

rm(CNV.FC, sample.info, anno.col, anno.row, anno.color)

clusters.res$CNV_HC <- cutree(hc.tumor.CNV, 2)[sample.selected]

aricode::AMI(clusters.res$Peak_NMF, clusters.res$CNV_HC)
aricode::AMI(clusters.res$Peak_NMF, clusters.res$LSI_HC)
aricode::AMI(clusters.res$LSI_HC, clusters.res$CNV_HC)



hc.tumor.peak <- hclust(
    d = dist(t(peatmat[feature.selcted, sample.selected])),
    method = "ward.D2"
)
pdf("Dendrogram.Tumor.peak.clusters.pdf", 6, 5)
plot(hc.tumor.peak)
dev.off()

cutree(hc.tumor, 2)[sample.selected]
cutree(hc.tumor.peak, 2)[sample.selected]
aricode::AMI(
    cutree(hc.tumor, 2)[sample.selected],
    cutree(hc.tumor.peak, 2)[sample.selected]
)
