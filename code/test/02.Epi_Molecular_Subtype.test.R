###################################################
# test different methods for subtype classification
# @liuzhenyu-yyy
# 2023-08-02
###################################################

setwd("E:/LabWork/Project/CRC_NGS_ATAC/CRC_Epi_scATAC/Results/02.Epi_Molecular_Subtype")

source("../../code/00.Requirements.R")
load("Epi_Molecular_Subtype.RData")

# 1. prepare ArchR project ----
## 1.1 load saved project ----
proj_Epi <- loadArchRProject(project.dir.epi, force = TRUE)
proj_all <- loadArchRProject(project.dir.all, force = TRUE)
colnames(proj_all@cellColData)
proj_all$Clusters_type %>% table()
proj_Epi$Sample %>% table()

proj_Epi$Epi_type <- "Malignant"
proj_Epi$Epi_type[proj_Epi$Clusters %in% c("C3", "C4")] <- "Normal"
proj_Epi$Epi_type[proj_Epi$Clusters %in% c("C28", "C9")] <- "Adenoma"

table(proj_Epi$Epi_type, proj_Epi$new_location)

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

## 1.2 pseudo-bulk repicates ----
# on WSL
# setwd("/mnt/e/LabWork/Project/CRC_NGS_ATAC/CRC_Epi_scATAC/Results/02.Epi_Molecular_Subtype")
# library(ArchR)
# proj_Epi <- loadArchRProject("../../Project_Dir_Epi", force = TRUE)
# table(proj_Epi$Clusters)

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
sePeaks@colData
sePeaks@assays@data$PeakMatrix
rowData(sePeaks)

dim(sePeaks)
sample.selected <- setdiff(rownames(sePeaks@colData), c("C3", "C4", "C28", "C9"))

# 2. clustering in LSI space ----
# subtype by LSI coordnates
# merge LSI coordinates
LSI.coord <- proj_Epi@reducedDims$IterativeLSI_merge
names(LSI.coord)
LSI.coord <- LSI.coord$matSV %>% as.data.frame()

identical(rownames(LSI.coord), rownames(sample.info.epi))
LSI.coord <- aggregate(LSI.coord,
    by = list(sample.info.epi$Clusters),
    FUN = mean
)
rownames(LSI.coord) <- LSI.coord$Group.1
LSI.coord$Group.1 <- NULL

# hclust of malignant clusters
hc.tumor <- hclust(
    d = dist(LSI.coord[sample.selected, ]),
    method = "ward.D2"
)
rm(LSI.coord)
plot(hc.tumor)

pdf("Dendrogram.Tumor.clusters.pdf", 6, 5)
plot(hc.tumor, main = "Hierarchical clustering of malignant clusters")
dev.off()

# groups <- list(
#     "D1" = cluster.tumor[cluster.tumor == 2] %>% names(),
#     "D2" = cluster.tumor[cluster.tumor == 1] %>% names()
# )

clusters.res <- data.frame(
    row.names = sample.selected,
    LSI_HC = cutree(hc.tumor, 2)[sample.selected]
)

# 3. clustering in peak space ----
library(NMF)
library(diptest)

## 3.1. feature selection by DIP test ----
dip.pvalue <- apply(sePeaks@assays@data$PeakMatrix, 1, function(x) {
    return(dip.test(x)$p.value)
})
quantile(dip.pvalue)
dip.pvalue[1:5]
feature.selcted <- dip.pvalue[dip.pvalue < 0.05] %>% names()
table(dip.pvalue < 0.05)
input.data <- sePeaks@assays@data$PeakMatrix[feature.selcted, sample.selected]

## 3.2. run NMF ----
res.ranks <- nmf(peatmat[feature.selcted, sample.selected], 2:10, nrun = 50)

pdf("NMF.rank.res.pdf", 8, 6)
plot(res.ranks)
dev.off()
rm(res.ranks)

peak.num.2c <- nmf(input.data,
    rank = 2
)
predict(peak.num.2c)
clusters.res$Peak_NMF <- predict(peak.num.2c)[sample.selected] %>%
    as.numeric()
pdf("NMF.2clusters.pdf", 5,4 )
NMF::consensusmap(peak.num.2c)
dev.off()
aricode::AMI(
    clusters.res$LSI_HC %>% as.character(),
    clusters.res$Peak_NMF %>% as.character()
)

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

save.image("Epi_Molecular_Subtype.RData")

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
