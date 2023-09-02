setwd("E:/LabWork/Project/CRC_NGS_ATAC/CRC_Epi_scATAC/Results/04.Epi_CIMP")

source("../../code/00.Requirements.R")

# 1. Calculate ATAC-CGI matrix ----
## 1.1 load saved project ----
proj_Epi <- loadArchRProject(project.dir.epi, force = TRUE)
table(proj_Epi$Epi_type)
getAvailableMatrices(proj_Epi)

## 1.2. add CGI matrix ----
# load CGI bed file
CGI.bed <- read.table("E:/LabWork/genome/hg38/hg38.CGI.bed")
colnames(CGI.bed) <- c("seqnames", "start", "end", "id", "v5")
CGI.bed <- GRanges(CGI.bed)
names(CGI.bed) <- CGI.bed$id
CGI.bed <- sort(CGI.bed)

# add Feature matrix
proj_Epi <- addFeatureMatrix(proj_Epi,
    features = CGI.bed,
    matrixName = "CGI_Matrix",
    binarize = FALSE,
    force = TRUE
)

## 1.3. group by clusters ----
seCGI <- getGroupSE(
    ArchRProj = proj_Epi,
    useMatrix = "CGI_Matrix",
    groupBy = "Clusters",
    divideN = TRUE,
    scaleTo = NULL
)

rownames(seCGI) <- paste(rowData(seCGI)$seqnames,
    rowData(seCGI)$start,
    rowData(seCGI)$end,
    sep = "_"
)
rownames(seCGI@assays@data$CGI_Matrix) <- rownames(seCGI)

CGImat.cluster <- seCGI@assays@data$CGI_Matrix
CGImat.cluster[1:5, 1:5]

cluster.info <- readRDS("../02.Epi_Molecular_Subtype/cluster.info.rds")
cluster.info$Epi_Group <- factor(cluster.info$Epi_Group,
    levels = c("Normal", "Adenoma", "Group_1", "Group_2")
)
rm(temp, seCGI)

CGI.bed$peak_id <- paste(seqnames(CGI.bed),
    start(CGI.bed),
    end(CGI.bed),
    sep = "_"
)
table(CGI.bed$peak_id %in% rownames(CGImat.cluster))
CGImat.cluster <- CGImat.cluster[CGI.bed$peak_id, ]
rownames(CGImat.cluster) <- CGI.bed$id

# 2. identify CIMP status ----
mycolor <- list(
    "Gender" = c("Female" = "#cab2d6", "Male" = "#6a3d9a"),
    "MSI" = c("MSS" = "#ffff99", "MSI-H" = "#b15928", "NA" = "gray70"),
    "Side" = c("Right" = "#b2df8a", "Left" = "#33a02c"),
    "Group" = c("Group_1" = "#62b7e6", "Group_2" = "#283891"),
    "Epi_Group" = c(
        "Normal" = "#208a42", "Adenoma" = "#d51f26",
        "Group_1" = "#62b7e6", "Group_2" = "#283891"
    ),
    "CIMP_Group" = c("CIMP_High" = "#59b891", "CIMP_Low" = "#fa774f", "CIMP_Negative" = "#7a8dbf")
)

## 2.1. cluster by marker signal ----
marker.CIMP <- c(
    "CGI_676", "CGI_10999", "CGI_3874", "CGI_18521", "CGI_21212",
    "CGI_8237", "CGI_8944", "CGI_25736"
)
names(marker.CIMP) <- c(
    "RUNX3", "CACNA1G", "IGF2", "MLH1", "NEUROG1",
    "CRABP1", "SOCS1", "CDKN2A"
)

marker.sig <- CGImat.cluster[marker.CIMP, rownames(cluster.info)]
marker.sig <- apply(marker.sig, 1, function(x) {
    (x - mean(x)) / sd(x)
}) %>%
    t() %>%
    as.data.frame()

rownames(marker.sig) <- paste(names(marker.CIMP),
    marker.CIMP,
    sep = "_"
)

anno.col <- cluster.info %>%
    select(c(CIMP_Group, Side_Major, MSI_Status_Major, Gender_Major))
colnames(anno.col) <- c("CIMP_Group", "Side", "MSI", "Gender")

pdf("Heatmap_CIMP.marker.sig.pdf", 9, 3.5)
pheatmap(marker.sig,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_method = "ward.D2",
    cutree_cols = 3,
    annotation_col = anno.col,
    annotation_colors = mycolor,
    color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
    breaks = seq(-3, 3, 0.06)
)
dev.off()

pdf("Heatmap_CIMP.marker.sig.binary.pdf", 9, 3.5)
p <- pheatmap(marker.sig,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_method = "ward.D2",
    cutree_cols = 3,
    annotation_col = anno.col,
    annotation_colors = mycolor,
    color = c("#4575B4", "white", "#D73027"),
    breaks = c(-3, -0.2, 0.2, 3)
)
dev.off()

## 2.2. assign GIMP group ----
identical(rownames(cluster.info), colnames(marker.sig))
cluster.info$CIMP_nMarker <- colSums(marker.sig < -0.2)

cluster.info$CIMP_Group <-
    c("CIMP_High", "CIMP_Low", "CIMP_Negative")[cutree(p$tree_col, k = 3)]

pdf("Bar.nCIMP.Marker.pdf", 5, 2.5)
ggplot(cluster.info) +
    geom_bar(aes(x = CIMP_nMarker, fill = CIMP_Group), width = 0.6) +
    scale_fill_manual(values = mycolor$CIMP_Group) +
    theme_classic() +
    xlab("NO. of less accessible CIMP markers") +
    ylab("NO. of clusters") +
    scale_x_continuous(breaks = seq(0, 8, 1))
dev.off()

sample.info <- proj_Epi@cellColData %>% as.data.frame()

temp <- cluster.info$CIMP_Group
names(temp) <- rownames(cluster.info)
temp[c("C3", "C4")] <- "Normal"
temp[c("C9", "C28")] <- "Adenoma"

sample.info$CIMP_Group <- temp[sample.info$Clusters]
table(sample.info$CIMP_Group)

proj_Epi <- addCellColData(
    ArchRProj = proj_Epi,
    data = sample.info$CIMP_Group,
    name = "CIMP_Group",
    cells = rownames(sample.info),
    force = TRUE
)

p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    name = "CIMP_Group", embedding = "UMAP",
    pal = c(mycolor$Epi_Group, mycolor$CIMP_Group),
    size = 0.2, plotAs = "points",
    labelMeans = FALSE
)
pdf("UAMP.Epi.CIMP_Group.pdf", 7, 7)
plot(p)
dev.off()

rm(temp, anno.col, p)

## 2.3. compare with normal ----
temp <- cluster.info$CIMP_Group
names(temp) <- rownames(cluster.info)
temp[c("C3", "C4")] <- "Normal"
temp[c("C9", "C28")] <- "Adenoma"
temp <- temp[colnames(CGImat.cluster)]
temp <- aggregate(t(CGImat.cluster[marker.CIMP, ]),
    by = list(temp), FUN = mean
)
rownames(temp) <- temp$Group.1
colnames(temp) <- c("Group", paste(names(marker.CIMP),
    marker.CIMP,
    sep = "_"
))
temp$Group <- NULL
temp <- t(temp) %>% as.data.frame()
temp$CGI <- names(marker.CIMP)

library(ggpubr)
p.list <- list()
for (one in c("CIMP_High", "CIMP_Low", "CIMP_Negative")) {
    p.list[[one]] <- temp %>%
        select("Normal", one, "CGI") %>%
        melt(id.vars = "CGI") %>%
        ggpaired(
            x = "variable", y = "value", color = "CGI",
            line.color = "gray75", line.size = 0.4,
            palette = "jco", width = 0,
            xlab = "CIMP status",
            ylab = "Average ATAC signal",
        ) +
        stat_compare_means(
            paired = TRUE,
            label = "p.format", label.x.npc = "center"
        ) +
        scale_x_discrete(expand = c(0.2, 0.2))
}

pdf("Dot.CIMP.marker.sig.pdf", 9, 4.5)
wrap_plots(p.list, ncol = 3, nrow = 1)
dev.off()

pdf("Dot.CIMP.marker.sig.one.pdf", 6, 6)
for (i in 1:3) {
    print(p.list[[i]])
}
dev.off()
rm(p.list, i, one, temp)



save.image("04.Epi_CIMP.RData")
