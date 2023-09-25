setwd("E:/LabWork/Project/CRC_NGS_ATAC/CRC_Epi_scATAC/Results/01.All_scCNV")

source("../../code/00.Requirements.R")

# 1. load data ----
proj_CRC <- loadArchRProject(project.dir.all, force = TRUE)
sample.info <- proj_CRC@cellColData %>% as.data.frame()
getAvailableMatrices(proj_CRC)

colnames(proj_CRC@cellColData)
table(proj_CRC$Clusters_type)

mycolor <- list()

# 2. basic visulization of all cells----
## 2.1. cell meta data ----
table(proj_CRC$Clusters_type)
proj_CRC$Clusters_type <- gsub("_.+?$", "", proj_CRC$Clusters_type)
proj_CRC$Clusters_type <- gsub("Macrophages", "Myeloid", proj_CRC$Clusters_type)

mycolor$Clusters_type <- paletteDiscrete(proj_CRC$Clusters_type, set = "summerNight")
names(ArchRPalettes)
scales::show_col(ArchRPalettes$blueYellow)

proj_CRC$Sample_2 <- gsub("-nofacs", "", proj_CRC$Sample)
temp <- table(proj_CRC$Sample_2) %>% names()
rename.patient <- paste0("P", formatC(seq_along(temp), width = 2, flag = 0))
names(rename.patient) <- temp
proj_CRC$Sample_2 <- rename.patient[proj_CRC$Sample_2]
p <- plotEmbedding(
    ArchRProj = proj_CRC, colorBy = "cellColData",
    name = "Sample_2", embedding = "UMAP",
    size = 0.2, plotAs = "points",
    labelMeans = FALSE
)
pdf("UMAP.All.Sample1.pdf", 7, 7)
plot(p)
dev.off()

p <- plotEmbedding(
    ArchRProj = proj_CRC, colorBy = "cellColData",
    name = "Clusters_type", embedding = "UMAP",
    size = 0.2, plotAs = "points",
    pal = mycolor$Clusters_type,
    labelMeans = FALSE
)
pdf("UMAP.All.Clusters_type.pdf", 6, 6)
plot(p)
dev.off()

p <- plotEmbedding(
    ArchRProj = proj_CRC, colorBy = "cellColData",
    name = "Clusters", embedding = "UMAP",
    size = 0.2, plotAs = "points"
)
pdf("UMAP.All.Clusters.pdf", 7, 7)
plot(p)
dev.off()

## 2.2. marker genes ----
markers <- c("EPCAM", "KRT19", "COL2A1", "THY1", "PTPRC", "CD3D", "CD79A", "CD14")

proj_CRC <- addImputeWeights(proj_CRC, reducedDims = "IterativeLSI_merge")
p <- plotEmbedding(
    ArchRProj = proj_CRC, colorBy = "GeneScoreMatrix",
    name = markers, embedding = "UMAP",
    pal = ArchRPalettes$blueYellow,
    imputeWeights = getImputeWeights(proj_CRC),
    size = 0.2, plotAs = "points", rastr = TRUE
)

pdf("UMAP.All.markers.pdf", 16, 10)
patchwork::wrap_plots(plotlist = p, nrow = 2, ncol = 4, byrow = TRUE)
dev.off()

## 2.3. barplot composition----
mycolor$Location <- c(
    "Normal tissue" = "#208a42", "Adenoma" = "#d51f26",
    "Cancer" = "#272d6a", "Lymph" = "#89288f"
)
scales::show_col(ArchRPalettes[[1]])
table(proj_CRC$location)

plot.data <- table(proj_CRC$location, proj_CRC$Clusters_type) %>%
    as.data.frame()
temp <- c("AD" = "Adenoma", "C" = "Cancer", "LN" = "Lymph", "T" = "Normal tissue")

plot.data$Var1 <- temp[plot.data$Var1]
table(plot.data$Var1)
colnames(plot.data) <- c("Location", "Clusters_type", "Count")
plot.data$Location <- factor(plot.data$Location,
    levels = c("Normal tissue", "Adenoma", "Cancer", "Lymph")
)
plot.data$Clusters_type <- factor(plot.data$Clusters_type,
    levels = rev(c("Epithelial", "Fibroblast", "T", "B", "Myeloid"))
)

temp <- table(proj_CRC$location)[names(temp)] %>% as.data.frame()
temp$Var1 <- c("Adenoma", "Cancer", "Lymph", "Normal tissue")

pdf("Bar.Location.Cell_Type.pdf", 3.6, 3)
ggplot() +
    geom_bar(data = plot.data, aes(x = Location, y = Count, fill = Clusters_type),
        stat = "identity", width = 0.7) +
    scale_fill_manual(values = mycolor$Clusters_type) +
    geom_text(data = temp, aes(x = Var1, y = Freq, label = Freq), vjust = -0.3, size = 2) +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ylab("No. of cells")
dev.off()

table(proj_CRC$Sample)
plot.data <- table(proj_CRC$Sample, proj_CRC$Clusters_type) %>%
    as.data.frame()
colnames(plot.data) <- c("Sample", "Clusters_type", "Count")

plot.data$Percent <- plot.data$Count / table(proj_CRC$Sample)[plot.data$Sample] * 100 %>% as.numeric()

plot.data$MACS <- "MACS"
plot.data[grep("COAD13", plot.data$Sample), ]$MACS <- "no MACS"
plot.data$Sample <- gsub("-nofacs", "", plot.data$Sample)

# orders <- plot.data %>%
#     filter(Clusters_type == "Epithelial") %>%
#     arrange(desc(Percent)) %>%
#     pull(Sample)
# plot.data$Sample <- factor(plot.data$Sample, levels = orders)
plot.data$Clusters_type <- factor(plot.data$Clusters_type,
    levels = c("Epithelial", "Fibroblast", "T", "B", "Myeloid")
)

pdf("Bar.Sample.Cell_Type.pdf", 7, 3)
ggplot() +
    geom_bar(
        data = plot.data, aes(x = Sample, y = Count, fill = Clusters_type),
        stat = "identity", show.legend = FALSE
    ) +
    scale_fill_manual(values = mycolor$Clusters_type) +
        theme_classic() +
        facet_grid(cols = vars(MACS), scales = "free_x", space = "free_x") +
        xlab("Percent of cells") +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1)
        )
dev.off()

## 2.4 big wig for IGV visualization ----
table(proj_CRC$Clusters_type)
proj_CRC <- addGroupCoverages(
    ArchRProj = proj_CRC,
    groupBy = "Clusters_type"
)
getGroupBW(
    ArchRProj = proj_CRC,
    groupBy = "Clusters_type",
    normMethod = "ReadsInTSS",
    tileSize = 30,
    maxCells = 2000,
    threads = 4
)

proj_CRC <- saveArchRProject(ArchRProj = proj_CRC, load = TRUE)
# 3. scCNV analysis ----
## 3.1. call scCNV ----

# dir.create("cells")
# for (one in unique(proj_CRC$Sample)) {
#     keep <- proj_CRC@cellColData %>%
#         as.data.frame() %>%
#         filter(Sample == one) %>%
#         rownames()
#     saveRDS(keep, paste("data/cells/", one, ".rds", sep = ""))
# }
#
## Call scCNV from each arrow file with script: CNV_from_Arrow.R
## command: Rscript CNV_from_Arrow.R Arrows/${sample}.arrow cells/${sample}.rds ${sample}
## results save in data/rds
#
# obj.list <- dir("data/rds/", full.names = TRUE)
# length(obj.list)

## merge cna object of each sample (each arrow file)
# cnaobj <- readRDS(obj.list[1])
# for (one in 2:length(obj.list)) {
#   temp <- readRDS(obj.list[one])
#   cnaobj <- cbind(cnaobj, temp)
#   rm(temp)
#   print(ncol(cnaobj))
# }
# saveRDS(cnaobj, "data/cnaobj.rds")

cnaobj <- readRDS("data/cnaobj.rds")

# fold change
names(cnaobj@assays@data)
cnaobj@assays@data$CNA[1:5, 1:5]

CNV.FC <- as.data.frame(cnaobj@assays@data$log2FC)
sw <- cnaobj@rowRanges # 1190
sw.filtered <- sw[sw$percentEffectiveLength > 90] # 1067
length(sw.filtered)

colnames(CNV.FC) <- cnaobj@colData@rownames
rownames(CNV.FC) <- cnaobj@rowRanges$name
dim(CNV.FC)

CNV.FC <- as.data.frame(t(CNV.FC))
CNV.FC <- CNV.FC[rownames(sample.info), sw.filtered$name]

cutoff <- 1.5
CNV.FC[CNV.FC > cutoff] <- cutoff
CNV.FC[CNV.FC < (-cutoff)] <- (-cutoff)

table(rownames(sample.info) %in% rownames(CNV.FC))

## 3.2 Plot scCNV for all cells ----
samples <- sort(unique(proj_CRC$Sample))
color.patient <- paletteDiscrete(samples, set = "stallion", reverse = FALSE)
color.cell <- paletteDiscrete(proj_CRC$Clusters_type, set = "stallion", reverse = FALSE)

anno.row <- sample.info %>%
    select(Sample, Clusters_type)
anno.col <- data.frame(
    row.names = sw.filtered$name,
    "seqnames" = seqnames(sw.filtered)
)

anno.color <- list(
    Clusters_type = color.cell,
    Sample = color.patient,
    seqnames = rep(c("#969696", "#212121"), 11)
)
names(anno.color$seqnames) <- paste("chr", 1:22, sep = "")

#CNV by patients
dir.create("By_Sample")
for (one in unique(proj_CRC$Sample)) {
    print(one)
    cell.select <- rownames(sample.info)[sample.info$Sample == one]
    png(paste0("By_Sample/Heatmap.CNV.FC.", one, ".png"), 1000, 750)
    pheatmap::pheatmap(CNV.FC[cell.select, ],
        color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
        cluster_rows = TRUE, cluster_cols = FALSE,
        clustering_method = "ward.D2",
        annotation_row = anno.row[cell.select, ],
        annotation_col = anno.col,
        annotation_colors = anno.color,
        annotation_legend = FALSE,
        show_rownames = FALSE, show_colnames = FALSE,
        treeheight_row = 0,
        main = one
    )
    dev.off()
}
gc()

# CNV by Cell Type
dir.create("By_CellType")
for (one in unique(proj_CRC$Clusters_type)) {
    print(one)
    cell.select <- rownames(sample.info)[sample.info$Clusters_type == one]
    png(paste0("By_CellType/Heatmap.CNV.FC.", one, ".png"), 1000, 750)
    pheatmap::pheatmap(CNV.FC[cell.select, ],
        color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
        cluster_rows = TRUE, cluster_cols = FALSE,
        clustering_method = "ward.D2",
        annotation_row = anno.row[cell.select, ],
        annotation_col = anno.col,
        annotation_colors = anno.color,
        annotation_legend = FALSE,
        show_rownames = FALSE, show_colnames = FALSE,
        treeheight_row = 0,
        main = one
    )
    dev.off()
}
gc()

# 4. Epithelium cells subset ----
## 4.1. load epi project ----
proj_Epi <- loadArchRProject(project.dir.epi, force = TRUE)
table(proj_CRC$Clusters_type)
table(proj_Epi$Clusters_type)

colnames(proj_CRC@cellColData)
colnames(proj_Epi@cellColData)
table(proj_Epi$Clusters)

proj_Epi$Sample_2 <- gsub("-nofacs", "", proj_Epi$Sample)

p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    name = "Sample_2", embedding = "UMAP",
    size = 0.2, plotAs = "points",
    labelMeans = FALSE
)
pdf("UMAP.Epi.Sample.pdf", 7, 7)
plot(p)
dev.off()

p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    name = "Clusters_type", embedding = "UMAP",
    size = 0.2, plotAs = "points"
)
pdf("UMAP.Epi.Clusters_type.pdf", 7, 7)
plot(p)
dev.off()

p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    name = "Clusters", embedding = "UMAP",
    size = 0.2, plotAs = "points"
)

pdf("UMAP.Epi.Clusters.pdf", 7, 7)
plot(p)
dev.off()

p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    name = "new_location", embedding = "UMAP",
    size = 0.2, plotAs = "points"
)

pdf("UMAP.Epi.new_location.pdf", 7, 7)
plot(p)
dev.off()

p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    name = "Type_location", embedding = "UMAP",
    size = 0.2, plotAs = "points"
)

pdf("UMAP.Epi.Type_location.pdf", 7, 7)
plot(p)
dev.off()

## 4.2. correlation of epithlium clusters ---
mycolor$Cell_Type <- c("Normal" = "#208a42", "Adenoma" = "#d51f26", "Malignant" = "#272d6a")
sePeaks <- getGroupSE(
    ArchRProj = proj_Epi,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters",
    divideN = TRUE,
    scaleTo = NULL
)
ann.row <- data.frame(
    row.names = colnames(sePeaks@assays@data$PeakMatrix),
    "Cell_Type" = rep("Malignant", ncol(sePeaks@assays@data$PeakMatrix)),
    "test" = 0
)
ann.row[c("C3", "C4"), ]$Cell_Type <- "Normal"
ann.row[c("C28", "C9"), ]$Cell_Type <- "Adenoma"
ann.row$Cell_Type <- factor(ann.row$Cell_Type,
    levels = c("Normal", "Adenoma", "Malignant")
)
ann.row <- ann.row[order(ann.row$Cell_Type), ]
ann.row$test <- NULL

dist.clusters <- dist(t(sePeaks@assays@data$PeakMatrix))
dist.clusters <- as.matrix(dist.clusters)
dist.clusters <- dist.clusters[rownames(ann.row), rownames(ann.row)]
pdf("Heatmap.cluster.distace.pdf", 5, 4)
pheatmap(dist.clusters,
    color = colorRampPalette(c("red", "White", "blue"))(100),
    annotation_row = ann.row,
    annotation_col = ann.row,
    cluster_rows = FALSE, cluster_cols = FALSE,
    method = "ward.D2",
    annotation_color = mycolor
)

# 5. scCNV for Epi clusters ----
load("CRC_CNV.rda")
sample.info.epi <- proj_Epi@cellColData %>% as.data.frame()

anno.row <- sample.info.epi %>%
    select(Sample, Clusters, new_location)

anno.color <- list(
    Clusters = paletteDiscrete(proj_Epi$Clusters),
    new_location = paletteDiscrete(proj_Epi$new_location),
    Sample = paletteDiscrete(proj_Epi$Sample),
    seqnames = rep(c("#969696", "#212121"), 11)
)
names(anno.color$seqnames) <- paste("chr", 1:22, sep = "")
anno.color$Clusters

dir.create("Epi_By_Cluster")
one <- "C1"
for (one in unique(proj_Epi$Clusters)) {
    print(one)
    cell.select <- rownames(sample.info.epi)[sample.info.epi$Clusters == one]
    png(paste0("Epi_By_Cluster/Heatmap.CNV.FC.", one, ".png"), 1000, 750)
    pheatmap::pheatmap(CNV.FC[cell.select, ],
        color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
        cluster_rows = TRUE, cluster_cols = FALSE,
        clustering_method = "ward.D2",
        annotation_row = anno.row[cell.select, ],
        annotation_col = anno.col,
        annotation_colors = anno.color,
        annotation_legend = FALSE,
        show_rownames = FALSE, show_colnames = FALSE,
        treeheight_row = 0,
        main = one
    )
    dev.off()
}
gc()

dir.create("Epi_By_Sample")
for (one in unique(proj_Epi$Sample)) {
    print(one)
    cell.select <- rownames(sample.info.epi)[sample.info.epi$Sample == one]
    png(paste0("Epi_By_Sample/Heatmap.CNV.FC.", one, ".png"), 1000, 750)
    pheatmap::pheatmap(CNV.FC[cell.select, ],
        color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
        cluster_rows = TRUE, cluster_cols = FALSE,
        clustering_method = "ward.D2",
        annotation_row = anno.row[cell.select, ],
        annotation_col = anno.col,
        annotation_colors = anno.color,
        annotation_legend = FALSE,
        show_rownames = FALSE, show_colnames = FALSE,
        treeheight_row = 0,
        main = one
    )
    dev.off()
}
gc()

save(CNV.FC, sample.info, anno.col, anno.row, anno.color,
    file = "CRC_CNV.rda"
)

rm(CNV.FC, sample.info, anno.col, anno.row, anno.color)
rm(proj_CRC, proj_Epi)
