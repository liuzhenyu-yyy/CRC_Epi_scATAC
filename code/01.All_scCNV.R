setwd("E:/LabWork/Project/CRC_NGS_ATAC/CRC_Epi_scATAC/Results/01.All_scCNV")

source("../../code/00.Requirements.R")

# 1. load data ----
proj_CRC <- loadArchRProject(project.dir.all, force = TRUE)
sample.info <- proj_CRC@cellColData %>% as.data.frame()
getAvailableMatrices(proj_CRC)

colnames(proj_CRC@cellColData)
table(proj_CRC$Clusters_type)

p <- plotEmbedding(
    ArchRProj = proj_CRC, colorBy = "cellColData",
    name = "Sample", embedding = "UMAP",
    size = 0.2, plotAs = "points"
)
pdf("UAMP.All.Sample.pdf", 7, 7)
plot(p)
dev.off()

p <- plotEmbedding(
    ArchRProj = proj_CRC, colorBy = "cellColData",
    name = "Clusters_type", embedding = "UMAP",
    size = 0.2, plotAs = "points"
)
pdf("UAMP.All.Clusters_type.pdf", 7, 7)
plot(p)
dev.off()

p <- plotEmbedding(
    ArchRProj = proj_CRC, colorBy = "cellColData",
    name = "Clusters", embedding = "UMAP",
    size = 0.2, plotAs = "points"
)
pdf("UAMP.All.Clusters.pdf", 7, 7)
plot(p)
dev.off()

# 2.Call scCNV ----

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

# 3. Plot scCNV for all cells ----
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

## 3.1. CNV by patients ----
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

## 3.2. CNV by Cell Type ----
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

# 4. load epi project ----
proj_Epi <- loadArchRProject(project.dir.epi, force = TRUE)
table(proj_CRC$Clusters_type)
table(proj_Epi$Clusters_type)

colnames(proj_CRC@cellColData)
colnames(proj_Epi@cellColData)
table(proj_Epi$Clusters)

p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    name = "Sample", embedding = "UMAP",
    size = 0.2, plotAs = "points"
)
pdf("UAMP.Epi.Sample.pdf", 7, 7)
plot(p)
dev.off()

p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    name = "Clusters_type", embedding = "UMAP",
    size = 0.2, plotAs = "points"
)
pdf("UAMP.Epi.Clusters_type.pdf", 7, 7)
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

p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    name = "new_location", embedding = "UMAP",
    size = 0.2, plotAs = "points"
)

pdf("UAMP.Epi.new_location.pdf", 7, 7)
plot(p)
dev.off()

p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    name = "Type_location", embedding = "UMAP",
    size = 0.2, plotAs = "points"
)

pdf("UAMP.Epi.Type_location.pdf", 7, 7)
plot(p)
dev.off()

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
