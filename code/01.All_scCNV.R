setwd("E:/LabWork/Project/CRC_NGS_ATAC/CRC_Epi_scATAC/Results/01.All_scCNV")

source("../../code/00.Requirements.R")

# 1. load data ----
proj_CRC <- loadArchRProject(project.dir.all, force = TRUE) # 16778 cells
sample.info <- proj_CRC@cellColData %>% as.data.frame()
getAvailableMatrices(proj_CRC)
mean(proj_CRC$TSSEnrichment)
colnames(proj_CRC@cellColData)
table(proj_CRC$Clusters_type)

mycolor <- list()

# 2. basic visulization of all cells----
## 2.1. cell meta data ----
table(proj_CRC$Clusters_type)
proj_CRC$Clusters_type <- gsub("_.+?$", "", proj_CRC$Clusters_type)
proj_CRC$Clusters_type <- gsub("Macrophages", "Myeloid", proj_CRC$Clusters_type)
View(proj_CRC@cellColData %>% as.data.frame())

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

temp <- proj_CRC$cell_type
temp[temp == "C11"] <- "C1"
temp[temp == "C12"] <- "C3"
temp[temp == "LN"] <- "C2"

temp[!temp %in% c("C1", "C2", "C3")] <- "None"

temp <- paste(proj_CRC$Sample_2, temp, sep = "_")
temp[grep("None", temp)] <- "None"

proj_CRC$Replication <- temp

# colors <- paletteDiscrete(proj_Epi$Replication)
# colors["None"] <- "gray75"
p <- plotEmbedding(
    ArchRProj = proj_CRC, colorBy = "cellColData",
    name = "Replication", embedding = "UMAP",
    size = 0.2, plotAs = "points",
    pal = colors,
    labelMeans = FALSE
)
pdf("UMAP.All.Replication.pdf", 7, 7)
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
proj_CRC$location_new <- gsub("LN", "C", proj_CRC$location)

plot.data <- table(proj_CRC$location_new, proj_CRC$Clusters_type) %>%
    as.data.frame()
temp <- c("AD" = "Adenoma", "C" = "Cancer", "T" = "Normal tissue")

plot.data$Var1 <- temp[plot.data$Var1]
table(plot.data$Var1)
colnames(plot.data) <- c("Location", "Clusters_type", "Count")
plot.data$Location <- factor(plot.data$Location,
    levels = c("Normal tissue", "Adenoma", "Cancer")
)
plot.data$Clusters_type <- factor(plot.data$Clusters_type,
    levels = rev(c("Epithelial", "Fibroblast", "T", "B", "Myeloid"))
)

temp <- table(proj_CRC$location_new)[names(temp)] %>% as.data.frame()
temp$Var1 <- c("Adenoma", "Cancer", "Normal tissue")

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

table(proj_CRC$Sample_2, proj_CRC$Sample)
table(proj_CRC$Sample_2)
plot.data <- table(proj_CRC$Sample_2, proj_CRC$Clusters_type) %>%
    as.data.frame()
colnames(plot.data) <- c("Sample", "Clusters_type", "Count")

plot.data$Percent <- plot.data$Count / table(proj_CRC$Sample)[plot.data$Sample] * 100 %>% as.numeric()

plot.data$MACS <- "MACS"
plot.data[grep("P09", plot.data$Sample), ]$MACS <- "no MACS"
# plot.data$Sample <- gsub("-nofacs", "", plot.data$Sample)

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

## 2.4. big wig for IGV visualization ----
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
# 3. scCNV analysis of all cells ----
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
saveRDS(sw, "sw.rds")
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

cluster.rename <- read.table("cluster_rename.txt", header = 1)
temp <- cluster.rename$manual
names(temp) <- cluster.rename$Previous
proj_Epi$Clusters_2 <- temp[proj_Epi$Clusters]
mycolor <- ArchR::paletteDiscrete(proj_Epi$Clusters)
names(mycolor) <- temp[names(mycolor)]

p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    name = "Clusters_2", embedding = "UMAP",
    pal = mycolor,
    size = 0.2, plotAs = "points"
)

pdf("UMAP.Epi.Clusters_2.pdf", 7, 7)
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

proj_Epi$Sample_2 <- gsub("-nofacs", "", proj_Epi$Sample)
proj_Epi$Sample_2 <- rename.patient[proj_Epi$Sample_2]

temp <- proj_Epi$cell_type
table(proj_Epi$Sample_2, proj_Epi$cell_type)
temp[temp == "C11"] <- "C1"
temp[temp == "C12"] <- "C3"
temp[temp == "LN"] <- "C2"

temp[!temp %in% c("C1", "C2", "C3")] <- "None"

temp <- paste(proj_Epi$Sample_2, temp, sep = "_")
temp[grep("None", temp)] <- "None"

proj_Epi$Replication <- temp

colors <- paletteDiscrete(proj_Epi$Replication)
colors["None"] <- "gray75"
p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    name = "Replication", embedding = "UMAP",
    size = 0.2, plotAs = "points",
    pal = colors,
    labelMeans = FALSE
)
pdf("UMAP.Epi.Replication.pdf", 7, 7)
plot(p)
dev.off()

## 4.2. batch effect -----
AMI.cluster <- c(
    Normal = aricode::AMI(
        sample.info.epi %>% filter(Epi_type == "Normal") %>% pull(Clusters),
        sample.info.epi %>% filter(Epi_type == "Normal") %>% pull(Sample)
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

table(proj_Epi$Sample)
proj_Epi@reducedDims$Harmony$params[3]
proj_Epi@reducedDims$IterativeLSI_merge$useMatrix
proj_Epi@embeddings$UMAP$params

pdf("Harmony.converge.pdf", 5, 4)
proj_Epi <- addHarmony(
    ArchRProj = proj_Epi,
    reducedDims = "IterativeLSI_merge",
    name = "Harmony",
    groupBy = "Sample",
    plot_convergence = TRUE,
    force = TRUE
)
dev.off()

proj_Epi <- addUMAP(
    ArchRProj = proj_Epi, reducedDims = "Harmony",
    name = "UMAP_Harmony",
    nNeighbors = 30, minDist = 0.4,
    metric = "cosine",
    force = TRUE
)

proj_Epi <- addClusters(
    input = proj_Epi,
    reducedDims = "Harmony",
    dimsToUse = 1:25,
    knnAssign = 6,
    maxClusters = 50,
    method = "Seurat",
    name = "Clusters_Harmony",
    force = TRUE
)
sample.info.epi <- proj_Epi@cellColData %>% as.data.frame()

p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    name = "Clusters_Harmony", embedding = "UMAP_Harmony",
    # pal = mycolor,
    size = 0.2, plotAs = "points"
)

pdf("UMAP.Harmony.Epi.Clusters_Harmony.pdf", 7, 7)
plot(p)
dev.off()

p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    name = "Clusters_2", embedding = "UMAP_Harmony",
    pal = mycolor,
    size = 0.2, plotAs = "points"
)

pdf("UMAP.Harmony.Epi.Clusters_2.pdf", 7, 7)
plot(p)
dev.off()

p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    name = "Clusters_type", embedding = "UMAP_Harmony",
    size = 0.2, plotAs = "points"
)
pdf("UMAP.Harmony.Epi.Clusters_type.pdf", 7, 7)
plot(p)
dev.off()

p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    name = "Sample_2", embedding = "UMAP_Harmony",
    size = 0.2, plotAs = "points",
    labelMeans = FALSE
)
pdf("UMAP.Harmony.Epi.Sample.pdf", 7, 7)
plot(p)
dev.off()

p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    name = "new_location", embedding = "UMAP_Harmony",
    size = 0.2, plotAs = "points"
)

pdf("UMAP.Harmony.Epi.new_location.pdf", 7, 7)
plot(p)
dev.off()

p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    name = "Type_location", embedding = "UMAP_Harmony",
    size = 0.2, plotAs = "points"
)

pdf("UMAP.Harmony.Epi.Type_location.pdf", 7, 7)
plot(p)
dev.off()

AMI.cluster.harmony <- c(
    Normal = aricode::AMI(
        sample.info.epi %>% filter(Epi_type == "Normal") %>% pull(Clusters_Harmony),
        sample.info.epi %>% filter(Epi_type == "Normal") %>% pull(Sample)
    ),
    Adenoma = aricode::AMI(
        sample.info.epi %>% filter(Epi_type == "Adenoma") %>% pull(Clusters_Harmony),
        sample.info.epi %>%
            filter(Epi_type == "Adenoma") %>%
            pull(Sample)
    ),
    Malignant = aricode::AMI(
        sample.info.epi %>% filter(Epi_type == "Malignant") %>% pull(Clusters_Harmony),
        sample.info.epi %>% filter(Epi_type == "Malignant") %>% pull(Sample)
    )
)
AMI.cluster
AMI.cluster.harmony

## 4.3. nClusters and resolution ----
resolutions <- seq(0.1, 2, 0.05)
nClusters <- data.frame(
    "resolutions" = resolutions,
    "Total" = 0,
    "Normal" = 0,
    "Adenoma" = 0,
    "Malignant" = 0
)

one <- 1
for (one in seq_len(nrow(nClusters))) {
    proj_Epi <- addClusters(
        input = proj_Epi,
        reducedDims = "IterativeLSI_merge",
        dimsToUse = 1:25,
        knnAssign = 6,
        maxClusters = 50,
        method = "Seurat",
        name = "Clusters_test",
        resolution = nClusters$resolutions[one],
        force = TRUE
    )
    temp <- table(proj_Epi$Clusters_test, proj_Epi$Epi_type)
    temp <- colnames(temp)[apply(temp, 1, which.max)]

    nClusters[one, "Total"] <- length(temp)
    nClusters[one, "Normal"] <- sum(temp == "Normal")
    nClusters[one, "Adenoma"] <- sum(temp == "Adenoma")
    nClusters[one, "Malignant"] <- sum(temp == "Malignant")
}
nClusters <- nClusters %>% mutate(Total = Normal + Adenoma + Malignant)
saveRDS(nClusters, "nClusters.rds")

pdf("Dot.resolution.nClusters.pdf", 5, 3)
ggplot(nClusters %>% filter(resolutions <= 1.15), aes(x = resolutions)) +
    geom_line(aes(y = Total, color = "Total")) +
    geom_line(aes(y = Normal, color = "Normal")) +
    geom_line(aes(y = Adenoma, color = "Adenoma")) +
    geom_line(aes(y = Malignant, color = "Malignant")) +
    geom_point(aes(y = Total, color = "Total")) +
    geom_point(aes(y = Normal, color = "Normal")) +
    geom_point(aes(y = Adenoma, color = "Adenoma")) +
    geom_point(aes(y = Malignant, color = "Malignant")) +
    scale_color_manual(values = c(
        "Total" = "#5e2953", "Normal" = "#208a42",
        "Adenoma" = "#d51f26", "Malignant" = "#272d6a"
    )) +
    geom_hline(yintercept = c(25, 2), linetype = "dashed") +
    geom_vline(xintercept = 0.5, linetype = "dashed") +
    ylab("Number of clusters") +
    scale_y_continuous(breaks = c(0, 2, 10, 20, 25, 30), labels = c(0, 2, 10, 20, 25, 30)) +
    scale_x_continuous(breaks = seq(0, 1.1, 0.1), labels = seq(0, 1.1, 0.1)) +
    theme_classic()
dev.off()

## 4.5. barplot composition ----

table(proj_Epi$location)
proj_Epi$location_new <- gsub("LN", "C", proj_Epi$location)

proj_Epi$Sample_2 <- gsub("-nofacs", "", proj_Epi$Sample)
proj_Epi$Sample_2 <- rename.patient[proj_Epi$Sample_2]

plot.data <- table(proj_Epi$location_new, proj_Epi$Epi_type, proj_Epi$Sample_2) %>%
    as.data.frame() %>%
    filter(Var1 == "C") %>%
    mutate(Var1 = "Cancer")

colnames(plot.data) <- c("Location", "Epi_type", "Patient", "Count")

plot.data$Epi_type <- factor(plot.data$Epi_type,
    levels = c("Normal", "Adenoma", "Malignant")
)
plot.data <- plot.data %>%
    filter(Patient != "P20")

pdf("Bar.Cancer.Cell_Type.pdf", 7, 3)
ggplot() +
    geom_bar(
        data = plot.data, aes(x = Patient, y = Count, fill = Epi_type),
        stat = "identity", width = 0.7, position = "fill"
    ) +
    scale_fill_manual(values = c(
        "Normal" = "#208a42", "Adenoma" = "#d51f26", "Malignant" = "#272d6a"
    )) +
    facet_wrap(~Location) +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ylab("Fraction of cells")
dev.off()

# 5. scCNV for Epi clusters ----
load("CRC_CNV.rda")
proj_Epi$Epi_type <- "Malignant"
proj_Epi$Epi_type[proj_Epi$Clusters %in% c("C3", "C4")] <- "Normal"
proj_Epi$Epi_type[proj_Epi$Clusters %in% c("C28", "C9")] <- "Adenoma"

table(proj_Epi$Epi_type)

anno.row <- sample.info.epi %>%
    select(Sample, Clusters, Epi_type)

anno.color <- list(
    Clusters = paletteDiscrete(proj_Epi$Clusters),
    Epi_type = paletteDiscrete(proj_Epi$Epi_type),
    Sample = paletteDiscrete(proj_Epi$Sample),
    seqnames = rep(c("#969696", "#212121"), 11)
)
names(anno.color$seqnames) <- paste("chr", 1:22, sep = "")
anno.color$Clusters

dir.create("Epi_By_Cluster")
for (one in unique(proj_Epi$Clusters)) {
    print(one)
    cell.select <- rownames(sample.info.epi)[sample.info.epi$Clusters == one]
    png(paste0("Epi_By_Cluster/Heatmap.CNV.FC.", one, ".png"), 800, 600)
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
sample.info.epi$Epi_type <- factor(sample.info.epi$Epi_type,
    levels = c("Normal", "Adenoma", "Malignant")
)

for (one in unique(proj_Epi$Sample)) {
    print(one)
    cell.select <- rownames(sample.info.epi)[sample.info.epi$Sample == one]
    png(paste0("Epi_By_Sample/Heatmap.CNV.FC.png"), 800, 600)
    p <- pheatmap::pheatmap(CNV.FC[cell.select, ],
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
    cell.select <- cell.select[p$tree_row$order]
    cell.select <- cell.select[order(sample.info.epi[cell.select, ]$Epi_type)]
    while (dev.cur() != 1) {
        dev.off()
    }
    png(paste0("Epi_By_Sample/Heatmap.CNV.FC.", one, ".png"), 800, 600)
    pheatmap::pheatmap(CNV.FC[cell.select, ],
        color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
        cluster_rows = FALSE, cluster_cols = FALSE,
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

pdf("CNV.legend.pdf", 3, 3)
ggplot(data = data.frame()) +
    geom_point(aes(x = 1:100, y = 1:100, color = 1:100)) +
    scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 50)
dev.off()

plot(colorRampPalette(rev(brewer.pal(9, "RdBu")))(100))

# 6. marker peak for each type ----
## 6.1 identify marker peaks ----
table(proj_Epi$Clusters)
table(proj_Epi$Epi_type)

markersPeaks <- getMarkerFeatures(
    ArchRProj = proj_Epi,
    useMatrix = "PeakMatrix",
    groupBy = "Epi_type",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.001 & Log2FC >= 1")
marker.list <- sapply(markerList, function(x) {
    x <- x %>%
        as.data.frame() %>%
        mutate(name = paste(seqnames, start, end, sep = "_")) %>%
        pull(name)
    return(x)
})
sapply(marker.list, length)

temp <- unlist(marker.list)  %>% table()
dup <- names(temp)[temp > 1]
marker.list <- sapply(marker.list, function(x) {
    x <- x[!(x %in% dup)]
    return(x)
})

## 6.2 heatmap of each clusters ----
proj_Epi$Sample_Type <- paste(proj_Epi$Sample, proj_Epi$Epi_type, sep = "_")

sePeaks <- getGroupSE(
    ArchRProj = proj_Epi,
    useMatrix = "PeakMatrix",
    groupBy = "Sample_Type",
    divideN = TRUE,
    scaleTo = NULL
)
saveRDS(sePeaks, "sePeaks.Sample_Type.rds")
sePeaks <- readRDS("sePeaks.Sample_Type.rds")

selected <- table(proj_Epi$Sample_Type)
selected <- names(selected)[selected > 10]
selected <- setdiff(selected, "COAD09_Adenoma")

# sePeaks <- readRDS("../03.Epi_Molecular_Subtype/sePeaks.cluster.rds")
PeakMatrix <- sePeaks@assays@data$PeakMatrix

rownames(PeakMatrix) <- sePeaks@elementMetadata %>%
    as.data.frame() %>%
    mutate(name = paste(seqnames, start, end, sep = "_")) %>%
    pull(name)

table(c(marker.list$Normal, marker.list$Adenoma, marker.list$Malignant) %in% rownames(PeakMatrix))

plot.data <- PeakMatrix[
    c(marker.list$Normal, marker.list$Adenoma, marker.list$Malignant),
    selected
]

plot.data <- apply(plot.data, 1, function(x) {
    x <- (x - mean(x)) / sd(x)
    return(x)
}) %>% t()

plot.data[plot.data > 1.5] <- 1.5
plot.data[plot.data < (-1.5)] <- -1.5

names(marker.list)
anno.row <- data.frame(
    row.names = rownames(plot.data),
    Group = factor(rep(c("Normal", "Adenoma", "Malignant"),
        times = c(length(marker.list$Normal), length(marker.list$Adenoma), length(marker.list$Malignant))
    ), levels = c("Normal", "Adenoma", "Malignant")),
    Peak = "Up"
)
# anno.row <- anno.row[order(anno.row$Group, rowMeans(plot.data[, 18:27]) - rowMeans(plot.data[, 3:17])), ]

ann.col <- data.frame(
    row.names = colnames(plot.data),
    Epi_Type = gsub("^.+?_", "", colnames(plot.data)),
    temp = rep("temp", length(colnames(plot.data)))
)
ann.col$Epi_Type <- factor(ann.col$Epi_Type, levels = c("Normal", "Adenoma", "Malignant"))

ann.col <- ann.col[order(ann.col$Epi_Type, rnorm(ncol(plot.data))), ]
# ann.col <- ann.col[order(ann.col$Epi_Type, colSums(plot.data)), ]
plot.data <- plot.data[row.names(anno.row), rownames(ann.col)]

png("Heatmap.marker.Epi_Type.patient.Up.png", 800, 800)
pheatmap(plot.data,
    # color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")[2:8]))(100),
    scale = "none",
    annotation_row = anno.row %>% select(Group),
    annotation_col = ann.col %>% select(Epi_Type),
    annotation_colors = list(
        Group = c(Normal = "#208a42", Adenoma = "#d51f26", Malignant = "#272d6a"),
        Epi_Type = c(Normal = "#208a42", Adenoma = "#d51f26", Malignant = "#272d6a")
    ),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    gaps_row = c(length(marker.list$Normal), length(marker.list$Normal) + length(marker.list$Adenoma)),
    gaps_col = c(28, 35)
)
dev.off()

# 7. TF deviation heatmap ----
plot.data <- data.table::fread("T-AD-C-specific_heatmap_data.csv") %>%
    as.data.frame()

rownames(plot.data) <- plot.data$V1
plot.data$V1 <- NULL
plot.data[1:5, 1:5]

table(colnames(plot.data) %in% rownames(proj_Epi@cellColData))

sample.info.epi <- proj_Epi@cellColData %>%
    as.data.frame() %>%
    .[colnames(plot.data), ]

anno.col <- sample.info.epi %>%
    select(Sample, Epi_type)

cutoff <- 2
plot.data[plot.data > cutoff] <- cutoff
plot.data[plot.data < -cutoff] <- -cutoff

# scales::show_col(ArchRPalettes$solarExtra)
# table(anno.col$Epi_type)
dev.off()

pdf("Heatmap.TF.Epi.type.pdf", 12, 7)
pheatmap::pheatmap(plot.data,
    scale = "none",
    color = colorRampPalette(ArchRPalettes$solarExtra[3:7])(100),
    cluster_rows = FALSE, cluster_cols = FALSE,
    annotation_col = anno.col,
    annotation_colors = list(
        Epi_type = paletteDiscrete(proj_Epi$Epi_type),
        Sample = paletteDiscrete(proj_Epi$Sample)
    ),
    show_rownames = FALSE, show_colnames = FALSE,
    gaps_col = c(1910, 1910 + 800),
    gaps_row = c(45, 53)
)
dev.off()

save(CNV.FC, sample.info, anno.col, anno.row, anno.color,
    file = "CRC_CNV.rda"
)

rm(CNV.FC, sample.info, anno.col, anno.row, anno.color)
rm(proj_CRC, proj_Epi)
