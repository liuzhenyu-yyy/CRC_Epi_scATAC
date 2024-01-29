setwd("E:/LabWork/Project/CRC_NGS_ATAC/CRC_Epi_scATAC/Results/02.Epi_AD_Methylation")
load("02.Epi_AD_Methylation.RData")
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

write.table(
    as.data.frame(diff_peaks$AD_vs_NA$Up)[, 1:3],
    file = "T_vs_AD_ADhigh.bed",
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)
write.table(
    as.data.frame(diff_peaks$AD_vs_NA$Down)[, 1:3],
    file = "T_vs_AD_Thigh.bed",
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)

table(duplicated(diff_peaks$AD_vs_NA$Up$idx))
table(duplicated(diff_peaks$AD_vs_NA$Up$peak_id))

## 1.3. enrichment with genome annotations ----
# genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene, single.strand.genes.only = FALSE) %>%
#     unlist()
# strand(genes) <- "*"
# genes <- reduce(genes)
# intergenic <- gaps(genes)
# intergenic <- intergenic[strand(intergenic) == "*"]
# sum(width(intergenic)) / 1e9

# repeats <- readRDS("E:/LabWork/genome/hg38/repeats/repeats.all.gr.RDS")

# anno.list <- list(
#     intergenic = intergenic,
#     intragenic = genes,
#     promoter = GenomicFeatures::promoters(TxDb.Hsapiens.UCSC.hg38.knownGene),
#     exon = GenomicFeatures::exonicParts(TxDb.Hsapiens.UCSC.hg38.knownGene),
#     intron = GenomicFeatures::intronicParts(TxDb.Hsapiens.UCSC.hg38.knownGene),
#     CGI = ChIPseeker::readPeakFile("E:/LabWork/genome/hg38/hg38.CGI.bed")
# )
# anno.list <- c(anno.list, repeats[c("LINE", "LTR", "Satellite", "Simple_repeat", "SINE")])
# rm(intergenic, genes, repeats)

# anno.list <- sapply(anno.list, function(x) {
#     x <- x[seqnames(x) %in% paste("chr", c(1:22, "X"), sep = "")]
#     x <- reduce(x, ignore.strand = TRUE)
#     return(x)
# })

# sapply(anno.list, function(x) {
#     return(sum(width(x)))
# }) / 3e9

# dir.create("LOLA_Core/LOLA_Core/regions")
# for (one in names(anno.list)) {
#     temp <- as.data.frame(anno.list[[one]])[, 1:3]
#     write.table(temp, file.path("LOLA_Core/regions", paste(one, "bed", sep = ".")),
#         sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
#     )
# }

# run LOLA test
library(LOLA)
regionDB <- loadRegionDB("E:/LabWork/genome/hg38/LOLA_Core")
sapply(diff_peaks$AD_vs_NA, length)
locResults <- runLOLA(diff_peaks$AD_vs_NA, proj_Epi@peakSet, regionDB, cores = 1)
plot.data <- locResults %>%
    as.data.frame() %>%
    select(userSet, collection, pValueLog, oddsRatio, qValue, filename) %>%
    mutate(
        direction = ifelse(userSet == "Up", 1, -1),
        filename = gsub(".bed", "", filename)
    )

pdf("AD_T_diffpeak_enrichment2.pdf", 6, 4)
ggplot(plot.data) +
    geom_bar(aes(x = pValueLog * direction, y = filename, fill = oddsRatio * direction), stat = "identity") +
    facet_grid(rows = vars(collection), scales = "free_y", space = "free") +
    ggbreak::scale_x_break(breaks = c(-160, -55)) +
    ggbreak::scale_x_break(breaks = c(55, 315)) +
    scale_fill_gradient2(low = "#1774cd", high = "#cd2525", mid = "#f7f7f7", midpoint = 0) +
    scale_x_continuous(breaks = c(-170, -45, -30, -15, 0, 15, 30, 45, 320), limits = c(-175, 330)) +
    theme(axis.text.x.top = element_blank()) +
    theme_classic() +
    xlab("-log10(qValue)") +
    ylab("Genome Annotation")
dev.off()
rm(regionDB, plot.data)

## 1.4. motif enrichment of AD peaks ----
dir.create("Homer")
write.table(
    diff_peaks$AD_vs_NA$Up$peak_id %>%
        gsub("_", "\t", .) %>%
        as.data.frame() %>%
        mutate(
            "Name" = diff_peaks$AD_vs_NA$Up$peak_id,
            "length" = 500,
            "strand" = "."
        ),
    "Homer/marker.peak.AD_vs_NA.Up.bed",
    col.names = FALSE,
    sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
    diff_peaks$AD_vs_NA$Down$peak_id %>%
        gsub("_", "\t", .) %>%
        as.data.frame() %>%
        mutate(
            "Name" = diff_peaks$AD_vs_NA$Down$peak_id,
            "length" = 500,
            "strand" = "."
        ),
    "Homer/marker.peak.AD_vs_NA.Down.bed",
    col.names = FALSE,
    sep = "\t", quote = FALSE, row.names = FALSE
)

# /mnt/d/WSL2/homer/bin/findMotifsGenome.pl ./marker.peak.AD_vs_NA.Up.bed hg38 ./AD_Up/ -size 200
homer.res.ad <- list(
    "Up" = homer.parser(paste0("homer/AD_Up/knownResults.txt"),
        log.p.value = 50, log2.enrichment = 0.5
    ),
    "Down" = homer.parser(paste0("homer/AD_Down/knownResults.txt"),
        log.p.value = 50, log2.enrichment = 0.5
    )
)
pdf("homer/Dot.motif.AD.Up.pdf", 5, 4)
ggplot(homer.res.ad$Up, aes(x = Log2_Enrichment, y = log.p.value)) +
    geom_point(aes(color = Diff), size = 1) +
    geom_vline(xintercept = 0.5, linetype = "dashed") +
    geom_hline(yintercept = 50, linetype = "dashed") +
    scale_color_manual(values = c("none" = "grey", "up" = "red")) +
    ggrepel::geom_text_repel(aes(label = TF), size = 2, max.overlaps = 30) +
    theme_classic()
dev.off()
pdf("homer/Dot.motif.AD.Down.pdf", 5, 4)
ggplot(homer.res.ad$Down, aes(x = Log2_Enrichment, y = log.p.value)) +
    geom_point(aes(color = Diff), size = 1) +
    geom_vline(xintercept = 0.5, linetype = "dashed") +
    geom_hline(yintercept = 50, linetype = "dashed") +
    scale_color_manual(values = c("none" = "grey", "up" = "red")) +
    ggrepel::geom_text_repel(aes(label = TF), size = 2, max.overlaps = 30) +
    theme_classic()
dev.off()

## 1.5. serrated vs conventional adenoma ----
ad.markers <- list(
    "Serrated" = c(
        "ALDOB", "MUC5AC", "AQP5", "TACSTD2", "FSCN1", "TFF2", "ANXA10", "REG4",
        "MUC17", "S100P", "GSDMB", "GSDMD", "IL18", "MKD", "RARA", "RXRA", "AGRN"
    ),
    "Conventional" = c(
        "CLDN2", "CD44", "AXIN2", "RNF43", "TGFBI", "EPHB2", "TEAD2", "CDX2", "LGR5", "ASCL2"
    )
)
proj_Epi$Sample_Type <- paste(proj_Epi$Sample, proj_Epi$Epi_type, sep = "_")
selected <- table(proj_Epi$Sample_Type)
selected <- names(selected)[selected > 10]
selected <- setdiff(selected, "COAD09_Adenoma") %>%
    grep("Adenoma", ., value = TRUE)

p.list <- plotGroups(proj_Epi,
    groupBy = "Sample_Type",
    colorBy = "GeneScoreMatrix",
    name = unlist(ad.markers) %>% intersect(getFeatures(proj_Epi)),
    plotAs = "violin"
)
length(p.list)

patient.rename <- c("P04", "P11", "P12", "P13", "P19", "P20", "P28")
names(patient.rename) <- c("COAD06", "COAD16", "COAD17", "COAD18", "COAD24", "COAD25", "COAD34")

p.list.2 <- lapply(p.list, function(x) {
    temp <- x$data %>%
        filter(x %in% selected) %>%
        mutate(Patient = gsub("-nofacs|_Adenoma", "", x)) %>%
        mutate("Type" = ifelse(Patient == ("COAD16"), "Serrated", "Conventional"))
    temp$Patient <- patient.rename[temp$Patient]
    p <- ggplot(temp) +
        geom_violin(aes(x = Patient, y = y, fill = Type), alpha = 0.5, show.legend = FALSE) +
        facet_grid(cols = vars(Type), scales = "free_x", space = "free_x") +
        ylab("Gene acvitiy") +
        theme_classic()
    return(p)
})

for (one in names(p.list.2)) {
    p.list.2[[one]] <- p.list.2[[one]] + ggtitle(one)
}

pdf("violin.AD.Serrated.markers.pdf", 21, 12)
patchwork::wrap_plots(p.list.2, ncol = 7)
dev.off()

pdf("violin.AD.Serrated.markers.selected.pdf", 10, 5)
patchwork::wrap_plots(p.list.2[c("MUC5AC", "TFF2", "AGRN", "RARA", "CLDN2", "LGR5", "ASCL2", "CD44")], ncol = 4)
dev.off()

# 2. AD diff peak keep in cancer ----
## 2.1. overlap with cancer peaks & CGI ----
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

# with CGI
CGI.bed <- ChIPseeker::readPeakFile("E:/LabWork/genome/hg38/hg38.CGI.bed")
length(CGI.bed)

temp <- diff_peaks$AD_vs_NA$Down
v <- Venn(Weight = c(
    "00" = 0,
    "10" = sum(countOverlaps(CGI.bed, temp) == 0),
    "01" = sum(countOverlaps(temp, CGI.bed) == 0),
    "11" = sum(countOverlaps(temp, CGI.bed) > 0)
), numberOfSets = 2, SetNames = c("CGI", "AD_down"))
pdf("Venn.Down.peak.all.vs.CGI.pdf", 6, 6)
plot(v, doWeights = FALSE, show = list(Faces = FALSE))
dev.off()

temp <- peakset[temp$peak_id] %>% subset(peakType == "Promoter")
length(temp)
v <- Venn(Weight = c(
    "00" = 0,
    "10" = sum(countOverlaps(CGI.bed, temp) == 0),
    "01" = sum(countOverlaps(temp, CGI.bed) == 0),
    "11" = sum(countOverlaps(temp, CGI.bed) > 0)
), numberOfSets = 2, SetNames = c("CGI", "AD_down_promoter"))

pdf("Venn.Down.peak.promoter.vs.CGI.pdf", 6, 6)
plot(v, doWeights = FALSE, show = list(Faces = FALSE))
dev.off()

rm(v, v1, v2)

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
    peak = rep("Down", length(rownames(peakmat.type)))
    #peakType = peakset[rownames(peakmat.type)]$peakType
)

anno.row[diff_peaks$AD_vs_NA$Up$peak_id, "peak"] <- "Up"
table(anno.row$peak)

anno.col <- data.frame(
    row.names = c("Normal", "Adenoma", "Malignant"),
    Epi_type = c("Normal", "Adenoma", "Malignant")
)
mycolor$peak <- c("Up" = "#cd2525", "Down" = "#1774cd")

pdf("Heatmap.AD_vs_NA.all.type1.pdf", 5.5, 6)
heatmap.merge <- pheatmap(peakmat.type,
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

cluster.info$Epi_Group <- cluster.info$Epi_type %>% as.character()
cluster.info[rownames(temp), ]$Epi_Group <- temp$Epi_Group
cluster.info$Epi_Group <- factor(cluster.info$Epi_Group,
    levels = c("Normal", "Adenoma", "Group_1", "Group_2")
)
rm(temp, sePeaks)

anno.col <- cluster.info %>% select(c("Epi_type", "Epi_Group"))

plot.data <- peakmat.cluster[diff_peaks$AD_vs_NA$Up$peak_id, ]
plot.data <- plot.data[, order(cluster.info$Epi_Group)]

pdf("Heatmap.AD_vs_NA.cluster.Up.pdf", 6, 6)
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

pdf("Heatmap.AD_vs_NA.cluster.Down.pdf", 6, 6)
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

# patient level
proj_Epi$Sample_Type <- paste(proj_Epi$Sample, proj_Epi$Epi_type, sep = "_")

sePeaks <- readRDS("../01.All_scCNV/sePeaks.Sample_Type.rds")

selected <- table(proj_Epi$Sample_Type)
selected <- names(selected)[selected > 10]
selected <- setdiff(selected, "COAD09_Adenoma")

rownames(sePeaks) <- paste(rowData(sePeaks)$seqnames,
    rowData(sePeaks)$start,
    rowData(sePeaks)$end,
    sep = "_"
)
rownames(sePeaks@assays@data$PeakMatrix) <- rownames(sePeaks)

peakmat.patient <- sePeaks@assays@data$PeakMatrix
peakmat.patient <- peakmat.patient[c(diff_peaks$AD_vs_NA$Up$peak_id, diff_peaks$AD_vs_NA$Down$peak_id), ]
peakmat.patient <- peakmat.patient[!duplicated(rownames(peakmat.patient)), selected]
rm(sePeaks)

plot.data <- apply(peakmat.patient, 1, function(x) {
    x <- (x - mean(x)) / sd(x)
    return(x)
}) %>% t()
plot.data[plot.data > 1.5] <- 1.5
plot.data[plot.data < (-1.5)] <- -1.5

anno.col <- data.frame(
    row.names = colnames(plot.data),
    Epi_type = gsub("^.+?_", "", colnames(plot.data)) %>%
        factor(levels = c("Normal", "Adenoma", "Malignant")),
    temp = "none"
)
anno.col <- anno.col[order(anno.col$Epi_type, rnorm(nrow(anno.col))), ]

anno.row <- data.frame(
    row.names = rownames(plot.data),
    peak = ifelse(rownames(plot.data) %in% diff_peaks$AD_vs_NA$Up$peak_id, "Up", "Down"),
    change = "none"
)
anno.row$change <- log2(rowMeans(peakmat.patient[, grep("Malignant", colnames(peakmat.patient))]) /
    rowMeans(peakmat.patient[, grep("Normal", colnames(peakmat.patient))]))
anno.row$status <- "none"
anno.row$status[anno.row$change > 1] <- "Up"
anno.row$status[anno.row$change < -1] <- "Down"
table(anno.row$status)
anno.row$status <- ifelse(anno.row$status == anno.row$peak, "Keep", "Revert")

anno.row <- anno.row[order(anno.row$peak, anno.row$status, rowSums(plot.data[, anno.col$Epi_type != "Malignant"])), ]
plot.data <- plot.data[rownames(anno.row), rownames(anno.col)]

png("Heatmap.AD_vs_NA.patient.png", 600, 600)
pheatmap(plot.data,
    scale = "none",
    show_rownames = FALSE,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    annotation_col = anno.col %>% select(Epi_type),
    annotation_row = anno.row %>% select(peak, status),
    annotation_colors = c(mycolor, list(status = c("Keep" = "#81c6a0", "Revert" = "#f0bd3b")))
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

wilcox.test(
    cluster.info[cluster.info$Epi_type == "Normal", "Mean_AD_Up"],
    cluster.info[cluster.info$Epi_type == "Adenoma", "Mean_AD_Up"]
)
wilcox.test(
    cluster.info[cluster.info$Epi_type == "Adenoma", "Mean_AD_Up"],
    cluster.info[cluster.info$Epi_type == "Malignant", "Mean_AD_Up"]
)
pdf("Line.Cluster.AD_vs_NA.Up.pdf", 4, 3)
ggplot(plot.data, aes(x = Epi_type, y = Mean)) +
    geom_line(aes(group = 1), size = 1, color = "#cd2525") +
    geom_ribbon(aes(x = 1:3, ymin = Mean - SD, ymax = Mean + SD),
        alpha = 0.2, fill = "#cd2525"
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
pdf("Line.Cluster.AD_vs_NA.Down.pdf", 4, 3)
ggplot(plot.data, aes(x = Epi_type, y = Mean)) +
    geom_line(aes(group = 1), size = 1, color = "#1774cd") +
    geom_ribbon(aes(x = 1:3, ymin = Mean - SD, ymax = Mean + SD),
        alpha = 0.2, fill = "#1774cd"
    ) +
    geom_point(aes(color = Epi_type), size = 2, show.legend = FALSE) +
    theme_classic() +
    scale_color_manual(values = mycolor$Epi_type) +
    scale_x_discrete(expand = c(0.02, 0.05)) +
    scale_y_continuous(labels = c(0.1, 0.2, 0.3), breaks = c(0.1, 0.2, 0.3)) +
    ylab("Mean ATAC signal")
dev.off()
rm(plot.data)

## 2.4. percent of tumor peaks in adenoma ----
cluster.info <- readRDS("../03.Epi_Molecular_Subtype/cluster.info.rds")

peaks.clusters.up <- list()
peaks.clusters.down <- list()
sapply(rownames(cluster.info), function(one) {
    temp <- read.table(paste0("../03.Epi_Molecular_Subtype/diff_peak_cluster/marker.peak.tumor.", one, ".tsv"),
        header = TRUE
    ) %>%
        mutate(
            id = paste(seqnames, start, end, sep = "_")
        )
    peaks.clusters.up[[one]] <<- temp$id
    return(1)
})
peaks.clusters.down <- sapply(rownames(cluster.info), function(one) {
    temp <- read.table(paste0("../03.Epi_Molecular_Subtype/diff_peak_cluster/down/marker.peak.tumor.", one, ".tsv"),
        header = TRUE, sep = "\t", stringsAsFactors = FALSE
    )
    return(paste(temp$seqnames, temp$start, temp$end, sep = "_"))
})

plot.data <- data.frame(
    row.names = rownames(cluster.info),
    "Cluster" = rownames(cluster.info),
    "n_Up" = sapply(peaks.clusters.up, length),
    "n_Down" = sapply(peaks.clusters.down, length),
    "n_Up_Adenoma" = sapply(peaks.clusters.up, function(x) {
        return(sum(x %in% diff_peaks$AD_vs_NA$Up$peak_id))
    }),
    "n_Down_Adenoma" = sapply(peaks.clusters.up, function(x) {
        return(sum(x %in% diff_peaks$AD_vs_NA$Down$peak_id))
    })
)
(plot.data$n_Up / plot.data$n_Down) %>% mean() # 2.474568
(plot.data$n_Up > plot.data$n_Down) %>% sum() # 22

plot.data$n_Up_Malignant <- plot.data$n_Up - plot.data$n_Up_Adenoma
plot.data$n_Down_Malignant <- plot.data$n_Down - plot.data$n_Down_Adenoma
plot.data$n_Up <- NULL
plot.data$n_Down <- NULL

plot.data <- reshape2::melt(plot.data, id.vars = "Cluster")
plot.data$value[grep("Down", plot.data$variable)] <- plot.data$value[grep("Down", plot.data$variable)] * -1
plot.data$variable <- factor(plot.data$variable,
    levels = c("n_Up_Malignant", "n_Up_Adenoma", "n_Down_Malignant", "n_Down_Adenoma")
)

cluster.rename <- read.table("../01.All_scCNV/cluster_rename.txt", header = TRUE)
rownames(cluster.rename) <- cluster.rename$Cluster
plot.data$Cluster_rename <- cluster.rename[plot.data$Cluster, "manual"]
plot.data$Cluster_rename <- factor(plot.data$Cluster_rename,
    levels = paste("C", 5:29, sep = "")
)
table(plot.data$Cluster_rename)

pdf("Barplot.cluster.tumor.peak.in.adenoma.pdf", 6, 3)
ggplot(plot.data, aes(x = Cluster_rename, y = value / 1000, fill = variable)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_classic() +
    scale_fill_manual(values = c("#cd2525", "#f89797", "#1774cd", "#75b6f4")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Number of peaks (10^3)") +
    xlab("Cluster_rename")
dev.off()

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
    Methylation = "not significant",
    peak = rep("Down", length(rownames(betamat.peak)))
)

anno.row[rownames(anno.row) %in% diff_peaks$AD_vs_NA$Up$peak_id, "peak"] <- "Up"
table(anno.row$peak)
anno.row$peak <- factor(anno.row$peak,
    levels = c("Up", "Down")
)

temp <- data.frame(
    row.names = rownames(betamat.peak),
    beta_normal = rowMeans(betamat.peak[, sample.info.array$ID[sample.info.array$Location == "Normal"]], na.rm = TRUE),
    beta_adenoma = rowMeans(betamat.peak[, sample.info.array$ID[sample.info.array$Location == "Adenoma"]], na.rm = TRUE)
)
anno.row$Methylation[temp$beta_adenoma - temp$beta_normal > 0.1] <- "Increased"
anno.row$Methylation[temp$beta_adenoma - temp$beta_normal < -0.1] <- "Decreased"
table(anno.row$Methylation)
anno.row$Methylation <- factor(anno.row$Methylation,
    levels = c("Increased", "Decreased", "not significant")
)
table(anno.row$peak, anno.row$Methylation)

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
        "Methylation" = c("Increased" = "#cd2525", "Decreased" = "#1774cd"),
        "peak" = mycolor$peak,
    )
)

# plot.data <- betamat.peak[p$tree_row$order, p$tree_col$order]
plot.data <- betamat.peak[
    order(
        anno.row$peak,
        # anno.row$Methylation,
        rowMeans(betamat.peak[, sample.info.array[colnames(betamat.peak), ]$Location == "Normal"], na.rm = TRUE) *
            ifelse(anno.row$peak == "Up", -1, 1)
    ),
    p$tree_col$order
]
plot.data <- plot.data[order(rownames(plot.data) %in% diff_peaks$AD_vs_NA$Up$peak_id), ]
plot.data <- plot.data[, order(sample.info.array[colnames(plot.data), ]$Location)]
table(rownames(plot.data) %in% diff_peaks$AD_vs_NA$Up$peak_id)
# plot.data <- plot.data[c(1:1257, rev(1258:4736)), ]

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
        "Methylation" = c("Increased" = "#fbb025", "Decreased" = "#44aaee", "not significant" = "gray50")
    ),
    gaps_row = c(1257)
)
dev.off()

rm(anno.col, anno.row, plot.data, p)

## 3.3. Line chart of aggregated peaks ----
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

wilcox.test(
    sample.info.array[sample.info.array$Location == "Normal", "Mean_AD_Up"],
    sample.info.array[sample.info.array$Location == "Adenoma", "Mean_AD_Up"]
)
wilcox.test(
    sample.info.array[sample.info.array$Location == "Adenoma", "Mean_AD_Up"],
    sample.info.array[sample.info.array$Location == "Malignant", "Mean_AD_Up"]
)

pdf("Line.beta.Cluster.AD_vs_NA.Up.pdf", 4, 3)
ggplot(plot.data, aes(x = Location, y = Mean)) +
    geom_line(aes(group = 1), size = 1, lty = 2, color = "#cd2525") +
    geom_ribbon(aes(x = 1:3, ymin = Mean - SD, ymax = Mean + SD),
        alpha = 0.2, fill = "#cd2525"
    ) +
    geom_point(aes(color = Location), size = 2, show.legend = FALSE) +
    theme_classic() +
    scale_color_manual(values = mycolor$Epi_type) +
    scale_x_discrete(expand = c(0.02, 0.05)) +
    scale_y_continuous(labels = c(0.4, 0.5, 0.6), breaks = c(0.4, 0.5, 0.6)) +
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
wilcox.test(
    sample.info.array[sample.info.array$Location == "Normal", "Mean_AD_Down"],
    sample.info.array[sample.info.array$Location == "Adenoma", "Mean_AD_Down"]
)
wilcox.test(
    sample.info.array[sample.info.array$Location == "Adenoma", "Mean_AD_Down"],
    sample.info.array[sample.info.array$Location == "Malignant", "Mean_AD_Down"]
)

pdf("Line.beta.Cluster.AD_vs_NA.Down.pdf", 4, 3)
ggplot(plot.data, aes(x = Location, y = Mean)) +
    geom_line(aes(group = 1), size = 1, lty = 2, color = "#1774cd") +
    geom_ribbon(aes(x = 1:3, ymin = Mean - SD, ymax = Mean + SD),
        alpha = 0.2, fill = "#1774cd"
    ) +
    geom_point(aes(color = Location), size = 2, show.legend = FALSE) +
    theme_classic() +
    scale_color_manual(values = mycolor$Epi_type) +
    scale_x_discrete(expand = c(0.02, 0.05)) +
    ylab("Mean beta value")
dev.off()
rm(plot.data)

# pdf("Beasworm.AD_vs_NA.Up.pdf", 4, 3.5)
# ggplot(sample.info.array, aes(x = Location, y = Mean_AD_Up)) +
#     geom_beeswarm(aes(color = Location), cex = 2.5, show.legend = FALSE) +
#     geom_violin(aes(fill = Location), alpha = 0.2, show.legend = FALSE) +
#     scale_color_manual(values = mycolor$Epi_type) +
#     ggpubr::stat_compare_means(
#         comparisons = list(
#             c("Normal", "Adenoma"),
#             c("Normal", "Malignant"),
#             c("Adenoma", "Malignant")
#         )
#     ) +
#     ylab("Mean beta value") +
#     scale_fill_manual(values = mycolor$Epi_type) +
#     theme_classic()
# dev.off()
# pdf("Beasworm.AD_vs_NA.Down.pdf", 4, 3.5)
# ggplot(sample.info.array, aes(x = Location, y = Mean_AD_Down)) +
#     geom_beeswarm(aes(color = Location), cex = 2.5, show.legend = FALSE) +
#     geom_violin(aes(fill = Location), alpha = 0.2, show.legend = FALSE) +
#     scale_color_manual(values = mycolor$Epi_type) +
#     scale_fill_manual(values = mycolor$Epi_type) +
#     ggpubr::stat_compare_means(
#         comparisons = list(
#             c("Normal", "Adenoma"),
#             c("Normal", "Malignant"),
#             c("Adenoma", "Malignant")
#         )
#     ) +
#     ylab("Mean beta value") +
#     theme_classic()
# dev.off()

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
        geom_quasirandom(aes(color = Location), cex = 1.5, show.legend = FALSE) +
        # geom_violin(aes(fill = Location), alpha = 0.2, show.legend = FALSE) +
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

pdf("Beasworm.AD_vs_NA.selected.gene1.pdf", 8, 10)
patchwork::wrap_plots(p.list, ncol = 3)
dev.off()
rm(p.list, plot.data, p, one)

# 4. Patient paired analysis ----
dir.create("Paired_Sample")
patient.selected <- c("COAD06", "COAD16", "COAD17", "COAD18", "COAD24", "COAD34")

## 4.1. get peak * type_sample matrix ----
table(proj_Epi$Sample, proj_Epi$Epi_type)

proj_Epi$Epi_type_Sample <- paste(proj_Epi$Epi_type,
    gsub("-nofacs", "", proj_Epi$Sample),
    sep = "_"
)
proj_Epi$Epi_type_Sample <- gsub("Normal_.*", "Normal", proj_Epi$Epi_type_Sample)
table(proj_Epi$Epi_type_Sample)

# patient * cell type level
sePeaks.epi.type <- getGroupSE(
    ArchRProj = proj_Epi,
    useMatrix = "PeakMatrix",
    groupBy = "Epi_type_Sample",
    divideN = TRUE,
    scaleTo = NULL
)
rownames(sePeaks.epi.type) <- paste(rowData(sePeaks.epi.type)$seqnames,
    rowData(sePeaks.epi.type)$start,
    rowData(sePeaks.epi.type)$end,
    sep = "_"
)
rownames(sePeaks.epi.type@assays@data$PeakMatrix) <- rownames(sePeaks.epi.type)

peakmat.type.sample <- sePeaks.epi.type@assays@data$PeakMatrix
peakmat.type.sample <- peakmat.type.sample[, grep(
    "Normal|COAD06|COAD16|COAD17|COAD18|COAD24|COAD34",
    colnames(peakmat.type.sample)
)]

## 4.2. identify Ad diff peak by patients ----
diff_peaks_patient_AD <- list()
for (one in patient.selected) {
    marker.peak.one <- getMarkerFeatures(
        ArchRProj = proj_Epi,
        useMatrix = "PeakMatrix",
        groupBy = "Epi_type_Sample",
        testMethod = "wilcoxon",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        useGroups = paste("Adenoma", one, sep = "_"),
        bgdGroups = "Normal"
    )

    pv <- plotMarkers(
        seMarker = marker.peak.one,
        name = paste("Adenoma", one, sep = "_"),
        cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1",
        plotAs = "Volcano"
    )
    pdf(paste("Paired_Sample/Volcano.markers.AD_vs_NA.", one, ".pdf", sep = ""), 5, 4)
    plot(pv)
    dev.off()

    marker.peak.one.df <- getMarkers(marker.peak.one,
        cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1"
    )
    marker.peak.one.df <- marker.peak.one.df[[1]] %>%
        as.data.frame() %>%
        mutate(
            id = paste(seqnames, start, end, sep = "_"),
            peak = ifelse(Log2FC > 0, "Up", "Down")
        )
    diff_peaks_patient_AD[[one]] <- split(marker.peak.one.df$id, marker.peak.one.df$peak)
}
rm(one, marker.peak.one, marker.peak.one.df, pv, sePeaks.epi.type)

temp <- unlist(diff_peaks_patient_AD, recursive = TRUE)
peakmat.type.sample <- peakmat.type.sample[temp, ]

temp[grep("Up", names(temp))] %>% table() %>% table()
temp[grep("Down", names(temp))] %>% table() %>% table()
rm(temp)

## 4.3. identify Ca diff peak by patients ----
diff_peaks_patient_Ca <- list()
for (one in patient.selected) {
    marker.peak.one <- getMarkerFeatures(
        ArchRProj = proj_Epi,
        useMatrix = "PeakMatrix",
        groupBy = "Epi_type_Sample",
        testMethod = "wilcoxon",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        useGroups = paste("Malignant", one, sep = "_"),
        bgdGroups = "Normal"
    )

    pv <- plotMarkers(
        seMarker = marker.peak.one,
        name = paste("Malignant", one, sep = "_"),
        cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1",
        plotAs = "Volcano"
    )
    pdf(paste("Paired_Sample/Volcano.markers.Ca_vs_NA.", one, ".pdf", sep = ""), 5, 4)
    plot(pv)
    dev.off()

    marker.peak.one.df <- getMarkers(marker.peak.one,
        cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1"
    )
    marker.peak.one.df <- marker.peak.one.df[[1]] %>%
        as.data.frame() %>%
        mutate(
            id = paste(seqnames, start, end, sep = "_"),
            peak = ifelse(Log2FC > 0, "Up", "Down")
        )
    diff_peaks_patient_Ca[[one]] <- split(marker.peak.one.df$id, marker.peak.one.df$peak)
}
rm(one, marker.peak.one, marker.peak.one.df, pv)

# Venn of AD / Ca peaks by patient
one <- patient.selected[2]
for (one in patient.selected[2:6]) {
    v <- Venn(list(
        Adenoma = unlist(diff_peaks_patient_AD[[one]]$Up),
        Malignant = unlist(diff_peaks_patient_Ca[[one]]$Up)
    ))
    pdf(paste("Paired_Sample/Venn.markers.Up.", one, ".pdf", sep = ""), 5, 4)
    gridExtra::grid.arrange(grid::grid.grabExpr(plot(v,
        doWeights = FALSE,
        show = list(Faces = FALSE)
    )), top = paste(one, "Up", sep = " "))
    dev.off()

    v <- Venn(list(
        Adenoma = unlist(diff_peaks_patient_AD[[one]]$Down),
        Malignant = unlist(diff_peaks_patient_Ca[[one]]$Down)
    ))
    pdf(paste("Paired_Sample/Venn.markers.Down.", one, ".pdf", sep = ""), 5, 4)
    gridExtra::grid.arrange(grid::grid.grabExpr(plot(v,
        doWeights = FALSE,
        show = list(Faces = FALSE)
    )), top = paste(one, "Down", sep = " "))
    dev.off()
}
rm(v, one)

## 4.4. heatmap of AD peaks by patient ----
anno.col <- data.frame(
    row.names = colnames(peakmat.type.sample),
    Epi_type = gsub("_.*", "", colnames(peakmat.type.sample))
)
diff_peaks_patient_AD_state <- list()

for (one in patient.selected[c(2, 4:6)]) {
    one <- patient.selected[c(6)]
    selected <- grep(paste("Normal", one, sep = "|"),
        colnames(peakmat.type.sample),
        value = TRUE
    ) %>%
        sort() %>%
        .[c(3, 1, 2)]
    plot.data <- peakmat.type.sample[unlist(diff_peaks_patient_AD[[one]]), selected]

    anno.row <- data.frame(
        row.names = rownames(plot.data),
        change = log2(plot.data[, selected[3]] / plot.data[, selected[1]]),
        peak = gsub("[0-9]+", "", names(unlist(diff_peaks_patient_AD[[one]])))
    )

    anno.row$status <- "none"
    anno.row$status[anno.row$change > 1] <- "Up"
    anno.row$status[anno.row$change < -1] <- "Down"
    anno.row$status <- ifelse(anno.row$status == anno.row$peak, "Keep", "Revert")
    anno.row$change <- NULL

    # diff_peaks_patient_AD_state[[one]] <- table(anno.row$status, anno.row$peak)
    p <- pheatmap(plot.data,
        color = colorRampPalette(rev(brewer.pal(7, "RdYlBu"))[2:6])(100),
        scale = "row",
        show_rownames = FALSE,
        clustering_method = "ward.D2",
        cluster_cols = FALSE,
        cluster_rows = TRUE,
        annotation_col = anno.col,
        annotation_row = anno.row[, c(2, 1)],
        annotation_colors = c(mycolor, list(status = c("Keep" = "#81c6a0", "Revert" = "#f0bd3b")))
    )
    plot.data <- plot.data[p$tree_row$order, ]
    plot.data <- plot.data[order(
        anno.row[rownames(plot.data), "peak"],
        anno.row[rownames(plot.data), "status"]
    ), ]

    while (dev.cur() != 1) {
        dev.off()
    }
    pdf(paste("Paired_Sample/Heatmap.AD_vs_NA.", one, ".pdf", sep = ""), 5.5, 6)
    pheatmap(plot.data,
        color = colorRampPalette(rev(brewer.pal(7, "RdYlBu"))[2:6])(100),
        scale = "row",
        show_rownames = FALSE,
        clustering_method = "ward.D2",
        cluster_cols = FALSE,
        cluster_rows = FALSE,
        annotation_col = anno.col,
        annotation_row = anno.row[, c(2, 1)],
        annotation_colors = c(mycolor, list(status = c("Keep" = "#81c6a0", "Revert" = "#f0bd3b")))
    )
    dev.off()
}
rm(selected, plot.data, anno.col, anno.row, one)

temp <- diff_peaks_patient_AD_state$COAD16 +
    diff_peaks_patient_AD_state$COAD18 +
    diff_peaks_patient_AD_state$COAD24 +
    diff_peaks_patient_AD_state$COAD34
temp[1] / (temp[1] + temp[2]) # 0.5586728 of down peaks
temp[3] / (temp[3] + temp[4]) # 0.7741892 of Up peaks
sapply(diff_peaks_patient_AD_state[c("COAD16", "COAD18", "COAD24", "COAD34")], function(x) {
    return(c(x[1] / (x[1] + x[2]), x[3] / (x[3] + x[4])))
})
rm(temp)

## 4.5. Dot visulization of methlation-related genes ----
p <- plotGroups(
    ArchRProj = proj_Epi,
    groupBy = "Epi_type",
    colorBy = "GeneScoreMatrix",
    name = gene.selected$gene,
    plotAs = "violin",
    maxCells = 8000
)

p.list <- list()
for (one in names(p)) {
    plot.data <- p[[one]]$data
    colnames(plot.data) <- c("Epi_type", "Value")
    plot.data$Sample <- as.data.frame(proj_Epi@cellColData)[rownames(plot.data), ]$Sample %>%
        gsub("-nofacs", "", .)
    plot.data$Gene <- one

    plot.data$Epi_type <- factor(plot.data$Epi_type,
        levels = c("Normal", "Adenoma", "Malignant")
    )
    plot.data <- plot.data[plot.data$Sample %in% patient.selected[c(2, 4:6)], ]

    p.list[[one]] <- ggplot(plot.data, aes(x = Epi_type, y = Value)) +
        geom_quasirandom(aes(color = Epi_type), cex = 0.5, show.legend = FALSE) +
        scale_color_manual(values = mycolor$Epi_type) +
        scale_fill_manual(values = mycolor$Epi_type) +
        ggpubr::stat_compare_means(
            comparisons = list(
                c("Normal", "Adenoma"),
                c("Normal", "Malignant"),
                c("Adenoma", "Malignant")
            ),
            size = 2, label.y.npc = 0.5,
        ) +
        ylab("Gene Activity Score") +
        # ggtitle(one) +
        theme_classic() +
        facet_grid(cols = vars(Sample), rows = vars(Gene), scales = "free_y") +
        theme(axis.title.x = element_blank(), axis.text.x = element_blank())
}
pdf("Paired_Sample/Beasworm.AD_vs_NA.Up.patient.pdf", 6, 10)
patchwork::wrap_plots(p.list[1:6], ncol = 1)
dev.off()
pdf("Paired_Sample/Beasworm.AD_vs_NA.Down.patient.pdf", 6, 10)
patchwork::wrap_plots(p.list[7:12], ncol = 1)
dev.off()
rm(p, plot.data, p.list, one)

gc()
save.image("03.Epi_AD_Methylation.RData")
