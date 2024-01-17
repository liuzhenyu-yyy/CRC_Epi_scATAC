
## 6.3. co-binding of TFs ----
pdf("TF_motif_cluster/Upset.Consensus.Peaks.pdf", 7, 4.5)
upset(
    fromList(
        peaks.common
    ),
    nset = 4,
    order.by = "freq"
)
dev.off()

for (one in c("Group_1", "Group_2")) {
    cluster.selected <- cluster.info %>%
        filter(Epi_Group == one) %>%
        rownames()

    for (genes in list(c("HNF4A", "PPARA"), c("SOX2", "FOXA3"))) {
        peak.selected1 <- motif.match[, genes[1]] %>%
            .[.] %>%
            names() %>%
            intersect(., unlist(peaks.clusters[cluster.selected]) %>% unique())
        peak.selected2 <- motif.match[, genes[2]] %>%
            .[.] %>%
            names() %>%
            intersect(., unlist(peaks.clusters[cluster.selected]) %>% unique())
        temp <- list(gene1 = peak.selected1, genes2 = peak.selected2)
        names(temp) <- c(genes[1], genes[2])
        v <- Venn(temp)

        pdf(paste("TF_motif_cluster/Venn.", one, "_", genes[1], "_", genes[2], ".pdf", sep = ""), 5, 4)
        gridExtra::grid.arrange(grid::grid.grabExpr(plot(v,
            doWeights = TRUE,
            show = list(Faces = FALSE)
        )), top = paste(genes[1], genes[2], one, sep = " "))
        dev.off()
    }
}
rm(one, genes, cluster.selected, peak.selected1, peak.selected2, temp, v, motif.match)

## 6.4. consensus peak genes ----
# get consensus peaks to genes
sapply(peaks.common, length)
peakset <- proj_Epi@peakSet
names(peakset) <- paste(seqnames(peakset), start(peakset), end(peakset), sep = "_")
table(peakset$peakType)

peaks.common.2gene <- sapply(peaks.common, function(one) {
    temp <- peakset[one]
    temp <- temp[temp$peakType != "Distal"]
    temp <- unique(temp$nearestGene)
    temp <- temp[!is.na(temp)]
    return(temp)
})
sapply(peaks.common.2gene, length)

length(unlist(peaks.common.2gene) %>% unique()) # 1582 geness
rm(peakset)
saveRDS(peaks.common.2gene, "peaks.common.2gene.rds")

# overlap with iCMS markers
markers.iCMS <- read.table("E:/LabWork/Project/CRC_NGS_ATAC/iCMS markers.txt",
    stringsAsFactors = FALSE,
    sep = "\t", header = TRUE
)
markers.iCMS <- base::as.list(markers.iCMS)
genes <- getFeatures(proj_Epi, useMatrix = "GeneScoreMatrix")
markers.iCMS <- lapply(markers.iCMS, function(x) intersect(x, genes))
lapply(markers.iCMS, length)

for (TF in names(peaks.common.2gene)) {
    pdf(paste0("TF_motif_cluster/Venn.peaks.common.", TF, ".iCMS.pdf"), 5, 4)
    for (genes in names(markers.iCMS)) {
        temp <- list(target = peaks.common.2gene[[TF]], iCMS = markers.iCMS[[genes]])
        names(temp) <- c(TF, genes)
        v <- Vennerable::Venn(temp)
        plot(v, doWeights = TRUE, show = list(Faces = FALSE))
    }
    dev.off()
}
rm(TF, genes, temp, v, markers.iCMS)

# get cell by concensus gene matrix
cGSMat <- getMatrixFromProject(
    ArchRProj = proj_Epi,
    useMatrix = "GeneScoreMatrix",
    binarize = FALSE
)
dim(cGSMat)
rownames(cGSMat@assays@data$GeneScoreMatrix) <- rowData(cGSMat)$name
cGSMat <- cGSMat@assays@data$GeneScoreMatrix
cGSMat <- as.data.frame(cGSMat)
cGSMat <- cGSMat[unlist(peaks.common.2gene) %>% unique(), ]

# reduction on consensus peak-genes
pca.cGSMat <- prcomp(t(cGSMat)[rownames(sample.info.tumor), ],
    center = TRUE, scale. = TRUE
)

pdf("TF_motif_cluster/Elbow.pca.cGSMat.pdf", 5, 5)
plot(pca.cGSMat$sdev[1:50]^2)
dev.off()

sample.info.tumor$PCA_cGS_1 <- pca.cGSMat$x[rownames(sample.info.tumor), 1]
sample.info.tumor$PCA_cGS_2 <- pca.cGSMat$x[rownames(sample.info.tumor), 2]
sample.info.tumor$PCA_cGS_3 <- pca.cGSMat$x[rownames(sample.info.tumor), 3]

pdf("TF_motif_cluster/PCA_cGSMat.Epi_Group.pdf", 12, 3.5)
ggplot(sample.info.tumor, aes(x = PCA_cGS_1, y = PCA_cGS_2)) +
    geom_point(aes(color = Epi_Group), size = 0.2) +
    scale_color_manual(values = mycolor$Epi_Group) +
    theme_classic() +
    NoLegend() +
    ggplot(sample.info.tumor, aes(x = PCA_cGS_1, y = PCA_cGS_3)) +
    geom_point(aes(color = Epi_Group), size = 0.2) +
    scale_color_manual(values = mycolor$Epi_Group) +
    theme_classic() +
    NoLegend() +
    ggplot(sample.info.tumor, aes(x = PCA_cGS_2, y = PCA_cGS_3)) +
    geom_point(aes(color = Epi_Group), size = 0.2) +
    scale_color_manual(values = mycolor$Epi_Group) +
    theme_classic()
dev.off()
pdf("TF_motif_cluster/PCA_cGS.CIMP_Group.pdf", 12, 3.5)
ggplot(sample.info.tumor, aes(x = PCA_cGS_1, y = PCA_cGS_2)) +
    geom_point(aes(color = CIMP_Group), size = 0.2) +
    scale_color_manual(values = mycolor$CIMP_Group) +
    theme_classic() +
    NoLegend() +
    ggplot(sample.info.tumor, aes(x = PCA_cGS_1, y = PCA_cGS_3)) +
    geom_point(aes(color = CIMP_Group), size = 0.2) +
    scale_color_manual(values = mycolor$CIMP_Group) +
    theme_classic() +
    NoLegend() +
    ggplot(sample.info.tumor, aes(x = PCA_cGS_2, y = PCA_cGS_3)) +
    geom_point(aes(color = CIMP_Group), size = 0.2) +
    scale_color_manual(values = mycolor$CIMP_Group) +
    theme_classic()
dev.off()

umap.cGSMat <- uwot::umap(pca.cGSMat$x[, 1:15],
    n_neighbors = 30,
    metric = "cosine",
    min_dist = 1
)
sample.info.tumor$UMAP_cGS_1 <- umap.cGSMat[rownames(sample.info.tumor), 1]
sample.info.tumor$UMAP_cGS_2 <- umap.cGSMat[rownames(sample.info.tumor), 2]

pdf("TF_motif_cluster/UMAP_cGS.Epi_Group.pdf", 4.5, 4)
ggplot(sample.info.tumor, aes(x = UMAP_cGS_1, y = UMAP_cGS_2)) +
    geom_point(aes(color = Epi_Group), size = 0.2) +
    scale_color_manual(values = mycolor$Epi_Group) +
    theme_bw() +
    theme(
        axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    ) +
    coord_fixed()
dev.off()
pdf("TF_motif_cluster/UMAP_cGS.CIMP_Group.pdf", 5, 4)
ggplot(sample.info.tumor, aes(x = UMAP_cGS_1, y = UMAP_cGS_2)) +
    geom_point(aes(color = CIMP_Group), size = 0.2) +
    scale_color_manual(values = mycolor$CIMP_Group) +
    theme_bw() +
    theme(
        axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    ) +
    coord_fixed()
dev.off()
pdf("TF_motif_cluster/UMAP_cGS.Sample.pdf", 7, 6)
ggplot(sample.info.tumor, aes(x = UMAP_cGS_1, y = UMAP_cGS_2)) +
    geom_point(aes(color = Sample), size = 0.2) +
    scale_color_manual(values = paletteDiscrete(sample.info.tumor$Sample)) +
    theme_bw() +
    theme(
        axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    ) +
    coord_fixed()
dev.off()
rm(cGSMat, pca.cGSMat, umap.cGSMat)

# 7. TF target in adenema? paired sample ----
## 7.1. get data: paired diff peak & motif match ----
proj_Epi$Patient_Epi <- paste(proj_Epi$Patient, proj_Epi$Epi_Group, sep = "_")
proj_Epi$Patient_Epi[grep("Normal", proj_Epi$Epi_Group)] <- "Normal"
table(proj_Epi$Epi_Group[proj_Epi$Patient == "COAD24"])

# only paired samples
marker.paired.list <- list()
for (patient in c("COAD16", "COAD18", "COAD24", "COAD34")) {
    for (one in c("Group_2", "Adenoma")) {
        if (patient == "COAD34" && one == "Group_2") {
            one <- "Group_1"
        }
        marker.peak.one <- getMarkerFeatures(
            ArchRProj = proj_Epi,
            useMatrix = "PeakMatrix",
            groupBy = "Patient_Epi",
            testMethod = "wilcoxon",
            bias = c("TSSEnrichment", "log10(nFrags)"),
            useGroups = paste(patient, one, sep = "_"),
            bgdGroups = c("Normal")
        )

        marker.peak.one.df <- getMarkers(marker.peak.one,
            cutOff = "FDR <= 0.01 & Log2FC >= 1"
        )[[1]] %>% as.data.frame()

        marker.paired.list[[paste(patient, one, sep = "_")]] <-
            paste(marker.peak.one.df$seqnames,
                marker.peak.one.df$start,
                marker.peak.one.df$end,
                sep = "_"
            )
    }
}
saveRDS(marker.paired.list, "marker.paired.list.RDS")
marker.paired.list <- readRDS("marker.paired.list.RDS")
sapply(marker.paired.list, length)

# motif match matrix
motif.match <- getMatches(ArchRProj = proj_Epi, name = "Motif")
temp <- rowRanges(motif.match)
motif.match <- motif.match@assays@data$matches

rownames(motif.match) <- paste(seqnames(temp), start(temp), end(temp), sep = "_")
colnames(motif.match) <- gsub("_.+?$", "", colnames(motif.match))
motif.match[1:5, 1:5]
rm(temp)

temp <- peak.selected <- motif.match[, "PPARA"] %>%
    .[.] %>%
    names()

one <- intersect(marker.paired.list$COAD24_Adenoma, temp)
sum(one %in% peaks.common[["PPARA"]]) / length(one)

## 7.2. iCMS2 TFs ----
cluster.selected <- cluster.info %>%
    filter(Epi_Group == "Group_2") %>%
    rownames()
marker.paired.list.sub <- marker.paired.list[c(1:6)]

# HNF4A
peak.selected <- motif.match[, "HNF4A"] %>%
    .[.] %>%
    names() %>%
    intersect(., c(
        unlist(peaks.clusters[cluster.selected]),
        unlist(marker.paired.list.sub)
    ) %>% unique())
length(peak.selected) # 20724 HNF4A peaks,AD + CRC

plot.data1 <- matrix(0, ncol = length(peak.selected), nrow = length(marker.paired.list.sub))
colnames(plot.data1) <- peak.selected
rownames(plot.data1) <- names(marker.paired.list.sub)
plot.data1 <- as.data.frame(plot.data1)

for (patient in c("COAD16", "COAD18", "COAD24")) {
    for (one in c("Group_2", "Adenoma")) {
        temp <- paste(patient, one, sep = "_")
        plot.data1[temp, colnames(plot.data1) %in% marker.paired.list.sub[[temp]]] <- 1
    }
}
plot.data1 <- plot.data1[c(2, 4, 6, 1, 3, 5), ]

plot.data2 <- matrix(0, ncol = length(peak.selected), nrow = length(cluster.selected))
colnames(plot.data2) <- peak.selected
rownames(plot.data2) <- cluster.selected
plot.data2 <- as.data.frame(plot.data2)

for (one in cluster.selected) {
    plot.data2[one, peak.selected %in% peaks.clusters[[one]]] <- 1
}

rownames(plot.data2) <- paste(cluster.selected, cluster.info[cluster.selected, "Patient_Major"], sep = "_")
plot.data <- rbind(plot.data1, plot.data2)

pdf("TF_motif_cluster/Heatmap.Paired.HNF4A.PatientPeaks.pdf", 10, 6)
pheatmap(
    plot.data,
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_color = mycolor,
    color = colorRampPalette(c("gray80", "#cd2525"))(2),
    clustering_method = "ward.D2",
    cutree_cols = 12
)
dev.off()

# PPARA
peak.selected <- motif.match[, "PPARA"] %>%
    .[.] %>%
    names() %>%
    intersect(., c(
        unlist(peaks.clusters[cluster.selected]),
        unlist(marker.paired.list.sub)
    ) %>% unique())
length(peak.selected) # 10081 HNF4A peaks,AD + CRC

plot.data1 <- matrix(0, ncol = length(peak.selected), nrow = length(marker.paired.list.sub))
colnames(plot.data1) <- peak.selected
rownames(plot.data1) <- names(marker.paired.list.sub)
plot.data1 <- as.data.frame(plot.data1)

for (patient in c("COAD16", "COAD18", "COAD24")) {
    for (one in c("Group_2", "Adenoma")) {
        temp <- paste(patient, one, sep = "_")
        plot.data1[temp, colnames(plot.data1) %in% marker.paired.list.sub[[temp]]] <- 1
    }
}
rowSums(plot.data1)
plot.data1 <- plot.data1[c(2, 4, 6, 1, 3, 5), ]

plot.data2 <- matrix(0, ncol = length(peak.selected), nrow = length(cluster.selected))
colnames(plot.data2) <- peak.selected
rownames(plot.data2) <- cluster.selected
plot.data2 <- as.data.frame(plot.data2)

for (one in cluster.selected) {
    plot.data2[one, peak.selected %in% peaks.clusters[[one]]] <- 1
}

rownames(plot.data2) <- paste(cluster.selected, cluster.info[cluster.selected, "Patient_Major"], sep = "_")
plot.data <- rbind(plot.data1, plot.data2)

pdf("TF_motif_cluster/Heatmap.Paired.PPARA.PatientPeaks.pdf", 10, 6)
pheatmap(
    plot.data,
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_color = mycolor,
    color = colorRampPalette(c("gray80", "#cd2525"))(2),
    clustering_method = "ward.D2",
    cutree_cols = 12
)
dev.off()

## 7.3. iCMS3 TFs ----
cluster.selected <- cluster.info %>%
    filter(Epi_Group == "Group_1") %>%
    rownames()
marker.paired.list.sub <- marker.paired.list[c(7:8)]

# SOX2
peak.selected <- motif.match[, "SOX2"] %>%
    .[.] %>%
    names() %>%
    intersect(., c(
        unlist(peaks.clusters[cluster.selected]),
        unlist(marker.paired.list.sub)
    ) %>% unique())
length(peak.selected) # 6175 SOX2 peaks,AD + CRC

plot.data1 <- matrix(0, ncol = length(peak.selected), nrow = length(marker.paired.list.sub))
colnames(plot.data1) <- peak.selected
rownames(plot.data1) <- names(marker.paired.list.sub)
plot.data1 <- as.data.frame(plot.data1)

for (one in c("Group_1", "Adenoma")) {
    temp <- paste("COAD34", one, sep = "_")
    plot.data1[temp, colnames(plot.data1) %in% marker.paired.list.sub[[temp]]] <- 1
}

rowSums(plot.data1)
plot.data1 <- plot.data1[c(2, 1), ]

plot.data2 <- matrix(0, ncol = length(peak.selected), nrow = length(cluster.selected))
colnames(plot.data2) <- peak.selected
rownames(plot.data2) <- cluster.selected
plot.data2 <- as.data.frame(plot.data2)

for (one in cluster.selected) {
    plot.data2[one, peak.selected %in% peaks.clusters[[one]]] <- 1
}

rownames(plot.data2) <- paste(cluster.selected, cluster.info[cluster.selected, "Patient_Major"], sep = "_")
plot.data <- rbind(plot.data1, plot.data2)

pdf("TF_motif_cluster/Heatmap.Paired.SOX2.PatientPeaks.pdf", 10, 5)
pheatmap(
    plot.data,
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_color = mycolor,
    color = colorRampPalette(c("gray80", "#cd2525"))(2),
    clustering_method = "ward.D2",
    cutree_cols = 12
)
dev.off()

# FOXA3
peak.selected <- motif.match[, "FOXA3"] %>%
    .[.] %>%
    names() %>%
    intersect(., c(
        unlist(peaks.clusters[cluster.selected]),
        unlist(marker.paired.list.sub)
    ) %>% unique())
length(peak.selected) # 18686 FOXA3 peaks,AD + CRC

plot.data1 <- matrix(0, ncol = length(peak.selected), nrow = length(marker.paired.list.sub))
colnames(plot.data1) <- peak.selected
rownames(plot.data1) <- names(marker.paired.list.sub)
plot.data1 <- as.data.frame(plot.data1)

for (one in c("Group_1", "Adenoma")) {
    temp <- paste("COAD34", one, sep = "_")
    plot.data1[temp, colnames(plot.data1) %in% marker.paired.list.sub[[temp]]] <- 1
}

rowSums(plot.data1)
plot.data1 <- plot.data1[c(2, 1), ]

plot.data2 <- matrix(0, ncol = length(peak.selected), nrow = length(cluster.selected))
colnames(plot.data2) <- peak.selected
rownames(plot.data2) <- cluster.selected
plot.data2 <- as.data.frame(plot.data2)

for (one in cluster.selected) {
    plot.data2[one, peak.selected %in% peaks.clusters[[one]]] <- 1
}

rownames(plot.data2) <- paste(cluster.selected, cluster.info[cluster.selected, "Patient_Major"], sep = "_")
plot.data <- rbind(plot.data1, plot.data2)

pdf("TF_motif_cluster/Heatmap.Paired.FOXA3.PatientPeaks.pdf", 10, 5)
pheatmap(
    plot.data,
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_color = mycolor,
    color = colorRampPalette(c("gray80", "#cd2525"))(2),
    clustering_method = "ward.D2",
    cutree_cols = 12
)
dev.off()
