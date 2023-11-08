setwd("E:/LabWork/Project/CRC_NGS_ATAC/CRC_Epi_scATAC/Results/04.Epi_CIMP")
load("04.Epi_CIMP.RData")
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
    "CIMP_Group" = c(
        "Normal" = "#208a42", "Adenoma" = "#d51f26",
        "CIMP_High" = "#59b891", "CIMP_Low" = "#fa774f", "CIMP_Negative" = "#7a8dbf"
    )
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
write.csv(marker.sig, "marker.sig.csv")
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
write.csv(cluster.info, "cluster.info.csv")
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

## 2.4. compare with molecular subtypes ----
plot.data <- table(proj_Epi$CIMP_Group, proj_Epi$Epi_Group) %>%
    as.data.frame()
plot.data <- plot.data %>%
    filter(Var1 != "Normal" & Var1 != "Adenoma") %>%
    filter(Var2 != "Normal" & Var2 != "Adenoma")

plot.data$Pct <- plot.data$Freq / table(proj_Epi$CIMP_Group)[plot.data$Var1] * 100
plot.data$Pct <- round(plot.data$Pct, 1) %>%
    paste0("%")

pdf("Bar.CMS.vs.CIMP.pdf", 3, 3)
ggplot(plot.data) +
    geom_bar(aes(x = Var1, y = Freq, fill = Var2),
        stat = "identity", position = position_fill(reverse = TRUE)
    ) +
    geom_text(aes(x = Var1, y = Freq, label = Pct),
        size = 3,
        position = position_fill(vjust = 0.5)
    ) +
    xlab("Peak type") +
    ylab("Percent of peaks") +
    scale_fill_manual(values = mycolor$Epi_Group) +
    theme_classic() +
    theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

plot.data <- table(ifelse(proj_Epi$CIMP_Group == "CIMP_Negative", "yes", "no"), proj_Epi$Epi_Group) %>%
    as.data.frame() %>%
        filter(Var1 != "Normal" & Var1 != "Adenoma") %>%
        filter(Var2 != "Normal" & Var2 != "Adenoma")
test.data <- matrix(plot.data$Freq, 2)
colnames(test.data) <- c("Group_1", "Group_2")
rownames(test.data) <- c("no", "yes")
chisq.test(test.data)
fisher.test(test.data)

# 3. find markers for all groups ----
## 3.1. identify diff peaks of each cluster ----
table(proj_Epi$CIMP_Group, proj_Epi$Epi_Group)

markersPeaks.CIMP <- getMarkerFeatures(
    ArchRProj = proj_Epi,
    useMatrix = "PeakMatrix",
    groupBy = "CIMP_Group",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = c("CIMP_High", "CIMP_Low", "CIMP_Negative"),
    bgdGroups = "Normal"
)
saveRDS(markersPeaks.CIMP, "markersPeaks.CIMP.rds")
# markersPeaks.CIMP <- readRDS("markersPeaks.CIMP.rds")

heatmapPeaks <- plotMarkerHeatmap(
    seMarker = markersPeaks.CIMP,
    pal = colorRampPalette(rev(brewer.pal(7, "RdYlBu")))(100),
    cutOff = "FDR <= 0.01 & Log2FC>= 1",
    clusterCols = FALSE,
    labelRows = FALSE,
    limits = c(-1.1, 1.1)
)

pdf("Heatmap.CIMP.marker.peaks1.pdf", 6, 5)
plot(heatmapPeaks)
dev.off()
rm(heatmapPeaks)

# volcano plot
for (one in c("CIMP_High", "CIMP_Low", "CIMP_Negative")) {
    message(paste("plot markers Volcano for ", one, " ...", sep = ""))
    pv <- plotMarkers(
        seMarker = markersPeaks.CIMP,
        name = c(one),
        cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1",
        plotAs = "Volcano"
    )
    pdf(paste("Volcano.markers.", one, ".pdf", sep = ""), 5, 4)
    plot(pv)
    dev.off()
    rm(pv, one)
}

## 3.2. location distribution ----
peakset <- proj_Epi@peakSet
names(peakset) <- paste(
    seqnames(peakset), start(peakset), end(peakset),
    sep = "_"
)
quantile(peakset$distToTSS, probs = seq(0, 1, 0.1))

marker.peak.up.list <- getMarkers(markersPeaks.CIMP,
    cutOff = "FDR <= 0.01 & Log2FC >= 1"
) %>% lapply(function(x) {
    x <- x %>%
        as.data.frame() %>%
        mutate(
            peak_id = paste(seqnames, start, end, sep = "_")
        )
    return(peakset[x$peak_id])
})

marker.peak.down.list <- getMarkers(markersPeaks.CIMP,
    cutOff = "FDR <= 0.01 & Log2FC <= -1"
) %>% lapply(function(x) {
    x <- x %>%
        as.data.frame() %>%
        mutate(
            peak_id = paste(seqnames, start, end, sep = "_")
        )
    return(peakset[x$peak_id])
})

lapply(marker.peak.down.list, length)

for (one in c("CIMP_High", "CIMP_Low", "CIMP_Negative")) {
    plot.data1 <- table(marker.peak.up.list[[one]]$peakType) %>%
        as.data.frame() %>%
        mutate("Group" = "Up")
    plot.data1$Pct <- plot.data1$Freq / sum(plot.data1$Freq) * 100
    plot.data2 <- table(marker.peak.down.list[[one]]$peakType) %>%
        as.data.frame() %>%
        mutate("Group" = "Down")
    plot.data2$Pct <- plot.data2$Freq / sum(plot.data2$Freq) * 100

    plot.data <- rbind(plot.data1, plot.data2)
    plot.data$Pct <- round(plot.data$Pct, 1) %>%
        paste0("%")

    pdf(paste0("Bar.marker.peakType", one, ".pdf"), 5, 2)
    p <- ggplot(plot.data) +
        geom_bar(aes(x = Group, y = Freq, fill = Var1),
            stat = "identity", position = position_fill(reverse = TRUE)
        ) +
        geom_text(aes(x = Group, y = Freq, label = Pct),
            size = 3,
            position = position_fill(vjust = 0.5)
        ) +
        xlab("Peak type") +
        ylab("Percent of peaks") +
        coord_flip() +
        scale_fill_manual(values = c(
            "Distal" = "#59b795", "Exonic" = "#f8784f",
            "Intronic" = "#7a8dbf", "Promoter" = "#df73b7"
        )) +
        theme_classic() +
        theme(axis.title = element_blank())
    plot(p)
    dev.off()
}
rm(one, p, plot.data, plot.data1, plot.data2, peakset)

## 3.3. functional enrichment ----
marker.peak2gene.list <- lapply(
    marker.peak.list,
    function(x) {
        peakset[x$peak_id] %>%
            subset(distToTSS < 5e4) %>%
            .$nearestGene %>%
            unique()
    }
)
lapply(marker.peak2gene.list, length) # 6422 9698 8446

pdf("Upset.markers.group_vs_normal.pdf", 7, 4.5)
upset(
    fromList(marker.peak2gene.list),
    nset = 3,
    order.by = "freq"
)
dev.off()

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
write.csv(GO.marker.list$CIMP_High, "GO.marker.CIMP_High.csv")
write.csv(GO.marker.list$CIMP_Low, "GO.marker.CIMP_Low.csv")
write.csv(GO.marker.list$CIMP_Negative, "GO.marker.CIMP_Negative.csv")

# no good results
rm(GO.marker.list, marker.peak2gene.list)

# 4. motif analysis of each subtype ----
dir.create("homer")
dir.create("bed")

write.table(
    names(marker.peak.up.list$CIMP_High) %>%
        gsub("_", "\t", .) %>%
        as.data.frame() %>%
        mutate(
            "Name" = names(marker.peak.up.list$CIMP_High),
            "length" = 500,
            "strand" = "."
        ),
    "bed/marker.peak.tumor.CIMP_High.bed",
    col.names = FALSE,
    sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
    names(marker.peak.up.list$CIMP_Low) %>%
        gsub("_", "\t", .) %>%
        as.data.frame() %>%
        mutate(
            "Name" = names(marker.peak.up.list$CIMP_Low),
            "length" = 500,
            "strand" = "."
        ),
    "bed/marker.peak.tumor.CIMP_Low.bed",
    col.names = FALSE,
    sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
    names(marker.peak.up.list$CIMP_Negative) %>%
        gsub("_", "\t", .) %>%
        as.data.frame() %>%
        mutate(
            "Name" = names(marker.peak.up.list$CIMP_Negative),
            "length" = 500,
            "strand" = "."
        ),
    "bed/marker.peak.tumor.CIMP_Negative.bed",
    col.names = FALSE,
    sep = "\t", quote = FALSE, row.names = FALSE
)

## 4.1. run HOMER ----
# nohup findMotifsGenome.pl bed/marker.peak.tumor.CIMP_High.bed hg38 homer/CIMP_High -size 200
# nohup findMotifsGenome.pl bed/marker.peak.tumor.CIMP_Low.bed hg38 homer/CIMP_Low -size 200
# nohup findMotifsGenome.pl bed/marker.peak.tumor.CIMP_Negative.bed hg38 homer/CIMP_Negative -size 200

homer.res <- list(
    "CIMP_High" = homer.parser("homer/CIMP_High/knownResults.txt") %>%
        mutate(Group = "CIMP_High"),
    "CIMP_Low" = homer.parser("homer/CIMP_Low/knownResults.txt") %>%
        mutate(Group = "CIMP_Low"),
    "CIMP_Negative" = homer.parser("homer/CIMP_Negative/knownResults.txt") %>%
        mutate(Group = "CIMP_Negative")
)

sapply(homer.res, dim) # 419 TFs
identical(homer.res$CIMP_High$TF, homer.res$CIMP_Low$TF)

## 4.2. identify significant TFs ----
dir.create("TF_motif")
pdf("TF_motif/Dot.motif.CIMP_High.pdf", 5, 4)
ggplot(homer.res$CIMP_High, aes(x = Log2_Enrichment, y = log.p.value)) +
    geom_point(aes(color = Diff), size = 1) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    geom_hline(yintercept = 50, linetype = "dashed") +
    scale_color_manual(values = c("none" = "grey", "up" = "red")) +
    ggrepel::geom_text_repel(aes(label = TF), size = 2, max.overlaps = 30) +
    theme_classic()
dev.off()
pdf("TF_motif/Dot.motif.CIMP_Low.pdf", 5, 4)
ggplot(homer.res$CIMP_Low, aes(x = Log2_Enrichment, y = log.p.value)) +
    geom_point(aes(color = Diff), size = 1) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    geom_hline(yintercept = 50, linetype = "dashed") +
    scale_color_manual(values = c("none" = "grey", "up" = "red")) +
    ggrepel::geom_text_repel(aes(label = TF), size = 2, max.overlaps = 30) +
    theme_classic()
dev.off()
pdf("TF_motif/Dot.motif.CIMP_Negative.pdf", 5, 4)
ggplot(homer.res$CIMP_Negative, aes(x = Log2_Enrichment, y = log.p.value)) +
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
    row.names = homer.res$CIMP_High$TF,
    "CIMP_High" = homer.res$CIMP_High$Log2_Enrichment,
    "CIMP_Low" = homer.res$CIMP_Low$Log2_Enrichment,
    "CIMP_Negative" = homer.res$CIMP_Negative$Log2_Enrichment
)

pdf("TF_motif/Heatmap.motif.sig.pdf", 10, 4)
p <- pheatmap::pheatmap(t(FC.mat[TF.sig, ]),
    color = colorRampPalette(ArchRPalettes$comet)(100),
    scale = "column",
    cluster_rows = FALSE,
    clustering_method = "ward.D2",
    cutree_cols = 4
)
dev.off()

## 4.3. visulize CIMP_High specific ----
# dot plot
TF.selected <- cutree(p$tree_col, 4) %>%
    .[. %in% c(2)] %>%
    names()

plot.data <- lapply(homer.res, function(x) {
    x <- x[x$TF %in% TF.selected, ]
    return(x)
}) %>% do.call(rbind, .)

plot.data$Group <- factor(plot.data$Group, levels = rev(c("CIMP_High", "CIMP_Low", "CIMP_Negative")))
plot.data$TF <- factor(plot.data$TF, levels = rev(TF.selected))
plot.data[plot.data$log.p.value > 1000, ]$log.p.value <- 1000

pdf("TF_motif/Dot.motif.sig.CIMP_High1.pdf", 6, 2)
ggplot(plot.data) +
    geom_point(aes(
        x = Group, y = TF,
        fill = log.p.value, size = Log2_Enrichment
    ), pch = 21) +
    scale_fill_viridis_c() +
        theme_bw() +
        coord_flip() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

rm(FC.mat, plot.data, p)

# umap
TF.CISBP <- getFeatures(proj_Epi, useMatrix = "MotifMatrix")
grep("TEAD", TF.CISBP, value = TRUE)
# align CISBP and HOMER
TF.sig.align <- sapply(TF.sig, function(x) {
    res <- paste("z:", x, sep = "")
    res <- grep(res, TF.CISBP, value = TRUE)
    return(res)
}) %>%
    unlist() %>%
    sort()
names(TF.sig.align) <- TF.sig.align %>%
    gsub("z:", "", .) %>%
    gsub("_.+?$", "", .)

TF.selected <- TF.sig.align[c(TF.selected, "TEAD1")]
TF.selected <- TF.selected[!is.na(TF.selected)]
TF.selected <- TF.selected[c("LEF1", "TCF3", "TCF7", "TEAD1")]
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
pdf("TF_motif/UMAP.TF.sig.CIMP_High1.pdf", 15, 8)
wrap_plots(c(p1, p2), ncol = 4)
dev.off()
rm(p1, p2, p)

# footprints
motifPositions <- getPositions(proj_Epi)

proj_Epi <- addGroupCoverages(
    ArchRProj = proj_Epi,
    groupBy = "CIMP_Group"
)

seFoot <- getFootprints(
    ArchRProj = proj_Epi,
    positions = motifPositions[gsub("z:", "", TF.selected)],
    groupBy = "CIMP_Group",
    useGroups = c("CIMP_High", "CIMP_Low", "CIMP_Negative"),
)

plotFootprints(
    seFoot = seFoot,
    ArchRProj = proj_Epi,
    normMethod = "Subtract",
    addDOC = FALSE,
    plotName = "Footprints_CIMP_Group_Subtract.sw10.pdf",
    pal = c(mycolor$CIMP_Group, "Normal" = "#208a42"),
    smoothWindow = 10
)

## 4.4. positive feedback loop ----
write.table(
    names(marker.peak.down.list$CIMP_High) %>%
        gsub("_", "\t", .) %>%
        as.data.frame() %>%
        mutate(
            "Name" = names(marker.peak.down.list$CIMP_High),
            "length" = 500,
            "strand" = "."
        ),
    "bed/marker.peak.tumor.CIMP_High.down.bed",
    col.names = FALSE,
    sep = "\t", quote = FALSE, row.names = FALSE
)
# findMotifsGenome.pl bed/marker.peak.tumor.CIMP_High.down.bed hg38 homer/CIMP_High_down -size 200
CIMP.H.down.homer <- homer.parser(
    "homer/CIMP_High_down/knownResults.txt",
    log.p.value = 50, log2.enrichment = 1
)
CIMP.H.down.homer <- CIMP.H.down.homer[CIMP.H.down.homer$TF != "UNKNOWN", ]

CIMP.H.down.homer[CIMP.H.down.homer$Diff != "up", ]$TF <- NA
CIMP.H.down.homer[CIMP.H.down.homer$Diff != "up", ]$Family <- NA

pdf("TF_motif/Dot.motif.CIMP_High.down.pdf", 5, 4)
ggplot(CIMP.H.down.homer, aes(x = Log2_Enrichment, y = log.p.value)) +
    geom_point(aes(color = Family), size = 1) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    geom_hline(yintercept = 50, linetype = "dashed") +
    scale_color_manual(values = brewer.pal(9, "Set3")) +
    ggrepel::geom_text_repel(aes(label = TF), size = 2, max.overlaps = 30) +
    theme_classic()
dev.off()

CIMP.H.down.motif <- CIMP.H.down.homer %>%
    filter(Diff == "up") %>%
    pull(TF)
CIMP.H.down.motif <- CIMP.H.down.motif[!CIMP.H.down.motif == "UNKNOWN"]

CIMP.H.down.gene <- marker.peak.down.list$CIMP_High %>%
    subset(peakType == "Promoter") %>%
    # subset(peakType != "Distal") %>%
    .$nearestGene %>%
    unique()

intersect(CIMP.H.down.motif, CIMP.H.down.gene) # "HOXB13" "KLF3"

v <- Venn(list(
    "PLS" = CIMP.H.down.gene,
    "TF" = CIMP.H.down.motif
))

pdf("TF_motif/Venn.PLS.TF.pdf", 4, 4)
plot(v, doWeights = TRUE, show = list(Faces = FALSE))
dev.off()

rm(v, motifPositions, TF.CISBP)

proj_Epi <- saveArchRProject(ArchRProj = proj_Epi, load = TRUE)
saveRDS(cluster.info, "cluster.info.rds")
save.image("04.Epi_CIMP.RData")
