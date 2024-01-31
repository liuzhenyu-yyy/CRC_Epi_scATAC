setwd("E:/LabWork/Project/CRC_NGS_ATAC/CRC_Epi_scATAC/Results/06.Epi_Intratumor")
load("Epi_intratumor.RData")

source("../../code/00.Requirements.R")
rm(foo)

proj_Epi <- loadArchRProject(project.dir.epi, force = TRUE)
sample.info.epi <- proj_Epi@cellColData %>% as.data.frame()

mycolor <- list(
    "Gender" = c("Female" = "#cab2d6", "Male" = "#6a3d9a"),
    "MSI" = c("MSS" = "#ffff99", "MSI-H" = "#b15928", "NA" = "gray70"),
    "Side" = c("Right" = "#b2df8a", "Left" = "#33a02c"),
    "Group" = c("Group_1" = "#62b7e6", "Group_2" = "#283891"),
    "Epi_Group" = c(
        "Normal" = "#208a42", "Adenoma" = "#d51f26",
        "Group_1" = "#62b7e6", "Group_2" = "#283891"
    ),
    "iCMS" = c(
        "Normal" = "#208a42", "Adenoma" = "#d51f26",
        "iCMS2" = "#283891", "iCMS3" = "#62b7e6 "
    ),
    "CIMP_Group" = c(
        "Normal" = "gray50", "Adenoma" = "gray51",
        "CIMP_High" = "#f76960", "CIMP_Low" = "#fbe625", "CIMP_Negative" = "#89dc0e"
    ),
    "Clusters" = ArchR::paletteDiscrete(proj_Epi$Clusters),
    "tree" = ArchR::paletteDiscrete(factor(1:8), set = "bear")
)

# 1. identify ITH of COAD26----
## 1.1. CNV of specific patient ----
dir.create("CNV_raw")
dir.create("CNV_qc")
load("../01.All_scCNV/CRC_CNV.rda")

patient.selected <- c("COAD12", "COAD17", "COAD24", "COAD30", "COAD32", "COAD33", "COAD35")
scales::show_col(mycolor$tree)

close.dev <- function() {
    while (dev.cur() != 1) {
        dev.off()
    }
}
for (one in patient.selected) {
    cell.selected <- sample.info.epi %>%
        filter(Epi_type %in% c("Malignant")) %>%
        filter(Patient == one) %>%
        rownames()

    anno.row <- sample.info.epi[cell.selected, ] %>%
        select(c("Clusters"))

    close.dev()
    png(paste0("CNV_raw/Heatmap.CNV.ITH.", one, ".raw.png"), 1000, 750)
    tree.all <- pheatmap::pheatmap(CNV.FC[cell.selected, ],
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
        cutree_rows = 5,
        gaps_col = which(!duplicated(anno.col$seqnames)) - 1
    )
    dev.off()

    temp <- cutree(tree.all$tree_row, 5)
    anno.row$tree <- factor(temp[rownames(anno.row)], levels = 1:5)

    close.dev()
    png(paste0("CNV_raw/Heatmap.CNV.ITH.", one, ".clustered.png"), 1000, 750)
    pheatmap::pheatmap(CNV.FC[cell.selected, ],
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
        cutree_rows = 5,
        gaps_col = which(!duplicated(anno.col$seqnames)) - 1
    )
    dev.off()
}

cluster.remove <- list(
    "COAD12" = c(1, 5),
    "COAD17" = c(3:5),
    "COAD24" = c(3, 5),
    "COAD30" = c(5),
    "COAD32" = c(1),
    "COAD33" = c(4, 5),
    "COAD35" = c(1:3)
)
ncluster <- c(
    "COAD12" = 2,
    "COAD17" = 2,
    "COAD24" = 2,
    "COAD30" = 4,
    "COAD32" = 5,
    "COAD33" = 4,
    "COAD35" = 3
)
dir.create("CNV_qc")
subclone.list <- list()
for (one in patient.selected) {
    cell.selected <- sample.info.epi %>%
        filter(Epi_type %in% c("Malignant")) %>%
        filter(Patient == one) %>%
        rownames()

    anno.row <- sample.info.epi[cell.selected, ] %>%
        select(c("Clusters"))

    close.dev()
    png(paste0("CNV_qc/Heatmap.CNV.ITH.", one, ".raw.png"), 1000, 750)
    tree.all <- pheatmap::pheatmap(CNV.FC[cell.selected, ],
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
        cutree_rows = 5,
        gaps_col = which(!duplicated(anno.col$seqnames)) - 1
    )
    dev.off()

    temp <- cutree(tree.all$tree_row, 5)
    anno.row$tree <- factor(temp[rownames(anno.row)], levels = 1:5)
    cell.selected <- cell.selected[!anno.row[cell.selected, ]$tree %in% cluster.remove[[one]]]
    anno.row <- anno.row[rownames(anno.row) %in% cell.selected, ]

    close.dev()
    png(paste0("CNV_qc/Heatmap.CNV.ITH.", one, ".qc.png"), 1000, 750)
    tree.all <- pheatmap::pheatmap(CNV.FC[cell.selected, ],
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
        cutree_rows = ncluster[one],
        gaps_col = which(!duplicated(anno.col$seqnames)) - 1
    )
    dev.off()

    temp <- cutree(tree.all$tree_row, ncluster[one])
    anno.row$tree <- factor(temp[rownames(anno.row)], levels = 1:ncluster[one])

    subclone.list[[one]] <- split(rownames(anno.row), anno.row$tree)
    names(subclone.list[[one]]) <- paste("subclone_", 1:ncluster[one], sep = "")

    close.dev()
    png(paste0("CNV_qc/Heatmap.CNV.ITH.", one, ".qc.png"), 1000, 750)
    pheatmap::pheatmap(CNV.FC[cell.selected, ],
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
        cutree_rows = ncluster[one],
        gaps_col = which(!duplicated(anno.col$seqnames)) - 1
    )
    dev.off()
}
saveRDS(subclone.list, "subclone.list.rds")

## 1.2. cross patient ----
plot.data <- sample.info.epi %>%
    filter(Epi_type == "Malignant") %>%
    filter(Patient %in% patient.selected)

plot.data <- table(plot.data$Clusters, plot.data$Patient) %>%
    as.data.frame()
colnames(plot.data) <- c("Clusters", "Patient", "Freq")
plot.data$iCMS <- ifelse(plot.data$Patient %in% c("COAD30", "COAD35"), "iCMS3", "iCMS2")

pdf("Bar.Pct.Patient.Cluster.pdf", 4, 3)
ggplot(plot.data, aes(x = Patient, y = Freq, fill = Clusters)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(x = "Patient", y = "Number of cells", fill = "Patient") +
    facet_grid(cols = vars(iCMS), scales = "free", space = "free") +
    scale_fill_manual(values = ArchR::paletteDiscrete(proj_Epi$Clusters))
dev.off()

plot.data <- as.data.frame(ncluster) %>%
    mutate(Patient = rownames(.)) %>%
    mutate(iCMS = ifelse(Patient %in% c("COAD30", "COAD35"), "iCMS3", "iCMS2"))
pdf("Bar.Pct.Patient.nSubclone.pdf", 4, 3)
ggplot(plot.data, aes(x = Patient, y = ncluster, fill = Patient)) +
    geom_bar(stat = "identity") +
    facet_grid(cols = vars(iCMS), scales = "free", space = "free") +
    theme_classic() +
    ylab("No. of subclones") +
    scale_fill_manual(values = ArchR::paletteDiscrete(proj_Epi$Patient))
dev.off()

rm(plot.data)

events.list <- list(
    "COAD12" = c("9p", "9q"),
    "COAD17" = c("10p", "10q"),
    "COAD24" = c("10p", "10q"),
    "COAD30" = c("2p", "7p", "7q", "18", "20", "21"),
    "COAD32" = c(),
    "COAD33" = c("3q", "4", "5", "7"),
    "COAD35" = c("7", "20", "21")
)

# 2. copy-number phylogenies ----
## 2.1. all cells ----
load("../01.All_scCNV/CRC_CNV.rda")
rm(sample.info, anno.row)

cell.selected <- sample.info.epi %>%
    filter(Clusters %in% c("C17", "C22")) %>%
    filter(Patient == "COAD32") %>%
    rownames()

anno.row <- sample.info.epi[cell.selected, ] %>%
    select(c("Clusters"))
CNV.FC <- CNV.FC[cell.selected, ]

dev.off()
png("Heatmap.CNV.ITH.all.png", 1000, 750)
tree.all <- pheatmap::pheatmap(CNV.FC,
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
    cutree_rows = 5,
    gaps_col = which(!duplicated(anno.col$seqnames)) - 1
)
dev.off()

temp <- cutree(tree.all$tree_row, 5)
anno.row$tree <- factor(temp[rownames(anno.row)], levels = 1:5)

pdf("legents.pdf", 6, 3)
scales::show_col(mycolor$tree)
dev.off()

png("Heatmap.CNV.ITH.all.png", 1000, 750)
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
    cutree_rows = 5,
    gaps_col = which(!duplicated(anno.col$seqnames)) - 1
)
dev.off()

## 2.2. QC cells ----
cell.selected <- names(temp)[temp != 1]
anno.row <- anno.row[cell.selected, ]
CNV.FC <- CNV.FC[cell.selected, ]

png("Heatmap.CNV.ITH.QC.png", 1000, 750)
tree.all <- pheatmap::pheatmap(CNV.FC,
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
    cutree_rows = 5,
    gaps_col = which(!duplicated(anno.col$seqnames)) - 1
)
dev.off()

temp <- cutree(tree.all$tree_row, 5)
anno.row$tree <- factor(temp[rownames(anno.row)], levels = 1:5)
png("Heatmap.CNV.ITH.QC.png", 1000, 750)
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
    cutree_rows = 5,
    gaps_col = which(!duplicated(anno.col$seqnames)) - 1
)
dev.off()

anno.row$subclone <- c("subclone_1", "subclone_3", "subclone_4", "subclone_5", "subclone_2")[anno.row$tree]

sample.info.ITH <- sample.info.epi[cell.selected, ] %>%
    mutate(Subclone = anno.row[cell.selected, ]$subclone)
table(sample.info.ITH$Subclone, sample.info.ITH$Clusters)

sample.info.ITH <- sample.info.ITH %>%
    filter(!(Subclone == "subclone_1" & Clusters == "C22")) %>%
    filter(!(Subclone == "subclone_5" & Clusters == "C17"))

mycolor$Subclone <- mycolor$tree[c(1, 5, 2, 3, 4)]
names(mycolor$Subclone) <- paste("subclone_", 1:5, sep = "")
mycolor$tree <- NULL

pdf("legents.subclone.pdf", 6, 3)
scales::show_col(mycolor$Subclone)
dev.off()

## 2.3. difference on iCMS modules ----
wilcox.test(
    sample.info.ITH[sample.info.ITH$Clusters == "C17", ]$Module.iCMS2,
    sample.info.ITH[sample.info.ITH$Clusters == "C22", ]$Module.iCMS2
)$p.value

pdf("Violin.Module.iCMS2.pdf", 5, 3)
ggplot(sample.info.ITH, aes(x = Module.iCMS2, y = factor(Subclone,
    levels = rev(c("subclone_2", "subclone_1", "subclone_3", "subclone_5", "subclone_4"))
))) +
    geom_violin(aes(fill = Subclone)) +
    theme_bw() +
    # ggpubr::stat_compare_means(
    #     method = "wilcox.test",
    #     comparisons = list(
    #         c("subclone_1", "subclone_4"), c("subclone_1", "subclone_5"),
    #         c("subclone_2", "subclone_4"), c("subclone_2", "subclone_5"),
    #          c("subclone_3", "subclone_4"), c("subclone_3", "subclone_5")
    #     ), label = "p.signif"
    # ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    facet_grid(rows = vars(Clusters), scales = "free", space = "free") +
    labs(x = "Module.iCMS2", y = "Subclone") +
    ggtitle("Inter-cluster  p-value < 2.2e-16") +
    scale_fill_manual(values = mycolor$Subclone)
dev.off()

## 2.4 subclone-by-arm copy number matrix ----
# sliding window
sw <- readRDS("../01.All_scCNV/sw.rds")
names(sw) <- sw$name # 1190
sw <- sw[colnames(CNV.FC)] # 1067

chr.anno <- fread("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz",
    col.names = c("chrom", "chromStart", "chromEnd", "name", "gieStain")
)

# hg38 chromosome arms
chr.anno <- chr.anno[chr.anno$chrom %in% unique(seqnames(sw)), ]
colnames(chr.anno) <- c("seqnames", "start", "end", "band", "gieStain")
chr.anno <- GenomicRanges::GRanges(chr.anno)

chr.anno$arm <- chr.anno$band
chr.anno$arm <- sub("p.*", "p", chr.anno$arm)
chr.anno$arm <- sub("q.*", "q", chr.anno$arm)
chr.anno$arm <- paste(seqnames(chr.anno), chr.anno$arm, sep = "")
table(chr.anno$arm) %>% length()
start(chr.anno) <- start(chr.anno) + 1
end(chr.anno) <- end(chr.anno) + 1

# align
sw$arm <- "none"
for (one in names(sw)) {
    temp <- sw[one]
    temp <- countOverlaps(chr.anno, temp) > 0
    temp <- chr.anno$arm[temp]
    temp <- table(temp) %>%
        sort(decreasing = TRUE) %>%
        names() %>%
        head(1)
    sw[one]$arm <- temp
}
table(sw$arm)
setdiff(paste("chr", rep(1:22, each = 2), c("p", "q"), sep = ""), unique(sw$arm))

CNV.FC <- CNV.FC[rownames(sample.info.ITH), ]

CNV.subclone.sw <- CNV.FC %>%
    group_by(sample.info.ITH$Subclone) %>%
    summarise(across(everything(), list(mean))) %>%
    as.data.frame()
CNV.subclone.sw$`sample.info.ITH$Subclone` <- NULL
rownames(CNV.subclone.sw) <- paste("subclone_", 1:5, sep = "")
colnames(CNV.subclone.sw) <- gsub("_1", "", colnames(CNV.subclone.sw))

identical(colnames(CNV.subclone.sw), names(sw))
CNV.subclone.arm <- CNV.subclone.sw %>%
    t() %>%
    as.data.frame() %>%
    group_by(sw$arm) %>%
    summarise(across(everything(), list(mean))) %>%
    as.data.frame()
rownames(CNV.subclone.arm) <- CNV.subclone.arm$`sw$arm`
CNV.subclone.arm$`sw$arm` <- NULL
CNV.subclone.arm <- t(CNV.subclone.arm) %>% as.data.frame()
rownames(CNV.subclone.arm) <- paste("subclone_", 1:5, sep = "")

temp <- paste("chr", rep(1:22, each = 2), c("p", "q"), sep = "")
temp <- temp[temp %in% colnames(CNV.subclone.arm)]
CNV.subclone.arm <- CNV.subclone.arm[, temp]

CNV.subclone <- matrix(0, nrow = 5, ncol = 44)
rownames(CNV.subclone) <- paste("subclone_", 1:5, sep = "")
colnames(CNV.subclone) <- paste("chr", rep(1:22, each = 2), c("p", "q"), sep = "")
CNV.subclone <- as.data.frame(CNV.subclone)
CNV.subclone[rownames(CNV.subclone.arm), colnames(CNV.subclone.arm)] <- CNV.subclone.arm

pheatmap::pheatmap(CNV.subclone,
    color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
    border_color = FALSE,
    cluster_rows = FALSE, cluster_cols = FALSE
)

CN.subclone <- CNV.subclone
CN.subclone[CN.subclone > 0.25] <- 100
CN.subclone[CN.subclone < (-0.25)] <- (-100)

pdf("Heatmap.CNV.subclone.pdf", 8, 3)
pheatmap::pheatmap(CN.subclone[5:1, ],
    color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
    border_color = FALSE,
    cluster_rows = FALSE, cluster_cols = FALSE,
    gaps_col = (1:10) * 2
)
dev.off()

CN.subclone.sw <- CNV.subclone.sw
CN.subclone.sw[] <- 2

for (one in colnames(CN.subclone)) {
    for (sb in rownames(CN.subclone)) {
        sws <- sw[sw$arm == one] %>% names()
        if (length(sws) == 0) {
            next
        }
        if (CN.subclone[sb, one] == 100) {
            CN.subclone.sw[sb, sws] <- 3
        } else if (CN.subclone[sb, one] == (-100)) {
            CN.subclone.sw[sb, sws] <- 1
        }
    }
}

pdf("Heatmap.CNV.subclone.sw.pdf", 8, 3)
pheatmap::pheatmap(CN.subclone.sw[5:1, ],
    color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
    border_color = FALSE,
    cluster_rows = FALSE, cluster_cols = FALSE
)
dev.off()

out.data <- CN.subclone %>%
    mutate(Subclone = rownames(CN.subclone)) %>%
    reshape2::melt(., id.vars = "Subclone")
out.data$CN <- out.data$value
out.data$CN[out.data$value == 100] <- 3
out.data$CN[out.data$value == (-100)] <- 1
out.data$CN[out.data$value != 100 & out.data$value != (-100)] <- 2

out.data$chrom <- gsub("p|q", "", out.data$variable)
out.data$start <- 0
out.data$start[grep("q", out.data$variable)] <- 1
out.data$end <- 1
out.data$end[grep("q", out.data$variable)] <- 2

out.data <- out.data[, c(1, 5:7, 4)]
colnames(out.data) <- c("sample_id", "chrom", "start", "end", "cn_total")
write.table(out.data, "CNV.subclone.arm.txt", sep = "\t", quote = FALSE, row.names = FALSE)

rm(CNV.FC, anno.col, anno.row, anno.color, cell.selected, temp, tree.all)
rm(CNV.subclone, CNV.subclone.arm, CNV.subclone.sw, CN.subclone, out.data, CN.subclone.sw)
rm(sw, sb, sws, one, chr.anno, plot.data, temp)

# run medicc2
# medicc2 CNV.subclone.arm.txt ./MEDICC2 --total-copy-numbers --input-allele-columns cn_total

# 3. differential analysis ----
## 3.1. identify sobclone peaks  ----
sample.info.epi$is_COAD32 <- "None"
sample.info.epi[sample.info.epi$Patient == "COAD32", ]$is_COAD32 <- "COAD32"
proj_Epi$is_COAD32 <- sample.info.epi$is_COAD32

p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    pal = c("None" = "gray75", "COAD32" = "#791618"),
    name = "is_COAD32", embedding = "UMAP",
    size = 0.2, plotAs = "points",
    labelMeans = FALSE
)
pdf("UAMP.Epi.is_COAD32.pdf", 5, 5)
plot(p)
dev.off()

sample.info.epi$Subclone <- "None"
sample.info.epi[sample.info.epi$Epi_type == "Normal", ]$Subclone <- "Normal"
sample.info.epi[rownames(sample.info.ITH), ]$Subclone <- sample.info.ITH$Subclone
proj_Epi$Subclone <- sample.info.epi$Subclone

p <- plotEmbedding(
    ArchRProj = proj_Epi, colorBy = "cellColData",
    pal = c(mycolor$Subclone, "None" = "gray75", "Normal" = "gray75"),
    name = "Subclone", embedding = "UMAP",
    size = 0.2, plotAs = "points",
    labelMeans = FALSE
)
pdf("UAMP.Epi.Subclone.pdf", 5, 5)
plot(p)
dev.off()

marker.subclone <- getMarkerFeatures(
    ArchRProj = proj_Epi,
    useMatrix = "PeakMatrix",
    groupBy = "Subclone",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = c("subclone_2", "subclone_1", "subclone_3", "subclone_5", "subclone_4"),
    bgdGroups = c("Normal")
)

saveRDS(marker.subclone, "markersPeaks.subclone.rds")
marker.subclone <- readRDS("markersPeaks.subclone.rds")
# markersPeaks.CIMP <- readRDS("markersPeaks.CIMP.rds")

# volcano plot
for (one in names(marker.peak.list)) {
    message(paste("plot markers Volcano for ", one, " ...", sep = ""))
    pv <- plotMarkers(
        seMarker = marker.subclone,
        name = c(one),
        cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1",
        plotAs = "Volcano"
    )
    pdf(paste("Volcano.markers.", one, ".pdf", sep = ""), 5, 4)
    plot(pv)
    dev.off()
    rm(pv, one)
}

heatmapPeaks <- plotMarkerHeatmap(
    seMarker = marker.subclone,
    pal = colorRampPalette(rev(brewer.pal(7, "RdYlBu")))(100),
    cutOff = "FDR <= 0.01 & Log2FC>= 1",
    clusterCols = FALSE,
    labelRows = FALSE,
    limits = c(-1.1, 1.1)
)

pdf("Heatmap.marker.peaks.subclone.pdf", 6, 5)
plot(heatmapPeaks)
dev.off()
rm(heatmapPeaks)

marker.peak.list <- getMarkers(marker.subclone,
    cutOff = "FDR <= 0.01 & Log2FC >= 1"
) %>% lapply(function(x) {
    x <- x %>%
        as.data.frame() %>%
        mutate(
            peak_id = paste(seqnames, start, end, sep = "_")
        )
    return(x$peak_id)
})

dir.create("bed")
for (one in names(marker.peak.list)) {
    write.table(
        marker.peak.list[[one]] %>%
            gsub("_", "\t", .) %>%
            as.data.frame() %>%
            mutate(
                "Name" = marker.peak.list[[one]],
                "length" = 500,
                "strand" = "."
            ),
        paste0("bed/marker.peak.", one, ".bed"),
        col.names = FALSE,
        sep = "\t", quote = FALSE, row.names = FALSE
    )
}

## 3.2. run Homer, identify significant TFs ----
# nohup findMotifsGenome.pl bed/marker.peak.subclone_1.bed hg38 homer/subclone_1 -size 200 > homer/subclone1.txt 2>&1 &
homer.res <- list(
    "subclone_1" = homer.parser("homer/subclone_1/knownResults.txt") %>%
        mutate(Group = "subclone_1"),
    "subclone_3" = homer.parser("homer/subclone_3/knownResults.txt") %>%
        mutate(Group = "subclone_3"),
    "subclone_4" = homer.parser("homer/subclone_4/knownResults.txt") %>%
        mutate(Group = "subclone_4"),
    "subclone_5" = homer.parser("homer/subclone_5/knownResults.txt") %>%
        mutate(Group = "subclone_5")
)

sapply(homer.res, dim) # 419 TFs
identical(homer.res$subclone_1$TF, homer.res$subclone_3$TF)

dir.create("TF_motif")
for (one in names(homer.res)) {
    p <- ggplot(homer.res[[one]], aes(x = Log2_Enrichment, y = log.p.value)) +
        geom_point(aes(color = Diff), size = 1) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        geom_hline(yintercept = 50, linetype = "dashed") +
        scale_color_manual(values = c("none" = "grey", "up" = "red")) +
        ggrepel::geom_text_repel(aes(label = TF), size = 2, max.overlaps = 30) +
        theme_classic()
    pdf(paste0("TF_motif/Dot.motif.", one, ".pdf"), 5, 4)
    plot(p)
    dev.off()
}
rm(p)

## 3.3 compare TFs ----
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
    row.names = homer.res$subclone_1$TF,
    "subclone_1" = homer.res$subclone_1$Log2_Enrichment,
    "subclone_3" = homer.res$subclone_3$Log2_Enrichment,
    "subclone_4" = homer.res$subclone_4$Log2_Enrichment,
    "subclone_5" = homer.res$subclone_5$Log2_Enrichment
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

## 3.4. visulize subclone specific ----
# dot plot
TF.selected <- c("FOXA3", "FOXA2", "CDX2", "LEF1", "TCF3", "TR4", "PPARA", "HNF4A")

plot.data <- lapply(homer.res, function(x) {
    x <- x[x$TF %in% TF.selected, ]
    return(x)
}) %>% do.call(rbind, .)

plot.data$Group <- factor(plot.data$Group, levels = c("subclone_4", "subclone_5", "subclone_3", "subclone_1"))
plot.data$TF <- factor(plot.data$TF, levels = rev(TF.selected))
plot.data[plot.data$log.p.value > 1000, ]$log.p.value <- 1000

plot.data$Clusters <- ifelse(plot.data$Group %in% c("subclone_4", "subclone_5"), "C22", "C17")
plot.data$Clusters <- factor(plot.data$Clusters, levels = c("C22", "C17"))

pdf("TF_motif/Dot.motif.sig.subclone.pdf", 3.5, 3.6)
ggplot(plot.data) +
    geom_point(aes(
        x = Group, y = TF,
        fill = log.p.value, size = Log2_Enrichment
    ), pch = 21) +
    scale_fill_viridis_c() +
        theme_bw() +
        facet_wrap(~Clusters, scales = "free_x") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# violin plot
TF.CISBP <- getFeatures(proj_Epi, useMatrix = "MotifMatrix") %>%
    grep("z:", ., value = TRUE)

TF.selected.name <- sapply(TF.selected, function(x) {
    grep(x, TF.CISBP, value = TRUE)
}) %>%
    do.call(c, .)
TF.selected.name <- TF.selected.name[c(1, 2, 5, 4, 6, 7)]
p <- plotGroups(
    ArchRProj = proj_Epi,
    groupBy = "Subclone",
    colorBy = "MotifMatrix",
    name = TF.selected.name,
    plotAs = "violin",
)

plot.data <- data.frame()

for (one in names(p)) {
    temp <- p[[one]]$data %>%
        filter(!x %in% c("None", "Normal"))
    colnames(temp) <- c("Subclone", "Deviation")
    temp$TF <- one
    plot.data <- rbind(plot.data, temp)
}

plot.data$Subclone <- factor(plot.data$Subclone,
    levels = c("subclone_4", "subclone_5", "subclone_3", "subclone_1", "subclone_2")
)
plot.data$TF <- factor(plot.data$TF, levels = (TF.selected.name))

pdf("TF_motif/Violin.motif.sig.subclone.pdf", 8, 3.5)
ggplot(plot.data, aes(x = Subclone, y = Deviation)) +
    geom_violin(aes(fill = Subclone), show.legend = FALSE) +
    theme_classic() +
    ggpubr::stat_compare_means(
        method = "wilcox.test",
        comparisons = list(
            c("subclone_1", "subclone_3"), c("subclone_3", "subclone_5"),
             c("subclone_4", "subclone_5")
        ), size = 2, label = "p.signif"
    ) +
    facet_wrap(~TF, scales = "free") +
    labs(x = "Subclone", y = "Deviation") +
    scale_fill_manual(values = mycolor$Subclone) +
    coord_flip()
dev.off()

rm(plot.data, p, p.list, one)

# footprints
motifPositions <- getPositions(proj_Epi)

proj_Epi <- addGroupCoverages(
    ArchRProj = proj_Epi,
    groupBy = "Subclone"
)

seFoot <- getFootprints(
    ArchRProj = proj_Epi,
    positions = motifPositions[gsub("z:", "", TF.selected.name)],
    groupBy = "Subclone",
    useGroups = c("subclone_1", "subclone_3", "subclone_4", "subclone_5")
)

plotFootprints(
    seFoot = seFoot,
    ArchRProj = proj_Epi,
    normMethod = "Subtract",
    addDOC = FALSE,
    plotName = "Footprints_subclone_Subtract.sw10.pdf",
    pal = c(mycolor$Subclone),
    smoothWindow = 10
)

# 4. Validate in other patients ----
dir.create("All_patient")
subclone.list <- readRDS("subclone.list.rds")
sample.info.epi$Subclone_all <- "None"

for (patient in names(subclone.list)) {
    for (subclone in names(subclone.list[[patient]])) {
        sample.info.epi[rownames(sample.info.epi) %in%
            subclone.list[[patient]][[subclone]], ]$Subclone_all <- paste(patient, subclone, sep = "_")
    }
}
sample.info.epi[sample.info.epi$Epi_type == "Normal", ]$Subclone_all <- "Normal"
table(sample.info.epi$Subclone_all)

proj_Epi$Subclone_all <- sample.info.epi$Subclone_all
colnames(sample.info.epi)

## 4.1. iCMS modules and TFs ----
# iCMS modules
p <- plotGroups(
    ArchRProj = proj_Epi,
    groupBy = "Subclone_all",
    colorBy = "cellColdata",
    name = "Module.iCMS2",
    plotAs = "violin",
)
plot.data <- p$data
colnames(plot.data) <- c("Subclone_all", "Module.iCMS2")
plot.data <- plot.data %>%
    mutate(Patient = gsub("_.*", "", Subclone_all),
    Subclone = gsub("^.+?_", "", Subclone_all)) %>%
    filter(Patient %in% c("COAD12", "COAD17", "COAD24", "COAD33"))

pdf("All_patient/Violin.all.Module.iCMS2.pdf", 5, 3)
ggplot(plot.data, aes(x = Subclone, y = Module.iCMS2, fill = Subclone)) +
    geom_violin() +
    ggpubr::stat_compare_means(
        comparisons = list(
            c("subclone_1", "subclone_2")
        ), size = 2
    ) +
    facet_grid(cols = vars(Patient), scales = "free", space = "free") +
    theme_classic()
dev.off()

# TF activity
p <- plotGroups(
    ArchRProj = proj_Epi,
    groupBy = "Subclone_all",
    colorBy = "MotifMatrix",
    name = TF.selected.name[5:6],
    plotAs = "violin",
)

plot.data <- p[[1]]$data
colnames(plot.data) <- c("Subclone_all", "Activity")
plot.data <- plot.data %>%
    mutate(Patient = gsub("_.*", "", Subclone_all),
    Subclone = gsub("^.+?_", "", Subclone_all)) %>%
    filter(Patient %in% c("COAD12", "COAD17", "COAD24", "COAD33"))

pdf("All_patient/Violin.all.PPARA.pdf", 5, 3)
ggplot(plot.data, aes(x = Subclone, y = Activity, fill = Subclone)) +
    geom_violin() +
    ggpubr::stat_compare_means(
        comparisons = list(
            c("subclone_1", "subclone_2")
        ), size = 2
    ) +
    facet_grid(cols = vars(Patient), scales = "free", space = "free") +
    theme_classic()
dev.off()

plot.data <- p[[2]]$data
colnames(plot.data) <- c("Subclone_all", "Activity")
plot.data <- plot.data %>%
    mutate(Patient = gsub("_.*", "", Subclone_all),
    Subclone = gsub("^.+?_", "", Subclone_all)) %>%
    filter(Patient %in% c("COAD12", "COAD17", "COAD24", "COAD33"))

pdf("All_patient/Violin.all.HNF4A.pdf", 5, 3)
ggplot(plot.data, aes(x = Subclone, y = Activity, fill = Subclone)) +
    geom_violin() +
    ggpubr::stat_compare_means(
        comparisons = list(
            c("subclone_1", "subclone_2")
        ), size = 2
    ) +
    facet_grid(cols = vars(Patient), scales = "free", space = "free") +
    theme_classic()
dev.off()

## 4.2. heatmap for each patient ----
sePeaks <- getGroupSE(
    ArchRProj = proj_Epi,
    useMatrix = "PeakMatrix",
    groupBy = "Subclone_all",
    divideN = TRUE,
    scaleTo = NULL
)
PeakMatrix.subclone <- sePeaks@assays@data$PeakMatrix
rownames(PeakMatrix.subclone) <- rowData(sePeaks) %>%
    as.data.frame() %>%
    mutate("id" = paste(seqnames, start, end, sep = "_")) %>%
    pull(id)

for (patient in patient.selected) {
    message(paste("identify subclone markers for ", patient, " ...", sep = ""))
    marker.subclone.one <- getMarkerFeatures(
        ArchRProj = proj_Epi,
        useMatrix = "PeakMatrix",
        groupBy = "Subclone_all",
        testMethod = "wilcoxon",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        useGroups = grep(patient, unique(proj_Epi$Subclone_all), value = TRUE),
        bgdGroups = c("Normal")
    )

    marker.subclone <- getMarkers(marker.subclone.one,
        cutOff = "FDR <= 0.01 & Log2FC >= 1"
    ) %>%
        lapply(., as.data.frame)

    marker.subclone.id <- lapply(marker.subclone, function(x) {
        x <- x %>%
            mutate(
                peak_id = paste(seqnames, start, end, sep = "_")
            )
        return(x$peak_id)
    }) %>%
        unlist() %>%
            unique()

    plot.data <- PeakMatrix.subclone[
        marker.subclone.id,
        c("Normal", grep(patient, colnames(PeakMatrix.subclone), value = TRUE))
    ]
    colnames(plot.data) <- gsub("COAD.+?_", "", colnames(plot.data))
    plot.data <- apply(plot.data, 1, function(x) {
        x <- x - mean(x)
        x <- x / sd(x)
        return(x)
    })
    plot.data[plot.data > 1.5] <- 1.5
    plot.data[plot.data < (-1.5)] <- (-1.5)

    pdf(paste0("All_patient/Heatmap.marker.", patient, ".peaks.subclone.pdf"), 7, 3.5)
    pheatmap::pheatmap(plot.data,
        scale = "column",
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        show_colnames = FALSE,
        main = length(marker.subclone.id)
    )
    dev.off()
}

gc()
save.image("Epi_intratumor.RData")
