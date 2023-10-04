setwd("E:/LabWork/Project/CRC_NGS_ATAC/CRC_Epi_scATAC/Results/05.Epi_TF_Clustering")
load("05.Epi_TF_Clustering.RData")

source("../../code/00.Requirements.R")
rm(foo)

# 1. Reduction on Motif matrix ----
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

## 1.1. get motif devation matrix ----
proj_Epi <- loadArchRProject(project.dir.epi, force = TRUE)
table(proj_Epi$Epi_type)
table(proj_Epi$Sample)
getAvailableMatrices(proj_Epi)

MotifMat <- getMatrixFromProject(
    ArchRProj = proj_Epi,
    useMatrix = "MotifMatrix",
    useSeqnames = "deviations",
    binarize = FALSE
)
dim(MotifMat)

MotifMat <- MotifMat@assays@data$deviations
apply(MotifMat, 1, function(x) {
    sum(is.na(x))
}) %>% sort()
MotifMat <- MotifMat[rownames(MotifMat) != "ENSG00000250542_156", ]
MotifMat <- t(MotifMat)

sample.info.epi <- proj_Epi@cellColData %>% as.data.frame()
sample.info.tumor <- sample.info.epi %>%
    filter(Epi_type == "Malignant")

## 1.2. Reduction on Motif matrix ----
pca.MotifMat <- prcomp(MotifMat[rownames(sample.info.tumor), ],
    center = TRUE, scale. = TRUE
)
sum(pca.MotifMat$sdev^2)

pdf("Elbow.pca.MotifMat.pdf", 5, 5)
plot(pca.MotifMat$sdev[1:50]^2)
dev.off()

sample.info.tumor$PCA_Motif_1 <- pca.MotifMat$x[rownames(sample.info.tumor), 1]
sample.info.tumor$PCA_Motif_2 <- pca.MotifMat$x[rownames(sample.info.tumor), 2]
sample.info.tumor$PCA_Motif_3 <- pca.MotifMat$x[rownames(sample.info.tumor), 3]

pdf("PCA_MotifMat.Epi_Group.pdf", 12, 3.5)
ggplot(sample.info.tumor, aes(x = PCA_Motif_1, y = PCA_Motif_2)) +
    geom_point(aes(color = Epi_Group), size = 0.2) +
    scale_color_manual(values = mycolor$Epi_Group) +
    theme_classic() +
    NoLegend() +
    ggplot(sample.info.tumor, aes(x = PCA_Motif_1, y = PCA_Motif_3)) +
    geom_point(aes(color = Epi_Group), size = 0.2) +
    scale_color_manual(values = mycolor$Epi_Group) +
    theme_classic() +
    NoLegend() +
    ggplot(sample.info.tumor, aes(x = PCA_Motif_2, y = PCA_Motif_3)) +
    geom_point(aes(color = Epi_Group), size = 0.2) +
    scale_color_manual(values = mycolor$Epi_Group) +
    theme_classic()
dev.off()
pdf("PCA_MotifMat.CIMP_Group.pdf", 12, 3.5)
ggplot(sample.info.tumor, aes(x = PCA_Motif_1, y = PCA_Motif_2)) +
    geom_point(aes(color = CIMP_Group), size = 0.2) +
    scale_color_manual(values = mycolor$CIMP_Group) +
    theme_classic() +
    NoLegend() +
    ggplot(sample.info.tumor, aes(x = PCA_Motif_1, y = PCA_Motif_3)) +
    geom_point(aes(color = CIMP_Group), size = 0.2) +
    scale_color_manual(values = mycolor$CIMP_Group) +
    theme_classic() +
    NoLegend() +
    ggplot(sample.info.tumor, aes(x = PCA_Motif_2, y = PCA_Motif_3)) +
    geom_point(aes(color = CIMP_Group), size = 0.2) +
    scale_color_manual(values = mycolor$CIMP_Group) +
    theme_classic()
dev.off()

umap.MotifMat <- uwot::umap(pca.MotifMat$x[, 1:15],
    n_neighbors = 30,
    metric = "cosine",
    min_dist = 1
)
sample.info.tumor$UMAP_Motif_1 <- umap.MotifMat[rownames(sample.info.tumor), 1]
sample.info.tumor$UMAP_Motif_2 <- umap.MotifMat[rownames(sample.info.tumor), 2]

pdf("UMAP_MotifMat.Epi_Group.pdf", 4.5, 4)
ggplot(sample.info.tumor, aes(x = UMAP_Motif_1, y = UMAP_Motif_2)) +
    geom_point(aes(color = Epi_Group), size = 0.2) +
    scale_color_manual(values = mycolor$Epi_Group) +
    theme_bw() +
    theme(
        axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    ) +
    coord_fixed()
dev.off()
pdf("UMAP_MotifMat.CIMP_Group.pdf", 5, 4)
ggplot(sample.info.tumor, aes(x = UMAP_Motif_1, y = UMAP_Motif_2)) +
    geom_point(aes(color = CIMP_Group), size = 0.2) +
    scale_color_manual(values = mycolor$CIMP_Group) +
    theme_bw() +
    theme(
        axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    ) +
    coord_fixed()
dev.off()
pdf("UMAP_MotifMat.Sample.pdf", 7, 6)
ggplot(sample.info.tumor, aes(x = UMAP_Motif_1, y = UMAP_Motif_2)) +
    geom_point(aes(color = Sample), size = 0.2) +
    scale_color_manual(values = paletteDiscrete(sample.info.tumor$Sample)) +
    theme_bw() +
    theme(
        axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    ) +
    coord_fixed()
dev.off()
rm(umap.MotifMat)

## 1.3. clustering on Motif matrix ----
kmeans.MotifMat <- kmeans(pca.MotifMat$x[, 1:15],
    nstart = 100, centers = 6
)
sample.info.tumor$Kmeans_Cluster <- kmeans.MotifMat$cluster[rownames(sample.info.tumor)]
rm(kmeans.MotifMat)

pdf("UMAP_MotifMat.Kmeans_Cluster.pdf", 5, 4)
ggplot(sample.info.tumor, aes(x = UMAP_Motif_1, y = UMAP_Motif_2)) +
    geom_point(aes(color = as.character(Kmeans_Cluster)),
        size = 0.2
    ) +
    theme_classic() +
    coord_fixed()
dev.off()

aricode::AMI(sample.info.tumor$Kmeans_Cluster, sample.info.tumor$Epi_Group)
aricode::AMI(sample.info.tumor$Kmeans_Cluster, sample.info.tumor$CIMP_Group)

saveRDS(MotifMat, "MotifMat.RDS")
rm(MotifMat)

# 2. WGCNA on motifX-by-cluster matrix ----
## 2.1. get motif enrichment matrix ----
seMotif.cluster <- getGroupSE(
    ArchRProj = proj_Epi,
    useMatrix = "MotifMatrix",
    groupBy = "Clusters",
    divideN = TRUE,
    scaleTo = NULL
)
dim(seMotif.cluster)
rowData(seMotif.cluster)
seMotif.cluster <- seMotif.cluster[rowData(seMotif.cluster)$seqnames == "deviations", ]
rownames(seMotif.cluster@assays@data$MotifMatrix) <- rowData(seMotif.cluster)$name

MotifMat.cluster <- seMotif.cluster@assays@data$MotifMatrix
rownames(MotifMat.cluster)

MotifMat.cluster["ENSG00000250542_156", ]
MotifMat.cluster <- MotifMat.cluster[rownames(MotifMat.cluster) != "ENSG00000250542_156", ]
MotifMat.cluster <- t(MotifMat.cluster)

cluster.info <- readRDS("../04.Epi_CIMP/cluster.info.rds")
MotifMat.cluster <- MotifMat.cluster[rownames(cluster.info), ]
rm(seMotif.cluster)

## 2.2. motif meta data ----
# TF symbol
TF.info <- data.frame(
    row.names = colnames(MotifMat.cluster),
    "DBID" = colnames(MotifMat.cluster),
    "Symbol" = colnames(MotifMat.cluster) %>% gsub("\\_.+?$", "", .)
)
rownames(TF.info) <- TF.info$DBID

temp <- grep("NKX", TF.info$Symbol)
names(temp) <- TF.info$Symbol[temp]
names(temp) <- sapply(names(temp), function(x) {
    paste0(substr(x, 1, 4), "-", substr(x, 5, 6))
})
TF.info$Symbol[temp] <- names(temp)
rm(temp)

setdiff(TF.info$Symbol, rownames(seGene))
TF.info["MYCL1", ]$Symbol <- "MYCL"
TF.info["MLL", ]$Symbol <- "KMT2A"
TF.info["CGBP", ]$Symbol <- "CXXC1"
TF.info["PRKRIR", ]$Symbol <- "THAP12"
TF.info["T", ]$Symbol <- "TBXT"
TF.info["ZFP161", ]$Symbol <- "ZBTB14"
TF.info["CSDA", ]$Symbol <- "YBX3"
TF.info["ENSG00000229544", ]$Symbol <- "NKX1-2"
TF.info["C11orf9", ]$Symbol <- "MYRF"
TF.info["ZNF187", ]$Symbol <- "ZSCAN26"
TF.info["ZNF238", ]$Symbol <- "ZBTB18"

# TF Family
TF.CISBP <- fread("E:/LabWork/genome/hg38/CISBP/TF_Information.txt",
    sep = "\t", header = TRUE
) %>% as.data.frame()
TF.CISBP <- TF.CISBP[!duplicated(TF.CISBP$TF_Name), ]
rownames(TF.CISBP) <- gsub("_.+?$", "", TF.CISBP$TF_Name)

table(TF.info$Symbol %in% rownames(TF.CISBP)) # 79 missing

TF.info$Family <- TF.CISBP[TF.info$Symbol, ]$Family_Name
table(is.na(TF.info$Family)) # 79 missing
rm(TF.CISBP, temp)

# expression
proj_Epi$Epi_type
seGene <- getGroupSE(
    ArchRProj = proj_Epi,
    useMatrix = "GeneScoreMatrix",
    groupBy = "Epi_type",
    divideN = TRUE,
    scaleTo = NULL
)
rownames(seGene@assays@data$GeneScoreMatrix) <- rowData(seGene)$name
seGene <- seGene@assays@data$GeneScoreMatrix
seGene <- as.data.frame(seGene)

table(TF.info$Symbol %in% rownames(seGene)) # 73 missing
setdiff(TF.info$Symbol, rownames(seGene))

temp <- seGene[TF.info$Symbol, ]
colnames(temp) <- paste("GeneScore", colnames(temp), sep = "_")
TF.info <- cbind(TF.info, temp)
rm(seGene, temp)
table(is.na(TF.info$GeneScore_Normal), is.na(TF.info$Family))

## 2.3. network construction & module detection ----
library(WGCNA)
dir.create("WGCNA")

# analysis of network topology
sft <- pickSoftThreshold(MotifMat.cluster,
    powerVector = c(c(1:18)),
    verbose = 5
)
pdf("WGCNA/Dot.SoftThreshold.pdf", 8, 4)
ggplot(sft$fitIndices, aes(x = Power, y = -sign(slope) * SFT.R.sq)) +
    geom_text(aes(label = Power), size = 3) +
    xlab("Soft Threshold (power)") +
    ylab("Scale Free Topology Model Fit,signed R^2") +
    ggtitle("Scale independence") +
    geom_hline(yintercept = 0.80, linetype = 2, color = "#f87575") +
    theme_classic() +
    ggplot(sft$fitIndices, aes(x = Power, y = mean.k.)) +
    geom_text(aes(label = Power), size = 3) +
    xlab("Soft Threshold (power)") +
    ylab("Mean Connectivity") +
    ggtitle("Mean connectivity") +
    theme_classic()
dev.off()
rm(sft)

# network construction and module detection
net <- blockwiseModules(MotifMat.cluster,
    power = 7,
    networkType = "unsigned",
    TOMType = "signed",
    minModuleSize = 20,
    mergeCutHeight = 0.25,
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs = TRUE,
    saveTOMFileBase = "WGCNA/Motif.Cluster.TOM",
    verbose = 3
)
table(net$colors)
TF.info$Module <- net$colors[TF.info$DBID]

# Plot the dendrogram and the module colors
scales::show_col(ArchRPalettes$calm)
mycolor$TF_Module <- c("0" = "gray", ArchRPalettes$calm[1:13])

moduleColors <- labels2colors(net$colors,
    colorSeq = ArchR::ArchRPalettes$calm
)
table(moduleColors, net$colors)
pdf("WGCNA/Dendrogram.pdf", 6, 4)
plotDendroAndColors(net$dendrograms[[1]],
    moduleColors[net$blockGenes[[1]]],
    "TF Module colors",
    dendroLabels = FALSE, hang = 0.03,
    addGuide = TRUE, guideHang = 0.05
)
dev.off()

# plot correlation of all genes
dissTOM <- 1 - TOMsimilarityFromExpr(MotifMat.cluster, power = 7)
plotTOM <- dissTOM^10
diag(plotTOM) <- NA

pdf("WGCNA/Heatmap.network.all.pdf", 8, 8)
TOMplot(plotTOM,
    dendro = net$dendrograms[[1]],
    Colors = moduleColors,
    main = "Network heatmap plot, all genes",
    col = colorRampPalette(brewer.pal(9, "RdBu"))(100)
    # col = gplots::colorpanel(250, "red", "orange", "lemonchiffon")
)
dev.off()
rm(dissTOM, plotTOM)

## 2.4. clustering results and family ----
TF.info.sub <- TF.info %>% subset(!is.na(Family))
table(TF.info.sub$Family) %>% sort()

cM <- confusionMatrix(
    TF.info.sub$Family,
    TF.info.sub$Module
) %>%
    as.data.frame()
cM <- cM[rowSums(cM) >= 5, ]
cM <- cM / rowSums(cM)

pdf("WGCNA/Heatmap_cM_Family.pdf", 6, 6)
pheatmap(cM, scale = "none", clustering_method = "ward.D2")
dev.off()

cM <- confusionMatrix(
    TF.info.sub$Module,
    TF.info.sub$Family
) %>%
    as.data.frame()
cM <- cM[rowSums(cM) >= 5, ]
cM <- cM / rowSums(cM)

pdf("WGCNA/Heatmap_cM_Module.pdf", 8, 5)
pheatmap(cM, scale = "none", clustering_method = "ward.D2")
dev.off()
rm(TF.info.sub, cM)

## 2.5. output network and meta data ----
TOM <- TOMsimilarityFromExpr(MotifMat.cluster, power = 7)

cyt <- exportNetworkToCytoscape(TOM,
    edgeFile = paste("WGCNA/CytoscapeInput-edges-all.txt", sep = ""),
    nodeFile = paste("WGCNA/CytoscapeInput-nodes-all.txt", sep = ""),
    weighted = TRUE,
    threshold = 0.02,
    nodeNames = colnames(MotifMat.cluster),
    nodeAttr = paste0("ME", net$colors)
)
rm(cyt, TOM)
write.table(TF.info, "WGCNA/TF.info.csv",
    sep = ",", quote = FALSE, row.names = FALSE
)

# 3. TF module and clinical trait & cancer subtype ----
## 3.1. Module & Trait correlation ----
MEs <- orderMEs(net$MEs, greyName = "ME0")

colnames(cluster.info)
datTraits <- cluster.info %>%
    select(c("Gender_Major", "MSI_Status_Major", "Side_Major", "Epi_Group", "CIMP_Group"))
datTraits$Gender_Major <- ifelse(datTraits$Gender_Major == "Male", 1, 0)
datTraits$MSI_Status_Major <- ifelse(datTraits$MSI_Status_Major == "MSS", 1, 0)
datTraits$MSI_Status_Major[cluster.info$MSI_Status_Major == "NA"] <- NA
datTraits$Side_Major <- ifelse(datTraits$Side_Major == "Right", 1, 0)
datTraits$Epi_Group <- ifelse(datTraits$Epi_Group == "Group_1", 1, 0)
datTraits$CIMP_Group <- factor(datTraits$CIMP_Group,
    levels = c("CIMP_Negative", "CIMP_Low", "CIMP_High")
) %>% as.numeric()
colnames(datTraits) <- gsub("_Major", "", colnames(datTraits))

identical(rownames(datTraits), rownames(MEs))
cor.module.trait <- cor(MEs, datTraits, use = "p")
p.module.trait <- corPvalueStudent(cor.module.trait, nrow(MotifMat.cluster))
cor.module.trait <- t(cor.module.trait) %>% as.data.frame()
p.module.trait <- t(p.module.trait) %>% as.data.frame()

orders <- gsub("ME", "", colnames(cor.module.trait)) %>%
    as.numeric() %>%
    order()
pdf("WGCNA/Heatmap.cor.Module_Trait1.pdf", 3, 5)
corrplot(
    t(as.matrix(cor.module.trait[, orders])),
    method = "square",
    # addCoef.col = "black",
    tl.col = "black",
    tl.cex = 1,
    p.mat = t(as.matrix(p.module.trait[, orders])),
    sig.level = c(0.001, 0.01, 0.05),
    pch.cex = 1.3,
    insig = "label_sig",
    col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
)
dev.off()
rm(orders)
identical(rownames(cluster.info), rownames(MEs))
plot.data <- cbind(cluster.info %>%
    select(c("MSI_Status_Major", "Epi_Group", "CIMP_Group")), MEs)
table(plot.data$MSI_Status_Major)

pdf("WGCNA/Box.Module_MSI.pdf", 2.5, 6)
plot.data %>%
    select(c("MSI_Status_Major", "ME9", "ME12")) %>%
    filter(!MSI_Status_Major == "NA") %>%
    melt() %>%
    mutate("Module" = gsub("ME", "Module ", variable)) %>%
    ggplot(aes(x = MSI_Status_Major, y = value)) +
    geom_boxplot(aes(fill = MSI_Status_Major), show.legend = FALSE) +
    scale_fill_manual(values = mycolor$MSI) +
    ggpubr::stat_compare_means(
        label = "p.format", label.x.npc = "center"
    ) +
    theme_classic() +
    facet_wrap(~Module, nrow = 2) +
    ylab("Module Eigenvalue")
dev.off()

pdf("WGCNA/Box.Module_Epi_Group.pdf", 2.5, 6)
plot.data %>%
    select(c("Epi_Group", "ME5", "ME8")) %>%
    # filter(!MSI_Status_Major == "NA") %>%
    melt() %>%
    mutate("Module" = gsub("ME", "Module ", variable)) %>%
    ggplot(aes(x = Epi_Group, y = value)) +
    geom_boxplot(aes(fill = Epi_Group), show.legend = FALSE) +
    scale_fill_manual(values = mycolor$Epi_Group) +
    ggpubr::stat_compare_means(
        label = "p.format", label.x.npc = "center"
    ) +
    theme_classic() +
    facet_wrap(~Module, nrow = 2) +
    ylab("Module Eigenvalue")
dev.off()

pdf("WGCNA/Box.Module_CIMP_Group.pdf", 3.5, 6)
plot.data %>%
    select(c("CIMP_Group", "ME11", "ME8")) %>%
    # filter(!MSI_Status_Major == "NA") %>%
    melt() %>%
    mutate("Module" = gsub("ME", "Module ", variable)) %>%
    ggplot(aes(x = CIMP_Group, y = value)) +
    geom_boxplot(aes(fill = CIMP_Group), show.legend = FALSE) +
    scale_fill_manual(values = mycolor$CIMP_Group) +
    ggpubr::stat_compare_means(
        comparisons = list(c("CIMP_Negative", "CIMP_High")),
        label = "p.format", label.x.npc = "center"
    ) +
    theme_classic() +
    facet_wrap(~Module, nrow = 2) +
    ylab("Module Eigenvalue")
dev.off()

## 3.2. Calculate Gene Trait Significance and Module Membership ----
identical(rownames(MotifMat.cluster), rownames(MEs))

# Module Membership
cor.gene.module <- cor(MotifMat.cluster, MEs, use = "p") %>% as.data.frame()
p.gene.module <- corPvalueStudent(
    as.matrix(cor.gene.module),
    nrow(MotifMat.cluster)
) %>% as.data.frame()

# Gene Trait Significance
identical(rownames(MotifMat.cluster), rownames(datTraits))
cor.gene.traits <- cor(MotifMat.cluster, datTraits, use = "p") %>% as.data.frame()
p.gene.traits <- corPvalueStudent(
    as.matrix(cor.gene.traits),
    nrow(MotifMat.cluster)
) %>% as.data.frame()

## 3.3. GS MM correlation ----
dir.create("GS_MM")

ME.selected <- c(9, 8, 8, 11, 12, 5)
Trait.selected <- c("MSI_Status", "Epi_Group", "CIMP_Group", "CIMP_Group", "MSI_Status", "Epi_Group")

temp <- table(TF.info$Family) %>%
    sort(decreasing = TRUE) %>%
    head(20) %>%
    names() %>%
    gsub("Unknown", "Others", .)

mycolor$TF_Family <- ArchRPalettes$calm[c(1:4, 7, 6, 8:20, 5)]
names(mycolor$TF_Family) <- temp
mycolor$TF_Family["Others"] <- "gray50"
scales::show_col(mycolor$TF_Family)

for (i in seq_along(ME.selected)) {
    gene.selected <- net$colors[net$colors == ME.selected[i]] %>% names()
    plot.data <- data.frame(
        row.names = gene.selected,
        "MM" = cor.gene.module[gene.selected, paste("ME", ME.selected[i], sep = "")],
        "GS" = cor.gene.traits[gene.selected, Trait.selected[i]]
    )

    plot.data <- cbind(plot.data, TF.info[gene.selected, ])
    if (sum(is.na(plot.data$GeneScore_Malignant)) > 0) {
        plot.data[is.na(plot.data$GeneScore_Malignant), ]$GeneScore_Malignant <- 0
    }
    plot.data[plot.data$GeneScore_Malignant < 30, ]$Symbol <- NA

    if (sum(is.na(plot.data$Family)) > 0) {
        plot.data[is.na(plot.data$Family), ]$Family <- "Others"
    }
    if (sum(!plot.data$Family %in% names(mycolor$TF_Family)) > 0) {
        plot.data[!plot.data$Family %in% names(mycolor$TF_Family), ]$Family <- "Others"
    }
    plot.data$Family <- factor(plot.data$Family, levels = names(mycolor$TF_Family))
    p <- ggplot(plot.data, aes(x = MM, y = GS)) +
        geom_point(aes(color = Family, size = log10(GeneScore_Malignant))) +
        ggrepel::geom_text_repel(aes(label = Symbol, color = Family),
            size = 2, max.overlaps = 20, point.padding = 1
        ) +
        scale_color_manual(values = mycolor$TF_Family) +
        ggpubr::stat_cor() +
        xlab(paste("Gene Module Membership in Module ", ME.selected[i], sep = "")) +
        ylab(paste("Gene Significance for ", Trait.selected[i])) +
        xlim(c(min(plot.data$MM) - 0.2, max(plot.data$MM) + 0.2)) +
        ylim(c(min(plot.data$GS) - 0.2, max(plot.data$GS) + 0.2)) +
        theme_classic() +
        scale_size_continuous(range = c(0, 5), limits = c(0, 2.5))
    pdf(paste("GS_MM/Dot.Module", ME.selected[i], "_", Trait.selected[i], ".pdf", sep = ""), 6, 4)
    plot(p)
    dev.off()
}

rm(plot.data, p, i, gene.selected, temp)

# 4. Vlsualization of the network ----
library("igraph")
dir.create("Network")

net.meta <- list(
    nodes = read.table("WGCNA/CytoscapeInput-nodes-all.txt",
        header = TRUE, as.is = TRUE, sep = "\t"
    )[, c(1, 3)],
    links = read.csv("WGCNA/CytoscapeInput-edges-all.txt",
        header = TRUE, as.is = TRUE, sep = "\t"
    )[, 1:3]
)
colnames(net.meta$nodes) <- c("ID", "ME")
net.meta$nodes <- cbind(net.meta$nodes, TF.info[net.meta$nodes$ID, ])
net.meta$nodes[!net.meta$nodes$Family %in% names(mycolor$TF_Family), ]$Family <- "Others"
table(net.meta$nodes$Family)

net.all <- graph_from_data_frame(
    d = net.meta$links,
    vertices = net.meta$nodes,
    directed = FALSE
)
table(V(net.all)$Family)
net.all <- simplify(net.all, remove.multiple = FALSE, remove.loops = TRUE)

V(net.all)$degree <- degree(net.all, mode = "all")
net.all <- induced_subgraph(net.all,
    vids = V(net.all)[V(net.all)$degree > 5]
)
net.all <- induced_subgraph(net.all,
    vids = V(net.all)[!is.na(V(net.all)$GeneScore_Malignant)]
)

l <- layout_with_fr(net.all)
pdf("Network/Net.all.module.pdf", 10, 10)
plot(net.all,
    edge.color = "gray80",
    edge.width = E(net.all)$weight * 0.1,
    vertex.size = log1p(V(net.all)$GeneScore_Malignant) * 1,
    vertex.color = mycolor$TF_Module[as.character(V(net.all)$Module)],
    vertex.label = NA,
    layout = l,
    main = "All TFs"
)
dev.off()
pdf("Network/Net.all.Family.pdf", 10, 10)
plot(net.all,
    edge.color = "gray80",
    edge.width = E(net.all)$weight * 0.1,
    vertex.size = log1p(V(net.all)$GeneScore_Malignant) * 1,
    vertex.color = mycolor$TF_Family[as.character(V(net.all)$Family)],
    vertex.label = NA,
    layout = l,
    main = "All TFs"
)
dev.off()

extrafont::loadfonts()
for (one in unique(ME.selected)) {
    net.sub <- induced_subgraph(net.all,
        vids = V(net.all)[V(net.all)$Module == one]
    )
    V(net.sub)$degree <- degree(net.sub, mode = "all")
    net.sub <- induced_subgraph(net.sub,
        vids = V(net.sub)[V(net.sub)$degree > 5]
    )
    pdf(paste("Network/Net.Module", one, ".Family.pdf", sep = ""), 7, 7)
    plot(net.sub,
        edge.color = "gray80",
        edge.width = E(net.sub)$weight * 0.1,
        vertex.size = log1p(V(net.sub)$GeneScore_Malignant) * 3,
        vertex.color = mycolor$TF_Family[as.character(V(net.sub)$Family)],
        vertex.label.family = "Arial",
        vertex.label.font = 4,
        vertex.frame.color = "gray",
        vertex.label.color = "black",
        vertex.label = V(net.sub)$Symbol,
        vertex.label.cex = log10(V(net.sub)$GeneScore_Malignant) * 0.4,
        layout = layout_with_fr(net.sub) * 0.6,
        main = paste("Module ", one, sep = "")
    )
    dev.off()
}
rm(net.sub, one)

# 5. TF enrichment in each patient ----
dir.create("TF_motif_cluster")

## 5.1 identify diff peaks in each cluster & run Homers ----
table(proj_Epi$Clusters)
dir.create("diff_peak_cluster")
dir.create("diff_peak_cluster/bed")

for (one in rownames(cluster.info)) {
    marker.peak.one <- getMarkerFeatures(
        ArchRProj = proj_Epi,
        useMatrix = "PeakMatrix",
        groupBy = "Clusters",
        testMethod = "wilcoxon",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        useGroups = one,
        bgdGroups = c("C3", "C4")
    )

    marker.peak.one.df <- getMarkers(marker.peak.one,
        cutOff = "FDR <= 0.01 & Log2FC >= 1"
    )[[1]] %>% as.data.frame()
    write.table(
        marker.peak.one.df,
        paste0("diff_peak_cluster/marker.peak.tumor.", one, ".tsv"),
        col.names = TRUE,
        sep = "\t", quote = FALSE, row.names = FALSE
    )

    marker.peak.one.df <- marker.peak.one.df  %>%
        select(seqnames, start, end) %>%
        mutate(
            id = paste(seqnames, start, end, sep = "_"),
            length = 500,
            strand = "."
        )
    write.table(
        marker.peak.one.df,
        paste0("diff_peak_cluster/bed/marker.peak.tumor.", one, ".bed"),
        col.names = FALSE,
        sep = "\t", quote = FALSE, row.names = FALSE
    )
}
rm(one, marker.peak.one.df, marker.peak.one, test)

# find motifs in each patient
# sh Run.Homoer.Motif.sh

## 5.2. cluster-wise analysis of module genes ----
homer.res.cluster <- list()

sapply(rownames(cluster.info), function(one) {
    homer.res <- homer.parser(paste0("homer/", one, "/knownResults.txt"))
    homer.res$Cluster <- one
    homer.res.cluster[[one]] <<- homer.res
    return(1)
})
sapply(homer.res.cluster, dim)

homer.res.cluster <- do.call(rbind, homer.res.cluster)
table(homer.res.cluster$Cluster)

for (i in seq_along(ME.selected)) {
    gene.selected <- net$colors[net$colors == ME.selected[i]] %>%
        names() %>%
        gsub("_.+?$", "", .)
    plot.data <- homer.res.cluster %>%
        filter(TF %in% gene.selected) %>%
        mutate(
            Epi_Group = cluster.info[.$Cluster, "Epi_Group"],
            CIMP_Group = cluster.info[.$Cluster, "CIMP_Group"],
            MSI_Status = cluster.info[.$Cluster, "MSI_Status_Major"]
        )

    # plot.data <- plot.data[order(plot.data$Epi_Group, plot.data$Cluster), ]
    # plot.data$Cluster <- factor(plot.data$Cluster, levels = unique(plot.data$Cluster))
    # plot.data$TF <- factor(plot.data$TF, levels = rev(TF.selected))

    if (max(plot.data$log.p.value) > 500) {
        plot.data[plot.data$log.p.value > 500, ]$log.p.value <- 500
    }
    p <- ggplot(plot.data) +
        geom_point(aes(
            x = Cluster, y = TF,
            fill = log.p.value, size = Log2_Enrichment
        ), pch = 21) +
        scale_fill_viridis_c() +
        facet_grid(
            cols = vars(get(Trait.selected[i])),
            rows = vars(paste0("Module", ME.selected[i])),
            scales = "free", space = "free"
        ) +
        theme_bw()
    pdf(paste("TF_motif_cluster/Dot.", Trait.selected[i], "_Module", ME.selected[i], ".pdf", sep = ""), 9, 4)
    plot(p)
    dev.off()
}
rm(p, plot.data, i, gene.selected)

## 5.3. detialed visulization of specific module ----
# Module 8: Epi_type
net$colors[grep(paste(c("ELF1", "EHF", "ETS1", "FOXM1", "MAFK", "FOXA3"), collapse = "|"),
    names(net$colors),
    value = TRUE
)]
TF.info[grep("SOX", TF.info$Symbol), ] %>%
    .[order(.$GeneScore_Malignant), c(1:4, 6:7)]

gene.selected <- net$colors[net$colors == 8] %>%
    names() %>%
    gsub("_.+?$", "", .) %>%
    c(., "FOXM1", "MAFK", "FOXA3") %>%
    setdiff(., c(
        "GSC", "FOXH1", "ESRRB", "CRX", # low enrichment
        "SOX10", "SOX15", "SOX3", "SOX7" # low expression
    ))
plot.data <- homer.res.cluster %>%
    filter(TF %in% gene.selected) %>%
    mutate(
        Epi_Group = cluster.info[.$Cluster, "Epi_Group"],
        Module = ifelse(TF %in% (net$colors[net$colors == 8] %>%
            names() %>%
            gsub("_.+?$", "", .)), "Module8", "Others") %>%
            factor(., levels = rev(c("Module8", "Others")))
    )
plot.data[plot.data$log.p.value > 500, ]$log.p.value <- 500

pdf("TF_motif_cluster/Dot.Selected.Module8_EpiGroup1.pdf", 9, 3.5)
ggplot(plot.data) +
    geom_point(aes(
        x = Cluster, y = TF,
        fill = log.p.value, size = Log2_Enrichment
    ), pch = 21) +
    scale_fill_viridis_c() +
    facet_grid(
        cols = vars(Epi_Group),
        rows = vars(Module),
        scales = "free", space = "free"
    ) +
    theme_bw()
dev.off()

# Module 11: CIMP_Group
gene.selected <- net$colors[net$colors == 11] %>%
    names() %>%
    gsub("_.+?$", "", .) %>%
    setdiff(., c("KLF1"))
plot.data <- homer.res.cluster %>%
    filter(TF %in% gene.selected) %>%
    mutate(
        CIMP_Group = cluster.info[.$Cluster, "CIMP_Group"]
    )
plot.data[plot.data$log.p.value > 300, ]$log.p.value <- 300

pdf("TF_motif_cluster/Dot.Selected.Module11_CIMP_Group.pdf", 9, 3.5)
ggplot(plot.data) +
    geom_point(aes(
        x = Cluster, y = TF,
        fill = log.p.value, size = Log2_Enrichment
    ), pch = 21) +
    scale_fill_viridis_c() +
    facet_grid(
        cols = vars(CIMP_Group),
        rows = vars("Module8"),
        scales = "free", space = "free"
    ) +
    theme_bw()
dev.off()

# Module 9: MSI
gene.selected <- net$colors[net$colors == 9] %>%
    names() %>%
    gsub("_.+?$", "", .) %>%
    setdiff(., c("ZKSCAN1", "NANOG", "MYB", "MAZ", "HOXD12", "HOXC9"))
plot.data <- homer.res.cluster %>%
    filter(TF %in% gene.selected) %>%
    mutate(
        MSI_Status = cluster.info[.$Cluster, "MSI_Status_Major"]
    ) %>%
    filter(MSI_Status != "NA")
plot.data[plot.data$log.p.value > 500, ]$log.p.value <- 500

pdf("TF_motif_cluster/Dot.Selected.Module9_MSI.pdf", 8, 3)
ggplot(plot.data) +
    geom_point(aes(
        x = Cluster, y = TF,
        fill = log.p.value, size = Log2_Enrichment
    ), pch = 21) +
    scale_fill_viridis_c() +
    facet_grid(
        cols = vars(MSI_Status),
        rows = vars("Module9"),
        scales = "free", space = "free"
    ) +
    theme_bw()
dev.off()
rm(plot.data, gene.selected, p)

# 6. same TF, different peaks in each cluster ----
## 6.1. get data ----
# significant peaks of each clusters
peaks.clusters <- list()
sapply(rownames(cluster.info), function(one) {
    temp <- read.table(paste0("diff_peak_cluster/marker.peak.tumor.", one, ".tsv"), header = TRUE) %>%
        mutate(
            id = paste(seqnames, start, end, sep = "_")
        )
    peaks.clusters[[one]] <<- temp$id
    return(1)
})
sapply(peaks.clusters, length)

# motif match matrix
motif.match <- getMatches(ArchRProj = proj_Epi, name = "Motif")
temp <- rowRanges(motif.match)
motif.match <- motif.match@assays@data$matches

rownames(motif.match) <- paste(seqnames(temp), start(temp), end(temp), sep = "_")
colnames(motif.match) <- gsub("_.+?$", "", colnames(motif.match))
motif.match[1:5, 1:5]
rm(temp)

## 6.2. plot diff peak with motif in each cluster ----
peaks.common <- list()

# Group2 TFs: HNF4A, PPARA
cluster.selected <- cluster.info %>%
    filter(Epi_Group == "Group_2") %>%
    rownames()

# HNF4A
peak.selected <- motif.match[, "HNF4A"] %>%
    .[.] %>%
    names() %>%
    intersect(., unlist(peaks.clusters[cluster.selected]) %>% unique())
length(peak.selected) # 20059 HNF4A peaks

plot.data <- matrix(0, nrow = length(peak.selected), ncol = length(cluster.selected))
rownames(plot.data) <- peak.selected
colnames(plot.data) <- cluster.selected
for (one in cluster.selected) {
    plot.data[peak.selected %in% peaks.clusters[[one]], one] <- 1
}

pdf("TF_motif_cluster/Heatmap.HNF4A.clusterPeaks.pdf", 10, 6)
p <- pheatmap(
    t(plot.data)[c(
        "C18", "C2", "C16", "C13", "C22", "C17", "C20",
        "C12", "C21", "C19", "C7", "C23", "C14", "C15", "C8"
    ), ],
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_row = cluster.info %>% select(c("Epi_Group")),
    annotation_color = mycolor,
    color = colorRampPalette(c("gray80", "#cd2525"))(2),
    clustering_method = "ward.D2",
    cutree_cols = 17
)
dev.off()

temp <- cutree(p$tree_col, 17)
rowSums(plot.data[names(temp), ]) %>%
    aggregate(., by = list(temp), FUN = mean) # 4 and 15
peaks.common[["HNF4A"]] <- names(temp)[temp %in% c(4, 15)]

anno.col <- data.frame(
    row.names = rownames(plot.data),
    "con" = ifelse(rownames(plot.data) %in% peaks.common[["HNF4A"]], "Consensus", "Specific")
)

pdf("TF_motif_cluster/Heatmap.HNF4A.clusterPeaks.pdf", 10, 6)
p <- pheatmap(
    t(plot.data)[c(
        "C18", "C2", "C16", "C13", "C22", "C17", "C20",
        "C12", "C21", "C19", "C7", "C23", "C14", "C15", "C8"
    ), ],
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_row = cluster.info %>% select(c("Epi_Group")),
    annotation_col = anno.col,
    annotation_color = c(mycolor, list(con = c("Consensus" = "#86d786", "Specific" = "#f6be43"))),
    color = colorRampPalette(c("gray80", "#cd2525"))(2),
    clustering_method = "ward.D2",
    cutree_cols = 17
)
dev.off()

# PPARA
peak.selected <- motif.match[, "PPARA"] %>%
    .[.] %>%
    names() %>%
    intersect(., unlist(peaks.clusters[cluster.selected]) %>% unique())
length(peak.selected) # 9720 HNF4A peaks

plot.data <- matrix(0, nrow = length(peak.selected), ncol = length(cluster.selected))
rownames(plot.data) <- peak.selected
colnames(plot.data) <- cluster.selected
for (one in cluster.selected) {
    plot.data[peak.selected %in% peaks.clusters[[one]], one] <- 1
}

pdf("TF_motif_cluster/Heatmap.PPARA.clusterPeaks.pdf", 10, 6)
p <- pheatmap(
    t(plot.data)[c(
        "C18", "C2", "C16", "C13", "C22", "C17", "C20",
        "C12", "C21", "C19", "C7", "C23", "C14", "C15", "C8"
    ), ],
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_row = cluster.info %>% select(c("Epi_Group")),
    annotation_color = mycolor,
    color = colorRampPalette(c("gray80", "#cd2525"))(2),
    clustering_method = "ward.D2",
    cutree_cols = 15
)
dev.off()

temp <- cutree(p$tree_col, 15)
rowSums(plot.data[names(temp), ]) %>%
    aggregate(., by = list(temp), FUN = mean) # 1 and 4
peaks.common[["PPARA"]] <- names(temp)[temp %in% c(1, 4)]

anno.col <- data.frame(
    row.names = rownames(plot.data),
    "con" = ifelse(rownames(plot.data) %in% peaks.common[["PPARA"]], "Consensus", "Specific")
)

pdf("TF_motif_cluster/Heatmap.PPARA.clusterPeaks.pdf", 10, 6)
p <- pheatmap(
    t(plot.data)[c(
        "C18", "C2", "C16", "C13", "C22", "C17", "C20",
        "C12", "C21", "C19", "C7", "C23", "C14", "C15", "C8"
    ), ],
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_row = cluster.info %>% select(c("Epi_Group")),
    annotation_col = anno.col,
    annotation_color = c(mycolor, list(con = c("Consensus" = "#86d786", "Specific" = "#f6be43"))),
    color = colorRampPalette(c("gray80", "#cd2525"))(2),
    clustering_method = "ward.D2",
    cutree_cols = 15
)
dev.off()

# Group1 TFs: SOX4, FOXA3
# SOX4
cluster.selected <- cluster.info %>%
    filter(Epi_Group == "Group_1") %>%
    rownames()

peak.selected <- motif.match[, "SOX4"] %>%
    .[.] %>%
    names() %>%
    intersect(., unlist(peaks.clusters[cluster.selected]) %>% unique())
length(peak.selected) # 8540 SOX4 peaks

plot.data <- matrix(0, nrow = length(peak.selected), ncol = length(cluster.selected))
rownames(plot.data) <- peak.selected
colnames(plot.data) <- cluster.selected
for (one in cluster.selected) {
    plot.data[peak.selected %in% peaks.clusters[[one]], one] <- 1
}

pdf("TF_motif_cluster/Heatmap.SOX4.clusterPeaks.pdf", 10, 4.5)
p <- pheatmap(
    t(plot.data)[c(
        "C10", "C11", "C26", "C29", "C27", "C1", "C6",
        "C24", "C25", "C5"
    ), ],
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_row = cluster.info %>% select(c("Epi_Group")),
    annotation_color = mycolor,
    color = colorRampPalette(c("gray80", "#cd2525"))(2),
    clustering_method = "ward.D2",
    cutree_cols = 11
)
dev.off()

temp <- cutree(p$tree_col, 11)
rowSums(plot.data[names(temp), ]) %>%
    aggregate(., by = list(temp), FUN = mean) # 8
peaks.common[["SOX4"]] <- names(temp)[temp %in% c(8)]

anno.col <- data.frame(
    row.names = rownames(plot.data),
    "con" = ifelse(rownames(plot.data) %in% peaks.common[["SOX4"]], "Consensus", "Specific")
)

pdf("TF_motif_cluster/Heatmap.SOX4.clusterPeaks.pdf", 10, 4.5)
p <- pheatmap(
    t(plot.data)[c(
        "C10", "C11", "C26", "C29", "C27", "C1", "C6",
        "C24", "C25", "C5"
    ), ],
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_row = cluster.info %>% select(c("Epi_Group")),
    annotation_col = anno.col,
    annotation_color = c(mycolor, list(con = c("Consensus" = "#86d786", "Specific" = "#f6be43"))),
    color = colorRampPalette(c("gray80", "#cd2525"))(2),
    clustering_method = "ward.D2",
    cutree_cols = 11
)
dev.off()

# FOXA3
peak.selected <- motif.match[, "FOXA3"] %>%
    .[.] %>%
    names() %>%
    intersect(., unlist(peaks.clusters[cluster.selected]) %>% unique())
length(peak.selected) # 18275 SOX4 peaks

plot.data <- matrix(0, nrow = length(peak.selected), ncol = length(cluster.selected))
rownames(plot.data) <- peak.selected
colnames(plot.data) <- cluster.selected
for (one in cluster.selected) {
    plot.data[peak.selected %in% peaks.clusters[[one]], one] <- 1
}

pdf("TF_motif_cluster/Heatmap.FOXA3.clusterPeaks.pdf", 10, 4.5)
p <- pheatmap(
    t(plot.data)[c(
        "C10", "C11", "C26", "C29", "C27", "C1", "C6",
        "C24", "C25", "C5"
    ), ],
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_row = cluster.info %>% select(c("Epi_Group")),
    annotation_color = mycolor,
    color = colorRampPalette(c("gray80", "#cd2525"))(2),
    clustering_method = "ward.D2",
    cutree_cols = 11
)
dev.off()

temp <- cutree(p$tree_col, 11)
rowSums(plot.data[names(temp), ]) %>%
    aggregate(., by = list(temp), FUN = mean) # 8
peaks.common[["FOXA3"]] <- names(temp)[temp %in% c(8)]

anno.col <- data.frame(
    row.names = rownames(plot.data),
    "con" = ifelse(rownames(plot.data) %in% peaks.common[["FOXA3"]], "Consensus", "Specific")
)

pdf("TF_motif_cluster/Heatmap.FOXA3.clusterPeaks.pdf", 10, 4.5)
p <- pheatmap(
    t(plot.data)[c(
        "C10", "C11", "C26", "C29", "C27", "C1", "C6",
        "C24", "C25", "C5"
    ), ],
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_row = cluster.info %>% select(c("Epi_Group")),
    annotation_col = anno.col,
    annotation_color = c(mycolor, list(con = c("Consensus" = "#86d786", "Specific" = "#f6be43"))),
    color = colorRampPalette(c("gray80", "#cd2525"))(2),
    clustering_method = "ward.D2",
    cutree_cols = 11
)
dev.off()

rm(p, peak.selected, cluster.selected, plot.data, motif.match, one, temp, anno.col)

## 6.3. clustering on consensus peak genes ----
# get consensus peaks to genes
sapply(peaks.common, length)
peakset <- proj_Epi@peakSet
names(peakset) <- paste(seqnames(peakset), start(peakset), end(peakset), sep = "_")
table(peakset$peakType)

peaks.common.2gene <- sapply(peaks.common, function(one) {
    temp <- peakset[one]
    temp <- temp[temp$peakType != "Distal"]
    return(unique(temp$nearestGene))
})
sapply(peaks.common.2gene, length)
peaks.common.2gene <- unlist(peaks.common.2gene) %>% unique()
peaks.common.2gene <- peaks.common.2gene[!is.na(peaks.common.2gene)]
length(peaks.common.2gene) # 1582 geness
rm(peakset)
saveRDS(peaks.common.2gene, "peaks.common.2gene.rds")

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
cGSMat <- cGSMat[peaks.common.2gene, ]

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
pdf("TF_motif_cluster/PCA_MotifMat.CIMP_Group.pdf", 12, 3.5)
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

pdf("TF_motif_cluster/UMAP_MotifMat.Epi_Group.pdf", 4.5, 4)
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
pdf("TF_motif_cluster/UMAP_MotifMat.CIMP_Group.pdf", 5, 4)
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
pdf("TF_motif_cluster/UMAP_MotifMat.Sample.pdf", 7, 6)
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

gc()
save.image("05.Epi_TF_Clustering.RData")
