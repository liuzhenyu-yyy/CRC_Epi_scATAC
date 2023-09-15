setwd("E:/LabWork/Project/CRC_NGS_ATAC/CRC_Epi_scATAC/Results/05.Epi_TF_Clustering")
# load("05.Epi_TF_Clustering.RData")

source("../../code/00.Requirements.R")
rm(foo, homer.parser)

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
    theme_classic() +
    coord_fixed()
dev.off()
pdf("UMAP_MotifMat.CIMP_Group.pdf", 5, 4)
ggplot(sample.info.tumor, aes(x = UMAP_Motif_1, y = UMAP_Motif_2)) +
    geom_point(aes(color = CIMP_Group), size = 0.2) +
    scale_color_manual(values = mycolor$CIMP_Group) +
    theme_classic() +
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
plotTOM <- dissTOM ^ 10
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

pdf("Heatmap.cor.Module_Trait.pdf", 7, 4)
corrplot(
    as.matrix(cor.module.trait),
    method = "square",
    # addCoef.col = "black",
    tl.col = "black",
    tl.cex = 1,
    p.mat = as.matrix(p.module.trait),
    sig.level = c(0.001, 0.01, 0.05),
    pch.cex = 1.3,
    insig = "label_sig",
    col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
)
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
            scale_size_continuous(range = c(0, 5), limits = c(0, 2.5)) +
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
        vertex.label.font = 3,
        vertex.frame.color = "gray",
        vertex.label.color = "gray20",
        vertex.label = V(net.sub)$Symbol,
        vertex.label.cex = log10(V(net.sub)$GeneScore_Malignant) * 0.3,
        layout = layout_with_fr(net.sub),
        main = paste("Module ", one, sep = "")
    )
    dev.off()
}
rm(net.sub, one)

save.image("05.Epi_TF_Clustering.RData")
