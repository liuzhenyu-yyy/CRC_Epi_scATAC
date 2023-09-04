if (TRUE) {
    library(ArchR)
    library(BSgenome)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(GenomicRanges)
    library(RColorBrewer)
    library(Seurat)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(UpSetR)
    library(corrplot)
    library(dplyr)
    library(genomation)
    library(ggplot2)
    library(org.Hs.eg.db)
    library(patchwork)
    library(pheatmap)
    library(viridis)
    library(Vennerable)
}

if (TRUE) {
    addArchRGenome("hg38")
    source("E:/LabWork/code/archr_track.R")
    # source("E:/LabWork/code/MyFunction.R")
    # ReadGeneInfo("hg38", classify.protein = FALSE)
}

project.dir.all <- "E:/LabWork/Project/CRC_NGS_ATAC/CRC_Epi_scATAC/Project_Dir_All"
project.dir.epi <- "E:/LabWork/Project/CRC_NGS_ATAC/CRC_Epi_scATAC/Project_Dir_Epi"

homer.parser <- function(res, log.p.value = 50, log2.enrichment = 1) {
    res <- fread(res,
        header = TRUE
    ) %>%
        as.data.frame() %>%
        .[order(.[, 1]), ]
    colnames(res) <- c(
        "Name", "Consensus", "p.value", "log.p.value", "q.value",
        "n.Targets", "Perc.Targets", "n.Background", "Perc.Background"
    )
    res$Perc.Targets <- res$Perc.Targets %>%
        gsub("%", "", .) %>%
        as.numeric()
    res$Perc.Background <- res$Perc.Background %>%
        gsub("%", "", .) %>%
        as.numeric()
    res$Log2_Enrichment <- log2(res$Perc.Targets / res$Perc.Background)
    res$log.p.value <- 0 - res$log.p.value
    res$TF <- gsub("\\(.+?$", "", res$Name) %>% toupper()
    res <- res[!duplicated(res$TF), ]
    rownames(res) <- res$TF

    anno <- read.csv("homer/motif.anno.csv", header = TRUE)
    rownames(anno) <- anno$TF
    res$Family <- anno[res$TF, "Anno"]

    # res <- res[order(res$log.p.value, decreasing = TRUE), ]
    # res <- res[!duplicated(res$Family), ]

    res$Diff <- "none"
    res$Diff[res$Log2_Enrichment >= log2.enrichment & res$log.p.value >= log.p.value] <- "up"
    return(res)
}