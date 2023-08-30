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
