setwd("E:/LabWork/Project/CRC_NGS_ATAC/CRC_Epi_scATAC/Results/03.Epi_AD_Methylation")
# load("03.Epi_AD_Methylation.RData")
source("../../code/00.Requirements.R")

if (TRUE) {
    library(dplyr)
    library(minfi)
    library(rtracklayer)
    library(ChIPseeker)
}

# 1. load data ----
## 1.1. sample meta data ----
sample.info.array <- read.table("data/sample.info.array.txt", header = FALSE, sep = "\t")
colnames(sample.info.array) <- c("ID", "Sample", "Number")
rownames(sample.info.array) <- sample.info.array$ID
table(sample.info.array$Sample)

sample.info.array$Location <- gsub("CRC normal", "normal", sample.info.array$Sample)
sample.info.array$Location <- gsub("normal", "Normal", sample.info.array$Location)
sample.info.array$Location <- gsub("adenoma", "Adenoma", sample.info.array$Location)
table(sample.info.array$Location)
sample.info.array$Location <- factor(sample.info.array$Location,
    levels = c("Normal", "Adenoma", "CRC")
)

mycolor <- c("#208a42", "#d51f26", "#272e6a")
names(mycolor) <- c("Normal", "Adenoma", "CRC")

## 1.2. methylation beta value ----
beta.mat <- data.table::fread("data/GSE48684_series_matrix.txt",
    sep = "\t",
    header = TRUE, blank.lines.skip = FALSE
)
beta.mat <- as.data.frame(beta.mat)
rownames(beta.mat) <- beta.mat$ID_REF
beta.mat$ID_REF <- NULL
beta.mat[1:5, 1:5]

beta.mat$GSM1183439 %>% quantile(na.rm = TRUE)
quantile(beta.mat, na.rm = TRUE)
hist(beta.mat$GSM1235158, breaks = 100)

identical(colnames(beta.mat), rownames(sample.info.array))

pdf("DensityBin.pdf", 5, 10)
densityBeanPlot(as.matrix(beta.mat[, 1:10]),
    sampGroups = sample.info.array$Sample[1:10],
    sampNames = sample.info.array$ID[1:10]
)
dev.off()

pdf("MDS.pdf", 7, 7)
mdsPlot(as.matrix(beta.mat),
    numPositions = 2000,
    sampGroups = sample.info.array$Sample,
    sampNames = sample.info.array$ID
)
dev.off()

## 1.3. probe info ----
probe.info <- read.csv("data/HumanMethylation450_15017482_v1-2.csv",
    header = TRUE
)
rownames(probe.info) <- probe.info$IlmnID
probe.info <- probe.info[!is.na(probe.info$MAPINFO), ]
table(rownames(beta.mat) %in% probe.info$IlmnID)

beta.mat <- beta.mat[rownames(beta.mat) %in% probe.info$IlmnID, ]

quantile(nchar(probe.info$SourceSeq)) # 50

probe.info <- probe.info[rownames(beta.mat), ] %>%
    dplyr::select(c(
        "IlmnID", "CHR", "MAPINFO", "Strand",
        "UCSC_RefGene_Name", "UCSC_RefGene_Group", "UCSC_CpG_Islands_Name"
    ))
probe.info[1:5, ]

probe.info$Strand <- ifelse(probe.info$Strand == "F", "+", "-")
probe.info$Start <- probe.info$MAPINFO
probe.info$End <- probe.info$MAPINFO
probe.info$CHR <- paste0("chr", probe.info$CHR)

probe.bed <- GRanges(probe.info)

probe.bed.hg38 <- liftOver(
    probe.bed,
    import.chain("E:/LabWork/genome/hg19ToHg38.over.chain")
) %>% unlist()

length(probe.bed.hg38) - length(probe.bed) # loss 168
end(probe.bed.hg38) <- start(probe.bed.hg38) + 50

beta.mat <- beta.mat[rownames(beta.mat) %in% names(probe.bed.hg38), ]
probe.bed.hg38 <- probe.bed.hg38[rownames(beta.mat), ]
rm(probe.bed, probe.info)

save(beta.mat, sample.info.array, probe.bed.hg38, file = "data/methylation.rda")
