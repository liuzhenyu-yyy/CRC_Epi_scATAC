##########################################################
###
### Generate fragments file from arrow files
### Input: arrow file (for ArchR)
### Output: fragments file
### Author: Zhenyu Liu
###
##########################################################

setwd("E:/LabWork/Project/CRC_NGS_ATAC/CRC_Epi_scATAC/Results/Export")

source("../../code/00.Requirements.R")
rm(foo)

arrows <- dir("E:/LabWork/Project/CRC_NGS_ATAC/CRC_Epi_scATAC/Project_Dir_All/ArrowFiles", full.names = TRUE)
names(arrows) <- basename(arrows) %>% gsub("\\.arrow$|-nofacs", "", .)

patient.info <- read.table("patient.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
patient.rename <- patient.info$Patient
names(patient.rename) <- patient.info$Old_ID
options(scipen = 999)

for (one in names(arrows)[29]) {
    message(Sys.time(), " : ", "Extracting fragment for ", one)
    gr <- getFragmentsFromArrow(
        ArrowFile = arrows[one]
    )
    # gr$RG %>% as.character()

    gr <- sortSeqlevels(gr)
    gr <- sort(gr)

    df <- data.frame(
        chrom = seqnames(gr),
        chromStart = start(gr) - 1,
        chromEnd = end(gr),
        barcode = gr$RG %>% as.character() %>% gsub("^.+?#", "", .) %>% gsub(one, patient.rename[one], .),
        readSupport = 1
    )

    message(Sys.time(), " : ", "Writing into files for ", one)
    data.table::fwrite(df,
        file = paste0("./", patient.rename[one], ".sort.fragments.tsv.gz"),
        quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
    )
    rm(gr, df, one)
    gc()
}
