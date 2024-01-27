setwd("E:/LabWork/Project/CRC_NGS_ATAC/CRC_Epi_scATAC/Results/Export")

source("../../code/00.Requirements.R")
rm(foo)

arrows <- dir("E:/LabWork/Project/CRC_NGS_ATAC/CRC_Epi_scATAC/Project_Dir_All/ArrowFiles", full.names = TRUE)
names(arrows) <- basename(arrows) %>% gsub("\\.arrow$|-nofacs", "", .)

patient.info <- read.table("patient.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
patient.rename <- patient.info$Patient
names(patient.rename) <- patient.info$Old_ID

for (one in names(arrows)) {
    message(Sys.time(), " : ", "Extracting fragment for ", one)
    gr <- getFragmentsFromArrow(
        ArrowFile = arrows[one]
    )
    gr$RG %>% as.character()

    df <- data.frame(
        chrom = seqnames(gr),
        chromStart = start(gr) - 1,
        chromEnd = end(gr),
        barcode = gr$RG %>% as.character() %>% gsub("^.+?#", "", .) %>% gsub(one, patient.rename[one], .),
        readSupport = 1
    )

    message(Sys.time(), " : ", "Writing into files for ", one)
    write.table(df,
        file = paste0("F:/temp/", patient.rename[one], "_fragments.tsv"),
        quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
    )
}
