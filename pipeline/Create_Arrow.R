args <- commandArgs(trailingOnly = TRUE)

ref <- args[1]
files <- args[2]

library(ArchR)

addArchRGenome(ref)

names <- gsub(".fragments.gz", "", files, fixed = TRUE)
names <- gsub(".sort", "", names, fixed = TRUE)
names <- gsub(".flank", "", names, fixed = TRUE)
names <- gsub(".flank", "", names, fixed = TRUE)
names <- gsub("_20200917_fragments.tsv.gz.tbi.gz", "", names, fixed = TRUE)

names(files) <- names

ArrowFiles <- createArrowFiles(
    inputFiles = files,
    sampleNames = names(files),
    minTSS = 0,
    minFrags = 0,
    maxFrags = 1e+09,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE
)