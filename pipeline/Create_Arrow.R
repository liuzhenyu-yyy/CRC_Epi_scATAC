args <- commandArgs(trailingOnly = TRUE)

ref <- args[1]
files <- args[2]
names <- args[3]

library(ArchR)

addArchRGenome(ref)

names(files) <- names

ArrowFiles <- createArrowFiles(
    inputFiles = files,
    sampleNames = names,
    minTSS = 4,
    minFrags = 3000,
    maxFrags = 1e+09,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE
)