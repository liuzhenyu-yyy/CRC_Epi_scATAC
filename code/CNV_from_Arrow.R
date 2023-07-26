args <- commandArgs(trailingOnly = TRUE)

# functions ----
library(ArchR)
library(Matrix)
library(SummarizedExperiment)
library(matrixStats)
library(readr)
library(GenomicRanges)
library(magrittr)
library(edgeR)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(1)

"%ni%" <- Negate("%in%")

countInsertions <- function(query, fragments, by = "RG") {
    # Count By Fragments Insertions
    inserts <- c(
        GRanges(
            seqnames = seqnames(fragments),
            ranges = IRanges(
                start(fragments),
                start(fragments)
            ),
            RG = mcols(fragments)[, by]
        ),
        GRanges(
            seqnames = seqnames(fragments),
            ranges = IRanges(
                end(fragments),
                end(fragments)
            ),
            RG = mcols(fragments)[, by]
        )
    )
    by <- "RG"
    overlapDF <- DataFrame(findOverlaps(query, inserts,
        ignore.strand = TRUE,
        maxgap = -1L,
        minoverlap = 0L,
        type = "any"
    ))
    overlapDF$name <- mcols(inserts)[overlapDF[, 2], by]
    overlapTDF <- transform(overlapDF, id = match(name, unique(name)))
    # Calculate Overlap Stats
    inPeaks <- table(overlapDF$name)
    total <- table(mcols(inserts)[, by])
    total <- total[names(inPeaks)]
    frip <- inPeaks / total
    # Summarize
    sparseM <- Matrix::sparseMatrix(
        i = overlapTDF[, 1],
        j = overlapTDF[, 4],
        x = rep(1, nrow(overlapTDF)),
        dims = c(length(query), length(unique(overlapDF$name)))
    )
    colnames(sparseM) <- unique(overlapDF$name)
    total <- total[colnames(sparseM)]
    frip <- frip[colnames(sparseM)]
    out <- list(counts = sparseM, frip = frip, total = total)
    return(out)
}

makeWindows <- function(genome, blacklist, windowSize = 10e6, slidingSize = 2e6) {
    chromSizes <- GRanges(names(seqlengths(genome)), IRanges(1, seqlengths(genome)))
    chromSizes <- GenomeInfoDb::keepStandardChromosomes(chromSizes, pruning.mode = "coarse")
    windows <- slidingWindows(x = chromSizes, width = windowSize, step = slidingSize) %>%
        unlist() %>%
        .[which(width(.) == windowSize), ]
    mcols(windows)$wSeq <- as.character(seqnames(windows))
    mcols(windows)$wStart <- start(windows)
    mcols(windows)$wEnd <- end(windows)
    message("Subtracting Blacklist...")
    windowsBL <- lapply(seq_along(windows), function(x) {
        if (x %% 100 == 0) {
            message(sprintf("%s of %s", x, length(windows)))
        }
        gr <- GenomicRanges::setdiff(windows[x, ], blacklist)
        mcols(gr) <- mcols(windows[x, ])
        return(gr)
    })
    names(windowsBL) <- paste0("w", seq_along(windowsBL))
    windowsBL <- unlist(GenomicRangesList(windowsBL), use.names = TRUE)
    mcols(windowsBL)$name <- names(windowsBL)
    message("Adding Nucleotide Information...")
    windowSplit <- split(windowsBL, as.character(seqnames(windowsBL)))
    windowNuc <- lapply(seq_along(windowSplit), function(x) {
        message(sprintf("%s of %s", x, length(windowSplit)))
        chrSeq <- Biostrings::getSeq(
            genome,
            chromSizes[which(as.character(seqnames(chromSizes)) == names(windowSplit)[x])]
        )
        grx <- windowSplit[[x]]
        aFreq <- alphabetFrequency(Biostrings::Views(chrSeq[[1]], ranges(grx)))
        mcols(grx)$GC <- rowSums(aFreq[, c("G", "C")]) / rowSums(aFreq)
        mcols(grx)$AT <- rowSums(aFreq[, c("A", "T")]) / rowSums(aFreq)
        return(grx)
    }) %>%
        GenomicRangesList() %>%
        unlist() %>%
        sortSeqlevels() %>%
        sort()
    windowNuc$N <- 1 - (windowNuc$GC + windowNuc$AT)
    windowNuc
}

scCNA <- function(windows, fragments, neighbors = 100, LFC = 1.5, FDR = 0.1, force = FALSE, remove = c("chrM", "chrX", "chrY")) {
    # Keep only regions in filtered chromosomes
    windows <- GenomeInfoDb::keepStandardChromosomes(windows, pruning.mode = "coarse")
    fragments <- GenomeInfoDb::keepStandardChromosomes(fragments, pruning.mode = "coarse")
    windows <- windows[seqnames(windows) %ni% remove]
    fragments <- fragments[seqnames(fragments) %ni% remove]

    # Count Insertions in windows
    message("Getting Counts...")
    counts <- countInsertions(windows, fragments, by = "RG")[[1]]
    message("Summarizing...")
    windowSummary <- GenomicRangesList()
    countSummary <- matrix(nrow = length(unique(windows$name)), ncol = ncol(counts))
    for (x in seq_along(unique(mcols(windows)$name))) {
        if (x %% 100 == 0) {
            message(sprintf("%s of %s", x, length(unique(mcols(windows)$name))))
        }
        idx <- which(mcols(windows)$name == unique(mcols(windows)$name)[x])
        wx <- windows[idx, ]
        wo <- GRanges(mcols(wx)$wSeq, ranges = IRanges(mcols(wx)$wStart, mcols(wx)$wEnd))[1, ]
        mcols(wo)$name <- mcols(wx)$name[1]
        mcols(wo)$effectiveLength <- sum(width(wx))
        mcols(wo)$percentEffectiveLength <- 100 * sum(width(wx)) / width(wo)
        mcols(wo)$GC <- sum(mcols(wx)$GC * width(wx)) / width(wo)
        mcols(wo)$AT <- sum(mcols(wx)$AT * width(wx)) / width(wo)
        mcols(wo)$N <- sum(mcols(wx)$N * width(wx)) / width(wo)
        countSummary[x, ] <- Matrix::colSums(counts[idx, , drop = FALSE])
        windowSummary[[x]] <- wo
    }
    windowSummary <- unlist(windowSummary)

    # Keep only regions with less than 0.1% N
    keep <- which(windowSummary$N < 0.001)
    windowSummary <- windowSummary[keep, ]
    countSummary <- countSummary[keep, ]

    # Now determine the nearest neighbors by GC content
    message("Computing Background...")
    bdgMean <- matrix(nrow = nrow(countSummary), ncol = ncol(countSummary))
    bdgSd <- matrix(nrow = nrow(countSummary), ncol = ncol(countSummary))
    log2FC <- matrix(nrow = nrow(countSummary), ncol = ncol(countSummary))
    z <- matrix(nrow = nrow(countSummary), ncol = ncol(countSummary))
    pval <- matrix(nrow = nrow(countSummary), ncol = ncol(countSummary))

    for (x in seq_len(nrow(countSummary))) {
        if (x %% 100 == 0) {
            message(sprintf("%s of %s", x, nrow(countSummary)))
        }
        # Get Nearest Indices
        idxNN <- head(order(abs(windowSummary$GC[x] - windowSummary$GC)), neighbors + 1)
        idxNN <- idxNN[idxNN %ni% x]
        # Background
        if (any(colMeans(countSummary[idxNN, ]) == 0)) {
            if (force) {
                message("Warning! Background Mean = 0 Try a higher neighbor count or remove cells with 0 in colMins")
            } else {
                stop("Background Mean = 0!")
            }
        }
        bdgMean[x, ] <- colMeans(countSummary[idxNN, ])
        bdgSd[x, ] <- matrixStats::colSds(countSummary[idxNN, ])
        log2FC[x, ] <- log2((countSummary[x, ] + 1e-5) / (bdgMean[x, ] + 1e-5))
        z[x, ] <- (countSummary[x, ] - bdgMean[x, ]) / bdgSd[x, ]
        pval[x, ] <- 2 * pnorm(-abs(z[x, ]))
    }
    padj <- apply(pval, 2, function(x) p.adjust(x, method = "fdr"))
    CNA <- matrix(0, nrow = nrow(countSummary), ncol = ncol(countSummary))
    CNA[which(log2FC >= LFC & padj <= FDR)] <- 1

    se <- SummarizedExperiment(
        assays = SimpleList(
            CNA = CNA,
            counts = countSummary,
            log2FC = log2FC,
            padj = padj,
            pval = pval,
            z = z,
            bdgMean = bdgMean,
            bdgSd = bdgSd
        ),
        rowRanges = windowSummary
    )
    colnames(se) <- colnames(counts)

    return(se)
}

# get fragements from arrow ----
file_arrow <- args[1]
keep <- readRDS(args[2])
prefix <- args[3]

fragments <- getFragmentsFromArrow(
    ArrowFile = file_arrow,
    cellNames = keep
)

blacklist <- import.bed("/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/database/hg38/hg38-blacklist.v2.bed")
windows <- makeWindows(
    genome = BSgenome.Hsapiens.UCSC.hg38,
    blacklist = blacklist
)

cnaObj <- scCNA(windows, fragments,
    neighbors = 100, LFC = 1.5, FDR = 0.1,
    force = FALSE, remove = c("chrM", "chrX", "chrY")
)

out_cna <- paste(prefix, "_CNV_LFC_GC.rds", sep = "")
saveRDS(cnaObj, out_cna)