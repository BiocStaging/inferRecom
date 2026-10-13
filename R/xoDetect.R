#' Detect Meiotic Recombinations from PLINK Genotype Data
#'
#' This function identifies recombination events in family-based genotype data,
#' using PLINK-formatted input files and a genetic map. It supports 2- or
#' 3-child family structures, and can optionally write results to disk if an
#' output path is provided.
#'
#' @param plinkFile Character string; path to the PLINK prefix (without
#'     extensions).
#' @param mapFile Character string; path to a tab-delimited map file containing
#'     columns: \code{chr}, \code{pos}, \code{rate}, \code{cM}.
#' @param familySize Integer; either 2 or 3, specifying the family pedigree size.
#' @param parent Character; either \code{"mother"} or \code{"father"},
#'     specifying the parent to trace crossovers from.
#' @param rsOnly Logical; whether to restrict to SNPs with names beginning
#'     with "rs". Default is \code{FALSE}.
#' @param snpFilter Integer; minimum number of informative SNPs separating
#'     two putative crossovers in the same child. A crossover that is
#'     reversed within fewer SNPs, typically because of genotyping error, is
#'     removed together with its reversal. Default is 5.
#' @param cmFilter Numeric; minimum genetic distance (cM) separating putative
#'     crossovers. Default is 1.
#' @param caseControl Logical; if TRUE, function returns a
#'    \link[GenomicRanges:GRangesList-class]{GRangesList} of case and control
#'    crossovers for 3+ child families. Default is \code{FALSE}.
#' @param out Character or NULL; if provided, output is written as a CSV at this
#'     path.
#' @param BPPARAM A BiocParallelParam object specifying parallel execution.
#'     Default is \code{SerialParam()} for serial execution. Use
#'     \code{MulticoreParam(workers = n)} for parallel processing across n cores.
#'
#' @return A \link[GenomicRanges:GRanges-class]{GRanges} listing detected
#'     crossover intervals with metadata columns: \code{childId} (in 3-child
#'     families), \code{familyId}, \code{startSnp}, \code{finishSnp},
#'     \code{startPos}, \code{finishPos}, \code{startCm}, \code{finishCm}.
#'     Chromosome names use the UCSC style (e.g. \code{chr4}).
#'
#' @details
#' This function reads PLINK genotype data (via \pkg{snpStats}) and infers
#' loci where meiotic recombination has occurred, resolved to the nearest
#' informative SNPs. The algorithm identifies informative marker configurations
#' where parents are heterozygous/homozygous and filters for Mendelian
#' consistency before detecting inheritance state changes that indicate
#' crossovers.
#'
#' For 3-child families, children are analyzed in triples to allow
#' identification of the specific recombinant child. For 2-child families, state
#' changes between the two siblings are used to infer crossover locations, but
#' the individual cannot be determined.
#'
#' If \code{out} is provided, the resulting table is written to disk as
#' UTF-8 encoded CSV without row names.
#'
#' Parallel processing is performed by family, improving efficiency for datasets
#' with many families.
#'
#' @examples
#' # Serial execution (default)
#'
#' # Load example data
#' dataPath <- system.file("extdata", package = "inferRecom")
#' plinkFile <- file.path(dataPath, "simCEU")
#'
#' # Basic maternal crossover detection
#' mapFemale <- file.path(dataPath, "female_chr4.txt")
#' xoMat3 <- xoDetect(
#'     plinkFile = plinkFile,
#'     mapFile = mapFemale,
#'     familySize = 3,
#'     parent = "mother"
#' )
#'
#' # View results
#' xoMat3
#'
#' # Basic 2-child family
#' xoMat2 <- xoDetect(
#'     plinkFile = plinkFile,
#'     mapFile = mapFemale,
#'     familySize = 2,
#'     parent = "mother"
#' )
#'
#' # View results
#' xoMat2
#'
#' # Maternal crossover detection with case/control separation
#' mapFemale <- file.path(dataPath, "female_chr4.txt")
#' xoMat3CC <- xoDetect(
#'     plinkFile = plinkFile,
#'     mapFile = mapFemale,
#'     familySize = 3,
#'     parent = "mother",
#'     caseControl = TRUE
#' )
#'
#' # View results
#' xoMat3CC$case
#' xoMat3CC$control
#'
#' @export
xoDetect <- function(plinkFile,
                     mapFile,
                     familySize = c(2, 3),
                     parent = c("mother", "father"),
                     rsOnly = FALSE,
                     snpFilter = 5,
                     cmFilter = 1,
                     caseControl = FALSE,
                     out = NULL,
                     BPPARAM = SerialParam()) {
    ## Argument validation
    parent <- match.arg(parent)
    familySize <- as.integer(match.arg(as.character(familySize), c("2", "3")))
    .checkString(plinkFile)
    .checkString(mapFile)
    .checkFlag(rsOnly)
    .checkNumber(snpFilter, allowZero = TRUE, type = "integer")
    .checkNumber(cmFilter, allowZero = TRUE)
    .checkFlag(caseControl)
    .checkString(out, allowNull = TRUE)

    if (familySize == 2L && caseControl) {
        stop("No case-control separation for 2-child families")
    }

    ## Read PLINK data
    pedData <- .readPlinkGenotypes(plinkFile, rsOnly)
    mapDf <- .readGeneticMap(mapFile)

    ## Recode genotypes as 0 (missing), 1, 2, 3 for crossover detection
    genoMatrix <- pedData$geno + 1
    genoMatrix[is.na(genoMatrix)] <- 0

    famData <- pedData$fam

    ## Detect crossovers in all families of the requested size
    resDf <- .xoDetectFamilies(
        pedData, genoMatrix, famData, mapDf, parent,
        snpFilter, cmFilter, familySize, BPPARAM
    )

    ## Write output if requested
    if (!is.null(out) && nrow(resDf) > 0L) {
        write.csv(as.data.frame(resDf),
            file = out,
            row.names = FALSE, fileEncoding = "UTF-8"
        )
        message("Results written to: ", out)
    }

    ## Convert to GRanges, split by case/control status if requested
    .asGRangesResult(resDf, mapDf$chr[1], caseControl, famData,
        idCol = "childId"
    )
}


## ========================== INTERNAL HELPERS ===========================

## Identify informative SNP types (heterozygous/homozygous parent pairs)
.getInformativeSNPs <- function(genoSubset) {
    ## Remove missing data (genotype == 0)
    genoNoMissing <- genoSubset[apply(genoSubset, 1, function(row) {
        all(row != 0)
    }), , drop = FALSE]

    if (nrow(genoNoMissing) == 0L) {
        return(NULL)
    }

    ## Informative SNP configurations: (1,2), (2,1), (3,2), (2,3)
    infSnpTypes <- list(c(1, 2), c(2, 1), c(3, 2), c(2, 3))

    infIndex <- vapply(seq_len(nrow(genoNoMissing)), function(i) {
        parentGeno <- as.vector(genoNoMissing[i, seq_len(2)])
        any(vapply(infSnpTypes, function(x) {
            identical(x, parentGeno)
        }, logical(1)))
    }, logical(1))

    if (!any(infIndex)) {
        return(NULL)
    }

    genoInformative <- genoNoMissing[infIndex, , drop = FALSE]

    if (nrow(genoInformative) == 0L || is.vector(genoInformative)) {
        return(NULL)
    }

    ## Filter for Mendelian consistency
    mendelConsistent <- vapply(seq_len(nrow(genoInformative)), function(i) {
        all(genoInformative[i, 3:ncol(genoInformative)] %in%
            as.vector(genoInformative[i, seq_len(2)]))
    }, logical(1))

    if (!any(mendelConsistent)) {
        return(NULL)
    }

    genoFiltered <- genoInformative[mendelConsistent, , drop = FALSE]

    if (nrow(genoFiltered) == 0L || is.vector(genoFiltered)) {
        return(NULL)
    }

    genoFiltered
}

## Create binary inheritance matrix for parent-specific haplotypes: 1 where a
## child's genotype matches the mother's (column 2), 0 otherwise
.makeBinaryInheritance <- function(parentSubset) {
    if (nrow(parentSubset) <= 1L) {
        return(NULL)
    }

    binaryMatrix <- (parentSubset[, -(1:2), drop = FALSE] == parentSubset[, 2]) * 1
    dimnames(binaryMatrix) <- list(rownames(parentSubset), NULL)
    binaryMatrix
}

## SNP filter: remove runs of an inheritance state that are shorter than
## snpFilter SNPs and have the same state on both sides. Such a run is a
## crossover followed by its reversal in the same child within too few SNPs,
## typically caused by genotyping error; both state changes are removed.
.smoothShortRuns <- function(state, snpFilter) {
    runs <- rle(unname(state))
    nRuns <- length(runs$lengths)

    if (nRuns < 3L) {
        return(state)
    }

    inner <- seq(2L, nRuns - 1L)
    short <- inner[runs$lengths[inner] < snpFilter &
        runs$values[inner - 1L] == runs$values[inner + 1L]]
    runs$values[short] <- runs$values[short - 1L]

    smoothed <- inverse.rle(runs)
    names(smoothed) <- names(state)
    smoothed
}

## Apply cM filter to remove closely-spaced crossovers
.applyCMFilter <- function(xoData, cmFilter) {
    if (is.na(cmFilter) || nrow(xoData) <= 1L) {
        return(xoData)
    }

    ## Compare crossovers within each child, or within each family when the
    ## child is unknown (2-child families), sorted by position
    groupCol <- if ("childId" %in% colnames(xoData)) "childId" else "familyId"
    xoSorted <- xoData[order(xoData[[groupCol]], xoData$startPos), ]

    ## Drop both crossovers of each consecutive pair in the same group that
    ## are less than cmFilter cM apart
    nXo <- nrow(xoSorted)
    group <- xoSorted[[groupCol]]
    tooClose <- group[-1L] == group[-nXo] &
        xoSorted$startCm[-1L] - xoSorted$finishCm[-nXo] < cmFilter
    pairStart <- which(tooClose)
    dropIndex <- unique(c(pairStart, pairStart + 1L))

    if (length(dropIndex) > 0L) {
        xoSorted <- xoSorted[-dropIndex, ]
    }

    xoSorted
}

## Add physical and genetic positions to crossover results
.addPositions <- function(xoData, pedMap, geneticMap) {
    if (nrow(xoData) == 0L) {
        return(xoData)
    }

    startPosList <- numeric(nrow(xoData))
    finishPosList <- numeric(nrow(xoData))
    startCmList <- numeric(nrow(xoData))
    finishCmList <- numeric(nrow(xoData))

    for (rowIdx in seq_len(nrow(xoData))) {
        ## Physical positions
        startSnpIdx <- which(pedMap$snp.name == xoData$startSnp[rowIdx])
        finishSnpIdx <- which(pedMap$snp.name == xoData$finishSnp[rowIdx])

        if (length(startSnpIdx) > 0L) {
            startPosList[rowIdx] <- pedMap$position[startSnpIdx[1]]
        } else {
            startPosList[rowIdx] <- NA_real_
        }

        if (length(finishSnpIdx) > 0L) {
            finishPosList[rowIdx] <- pedMap$position[finishSnpIdx[1]]
        } else {
            finishPosList[rowIdx] <- NA_real_
        }

        ## Genetic positions (cM)
        if (!is.na(startPosList[rowIdx])) {
            startCmList[rowIdx] <- geneticMap$cM[which.min(abs(startPosList[rowIdx] -
                geneticMap$pos))]
        } else {
            startCmList[rowIdx] <- NA_real_
        }

        if (!is.na(finishPosList[rowIdx])) {
            finishCmList[rowIdx] <- geneticMap$cM[which.min(abs(finishPosList[rowIdx] -
                geneticMap$pos))]
        } else {
            finishCmList[rowIdx] <- NA_real_
        }
    }

    xoData$startPos <- startPosList
    xoData$finishPos <- finishPosList
    xoData$startCm <- startCmList
    xoData$finishCm <- finishCmList

    xoData
}


## ---------------- Process single family (for parallel execution) ------------

.processSingleFamily <- function(currentFamId, famSubset, genoMatrix, parent,
                                 snpFilter, familySize) {
    ## Genotypes ordered father, mother, children
    family <- .orderedFamily(currentFamId, famSubset, genoMatrix,
        minMembers = if (familySize == 2L) 4L else 5L
    )
    if (is.null(family)) {
        return(NULL)
    }

    ## Get informative SNPs
    informativeGeno <- .getInformativeSNPs(family$geno)
    if (is.null(informativeGeno)) {
        return(NULL)
    }

    ## Select parent-specific informative SNPs
    parentColIdx <- if (parent == "father") 1L else 2L
    parentSpecificGeno <- informativeGeno[informativeGeno[, parentColIdx] == 2, , drop = FALSE]

    ## Create binary inheritance matrix
    binaryInheritance <- .makeBinaryInheritance(parentSpecificGeno)
    if (is.null(binaryInheritance)) {
        return(NULL)
    }

    if (familySize == 2L) {
        return(.processTwoChildFamily(currentFamId, binaryInheritance, snpFilter))
    }

    .processThreeChildFamily(
        currentFamId, binaryInheritance, family$fam, snpFilter
    )
}


## ---------------- Two-child family processing ----------------

.processTwoChildFamily <- function(currentFamId, binaryInheritance, snpFilter) {
    if (nrow(binaryInheritance) < 2L) {
        return(NULL)
    }

    ## Inheritance state: 1 if the children share the parental allele, else 2
    inheritanceState <- ifelse(
        binaryInheritance[, 1] == binaryInheritance[, 2], 1, 2
    )
    names(inheritanceState) <- rownames(binaryInheritance)

    ## Apply SNP filter
    inheritanceState <- .smoothShortRuns(inheritanceState, snpFilter)

    ## Each state change is a crossover in one of the two children
    changeIdx <- which(diff(inheritanceState) != 0)

    if (length(changeIdx) == 0L) {
        return(NULL)
    }

    snpNames <- names(inheritanceState)
    cbind(
        familyId = as.character(currentFamId),
        startSnp = snpNames[changeIdx],
        finishSnp = snpNames[changeIdx + 1L]
    )
}


## ---------------- Three-child family processing ----------------

## Inheritance state of a triple of children at a SNP, indexed by the 0/1
## pattern of the three children read as a 3-bit number (plus 1). A pattern
## and its complement describe the same state. The codes are chosen so that
## the size of a change between two states identifies the one child that
## recombined (see .getChildIdx()).
.tripleStateCodes <- c(1, 3, 14, 7, 7, 14, 3, 1)

## Child (1 to 3 within the triple) whose crossover produces a change of this
## size between .tripleStateCodes, or NA if the change does not identify one
.getChildIdx <- function(stateChange) {
    childByChange <- c(
        "2" = 3L, "7" = 3L, "4" = 2L, "13" = 2L, "6" = 1L, "11" = 1L
    )
    unname(childByChange[as.character(abs(stateChange))])
}

.processThreeChildFamily <- function(currentFamId, binaryInheritance,
                                     currentFamilyOrdered, snpFilter) {
    numChildren <- ncol(binaryInheritance)

    ## Analyse children in triples; an incomplete final triple is completed
    ## with the preceding children
    tripleEvents <- lapply(seq(1L, numChildren, by = 3L), function(start) {
        tripleCols <- seq(min(start, numChildren - 2L), length.out = 3L)
        .tripleCrossovers(
            binaryInheritance[, tripleCols, drop = FALSE],
            currentFamilyOrdered$member[tripleCols + 2L],
            currentFamId, snpFilter
        )
    })

    do.call(rbind, tripleEvents)
}

## Crossovers in one triple of children: a character matrix with one row per
## crossover, or NULL if there are none
.tripleCrossovers <- function(tripleGeno, tripleChildIds, currentFamId,
                              snpFilter) {
    ## Too few informative SNPs to call crossovers
    if (nrow(tripleGeno) < 10L) {
        return(NULL)
    }

    ## Inheritance state at each SNP, with 0 for an unrecognised pattern
    inheritanceState <- .tripleStateCodes[drop(tripleGeno %*% c(4, 2, 1)) + 1]
    inheritanceState[is.na(inheritanceState)] <- 0
    names(inheritanceState) <- rownames(tripleGeno)

    ## Apply SNP filter
    inheritanceState <- .smoothShortRuns(inheritanceState, snpFilter)

    ## Each state change that identifies a child is a crossover in that child
    stateChange <- diff(inheritanceState)
    changeIdx <- which(stateChange != 0)
    childIdx <- .getChildIdx(stateChange[changeIdx])
    changeIdx <- changeIdx[!is.na(childIdx)]
    childIdx <- childIdx[!is.na(childIdx)]

    if (length(changeIdx) == 0L) {
        return(NULL)
    }

    snpNames <- names(inheritanceState)
    cbind(
        childId = as.character(tripleChildIds[childIdx]),
        familyId = as.character(currentFamId),
        startSnp = snpNames[changeIdx],
        finishSnp = snpNames[changeIdx + 1L]
    )
}


## ---------------- All families of one size (with parallel support) ----------

## Detect crossovers in all families with exactly two children (familySize 2,
## 4 members) or with three or more children (familySize 3, 5+ members)
.xoDetectFamilies <- function(pedData, genoMatrix, famData, geneticMap, parent,
                              snpFilter, cmFilter, familySize, BPPARAM) {
    ## Identify families of the requested size (2 parents + children)
    famCounts <- table(famData[, 1])
    famKeep <- if (familySize == 2L) {
        names(famCounts)[famCounts == 4]
    } else {
        names(famCounts)[famCounts >= 5]
    }
    famSubset <- famData[famData[, 1] %in% famKeep, , drop = FALSE]
    famIds <- unique(famSubset[, 1])
    pedMap <- pedData$map

    ## Process families in parallel
    xoList <- bplapply(famIds, function(currentFamId) {
        .processSingleFamily(currentFamId, famSubset, genoMatrix, parent,
            snpFilter,
            familySize = familySize
        )
    }, BPPARAM = BPPARAM)

    ## Remove NULL results
    xoList <- xoList[!vapply(xoList, is.null, logical(1))]

    if (length(xoList) == 0L) {
        return(.xoResultDataFrame(data.frame(), familySize))
    }

    xoAllEvents <- do.call(rbind, xoList)
    xoAllEvents <- as.data.frame(xoAllEvents, stringsAsFactors = FALSE)

    ## Add positions, remove duplicates and apply cM filter
    xoAllEvents <- .addPositions(xoAllEvents, pedMap, geneticMap)
    xoAllEvents <- distinct(xoAllEvents)
    xoAllEvents <- .applyCMFilter(xoAllEvents, cmFilter)

    .xoResultDataFrame(xoAllEvents, familySize)
}

## Crossover results as a DataFrame, one row per event. The recombinant child
## (childId) can only be identified in families with three or more children.
.xoResultDataFrame <- function(xoEvents, familySize) {
    resultDf <- DataFrame(
        familyId = as.character(xoEvents$familyId),
        startSnp = as.character(xoEvents$startSnp),
        finishSnp = as.character(xoEvents$finishSnp),
        startPos = as.numeric(xoEvents$startPos),
        finishPos = as.numeric(xoEvents$finishPos),
        startCm = as.numeric(xoEvents$startCm),
        finishCm = as.numeric(xoEvents$finishCm)
    )

    if (familySize == 3L) {
        resultDf <- cbind(
            DataFrame(childId = as.character(xoEvents$childId)), resultDf
        )
    }

    resultDf
}
