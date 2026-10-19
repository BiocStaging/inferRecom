#' Detect Runs of Homozygosity in PLINK Genotype Data
#'
#' This function identifies runs of homozygosity (ROH) in PLINK genotype data
#' based on minimum physical length and SNP count thresholds. It can optionally
#' filter for specific samples, apply a minimum genetic distance criterion, and
#' separate results by case/control status.
#'
#' @param plinkFile Character string; path to the PLINK prefix (without
#'     extensions).
#' @param mapFile Character string or NULL; path to a tab-delimited genetic map
#'     file containing columns: \code{chr}, \code{pos}, \code{cM}. If NULL, cM
#'     filtering is not available. Default is NULL.
#' @param minMb Numeric; minimum physical length of ROH in megabases.
#'     Default is 1.
#' @param minSnps Integer; minimum number of consecutive homozygous SNPs
#'     required for an ROH. Default is 5.
#' @param minCm Numeric or NULL; minimum genetic length of ROH in centiMorgans.
#'     If NULL, no genetic length filtering is applied. Requires mapFile to be
#'     specified. Default is NULL.
#' @param rsOnly Logical; whether to restrict analysis to SNPs with names
#'     beginning with "rs". Default is FALSE.
#' @param sampleIds Character vector or NULL; specific sample IDs to analyze. If
#'     NULL, all samples are analyzed. Default is NULL.
#' @param caseControl Logical; if TRUE, function returns a
#'     \link[GenomicRanges:GRangesList-class]{GRangesList} of case and
#'     control ROH, using the affected status in the FAM file with PLINK
#'     coding: 2 = affected (case); 1 = unaffected or missing (control).
#'     Default is FALSE.
#' @param BPPARAM A BiocParallelParam object specifying parallel execution.
#'     Default is \code{SerialParam()} for serial execution. Use
#'     \code{MulticoreParam(workers = n)} for parallel processing across n
#'     cores.
#'
#' @return A \link[GenomicRanges:GRanges-class]{GRanges} object containing
#'     detected ROH with metadata columns: \code{sampleId}, \code{startSnp},
#'     \code{finishSnp}, \code{numSnps}. If genetic map is provided, also
#'     includes \code{startCm}, \code{finishCm}. If
#'     \code{caseControl = TRUE}, returns a
#'     \link[GenomicRanges:GRangesList-class]{GRangesList} with
#'     two elements: \code{case} and \code{control}. Returns an empty
#'     \link[GenomicRanges:GRanges-class]{GRanges} if no ROH are detected.
#'     Chromosome names use the UCSC style (e.g. \code{chr4}).
#'
#' @details
#' Runs of homozygosity (ROH) are continuous genomic segments where an
#' individual is homozygous at all marker positions. ROH can indicate
#' autozygosity (inheritance of identical haplotypes from a common ancestor)
#' and are used to estimate inbreeding, identify disease-associated loci,
#' and study population history.
#'
#' This function identifies ROH by:
#' \enumerate{
#'   \item Converting genotypes to binary (1 = homozygous,
#'         0 = heterozygous/missing)
#'   \item Identifying runs using run-length encoding
#'   \item Filtering runs by minimum SNP count and physical length
#'   \item Optionally filtering by minimum genetic length (cM)
#' }
#'
#' Genotypes from \pkg{snpStats} are coded as:
#' \itemize{
#'   \item 0 = homozygous reference (AA) - counted as homozygous
#'   \item 1 = heterozygous (AB) - breaks ROH
#'   \item 2 = homozygous alternate (BB) - counted as homozygous
#'   \item NA = missing - breaks ROH
#' }
#'
#' Filtering strategy:
#' \itemize{
#'   \item Physical length (minMb): Always applied
#'   \item SNP count (minSnps): Always applied
#'   \item Genetic length (minCm): Applied only if mapFile is provided and
#'         minCm is not NULL
#' }
#'
#' Parallel processing is performed by individual, improving efficiency for
#' datasets with many samples.
#'
#' @examples
#' # Serial execution (default)
#'
#' # Get path to example data
#' dataPath <- system.file("extdata", package = "inferRecom")
#' plinkFile <- file.path(dataPath, "simCEU")
#'
#' # Detect ROH
#' rohData <- hzRun(
#'     plinkFile = plinkFile,
#'     minMb = 1,
#'     minSnps = 50
#' )
#'
#' # View results
#' head(rohData)
#'
#' # Analyze specific individuals
#' sampleIds <- c("F3C1", "F4C3", "F5P2")
#' rohSubset <- hzRun(
#'     plinkFile = plinkFile,
#'     minMb = .5,
#'     minSnps = 50,
#'     sampleIds = sampleIds
#' )
#'
#' rohSubset
#'
#' # Separate by case/control status
#' rohCc <- hzRun(
#'     plinkFile = plinkFile,
#'     minMb = .5,
#'     minSnps = 50,
#'     caseControl = TRUE
#' )
#'
#' rohCases <- rohCc$case
#' rohControls <- rohCc$control
#'
#' rohCases
#' rohControls
#'
#' @export
hzRun <- function(plinkFile,
                  mapFile = NULL,
                  minMb = 1,
                  minSnps = 5,
                  minCm = NULL,
                  rsOnly = FALSE,
                  sampleIds = NULL,
                  caseControl = FALSE,
                  BPPARAM = SerialParam()) {
    ## Argument validation
    .checkString(plinkFile)
    .checkString(mapFile, allowNull = TRUE)
    .checkNumber(minMb)
    .checkNumber(minSnps, type = "integer")
    .checkNumber(minCm, allowNull = TRUE)
    .checkFlag(rsOnly)
    .checkCharacterOrNull(sampleIds)
    .checkFlag(caseControl)

    if (!is.null(minCm) && is.null(mapFile)) {
        stop("'mapFile' must be provided when 'minCm' is specified")
    }

    minSnps <- as.integer(minSnps)

    ## Read PLINK data
    plinkData <- .readPlinkGenotypes(plinkFile, rsOnly)
    genoMatrix <- plinkData$geno
    mapData <- plinkData$map
    famData <- plinkData$fam

    ## Read genetic map if provided
    geneticMap <- NULL
    if (!is.null(mapFile)) {
        geneticMap <- .readGeneticMap(mapFile)
    }

    ## Filter for specific IDs if provided
    if (!is.null(sampleIds)) {
        sampleIdx <- which(colnames(genoMatrix) %in% sampleIds)
        if (length(sampleIdx) == 0L) {
            stop("None of the specified IDs found in genotype data")
        }
        genoMatrix <- genoMatrix[, sampleIdx, drop = FALSE]
    }

    ## Get physical positions for SNPs
    snpNames <- rownames(genoMatrix)
    locationVec <- .getSnpLocations(snpNames, mapData)

    ## Convert threshold to base pairs
    minBp <- minMb * 1e6

    ## Detect ROH for all samples (with parallel support)
    rohList <- .detectRohAllSamples(
        genoMatrix, locationVec, snpNames,
        minSnps, minBp, geneticMap, minCm,
        BPPARAM
    )

    ## Return empty GRanges or GRangesList if no ROH were detected
    if (length(rohList) == 0L) {
        return(.emptyGRangesResult(caseControl))
    }

    ## Build DataFrame with the columns common to all outputs
    rohMat <- do.call(rbind, rohList)
    resultDf <- DataFrame(
        sampleId = rohMat[, 1],
        startSnp = rohMat[, 2],
        finishSnp = rohMat[, 3],
        numSnps = as.integer(rohMat[, 4]),
        startPos = as.numeric(rohMat[, 5]),
        finishPos = as.numeric(rohMat[, 6])
    )

    ## Add genetic positions if a map was provided
    if (!is.null(geneticMap)) {
        resultDf$startCm <- as.numeric(rohMat[, 7])
        resultDf$finishCm <- as.numeric(rohMat[, 8])
    }

    ## Convert to GRanges, split by case/control status if requested
    .asGRangesResult(resultDf, mapData$chromosome[1], caseControl, famData,
        idCol = "sampleId"
    )
}


## ========================== INTERNAL HELPERS ===========================

#' Get Physical Positions for SNPs
#'
#' Map SNP names to physical positions
#'
#' @param snpNames Character vector of SNP names
#' @param mapData Data frame with SNP map information
#'
#' @return Numeric vector of positions
#' @keywords internal
.getSnpLocations <- function(snpNames, mapData) {
    ## Create lookup table for matching
    lookupTable <- setNames(mapData$position, mapData$snp.name)

    positions <- lookupTable[snpNames]

    positions[is.na(positions)] <- 0

    as.numeric(positions)
}

#' Classify Genotype as Homozygous or Heterozygous
#'
#' Converts numeric genotype codes to binary (1 = homozygous, 0 = heterozygous)
#'
#' @param genoVec Numeric vector of genotype codes (0=AA, 1=AB, 2=BB,
#'     NA=missing)
#'
#' @return Integer vector: 1 for homozygous (0 or 2), 0 for heterozygous (1)
#'     or missing (NA)
#' @keywords internal
.isHomozygous <- function(genoVec) {
    ## Homozygous: 0 (AA) or 2 (BB)
    ## Heterozygous: 1 (AB)
    ## Missing: NA (treated as 0 to break runs)
    resultVec <- as.integer(genoVec == 0 | genoVec == 2)
    resultVec[is.na(resultVec)] <- 0L
    resultVec
}

#' Detect ROH for a Single Sample
#'
#' Identifies runs of homozygosity in one individual
#'
#' @param genoVec Numeric vector of genotypes for one sample (0, 1, 2, or NA)
#' @param locationVec Physical positions of SNPs in base pairs
#' @param snpNames SNP identifiers
#' @param sampleId Sample identifier
#' @param minSnps Minimum number of SNPs
#' @param minBp Minimum physical length in base pairs
#' @param geneticMap Data frame with genetic map (chr, pos, cM), or NULL
#' @param minCm Minimum genetic length in cM, or NULL
#'
#' @return Matrix with ROH information, or NULL if no ROH detected
#' @keywords internal
.detectRohSample <- function(genoVec, locationVec, snpNames, sampleId,
                             minSnps, minBp, geneticMap = NULL, minCm = NULL) {
    ## Runs of consecutive homozygous SNPs (heterozygous or missing SNPs end
    ## a run), with their first and last SNP and physical positions
    runs <- rle(.isHomozygous(genoVec))
    endIdx <- cumsum(runs$lengths)
    roh <- data.frame(start = endIdx - runs$lengths + 1L, end = endIdx)
    roh <- roh[runs$values == 1L, , drop = FALSE]
    roh$startBp <- locationVec[roh$start]
    roh$endBp <- locationVec[roh$end]

    ## Filter by number of SNPs and physical length
    roh <- roh[roh$end - roh$start + 1L >= minSnps &
        roh$endBp - roh$startBp >= minBp, , drop = FALSE]

    ## Add genetic positions, and filter by genetic length if requested
    if (!is.null(geneticMap)) {
        roh$startCm <- .getGeneticPosition(roh$startBp, geneticMap)
        roh$endCm <- .getGeneticPosition(roh$endBp, geneticMap)

        if (!is.null(minCm)) {
            lengthCm <- roh$endCm - roh$startCm
            roh <- roh[!is.na(lengthCm) & lengthCm >= minCm, , drop = FALSE]
        }
    }

    if (nrow(roh) == 0L) {
        return(NULL)
    }

    ## Character matrix with one row per ROH; cM columns only with a map
    resultMat <- cbind(
        sampleId, snpNames[roh$start], snpNames[roh$end],
        as.character(roh$end - roh$start + 1L),
        as.character(roh$startBp), as.character(roh$endBp)
    )
    if (!is.null(geneticMap)) {
        resultMat <- cbind(
            resultMat, as.character(roh$startCm), as.character(roh$endCm)
        )
    }

    unname(resultMat)
}

#' Process Single Individual for ROH Detection
#'
#' Wrapper function for parallel processing of a single individual
#'
#' @param sampleIdx Integer index of sample in genotype matrix
#' @param genoMatrix Matrix of genotypes (SNPs x samples)
#' @param locationVec Physical positions of SNPs
#' @param snpNames SNP identifiers
#' @param sampleIds Sample identifiers
#' @param minSnps Minimum number of SNPs
#' @param minBp Minimum physical length
#' @param geneticMap Genetic map data frame or NULL
#' @param minCm Minimum genetic length or NULL
#'
#' @return Matrix with ROH information, or NULL
#' @keywords internal
.processSingleIndividual <- function(sampleIdx, genoMatrix, locationVec,
                                     snpNames, sampleIds, minSnps, minBp,
                                     geneticMap, minCm) {
    .detectRohSample(
        genoVec = genoMatrix[, sampleIdx],
        locationVec = locationVec,
        snpNames = snpNames,
        sampleId = sampleIds[sampleIdx],
        minSnps = minSnps,
        minBp = minBp,
        geneticMap = geneticMap,
        minCm = minCm
    )
}

#' Detect ROH Across All Samples
#'
#' Efficiently processes all samples to detect runs of homozygosity
#' with optional parallel processing
#'
#' @param genoMatrix Matrix of genotypes (SNPs x samples)
#' @param locationVec Physical positions of SNPs in base pairs
#' @param snpNames SNP identifiers
#' @param minSnps Minimum number of SNPs
#' @param minBp Minimum physical length in base pairs
#' @param geneticMap Data frame with genetic map, or NULL
#' @param minCm Minimum genetic length in cM, or NULL
#' @param BPPARAM BiocParallelParam object for parallel execution
#'
#' @return List of matrices, one per sample with detected ROH
#' @keywords internal
.detectRohAllSamples <- function(genoMatrix, locationVec, snpNames,
                                 minSnps, minBp, geneticMap = NULL,
                                 minCm = NULL, BPPARAM = SerialParam()) {
    numSamples <- ncol(genoMatrix)
    sampleIds <- colnames(genoMatrix)

    ## Process all samples in parallel
    rohList <- bplapply(seq_len(numSamples), function(i) {
        .processSingleIndividual(
            sampleIdx = i,
            genoMatrix = genoMatrix,
            locationVec = locationVec,
            snpNames = snpNames,
            sampleIds = sampleIds,
            minSnps = minSnps,
            minBp = minBp,
            geneticMap = geneticMap,
            minCm = minCm
        )
    }, BPPARAM = BPPARAM)

    ## Remove NULL results
    rohList <- rohList[!vapply(rohList, is.null, logical(1))]

    rohList
}
