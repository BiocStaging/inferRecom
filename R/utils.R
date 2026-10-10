## ========================== SHARED INTERNAL HELPERS ===========================

## ---------------- Argument checks ----------------
## Each stops with a message naming the argument, taken from the call, e.g.
## .checkFlag(rsOnly) gives "'rsOnly' must be a single logical value"

.checkString <- function(x, allowNull = FALSE, name = deparse(substitute(x))) {
    if (allowNull && is.null(x)) {
        return(invisible(TRUE))
    }
    if (!is.character(x) || length(x) != 1L || is.na(x)) {
        stop(sprintf(
            "'%s' must be a single character string%s",
            name, if (allowNull) " or NULL" else ""
        ), call. = FALSE)
    }
    invisible(TRUE)
}

.checkFlag <- function(x, name = deparse(substitute(x))) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
        stop(sprintf("'%s' must be a single logical value", name),
            call. = FALSE
        )
    }
    invisible(TRUE)
}

.checkCharacterOrNull <- function(x, name = deparse(substitute(x))) {
    if (!is.null(x) && !is.character(x)) {
        stop(sprintf("'%s' must be a character vector or NULL", name),
            call. = FALSE
        )
    }
    invisible(TRUE)
}

## A single number greater than zero, or at least zero if allowZero is TRUE.
## type names the expected value in the message ("numeric value", "integer")
.checkNumber <- function(x, allowZero = FALSE, allowNull = FALSE,
                         type = "numeric value",
                         name = deparse(substitute(x))) {
    if (allowNull && is.null(x)) {
        return(invisible(TRUE))
    }
    if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
        x < 0 || (!allowZero && x == 0)) {
        stop(sprintf(
            "'%s' must be a %s %s%s", name,
            if (allowZero) "non-negative" else "positive", type,
            if (allowNull) " or NULL" else ""
        ), call. = FALSE)
    }
    invisible(TRUE)
}

## ---------------- Data input and output ----------------

## Remove a .bed/.bim/.fam extension so plinkFile is the PLINK prefix. Other
## extensions are kept because PLINK prefixes may contain dots (e.g. "cohort.v2")
.stripPlinkExtension <- function(plinkFile) {
    if (tolower(file_ext(plinkFile)) %in% c("bed", "bim", "fam")) {
        plinkFile <- file_path_sans_ext(plinkFile)
    }
    plinkFile
}

## Read PLINK binary files. Returns a list with the genotypes as a numeric
## matrix (SNPs in rows, samples in columns; 0/1/2 copies of allele 2, NA if
## missing), plus the SNP map and FAM data from read.plink(). If rsOnly is
## TRUE, only SNPs with an "rs" prefix are kept.
.readPlinkGenotypes <- function(plinkFile, rsOnly = FALSE) {
    plinkFile <- .stripPlinkExtension(plinkFile)

    searchFiles <- sprintf("%s.%s", plinkFile, c("bed", "bim", "fam"))
    missingFiles <- searchFiles[!file.exists(searchFiles)]

    if (length(missingFiles) > 0L) {
        stop("PLINK file(s) not found: ", paste(missingFiles, collapse = ", "))
    }

    pedData <- read.plink(plinkFile)
    geno <- t(as(pedData$genotypes, "numeric"))

    if (isTRUE(rsOnly)) {
        rsIdx <- grep("^rs", rownames(geno), ignore.case = TRUE)
        if (length(rsIdx) == 0L) {
            stop("No SNPs with 'rs' prefix found")
        }
        geno <- geno[rsIdx, , drop = FALSE]
    }

    list(geno = geno, map = pedData$map, fam = pedData$fam)
}

## Read a tab-delimited genetic map and return its chr, pos and cM columns
.readGeneticMap <- function(mapFile) {
    if (!file.exists(mapFile)) {
        stop("Genetic map file not found: ", mapFile)
    }

    mapDf <- read.delim(mapFile, stringsAsFactors = FALSE)

    if (!all(c("chr", "pos", "cM") %in% colnames(mapDf))) {
        stop("mapFile must have columns: chr, pos, cM")
    }

    mapDf[, c("chr", "pos", "cM")]
}

## Convert a results DataFrame with startPos and finishPos columns to GRanges
## on the given chromosome. If caseControl is TRUE, split into a GRangesList
## with "case" and "control" elements, using the affected status in famData of
## the samples named in column idCol.
.asGRangesResult <- function(resultDf, chromosome, caseControl, famData,
                             idCol) {
    if (nrow(resultDf) == 0L) {
        return(.emptyGRangesResult(caseControl))
    }

    resultDf$chromosome <- .ucscChromosome(chromosome)

    if (caseControl) {
        affectedIds <- famData$member[famData$affected == 1]
        resultDf$caseStatus <- ifelse(resultDf[[idCol]] %in% affectedIds,
            "case",
            "control"
        )
    }

    grResult <- makeGRangesFromDataFrame(
        df = resultDf,
        seqnames.field = "chromosome",
        start.field = "startPos",
        end.field = "finishPos",
        keep.extra.columns = TRUE
    )

    if (caseControl) {
        grResult <- split(grResult, grResult$caseStatus)
    }

    grResult
}

## UCSC-style chromosome names ("chr4", "chrX", "chrM") from PLINK or UCSC
## names. PLINK codes chromosomes X, Y and MT as 23, 24 and 26.
.ucscChromosome <- function(chromosome) {
    chromosome <- sub("^chr", "", as.character(chromosome), ignore.case = TRUE)
    plinkCodes <- c("23" = "X", "24" = "Y", "26" = "M", "MT" = "M")
    recode <- chromosome %in% names(plinkCodes)
    chromosome[recode] <- plinkCodes[chromosome[recode]]
    paste0("chr", chromosome, recycle0 = TRUE)
}

## Empty result in the same form as .asGRangesResult()
.emptyGRangesResult <- function(caseControl) {
    if (caseControl) {
        return(GRangesList(case = GRanges(), control = GRanges()))
    }
    GRanges()
}
