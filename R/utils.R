## ========================== SHARED INTERNAL HELPERS ===========================

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
