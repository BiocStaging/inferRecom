#' Phase Haplotypes Using Crossover Information
#'
#' This function phases parental and child haplotypes based on detected
#' crossover events from maternal and paternal meioses.
#'
#' @param plinkFile Character string; path to the PLINK prefix (without
#'     extensions).
#' @param xoDetectPaternal \link[GenomicRanges:GRanges-class]{GRanges}; output
#'     from \code{xoDetect()} for paternal crossovers (parent = "father").
#' @param xoDetectMaternal \link[GenomicRanges:GRanges-class]{GRanges}; output
#'     from \code{xoDetect()} for maternal crossovers (parent = "mother").
#' @param famIds Character vector or NULL; specific family IDs to phase. If
#'     NULL, all families with 5+ members are analyzed. Default is NULL.
#' @param rsOnly Logical; whether to restrict to SNPs with names beginning
#'     with "rs". Default is TRUE.
#' @param BPPARAM A BiocParallelParam object specifying parallel execution.
#'     Default is \code{SerialParam()} for serial execution. Use
#'     \code{MulticoreParam(workers = n)} for parallel processing across n
#'     cores.
#'
#' @return A \link[VariantAnnotation:VCF-class]{VCF} object with one row per
#'     SNP and one column per family member, for all phased families. The
#'     \code{GT} genotype field holds phased genotypes (e.g. \code{"0|1"}),
#'     coded against \code{REF} (PLINK allele 1) and \code{ALT} (PLINK
#'     allele 2); \code{"."} marks an allele that could not be phased. For
#'     children the first allele is paternal and the second maternal; for
#'     parents the alleles are the two phased haplotypes. \code{colData}
#'     holds \code{familyId}, \code{fatherId}, \code{motherId}, \code{sex}
#'     and \code{affected} from the FAM file. Families without phase
#'     information are omitted with a message. Chromosome names use the UCSC
#'     style (e.g. \code{chr4}). Use \code{VariantAnnotation::writeVcf()} to
#'     write the result to a file.
#'
#' @details
#' This function performs haplotype phasing using family-based genetic data
#' and detected crossover events. The phasing process:
#' \enumerate{
#'   \item Identifies informative SNPs where parents are heterozygous/homozygous
#'   \item Assigns alleles (1 or 2) to children based on inheritance
#'   \item Converts allele codes to nucleotides (A, C, G, T)
#'   \item Imputes homozygous sites from the other parent
#'   \item Phases parental haplotypes using crossover breakpoints
#'   \item Adds back non-informative SNPs with inferred genotypes where possible
#'   \item Outputs diploid genotypes for all family members across all SNPs
#' }
#'
#' Non-informative SNPs are handled as follows. In every case, a child
#' whose genotype is missing or inconsistent with its parents' genotypes is
#' left unphased.
#' \itemize{
#'   \item Both parents homozygous: each parent passes on its only allele.
#'   \item Both parents heterozygous: phase uncertain, set to unknown.
#'   \item One parent heterozygous, one homozygous (informative SNPs left out
#'         because of a missing or inconsistent genotype): each child's alleles
#'         follow from the genotypes. The heterozygous parent is phased from
#'         which of its haplotypes the children inherited at the nearest
#'         phased SNPs on both sides, and left unknown if that evidence is
#'         missing or conflicting.
#' }
#'
#' The function requires crossover detection results from both parents as
#' GRanges objects. Run \code{xoDetect()} separately for maternal and paternal
#' crossovers before using this function.
#'
#' Parallel processing is performed by family, improving efficiency for datasets
#' with many families.
#'
#' @examples
#' # Load example data
#' dataPath <- system.file("extdata", package = "inferRecom")
#' plinkFile <- file.path(dataPath, "simCEU")
#' mapFemale <- file.path(dataPath, "female_chr4.txt")
#' mapMale <- file.path(dataPath, "male_chr4.txt")
#'
#' # Detect crossovers for both parents (returns GRanges objects)
#' xoPat <- xoDetect(
#'     plinkFile = plinkFile,
#'     mapFile = mapMale,
#'     familySize = 3,
#'     parent = "father"
#' )
#' xoMat <- xoDetect(
#'     plinkFile = plinkFile,
#'     mapFile = mapFemale,
#'     familySize = 3,
#'     parent = "mother"
#' )
#'
#' # Phase haplotypes (serial execution, default)
#' phased <- xoPhase(
#'     plinkFile = plinkFile,
#'     xoDetectPaternal = xoPat,
#'     xoDetectMaternal = xoMat
#' )
#'
#' phased
#'
#' # Phased genotypes for the members of family F6
#' f6 <- colData(phased)$familyId == "F6"
#' head(geno(phased)$GT[, f6])
#'
#' # Phase haplotypes for subset of families
#' phasedSubset <- xoPhase(
#'     plinkFile = plinkFile,
#'     xoDetectPaternal = xoPat,
#'     xoDetectMaternal = xoMat,
#'     famIds = c("F4", "F6")
#' )
#'
#' phasedSubset
#'
#' # Write to a VCF file
#' vcfFile <- tempfile(fileext = ".vcf")
#' writeVcf(phasedSubset, vcfFile)
#'
#' @export
xoPhase <- function(plinkFile,
    xoDetectPaternal,
    xoDetectMaternal,
    famIds = NULL,
    rsOnly = TRUE,
    BPPARAM = SerialParam()) {
    ## Argument validation
    .checkString(plinkFile)
    .checkCrossovers(xoDetectPaternal)
    .checkCrossovers(xoDetectMaternal)
    .checkCharacterOrNull(famIds)
    .checkFlag(rsOnly)

    ## Read PLINK data, and crossovers as data frames
    plinkData <- .readPlinkGenotypes(plinkFile, rsOnly)
    bim <- .bimFromMap(plinkData$map)
    xoPat <- as.data.frame(xoDetectPaternal)
    xoMat <- as.data.frame(xoDetectMaternal)

    ## Phase each family using parallel processing
    familyIds <- .familiesToPhase(plinkData$fam, famIds)
    phaseList <- bplapply(familyIds, function(id) {
        .phaseSingleFamily(
            familyId = id,
            geno = plinkData$geno,
            fam = plinkData$fam,
            bim = bim,
            xoPat = xoPat,
            xoMat = xoMat
        )
    }, BPPARAM = BPPARAM)
    names(phaseList) <- familyIds

    ## Drop families without phase information
    noPhase <- vapply(phaseList, is.null, logical(1))
    if (any(noPhase)) {
        message(
            "No phase information for families: ",
            paste(familyIds[noPhase], collapse = ", ")
        )
    }

    .phaseListToVcf(phaseList[!noPhase], bim, plinkData$fam)
}


## ========================== INTERNAL HELPERS ===========================

#' Check Crossover Input
#'
#' Checks that a crossover object is a GRanges on a single chromosome, and
#' warns if it is empty.
#'
#' @param x Crossover results from \code{xoDetect()}
#' @param name Argument name used in messages
#'
#' @return Invisibly TRUE; stops with an error otherwise
#' @keywords internal
.checkCrossovers <- function(x, name = deparse(substitute(x))) {
    if (!is(x, "GRanges")) {
        stop(sprintf("'%s' must be a GRanges object", name), call. = FALSE)
    }

    if (!.checkSameChromosome(x)) {
        stop("GRanges object can only contain ranges from a single chromosome",
            call. = FALSE
        )
    }

    if (length(x) == 0L) {
        warning(sprintf("Empty GRanges object provided for '%s'", name),
            call. = FALSE
        )
    }

    invisible(TRUE)
}

#' SNP Information from the PLINK Map
#'
#' @param map SNP map from \code{read.plink()}
#'
#' @return Data frame with columns chr, snp, cm, pos, a1 and a2, in BIM file
#'     column order
#' @keywords internal
.bimFromMap <- function(map) {
    data.frame(
        chr = map$chromosome,
        snp = map$snp.name,
        cm = map$cM,
        pos = map$position,
        a1 = map$allele.1,
        a2 = map$allele.2,
        row.names = NULL
    )
}

#' Families to Phase
#'
#' @param fam FAM data frame
#' @param famIds Character vector of family IDs, or NULL
#'
#' @return famIds if given, otherwise the IDs of all families with 5 or more
#'     members (two parents and three or more children), in file order
#' @keywords internal
.familiesToPhase <- function(fam, famIds) {
    if (!is.null(famIds)) {
        return(famIds)
    }

    famCounts <- table(fam[, 1])
    unique(fam[fam[, 1] %in% names(famCounts)[famCounts >= 5], 1])
}

#' Check if all chromosomes in GRanges are the same
#'
#' @param gr GRanges object
#'
#' @return Logical; TRUE if all chromosomes are identical
#' @keywords internal
.checkSameChromosome <- function(gr) {
    if (length(gr) == 0) {
        return(TRUE)
    }

    chroms <- as.character(seqnames(gr))
    length(unique(chroms)) == 1
}

#' Phase a Single Family
#'
#' @param familyId Family identifier
#' @param geno Numeric genotype matrix
#' @param fam FAM data frame
#' @param bim BIM data frame
#' @param xoPat Paternal crossover data frame
#' @param xoMat Maternal crossover data frame
#'
#' @return DataFrame with phased haplotypes, or NULL if the family cannot be
#'     phased
#' @keywords internal
.phaseSingleFamily <- function(familyId, geno, fam, bim, xoPat, xoMat) {
    ## Genotypes ordered father, mother, children
    family <- .orderedFamily(familyId, fam, geno, minMembers = 3L)
    if (is.null(family)) {
        return(NULL)
    }
    genoSubset <- family$geno
    fatherId <- colnames(genoSubset)[1]
    motherId <- colnames(genoSubset)[2]

    ## Get informative genotypes
    genoAll <- .getInformativeGenotypes(genoSubset)

    if (is.null(genoAll)) {
        return(NULL)
    }

    ## Alleles each parent passed to each child, as nucleotides
    alleles <- .parentAlleleTables(genoAll, bim)

    ## Phase children, and parents using crossover information
    phaseChildren <- .phaseChildren(alleles$pat, alleles$mat)
    childSnps <- as.character(phaseChildren$rsID)
    patPhase <- .phaseParent(alleles$pat, xoPat, bim, fatherId)
    matPhase <- .phaseParent(alleles$mat, xoMat, bim, motherId)

    ## Combine phased informative SNPs, in the children's SNP order
    famPhaseInformative <- cbind(
        phaseChildren,
        .alignToSnps(patPhase, childSnps),
        .alignToSnps(matPhase, childSnps)
    )

    ## Add non-informative SNPs back to the output
    famPhaseFull <- .addNonInformativeSNPs(
        phasedData = famPhaseInformative,
        genoSubset = genoSubset,
        bim = bim,
        fatherId = fatherId,
        motherId = motherId,
        fatherCol = 1L,
        motherCol = 2L
    )

    DataFrame(famPhaseFull)
}

#' Alleles Passed on by Each Parent
#'
#' @param genoAll Informative genotypes, ordered father, mother, children
#' @param bim BIM data frame
#'
#' @return List with \code{pat} and \code{mat}: for each parent, the
#'     nucleotide each child received at the SNPs where that parent is
#'     heterozygous, plus the SNPs where it is homozygous (marked
#'     \code{hom = 0})
#' @keywords internal
.parentAlleleTables <- function(genoAll, bim) {
    patActg <- .convertToNucleotides(.assignAlleleCodes(genoAll,
        parent = "father", fatherCol = 1L, motherCol = 2L
    ), bim)
    matActg <- .convertToNucleotides(.assignAlleleCodes(genoAll,
        parent = "mother", fatherCol = 1L, motherCol = 2L
    ), bim)

    list(
        pat = .addHomozygousSNPs(patActg, matActg, genoAll, bim, "father"),
        mat = .addHomozygousSNPs(matActg, patActg, genoAll, bim, "mother")
    )
}

#' Align Phased Haplotypes to SNPs
#'
#' @param phase Data frame of phased haplotypes with SNP IDs as row names
#' @param snps SNP IDs in the required order
#'
#' @return The rows of phase for snps, in that order, without row names
#' @keywords internal
.alignToSnps <- function(phase, snps) {
    aligned <- phase[match(snps, rownames(phase)), , drop = FALSE]
    rownames(aligned) <- NULL
    aligned
}

#' Add Non-Informative SNPs to Phased Output
#'
#' Adds SNPs that were not informative for phasing back into the output,
#' with alleles inferred from the parental genotypes where possible (see
#' \code{.phaseBothHomozygous()} and \code{.phaseOneHeterozygous()}).
#'
#' @param phasedData Data frame with phased informative SNPs
#' @param genoSubset Full genotype matrix for the family
#' @param bim BIM data frame with SNP information
#' @param fatherId Father identifier
#' @param motherId Mother identifier
#' @param fatherCol Column index for father in genoSubset
#' @param motherCol Column index for mother in genoSubset
#'
#' @return Data frame with all SNPs (informative and non-informative)
#' @keywords internal
.addNonInformativeSNPs <- function(phasedData, genoSubset, bim,
    fatherId, motherId,
    fatherCol, motherCol) {
    ## Identify non-informative SNPs
    nonInfSnps <- setdiff(rownames(genoSubset), as.character(phasedData$rsID))

    if (length(nonInfSnps) == 0) {
        ## All SNPs were informative, return as is
        return(phasedData)
    }

    ## Family member IDs; child IDs from the Pat/Mat haplotype columns
    parentCols <- paste0(rep(c(fatherId, motherId), each = 2), c("_1", "_2"))
    childCols <- setdiff(
        colnames(phasedData), c("rsID", "location", parentCols)
    )
    ids <- list(
        father = fatherId, mother = motherId,
        children = unique(gsub("(Pat|Mat)$", "", childCols))
    )

    ## Phase each non-informative SNP; one character row per SNP
    nonInfMatrix <- t(vapply(nonInfSnps, function(snpId) {
        .phaseNonInformativeSnp(
            snpId, genoSubset, bim, ids, fatherCol, motherCol, phasedData
        )
    }, character(ncol(phasedData)), USE.NAMES = FALSE))
    colnames(nonInfMatrix) <- colnames(phasedData)
    nonInfDf <- as.data.frame(nonInfMatrix, stringsAsFactors = FALSE)

    ## Combine informative and non-informative SNPs, sorted by chromosome
    ## then position
    combinedDf <- rbind(phasedData, nonInfDf)
    combinedDf$location <- as.numeric(combinedDf$location)
    snpChr <- bim$chr[match(as.character(combinedDf$rsID), bim$snp)]
    combinedDf <- combinedDf[order(snpChr, combinedDf$location), ]
    rownames(combinedDf) <- NULL

    combinedDf
}

#' Phase One Non-Informative SNP
#'
#' @param snpId SNP identifier
#' @param genoSubset Full genotype matrix for the family
#' @param bim BIM data frame with SNP information
#' @param ids List with \code{father}, \code{mother} and \code{children} IDs
#' @param fatherCol,motherCol Column indices for the parents in genoSubset
#' @param phasedData Data frame with phased informative SNPs
#'
#' @return Named character vector with one value per column of phasedData:
#'     the SNP ID, its location and the phased alleles, NA where unknown. All
#'     NA for a SNP that is not in the BIM file.
#' @keywords internal
.phaseNonInformativeSnp <- function(snpId, genoSubset, bim, ids, fatherCol,
                                    motherCol, phasedData) {
    row <- setNames(
        rep(NA_character_, ncol(phasedData)), colnames(phasedData)
    )
    bimIdx <- match(snpId, bim$snp)

    if (is.na(bimIdx)) {
        return(row)
    }

    row["rsID"] <- snpId
    row["location"] <- as.character(bim$pos[bimIdx])

    parentGeno <- genoSubset[snpId, c(fatherCol, motherCol)]
    childGeno <- genoSubset[snpId, ids$children]
    alleles <- c(bim$a1[bimIdx], bim$a2[bimIdx])

    ## Infer alleles from the parental genotypes; leave NA when a parent's
    ## genotype is missing or both parents are heterozygous
    hetParents <- parentGeno == 1
    phased <- if (anyNA(parentGeno) || all(hetParents)) {
        character()
    } else if (!any(hetParents)) {
        .phaseBothHomozygous(parentGeno, childGeno, alleles, ids)
    } else {
        .phaseOneHeterozygous(
            parentGeno, childGeno, alleles, ids,
            bim$pos[bimIdx], phasedData
        )
    }

    row[names(phased)] <- phased
    row
}

#' Phase a SNP Where Both Parents Are Homozygous
#'
#' Each parent passes on its only allele, so every child with the expected
#' genotype carries the father's allele on its paternal haplotype and the
#' mother's allele on its maternal haplotype. Children with any other
#' genotype, which indicates a genotyping error, are left as NA.
#'
#' @param parentGeno Father and mother genotypes (0 or 2)
#' @param childGeno Named vector of child genotypes
#' @param alleles Allele 1 and allele 2 of the SNP
#' @param ids List with \code{father}, \code{mother} and \code{children} IDs
#'
#' @return Named character vector of phased alleles
#' @keywords internal
.phaseBothHomozygous <- function(parentGeno, childGeno, alleles, ids) {
    fatherAllele <- alleles[parentGeno[[1]] / 2 + 1]
    motherAllele <- alleles[parentGeno[[2]] / 2 + 1]
    children <- ids$children[childGeno %in% (sum(parentGeno) / 2)]

    c(
        setNames(
            rep(fatherAllele, length(children)),
            paste0(children, "Pat", recycle0 = TRUE)
        ),
        setNames(
            rep(motherAllele, length(children)),
            paste0(children, "Mat", recycle0 = TRUE)
        ),
        setNames(rep(fatherAllele, 2), paste0(ids$father, c("_1", "_2"))),
        setNames(rep(motherAllele, 2), paste0(ids$mother, c("_1", "_2")))
    )
}

#' Phase a SNP Where One Parent Is Heterozygous
#'
#' The homozygous parent's haplotypes both carry its allele. The allele each
#' child received from the heterozygous parent follows from the genotypes
#' (see \code{.matchAlleleCode()}); children with a missing or inconsistent
#' genotype are left as NA. The heterozygous parent's phase follows from
#' which of its haplotypes each child inherited around the SNP (see
#' \code{.inheritedHaplotype()}): that haplotype carries the allele the child
#' received. The parent is left as NA if no child gives evidence or the
#' children disagree.
#'
#' @param parentGeno Father and mother genotypes, one of them 1
#' @param childGeno Named vector of child genotypes
#' @param alleles Allele 1 and allele 2 of the SNP
#' @param ids List with \code{father}, \code{mother} and \code{children} IDs
#' @param snpPos Physical position of the SNP
#' @param phasedData Data frame with phased informative SNPs
#'
#' @return Named character vector of phased alleles
#' @keywords internal
.phaseOneHeterozygous <- function(parentGeno, childGeno, alleles, ids,
    snpPos, phasedData) {
    hetParent <- if (parentGeno[[1]] == 1) "father" else "mother"
    hetId <- ids[[hetParent]]
    homId <- if (hetParent == "father") ids$mother else ids$father
    homGeno <- parentGeno[[if (hetParent == "father") 2 else 1]]
    homAllele <- alleles[homGeno / 2 + 1]
    hetSuffix <- if (hetParent == "father") "Pat" else "Mat"
    homSuffix <- if (hetParent == "father") "Mat" else "Pat"

    ## Allele each child received from the heterozygous parent
    received <- alleles[.matchAlleleCode(
        hetParent, parentGeno[[1]], parentGeno[[2]], childGeno
    )]
    names(received) <- ids$children
    children <- ids$children[!is.na(received)]

    ## Allele on the heterozygous parent's haplotype 1 implied by each child
    hap1 <- phasedData[[paste0(hetId, "_1")]]
    hap2 <- phasedData[[paste0(hetId, "_2")]]
    positions <- as.numeric(phasedData$location)
    impliedHap1 <- vapply(children, function(child) {
        inherited <- .inheritedHaplotype(
            phasedData[[paste0(child, hetSuffix)]], hap1, hap2,
            positions, snpPos
        )
        if (is.na(inherited)) {
            return(NA_character_)
        }
        allele <- received[[child]]
        if (inherited == 1L) allele else alleles[3L - match(allele, alleles)]
    }, character(1))
    impliedHap1 <- unique(impliedHap1[!is.na(impliedHap1)])

    phased <- c(
        setNames(
            received[children],
            paste0(children, hetSuffix, recycle0 = TRUE)
        ),
        setNames(
            rep(homAllele, length(children)),
            paste0(children, homSuffix, recycle0 = TRUE)
        ),
        setNames(rep(homAllele, 2), paste0(homId, c("_1", "_2")))
    )

    if (length(impliedHap1) == 1L) {
        hapAlleles <- c(impliedHap1, alleles[3L - match(impliedHap1, alleles)])
        phased <- c(phased, setNames(hapAlleles, paste0(hetId, c("_1", "_2"))))
    }

    phased
}

#' Haplotype a Child Inherited Around a SNP
#'
#' Which haplotype (1 or 2) of a heterozygous parent a child inherited around
#' a SNP, read from the nearest phased SNP on each side where the parent's
#' two haplotypes differ and the child's allele from that parent is known.
#'
#' @param childAlleles Child's alleles from the parent at the phased SNPs
#' @param hap1,hap2 Parent's haplotypes at the phased SNPs
#' @param positions Physical positions of the phased SNPs
#' @param snpPos Physical position of the SNP
#'
#' @return 1 or 2, or NA if the two sides disagree (the child has a crossover
#'     between them) or there is no usable SNP. Near a chromosome end, the
#'     one available side is used.
#' @keywords internal
.inheritedHaplotype <- function(childAlleles, hap1, hap2, positions, snpPos) {
    inherited <- ifelse(childAlleles == hap1, 1L,
        ifelse(childAlleles == hap2, 2L, NA_integer_)
    )
    usable <- !is.na(inherited) & hap1 != hap2
    left <- which(usable & positions < snpPos)
    right <- which(usable & positions > snpPos)

    sides <- unique(c(
        inherited[left][which.max(positions[left])],
        inherited[right][which.min(positions[right])]
    ))
    if (length(sides) == 1L) sides else NA_integer_
}

#' Get Informative Genotypes for Phasing
#'
#' Filters for informative SNPs and checks Mendelian consistency
#'
#' @param genoSubset Genotype matrix for one family
#'
#' @return Matrix of informative genotypes or NULL
#' @keywords internal
.getInformativeGenotypes <- function(genoSubset) {
    ## Remove SNPs with missing genotypes (NA; 0 is homozygous for allele 1)
    genoNoMissing <- genoSubset[apply(genoSubset, 1, function(row) {
        all(!is.na(row))
    }), , drop = FALSE]

    if (nrow(genoNoMissing) == 0L) {
        return(NULL)
    }

    ## Informative SNP configurations: (0,1), (1,0), (2,1), (1,2)
    ## In numeric coding: heterozygous/homozygous parent combinations
    infSnpTypes <- list(c(0, 1), c(1, 0), c(2, 1), c(1, 2))

    infIndex <- vapply(seq_len(nrow(genoNoMissing)), function(i) {
        parentGeno <- as.vector(genoNoMissing[i, seq_len(2)])
        any(vapply(infSnpTypes, function(x) {
            identical(x, parentGeno)
        }, logical(1)))
    }, logical(1))

    if (!any(infIndex)) {
        return(NULL)
    }

    genoInf <- genoNoMissing[infIndex, , drop = FALSE]

    if (nrow(genoInf) == 0L) {
        return(NULL)
    }

    ## Filter for Mendelian consistency
    mendel <- vapply(seq_len(nrow(genoInf)), function(i) {
        all(genoInf[i, 3:ncol(genoInf)] %in% genoInf[i, seq_len(2)])
    }, logical(1))

    if (!any(mendel)) {
        return(NULL)
    }

    genoAll <- genoInf[mendel, , drop = FALSE]

    if (nrow(genoAll) == 0L) {
        return(NULL)
    }

    genoAll
}

#' Assign Allele Codes for One Parent
#'
#' Determines which allele (1 or 2) was inherited from specified parent
#'
#' @param genoDf Informative genotype matrix
#' @param parent "father" or "mother"
#' @param fatherCol Column index for father
#' @param motherCol Column index for mother
#'
#' @return Data frame with allele codes
#' @keywords internal
.assignAlleleCodes <- function(genoDf, parent, fatherCol, motherCol) {
    parentCol <- if (parent == "father") fatherCol else motherCol

    ## Select SNPs where parent is heterozygous (genotype = 1)
    genoParent <- genoDf[genoDf[, parentCol] == 1, , drop = FALSE]

    if (nrow(genoParent) == 0L) {
        return(NULL)
    }

    ## Determine allele codes for children
    genoChildren <- genoParent[, -c(fatherCol, motherCol), drop = FALSE]
    codeMatrix <- .matchAlleleCode(
        parent, genoParent[, fatherCol], genoParent[, motherCol], genoChildren
    )

    ## Create data frame with SNP IDs
    result <- data.frame(
        rsID = rownames(genoParent),
        codeMatrix,
        stringsAsFactors = FALSE
    )
    colnames(result) <- c("rsID", colnames(genoChildren))
    rownames(result) <- rownames(genoParent)

    result
}

#' Match Allele Code
#'
#' Which allele (1 or 2) a child inherited from the heterozygous parent, given
#' the genotypes of both parents. The homozygous parent passes on allele 2
#' \code{homGt / 2} times, so the heterozygous parent passed on allele 2
#' \code{childGt - homGt / 2} times; the code is 1 more than this.
#'
#' @param parent "father" or "mother": the heterozygous parent
#' @param fatherGt Father genotype(s) (0, 1, or 2)
#' @param motherGt Mother genotype(s) (0, 1, or 2)
#' @param childGt Child genotype(s) (0, 1, or 2): a vector or a matrix with one
#'     row per SNP, recycled against \code{fatherGt} and \code{motherGt}
#'
#' @return Integer allele code(s) (1 or 2), with the dimensions of
#'     \code{childGt}; NA where the parent is not heterozygous, the other
#'     parent is not homozygous, or the child's genotype is inconsistent
#' @keywords internal
.matchAlleleCode <- function(parent, fatherGt, motherGt, childGt) {
    hetGt <- if (parent == "father") fatherGt else motherGt
    homGt <- if (parent == "father") motherGt else fatherGt

    code <- 1 + childGt - homGt / 2
    valid <- (hetGt == 1 & homGt %in% c(0, 2) & code %in% c(1, 2)) %in% TRUE
    code[!valid] <- NA

    storage.mode(code) <- "integer"
    code
}

#' Convert Allele Codes to Nucleotides
#'
#' @param alleleCodeDf Data frame with allele codes (1 or 2)
#' @param bim BIM data frame with allele information
#'
#' @return Data frame with nucleotide alleles: allele 1 or 2 from the BIM
#'     file, or "" for other codes and for SNPs that are not in the BIM file
#' @keywords internal
.convertToNucleotides <- function(alleleCodeDf, bim) {
    if (is.null(alleleCodeDf)) {
        return(NULL)
    }

    bimIdx <- match(alleleCodeDf$rsID, bim$snp)
    codes <- as.matrix(alleleCodeDf[, -1, drop = FALSE])

    ## Look up each SNP's alleles by row
    alleles <- matrix("", nrow(codes), ncol(codes))
    isA1 <- which(codes == 1 & !is.na(bimIdx))
    isA2 <- which(codes == 2 & !is.na(bimIdx))
    alleles[isA1] <- bim$a1[bimIdx[row(codes)[isA1]]]
    alleles[isA2] <- bim$a2[bimIdx[row(codes)[isA2]]]

    location <- as.character(bim$pos[bimIdx])
    location[is.na(bimIdx)] <- ""

    result <- data.frame(
        rsID = as.character(alleleCodeDf$rsID),
        location = location,
        alleles,
        stringsAsFactors = FALSE
    )
    colnames(result) <- c("rsID", "location", colnames(alleleCodeDf)[-1])
    rownames(result) <- rownames(alleleCodeDf)

    result
}

#' Add Homozygous SNPs from Other Parent
#'
#' @param actgMatrixHom Matrix for parent being imputed
#' @param actgMatrixHet Matrix from other parent (het SNPs)
#' @param genotypeDf Original genotype data frame
#' @param bim BIM data frame
#' @param parent "father" or "mother"
#'
#' @return Combined data frame with het and hom SNPs
#' @keywords internal
.addHomozygousSNPs <- function(actgMatrixHom, actgMatrixHet,
                               genotypeDf, bim, parent) {
    if (is.null(actgMatrixHom) || is.null(actgMatrixHet)) {
        return(actgMatrixHom)
    }

    pCol <- if (parent == "father") 1L else 2L
    nSamples <- ncol(actgMatrixHom) - 2L

    ## At the other parent's heterozygous SNPs this parent is homozygous, so
    ## every child received its only allele
    snps <- rownames(actgMatrixHet)
    bimIdx <- match(snps, bim$snp)
    genoIdx <- match(snps, rownames(genotypeDf))
    known <- !is.na(bimIdx) & !is.na(genoIdx)
    snps <- snps[known]
    bimIdx <- bimIdx[known]
    genoIdx <- genoIdx[known]

    actgMatrixHom$hom <- 1

    if (length(snps) == 0L) {
        return(actgMatrixHom)
    }

    homAllele <- ifelse(genotypeDf[genoIdx, pCol] == 0,
        bim$a1[bimIdx], bim$a2[bimIdx]
    )
    homDf <- data.frame(
        rsID = snps,
        location = as.character(bim$pos[bimIdx]),
        matrix(homAllele, nrow = length(snps), ncol = nSamples),
        hom = 0,
        stringsAsFactors = FALSE
    )
    colnames(homDf) <- colnames(actgMatrixHom)

    ## Combine and sort by location
    combinedDf <- rbind(actgMatrixHom, homDf)
    combinedDf$location <- as.numeric(combinedDf$location)
    combinedDf <- combinedDf[order(combinedDf$location), ]

    combinedDf
}

#' Phase Child Haplotypes
#'
#' @param actgCombPat Paternal allele data frame
#' @param actgCombMat Maternal allele data frame
#'
#' @return Data frame with phased children
#' @keywords internal
.phaseChildren <- function(actgCombPat, actgCombMat) {
    ## Ensure both data frames have the same SNPs
    ## Union of all SNPs
    allSnps <- unique(c(
        as.character(actgCombPat$rsID),
        as.character(actgCombMat$rsID)
    ))

    ## Create aligned data frames with all SNPs
    actgCombPatAligned <- .alignSnpData(actgCombPat, allSnps)
    actgCombMatAligned <- .alignSnpData(actgCombMat, allSnps)

    ## Child allele columns follow rsID and location; "hom" is last
    nChildren <- ncol(actgCombPatAligned) - 3L
    childCols <- seq_len(nChildren) + 2L
    childNames <- colnames(actgCombPatAligned)[childCols]

    ## Paternal and maternal alleles, interleaved so that each child's
    ## Pat column is followed by its Mat column
    haps <- cbind(
        as.matrix(actgCombPatAligned[, childCols, drop = FALSE]),
        as.matrix(actgCombMatAligned[, childCols, drop = FALSE])
    )
    haps <- haps[, c(rbind(childCols - 2L, childCols - 2L + nChildren)),
        drop = FALSE
    ]
    storage.mode(haps) <- "character"

    result <- data.frame(
        rsID = as.character(actgCombPatAligned$rsID),
        location = as.character(actgCombPatAligned$location),
        haps,
        stringsAsFactors = FALSE
    )
    colnames(result) <- c(
        "rsID", "location", paste0(rep(childNames, each = 2), c("Pat", "Mat"))
    )
    rownames(result) <- NULL

    result
}

#' Align SNP Data to Common Set
#'
#' Ensures a data frame contains all specified SNPs in order, filling missing
#' SNPs with NA values
#'
#' @param snpDf Data frame with SNP information (must have rsID column)
#' @param allSnps Character vector of all SNPs to include
#'
#' @return Data frame with all SNPs aligned
#' @keywords internal
.alignSnpData <- function(snpDf, allSnps) {
    if (nrow(snpDf) == 0 || length(allSnps) == 0) {
        return(snpDf)
    }

    ## Identify missing SNPs
    existingSnps <- as.character(snpDf$rsID)
    missingSnps <- setdiff(allSnps, existingSnps)

    if (length(missingSnps) == 0) {
        ## All SNPs present; reorder
        snpDf <- snpDf[match(allSnps, snpDf$rsID), , drop = FALSE]
        return(snpDf)
    }

    ## Template for missing SNPs
    nCols <- ncol(snpDf)
    nMissing <- length(missingSnps)

    missingDf <- as.data.frame(
        matrix(NA_character_,
            nrow = nMissing,
            ncol = nCols
        ),
        stringsAsFactors = FALSE
    )
    colnames(missingDf) <- colnames(snpDf)
    missingDf$rsID <- missingSnps

    ## Combine existing and missing data
    combinedDf <- rbind(snpDf, missingDf)

    ## Reorder to match allSnps
    combinedDf <- combinedDf[match(allSnps, combinedDf$rsID), , drop = FALSE]

    ## Ensure location is numeric for sorting
    if ("location" %in% colnames(combinedDf)) {
        combinedDf$location <- as.numeric(combinedDf$location)
    }

    combinedDf
}

#' Phase Parent Haplotypes Using Crossover Information
#'
#' @param actgMatrix Allele matrix for parent
#' @param xoDetectMatrix Crossover detection results
#' @param bim BIM data frame
#' @param parentId Parent identifier
#'
#' @return Data frame with two phased haplotypes
#' @keywords internal
.phaseParent <- function(actgMatrix, xoDetectMatrix, bim, parentId) {
    childCols <- setdiff(colnames(actgMatrix), c("rsID", "location", "hom"))
    childrenNoXo <- setdiff(childCols, xoDetectMatrix$childId)

    ## Haplotype 1: the alleles a child without crossovers received, or else
    ## those of the first child, switched back at its crossovers
    hap1 <- if (length(childrenNoXo) > 0) {
        as.character(actgMatrix[[childrenNoXo[1]]])
    } else {
        .undoCrossovers(actgMatrix, childCols[1], xoDetectMatrix, bim)
    }
    names(hap1) <- actgMatrix$rsID

    ## Haplotype 2 carries the other allele, except at homozygous sites
    hap2 <- vapply(seq_along(hap1), function(i) {
        .getOtherAllele(names(hap1)[i], hap1[i], bim)
    }, character(1))
    homSites <- which(actgMatrix$hom == 0)
    hap2[homSites] <- hap1[homSites]

    result <- data.frame(hap1, hap2, stringsAsFactors = FALSE)
    colnames(result) <- paste0(parentId, c("_1", "_2"))
    rownames(result) <- names(hap1)

    result
}

#' Undo a Child's Crossovers
#'
#' The alleles a child received from a parent, switched to the parent's other
#' allele wherever the child carries the parent's other haplotype, so that
#' they all lie on one parental haplotype.
#'
#' @param actgMatrix Allele table for the parent (see .phaseParent())
#' @param childName Child identifier (a column of actgMatrix)
#' @param xoDetectMatrix Crossover data frame for the parent
#' @param bim BIM data frame
#'
#' @return Character vector of alleles, one per row of actgMatrix
#' @keywords internal
.undoCrossovers <- function(actgMatrix, childName, xoDetectMatrix, bim) {
    alleles <- as.character(actgMatrix[[childName]])
    xoChild <- xoDetectMatrix[xoDetectMatrix$childId == childName, ]

    if (nrow(xoChild) == 0L) {
        return(alleles)
    }

    ## Count the crossovers completed at each SNP (end = position of the first
    ## SNP after the crossover). Where the count is odd the child carries the
    ## other parental haplotype, so switch the allele back. Homozygous sites
    ## are the same on both haplotypes.
    nXo <- findInterval(actgMatrix$location, sort(xoChild$end))
    switchIdx <- which(nXo %% 2 == 1 & actgMatrix$hom != 0)
    alleles[switchIdx] <- vapply(switchIdx, function(j) {
        .getOtherAllele(actgMatrix$rsID[j], alleles[j], bim)
    }, character(1))

    alleles
}

#' Get Complementary Allele
#'
#' @param snpId SNP identifier
#' @param allele Current allele (A, C, G, or T)
#' @param bim BIM data frame
#'
#' @return Complementary allele
#' @keywords internal
.getOtherAllele <- function(snpId, allele, bim) {
    bimIdx <- which(bim$snp == snpId)

    if (length(bimIdx) == 0L) {
        return(NA_character_)
    }

    if (bim$a1[bimIdx[1]] == allele) {
        return(bim$a2[bimIdx[1]])
    } else {
        return(bim$a1[bimIdx[1]])
    }
}

#' Convert Phased Families to a VCF Object
#'
#' @param phaseList Named list of phased DataFrames, one per family
#' @param bim BIM data frame
#' @param fam FAM data frame
#'
#' @return VCF object with phased GT for all family members
#' @keywords internal
.phaseListToVcf <- function(phaseList, bim, fam) {
    ## SNPs across all families, in genomic order
    snps <- unique(unlist(lapply(phaseList, function(x) as.character(x$rsID))))
    bimIdx <- match(snps, bim$snp)
    snpOrder <- order(bim$chr[bimIdx], bim$pos[bimIdx])
    snps <- snps[snpOrder]
    bimIdx <- bimIdx[snpOrder]
    ref <- bim$a1[bimIdx]
    alt <- bim$a2[bimIdx]

    ## Phased genotypes, with samples in FAM file order
    gtList <- lapply(phaseList, .familyGenotypes,
        snps = snps, ref = ref, alt = alt
    )
    gt <- do.call(cbind, c(
        list(matrix(character(), length(snps), 0L)),
        unname(gtList)
    ))
    sampleIds <- as.character(colnames(gt))
    sampleIds <- sampleIds[order(match(sampleIds, fam$member))]
    gt <- gt[, sampleIds, drop = FALSE]

    ## SNP positions and alleles; PLINK codes a missing allele as "0"
    snpRanges <- GRanges(
        seqnames = .ucscChromosome(bim$chr[bimIdx]),
        ranges = IRanges(bim$pos[bimIdx], width = 1L, names = snps)
    )
    altList <- DNAStringSetList(lapply(alt, function(a) {
        if (a == "0") character() else a
    }))

    VCF(
        rowRanges = snpRanges,
        colData = .vcfSampleData(sampleIds, fam),
        metadata = list(header = .phasedVcfHeader(sampleIds)),
        fixed = DataFrame(
            REF = DNAStringSet(ref),
            ALT = altList,
            QUAL = rep(NA_real_, length(snps)),
            FILTER = rep("PASS", length(snps))
        ),
        geno = SimpleList(GT = gt)
    )
}

#' Sample Information for the VCF
#'
#' @param sampleIds Character vector of sample IDs
#' @param fam FAM data frame
#'
#' @return DataFrame with familyId, fatherId, motherId, sex and affected
#'     from the FAM file, one row per sample
#' @keywords internal
.vcfSampleData <- function(sampleIds, fam) {
    famIdx <- match(sampleIds, fam$member)
    DataFrame(
        familyId = fam$pedigree[famIdx],
        fatherId = fam$father[famIdx],
        motherId = fam$mother[famIdx],
        sex = fam$sex[famIdx],
        affected = fam$affected[famIdx],
        row.names = sampleIds
    )
}

#' Phased Genotypes for One Family
#'
#' @param famPhase DataFrame of phased haplotypes for one family, with
#'     \code{<childID>Pat}/\code{<childID>Mat} and
#'     \code{<parentID>_1}/\code{<parentID>_2} columns
#' @param snps Character vector of SNP IDs (output rows)
#' @param ref Character vector of REF alleles, matching \code{snps}
#' @param alt Character vector of ALT alleles, matching \code{snps}
#'
#' @return Character matrix of phased GT values (SNPs x family members)
#' @keywords internal
.familyGenotypes <- function(famPhase, snps, ref, alt) {
    hapCols <- setdiff(colnames(famPhase), c("rsID", "location"))
    firstHap <- grep("(Pat|_1)$", hapCols, value = TRUE)
    secondHap <- sub("Pat$", "Mat", sub("_1$", "_2", firstHap))
    rows <- match(snps, famPhase$rsID)

    ## Code alleles as 0 (REF), 1 (ALT) or "." (unknown)
    alleleCodes <- function(cols) {
        alleles <- as.matrix(as.data.frame(famPhase)[rows, cols, drop = FALSE])
        codes <- matrix(".", nrow(alleles), ncol(alleles))
        codes[which(alleles == ref)] <- "0"
        codes[which(alleles == alt)] <- "1"
        codes
    }

    matrix(
        paste(alleleCodes(firstHap), alleleCodes(secondHap), sep = "|"),
        nrow = length(snps),
        dimnames = list(snps, sub("(Pat|_1)$", "", firstHap))
    )
}

#' VCF Header for Phased Genotypes
#'
#' @param sampleIds Character vector of sample IDs
#'
#' @return VCFHeader object
#' @keywords internal
.phasedVcfHeader <- function(sampleIds) {
    metaRow <- function(key, value) DataFrame(Value = value, row.names = key)

    VCFHeader(
        samples = sampleIds,
        header = DataFrameList(
            fileformat = metaRow("fileformat", "VCFv4.2"),
            fileDate = metaRow("fileDate", format(Sys.Date(), "%Y%m%d")),
            source = metaRow("source", "inferRecom::xoPhase"),
            FORMAT = DataFrame(
                Number = "1", Type = "String",
                Description = "Phased genotype", row.names = "GT"
            )
        )
    )
}
