## ============================================================================
## Setup and Helper Functions
## ============================================================================

test_that("xoPhase helper functions work correctly", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()

    # Test .checkSameChromosome with empty GRanges
    emptyGr <- GRanges()
    expect_true(inferRecom:::.checkSameChromosome(emptyGr))

    # Test .checkSameChromosome with single chromosome
    singleChrGr <- GRanges(
        seqnames = rep("chr4", 3),
        ranges = IRanges(start = c(1000, 2000, 3000), width = 100)
    )
    expect_true(inferRecom:::.checkSameChromosome(singleChrGr))

    # Test .checkSameChromosome with multiple chromosomes
    multiChrGr <- GRanges(
        seqnames = c("chr4", "chr5", "chr4"),
        ranges = IRanges(start = c(1000, 2000, 3000), width = 100)
    )
    expect_false(inferRecom:::.checkSameChromosome(multiChrGr))
})

test_that(".phaseParent uses crossover data to phase parental haplotypes", {
    # Parental haplotypes; the parent is homozygous at s5 (hom = 0)
    hap1 <- c("A", "G", "A", "G", "A", "G")
    hap2 <- c("G", "A", "G", "A", "A", "A")
    bim <- data.frame(snp = paste0("s", 1:6), a1 = "A", a2 = "G")
    makeActg <- function(...) {
        data.frame(
            rsID = paste0("s", 1:6), location = seq(100, 600, by = 100),
            ..., hom = c(1, 1, 1, 1, 0, 1), row.names = paste0("s", 1:6)
        )
    }
    expectHaps <- function(res) {
        expect_equal(unname(res$P1_1), hap1)
        expect_equal(unname(res$P1_2), hap2)
    }

    # C1 switches to hap2 after s3; C2 inherits hap1 intact. C2 must be used
    # even though C1 is the first child.
    actg <- makeActg(C1 = c(hap1[1:3], hap2[4:6]), C2 = hap1)
    xo <- data.frame(childId = "C1", start = 300, end = 400)
    expectHaps(inferRecom:::.phaseParent(actg, xo, bim, "P1"))

    # Every child has crossovers, so phase from C1's double crossover
    actg <- makeActg(
        C1 = c(hap1[1:3], hap2[4:5], hap1[6]),
        C2 = c(hap2[1:2], hap1[3:6])
    )
    xo <- data.frame(
        childId = c("C1", "C1", "C2"),
        start = c(300, 500, 200), end = c(400, 600, 300)
    )
    expectHaps(inferRecom:::.phaseParent(actg, xo, bim, "P1"))

    # Crossover boundary SNPs need not be rows of the phasing table
    actg <- makeActg(C1 = c(hap1[1:3], hap2[4:6]), C2 = c(hap2[1:2], hap1[3:6]))
    xo <- data.frame(
        childId = c("C1", "C2"), start = c(320, 220), end = c(380, 280)
    )
    expectHaps(inferRecom:::.phaseParent(actg, xo, bim, "P1"))
})

## ============================================================================
## Input Validation Tests
## ============================================================================

test_that("xoPhase validates input arguments", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()

    # Load crossover data
    xoPat <- readRDS(file.path(dirname(testData$plinkFile), "xoPat.rds"))
    xoMat <- readRDS(file.path(dirname(testData$plinkFile), "xoMat.rds"))

    # Test invalid plinkFile
    expect_error(
        xoPhase(
            plinkFile = c("file1", "file2"),
            xoDetectPaternal = xoPat,
            xoDetectMaternal = xoMat
        ),
        "'plinkFile' must be a single character string"
    )

    expect_error(
        xoPhase(
            plinkFile = 123,
            xoDetectPaternal = xoPat,
            xoDetectMaternal = xoMat
        ),
        "'plinkFile' must be a single character string"
    )

    # Test invalid xoDetect objects
    expect_error(
        xoPhase(
            plinkFile = testData$plinkFile,
            xoDetectPaternal = data.frame(),
            xoDetectMaternal = xoMat
        ),
        "'xoDetectPaternal' must be a GRanges object"
    )

    expect_error(
        xoPhase(
            plinkFile = testData$plinkFile,
            xoDetectPaternal = xoPat,
            xoDetectMaternal = list()
        ),
        "'xoDetectMaternal' must be a GRanges object"
    )

    # Test invalid famIds
    expect_error(
        xoPhase(
            plinkFile = testData$plinkFile,
            xoDetectPaternal = xoPat,
            xoDetectMaternal = xoMat,
            famIds = 123
        ),
        "'famIds' must be a character vector or NULL"
    )

    # Test invalid rsOnly
    expect_error(
        xoPhase(
            plinkFile = testData$plinkFile,
            xoDetectPaternal = xoPat,
            xoDetectMaternal = xoMat,
            rsOnly = "yes"
        ),
        "'rsOnly' must be a single logical value"
    )
})

test_that("xoPhase rejects mixed chromosomes", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()

    # Create GRanges with multiple chromosomes
    multiChrGr <- GRanges(
        seqnames = c("chr4", "chr5"),
        ranges = IRanges(start = c(1000, 2000), width = 100),
        childId = c("child1", "child1"),
        familyId = c("fam1", "fam1"),
        startSnp = c("rs1", "rs2"),
        finishSnp = c("rs3", "rs4")
    )

    singleChrGr <- GRanges(
        seqnames = "chr4",
        ranges = IRanges(start = 1000, width = 100),
        childId = "child1",
        familyId = "fam1",
        startSnp = "rs1",
        finishSnp = "rs2"
    )

    expect_error(
        xoPhase(
            plinkFile = testData$plinkFile,
            xoDetectPaternal = multiChrGr,
            xoDetectMaternal = singleChrGr
        ),
        "GRanges object can only contain ranges from a single chromosome"
    )
})

## ============================================================================
## Basic Functionality Tests
## ============================================================================

loadXo <- function(testData) {
    dataPath <- dirname(testData$plinkFile)
    list(
        pat = readRDS(file.path(dataPath, "xoPat.rds")),
        mat = readRDS(file.path(dataPath, "xoMat.rds"))
    )
}

test_that("xoPhase returns a VCF object", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()
    xo <- loadXo(testData)

    testResult <- xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xo$pat,
        xoDetectMaternal = xo$mat
    )

    expect_s4_class(testResult, "VCF")
    expect_true(nrow(testResult) > 0)
    expect_true(ncol(testResult) > 0)
    expect_true("GT" %in% names(geno(testResult)))
    expect_equal(
        colnames(colData(testResult)),
        c("familyId", "fatherId", "motherId", "sex", "affected")
    )

    # Every genotype is phased and coded 0 (REF), 1 (ALT) or . (unknown)
    expect_true(all(grepl("^[01.]\\|[01.]$", geno(testResult)$GT)))
})

test_that("xoPhase phases multiple families", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()
    xo <- loadXo(testData)

    testResult <- xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xo$pat,
        xoDetectMaternal = xo$mat
    )

    expect_true(length(unique(colData(testResult)$familyId)) > 1)
    expect_equal(anyDuplicated(colnames(testResult)), 0L)
})

test_that("xoPhase handles specific family IDs", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()
    xo <- loadXo(testData)

    testResult <- xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xo$pat,
        xoDetectMaternal = xo$mat,
        famIds = "F6"
    )

    expect_equal(unique(colData(testResult)$familyId), "F6")
})

test_that("xoPhase omits families without phase information", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()
    xo <- loadXo(testData)

    expect_message(
        testResult <- xoPhase(
            plinkFile = testData$plinkFile,
            xoDetectPaternal = xo$pat,
            xoDetectMaternal = xo$mat,
            famIds = c("F6", "F99")
        ),
        "No phase information for families: F99"
    )
    expect_equal(unique(colData(testResult)$familyId), "F6")

    # Still a VCF when no family can be phased
    testResult <- suppressMessages(xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xo$pat,
        xoDetectMaternal = xo$mat,
        famIds = "F99"
    ))
    expect_s4_class(testResult, "VCF")
    expect_equal(ncol(testResult), 0L)
})

test_that("xoPhase filters rs SNPs correctly", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()
    xo <- loadXo(testData)

    testResultRs <- xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xo$pat,
        xoDetectMaternal = xo$mat,
        rsOnly = TRUE
    )

    testResultAll <- xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xo$pat,
        xoDetectMaternal = xo$mat,
        rsOnly = FALSE
    )

    expect_true(all(grepl("^rs", rownames(testResultRs))))
    expect_true(nrow(testResultAll) >= nrow(testResultRs))
})

## ============================================================================
## Output Format Tests
## ============================================================================

test_that("xoPhase output is a RangedSummarizedExperiment", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()
    xo <- loadXo(testData)

    testResult <- xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xo$pat,
        xoDetectMaternal = xo$mat
    )

    expect_s4_class(testResult, "RangedSummarizedExperiment")
    expect_true(all(c("REF", "ALT") %in% colnames(mcols(rowRanges(testResult)))))
    expect_true(all(width(rowRanges(testResult)) == 1L))
})

test_that("xoPhase output can be written with writeVcf", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()
    xo <- loadXo(testData)

    testResult <- xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xo$pat,
        xoDetectMaternal = xo$mat,
        famIds = c("F4", "F6")
    )

    vcfFile <- tempfile(fileext = ".vcf")
    on.exit(unlink(vcfFile))
    writeVcf(testResult, vcfFile)

    expect_true(any(grepl("^##fileformat=VCF", readLines(vcfFile, n = 1))))

    roundTrip <- readVcf(vcfFile)
    expect_equal(colnames(roundTrip), colnames(testResult))
    expect_equal(unname(geno(roundTrip)$GT), unname(geno(testResult)$GT))
})

## ============================================================================
## Phasing Quality Tests
## ============================================================================

test_that("xoPhase produces SNPs in genomic order", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()
    xo <- loadXo(testData)

    testResult <- xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xo$pat,
        xoDetectMaternal = xo$mat
    )

    positions <- start(rowRanges(testResult))
    expect_false(anyNA(positions))
    expect_true(all(diff(positions) >= 0))
    expect_true(all(nchar(rownames(testResult)) > 0))
})

test_that("xoPhase maintains parent-child relationships", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()
    xo <- loadXo(testData)

    testResult <- xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xo$pat,
        xoDetectMaternal = xo$mat
    )

    # Each child's parents are samples from the same family
    samples <- colData(testResult)
    children <- samples[!is.na(samples$fatherId), ]
    expect_true(nrow(children) > 0)
    expect_true(all(children$fatherId %in% rownames(samples)))
    expect_true(all(children$motherId %in% rownames(samples)))
    expect_equal(
        samples[children$fatherId, "familyId"], children$familyId
    )
})

test_that("xoPhase handles homozygous sites correctly", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()
    xo <- loadXo(testData)

    testResult <- xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xo$pat,
        xoDetectMaternal = xo$mat
    )

    # Parents should be homozygous at some sites
    parents <- is.na(colData(testResult)$fatherId)
    parentGt <- geno(testResult)$GT[, parents]
    expect_true(any(parentGt %in% c("0|0", "1|1")))
})

## ============================================================================
## Edge Case Tests
## ============================================================================

test_that("xoPhase handles empty crossover data", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()

    # Create empty GRanges
    emptyGr <- GRanges()

    expect_warning(
        expect_warning(
            testResult <- xoPhase(
                plinkFile = testData$plinkFile,
                xoDetectPaternal = emptyGr,
                xoDetectMaternal = emptyGr
            ),
            "Empty GRanges object provided"
        ),
        "Empty GRanges object provided"
    )

    # Should still return a VCF
    expect_s4_class(testResult, "VCF")
})

test_that("xoPhase handles missing PLINK files", {
    xoPat <- GRanges(
        seqnames = "chr4",
        ranges = IRanges(start = 1000, width = 100),
        childId = "child1",
        familyId = "fam1"
    )

    expect_error(
        xoPhase(
            plinkFile = "nonexistent_file",
            xoDetectPaternal = xoPat,
            xoDetectMaternal = xoPat
        ),
        "PLINK file\\(s\\) not found"
    )
})

## ============================================================================
## Parallel Processing Tests
## ============================================================================

test_that("xoPhase works with parallel processing", {
    skip_if_not(checkTestData(), "Test data not available")
    skip_on_os("windows") # MulticoreParam doesn't work on Windows

    testData <- loadSimCeuData()
    xoPat <- readRDS(file.path(dirname(testData$plinkFile), "xoPat.rds"))
    xoMat <- readRDS(file.path(dirname(testData$plinkFile), "xoMat.rds"))

    # Serial execution
    resultSerial <- xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xoPat,
        xoDetectMaternal = xoMat,
        BPPARAM = SerialParam()
    )

    # Parallel execution
    resultParallel <- xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xoPat,
        xoDetectMaternal = xoMat,
        BPPARAM = MulticoreParam(workers = 2)
    )

    # Results should be equivalent
    expect_equal(dimnames(resultSerial), dimnames(resultParallel))
    expect_equal(geno(resultSerial)$GT, geno(resultParallel)$GT)
})

## ============================================================================
## Internal Helper Function Tests
## ============================================================================

test_that(".matchAlleleCode works correctly", {
    # Father heterozygous (1), mother homozygous (2)
    expect_equal(inferRecom:::.matchAlleleCode("father", 1, 2, 2), 2L)
    expect_equal(inferRecom:::.matchAlleleCode("father", 1, 2, 1), 1L)

    # Father heterozygous (1), mother homozygous (0)
    expect_equal(inferRecom:::.matchAlleleCode("father", 1, 0, 0), 1L)
    expect_equal(inferRecom:::.matchAlleleCode("father", 1, 0, 1), 2L)

    # Mother heterozygous (1), father homozygous (2)
    expect_equal(inferRecom:::.matchAlleleCode("mother", 2, 1, 2), 2L)
    expect_equal(inferRecom:::.matchAlleleCode("mother", 2, 1, 1), 1L)

    # Mother heterozygous (1), father homozygous (0)
    expect_equal(inferRecom:::.matchAlleleCode("mother", 0, 1, 0), 1L)
    expect_equal(inferRecom:::.matchAlleleCode("mother", 0, 1, 1), 2L)

    # Unexpected combination
    expect_true(is.na(inferRecom:::.matchAlleleCode("father", 0, 0, 0)))
})

test_that(".getInformativeGenotypes keeps 0/0 parents, drops missing data", {
    # Columns: father, mother, three children
    geno <- rbind(
        rs1 = c(0, 1, 0, 1, 1), # father 0/0, mother heterozygous
        rs2 = c(1, 2, 1, 2, 2), # father heterozygous, mother 2/2
        rs3 = c(1, 0, NA, 0, 1), # a child's genotype is missing
        rs4 = c(1, 1, 0, 1, 2), # both parents heterozygous
        rs5 = c(0, 1, 2, 1, 0) # child 2/2 is not possible
    )

    informative <- inferRecom:::.getInformativeGenotypes(geno)
    expect_equal(rownames(informative), c("rs1", "rs2"))
})

test_that(".addNonInformativeSNPs phases parents homozygous for different alleles", {
    # One informative SNP (rs1) has been phased already
    phasedData <- data.frame(
        rsID = "rs1", location = 100,
        C1Pat = "A", C1Mat = "G", C2Pat = "G", C2Mat = "G",
        F_1 = "A", F_2 = "G", M_1 = "G", M_2 = "G"
    )
    # rs2: father 0/0 (C/C), mother 2/2 (T/T); rs3: father 2/2 (A/A),
    # mother 0/0 (T/T), with child C2's genotype inconsistent
    geno <- rbind(rs1 = c(1, 2, 1, 2), rs2 = c(0, 2, 1, 1), rs3 = c(2, 0, 1, 0))
    colnames(geno) <- c("F", "M", "C1", "C2")
    bim <- data.frame(
        chr = 4, snp = c("rs1", "rs2", "rs3"), cm = 0, pos = c(100, 200, 300),
        a1 = c("A", "C", "T"), a2 = c("G", "T", "A")
    )

    phased <- inferRecom:::.addNonInformativeSNPs(
        phasedData, geno, bim, "F", "M", 1L, 2L
    )
    alleles <- function(snp) unlist(phased[phased$rsID == snp, -(1:2)])
    expect_equal(
        unname(alleles("rs2")), c("C", "T", "C", "T", "C", "C", "T", "T")
    )
    expect_equal(
        unname(alleles("rs3")), c("A", "T", NA, NA, "A", "A", "T", "T")
    )
})

test_that(".convertToNucleotides looks up alleles from the BIM file", {
    bim <- data.frame(
        snp = c("rs1", "rs2"), pos = c(1000L, 2000L),
        a1 = c("A", "C"), a2 = c("G", "T")
    )
    alleleCodes <- data.frame(
        rsID = c("rs1", "rs9", "rs2"),
        C1 = c(1L, 2L, 0L),
        C2 = c(2L, 1L, 1L),
        row.names = c("rs1", "rs9", "rs2")
    )

    nucleotides <- inferRecom:::.convertToNucleotides(alleleCodes, bim)
    expect_equal(colnames(nucleotides), c("rsID", "location", "C1", "C2"))
    expect_equal(nucleotides$location, c("1000", "", "2000"))
    expect_equal(nucleotides$C1, c("A", "", "")) # rs9 unknown; code 0 invalid
    expect_equal(nucleotides$C2, c("G", "", "C"))
    expect_null(inferRecom:::.convertToNucleotides(NULL, bim))
})

test_that(".getOtherAllele retrieves complementary alleles", {
    # Create mock BIM data
    bim <- data.frame(
        chr = c(4, 4, 4),
        snp = c("rs1", "rs2", "rs3"),
        cm = c(0, 0, 0),
        pos = c(1000, 2000, 3000),
        a1 = c("A", "C", "G"),
        a2 = c("T", "G", "A"),
        stringsAsFactors = FALSE
    )

    # Test getting other allele
    expect_equal(inferRecom:::.getOtherAllele("rs1", "A", bim), "T")
    expect_equal(inferRecom:::.getOtherAllele("rs1", "T", bim), "A")
    expect_equal(inferRecom:::.getOtherAllele("rs2", "C", bim), "G")
    expect_equal(inferRecom:::.getOtherAllele("rs3", "G", bim), "A")

    # Test missing SNP
    expect_true(is.na(inferRecom:::.getOtherAllele("rs999", "A", bim)))
})

test_that(".phaseChildren pairs each child's paternal and maternal alleles", {
    pat <- data.frame(
        rsID = c("rs1", "rs2"), location = c("100", "200"),
        C1 = c("A", "C"), C2 = c("G", "T"), hom = 1
    )
    mat <- data.frame(
        rsID = c("rs1", "rs2"), location = c("100", "200"),
        C1 = c("G", "T"), C2 = c("A", NA), hom = 1
    )

    phased <- inferRecom:::.phaseChildren(pat, mat)
    expect_equal(
        colnames(phased),
        c("rsID", "location", "C1Pat", "C1Mat", "C2Pat", "C2Mat")
    )
    expect_equal(phased$C1Pat, c("A", "C"))
    expect_equal(phased$C1Mat, c("G", "T"))
    expect_equal(phased$C2Mat, c("A", NA))
})

test_that(".alignSnpData aligns SNPs correctly", {
    # Create test data
    snpDf <- data.frame(
        rsID = c("rs1", "rs3", "rs5"),
        location = c(1000, 3000, 5000),
        value = c("A", "C", "G"),
        stringsAsFactors = FALSE
    )

    allSnps <- c("rs1", "rs2", "rs3", "rs4", "rs5")

    aligned <- inferRecom:::.alignSnpData(snpDf, allSnps)

    expect_equal(nrow(aligned), 5)
    expect_equal(as.character(aligned$rsID), allSnps)
    expect_equal(aligned$value[1], "A")
    expect_true(is.na(aligned$value[2])) # rs2 was missing
    expect_equal(aligned$value[3], "C")
    expect_true(is.na(aligned$value[4])) # rs4 was missing
    expect_equal(aligned$value[5], "G")
})

test_that(".inheritedHaplotype reads the haplotype from both sides", {
    hap1 <- c("A", "C", "G", "T")
    hap2 <- c("G", "T", "A", "C")
    positions <- c(100, 200, 400, 500)
    haplotype <- function(child, snpPos = 300) {
        inferRecom:::.inheritedHaplotype(child, hap1, hap2, positions, snpPos)
    }

    expect_equal(haplotype(c("A", "C", "G", "T")), 1L) # hap1 on both sides
    expect_equal(haplotype(c("G", "T", "A", "C")), 2L) # hap2 on both sides
    expect_true(is.na(haplotype(c("A", "C", "A", "C")))) # crossover between
    expect_equal(haplotype(c("A", "C", "A", "C"), snpPos = 50), 1L) # one side

    # Only SNPs where the parent's haplotypes differ and the child's allele
    # is known are used
    expect_true(is.na(inferRecom:::.inheritedHaplotype(
        c("A", NA), c("A", "C"), c("A", "T"), c(100, 400), 300
    )))
})

test_that(".phaseOneHeterozygous phases the parent from the children", {
    # Phased SNPs either side of position 300. The father is heterozygous
    # there; C1 inherited his haplotype 1, C2 his haplotype 2, and C3 has a
    # crossover between the sides.
    phasedData <- data.frame(
        rsID = paste0("rs", 1:4), location = c(100, 200, 400, 500),
        C1Pat = c("A", "C", "G", "T"), C1Mat = "A",
        C2Pat = c("G", "T", "A", "C"), C2Mat = "A",
        C3Pat = c("A", "C", "A", "C"), C3Mat = "A",
        F_1 = c("A", "C", "G", "T"), F_2 = c("G", "T", "A", "C"),
        M_1 = "A", M_2 = "A"
    )
    ids <- list(father = "F", mother = "M", children = c("C1", "C2", "C3"))
    phase <- function(childGeno) {
        inferRecom:::.phaseOneHeterozygous(
            c(F = 1, M = 0), childGeno, c("C", "T"), ids, 300, phasedData
        )
    }

    # Father C/T, mother C/C: C1 and C3 received T, C2 received C
    phased <- phase(c(C1 = 1, C2 = 0, C3 = 1))
    expect_equal(
        phased[c("C1Pat", "C1Mat", "C2Pat", "C2Mat", "C3Pat", "C3Mat")],
        c(
            C1Pat = "T", C1Mat = "C", C2Pat = "C", C2Mat = "C",
            C3Pat = "T", C3Mat = "C"
        )
    )
    # C1 (haplotype 1, T) and C2 (haplotype 2, C) both put T on haplotype 1
    expect_equal(phased[c("F_1", "F_2")], c(F_1 = "T", F_2 = "C"))
    expect_equal(phased[c("M_1", "M_2")], c(M_1 = "C", M_2 = "C"))

    # Children that disagree leave the father unphased
    conflict <- phase(c(C1 = 1, C2 = 1, C3 = 1))
    expect_false(any(c("F_1", "F_2") %in% names(conflict)))
    expect_equal(conflict[["C2Pat"]], "T")

    # A child that cannot have these parents is left unphased
    inconsistent <- phase(c(C1 = 1, C2 = 2, C3 = 1))
    expect_false(any(c("C2Pat", "C2Mat") %in% names(inconsistent)))
    expect_equal(inconsistent[c("F_1", "F_2")], c(F_1 = "T", F_2 = "C"))
})

test_that(".phaseBothHomozygous leaves inconsistent children unphased", {
    ids <- list(father = "F", mother = "M", children = c("C1", "C2"))
    phased <- inferRecom:::.phaseBothHomozygous(
        c(F = 0, M = 0), c(C1 = 0, C2 = 1), c("A", "G"), ids
    )

    expect_equal(phased[c("C1Pat", "C1Mat")], c(C1Pat = "A", C1Mat = "A"))
    expect_false(any(c("C2Pat", "C2Mat") %in% names(phased)))
    expect_equal(unname(phased[c("F_1", "F_2", "M_1", "M_2")]), rep("A", 4))
})

## ============================================================================
## Integration Tests
## ============================================================================

test_that("xoPhase integrates with xoDetect output", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()

    # Run xoDetect to generate crossovers
    xoPat <- xoDetect(
        plinkFile = testData$plinkFile,
        mapFile = testData$mapMale,
        familySize = 3,
        parent = "father"
    )

    xoMat <- xoDetect(
        plinkFile = testData$plinkFile,
        mapFile = testData$mapFemale,
        familySize = 3,
        parent = "mother"
    )

    # Phase using detected crossovers
    phased <- xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xoPat,
        xoDetectMaternal = xoMat
    )

    # Verify output
    expect_s4_class(phased, "VCF")
    expect_true(ncol(phased) > 0)
})

test_that("xoPhase output structure is consistent", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()
    xoPat <- readRDS(file.path(dirname(testData$plinkFile), "xoPat.rds"))
    xoMat <- readRDS(file.path(dirname(testData$plinkFile), "xoMat.rds"))

    # Run multiple times
    result1 <- xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xoPat,
        xoDetectMaternal = xoMat
    )

    result2 <- xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xoPat,
        xoDetectMaternal = xoMat
    )

    # Results should be identical
    expect_equal(dimnames(result1), dimnames(result2))
    expect_equal(geno(result1)$GT, geno(result2)$GT)
    expect_equal(colData(result1), colData(result2))
})
