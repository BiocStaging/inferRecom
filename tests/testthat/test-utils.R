## ============================================================================
## Shared Internal Helper Tests
## ============================================================================

test_that(".firstNonNA returns the first non-NA value", {
    expect_equal(inferRecom:::.firstNonNA(c(NA, NA, 5, 10)), 5)
    expect_equal(inferRecom:::.firstNonNA(c(1, 2, 3)), 1)
    expect_equal(inferRecom:::.firstNonNA(c(NA, "F1P1", "F1P2")), "F1P1")

    # NA of the same type when there is no non-NA value
    expect_identical(inferRecom:::.firstNonNA(c(NA, NA, NA)), NA)
    expect_identical(
        inferRecom:::.firstNonNA(c(NA_character_, NA)), NA_character_
    )
    expect_identical(inferRecom:::.firstNonNA(character()), NA_character_)
})

test_that(".orderedFamily orders genotypes father, mother, children", {
    # FAM rows listed in a different order from the genotype columns
    fam <- data.frame(
        pedigree = "F1",
        member = c("C2", "M", "C1", "F"),
        father = c("F", NA, "F", NA),
        mother = c("M", NA, "M", NA),
        row.names = c("C2", "M", "C1", "F")
    )
    geno <- matrix(1:12, nrow = 3, dimnames = list(
        paste0("rs", 1:3), c("C1", "F", "C2", "M")
    ))

    family <- inferRecom:::.orderedFamily("F1", fam, geno, minMembers = 4L)
    expect_equal(colnames(family$geno), c("F", "M", "C1", "C2"))
    expect_equal(family$geno[, "C2"], geno[, "C2"])
    expect_equal(family$fam$member, c("F", "M", "C1", "C2"))

    # Too few genotyped members, or a parent not genotyped
    expect_null(inferRecom:::.orderedFamily("F1", fam, geno, minMembers = 5L))
    expect_null(inferRecom:::.orderedFamily(
        "F1", fam, geno[, c("C1", "C2", "M")],
        minMembers = 3L
    ))
})

test_that("argument checks accept valid values", {
    plinkFile <- "simCEU"
    out <- NULL
    rsOnly <- TRUE
    famIds <- c("F1", "F2")
    minMb <- 0.5
    cmFilter <- 0

    expect_silent(inferRecom:::.checkString(plinkFile))
    expect_silent(inferRecom:::.checkString(out, allowNull = TRUE))
    expect_silent(inferRecom:::.checkFlag(rsOnly))
    expect_silent(inferRecom:::.checkCharacterOrNull(famIds))
    expect_silent(inferRecom:::.checkNumber(minMb))
    expect_silent(inferRecom:::.checkNumber(cmFilter, allowZero = TRUE))
})

test_that("argument checks name the argument in their errors", {
    plinkFile <- c("a", "b")
    rsOnly <- NA
    famIds <- 1
    minSnps <- 0
    snpFilter <- -1
    minCm <- "1"

    expect_error(
        inferRecom:::.checkString(plinkFile),
        "'plinkFile' must be a single character string"
    )
    expect_error(
        inferRecom:::.checkFlag(rsOnly),
        "'rsOnly' must be a single logical value"
    )
    expect_error(
        inferRecom:::.checkCharacterOrNull(famIds),
        "'famIds' must be a character vector or NULL"
    )
    expect_error(
        inferRecom:::.checkNumber(minSnps, type = "integer"),
        "'minSnps' must be a positive integer"
    )
    expect_error(
        inferRecom:::.checkNumber(snpFilter,
            allowZero = TRUE, type = "integer"
        ),
        "'snpFilter' must be a non-negative integer"
    )
    expect_error(
        inferRecom:::.checkNumber(minCm, allowNull = TRUE),
        "'minCm' must be a positive numeric value or NULL"
    )
})

test_that(".readGeneticMap returns the chr, pos and cM columns", {
    skip_if_not(checkMapFile("female"), "Test data not available")

    geneticMap <- inferRecom:::.readGeneticMap(loadSimCeuData()$mapFemale)

    expect_s3_class(geneticMap, "data.frame")
    expect_equal(colnames(geneticMap), c("chr", "pos", "cM"))
    expect_true(nrow(geneticMap) > 0)
})

test_that(".readGeneticMap rejects missing files and columns", {
    expect_error(
        inferRecom:::.readGeneticMap("nonexistent_map.txt"),
        "Genetic map file not found"
    )

    badMap <- tempfile(fileext = ".txt")
    on.exit(unlink(badMap))
    write.table(
        data.frame(chr = "chr4", position = 1, cM = 0), badMap,
        sep = "\t", quote = FALSE, row.names = FALSE
    )
    expect_error(
        inferRecom:::.readGeneticMap(badMap),
        "mapFile must have columns: chr, pos, cM"
    )
})

test_that(".asGRangesResult builds GRanges and splits by case/control", {
    resultDf <- DataFrame(
        childId = c("C1", "C2", "C3"),
        startPos = c(100, 200, 300),
        finishPos = c(150, 250, 350)
    )
    famData <- data.frame(member = c("C1", "C2", "C3"), affected = c(1, NA, 1))

    gr <- inferRecom:::.asGRangesResult(resultDf, "chr4", FALSE, famData,
        idCol = "childId"
    )
    expect_s4_class(gr, "GRanges")
    expect_equal(as.character(seqnames(gr)), rep("chr4", 3))
    expect_equal(start(gr), c(100, 200, 300))
    expect_equal(gr$childId, c("C1", "C2", "C3"))

    grl <- inferRecom:::.asGRangesResult(resultDf, "chr4", TRUE, famData,
        idCol = "childId"
    )
    expect_s4_class(grl, "GRangesList")
    expect_equal(grl$case$childId, c("C1", "C3"))
    expect_equal(grl$control$childId, "C2")
})

test_that(".asGRangesResult uses UCSC-style chromosome names", {
    resultDf <- DataFrame(sampleId = "S1", startPos = 100, finishPos = 150)

    seqNames <- vapply(list(4L, "4", "chr4", "CHR4"), function(chromosome) {
        gr <- inferRecom:::.asGRangesResult(resultDf, chromosome, FALSE, NULL,
            idCol = "sampleId"
        )
        as.character(seqnames(gr))
    }, character(1))

    expect_equal(seqNames, rep("chr4", 4))
})

test_that(".ucscChromosome converts PLINK chromosome codes", {
    expect_equal(
        inferRecom:::.ucscChromosome(
            c("1", "chr22", "23", "X", "24", "26", "MT")
        ),
        c("chr1", "chr22", "chrX", "chrX", "chrY", "chrM", "chrM")
    )
    expect_equal(inferRecom:::.ucscChromosome(character()), character())
})

test_that(".asGRangesResult returns empty results of the same type", {
    emptyDf <- DataFrame(startPos = numeric(), finishPos = numeric())

    expect_equal(
        inferRecom:::.asGRangesResult(emptyDf, "chr4", FALSE, NULL, "childId"),
        GRanges()
    )

    emptyCc <- inferRecom:::.asGRangesResult(
        emptyDf, "chr4", TRUE, NULL, "childId"
    )
    expect_s4_class(emptyCc, "GRangesList")
    expect_equal(names(emptyCc), c("case", "control"))
})
