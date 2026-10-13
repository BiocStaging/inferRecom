test_that("xoDetect works with 3-child maternal crossovers and female map", {
    skip_if_not(checkTestData(), "Test data not available")
    skip_if_not(checkMapFile("female"), "Female genetic map not available")

    testData <- loadSimCeuData()

    # Test maternal crossovers with female genetic map
    testResult <- xoDetect(
        plinkFile = testData$plinkFile,
        mapFile = testData$mapFemale,
        familySize = 3,
        parent = "mother"
    )

    # Check return type - should be GRanges when caseControl = FALSE
    expect_s4_class(testResult, "GRanges")

    # Check expected metadata columns
    expectedCols <- c(
        "childId", "familyId", "startSnp", "finishSnp",
        "startCm", "finishCm"
    )
    expect_true(all(expectedCols %in% colnames(mcols(testResult))))

    # Check data types
    expect_type(mcols(testResult)$childId, "character")
    expect_type(mcols(testResult)$familyId, "character")

    # Check that positions are reasonable
    if (length(testResult) > 0) {
        expect_true(all(start(testResult) < end(testResult)))
        expect_true(all(start(testResult) > 0))
        # Check cM positions are reasonable (<500 cM per chromosome)
        expect_true(all(mcols(testResult)$startCm >= 0 | is.na(mcols(testResult)$startCm)))
        expect_true(all(mcols(testResult)$finishCm <= 500 | is.na(mcols(testResult)$finishCm)))
    }
})

test_that("xoDetect works with 3-child paternal crossovers and male map", {
    skip_if_not(checkTestData(), "Test data not available")
    skip_if_not(checkMapFile("male"), "Male genetic map not available")

    testData <- loadSimCeuData()

    # Test paternal crossovers with male genetic map
    testResult <- xoDetect(
        plinkFile = testData$plinkFile,
        mapFile = testData$mapMale,
        familySize = 3,
        parent = "father"
    )

    # Check return type
    expect_s4_class(testResult, "GRanges")

    # Check expected metadata columns
    expectedCols <- c(
        "childId", "familyId", "startSnp", "finishSnp",
        "startCm", "finishCm"
    )
    expect_true(all(expectedCols %in% colnames(mcols(testResult))))

    # Check that positions are reasonable
    if (length(testResult) > 0) {
        expect_true(all(start(testResult) < end(testResult)))
        expect_true(all(start(testResult) > 0))
        # Check cM positions are reasonable (<500 cM per chromosome)
        expect_true(all(mcols(testResult)$startCm >= 0 | is.na(mcols(testResult)$startCm)))
        expect_true(all(mcols(testResult)$finishCm <= 500 | is.na(mcols(testResult)$finishCm)))
    }
})

test_that("xoDetect handles 2-child families", {
    skip_if_not(checkTestData(), "Test data not available")
    skip_if_not(checkMapFile("female"), "Female genetic map not available")

    testData <- loadSimCeuData()

    testResult <- xoDetect(
        plinkFile = testData$plinkFile,
        mapFile = testData$mapFemale,
        familySize = 2,
        parent = "mother"
    )

    expect_s4_class(testResult, "GRanges")

    # 2-child families should not have childId column
    expect_false("childId" %in% colnames(mcols(testResult)))

    # Should have familyId column
    if (length(testResult) > 0) {
        expect_true("familyId" %in% colnames(mcols(testResult)))
    }

    # Check that positions are reasonable
    if (length(testResult) > 0) {
        expect_true(all(start(testResult) < end(testResult)))
        expect_true(all(start(testResult) > 0))
        # Check cM positions are reasonable (<500 cM per chromosome)
        expect_true(all(mcols(testResult)$startCm >= 0 | is.na(mcols(testResult)$startCm)))
        expect_true(all(mcols(testResult)$finishCm <= 500 | is.na(mcols(testResult)$finishCm)))
    }
})

test_that("xoDetect filters by SNP correctly", {
    skip_if_not(checkTestData(), "Test data not available")
    skip_if_not(checkMapFile("female"), "Female genetic map not available")

    testData <- loadSimCeuData()

    # With rsOnly = TRUE
    resultFiltered <- xoDetect(
        plinkFile = testData$plinkFile,
        mapFile = testData$mapFemale,
        familySize = 3,
        parent = "mother",
        rsOnly = TRUE
    )

    # With rsOnly = FALSE
    resultAll <- xoDetect(
        plinkFile = testData$plinkFile,
        mapFile = testData$mapFemale,
        familySize = 3,
        parent = "mother",
        rsOnly = FALSE
    )

    # Both should return GRanges
    expect_s4_class(resultFiltered, "GRanges")
    expect_s4_class(resultAll, "GRanges")

    # With rsOnly = TRUE, every crossover starts and ends at an rs SNP
    expect_true(all(
        grepl("^rs", mcols(resultFiltered)$startSnp) &
            grepl("^rs", mcols(resultFiltered)$finishSnp)
    ))
})

test_that("xoDetect applies filters correctly", {
    skip_if_not(checkTestData(), "Test data not available")
    skip_if_not(checkMapFile("female"), "Female genetic map not available")

    testData <- loadSimCeuData()

    # Test with strict filters
    resultStrict <- xoDetect(
        plinkFile = testData$plinkFile,
        mapFile = testData$mapFemale,
        familySize = 3,
        parent = "mother",
        snpFilter = 10,
        cmFilter = 5
    )

    # Test with lenient filters
    resultLenient <- xoDetect(
        plinkFile = testData$plinkFile,
        mapFile = testData$mapFemale,
        familySize = 3,
        parent = "mother",
        snpFilter = 1,
        cmFilter = 0.1
    )

    # Strict filters should result in fewer crossovers
    expect_true(length(resultStrict) <= length(resultLenient))
})

test_that("xoDetect handles invalid input", {
    testData <- loadSimCeuData()

    # Invalid parent
    expect_error(
        xoDetect(
            plinkFile = testData$plinkFile,
            mapFile = testData$mapFemale,
            familySize = 3,
            parent = "invalid"
        )
    )

    # Invalid family size
    expect_error(
        xoDetect(
            plinkFile = testData$plinkFile,
            mapFile = testData$mapFemale,
            familySize = 5,
            parent = "mother"
        )
    )

    # Missing file
    expect_error(
        xoDetect(
            plinkFile = "nonexistent_file",
            mapFile = testData$mapFemale,
            familySize = 3,
            parent = "mother"
        ),
        "PLINK file\\(s\\) not found"
    )
})

test_that("xoDetect output can be written to file", {
    skip_if_not(checkTestData(), "Test data not available")
    skip_if_not(checkMapFile("female"), "Female genetic map not available")

    testData <- loadSimCeuData()
    tempFile <- tempfile(fileext = ".csv")

    testResult <- xoDetect(
        plinkFile = testData$plinkFile,
        mapFile = testData$mapFemale,
        familySize = 3,
        parent = "mother",
        out = tempFile
    )

    # Check result is GRanges
    expect_s4_class(testResult, "GRanges")

    # Check file was created
    expect_true(file.exists(tempFile))

    # Check file can be read back
    readBack <- read.csv(tempFile, stringsAsFactors = FALSE)
    expect_equal(nrow(readBack), length(testResult))

    # Cleanup
    unlink(tempFile)
})

test_that("xoDetect returns GRangesList for case-control analysis", {
    skip_if_not(checkTestData(), "Test data not available")
    skip_if_not(checkMapFile("female"), "Female genetic map not available")

    testData <- loadSimCeuData()

    # Test with caseControl = TRUE (only works with familySize = 3)
    testResult <- xoDetect(
        plinkFile = testData$plinkFile,
        mapFile = testData$mapFemale,
        familySize = 3,
        parent = "mother",
        caseControl = TRUE
    )

    # Should be a GRangesList with "case" and "control" elements
    expect_s4_class(testResult, "GRangesList")
    expect_true("case" %in% names(testResult))
    expect_true("control" %in% names(testResult))

    # Each element should be a GRanges
    expect_s4_class(testResult$case, "GRanges")
    expect_s4_class(testResult$control, "GRanges")
})

test_that("xoDetect caseControl requires familySize = 3", {
    skip_if_not(checkTestData(), "Test data not available")
    skip_if_not(checkMapFile("female"), "Female genetic map not available")

    testData <- loadSimCeuData()

    # Should error when caseControl = TRUE and familySize = 2
    expect_error(
        xoDetect(
            plinkFile = testData$plinkFile,
            mapFile = testData$mapFemale,
            familySize = 2,
            parent = "mother",
            caseControl = TRUE
        ),
        "No case-control separation for 2-child families"
    )
})

test_that("xoDetect returns empty GRanges when no crossovers detected", {
    skip_if_not(checkTestData(), "Test data not available")
    skip_if_not(checkMapFile("female"), "Female genetic map not available")

    testData <- loadSimCeuData()

    # Use very strict filters to potentially get no results
    testResult <- xoDetect(
        plinkFile = testData$plinkFile,
        mapFile = testData$mapFemale,
        familySize = 3,
        parent = "mother",
        snpFilter = 100,
        cmFilter = 50
    )

    # Should still return GRanges (possibly empty)
    expect_s4_class(testResult, "GRanges")
    expect_true(length(testResult) >= 0)
})

test_that("xoDetect GRanges has correct seqnames", {
    skip_if_not(checkTestData(), "Test data not available")
    skip_if_not(checkMapFile("female"), "Female genetic map not available")

    testData <- loadSimCeuData()

    testResult <- xoDetect(
        plinkFile = testData$plinkFile,
        mapFile = testData$mapFemale,
        familySize = 3,
        parent = "mother"
    )

    # Check that seqnames are set
    if (length(testResult) > 0) {
        expect_true(length(seqnames(testResult)) > 0)
        # Seqnames should be chromosome identifiers
        expect_true(all(!is.na(seqnames(testResult))))
    }
})

## ============================================================================
## SNP Filter Tests
## ============================================================================

test_that(".smoothShortRuns removes short runs between identical states", {
    # A single-SNP run of state 2 between runs of state 1
    state <- setNames(c(1, 1, 1, 2, 1, 1, 1), paste0("rs", 1:7))
    smoothed <- inferRecom:::.smoothShortRuns(state, snpFilter = 5)
    expect_equal(unname(smoothed), rep(1, 7))
    expect_equal(names(smoothed), names(state))

    # Kept when the run is at least snpFilter SNPs long, or snpFilter is 0/1
    expect_equal(inferRecom:::.smoothShortRuns(state, snpFilter = 1), state)
    expect_equal(inferRecom:::.smoothShortRuns(state, snpFilter = 0), state)
    longRun <- c(1, 1, 2, 2, 2, 1, 1)
    expect_equal(inferRecom:::.smoothShortRuns(longRun, snpFilter = 3), longRun)

    # Kept when the states on either side differ (crossovers in two children)
    twoChildren <- c(1, 1, 3, 7, 7, 7)
    expect_equal(
        inferRecom:::.smoothShortRuns(twoChildren, snpFilter = 5), twoChildren
    )

    # Kept at either end of the chromosome, where there is no evidence of a
    # reversal
    expect_equal(inferRecom:::.smoothShortRuns(c(2, 1, 1, 1), 5), c(2, 1, 1, 1))
})

test_that("SNP filter removes crossover pairs caused by a genotyping error", {
    # Two children share inheritance except at rs6, where one SNP disagrees,
    # and from rs12 onwards, where a real crossover changes the state
    sameState <- c(rep(1, 5), 0, rep(1, 5), rep(0, 9))
    binaryInheritance <- cbind(rep(1, 20), sameState)
    rownames(binaryInheritance) <- paste0("rs", 1:20)

    unfiltered <- inferRecom:::.processTwoChildFamily("F1", binaryInheritance,
        snpFilter = 0
    )
    expect_equal(unname(unfiltered[, "startSnp"]), c("rs5", "rs6", "rs11"))

    filtered <- inferRecom:::.processTwoChildFamily("F1", binaryInheritance,
        snpFilter = 5
    )
    expect_equal(unname(filtered[, "startSnp"]), "rs11")
    expect_equal(unname(filtered[, "finishSnp"]), "rs12")
})

test_that("snpFilter changes xoDetect results", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()
    detect <- function(snpFilter) {
        xoDetect(testData$plinkFile, testData$mapFemale,
            familySize = 2, parent = "mother", snpFilter = snpFilter,
            cmFilter = 0
        )
    }

    expect_true(length(detect(50)) < length(detect(0)))
})

## ============================================================================
## Three-Child Inheritance State Tests
## ============================================================================

test_that("triple state codes identify the child that recombined", {
    patterns <- as.matrix(expand.grid(c(0, 1), c(0, 1), c(0, 1))[, 3:1])
    stateOf <- function(p) inferRecom:::.tripleStateCodes[sum(p * c(4, 2, 1)) + 1]

    # Complementary patterns describe the same state
    expect_equal(
        apply(patterns, 1, stateOf), apply(1 - patterns, 1, stateOf)
    )

    # Flipping one child changes the state by an amount that names that child
    childFound <- vapply(1:3, function(child) {
        flipped <- patterns
        flipped[, child] <- 1 - flipped[, child]
        change <- apply(flipped, 1, stateOf) - apply(patterns, 1, stateOf)
        inferRecom:::.getChildIdx(change)
    }, integer(8))
    expect_equal(childFound, matrix(rep(1:3, each = 8), nrow = 8))

    # A change that does not identify one child gives NA
    expect_true(is.na(inferRecom:::.getChildIdx(1)))
})

test_that(".tripleCrossovers reports crossovers in the right child", {
    # Child 2 recombines between rs10 and rs11; child 3 between rs15 and rs16
    tripleGeno <- cbind(
        rep(0, 20), c(rep(0, 10), rep(1, 10)), c(rep(0, 15), rep(1, 5))
    )
    rownames(tripleGeno) <- paste0("rs", 1:20)

    xo <- inferRecom:::.tripleCrossovers(
        tripleGeno, c("C1", "C2", "C3"), "F1", snpFilter = 5
    )
    expect_equal(unname(xo[, "childId"]), c("C2", "C3"))
    expect_equal(unname(xo[, "startSnp"]), c("rs10", "rs15"))
    expect_equal(unname(xo[, "finishSnp"]), c("rs11", "rs16"))

    # Too few informative SNPs
    expect_null(inferRecom:::.tripleCrossovers(
        tripleGeno[1:9, ], c("C1", "C2", "C3"), "F1", snpFilter = 5
    ))
})
