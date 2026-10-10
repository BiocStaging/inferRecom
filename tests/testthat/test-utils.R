## ============================================================================
## Shared Internal Helper Tests
## ============================================================================

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

test_that(".asGRangesResult returns empty results of the same type", {
    emptyDf <- DataFrame(startPos = numeric(), finishPos = numeric())

    expect_equal(
        inferRecom:::.asGRangesResult(emptyDf, "chr4", FALSE, NULL, "childId"),
        GRanges()
    )

    emptyCc <- inferRecom:::.asGRangesResult(emptyDf, "chr4", TRUE, NULL, "childId")
    expect_s4_class(emptyCc, "GRangesList")
    expect_equal(names(emptyCc), c("case", "control"))
})
