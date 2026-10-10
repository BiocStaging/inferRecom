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
