## Create the example crossover results in inst/extdata/{xoPat,xoMat}.rds
##
## Paternal and maternal crossovers in the simCEU families with three or
## more children, detected by xoDetect() with the sex-specific genetic maps
## and default filters. They are used in the xoPhase() examples and tests.

library(inferRecom)

dataPath <- system.file("extdata", package = "inferRecom")
plinkFile <- file.path(dataPath, "simCEU")
outDir <- "."

xoPat <- xoDetect(
    plinkFile = plinkFile,
    mapFile = file.path(dataPath, "male_chr4.txt"),
    familySize = 3,
    parent = "father"
)
xoMat <- xoDetect(
    plinkFile = plinkFile,
    mapFile = file.path(dataPath, "female_chr4.txt"),
    familySize = 3,
    parent = "mother"
)

saveRDS(xoPat, file.path(outDir, "xoPat.rds"))
saveRDS(xoMat, file.path(outDir, "xoMat.rds"))
