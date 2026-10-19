## Download the genetic maps in inst/extdata/{female,male,sexavg}_chr4.txt
##
## Refined European sex-specific and sex-averaged recombination maps from
## Bherer, Campbell and Auton (2017), Nature Communications 8:14994, in
## genome build 37 (GRCh37):
## https://github.com/cbherer/Bherer_etal_SexualDimorphismRecombination
##
## The bundled files are the chromosome 4 files of the
## Refined_EUR_genetic_map_b37 archive, unchanged. Each is tab-delimited with
## columns chr, pos (bp), rate (cM/Mb) and cM (cumulative genetic position).

outDir <- "."
mapUrl <- paste0(
    "https://github.com/cbherer/Bherer_etal_SexualDimorphismRecombination/",
    "raw/master/Refined_EUR_genetic_map_b37.tar.gz"
)
mapFiles <- paste0(
    "Refined_EUR_genetic_map_b37/", c("female", "male", "sexavg"),
    "_chr4.txt"
)

archive <- tempfile(fileext = ".tar.gz")
download.file(mapUrl, archive, mode = "wb")

extractDir <- tempfile()
untar(archive, files = mapFiles, exdir = extractDir)
file.copy(file.path(extractDir, mapFiles), outDir, overwrite = TRUE)
