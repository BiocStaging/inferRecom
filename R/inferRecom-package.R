## Package-wide imports. All external functions are imported here so that
## code in the rest of the package calls them without a namespace prefix.

#' @importFrom BiocParallel bplapply MulticoreParam SerialParam
#' @importFrom dplyr distinct
#' @importFrom GenomicRanges as.data.frame GRanges GRangesList
#'     makeGRangesFromDataFrame seqnames
#' @importFrom methods as is
#' @importFrom S4Vectors DataFrame split
#' @importFrom snpStats read.plink
#' @importFrom stats setNames
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom utils read.delim write.csv write.table
NULL
