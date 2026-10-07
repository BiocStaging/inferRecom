## Package-wide imports. All external functions are imported here so that
## code in the rest of the package calls them without a namespace prefix.

#' @importFrom BiocParallel bplapply MulticoreParam SerialParam
#' @importFrom dplyr distinct
#' @importFrom GenomicRanges end GRanges GRangesList makeGRangesFromDataFrame
#'     mcols seqnames start strand width
#' @importFrom methods as is
#' @importFrom S4Vectors DataFrame
#' @importFrom snpStats read.plink
#' @importFrom stats setNames
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom utils read.delim write.csv write.table
NULL
