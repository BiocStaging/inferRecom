## Package-wide imports. All external functions are imported here so that
## code in the rest of the package calls them without a namespace prefix.

#' @importFrom BiocParallel bplapply MulticoreParam SerialParam
#' @importFrom Biostrings DNAStringSet DNAStringSetList
#' @importFrom dplyr distinct
#' @importFrom GenomicRanges as.data.frame GRanges GRangesList
#'     makeGRangesFromDataFrame seqnames
#' @importFrom IRanges DataFrameList IRanges
#' @importFrom methods as is
#' @importFrom S4Vectors DataFrame SimpleList split
#' @importFrom snpStats read.plink
#' @importFrom stats setNames
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom utils read.delim write.csv
#' @importFrom VariantAnnotation VCF VCFHeader
NULL
