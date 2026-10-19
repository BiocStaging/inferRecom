#' Example PLINK Dataset - simCEU
#'
#' Simulated pedigrees based on HapMap CEU (Utah residents with Northern and
#' Western European ancestry) genotypes, for demonstrating crossover
#' detection, phasing and ROH analysis. These data are for package
#' development and examples only and do not reflect any biological
#' principle: recombination was simulated with a uniform map of 100 cM
#' across the region, far above real recombination rates, so that the short
#' region contains crossovers.
#'
#' @format PLINK binary format files (.bed, .bim, .fam) containing:
#' \describe{
#'   \item{Families}{8 nuclear families, 40 individuals: three families with
#'     2 children, three with 3, one with 4 and one with 5}
#'   \item{SNPs}{1162 SNPs on chromosome 4 between 78.0 and 79.3 Mb (GRCh37)}
#'   \item{Genotypes}{Founder haplotypes sampled from HapMap CEU genotypes,
#'     with recombination simulated in each meiosis}
#' }
#'
#' @details
#' The dataset includes three PLINK binary files:
#' \itemize{
#'   \item \code{simCEU.bed} - Binary genotype data
#'   \item \code{simCEU.bim} - Variant information (SNP IDs, positions, alleles)
#'   \item \code{simCEU.fam} - Sample information (family IDs, relationships,
#'     sex and affected status). One child per family is affected (PLINK
#'     phenotype 2); everyone else is unaffected (1).
#' }
#'
#' Access the files using:
#' \preformatted{
#' dataPath <- system.file("extdata", package = "inferRecom")
#' plinkFile <- file.path(dataPath, "simCEU")
#' }
#'
#' @source Simulated with the R package \code{sim1000G} from HapMap CEU
#' genotypes (2010 release); see \code{inst/scripts/make-simCEU.R} and
#' \code{inst/scripts/README.md}.
#'
#' @return No R object is returned. Access files via
#' \code{system.file("extdata", package = "inferRecom")}.
#'
#' @references
#' The International HapMap 3 Consortium (2010). Integrating common and rare
#' genetic variation in diverse human populations. Nature, 467, 52-58.
#'
#' Dimitromanolakis, A., Xu, J., Krol, A., & Briollais, L. (2019). sim1000G:
#' a user-friendly genetic variant simulator in R for unrelated individuals and
#' family-based designs. BMC bioinformatics, 20(1), 26.
#' @name simCEU-data
NULL
#'
#' Genetic Maps for Chromosome 4
#'
#' Female, male and sex-averaged genetic maps for chromosome 4
#' (\code{female_chr4.txt}, \code{male_chr4.txt}, \code{sexavg_chr4.txt}).
#'
#' @format Tab-delimited text file with columns:
#' \describe{
#'   \item{chr}{Chromosome identifier (chr4)}
#'   \item{pos}{Physical position in base pairs (hg19/GRCh37)}
#'   \item{rate}{Recombination rate (cM/Mb) at this position}
#'   \item{cM}{Cumulative genetic distance in centiMorgans from chromosome
#'     start}
#' }
#'
#' Access the files using:
#' \preformatted{
#' femaleMap <- read.delim(file.path(dataPath, "female_chr4.txt"))
#' }
#'
#' @source The chromosome 4 files of the refined European maps of Bherer, et
#' al. (2017), unchanged; see \code{inst/scripts/make-genetic-maps.R}.
#'
#' @return No R object is returned. Access files via
#' \code{system.file("extdata", package = "inferRecom")}.
#'
#' @references
#' Bhérer, C., Campbell, C. L., & Auton, A. (2017). Refined genetic maps reveal
#' sexual dimorphism in human meiotic recombination at multiple scales. Nature
#' communications, 8(1), 14994.
#' @name geneticMaps-data
NULL
#' Example Crossover Data
#'
#' Pre-computed maternal and paternal crossover events from the simCEU dataset.
#' Used for testing and demonstrating haplotype phasing functionality.
#'
#' @format \link[GenomicRanges:GRanges-class]{GRanges} objects stored as RDS
#' files with metadata columns:
#' \describe{
#'   \item{childId}{Child identifier}
#'   \item{familyId}{Family identifier}
#'   \item{startSnp}{SNP name at crossover interval start}
#'   \item{finishSnp}{SNP name at crossover interval end}
#'   \item{startCm}{Genetic position at interval start (cM)}
#'   \item{finishCm}{Genetic position at interval end (cM)}
#' }
#'
#' @details
#' Two crossover datasets are provided:
#' \itemize{
#'   \item \code{xoMat.rds} - Maternal crossover events
#'   \item \code{xoPat.rds} - Paternal crossover events
#' }
#'
#' These objects contain detected crossover events from the families with
#' three or more children in the simCEU dataset. Crossovers were identified
#' using:
#'
#' \preformatted{
#' # Maternal crossovers
#' xoMat <- xoDetect(
#'   plinkFile = "simCEU",
#'   mapFile = "female_chr4.txt",
#'   familySize = 3,
#'   parent = "mother",
#'   snpFilter = 5,
#'   cmFilter = 1
#' )
#'
#' # Paternal crossovers
#' xoPat <- xoDetect(
#'   plinkFile = "simCEU",
#'   mapFile = "male_chr4.txt",
#'   familySize = 3,
#'   parent = "father",
#'   snpFilter = 5,
#'   cmFilter = 1
#' )
#' }
#'
#' Access the crossover data using:
#' \preformatted{
#' dataPath <- system.file("extdata", package = "inferRecom")
#' xoMat <- readRDS(file.path(dataPath, "xoMat.rds"))
#' xoPat <- readRDS(file.path(dataPath, "xoPat.rds"))
#' }
#'
#' @source Created by \code{inst/scripts/make-crossovers.R}.
#'
#' @name crossoverData
#'
#' @return No R object is returned. Access files via
#' \code{system.file("extdata", package = "inferRecom")}.
NULL
