## Simulate the example pedigrees in inst/extdata/simCEU.{bed,bim,fam}
##
## Eight nuclear families with simulated recombination, generated with the
## sim1000G package (Dimitromanolakis et al. 2019), following its vignette
## "Simulating family data":
## https://adimitromanolakis.github.io/sim1000G/inst/doc/SimulatingFamilyData.html
##
## These data are for package development and examples only. They do not
## reflect any biological principle: in particular, recombination is
## simulated with a uniform map of 100 cM across a ~1.3 Mb region, far above
## real recombination rates, so that the short region contains crossovers.
##
## The bundled files were made by an earlier run of this procedure without a
## fixed random seed, so a new run produces different genotypes with the same
## structure. See inst/scripts/README.md.
##
## Requires: sim1000G (CRAN) and snpStats (Bioconductor).

library(sim1000G)
library(snpStats)

## ---------------- Settings ----------------

## HapMap CEU genotypes (2010 release) for chromosome 4, GRCh37, as VCF
vcfFile <- "CEU_chr4.vcf"
outDir <- "."
outBase <- file.path(outDir, "simCEU")

## Number of children per family, for families 1 to 8
childrenPerFamily <- c(2, 2, 2, 3, 3, 3, 4, 5)

set.seed(1)

## ---------------- Functions ----------------

## Simulate families with sim1000G and write them as PLINK PED/MAP files.
## The simulation is reset first, with a uniform recombination map of
## 100 cM across the region. sim1000G's newFamilyWithOffspring() gives the
## parents IDs 1 (father) and 2 (mother), the children IDs 11, 12, ..., and
## each child a random sex. Returns the PED/MAP file prefix.
simulateFamilies <- function(vcf, familyIds, nOffspring, prefix) {
    SIM$reset()
    SIM$cm <- seq(0, 100, length = length(SIM$cm))

    fam <- do.call(rbind, lapply(familyIds, function(familyId) {
        newFamilyWithOffspring(familyId, noffspring = nOffspring)
    }))
    writePED(vcf, fam, prefix)

    prefix
}

## Read PED/MAP files with the same SNPs as one data set. The PED files are
## combined before reading, so that every individual's genotypes are coded
## against the same allele 1 and allele 2 (read.pedfile() chooses these per
## file)
readSimulatedPed <- function(prefixes) {
    combined <- tempfile(fileext = ".ped")
    on.exit(unlink(combined))
    writeLines(
        unlist(lapply(paste0(prefixes, ".ped"), readLines)),
        combined
    )

    read.pedfile(combined, snps = paste0(prefixes[1], ".map"))
}

## Give individuals readable IDs (F<family>P1 = father, F<family>P2 = mother,
## F<family>C<k> = child k) and mark one randomly chosen child per family as
## affected, using PLINK coding (2 = affected, 1 = unaffected)
relabelFamilies <- function(fam) {
    familyId <- paste0("F", fam$pedigree)
    newId <- function(id) {
        role <- ifelse(id == "1", "P1", ifelse(id == "2", "P2",
            paste0("C", suppressWarnings(as.integer(id)) - 10L)
        ))
        ifelse(is.na(id), NA_character_, paste0(familyId, role))
    }

    children <- which(!is.na(fam$father))
    cases <- vapply(split(children, familyId[children]), function(idx) {
        idx[sample.int(length(idx), 1L)]
    }, integer(1))
    affected <- rep(1L, nrow(fam))
    affected[cases] <- 2L

    data.frame(
        pedigree = familyId,
        member = newId(fam$member),
        father = newId(fam$father),
        mother = newId(fam$mother),
        sex = fam$sex,
        affected = affected
    )
}

## Nucleotides of allele 1 and allele 2 of each SNP. writePED() codes
## alleles as 1 (VCF REF) and 2 (VCF ALT), and read.pedfile() takes allele 1
## to be the first code it sees. A monomorphic SNP has no allele 2 in the
## data; it is given the other VCF allele.
snpAlleles <- function(map, vcf) {
    vcfIdx <- match(map$snp.names, vcf$vcf[, 3])
    ref <- as.character(vcf$vcf[vcfIdx, 4])
    alt <- as.character(vcf$vcf[vcfIdx, 5])
    nucleotide <- function(code) ifelse(code == "1", ref, alt)

    allele.1 <- nucleotide(map$allele.1)
    allele.2 <- ifelse(is.na(map$allele.2),
        ifelse(map$allele.1 == "1", alt, ref),
        nucleotide(map$allele.2)
    )
    data.frame(allele.1 = allele.1, allele.2 = allele.2)
}

## ---------------- Simulation ----------------

vcf <- readVCF(vcfFile,
    maxNumberOfVariants = 5000, min_maf = 0.02, max_maf = NA
)
readGeneticMap(chromosome = 4)
startSimulation(vcf, totalNumberOfIndividuals = 100)

## Families with the same number of children are simulated together, as in
## the original script; each group starts from a reset simulation
familyGroups <- split(seq_along(childrenPerFamily), childrenPerFamily)
prefixes <- vapply(names(familyGroups), function(nOffspring) {
    simulateFamilies(
        vcf, familyGroups[[nOffspring]], as.integer(nOffspring),
        file.path(outDir, paste0("sim", nOffspring))
    )
}, character(1))

## ---------------- PLINK binary files ----------------

## Equivalent to merging the PED/MAP files with plink --make-bed
sim <- readSimulatedPed(prefixes)
fam <- relabelFamilies(sim$fam)
alleles <- snpAlleles(sim$map, vcf)

write.plink(outBase,
    snps = sim$genotypes,
    pedigree = fam$pedigree, id = fam$member,
    father = fam$father, mother = fam$mother,
    sex = fam$sex, phenotype = fam$affected,
    chromosome = sim$map[, 1], genetic.distance = sim$map[, 3],
    position = sim$map[, 4],
    allele.1 = alleles$allele.1, allele.2 = alleles$allele.2
)
