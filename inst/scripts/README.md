# How the example data in `inst/extdata` were produced

The files in `inst/extdata` are used in the examples, tests and vignette of
inferRecom. The scripts in this directory document how each file was made.
They are not run when the package is built or checked.

The genotype data are **simulated for package development and examples
only**. They do not reflect any biological principle; see the note on the
recombination map below.

| File | Contents | Script |
|---|---|---|
| `simCEU.bed`, `simCEU.bim`, `simCEU.fam` | Simulated pedigrees in PLINK binary format | `make-simCEU.R` |
| `female_chr4.txt`, `male_chr4.txt`, `sexavg_chr4.txt` | Genetic maps for chromosome 4 | `make-genetic-maps.R` |
| `xoPat.rds`, `xoMat.rds` | Example `xoDetect()` results | `make-crossovers.R` |

## `simCEU`: simulated pedigrees

**Contents:** 8 nuclear families with 40 individuals in total: families F1–F3
have 2 children, F4–F6 have 3, F7 has 4 and F8 has 5. There are 1162 SNPs on
chromosome 4 between 78,002,267 and 79,298,781 bp (GRCh37), all with rs IDs.

**Source genotypes:** HapMap CEU genotypes (2010 release) for chromosome 4, as
a VCF file (`CEU_chr4.vcf`).

**Simulation:** with the R package sim1000G (Dimitromanolakis et al. 2019),
following its vignette
[*Simulating family data*](https://adimitromanolakis.github.io/sim1000G/inst/doc/SimulatingFamilyData.html):

1. `readVCF()` reads the VCF, keeping up to 5000 variants with minor allele
   frequency at least 0.02.
2. `startSimulation()` sets up the simulation for up to 100 individuals.
3. For each family size, the simulation is reset and the recombination map
   is set to **a uniform 100 cM across the region**
   (`SIM$cm <- seq(0, 100, length = length(SIM$cm))`).
4. `newFamilyWithOffspring()` creates each family. Two unrelated founders are
   sampled from the CEU haplotypes, and each child is produced by
   `SIM$mate()`, which simulates one meiosis per parent. Each child is given
   a random sex.
5. `writePED()` writes each group of families as PLINK PED/MAP files.
6. The PED/MAP files are combined and written in PLINK binary format with
   snpStats, which is equivalent to `plink --merge` and `plink --make-bed`.
   Individuals are given readable IDs: `F<n>P1` is the father, `F<n>P2` the
   mother and `F<n>C<k>` child *k*. One child per family is chosen at random
   to be affected (PLINK phenotype 2); everyone else is unaffected (1).
   Allele codes are converted to the VCF's nucleotides.

**About the recombination map:** 100 cM across about 1.3 Mb is far above
real recombination rates, which are roughly 1 cM/Mb on average. It was
chosen so that this short region contains many crossovers for the examples.
The map is uniform, so the simulated crossovers do not follow the
recombination hotspots or the female and male differences of the genetic
maps below, although those maps are used to analyse the data.

**Reproducibility:** the bundled files were made by an earlier run of this
procedure without a fixed random seed. Its relabelling step (step 6) was not
recorded and is reconstructed in `make-simCEU.R`. Running the script, which
sets a seed, produces data with the same structure but different genotypes.
Two other differences remain:

- In the bundled `simCEU.bim`, alleles are written in quotes and the genetic
  distance column is 0.
- In the bundled `simCEU.fam`, cases were originally coded 1 and everyone
  else 0. They were recoded to standard PLINK coding (cases 2, others 1) in
  inferRecom 0.99.3.

## Genetic maps

**Contents:** the female, male and sex-averaged refined European
recombination maps of Bhérer, Campbell and Auton (2017), for chromosome 4 in
genome build 37 (GRCh37). Each file is tab-delimited with columns `chr`,
`pos` (bp), `rate` (cM/Mb) and `cM` (cumulative genetic position).

**Source:** the files `female_chr4.txt`, `male_chr4.txt` and
`sexavg_chr4.txt` of `Refined_EUR_genetic_map_b37.tar.gz` at
<https://github.com/cbherer/Bherer_etal_SexualDimorphismRecombination>,
unchanged. `make-genetic-maps.R` downloads them, and the downloaded files
are identical to the bundled ones.

## Example crossovers

**Contents:** `xoPat.rds` and `xoMat.rds` hold the paternal and maternal
crossovers detected by `xoDetect()` in the simCEU families with three or
more children, using the male and female maps and the default filters. They
are used in the `xoPhase()` examples and tests.

**Source:** `make-crossovers.R`.

## References

- Bhérer C, Campbell CL, Auton A (2017). Refined genetic maps reveal sexual
  dimorphism in human meiotic recombination at multiple scales. *Nature
  Communications* 8, 14994.
- Dimitromanolakis A, Xu J, Krol A, Briollais L (2019). sim1000G: a
  user-friendly genetic variant simulator in R for unrelated individuals and
  family-based designs. *BMC Bioinformatics* 20, 26.
- The International HapMap 3 Consortium (2010). Integrating common and rare
  genetic variation in diverse human populations. *Nature* 467, 52–58.
