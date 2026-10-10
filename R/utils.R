## ========================== SHARED INTERNAL HELPERS ===========================

## Remove a .bed/.bim/.fam extension so plinkFile is the PLINK prefix. Other
## extensions are kept because PLINK prefixes may contain dots (e.g. "cohort.v2")
.stripPlinkExtension <- function(plinkFile) {
    if (tolower(file_ext(plinkFile)) %in% c("bed", "bim", "fam")) {
        plinkFile <- file_path_sans_ext(plinkFile)
    }
    plinkFile
}
