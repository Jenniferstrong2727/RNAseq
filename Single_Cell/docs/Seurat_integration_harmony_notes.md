R/4.2.0
(base) [stronj04@li04e02 scripts]$ module load R/4.2.0
(base) [stronj04@li04e02 scripts]$ R
R version 4.2.0 (2022-04-22) -- "Vigorous Calisthenics"
Copyright (C) 2022 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)
R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.
  Natural language support but running in an English locale
R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.
Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.
[Previously saved workspace restored]
> # 00_check_packages.R
required_pkgs <- c(
  "Seurat",
  "scDblFinder",
  "scater",
  "SingleCellExperiment",
  "rtracklayer",
  "harmony",
  "future",
  "cowplot",
  "readxl",
  "ggplot2"
)
cat("R version:\n")
print(R.version.string)
cat("\nChecking libraries...\n\n")
for (p in required_pkgs) {
  cat("----", p, "----\n")
  suppressPackageStartupMessages(
    library(p, character.only = TRUE)
  )
  cat("  version:", as.character(packageVersion(p)), "\n\n")
}
cat("All requested packages loaded successfully.\n")
R version:
[1] "R version 4.2.0 (2022-04-22)"
Checking libraries...
---- Seurat ----
  version: 5.1.0 
---- scDblFinder ----
  version: 1.12.0 
---- scater ----
  version: 1.32.1 
---- SingleCellExperiment ----
  version: 1.20.1 
---- rtracklayer ----
  version: 1.58.0 
---- harmony ----
  version: 1.2.1 
---- future ----
  version: 1.40.0 
---- cowplot ----
  version: 1.1.3 
---- readxl ----
  version: 1.4.3 
---- ggplot2 ----
  version: 3.5.1 


#if error try
#Seurat 5.3.1
#RcppAnnoy 0.0.22

