# Installation instructions for Bioconductor: https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html#installing-biocmanager
install.packages("BiocManager", repos = "https://cloud.r-project.org")

BiocManager::install()  # update Bioconductor packages
BiocManager::install("clusterProfiler", force=TRUE)
