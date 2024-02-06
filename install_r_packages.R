# Please note there is also a package via Anaconda: https://anaconda.org/bioconda/bioconductor-clusterprofiler
# Additionally, please be aware you may want to have a dedicated environment if
# you install via Anaconda because of conflicts with various dependencies

# Install igraph: https://github.com/igraph/rigraph
options(
  repos = c(
    igraph = 'https://igraph.r-universe.dev',
    CRAN = 'https://cloud.r-project.org'
  )
)
install.packages('igraph')

# Installation instructions for Bioconductor: https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html#installing-biocmanager
install.packages("BiocManager", repos = "https://cloud.r-project.org")

BiocManager::install()  # update Bioconductor packages
BiocManager::install("clusterProfiler", force=TRUE)
