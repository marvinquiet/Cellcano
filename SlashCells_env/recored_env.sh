python=3.6
conda install -c conda-forge r-base=4.1.1
conda install -c conda-forge r-devtools
conda install -c conda-forge r-cairo

## python package
pip install rpy2
pip install ipython

## R package
BiocManager
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
library(ArchR)
ArchR::installExtraPackages()

BiocManager::install("BSgenome.Mmusculus.UCSC.mm9")
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

