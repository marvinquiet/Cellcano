python=3.6
conda install -c conda-forge r-base=4.1.1
conda install -c conda-forge r-devtools
conda install -c conda-forge r-cairo

## python package
pip install rpy2

## R package
BiocManager
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
library(ArchR)
ArchR::installExtraPackages()

