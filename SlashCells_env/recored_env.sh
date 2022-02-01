python=3.6
conda install -c conda-forge r-base=4.1.1
conda install -c conda-forge r-devtools
conda install -c conda-forge r-cairo
conda install -c conda-forge tensorflow=1.15
#conda install -c conda-forge tensorflow-gpu=1.15

## python package
pip install rpy2
pip install ipython
pip install keras==2.1.0
pip install matplotlib==3.3.4
pip install anndata==0.7.4
pip install scanpy==1.5.0

## R package
BiocManager
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
library(ArchR)
ArchR::installExtraPackages()

BiocManager::install("BSgenome.Mmusculus.UCSC.mm9")
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

