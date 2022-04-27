conda create --prefix ./Pyramid python=3.8.12
## python package
pip install tensorflow==2.4.1
pip install anndata==0.7.4
pip install scanpy==1.8.2
pip install rpy2
pip install ipython
pip install keras
pip install matplotlib

## R packages
conda install -c conda-forge r-base=4.1.1
conda install -c conda-forge r-devtools
conda install -c conda-forge r-cairo
#conda install -c conda-forge r-xml

install.packages("BiocManager")
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
library(ArchR)
ArchR::installExtraPackages()

#BiocManager::install("BSgenome.Mmusculus.UCSC.mm9")
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")


## Other bioinfo tools
conda install -c bioconda samtools
conda install -c bioconda macs2


