## Cellcano Tutorial

**Authors:** Wenjing Ma (wenjing.ma@emory.edu), Dr. Hao Wu (hao.wu@emory.edu), Emory University

**Latest revision:** 05-30-2022

In this tutorial, we will guide you through Cellcano, a supervised cell type annotation (celltyping) tool in scATAC-seq. We will use human Peripheral Blood Mononuclear Cells (PBMC) datasets as examples. 

![img](workflow.png)

We are aiming at accurate and efficient celltyping in single cell genomics, including scRNA-seq, scATAC-seq, scBS-seq, etc. If you are interested, stay tuned!

### Installation

To use Cellcano, you need to first install Python and we recommend Python 3.8. In addition, Cellcano embeds [ArchR](https://www.archrproject.com/) to process raw scATAC-seq data (fragment files or bam files). If you want to use the `preprocess` option in Cellcano to convert the raw scATAC-seq data to a gene score matrix, an R environment along with the installation of ArchR package is a must.



If you already have your R and Python installed, you can simply use `pip install Cellcano` to install our package. 



Otherwise, we recommend [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands) to easily install our package. Once your conda is installed, the following commands can be used to create the environment.

```shell
## create a conda environment with Python installed
conda create -n Cellcano python=3.8.12
## activate the created Cellcano environment
conda activate Cellcano

## install R environment
## If you already have R and devtools installed, you can ignore the conda install
conda install -c conda-forge r-base=4.1.1
conda install -c conda-forge r-devtools
# conda install -c bioconda bioconductor-rhdf5lib  ## you might need this if installing ArchR has a rhdf5lib error

## install Cellcano package
pip install Cellcano

## Use Cellcano help package to check
Cellcano -h

## Then you will see the following console printout
usage: Cellcano [-h] {preprocess,train,predict} ...

Cellcano: a supervised celltyping pipeline for single-cell omics.

positional arguments:
  {preprocess,train,predict}
                        sub-command help.
    preprocess          Run ArchR to preprocess raw input data (*fragments.tsv.gz, *.bam) to gene score.
    train               Train a Cellcano model.
    predict             Use Cellcano model to predict cell types.

optional arguments:
  -h, --help            show this help message and exit
```



If you need ArchR to process your raw scATAC-seq data to gene score matrix, please also install the ArchR package. You can either enter in R console by typing `R` or in your Rstudio:

```R
## install ArchR according to: https://www.archrproject.com/
## install BiocManager
install.packages("BiocManager")
## install ArchR package
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
library(ArchR)
ArchR::installExtraPackages()
```

*****

### 1. Download and process sample scATAC-seq raw data

#### Reference dataset
First, we take a healthy donor `PBMC_Rep1` from [GSE129785](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129785) to be our reference dataset. This study measures the single-cell chromatin landscape of human immune cell development.

We start by downloading the raw data along with the processed metadata. The raw data is provided on GEO under accession number GSE129785. The metadata needs some preprocessing procedures. The original metadata can be first downloaded from GEO at [GSE129785 GSE129785_scATAC-Hematopoiesis-All.cell_barcodes.txt.gz](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129785/suppl/GSE129785_scATAC-Hematopoiesis-All.cell_barcodes.txt.gz). We extract out individuals related to PBMC using the `Group` column and extract out those cells lying in `PBMC_Rep1` category. Then, we annotate the cluster information according to the [WashU EpiGenome Browser](http://epigenomegateway.wustl.edu/legacy/?genome=hg19&session=HcbHMSgBCc&statusId=28207718) under the `SatPathy, Granja et al 2019 v1` track. For your convenience, we provide the processed metadata with curated major human PBMC cell types on [Dropbox](https://www.dropbox.com/s/d7e7wcb4nuy5l2q/GSM3722015_PBMC_Rep1_curated_metadata.csv). 

```shell
# create and enter directory
mkdir train_data
cd train_data
# Download sample fragment file
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3722nnn/GSM3722015/suppl/GSM3722015_PBMC_Rep1_fragments.tsv.gz
# Download processed cell metadata
wget https://www.dropbox.com/s/d7e7wcb4nuy5l2q/GSM3722015_PBMC_Rep1_curated_metadata.csv
```

Next, we use Cellcano to call ArchR to process the fragment file into gene score matrix. The `-i` option scans files ending with `fragments.tsv.gz` or `.bam` and load the files to ArchR.`-g` describes which genome will be used. `--threads` indicates how many threads will be used to run ArchR. ArchR would take around 30 mins to process the data with threads as 4. 

```shell
Cellcano preprocess -i train_data -o train_data -g hg19 --threads 4

## The output will be..
Output will be written to train_data
R[write to console]:
                                                   / |
                                                 /    \
            .                                  /      |.
            \\\                              /        |.
              \\\                          /           `|.
                \\\                      /              |.
                  \                    /                |\
                  \\#####\           /                  ||
                ==###########>      /                   ||
                 \\##==......\    /                     ||
            ______ =       =|__ /__                     ||      \\\
        ,--' ,----`-,__ ___/'  --,-`-===================##========>
       \               '        ##_______ _____ ,--,__,=##,__   ///
        ,    __==    ___,-,__,--'#'  ==='      `-'    | ##,-/
        -,____,---'       \\####\\________________,--\\_##,/
           ___      .______        ______  __    __  .______
          /   \     |   _  \      /      ||  |  |  | |   _  \
         /  ^  \    |  |_)  |    |  ,----'|  |__|  | |  |_)  |
        /  /_\  \   |      /     |  |     |   __   | |      /
       /  _____  \  |  |\  \\___ |  `----.|  |  |  | |  |\  \\___.
      /__/     \__\ | _| `._____| \______||__|  |__| | _| `._____|


R[write to console]: Setting default number of Parallel threads to 36.

R[write to console]: Setting default number of Parallel threads to 4.

R[write to console]: Setting default genome to Hg19.

...

R[write to console]: 2022-05-15 10:36:13 : (GSM3722015_PBMC_Rep1 : 1 of 1) Finished Creating Arrow File, 31.809 mins elapsed.

...

R[write to console]: 2022-05-15 10:37:39 : Organizing colData, 1.336 mins elapsed.

R[write to console]: 2022-05-15 10:37:39 : Organizing rowData, 1.336 mins elapsed.

R[write to console]: 2022-05-15 10:37:39 : Organizing rowRanges, 1.336 mins elapsed.

R[write to console]: 2022-05-15 10:37:39 : Organizing Assays (1 of 1), 1.336 mins elapsed.

R[write to console]: 2022-05-15 10:37:39 : Constructing SummarizedExperiment, 1.337 mins elapsed.

R[write to console]: 2022-05-15 10:37:40 : Finished Matrix Creation, 1.359 mins elapsed.
```

The gene score matrix will be in the `train_data` folder:

```
- ArchR_genescore.mtx.gz: stores the gene scores in a COO format.
- ArchR_genescore_barcodes.tsv: stores the cell barcodes.
- ArchR_genescore_genes.tsv: stores the gene symbols.
```

#### Target dataset

Here, as an example, we use a FACS-sorted human PBMC scATAC-seq dataset to be our target dataset. The data is provided on [GSE123578](GSE123578). The fragments files are having the suffixes as `fragments.tsv.gz`. 5 human PBMC cell types available: CD19 B cells, CD4 T cells, CD8 T cells, Monocytes and NK cells. In total, 21K cells are available. We downsample 5K cells out of them and provide the gene score matrix on [Dropbox](https://www.dropbox.com/s/e7g9vem3oxt096l/FACS5K_genescore.csv).

```shell
## create and enter directory
mkdir test_data
cd test_data

## download processed target data
wget https://www.dropbox.com/s/e7g9vem3oxt096l/FACS5K_genescore.csv
```

### 2. Pre-train Cellcano's first round

With the preprocessed reference gene score matrix, we can use it as input for training a first-round prediction. When the input data ends with `.mtx`, the `-i` option takes the common file prefix indicating the matrix, barcodes and genes.

```shell
## run the training process
Cellcano train -i train_data/ArchR_genescore -m train_data/GSM3722015_PBMC_Rep1_curated_metadata.csv -o output --prefix PBMC_Rep1_trained 

## The output will be
Skipping registering GPU devices...
_utils.<module>: INFO: Num GPUs Available: 0
train.load_train_adata: INFO: Loading data...
 This may take a while depending on your data size..
 
 ...
 
601/601 - 1s - loss: 0.3551 - accuracy: 0.8839
Epoch 49/100
601/601 - 1s - loss: 0.3499 - accuracy: 0.8833
Epoch 50/100
601/601 - 1s - loss: 0.3429 - accuracy: 0.8836
Epoch 51/100
601/601 - 1s - loss: 0.3457 - accuracy: 0.8830
Epoch 52/100
601/601 - 1s - loss: 0.3442 - accuracy: 0.8839
2022-05-15 17:26:51.405652: W tensorflow/python/util/util.cc:348] Sets are not currently considered sequences, but this may change in the future, so consider avoiding using them.
```

The trained model will be under the `output` directory with a folder named `PBMC_Rep1_trained`. Along with the trained model, the tSNE visualization and an [anndata](https://anndata.readthedocs.io/en/latest/) object storing the reference data are stored.

### 3. Predict on target dataset in second round

We use the trained model to predict on the target dataset. When the input data ends with `.csv`, we use the full path as input. 

```shell
## run the prediction process
Cellcano predict -i test_data/FACS5K_genescore.csv --trained_model output/PBMC_Rep1_trainedMLP_model -o output --prefix predict_FACS --tworound

## The output will be
Skipping registering GPU devices...
_utils.<module>: INFO: Num GPUs Available: 0
predict.predict: INFO: Loading data...
 This may take a while depending on your data size..
 
...

Epoch 28/30
63/63 - 0s - accuracy: 0.8685 - student_loss: 0.1856 - distillation_loss: 0.1795
Epoch 29/30
63/63 - 0s - accuracy: 0.8779 - student_loss: 0.6145 - distillation_loss: 0.2452
Epoch 30/30
63/63 - 0s - accuracy: 0.8685 - student_loss: 0.4724 - distillation_loss: 0.2203
=== Predicted celltypes:  {'CD8 T cells', 'CD4 T cells', 'B cells', 'NK cells', 'Monocytes'}
```

The result will be stored in a csv file under the output directory indicated by `-o` option with  specified prefix. In our case, it will be `predict_FACScelltypes.csv`. If using without `--tworound` option, the result will be directly predicted from the trained model.

### 4. Result analysis

We use Sankey plot to analyze the cell type assignment. Here, we use R to run the following analysis code.

```R
# We need ggsankey package to draw the Sankey plot
# install.packages("devtools")
# devtools::install_github("davidsjoberg/ggsankey")

library(ggsankey)

library(ggplot2)
library(dplyr)
## load predicted results
pred_df = read.csv("output/predict_FACScelltypes.csv", header=T, row.names=1)
## add cell type labels from FACS barcodes
pred_df$celltype = sapply(strsplit(pred_df$barcode, '_'), '[', 3) 
sankey_df = pred_df %>% make_long(celltype, pred_celltype)
## count how many cells in the link
count_df = sankey_df %>% dplyr::group_by(node)%>% tally()  
df = merge(sankey_df, count_df, by.x = 'node', by.y = 'node', all.x = TRUE)
g = ggplot(df, aes(x=x, next_x=next_x, node=node, next_node=next_node,
                          fill=factor(node), label=paste0(node," n=", n)))
g = g + geom_sankey(flow.alpha=0.5, node.color="black", show.legend=T)
g = g + geom_sankey_label(size=3, color="black", fill="white", hjust=-0.5)
g = g + theme_bw()
g = g + theme(legend.position="none")
g = g + theme(axis.title=element_blank(), axis.text.y = element_blank(),
              axis.ticks = element_blank(), panel.grid = element_blank())
g
```

![img](predict_FACS_sankey.png)

In the Sankey plot, the left nodes are the labels from FACS and the right nodes are predicted cell types. Most of the cells are correctly predicted while some cells in CD8 T cells are wrongly assigned as CD4 T cells and some cells in NK cells are wrongly assigned to CD8 T cells. The reason can be indicated from the tSNE visualization on FACS dataset below. The tSNE plots show the high-dimensional data in 2-dimension and labeld with FACS cell types. There are cell mixtures between NK cells and CD8 T cells as well as CD4 T cells and CD8 T cells. The mixture makes the prediction much harder. 

![img](FACS_tSNE.png)

