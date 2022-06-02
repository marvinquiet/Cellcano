## Cellcano Tutorial

**Authors:** [Wenjing Ma](https://marvinquiet.github.io/) (wenjing.ma@emory.edu), [Jiaying Lu](https://lujiaying.github.io/) (jiaying.lu@emory.edu), [Dr. Hao Wu](http://www.haowulab.org/) (hao.wu@emory.edu), Emory University

**Latest revision:** 06-01-2022

In this tutorial, we will guide you through Cellcano, a supervised cell type annotation (celltyping) tool in scATAC-seq. We will use human Peripheral Blood Mononuclear Cells (PBMC) datasets as examples. 

![img](workflow.png)

The above is a workflow figure illustrating how Cellcano works. Cellcano is a two-round algorithm which first uses a simple multi-layer perceptron (MLP) to predict confident cells in the target dataset and then in the second round, use those highly confident cells as guidance to predict remaining cells. This idea can be widely applicable to other single-cell genomics data, including scRNA-seq, scATAC-seq, scBS-seq, etc. If you are interested, stay tuned!

### Installation

Since Cellcano is a Python package, you will need to install Python and the recommended version is **Python 3.8**. In addition, Cellcano incorporates [ArchR](https://www.archrproject.com/) to do data preprocessing by turning raw scATAC-seq data (fragment files or bam files) into gene score matrices. Therefore, R environment and the [installation of ArchR](https://www.archrproject.com/) are needed.  



Once you have your R and Python installed, a simple command can be used to install Cellcano by typing:

```shell
pip install Cellcano
```

<u>Another installation option</u> which we recommend is to use [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands) to easily install our package or to manage your environment. When you get your conda installed, the command:

```shell
conda install -c marvinquiet cellcano-all
```

can be used to directly help you install all environments you need, including R, Python and Cellcano. Then you can follow the [ArchR installation](https://www.archrproject.com/) to install ArchR.



---

If you have successfully installed Cellcano, you could use the following command:

```shell
## Check the help documentation of Cellcano
Cellcano -h 
```

Then you will see the following console output:

```
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

*****

### 1. Download and process sample scATAC-seq raw data

#### Reference dataset
- First, we take one healthy donor labeled as `PBMC_Rep1` from [GSE129785](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129785) as our reference dataset. A brief introduction to this study is that it measures the single-cell chromatin landscape of human immune cell development.

  We start by downloading the raw data along with the processed metadata. The raw data is provided on GEO under accession number [GSE129785](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129785). The metadata needs some preprocessing procedures. The original metadata can be downloaded from GEO at [GSE129785 GSE129785_scATAC-Hematopoiesis-All.cell_barcodes.txt.gz](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129785/suppl/GSE129785_scATAC-Hematopoiesis-All.cell_barcodes.txt.gz). We preprocess the metadata by filtering the `Group `column to extract all cells from `PBMC_Rep1` donor. Then, we add the cell type labels based on the cluster annotation information provided by the [WashU EpiGenome Browser](http://epigenomegateway.wustl.edu/legacy/?genome=hg19&session=HcbHMSgBCc&statusId=28207718). The `SatPathy, Granja et al 2019 v1` track is used. For your convenience, we provide the processed metadata with curated major human PBMC cell types on [Dropbox](https://www.dropbox.com/s/d7e7wcb4nuy5l2q/GSM3722015_PBMC_Rep1_curated_metadata.csv). 

  Below are command lines which you could use to create the folder and download the data.

  ```shell
  # create and enter directory
  mkdir train_data
  cd train_data
  # Download sample fragment file
  wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3722nnn/GSM3722015/suppl/GSM3722015_PBMC_Rep1_fragments.tsv.gz
  # Download processed cell metadata
  wget https://www.dropbox.com/s/d7e7wcb4nuy5l2q/GSM3722015_PBMC_Rep1_curated_metadata.csv
  ```

- Next, we use ArchR to process the fragment file into gene score matrix. Since we embed it inside Cellcano, we can easily use the `Cellcano preprocess` to take in the raw data and turn it into the gene score matrix. Similarly, you can use `Cellcano preprocess -h` to check the usage.

  

  The `-i` option scans files ending with `fragments.tsv.gz` or `.bam` and load the files to ArchR. The`-g` option describes which genome will be used. `--threads` indicates how many threads are used to run ArchR. It would take around 30 mins for ArchR to process the data with threads as 4. 

  Here, our example command line will be:

  ```shell
  Cellcano preprocess -i train_data -o train_data -g hg19 --threads 4
  ```

  

  The output will become:

```shell
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
R[write to console]: 2022-05-15 10:37:39 : Constructing SummarizedExperiment, 1.337 mins elapsed.
R[write to console]: 2022-05-15 10:37:40 : Finished Matrix Creation, 1.359 mins elapsed.
```

The gene score matrix will be stored in the output `train_data` folder. 

1. ArchR_genescore.mtx.gz: stores the gene scores in a COO format.
2. ArchR_genescore_barcodes.tsv: stores the cell barcodes.
3. ArchR_genescore_genes.tsv: stores the gene symbols.

#### Target dataset

For the target dataset, we use a FACS-sorted human PBMC scATAC-seq dataset as example. The data can be downloaded from [GSE123578](GSE123578). The fragments files are having the suffixes as `fragments.tsv.gz`. Five cell types in human PBMC are provided by this study, including CD19 B cells, CD4 T cells, CD8 T cells, Monocytes and NK cells. In total, 21K cells are available. We then downsample 5K cells and provide the preprocessed gene score matrix on [Dropbox](https://www.dropbox.com/s/e7g9vem3oxt096l/FACS5K_genescore.csv).

```shell
## create and enter directory
mkdir test_data
cd test_data

## download processed target data
wget https://www.dropbox.com/s/e7g9vem3oxt096l/FACS5K_genescore.csv
```

### 2. Pre-train Cellcano's first round

With the preprocessed reference gene score matrix, we can use it as input for training a first-round prediction. When the input data ends with `.mtx`, the `-i` option takes the common file prefix indicating the matrix, barcodes and genes. `Cellcano train -h` can be used to check the usage.



The command for pre-training is:

```shell
Cellcano train -i train_data/ArchR_genescore -m train_data/GSM3722015_PBMC_Rep1_curated_metadata.csv -o output --prefix PBMC_Rep1_trained 
```

The output will be:

```shell
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

The pre-trained model will be stroed under the `output` directory with a folder named `PBMC_Rep1_trained`. Along with the trained model, the tSNE visualization and an [anndata](https://anndata.readthedocs.io/en/latest/) object storing the reference data are stored. The tSNE visualization color-labels cell type in the reference dataset which can be used to check the sepration of different cell types after summarizing to gene scores. The anndata object can be easily reloaded and combined with other anndata objects to form a larger reference dataset.

### 3. Predict on target dataset in Cellcano's second round

We use the pre-trained model to predict cell types on the target dataset. When the input data ends with `.csv`, we use the full path as input. 

The following command can be used for the prediction purpose. Since we are doing the second round which takes cells from the first round prediction as anchors to predict cell types for remaining cells, the `--tworound` option is needed. If this option is omitted, then the result will come from a direct MLP prediction. **We recommand using `--tworound`option because we discover the second-round prediction helps improve the prediction accuracy of non-anchor cells by ~3%. **By default, we use 40% cells from the target dataset as anchors to predict remaining cells. 

```shell
Cellcano predict -i test_data/FACS5K_genescore.csv --trained_model output/PBMC_Rep1_trainedMLP_model -o output --prefix predict_FACS --tworound
```

The output in the console will be:

```shell
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

The result will be stored in a csv file under the output directory indicated by `-o` option with  specified prefix. In our case, it will be `predict_FACScelltypes.csv`. 

### 4. Result analysis on our prediction

We use Sankey plot to analyze the cell type assignment. Here, we use R to run the following analysis.

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

In the Sankey plot, the left nodes are the true labels from FACS and the right nodes are predicted cell types. Most of the cells are correctly predicted. Although some cells in CD8 T cells are wrongly assigned as CD4 T cells and some cells in NK cells are wrongly assigned to CD8 T cells, the reason can be explained by the tSNE visualization on the FACS dataset below.  Cell label mixtures exist between NK cells and CD8 T cells as well as between CD4 T cells and CD8 T cells when the cells are colored by ground truth labels. This kind of cell mixture shows very high similarites among those cells even though they are assigned as different cell types. This is very hard for all classifiers to truly identify one from the other.  

![img](FACS_tSNE.png)

**Overall, Cellcano can achieve high prediction accuracy in supervised scATAC-seq celltyping and is fast and scalable to large number of cells.** If you have any questions using Cellcano, please let us know by emailing me or posting [Github issues](https://github.com/marvinquiet/Cellcano/issues).

