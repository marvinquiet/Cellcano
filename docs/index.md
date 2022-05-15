## Cellcano Tutorial

Authors: Wenjing Ma (wenjing.ma@emory.edu), Dr. Hao Wu (hao.wu@emory.edu), Emory University
Latest revision: 04-29-2022

In this tutorial, we will guide you through Cellcano, a supervised cell type annotation (celltyping) tool in scATAC-seq. We use human Peripheral Blood Mononuclear Cells (PBMC) datasets as examples. 

We are aiming at celltyping in single cell genomics, including scRNA-seq, scATAC-seq, scBS-seq, etc. If you are interested, stay tuned!

### 1. Download sample scATAC-seq raw data

#### Reference dataset
First, we take a healthy donor `PBMC_Rep1` from [GSE129785](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129785) to be our reference dataset. This study measures the single-cell chromatin landscape of human immune cell development.

We start by downloading the raw data along with the processed metadata. The raw data is provided on GEO under accession number GSE129785. The metadata needs some preprocessing procedures. The original metadata can be first downloaded from GEO at [GSE129785 GSE129785_scATAC-Hematopoiesis-All.cell_barcodes.txt.gz](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129785/suppl/GSE129785_scATAC-Hematopoiesis-All.cell_barcodes.txt.gz). We extract out individuals related to PBMC using the `Group` column and extract out those cells lying in `PBMC_Rep1, PBMC_Rep2, PBMC_Rep3, PBMC_Rep4` individuals. Then, we annotate the cluster information according to the [WashU EpiGenome Browser](http://epigenomegateway.wustl.edu/legacy/?genome=hg19&session=HcbHMSgBCc&statusId=28207718) under the `SatPathy, Granja et al 2019 v1` track. For your convenience, we provide the processed metadata with cell types in our Dropbox.

```shell
# create and enter directory
mkdir train_data
cd train_data
# Download sample fragment file
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3722nnn/GSM3722015/suppl/GSM3722015_PBMC_Rep1_fragments.tsv.gz
# Download processed cell metadata
wget https://www.dropbox.com/s/u4nba4rdqvljncg/GSE129785_scATAC-Hematopoiesis-PBMC.cell_barcodes.txt
```

Next, we use Cellcano to call ArchR to process the fragment file to gene score matrix.

```shell


```

#### Target dataset
Here, as an example, we use a FACS-sorted human PBMC scATAC-seq dataset to be our target dataset.





You can use the [editor on GitHub](https://github.com/marvinquiet/Pyramid/edit/master/docs/index.md) to maintain and preview the content for your website in Markdown files.

Whenever you commit to this repository, GitHub Pages will run [Jekyll](https://jekyllrb.com/) to rebuild the pages in your site, from the content in your Markdown files.

### Markdown

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [Basic writing and formatting syntax](https://docs.github.com/en/github/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/marvinquiet/Pyramid/settings/pages). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://support.github.com/contact) and we’ll help you sort it out.
