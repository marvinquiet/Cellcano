### Processing different scATAC-seq inputs

Utilize existing ArchR pipeline to process data to 5kb-bin count matrix
- Input: \*.bam, \*fragment.tsv.gz
- ArchR Output:
    - \*.mtx.gz: 5kb tile matrix
    - \*.barcodes.gz: barcode information
    - \*.5kbtile.gz: 5kb tile information
- Secondary Output:
    - \*.genescore.mtx.gz: ArchR best gene score matrix
    - \*.genes.gz: gene information
