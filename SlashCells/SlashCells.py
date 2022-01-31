import os, argparse

parser = argparse.ArgumentParser(description="SlashCells: a supervised celltyping pipeline.",
        prog='SlashCells')
subparsers = parser.add_subparsers(help='sub-command help.')

## create parser for preprocessing
preprocess_parser = subparsers.add_parser('preprocess', 
        help='Run ArchR to preprocess raw input data (*fragment.tsv.gz, *.bam) to gene score.',
        dest="")
preprocess_parser.add_argument('-i', '--input', dest='raw_input', required=True,
        help="Raw scATAC-seq input file.")
preprocess_parser.add_argument('--input_type', required=True,
        help="Indicate input type: fragment or bam",
        choices=['fragment', 'bam'])
preprocess_parser.add_argument('-g', '--genome', required=True,
        help="Indicate input genome: mm9, mm10, hg19 or hg38", 
        choices=['mm9', 'mm10', 'hg19', 'hg38'])

## create parser for training
train_parser = subparsers.add_parser('train', help='Train a SlashCells model.')
#train_parser.add_argument()

## create parser for prediction
predict_parser = subparsers.add_parser('predict', help='Use SlashCells model to predict cell types.')
#predict_parser.add_argument()

args = parser.parse_args()

print(args)


