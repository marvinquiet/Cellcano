import os, argparse
from Pyramid import *

parser = argparse.ArgumentParser(description="Pyramid: a supervised celltyping pipeline for single-cell omics.",
        prog='Pyramid')
subparsers = parser.add_subparsers(help='sub-command help.', dest='cmd_choice')

# ===========================
# Parser for preprocess data into gene score matrix
# -i, --input_dir
# -o, --output_dir
# -g, --genome
# --sample_names
# --tilemat
# --threads
# ============================

preprocess_parser = subparsers.add_parser('preprocess', 
        help='Run ArchR to preprocess raw input data (*fragments.tsv.gz, *.bam) to gene score.')
preprocess_parser.add_argument('-i', '--input_dir', dest='input_dir',
        required=True, type=str,
        help="Raw scATAC-seq input file. We use ArchR to process the input into a gene score matrix. Pyramid will automatically look for *fragment.tsv.gz or *.bam files under the folder.")
preprocess_parser.add_argument('-o', '--output_dir', dest='output_dir',
        type=str,
        help="Output directory to store results. Default: the same as input directory.")
preprocess_parser.add_argument('--sample_names', dest='sample_names',
        type=str, nargs='+',
        help="The corresponding sample names for the input files. If not provided, Pyramid will use the file name as the sample names.")
preprocess_parser.add_argument('-g', '--genome', required=True,
        type=str,
        help="Indicate input genome: mm9, mm10, hg19 or hg38", 
        choices=['mm9', 'mm10', 'hg19', 'hg38'])
preprocess_parser.add_argument('--threads',
        type=int,
        help="Threads used to run ArchR.", default=1)

# ===========================
# Parser for training gene score matrix to a model
# -i, --input
# -m, --metadata
# -o, --output_dir
# --fs
# --num_features
# --teacher_ns
# --student_ns
# ============================

train_parser = subparsers.add_parser('train', help='Train a Pyramid model.')
train_parser.add_argument('-i', '--input', dest='input', 
        type=str, required=True,
        help="A COO matrix prefix or a csv input.")
train_parser.add_argument('-m', '--metadata', dest='metadata', 
        type=str, required=True,
        help="The annotation dataframe with row as barcodes or cell ID and columns as celltype. Notice: please make sure that your cell type indicator be 'celltype'.")
train_parser.add_argument('-o', '--output_dir', dest='output_dir', type=str,
        help="Output directory.")
train_parser.add_argument('--fs', dest='fs', 
        help="Feature selection methods.", 
        choices=["F-test", "noFS", "seurat"], type=str, default="F-test")
train_parser.add_argument('--num_features', dest='num_features', 
        help="Feature selection methods.", type=int, default=1000)
train_parser.add_argument('--teacher_ns', dest='teacher_ns', 
        type=int, nargs='+',
        help="Teacher network structure.")
train_parser.add_argument('--student_ns', dest='student_ns', 
        type=int, nargs='+',
        help="Student network structure.")

## create parser for prediction
predict_parser = subparsers.add_parser('predict', help='Use Pyramid model to predict cell types.')
predict_parser.add_argument('--gs_prefix', dest='gs_prefix', required=True,
        help="Gene Score matrix generated for test data.")
predict_parser.add_argument('-o', '--output', dest='output', default="Pyramid_predict",
        help="A dataframe with row as barcodes and columns as predicted cell types.")

args = parser.parse_args()
print(args)

if "preprocess" == args.cmd_choice:
    preprocess.preprocess(args)

if "train" == args.cmd_choice:

    pass

if "predict" == args.cmd_choice:
    pass


