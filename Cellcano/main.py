import os, argparse
import logging
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(description="Cellcano: a supervised celltyping pipeline for single-cell omics.",
        prog='Cellcano')
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
            help="Raw scATAC-seq input file. We use ArchR to process the input into a gene score matrix. Cellcano will automatically look for *fragment.tsv.gz or *.bam files under the folder.")
    preprocess_parser.add_argument('-o', '--output_dir', dest='output_dir',
            type=str,
            help="Output directory to store results. Default: the same as input directory.")
    preprocess_parser.add_argument('--sample_names', dest='sample_names',
            type=str, nargs='+',
            help="The corresponding sample names for the input files. If not provided, Cellcano will use the file name as the sample names.")
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
    # --anndata
    # --model
    # -o, --output_dir
    # --prefix
    # --fs
    # --num_features
    # ============================
    train_parser = subparsers.add_parser('train', help='Train a Cellcano model.')
    train_parser.add_argument('-i', '--input', dest='input', 
            type=str,
            help="A COO matrix prefix or a csv input.")
    train_parser.add_argument('-m', '--metadata', dest='metadata', 
            type=str, 
            help="The annotation dataframe with row as barcodes or cell ID and columns as celltype. Notice: please make sure that your cell type indicator be 'celltype'.")
    train_parser.add_argument('--anndata', dest='anndata',
            type=str,
            help="Cellcano provides an option to load processed anndata object generated previously to reduce loading time.")
    train_parser.add_argument('--model', dest='model',
            type=str, default='MLP', choices=['MLP', 'KD', 'ADDA'],
            help="Model used to train the data.")
    train_parser.add_argument('-o', '--output_dir', dest='output_dir', 
            type=str, default="output",
            help="Output directory.")
    train_parser.add_argument('--prefix', dest='prefix', 
            type=str, default='train_',
            help="Output prefix.")
    train_parser.add_argument('--fs', dest='fs', 
            help="Feature selection methods.", 
            choices=["F-test", "noFS", "seurat"], type=str, default="F-test")
    train_parser.add_argument('--num_features', dest='num_features', 
            help="Feature selection methods.", type=int, default=1000)

    ## ===========================
    # Parser for predicting gene score matrix
    # -i, --input
    # --trained_model
    # --predict_type
    # -o, --output_dir
    # --prefix
    # ============================
    predict_parser = subparsers.add_parser('predict', help='Use Cellcano model to predict cell types.')
    predict_parser.add_argument('-i', '--input', dest='input', 
            type=str, required=True,
            help="A COO matrix prefix or a csv input.")
    predict_parser.add_argument('--trained_model', dest='trained_model',
            type=str, required=True,
            help="Path to the trained model.")
    predict_parser.add_argument('--tworound', dest='tworound',
            default=False, action="store_true",
            help="Whether predict with two rounds.")
    predict_parser.add_argument('-o', '--output_dir', dest='output_dir', 
            type=str, default="output",
            help="Output directory.")
    predict_parser.add_argument('--prefix', dest='prefix', 
            type=str, default='predict_',
            help="Output prefix.")
    
    args = parser.parse_args()
    logger.debug("User input arguments: ", args)
    return args


def main():
    args = parse_args()
    ## preprocessing
    if "preprocess" == args.cmd_choice:
        from Cellcano import preprocess
        preprocess.preprocess(args)
    ## training process 
    if "train" == args.cmd_choice:
        from Cellcano import train
        if not os.path.exists(args.output_dir):
            logger.info("Creating output directory: %s" % args.output_dir)
            os.makedirs(args.output_dir, exist_ok=True)
        if args.model == "KD":
            train.train_KD(args)
        if args.model == "ADDA":
            ## For future API
            pass
        if args.model == "MLP":
            train.train_MLP(args)
    ## prediction process
    if "predict" == args.cmd_choice:
        from Cellcano import predict
        if not os.path.exists(args.trained_model):
            logger.error("The input model does not exist.")
        predict.predict(args)
