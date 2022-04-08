'''
A celltyping pipeline integrating DFN, GEDFN, MLP, SVM and Random Forest

For scmap and CHETA, there are other Rscripts and then use bash to run
'''
import os, argparse

from ATACseq import method_utils, dataloading_utils
from ATACseq.preprocess.process_train_test_data import *

if __name__ == "__main__":
    data_dir = "/home/wma36/gpu/data"

    ## parse arguments
    import argparse
    parser = argparse.ArgumentParser(description="Celltyping pipeline.")
    parser.add_argument('data_source', help="Load which dataset",
        choices=[
            ## ==== immune cells Greenleaf 
            'immunePBMC_genescore', 'immunePBMC_genescore_major', 'immunePBMC_genescore_major_cv',
            'immunePBMC_ciceroGA', 
            'immunePBMC_SnapATAC_harmony', 'immunePBMC_SnapATAC_noharmony',
            'immunePBMC_ATACbin', 'immunePBMC_ATACbin_major',
            ## ==== AML study Greenleaf
            'AMLPBMC_genescore', 'AMLPBMC_genescore_major', 'AMLPBMC_genescore_major_cv',
            'AMLPBMC_ATACbin_major', 'AMLPBMC_ATACbin',
            ## ==== PBMC cross dataset, same modality
            'immunePBMC_to_AMLPBMC_ATACbin', 'AMLPBMC_to_immunePBMC_ATACbin',
            '10XPBMC_to_immunePBMC_ATACbin', '10XPBMC_to_AMLPBMC_ATACbin',
            'immunePBMC_to_10XPBMC_ATACbin', 'AMLPBMC_to_10XPBMC_ATACbin',
            'immunePBMC_to_AMLPBMC_genescore', 'AMLPBMC_to_immunePBMC_genescore',
            '10XPBMC_to_immunePBMC_genescore', '10XPBMC_to_AMLPBMC_genescore',
            'immunePBMC_to_10XPBMC_genescore', 'AMLPBMC_to_10XPBMC_genescore',
            ## === PBMC cross dataset, cross modality
            '10XPBMC_exprs_to_10XPBMC_genescore', '10XPBMC_exprs_to_immunePBMC_genescore',
            '10XPBMC_exprs_to_AMLPBMC_genescore',
            ## === quantile normalization
            '10XPBMC_exprs_to_10XPBMC_genescore_qnorm', '10XPBMC_exprs_to_immunePBMC_genescore_qnorm',
            '10XPBMC_exprs_to_AMLPBMC_genescore_qnorm',
            ## === combined individuals
            'combined_GreenleafPBMC_to_10X'
            ])

    parser.add_argument('-m', '--method', help="Run which method",
        choices=['MLP'], 
        required=True)
    parser.add_argument('--select_on', help="Feature selection on train or test, or None of them",
        choices=['train', 'test'])
    parser.add_argument('--select_method', help="Feature selection method, Seurat/FEAST or None",
            choices=['Seurat', 'FEAST', 'F-test'])
    parser.add_argument('--n_features', help="Number of features selected",
            default=1000, type=int)
    parser.add_argument('--train', help="Specify which as train")
    parser.add_argument('--test', help="Specify which as test")
    parser.add_argument('--sample_seed', help="Downsample seed in combined individual effect", 
            default=0, type=int)

    args = parser.parse_args()
    pipeline_dir = "/home/wma36/gpu/celltyping_methodTrials/results/result_PBMCs"

    result_prefix = pipeline_dir+os.sep+"result_"+args.data_source+'_'+\
        args.train+'_to_'+args.test
    os.makedirs(result_prefix, exist_ok=True)
    ## create file directory 
    if args.select_on is None and args.select_method is None:
        result_dir = result_prefix+os.sep+"no_feature"
    else:
        result_dir = result_prefix+os.sep+args.select_method+'_'+\
                str(args.n_features)+'_on_'+args.select_on
    os.makedirs(result_dir, exist_ok=True)

    load_ind, train_adata, test_adata = load_adata(result_dir)
    if not load_ind:
        if train_adata is None or test_adata is None:
            train_adata, test_adata = dataloading_utils.load_Greenleaf_adata(
                data_dir, args=args)
        if train_adata is None or test_adata is None:
            train_adata, test_adata = dataloading_utils.load_crossdataset_ATACbin_adata(
                data_dir, args=args)
        if train_adata is None or test_adata is None:
            train_adata, test_adata = dataloading_utils.load_crossdataset_genescore_adata(
                data_dir, args=args)

        if 'qnorm' in args.data_source:
            train_adata, test_adata = dataloading_utils.process_loaded_data(
                train_adata, test_adata, result_dir, args=args,
                preprocess=False, quantile_norm=True)
 
        elif 'genescore' in args.data_source or 'ciceroGA' in args.data_source:
            train_adata, test_adata = dataloading_utils.process_loaded_data(
                train_adata, test_adata, result_dir, args=args,
                preprocess=True)
                #preprocess=False)

        elif args.data_source in ["immunePBMC_SnapATAC_harmony", "immunePBMC_SnapATAC_noharmony"]:
            train_adata, test_adata = dataloading_utils.process_loaded_data(
                train_adata, test_adata, result_dir, args=args,
                preprocess=False, log=False)

        elif 'ATACbin' in args.data_source:
            train_adata, test_adata = dataloading_utils.process_loaded_data(
                train_adata, test_adata, result_dir, args=args, save_data=False,
                ATACseq=True)

        else:
            train_adata, test_adata = dataloading_utils.process_loaded_data(
                train_adata, test_adata, result_dir, args=args)
 
        print("Train anndata: \n", train_adata)
        print("Test anndata: \n", test_adata)
    method_utils.run_pipeline(args, train_adata, test_adata, data_dir, result_dir)

