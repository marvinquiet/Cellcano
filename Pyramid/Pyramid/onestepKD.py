'''
Method related utility pipelines
'''
import os, sys 
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scipy

## import my package
from ..utils import _utils, _logger

def load_train(inputfile):
    '''Load training data
    ---
    Input:
        - inputfile: input from argument parser
    ---
    Output:
        - adata: gene score anndata object
    '''
    if '.csv' in inputfile:
        adata = _csv_data_loader(inputfile)
    else:
        adata = _COOmtx_data_loader(inputfile)
    return adata


def train_onestepKD(args,
        teacher_MLP_DIMS=[128, 32, 8], student_MLP_DIMS=[64, 16], 
        pred_celltype_cols="pred_celltypes", prefix="Pyramid_onestep"):
    '''Train one step KD model
        1. Load train data and metadata
        2. Feature selection
        3. Log-Norm and scale data
    ---
    Input:
        - args: user's input arguments
    '''

    '''Two step prediction
    - Select low entropy cells from each cell type
    - Use low entropy cells as reference to predict the rest target cells
    '''

    ## load input data
    train_adata = load_train(args.input)
    metadata = _utils._metadata_loader(args.metadata)

    common_cells = set(train_adata.obs_names).intersection(set(metadata.index))
    print("%d common cells found between input data and metadata." % len(common_cells))
    if len(common_cells) < 100:
        print("There are too few cells. Pyramid might not be accurate.")
    train_adata = train_adata[list(common_cells)]
    train_adata.obs = train_adata.obs.merge(metadata, 
            left_on="barcode", right_index=True, how='left')

    ## Feature selection
    if args.fs == "noFS":
        print("Pyramid will not perform feature selection.\n")
        num_features = train_adata.shape[0]
    else:
        num_features = args.num_features
        if num_features < train_adata.shape[0]:
            print("Number of features is larger than data. Pyramid will not perform feature selection.\n")
            num_features = train_adata.shape[0]

    ## filter cells/genes, etc
    train_adata = _utils._process_adata(train_adata)
    print("Data shape after processing: %d*%d" % (train_adata.shape[0], train_adata.shape[1]))

    if args.fs == "F-test":
        print("Use F-test to select features.\n")
        if scipy.sparse.issparse(train_adata.X) or \
                isinstance(train_adata.X, pd.DataFrame):
            tmp_data = train_adata.X.toarray()
        else:
            tmp_data = train_adata.X
        ## write out original read count matrix
        tmp_df = pd.DataFrame(data=np.round(tmp_data, 3), 
                index=train_adata.obs_names, columns=train_adata.var_names).T
        tmp_df_path = args.output_dir+os.sep+"tmp_counts.csv"
        tmp_df.to_csv(tmp_df_path)
        ## write out cell annotations based on train
        cell_annots = tmp_adata.obs[celltype_label].tolist()
        cell_annots_path = args.output_dir+os.sep+"tmp_cell_annots.txt"
        with open(cell_annots_path, 'w') as f:
            for cell_annot in cell_annots:
                f.write("%s\n" % cell_annot)
        os.system("Rscript --vanilla " + FEAST_FTEST_RSCRIPT_PATH + " "+ tmp_df_path + " " + 
                cell_annots_path + " " + str(num_features))
        os.system("rm {}".format(cell_annots_path))  ## remove the temporaty cell annotations
        os.system("rm {}".format(tmp_df_path))  ## remove the temporaty counts

        ftest_file = args.output_dir+os.sep+'F-test_features.txt'
        with open(ftest_file) as f:
            features = f.read().splitlines()
        features.sort()
        train_adata = train_adata[:, features]

    if args.fs == "seurat":
        print("Use seurat in scanpy to select features.\n")
        sc.pp.highly_variable_genes(train_adata, n_top_genes=num_features, subset=True)


 






    ## prepare for data
    enc, x_train, y_train, x_test = \
            _utils._prepare_data(train_adata, test_adata, enc=None,
                               celltype_cols=celltype_cols)

    ## initialize teacher model
    teacher = _utils._init_MLP(x_train, y_train, dims=teacher_MLP_DIMS, seed=RANDOM_SEED)

    ## teacher model
    teacher = _utils._init_MLP(x_train, y_train, dims=teacher_MLP_DIMS, seed=RANDOM_SEED)
    teacher.compile()
    start = time.time()
    teacher.fit(x_train, y_train, batch_size=BATCH_SIZE)
    ## student model -> actually same model, just used the concept of distillation
    student = initialize_MLP(x_train, y_train, dims=student_MLP_DIMS, seed=RANDOM_SEED)
    # Initialize and compile distiller
    distiller = run_distiller(x_train, y_train, 
            student_model=student.model,
            teacher_model=teacher.model)
    y_pred = tf.nn.softmax(distiller.student.predict(x_test)).numpy()
    pred_celltypes = prob_to_label(enc, y_pred)
    test_adata.obs[pred_celltype_cols] = pred_celltypes
    ## analyze first step MLP
    analyze_MLP_prediction(train_adata, test_adata, result_dir=result_dir,  
                celltype_cols=celltype_cols, pred_celltype_cols=pred_celltype_cols,
                prefix=prefix)
    return test_adata

