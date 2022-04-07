'''
Method utils
'''
import os, sys 
import numpy as np
import pandas as pd

## import my packages
from utils import _utils

def run_MLP_pipeline(args, train_adata, test_adata, result_dir, enc=None,
        celltype_cols="cell.type", pred_celltype_cols="pred_celltypes",
        prefix=None):
    '''Directly run MLP to predict
    '''
    if prefix is None:
        prefix = args.method

    ## generate x_train, y_train, x_test
    enc, x_train, y_train, x_test = \
            prepare_data(train_adata, test_adata, enc=None, 
                    celltype_cols=celltype_cols)
    ### run MLP
    print("\n\n=== MLP\n")
    y_pred, run_time = run_MLP(x_train, y_train, x_test, 
            dims=MLP_DIMS, batch_size=BATCH_SIZE, seed=RANDOM_SEED,
            save_dir=result_dir+os.sep+prefix+'_checkpoint',
            sample_weight=sample_weight, class_weight=class_weight) ## save weights
    pred_celltypes = prob_to_label(enc, y_pred)
    test_adata.obs[pred_celltype_cols] = pred_celltypes
 
    print("\n\n=== Run time:", run_time)
    return y_pred, enc, run_time


def select_low_entropy_cells(args, test_adata, low_entropy_quantile,
        pred_celltype_cols):
    '''Select low entropy cells from each cell type
    '''
    low_entropy_cells = []
    for celltype in set(test_adata.obs[pred_celltype_cols]):
        celltype_df = test_adata.obs[test_adata.obs[pred_celltype_cols] == celltype]
        if args.random_anchor:
            seed(RANDOM_SEED)
            cells = sample(celltype_df.index.tolist(), math.ceil(low_entropy_quantile*celltype_df.shape[0]))
        else:
            entropy_cutoff = np.quantile(celltype_df['entropy'], q=low_entropy_quantile)
            ## change to < instead of <= to deal with ties
            cells = celltype_df.index[np.where(celltype_df['entropy'] < entropy_cutoff)[0]].tolist()
        low_entropy_cells.extend(cells)
    high_entropy_cells = list(set(test_adata.obs_names) - set(low_entropy_cells))
    test_adata.obs.loc[low_entropy_cells, 'entropy_status'] = "low"
    test_adata.obs.loc[high_entropy_cells, 'entropy_status'] = "high"
    return test_adata
 

def run_MLP_pipeline(args, train_adata, test_adata, result_dir, enc=None,
        celltype_cols="cell.type", pred_celltype_cols="pred_celltypes",
        prefix=None, sample_weight=None, class_weight=None):
    '''Run methods
    '''
    if prefix is None:
        prefix = args.method
    ## generate x_train, y_train, x_test
    enc, x_train, y_train, x_test = \
            prepare_data(train_adata, test_adata, enc=None, 
                    celltype_cols=celltype_cols)
    ### run MLP
    print("\n\n=== MLP\n")
    y_pred, run_time = run_MLP(x_train, y_train, x_test, 
            dims=MLP_DIMS, batch_size=BATCH_SIZE, seed=RANDOM_SEED,
            save_dir=result_dir+os.sep+prefix+'_checkpoint',
            sample_weight=sample_weight, class_weight=class_weight) ## save weights
    pred_celltypes = prob_to_label(enc, y_pred)
    test_adata.obs[pred_celltype_cols] = pred_celltypes
    analyze_MLP_prediction(train_adata, test_adata, result_dir=result_dir, run_time=run_time, 
                celltype_cols=celltype_cols, pred_celltype_cols=pred_celltype_cols,
                prefix=prefix)
 
    print("\n\n=== Run time:", run_time)
    return y_pred, enc, run_time

def run_Pyramid_pipeline(args, train_adata, test_adata, result_dir,
        teacher_MLP_DIMS=[128, 32, 8],
        celltype_cols="cell.type", pred_celltype_cols="pred_celltypes", 
        low_entropy_quantile=0.2, prefix="MLP_twostep_KD"):
    '''Two step prediction
    - Select low entropy cells from each cell type
    - Use low entropy cells as reference to predict the rest target cells
    '''
    if prefix is None:
        prefix = args.method

    ## save results to different quantile cutoff
    quantile_dir = result_dir+os.sep+str(low_entropy_quantile)
    os.makedirs(quantile_dir, exist_ok=True)

    enc, x_train, y_train, x_test = \
            prepare_data(train_adata, test_adata, enc=None,
                               celltype_cols=celltype_cols)
    ## teacher model
    teacher = initialize_MLP(x_train, y_train, dims=teacher_MLP_DIMS, seed=RANDOM_SEED)
    teacher.compile()
    start = time.time()
    teacher.fit(x_train, y_train, batch_size=BATCH_SIZE)
    ## student model -> actually same model, just used the concept of distillation
    student = initialize_MLP(x_train, y_train, dims=MLP_DIMS, seed=RANDOM_SEED)
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
                prefix='MLP_KD')
    y_test = enc.transform(test_adata.obs[[celltype_cols]]).toarray()

    ## === Step1: explore entropy
    test_adata = analyze_entropy_per_assigned(y_pred, y_test, test_adata, quantile_dir)
    test_adata = select_low_entropy_cells(args, test_adata,  
            low_entropy_quantile, pred_celltype_cols)
    analyze_entropy_per_celltype(test_adata, quantile_dir, celltype_cols)

    ## === Step2: take low entropy cells as reference
    test_adata.obs.rename(columns={pred_celltype_cols:'first_pred_celltypes'}, inplace=True)
    test_adata.obs[pred_celltype_cols] = test_adata.obs['first_pred_celltypes']

    low_entropy_cells = test_adata.obs_names[np.where(test_adata.obs['entropy_status'] == 'low')].tolist()
    high_entropy_cells = test_adata.obs_names[np.where(test_adata.obs['entropy_status'] == 'high')].tolist()
    test_ref_adata = test_adata[low_entropy_cells]
    test_tgt_adata = test_adata[high_entropy_cells]

    ## prepare data for next step KD
    enc, x_tgt_train, y_tgt_train, x_tgt_test = \
               prepare_data(test_ref_adata, test_tgt_adata, enc=enc, 
                       celltype_cols="first_pred_celltypes")
    teacher.fit(x_tgt_train, y_tgt_train, batch_size=BATCH_SIZE)
    distiller.fit(x_tgt_train, y_tgt_train, epochs=30, 
                       validation_split=0.0, verbose=2)
    end = time.time()
    traintime = end-start  ## record training time

    start = time.time()
    y_pred_tgt = tf.nn.softmax(distiller.student.predict(x_tgt_test)).numpy()
    end = time.time()
    predtime = end-start ## record prediction time

    pred_celltypes = prob_to_label(enc, y_pred_tgt)
    test_adata.obs.loc[high_entropy_cells, pred_celltype_cols] = pred_celltypes
    analyze_twostep_prediction(test_adata, result_dir=quantile_dir, prefix=prefix,
            traintime=traintime, predtime=predtime)
