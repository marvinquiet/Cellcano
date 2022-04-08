'''
A celltyping pipeline integrating DFN, GEDFN, MLP, SVM and Random Forest

For scmap and CHETA, there are other Rscripts and then use bash to run
'''

import os, sys, math, time
import anndata
import numpy as np
import pandas as pd
import scipy
import tensorflow as tf

from sklearn.preprocessing import OneHotEncoder
from sklearn import metrics

from random import seed
from numpy.random import seed
from tensorflow import set_random_seed

# GPU settings and reproducible
#os.environ["CUDA_VISIBLE_DEVICES"]="0"

#Set seeds
RANDOM_SEED=0
seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)
set_random_seed(RANDOM_SEED)

## every time reset all graphs and run something new
tf.compat.v1.reset_default_graph()

## import my package
from ATACseq import dataloading_utils

def run_MLP(x_train, y_train, x_test, dims=[128, 64, 32, 16], 
        batch_size=128, seed=0, save_dir=None):
    '''Pipeline for MLP
    '''
    from ATACseq.test_MLP import MLP

    mlp = MLP(dims)
    mlp.input_shape = (x_train.shape[1], )
    mlp.n_classes = len(set(y_train.argmax(1)))
    mlp.random_state = seed
 
    if save_dir is None:
        mlp.fit(x_train, y_train, batch_size=batch_size)
    else:
        if os.path.exists(save_dir):
            mlp.init_MLP()
            mlp.model.load_weights(save_dir)
        else:
            mlp.fit(x_train, y_train, batch_size=batch_size)
            mlp.model.save_weights(save_dir)
    mlp.model.summary()

    ## only calculate predicting time
    start = time.time()
    y_pred = mlp.predict(x_test)
    end = time.time()
    return y_pred, end-start

def run_SVM(x_train, y_train, x_test, kernel="rbf", seed=0):
    '''Fit SVM model with a RBF kernel
    SVM decision_function_shape: can be one-vs-one or one-vs-rest, in our case,
    according to our input, we should use one-vs-rest
    '''
    if "rbf" == kernel:
        from sklearn.svm import SVC
        model = SVC(decision_function_shape="ovr", kernel=kernel, random_state=seed)
        #model = SVC(decision_function_shape="ovo", kernel=kernel, random_state=seed)
    elif "linear" == kernel:
        from sklearn.svm import LinearSVC
        model = LinearSVC(multi_class='ovr', max_iter=1e4, random_state=seed)

    ## fit model
    model.fit(x_train, y_train)
    start = time.time()
    y_pred = model.predict(x_test)
    end = time.time()
    return y_pred, end-start

def analyze_prediction(y_pred, y_test, enc_obj, test_adata, n_clusters, 
        run_time=None, celltype_cols="cell.type", result_dir="./result", prefix="model"):
    ''' Analyze the celltyping prediction result
    @ y_pred: the result from predicting
    @ enc_obj: onehot encoder object to convert back to cell types
    @ n_clusters: number of clusters in train datasets
    @ thres: if thres, then we generate "unassigned" results
    '''
    ## evaluation metrics
    acc = metrics.accuracy_score(y_test, y_pred)
    ARI = metrics.cluster.adjusted_rand_score(y_test, y_pred)
    macroF1 = metrics.f1_score(y_test, y_pred, average="macro")
    print("*****=====", "Testing accuracy: ", acc, " Testing ARI: ", ARI, 
            "Testing MacroF1: ", macroF1, "=====*****")
    
    ## write metrics result to file
    with open(result_dir+os.sep+prefix+"_metrics.txt", 'w') as f:
        f.write("Acc:%s\n" % str(acc))
        f.write("ARI:%s\n" % str(ARI))
        f.write("macroF1:%s\n" % str(macroF1))
        if run_time is not None:
            f.write("runtime:%s\n" % str(run_time))

    ## analyze predicted cell types 
    pred_labels = y_pred
    pred_onehot = np.zeros((pred_labels.size, n_clusters))
    pred_onehot[np.arange(pred_labels.size), pred_labels] = 1
    pred_celltypes = enc_obj.inverse_transform(pred_onehot)
    print("=== Predicted celltypes: ", set(pred_celltypes.flatten()))

    test_adata.obs['pred_celltypes'] = pred_celltypes
    test_adata.obs.to_csv(result_dir+os.sep+prefix+"_predicted_obs.csv")
 
    ## visualization
    import matplotlib.pyplot as plt
    import scanpy.api as sc
    if len(test_adata.obsm) == 0:
        sc.pp.highly_variable_genes(test_adata, n_top_genes=1000, subset=True)
        sc.tl.pca(test_adata, random_state=RANDOM_SEED) ## PCA
        sc.tl.tsne(test_adata, n_pcs=50, use_rep="X_pca",  ## tSNE
                learning_rate=300,
                perplexity=30,
                n_jobs=4,
                random_state=RANDOM_SEED)
    sc.pl.tsne(test_adata, color=[celltype_cols, "pred_celltypes"], size=15)
    plt.savefig(result_dir+os.sep+prefix+"_prediction_result.png")
    print("=== Finish visualizing..")


def run_pipeline(args, train_adata, test_adata, data_dir, result_dir,
        batch_size = 128, celltype_cols="cell.type"):
    '''Run methods
    '''
    ## Hyperparameters for network
    MLP_dims = [128, 32, 8] if train_adata.shape[0] >=5000 else [64, 16]  ## set dimension

    ## OneHotEncoding the celltypes
    enc, x_train, y_train, x_test, y_test = \
            prepare_data(train_adata, test_adata, enc=None, 
                    train_celltype_cols=celltype_cols, test_celltype_cols=celltype_cols)

    n_clusters = len(set(train_adata.obs[celltype_cols]))
    if "MLP" == args.method:
        ### --- run MLP
        print("\n\n=== MLP\n")
        y_pred, exec_time = run_MLP(x_train, y_train, x_test, 
                dims=MLP_dims, batch_size=batch_size, seed=RANDOM_SEED,
                save_dir=result_dir+os.sep+'checkpoint') ## save weights
        print("\n\n=== Run time:", exec_time)
        analyze_prediction(y_pred.argmax(1), y_test.argmax(1), enc, test_adata, 
                n_clusters, exec_time, celltype_cols=celltype_cols, 
                result_dir=result_dir, prefix=args.method)

    if "SVM_RBF" == args.method:
        ## --- run SVM
        print("\n\n=== SVM RBF kernel\n")
        y_pred, exec_time = run_SVM(x_train, y_train.argmax(1), x_test)
        print("\n\n=== Run time:", exec_time)
        analyze_prediction(y_pred, y_test.argmax(1), enc, test_adata, 
                n_clusters, exec_time, celltype_cols=celltype_cols, 
                result_dir=result_dir, prefix=args.method)

    if "SVM_linear" == args.method:
        ## --- run SVM
        print("\n\n=== SVM linear kernel\n")
        y_pred, exec_time = run_SVM(x_train, y_train.argmax(1), x_test, kernel="linear")
        print("\n\n=== Run time:", exec_time)
        analyze_prediction(y_pred, y_test.argmax(1), enc, test_adata, 
                n_clusters, exec_time, celltype_cols=celltype_cols, 
                result_dir=result_dir, prefix=args.method)

    return y_pred, enc

def run_twostep_MLP_pipeline(args, train_adata, test_adata, data_dir, result_dir,
        batch_size=32, celltype_cols="cell.type", pred_celltype_cols="pred_celltypes", 
        entropy_selection="celltype", low_entropy_quantile=0.2, FS=False):
    '''Two step prediction
    - Select low entropy cells from each cell type
    - Use low entropy cells as reference to predict the rest target cells

    @entropy_selection: celltype - select in each cell type
                        all - select in all cells
    @FS: whether perform feature selection
    '''
    ## first step prediction
    train_fs_adata, test_fs_adata = dataloading_utils.process_loaded_data(
            train_adata.copy(), test_adata.copy(), result_dir, args=args,
            preprocess=True, save_data=False)
    y_pred, enc = run_pipeline(args, train_fs_adata, test_fs_adata, data_dir, result_dir)

    test_fs_cells = [cell.replace('-test', '') for cell in test_fs_adata.obs_names]
    test_adata = test_adata[test_fs_cells]
    test_adata.obs[pred_celltype_cols] = test_fs_adata.obs[pred_celltype_cols].tolist()
    test_adata.obsm = test_fs_adata.obsm
    y_test = enc.transform(test_adata.obs[[celltype_cols]]).toarray()

    ## save results to different quantile cutoff
    quantile_dir = result_dir+os.sep+'quantile_'+str(low_entropy_quantile)
    os.makedirs(quantile_dir, exist_ok=True)

    ## === Step1: explore entropy
    y_pred_label = y_pred.argmax(1)
    y_test_label = y_test.argmax(1)
    correct_idx = np.where(y_pred_label == y_test_label)[0]
    wrong_idx = np.where(y_pred_label != y_test_label)[0]
    y_pred_correct_probs = y_pred[correct_idx]
    y_pred_wrong_probs = y_pred[wrong_idx]

    y_pred_correct_entropy = [-np.nansum(y_pred_correct_probs[i]*np.log(y_pred_correct_probs[i])) for i in range(y_pred_correct_probs.shape[0])]
    y_pred_wrong_entropy = [-np.nansum(y_pred_wrong_probs[i]*np.log(y_pred_wrong_probs[i])) for i in range(y_pred_wrong_probs.shape[0])]

    test_adata.obs['assigned_status'] = "wrong"
    correct_cells = test_adata.obs_names[correct_idx]
    wrong_cells = test_adata.obs_names[wrong_idx]
    test_adata.obs.loc[correct_cells, 'assigned_status'] = "correct"

    test_adata.obs['entropy'] = None
    test_adata.obs.loc[correct_cells, 'entropy'] = y_pred_correct_entropy
    test_adata.obs.loc[wrong_cells, 'entropy'] = y_pred_wrong_entropy
    test_adata.obs['entropy'] = pd.to_numeric(test_adata.obs['entropy'])

    import matplotlib.pyplot as plt
    import scanpy.api as sc
    import seaborn as sns
    plt.rcParams.update({'font.size': 18})
    from textwrap import wrap

    fig = plt.figure(figsize=(10, 6), dpi=300)
    sns.boxplot(x='assigned_status', y='entropy', data=test_adata.obs)
    plt.savefig(quantile_dir+os.sep+'Entropy_per_assigned.png')
 
    if entropy_selection == 'celltype': 
        low_entropy_cells = []
        for celltype in set(test_adata.obs[pred_celltype_cols]):
            celltype_df = test_adata.obs[test_adata.obs[pred_celltype_cols] == celltype]
            entropy_cutoff = np.quantile(celltype_df['entropy'], q=low_entropy_quantile)
            cells = celltype_df.index[np.where(celltype_df['entropy'] <= entropy_cutoff)[0]].tolist()
            low_entropy_cells.extend(cells)
    if entropy_selection == 'all':
        entropy_cutoff = np.quantile(test_adata.obs['entropy'], q=low_entropy_quantile)
        low_entropy_cells = test_adata.obs.index[np.where(test_adata.obs['entropy'] <= entropy_cutoff)[0]].tolist()
    high_entropy_cells = list(set(test_adata.obs_names) - set(low_entropy_cells))

    test_adata.obs['entropy_status'] = "high"
    test_adata.obs.loc[low_entropy_cells, 'entropy_status'] = "low"

    low_entropy_df = test_adata[low_entropy_cells].obs
    high_entropy_df = test_adata[high_entropy_cells].obs

    fig = plt.figure(figsize=(10, 6), dpi=300)
    ax = sns.boxplot(x=celltype_cols, y='entropy', data=low_entropy_df)
    ax.set_xticklabels(ax.get_xticklabels(),rotation = 30)
    plt.tight_layout()
    plt.savefig(quantile_dir+os.sep+'Entropy_per_celltype_lowentropy.png')

    fig = plt.figure(figsize=(10, 6), dpi=300)
    ax = sns.boxplot(x=celltype_cols, y='entropy', data=high_entropy_df)
    ax.set_xticklabels(ax.get_xticklabels(),rotation = 30)
    plt.tight_layout()
    plt.savefig(quantile_dir+os.sep+'Entropy_per_celltype_highentropy.png')

    with open(quantile_dir+os.sep+"MLP_entropy_metrics.txt", 'w') as f:
        f.write('**Performance on low entropy cells**\n')
        f.write("Acc:%s\n" % str(metrics.accuracy_score(low_entropy_df[celltype_cols], low_entropy_df[pred_celltype_cols])))
        f.write("ARI:%s\n" % str(metrics.cluster.adjusted_rand_score(low_entropy_df[celltype_cols], low_entropy_df[pred_celltype_cols])))
        f.write("macroF1:%s\n\n" % str(metrics.f1_score(low_entropy_df[celltype_cols], low_entropy_df[pred_celltype_cols], average="macro")))
        f.write('**Performance on high entropy cells**\n')
        f.write("Acc:%s\n" % str(metrics.accuracy_score(high_entropy_df[celltype_cols], high_entropy_df[pred_celltype_cols])))
        f.write("ARI:%s\n" % str(metrics.cluster.adjusted_rand_score(high_entropy_df[celltype_cols], high_entropy_df[pred_celltype_cols])))
        f.write("macroF1:%s\n\n" % str(metrics.f1_score(high_entropy_df[celltype_cols], high_entropy_df[pred_celltype_cols], average="macro")))
        f.write('**Cell labels in low entropy cells**\n')
        f.write(str(low_entropy_df[pred_celltype_cols].value_counts()))
        f.write('\n\n')
        f.write('**Number of cells selected as low entropy cells: %d**\n\n' % low_entropy_df.shape[0])

    ## === Step2: take low entropy cells as reference
    test_adata.obs.rename(columns={pred_celltype_cols:'first_pred_celltypes'}, inplace=True)
    test_adata.obs[pred_celltype_cols] = test_adata.obs['first_pred_celltypes']

    if FS:
        test_ref_adata = test_adata[low_entropy_cells]
        test_tgt_adata = test_adata[high_entropy_cells]
        new_args = args
        new_args.n_features = 5000
        test_ref_adata, test_tgt_adata = dataloading_utils.process_loaded_data(
                test_ref_adata, test_tgt_adata, quantile_dir, args=new_args,
                preprocess=True, save_data=False)
    else:
        test_ref_adata = test_fs_adata[low_entropy_cells]
        test_tgt_adata = test_fs_adata[high_entropy_cells]
 
    enc, x_test_ref, y_test_ref, x_test_tgt, y_test_tgt = \
            prepare_data(test_ref_adata, test_tgt_adata, 
                    enc=enc, train_celltype_cols=pred_celltype_cols, test_celltype_cols=celltype_cols)

    #method_utils.run_pipeline(args, train_adata, test_adata, data_dir, result_dir)
    n_clusters = len(set(test_ref_adata.obs[pred_celltype_cols]))
    print("\n\n=== MLP\n")
    MLP_dims = [128, 32, 8] if test_ref_adata.shape[0] >=5000 else [64, 16]
    y_pred_tgt, exec_time = run_MLP(x_test_ref, y_test_ref, x_test_tgt, 
            dims=MLP_dims, batch_size=batch_size, seed=0, 
            save_dir=quantile_dir+os.sep+'secondround_checkpoint') ## run MLP
    print("\n\n=== Run time:", exec_time)

    y_pred_tgt_label = y_pred_tgt.argmax(1)
    y_test_tgt_label = y_test_tgt.argmax(1)

    pred_onehot = np.zeros((y_pred_tgt_label.size, n_clusters))
    pred_onehot[np.arange(y_pred_tgt_label.size), y_pred_tgt_label] = 1
    pred_celltypes = enc.inverse_transform(pred_onehot)
    print("=== Predicted celltypes: ", set(pred_celltypes.flatten()))
    test_adata.obs.loc[high_entropy_cells, pred_celltype_cols] = pred_celltypes

    with open(quantile_dir+os.sep+"MLP_entropy_metrics.txt", 'a') as f:
        f.write('**Performance on high entropy (second round)**\n')
        f.write("Acc:%s\n" % str(metrics.accuracy_score(y_test_tgt_label, y_pred_tgt_label)))
        f.write("ARI:%s\n" % str(metrics.cluster.adjusted_rand_score(y_test_tgt_label, y_pred_tgt_label)))
        f.write("macroF1:%s\n\n" % str(metrics.f1_score(y_test_tgt_label, y_pred_tgt_label, average="macro")))
        f.write('**Performance on all (second round)**\n')
        f.write("Acc:%s\n" % str(metrics.accuracy_score(test_adata.obs[celltype_cols], test_adata.obs[pred_celltype_cols])))
        f.write("ARI:%s\n" % str(metrics.cluster.adjusted_rand_score(test_adata.obs[celltype_cols], test_adata.obs[pred_celltype_cols])))
        f.write("macroF1:%s\n\n" % str(metrics.f1_score(test_adata.obs[celltype_cols], test_adata.obs[pred_celltype_cols], average="macro")))

    sc.pl.tsne(test_adata, 
            color=[celltype_cols, pred_celltype_cols, "first_pred_celltypes", "assigned_status", "entropy", "entropy_status"], 
            size=15, ncols=2, wspace=0.5)
    plt.savefig(quantile_dir+os.sep+'tSNE_testdata_entropy.png')

def run_iterative_MLP_pipeline(args, train_adata, test_adata, data_dir, result_dir,
        batch_size=32, iterative_rounds=10, celltype_cols="cell.type", pred_celltype_cols="pred_celltypes", 
        low_entropy_quantile=0.2):
    '''Iterative prediction
    - Select low entropy cells from each cell type
    - Use low entropy cells as reference to predict the rest target cells
    '''
    ## first step prediction
    y_pred, enc = run_pipeline(args, train_adata, test_adata, data_dir, result_dir)

    ## iteratively predict
    tmp_train_adata, tmp_test_adata = train_adata.copy(), test_adata.copy()
    tmp_test_adata.obs.rename(columns={pred_celltype_cols: pred_celltype_cols+'0'}, inplace=True)
    y_pred_entropy = [-np.nansum(y_pred[i]*np.log(y_pred[i])) for i in range(y_pred.shape[0])]
    tmp_test_adata.obs['entropy0'] = y_pred_entropy

    result_test_adata = tmp_test_adata.copy()
    for i in range(iterative_rounds):
        print("**Round %d" % i)
        low_entropy_cells = []
        for celltype in set(result_test_adata.obs[pred_celltype_cols+str(i)]):
            celltype_df = result_test_adata.obs[result_test_adata.obs[pred_celltype_cols+str(i)] == celltype]
            entropy_cutoff = np.quantile(celltype_df['entropy'+str(i)], q=low_entropy_quantile)
            cells = celltype_df.index[np.where(celltype_df['entropy'+str(i)] <= entropy_cutoff)[0]].tolist()
            low_entropy_cells.extend(cells)
        high_entropy_cells = list(set(result_test_adata.obs_names) - set(low_entropy_cells))

        tmp_train_adata = result_test_adata[low_entropy_cells]
        tmp_test_adata = result_test_adata[high_entropy_cells]
        enc, x_test_ref, y_test_ref, x_test_tgt, y_test_tgt = \
                prepare_data(tmp_train_adata, tmp_test_adata,
                        enc=enc, train_celltype_cols=pred_celltype_cols+str(i), 
                        test_celltype_cols=celltype_cols)
        n_clusters = len(set(tmp_train_adata.obs[pred_celltype_cols+str(i)]))
        MLP_dims = [128, 32, 8] if tmp_train_adata.shape[0] >=5000 else [64, 16]
        y_pred_new, exec_time = run_MLP(x_test_ref, y_test_ref, x_test_tgt,
                dims=MLP_dims, batch_size=batch_size, seed=0, 
                save_dir=result_dir+os.sep+str(low_entropy_quantile)+'_'+str(i)+'_checkpoint')  ## run MLP
        y_pred_new_entropy = [-np.nansum(y_pred_new[i]*np.log(y_pred_new[i])) for i in range(y_pred_new.shape[0])]

        ## set entropy
        result_test_adata.obs['entropy'+str(i+1)] = result_test_adata.obs['entropy'+str(i)]
        result_test_adata.obs.loc[high_entropy_cells, 'entropy'+str(i+1)] = y_pred_new_entropy
        ## set predicted cell types
        result_test_adata.obs[pred_celltype_cols+str(i+1)] = result_test_adata.obs[pred_celltype_cols+str(i)]
        y_pred_new_label = y_pred_new.argmax(1)
        pred_onehot = np.zeros((y_pred_new_label.size, n_clusters))
        pred_onehot[np.arange(y_pred_new_label.size), y_pred_new_label] = 1
        pred_celltypes = enc.inverse_transform(pred_onehot)
        result_test_adata.obs.loc[high_entropy_cells, pred_celltype_cols+str(i+1)] = pred_celltypes

    with open(result_dir+os.sep+str(low_entropy_quantile)+'_'+str(iterative_rounds)+"_MLP_entropy_metrics.txt", 'w') as f:
        for i in range(iterative_rounds+1):
            f.write("**Round %d \n" % i)
            f.write("Acc:%s\n" % str(metrics.accuracy_score(result_test_adata.obs[celltype_cols], 
                result_test_adata.obs[pred_celltype_cols+str(i)])))
            f.write("ARI:%s\n" % str(metrics.cluster.adjusted_rand_score(result_test_adata.obs[celltype_cols],
                result_test_adata.obs[pred_celltype_cols+str(i)])))
            f.write("macroF1:%s\n\n" % str(metrics.f1_score(result_test_adata.obs[celltype_cols], 
                result_test_adata.obs[pred_celltype_cols+str(i)], average="macro")))
    result_test_adata.obs.to_csv(result_dir+os.sep+str(low_entropy_quantile)+'_'+str(iterative_rounds)+'_iterative_obs.csv')


def prepare_data(train_adata, test_adata, enc=None, 
        train_celltype_cols="cell.type", test_celltype_cols="cell.type"):
    '''Prepare data for train/test
    '''
    if scipy.sparse.issparse(train_adata.X):
        x_train = train_adata.X.toarray()
    else:
        x_train = train_adata.X

    if scipy.sparse.issparse(test_adata.X):
        x_test = test_adata.X.toarray()
    else:
        x_test = test_adata.X
 
    if enc is None:
        enc = OneHotEncoder(handle_unknown='ignore')
        y_train = enc.fit_transform(train_adata.obs[[train_celltype_cols]]).toarray()
    else:
        y_train = enc.transform(train_adata.obs[[train_celltype_cols]]).toarray()
    y_test = enc.transform(test_adata.obs[[test_celltype_cols]]).toarray()
    return enc, x_train, y_train, x_test, y_test
