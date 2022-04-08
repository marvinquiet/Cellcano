## load library for reading Anndata H5AD file and selecting top features
#suppressMessages(library(FEAST))
suppressMessages(library(matrixStats))

## code from FEAST package: https://github.com/suke18/FEAST/blob/85102e6f9c91e88143ba2a1daebbe8422359d794/R/Fisher_F.R
cal_F2 = function(Y, classes){
    if (all(Y%%1 == 0)){
        L = colSums(Y) / median(colSums(Y))
        Ynorm = log2(sweep(Y, 2, L, FUN="/") + 1)
    }else{Ynorm = Y}
    cid = which(is.na(classes))
    if (length(cid) > 0){
        Ynorm = Ynorm[, -cid]
        classes = classes[-cid]
    }
    classes = as.factor(classes)
    unique_classes = unique(classes)
    k = length(unique(classes))
    row_class_mean = matrix(0, ncol = k, nrow = nrow(Y))

    row_mean = rowMeans(Ynorm)
    k = length(unique(classes))
    pb = txtProgressBar( min = 0, max = k, style = 3)
    for (i in seq_len(k)){
        setTxtProgressBar(pb, i)
        ix = which(classes == unique_classes[i])
        if (length(ix) > 1){
            tmp_mean =  rowMeans(Ynorm[, ix])
            row_class_mean[, i] = tmp_mean
        }else{
            row_class_mean[, i] = Ynorm[, ix]
        }
    }; close(pb)
    colnames(row_class_mean) = unique_classes
    ### make sure the classes are matched; otherwise, causing error ###
    table_class = table(classes)[unique_classes]
    BBS = table_class  %*% t((row_class_mean - row_mean)^2)
    BBS = BBS[1,]
    TSS = rowSums((Ynorm - row_mean) ^2)
    ESS = TSS - BBS
    df1 = k-1; df2 = ncol(Ynorm)-k
    F_scores = (BBS/df1) / (ESS/ df2)
    ps = pf(F_scores, df1, df2, lower.tail = FALSE)
    return(list(F_scores = as.numeric(F_scores), ps = ps))
}


args <- commandArgs(trailingOnly = TRUE)
df_dir <- args[1]
annotation_dir <- args[2]
n_features <- args[3]

processed_counts <- read.csv(df_dir, row.names=1)
#processed_counts <- process_Y(counts, thre=0)
#rm(counts) # release memory

cell_annotations <- readLines(annotation_dir)
F_res <- cal_F2(processed_counts, cell_annotations)$F_scores

names(F_res) <- rownames(processed_counts)
ixs <- order(F_res, decreasing=T)[1:n_features]
features <- names(F_res[ixs])
write(features, file.path(dirname(df_dir), "F-test_features.txt"))
