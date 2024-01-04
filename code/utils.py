import numpy as np
import pandas as pd

from sklearn.metrics import roc_auc_score, matthews_corrcoef

def round_pval(pval, stars=False, sig_level=0.05):
    if not stars:
        for i in range(2, 500):
            if np.round(pval, i) > 0:
                return np.round(pval, i)
    if stars:
        if pval > sig_level:
            return 'ns'
        elif pval > 0.01:
            return '*'
        elif pval > 0.001:
            return '**'
        elif pval > 0.0001:
            return '***'
        else:
            return '****'

    else:
        return pval


def get_avg_mcc(y_true, y_preds, y_probs):
    zero_dict = {x:0 for x in y_true.unique()}
    
    metric_per_class_list = []
    metric_weigted_per_class_list = []
    for idx in y_true.unique():
        num_class = len(y_true.loc[y_true == idx])

        y_true_replace_dict = zero_dict.copy()
        y_true_replace_dict[idx] = 1

        _y_true = y_true.replace(y_true_replace_dict) # replace current class with 1, others with 0
        _y_probs = y_probs[idx].T.loc[_y_true.index] 

        y_probs_rounded = np.rint(_y_probs) 

        mcc = matthews_corrcoef(y_true=_y_true, y_pred=y_probs_rounded)
        mcc_weight = mcc*(num_class/len(y_true))

        metric_per_class_list.append(mcc)
        metric_weigted_per_class_list.append(mcc_weight)

    return {'macro':np.mean(metric_per_class_list), 'micro':np.sum(metric_weigted_per_class_list)}

def get_avg_auc(y_true, y_probs):
    metric_per_class_list = []
    metric_weigted_per_class_list = []
    for target in y_true.unique():
        num_class = len(y_true.loc[y_true == target])

        replace_dict = {x:0 for x in y_true.unique() if x != target}
        replace_dict[target] = 1
        temp_y_true = y_true.replace(replace_dict)
        auc = roc_auc_score(temp_y_true, y_probs[target])
        metric_per_class_list.append(auc)
        
        auc_weight = auc*(num_class/len(y_true))
        metric_weigted_per_class_list.append(auc_weight)
    
    return {'macro':np.mean(metric_per_class_list), 'micro':np.sum(metric_weigted_per_class_list)}

def get_avg_acc(y_true, y_preds):
    metric_per_class_list = []
    metric_weigted_per_class_list = []
    acc_local_list = []
    for target in y_true.unique():
        num_class = len(y_true.loc[y_true == target])
        idxs = y_true.loc[y_true == target].index.tolist()

        replace_dict = {x:0 for x in y_true.unique() if x != target}
        replace_dict[target] = 1
        temp_y_true = y_true.replace(replace_dict)
        temp_y_preds = y_preds.replace(replace_dict)
        
        true_pred_df = pd.concat([temp_y_true, temp_y_preds], axis=1).loc[idxs]
        col1, col2 = true_pred_df.columns[0], true_pred_df.columns[1]
        
        acc = len(true_pred_df.loc[true_pred_df[col1] == true_pred_df[col2]])/len(true_pred_df)
        metric_per_class_list.append(acc)

        auc_weight = acc*(num_class/len(y_true))
        metric_weigted_per_class_list.append(auc_weight)
        
    return {'macro':np.mean(metric_per_class_list), 'micro':np.mean(y_true == y_preds)}