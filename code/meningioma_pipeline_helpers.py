import numpy as np
import pandas as pd
from tqdm import tqdm

from meningioma_utils import get_distance_df, get_centroid, score_centroid_pred

def step1_remove_low_freq_nns_fn(patient, scale_df, scale_df_dict, patient_mod_perc_common_nn_df, mod_prec_thresh):
    # update values
    patient_mod_perc_common_nn_df = patient_mod_perc_common_nn_df.loc[patient_mod_perc_common_nn_df 
                                                                  >= mod_prec_thresh]
    scale_df[scale_df < mod_prec_thresh] = 0
    scale_df_dict[patient]['1-low-feq'] = scale_df.copy()
    
    return scale_df, scale_df_dict

def step2_remove_outliers_fn(patient, umap_df, scale_df, scale_df_dict, outlier_quantile, patient_mod_perc_common_nn_df):
    base_x_cent, base_y_cent = get_centroid(umap_df=umap_df, scale_df=scale_df)

    # get distances
    _score_avg_df, _weighted_score_avg_df = get_distance_df(umap_df=umap_df, 
                                                            x_cent=base_x_cent, 
                                                            y_cent=base_y_cent, 
                                                            scale_df=scale_df)

    # find and remove outliers
    inliner_idxs = _score_avg_df[_score_avg_df.between(0, _score_avg_df.quantile(outlier_quantile))].index.tolist()
    outlier_idxs = [x for x in _score_avg_df.index if x not in inliner_idxs]

    # update values
    patient_mod_perc_common_nn_df = patient_mod_perc_common_nn_df.loc[patient_mod_perc_common_nn_df.index.isin(inliner_idxs)]
    scale_df.loc[outlier_idxs] = 0 
    scale_df_dict[patient]['2-outliers'] = scale_df.copy()
    
    return (scale_df, scale_df_dict, patient_mod_perc_common_nn_df, base_x_cent, base_y_cent)

def step3_remove_radius_fn(patient, radius, umap_df, scale_df, scale_df_dict, patient_mod_perc_common_nn_df): 
    _x_cent, _y_cent = get_centroid(umap_df=umap_df, scale_df=scale_df)

    values = np.sum((umap_df - np.array([_x_cent, _y_cent]))**2, axis=1)
    dist_df = pd.DataFrame(data=values, index=umap_df.index, columns=['dist']).sort_values('dist')

    radius_nn_idxs = dist_df.loc[dist_df['dist'] < radius].index.tolist()

    keep_idxs = [x for x in patient_mod_perc_common_nn_df.index if x in radius_nn_idxs]
    patient_mod_perc_common_nn_df = patient_mod_perc_common_nn_df.loc[keep_idxs]

    exclude_idxs = [x for x in scale_df.index if x not in radius_nn_idxs]
    scale_df.loc[exclude_idxs] = 0
    scale_df_dict[patient]['3-radius'] = scale_df.copy()
    
    return scale_df, scale_df_dict


def get_score_helper(patient,
                     scale_df_dict,
                     per_common_nn_df, 
                     umap_df, 
                     remove_low_freq_nns = True, 
                     remove_outliers = True, 
                     remove_radius = True,
                     mod_prec_thresh = 0.25,
                     outlier_quantile = 0.95,
                     radius = 0.75):
    
    scale_df_dict[patient]  = {}
    
    ### Stage 1 ###
    patient_mod_perc_common_nn_df = per_common_nn_df[patient]
    init_num_nn = len(patient_mod_perc_common_nn_df.dropna())
    init_nn_score = (1 - patient_mod_perc_common_nn_df).mean()
    
    scale_df = per_common_nn_df[patient].replace({np.nan:0})
    scale_df_dict[patient]['baseline'] = scale_df.copy()
    
    # 1) remove low fequency NNs
    if remove_low_freq_nns:
        scale_df, scale_df_dict = step1_remove_low_freq_nns_fn(patient=patient,
                                                         scale_df=scale_df, 
                                                         scale_df_dict=scale_df_dict, 
                                                         patient_mod_perc_common_nn_df=patient_mod_perc_common_nn_df, 
                                                         mod_prec_thresh=mod_prec_thresh)
    
    # 2) trim distant points 
    if remove_outliers:
        output = step2_remove_outliers_fn(patient=patient,
                                    umap_df=umap_df, 
                                    scale_df=scale_df, 
                                    scale_df_dict=scale_df_dict, 
                                    outlier_quantile=outlier_quantile, 
                                    patient_mod_perc_common_nn_df=patient_mod_perc_common_nn_df)
        
        scale_df, scale_df_dict, patient_mod_perc_common_nn_df, base_x_cent, base_y_cent = output
        
    # 3) remove points outside of radius
    if remove_radius:
        scale_df, scale_df_dict = step3_remove_radius_fn(patient=patient, 
                                                   radius=radius, 
                                                   umap_df=umap_df, 
                                                   scale_df_dict=scale_df_dict,
                                                   scale_df=scale_df, 
                                                   patient_mod_perc_common_nn_df=patient_mod_perc_common_nn_df)
    
    
    ## Score ##
    x_cent, y_cent = get_centroid(umap_df=umap_df, scale_df=scale_df)
    
    score_dict_for_patient = score_centroid_pred(patient=patient,
                                                 x_cent=x_cent, 
                                                 y_cent=y_cent, 
                                                 base_x_cent=base_x_cent,
                                                 base_y_cent=base_y_cent,
                                                 scale_df=scale_df, 
                                                 umap_df=umap_df, 
                                                 patient_mod_perc_common_nn_df=patient_mod_perc_common_nn_df)
    
    score_dict_for_patient['init_num_nn'] = init_num_nn
    score_dict_for_patient['init_nn_score'] = init_nn_score

    return scale_df, scale_df_dict, score_dict_for_patient