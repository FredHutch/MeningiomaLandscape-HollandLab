# Nicholas Nuechterlein
import pandas as pd
import numpy as np

def get_centroid(umap_df, scale_df):
    '''
    takes weighted average of umap_df points and returns the centroid
    umap_df: x,y coordiance of reference landscape
    scale_df: the weight each reference sample should contribue to centroid calculation

    '''
    
    temp_umap_df = pd.DataFrame(index=umap_df.index)
    temp_umap_df['umap_1'] = umap_df['umap_1'].multiply(scale_df)
    temp_umap_df['umap_2'] = umap_df['umap_2'].multiply(scale_df)

    x_cent, y_cent = temp_umap_df.sum()/scale_df.sum()
    
    return (x_cent, y_cent)

def load_reference_umap(umap_filepath='../data/reference_umap_coordinates.csv', normalize = False):
    umap_df = pd.read_csv(umap_filepath)
    umap_df = umap_df.set_index('coordinate_ID')

    rename_dict = {'UMAP1_2D':'umap_1', 'UMAP2_2D':'umap_2'}
    umap_df = umap_df.rename(columns=rename_dict)
    umap_df = umap_df[['umap_1', 'umap_2']]
    
    if normalize:
        umap_df = get_norm_umap(umap_df)
    
    return umap_df

def get_norm_umap(umap_df):
    # center and scale
    umap_df = umap_df - umap_df.mean()

    avg_dist = (umap_df**2).sum(axis=1).mean()
    umap_df = umap_df/np.sqrt(avg_dist)
    
    return umap_df

def normalize_training_umap(umap_df):
    # center and scale
    mean = umap_df.mean()
    umap_df = umap_df - mean

    avg_dist = (umap_df**2).sum(axis=1).mean()
    umap_df = umap_df/np.sqrt(avg_dist)
    
    return umap_df, mean, np.sqrt(avg_dist)

def get_umap_dist(umap_df):
    # center
    umap_df = umap_df - umap_df.mean()
    avg_dist = (umap_df**2).sum(axis=1).mean()
    return avg_dist
  

def get_distance_df(umap_df, x_cent, y_cent, scale_df):
    '''
    umap_df: x,y coordiance of reference landscape
    x_cent, y_cent: x,y cooridenates of centrid
    scale_df: the weight each reference sample should contribue to centroid calculation
    ---
    I think this function just scores the results (i.e., not part of the prediction)
    '''
    
    # find patients that were consitered NNs (were used to compute centroid)
    patient_nn_idxs = scale_df.loc[scale_df != 0].index.tolist()
    _umap_df = umap_df.loc[patient_nn_idxs]
    
    # demean UMAP
    _umap_df['umap_1'] = _umap_df['umap_1'] - x_cent
    _umap_df['umap_2'] = _umap_df['umap_2'] - y_cent
    
    # dist score 1
    score_avg_df = (_umap_df**2).sum(axis=1)
    
    # dist score 2 (weighted)
    weighted_score_avg_df = ((_umap_df**2).multiply(scale_df, axis=0)).sum(axis=1)
    
    return score_avg_df, weighted_score_avg_df

def load_meningioma_umap(umap_path, normalize = False):
    umap_df = pd.read_csv(umap_path)
    umap_df = umap_df.set_index('coordinate_ID')

    rename_dict = {'UMAP1_2D':'umap_1', 'UMAP2_2D':'umap_2'}
    umap_df = umap_df.rename(columns=rename_dict)
    umap_df = umap_df[['umap_1', 'umap_2']]
    
    if normalize:
        umap_df = get_norm_umap(umap_df)
    
    return umap_df

def score_centroid_pred(patient, x_cent, y_cent, scale_df, umap_df, patient_mod_perc_common_nn_df, base_x_cent, base_y_cent):
    score_avg_df, weighted_score_avg_df = get_distance_df(umap_df=umap_df, 
                                                          x_cent=x_cent, 
                                                          y_cent=y_cent, 
                                                          scale_df=scale_df)
    # dist 1
    score_avg_dist = score_avg_df.mean()
    
    # dist 2
    weighted_score_avg_dist = weighted_score_avg_df.mean()
    
    # dist 3
    x_gt, y_gt = umap_df.loc[patient]['umap_1'], umap_df.loc[patient]['umap_2']
    cent_dist = np.sqrt((x_gt-x_cent)**2 + (y_gt-y_cent)**2)
    
    base_cent_dist = np.sqrt((x_gt-base_x_cent)**2 + (y_gt-base_y_cent)**2)

    
    # nn score
    nn_score = (1 - patient_mod_perc_common_nn_df).mean()
    
    score_dict_for_patient = {'nn_score':nn_score,
                              'cent_dist':cent_dist,
                              'avg_dist':score_avg_dist, 
                              'avg_dist_weighted':weighted_score_avg_dist,
                              'x_cent':x_cent,
                              'y_cent':y_cent,
                              'base_cent_dist':base_cent_dist,
                              'num_nn_averaged':len(patient_mod_perc_common_nn_df)}

    return score_dict_for_patient