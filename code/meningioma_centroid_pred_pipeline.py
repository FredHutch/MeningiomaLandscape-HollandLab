# Nicholas Nuechterlein
import numpy as np
import pandas as pd
from tqdm import tqdm

from meningioma_utils import get_norm_umap,  get_umap_dist, score_centroid_pred, normalize_training_umap
from meningioma_pipeline_helpers import get_score_helper

def stage1_get_nn_bag_dict_new(umap_embeddings, umap_df, num_nn=100, radius=0.5, normalize=False, num_umaps=100):    
    '''
    nn_picked_dict = number of NNs per each patient on the static UMAP and each of the background UMAPs.
    nn_bag_dict
    '''
    avg_dist_list = []
    nn_picked_dict = {}
    nn_bag_dict = {}
    for patient in tqdm(umap_df.index):
        nns_list = []
        nn_picked_dict[patient] = {}
        for random_state in list(umap_embeddings.keys())[:num_umaps]:

            train_umap_df = umap_embeddings[random_state][patient]['train_umap_dict']
            test_umap_df = umap_embeddings[random_state][patient]['test_umap_dict']

            ##
            train_umap_df, mean, dist = normalize_training_umap(train_umap_df)
            test_umap_df = (test_umap_df - mean)/dist

            # data_umap_df = pd.concat([train_umap_df, test_umap_df])
            
            # if normalize:
            #     data_umap_df = get_norm_umap(umap_df=data_umap_df) # normalize UMAP

            # train_umap_df = data_umap_df.loc[train_umap_df.index]
            # test_umap_df = data_umap_df.loc[test_umap_df.index]

            # get nns
            x, y = test_umap_df.values[0]
            values = np.sum((train_umap_df - np.array([x, y]))**2, axis=1) # how far everything is from x,y
            sup_dist_df = pd.DataFrame(data=values, index=train_umap_df.index, columns=['dist']).sort_values('dist')

            # use radius to cut down nns (usually when sample is in a separate cluster)
            sup_dist_df = sup_dist_df.loc[sup_dist_df['dist'] < radius]

            # of those samples in the radius, take only the top num_nn of them
            pred_nn_idxs = sup_dist_df['dist'].sort_values().index.tolist()[:num_nn]
            nn_picked_dict[patient][random_state] = len(pred_nn_idxs)

            # add nns to nn bag
            nns_list = nns_list + pred_nn_idxs

        nn_bag_dict[patient] = nns_list
    
    return nn_bag_dict

def stage2_get_per_common_nn_df(nn_bag_dict, umap_df, num_background_itrs=100):
    df_list = [pd.DataFrame(index=umap_df.index)]
    for i, patient in enumerate(tqdm(nn_bag_dict.keys())):
        nn_bag_list = nn_bag_dict[patient]

        vals, cnts = np.unique(nn_bag_list, return_counts=True)
        value_counts = {vals[i]:cnts[i] for i in range(vals.shape[0])}
        patient_value_counts_df = pd.DataFrame.from_dict(data=value_counts, orient='index', columns=[patient])
        df_list.append(patient_value_counts_df)

    # the columns of this df are the colors of the nns on the plot
    per_common_nn_df = pd.concat(df_list, axis=1)/num_background_itrs
    
    return per_common_nn_df

def stage3_get_score(per_common_nn_df, umap_df, 
                     remove_low_freq_nns = True, 
                     remove_outliers = True, 
                     remove_radius = True,
                     mod_prec_thresh = 0.25,
                     outlier_quantile = 0.95,
                     radius = 0.75):

    score_dict = {}
    scale_df_dict = {}
    score_dict_for_patient = {}
    for patient in tqdm(per_common_nn_df.columns):

        score_dict[patient], scale_df_dict, score_dict_for_patient[patient] = get_score_helper(patient=patient,
                                               scale_df_dict=scale_df_dict,
                                               per_common_nn_df=per_common_nn_df, 
                                               umap_df=umap_df, 
                                               remove_low_freq_nns = remove_low_freq_nns, 
                                               remove_outliers = remove_outliers, 
                                               remove_radius = remove_radius,
                                               mod_prec_thresh = mod_prec_thresh,
                                               outlier_quantile = outlier_quantile,
                                               radius = radius)

    score_df = pd.DataFrame.from_dict(data=score_dict_for_patient, orient='index').sort_values('nn_score')
    score_df['per_diff_num_nn'] = (score_df['init_num_nn'] - score_df['num_nn_averaged'])/score_df['init_num_nn']

    return score_df, scale_df_dict