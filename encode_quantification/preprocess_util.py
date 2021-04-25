from scipy import stats
import numpy as np
import pandas as pd
from itertools import combinations
from sklearn.metrics import precision_score,recall_score,accuracy_score,roc_auc_score,f1_score,roc_curve,precision_recall_curve
from static_data import *

def normalize(arr):
    return arr
    # return arr/arr.sum()
def get_consistency_measures(abund_arr,K):
    condition_column_split_index = abund_arr.shape[1]//2
    arr_lst = []
    for abund_arr in [abund_arr[:,:condition_column_split_index],abund_arr[:,condition_column_split_index:]]:
        for i, j in combinations(range(abund_arr.shape[1]), 2):
            arr_lst.append((((abund_arr[:,i]<K) & (abund_arr[:,j]>=K)) | ((abund_arr[:,i]>=K) & (abund_arr[:,j]<K))))
    return np.vstack(arr_lst).mean()
def get_resolution_entropy(abund_arr,I):
    df = pd.DataFrame({'abund':abund_arr})
    df['group_range'] = pd.cut(df['abund'],I,right=False)
    P = df.groupby('group_range').count()['abund'].values / abund_arr.shape[0]
    P = P[P != 0]
    RE = -np.mean(np.log(P)*P)
    return RE
def get_single_sample_metric(metric,ground_truth, estimated,df=None,agg_gene=False):
    if metric == 'spearmanr':
        if (ground_truth.shape[0] == 1):
            return 1 if ground_truth.values[0] == estimated.values[0] else 0
        else:
            return stats.spearmanr(ground_truth, estimated).correlation
    elif metric == 'nrmse':
        if (np.std(ground_truth, axis=0) == 0):
            return np.square(estimated-ground_truth).mean(axis=0)
        else:
            return np.square(estimated-ground_truth).mean(axis=0) / np.std(ground_truth, axis=0)
    elif metric == 'mrd':
        if (np.count_nonzero(np.linalg.norm(ground_truth)) == 0):
            return np.mean(np.linalg.norm(estimated-ground_truth))
        if (agg_gene):
            df['mrd'] = np.linalg.norm(estimated-ground_truth) / np.linalg.norm(ground_truth)
            return np.mean(df[['mrd','gene']].groupby('gene').mean()).values[0]
        return  np.mean(np.linalg.norm(estimated-ground_truth) / np.linalg.norm(ground_truth))
    elif metric == 'RE':
        return get_resolution_entropy(estimated,100)
    elif metric == 'mean_arr':
        return df['arr'].mean()
def get_multi_sample_metric(metric,df,ground_truth, estimated,agg_gene=False):
    estimation_diff_expressed = df['alfc'] > df['alfc'].mean()
    if (ground_truth is not None):
        truth_diff_expressed = df['true_alfc'] > df['true_alfc'].mean()
    else:
        truth_diff_expressed = None
    if metric == 'precision':
        return precision_score(truth_diff_expressed,estimation_diff_expressed,zero_division=0)
    elif metric == 'recall':
        return recall_score(truth_diff_expressed,estimation_diff_expressed,zero_division=0)
    elif metric == 'accuracy':
        return accuracy_score(truth_diff_expressed,estimation_diff_expressed)
    elif metric == 'auc':
        try:
            return roc_auc_score(truth_diff_expressed,estimation_diff_expressed)
        except:
            return 0
    elif metric == 'f1':
        return f1_score(truth_diff_expressed,estimation_diff_expressed,zero_division=0)
    elif metric == 'spearmanr':
        if (ground_truth.shape[0] == 1):
            return 1 if ground_truth.values[0] == estimated.values[0] else 0
        else:
            return stats.spearmanr(ground_truth, estimated,axis=None).correlation
    elif metric == 'nrmse':
        return (np.square(estimated-ground_truth).mean(axis=0) / np.std(ground_truth, axis=0)).mean()
    elif metric == 'mrd':
        if (np.count_nonzero(np.linalg.norm(ground_truth)) == 0):
            return np.mean(np.linalg.norm(estimated-ground_truth))
        if (agg_gene):
            df['mrd'] = np.linalg.norm(estimated-ground_truth) / np.linalg.norm(ground_truth)
            return np.mean(df[['mrd','gene']].groupby('gene').mean()).values[0]
        return  np.mean(np.linalg.norm(estimated-ground_truth) / np.linalg.norm(ground_truth))
    elif metric =='CM':
        # return get_consistency_measures(estimated,np.median(estimated))
        return get_consistency_measures(estimated,1)
    elif metric == 'RM':
        std_1,std_2 = df['std#1'].values,df['std#2'].values
        return np.sqrt(np.mean([np.square(std_1),np.square(std_2)]))
    elif metric == 'RE':
        return np.median([get_resolution_entropy(estimated[:,i],100) for i in range(estimated.shape[1])])
    elif metric == 'mean_arr':
        arr_columns = [x for x in list(df.columns) if 'arr_' in x]
        return df[arr_columns].values.mean()
def preprocess_single_sample_util(df, kvalues_dict,num_exon_dict,isoform_length_dict,isoform_gene_dict):

    df['estimated_abund'] = normalize(df['estimated_abund'])
    df['true_abund'] = normalize(df['true_abund'])
    df['std'] = (df['true_abund'] - df['estimated_abund']) / 2
    df['arr'] = df['estimated_abund'] / df['true_abund']
    df.loc[np.isinf(df['arr']),'arr'] = 0
    df['log2_true_abund'] = np.log2(df['true_abund']+1)
    df = df.set_index('isoform').assign(K_value=pd.Series(kvalues_dict),num_exons=pd.Series(num_exon_dict),isoform_length=pd.Series(isoform_length_dict),gene=pd.Series(isoform_gene_dict)).reset_index()
    df = df.dropna()
    return df
def preprocess_multi_sample_diff_condition_util(estimated_df,true_expression_df, kvalues_dict,num_exon_dict,isoform_length_dict,isoform_gene_dict):
    condition_column_split_index = estimated_df.shape[1]//2
    estimated_arr = estimated_df.drop(0,axis=1).to_numpy()
    estimated_arr = normalize(estimated_arr)
    df = pd.DataFrame()
    df['isoform'] = estimated_df[0]
    df = df.set_index('isoform').assign(K_value=pd.Series(kvalues_dict),num_exons=pd.Series(num_exon_dict),isoform_length=pd.Series(isoform_length_dict),gene=pd.Series(isoform_gene_dict)).reset_index()
    df['std'] = np.std(np.log2(estimated_arr+1),axis=1)
    df['std#1'] = np.std(np.log2(estimated_arr[:,:condition_column_split_index]+1),axis=1)
    df['std#2'] = np.std(np.log2(estimated_arr[:,condition_column_split_index:]+1),axis=1)
    df['ave_estimated_abund'] = np.log2(estimated_arr+1).mean(axis=1)
    df['ave_estimated_abund#1'] = np.log2(estimated_arr[:,:condition_column_split_index]+1).mean(axis=1)
    df['ave_estimated_abund#2'] = np.log2(estimated_arr[:,condition_column_split_index:]+1).mean(axis=1)
    df['alfc'] = np.log2((estimated_arr[:,condition_column_split_index:] + 1).mean(axis=1)/(estimated_arr[:,:condition_column_split_index] + 1).mean(axis=1))
    df  = pd.concat((df,pd.DataFrame(estimated_arr).add_prefix('estimated_abund_')),axis=1)
    if (true_expression_df is not None):
        true_expression_arr = true_expression_df.drop(0,axis=1).to_numpy()
        true_expression_arr = normalize(true_expression_arr)
        df['ave_true_abund'] = np.log2(true_expression_arr+1).mean(axis=1)
        df['ave_true_abund#1'] = np.log2(true_expression_arr[:,:condition_column_split_index]+1).mean(axis=1)
        df['ave_true_abund#2'] = np.log2(true_expression_arr[:,condition_column_split_index:]+1).mean(axis=1)
        df['true_alfc'] = np.log2((true_expression_arr[:,condition_column_split_index:] + 1).mean(axis=1)/(true_expression_arr[:,:condition_column_split_index] + 1).mean(axis=1))
        df  = pd.concat((df,pd.DataFrame(true_expression_arr).add_prefix('true_abund_')),axis=1)
        estimated_columns = [x for x in list(df.columns) if 'estimated_abund_' in x]
        true_columns = [x for x in list(df.columns) if 'true_abund_' in x]
        arrs = []
        for i in range(len(estimated_columns)):
            arr = df[estimated_columns[i]] / df[true_columns[i]]
            arr[np.isinf(arr)] = 0
            arrs.append(arr)
        df  = pd.concat((df,pd.DataFrame(np.array(arrs).T).add_prefix('arr_')),axis=1)
    df = df.dropna()
    return df
def prepare_single_sample_table_metrics(ground_truth, estimated,df):
    metrics = {}
    for metric_dict in single_sample_table_metrics:
        metric = metric_dict['id']
        metrics[metric] = get_single_sample_metric(metric,ground_truth, estimated,df)
    return metrics
def prepare_multi_sample_diff_conditon_table_metrics(df,ground_truth_given):
    estimated_columns = [x for x in list(df.columns) if 'estimated_abund_' in x]
    estimated_arr = df[estimated_columns].to_numpy()
    metrics = {}
    if (ground_truth_given):
        true_columns = [x for x in list(df.columns) if 'true_abund_' in x]
        true_expression_arr = df[true_columns].to_numpy()
        for metric_dict in multi_sample_diff_condition_table_metrics:
            metric = metric_dict['id']
            metrics[metric] = get_multi_sample_metric(metric,df,true_expression_arr,estimated_arr)
    else:
        for metric_dict in multi_sample_diff_condition_without_ground_truth_table_metrics:
            metric = metric_dict['id']
            metrics[metric] = get_multi_sample_metric(metric,df,None,estimated_arr)
    return metrics
def prepare_grouped_violin_data(metric,df,agg_gene=False):
    estimated_columns = [x for x in list(df.columns) if 'estimated_abund_' in x]
    true_columns = [x for x in list(df.columns) if 'true_abund_' in x]
    if metric in ['precision','recall','accuracy','auc','f1']:
        return  get_multi_sample_metric(metric,df,df[true_columns].to_numpy(),df[estimated_columns].to_numpy())
    if metric in ['CM','RM','RE']:
        return  get_multi_sample_metric(metric,df,None,df[estimated_columns].to_numpy())
    all_samples_metric = []
    for i in range(len(estimated_columns)):
        if metric in ['nrmse','mrd','spearmanr']:
            temp = get_multi_sample_metric(metric,df,df[true_columns[i]], df[estimated_columns[i]],agg_gene)
            all_samples_metric.append(temp)
        elif metric in ['mean_arr']:
            arr_columns = [x for x in list(df.columns) if 'arr_' in x]
            temp = df[arr_columns[i]].mean()
            all_samples_metric.append(temp)
    return np.array(all_samples_metric)
def prepare_stats_box_plot_data(df,y_axis_names):
    estimated_columns = [x for x in list(df.columns) if 'estimated_abund_' in x]
    true_columns = [x for x in list(df.columns) if 'true_abund_' in x]
    all_cond_metric_dicts = []
    for start_col,end_col in zip([0,len(estimated_columns)//2],[len(estimated_columns)//2,len(estimated_columns)]):
        metric_dicts = []
        for metric in ['nrmse','mrd','spearmanr']:
            if metric in y_axis_names:
                metric_dict = {}
                vals = []
                for i in range(start_col,end_col):
                    true_column = true_columns[i]
                    estimated_column = estimated_columns[i]
                    vals.append(get_single_sample_metric(metric,df[true_column], df[estimated_column],df))
                metric_dict['Metric'] = metric
                metric_dict['Mean'] = np.mean(vals)
                metric_dict['Error'] = np.std(vals)
                metric_dicts.append(metric_dict)
        for metric in ['RE']:
            if metric in y_axis_names:
                metric_dict = {}
                vals = []
                for i in range(start_col,end_col):
                    estimated_column = estimated_columns[i]
                    vals.append(get_single_sample_metric(metric,None, df[estimated_column],df))
                metric_dict['Metric'] = metric
                metric_dict['Mean'] = np.mean(vals)
                metric_dict['Error'] = np.std(vals)
                metric_dicts.append(metric_dict)
        all_cond_metric_dicts.append(metric_dicts)

    return all_cond_metric_dicts
def prepare_consistency_measure_plot_data(df):
    estimated_columns = [x for x in list(df.columns) if 'estimated_abund_' in x]
    estimated_arr = df[estimated_columns].to_numpy()
    n_bins = 1000
    _,C_ranges = pd.cut(estimated_arr.flatten(), n_bins, retbins=True,include_lowest=True)
    C_ranges[0] = 1
    CM_list = [get_consistency_measures(estimated_arr,K) for K in C_ranges]
    return CM_list,C_ranges
def prepare_corr_box_plot_data(estimated_arr,true_expression_arr,shared_bins=None):
    n_bins = 4
    df = pd.DataFrame({'estimated_abund':estimated_arr})
    if (shared_bins is None):
        df['range'],shared_bins = pd.cut(true_expression_arr,n_bins,retbins=True)
    else:
        df['range'] = pd.cut(true_expression_arr,shared_bins)
    df['true_abund'] = [np.mean([interval.left,interval.right]) if interval is not np.nan else None for interval in df['range']]
    df = df.dropna()
    return df,shared_bins
# def get_multi_sample_same_conditon_metrics(estimated_df,true_expression_df,df):
#     estimated_arr = estimated_df.drop(0,axis=1).to_numpy()
#     spearmanr,NRMSE,MRD = '','',''
#     if (true_expression_df is not None):
#         true_expression_arr = true_expression_df.drop(0,axis=1).to_numpy()
#         spearmanr = stats.spearmanr(true_expression_arr, estimated_arr,axis=None).correlation
#         NRMSE = (np.square(estimated_arr-true_expression_arr).mean(axis=0) / (np.std(true_expression_arr, axis=0))).mean()
#         MRD = np.mean(np.linalg.norm(estimated_arr-true_expression_arr) /
#                     np.linalg.norm(true_expression_arr))
#     RM = np.sqrt(np.square(df['std']).sum())
#     return [{'spearmanr': spearmanr, 'nrmse': NRMSE, 'mrd': MRD,'rm':RM}]
# def preprocess_multi_sample_df_same_condition(estimated_df,true_expression_df,annotation):
#     kvalues_dict = get_kvalues_dict(annotation)
#     estimated_arr = estimated_df.drop(0,axis=1).to_numpy()
#     true_expression_arr = true_expression_df.drop(0,axis=1).to_numpy()
#     df = pd.DataFrame()
#     df['isoform'] = estimated_df[0]
#     df['std'] = np.std(np.log2(estimated_arr+1),axis=1)
#     df['ave_estimated_abund'] = np.log2(estimated_arr+1).mean(axis=1)
#     df = df.set_index('isoform').assign(K_value=pd.Series(kvalues_dict)).reset_index()
#     return df