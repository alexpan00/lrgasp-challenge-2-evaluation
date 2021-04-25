import base64
import io
import pandas as pd
import numpy as np
import time
import zipfile
from app import cache
from library.k_values.main import get_kvalues_dict
from preprocess_util import *

def load_data(contents):
    return pd.read_csv(contents[0], sep='\t',header=None,skiprows=1, comment='#')
    # content_type, content_string = contents.split(',')
    # decoded = base64.b64decode(content_string)
    # df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep='\t',header=None)
    # return df
    
def load_zipped_data(contents):
    if 'zip' in contents[0]:
        list_of_df = []
        method_names = []
        with zipfile.ZipFile(contents[0]) as myzip:
            list_of_files = myzip.namelist()
            for path in list_of_files:
                with myzip.open(path) as myfile:
                    list_of_df.append(pd.read_csv(myfile, sep='\t',skiprows=1,header=None))
                    method_names.append(path.split('.')[0])
        return list_of_df,method_names


def load_annotation(contents,is_long_read=True,K_value_selection='Condition_number'):
    with open(contents[0],'r') as f:
        return get_kvalues_dict(io.StringIO(f.read()),is_long_read,K_value_selection)
    # content_type, content_string = contents.split(',')
    # decoded = base64.b64decode(content_string)
    # return io.StringIO(decoded.decode('utf-8'))

def preprocess_single_sample(list_of_contents,replicate_column,is_long_read=True,K_value_selection='Condition_number'):
    estimated_df = load_data(list_of_contents[0]).set_index(0)[[replicate_column]]
    estimated_df.index.name = 'isoform'
    estimated_df.columns = ['estimated_abund']
    kvalues_dict,num_exon_dict,isoform_length_dict,isoform_gene_dict  = load_annotation(list_of_contents[1],is_long_read,K_value_selection)
    if (list_of_contents[2] is not None):
        true_expression_df = load_data(list_of_contents[2]).set_index(0)[[replicate_column]]
        true_expression_df.index.name = 'isoform'
        true_expression_df.columns = ['true_abund']
        df = estimated_df.join(true_expression_df,on='isoform',how='inner').reset_index()
    else:
        raise Exception('No ground truth data is given')
    anno_df = pd.DataFrame({'K_value':pd.Series(kvalues_dict),'num_exons':pd.Series(num_exon_dict),'isoform_length':pd.Series(isoform_length_dict),'gene':pd.Series(isoform_gene_dict)})
    anno_df.index.name = 'isoform'
    anno_df = anno_df.reset_index()
    df = preprocess_single_sample_util(df, kvalues_dict,num_exon_dict,isoform_length_dict,isoform_gene_dict)
    return df,anno_df


def preprocess_multi_sample_diff_condition(list_of_contents,ground_truth_given,is_long_read=True,K_value_selection='Condition_number'):

    estimated_df = load_data(list_of_contents[0]).set_index(0)
    kvalues_dict,num_exon_dict,isoform_length_dict,isoform_gene_dict  = load_annotation(list_of_contents[1],is_long_read,K_value_selection)

    if (ground_truth_given):
        true_expression_df = load_data(list_of_contents[2]).set_index(0)
        intersected_index = true_expression_df.index.intersection(estimated_df.index)
        estimated_df = estimated_df.loc[intersected_index,:].reset_index()
        true_expression_df = true_expression_df.loc[intersected_index,:].reset_index()
    else:
        estimated_df = estimated_df.reset_index()
        true_expression_df = None
    anno_df = pd.DataFrame({'K_value':pd.Series(kvalues_dict),'num_exons':pd.Series(num_exon_dict),'isoform_length':pd.Series(isoform_length_dict),'gene':pd.Series(isoform_gene_dict)})
    anno_df.index.name = 'isoform'
    anno_df = anno_df.reset_index()
    df = preprocess_multi_sample_diff_condition_util(estimated_df,true_expression_df, kvalues_dict,num_exon_dict,isoform_length_dict,isoform_gene_dict)
    
    return df,anno_df

def preprocess_single_sample_multi_method(list_of_contents,replicate_column,is_long_read=True,K_value_selection='Condition_number'):
    estimated_dfs,method_names = load_zipped_data(list_of_contents[0])
    kvalues_dict,num_exon_dict,isoform_length_dict,isoform_gene_dict  = load_annotation(list_of_contents[1],is_long_read,K_value_selection)
    if (list_of_contents[2] is not None):
        true_expression_df = load_data(list_of_contents[2]).set_index(0)[[replicate_column]]
        true_expression_df.index.name = 'isoform'
        true_expression_df.columns = ['true_abund']
    else:
        raise Exception('No ground truth data is given')
    dfs = []
    for estimated_df in estimated_dfs:
        estimated_df = estimated_df.set_index(0)[[replicate_column]]
        estimated_df.index.name = 'isoform'
        estimated_df.columns = ['estimated_abund']
        df = estimated_df.join(true_expression_df,on='isoform',how='inner').reset_index()
        dfs.append(preprocess_single_sample_util(df, kvalues_dict,num_exon_dict,isoform_length_dict,isoform_gene_dict))
    anno_df = pd.DataFrame({'K_value':pd.Series(kvalues_dict),'num_exons':pd.Series(num_exon_dict),'isoform_length':pd.Series(isoform_length_dict),'gene':pd.Series(isoform_gene_dict)})
    anno_df.index.name = 'isoform'
    anno_df = anno_df.reset_index()
    return dfs,anno_df,method_names

def preprocess_multi_sample_multi_method(list_of_contents,ground_truth_given,is_long_read=True,K_value_selection='Condition_number'):
    estimated_dfs,method_names = load_zipped_data(list_of_contents[0])
    kvalues_dict,num_exon_dict,isoform_length_dict,isoform_gene_dict  = load_annotation(list_of_contents[1],is_long_read,K_value_selection)
    if (ground_truth_given):
        true_expression_df = load_data(list_of_contents[2]).set_index(0)
    else:
        true_expression_df = None
    dfs = []
    for estimated_df in estimated_dfs:
        estimated_df = estimated_df.set_index(0)
        if (ground_truth_given):
            intersected_index = true_expression_df.index.intersection(estimated_df.index)
            estimated_df = estimated_df.loc[intersected_index,:].reset_index()
            temp_true_expression_df = true_expression_df.loc[intersected_index,:].reset_index()
            dfs.append(preprocess_multi_sample_diff_condition_util(estimated_df,temp_true_expression_df, kvalues_dict,num_exon_dict,isoform_length_dict,isoform_gene_dict))
        else:
            estimated_df = estimated_df.reset_index()
            dfs.append(preprocess_multi_sample_diff_condition_util(estimated_df,None, kvalues_dict,num_exon_dict,isoform_length_dict,isoform_gene_dict))

    anno_df = pd.DataFrame({'K_value':pd.Series(kvalues_dict),'num_exons':pd.Series(num_exon_dict),'isoform_length':pd.Series(isoform_length_dict),'gene':pd.Series(isoform_gene_dict)})
    anno_df.index.name = 'isoform'
    anno_df = anno_df.reset_index()
    return dfs,anno_df,method_names
# @cache.memoize()
# def calculate_statistics_multi_sample_same_condition(list_of_contents):
#     df = preprocess_files_multi_sample_same_condition(list_of_contents)
#     estimated_df = load_data(list_of_contents[0])
#     if (list_of_contents[2] is not None):
#         true_expression_df = load_data(list_of_contents[2])
#         return get_multi_sample_same_conditon_metrics(estimated_df,true_expression_df,df)
#     return []

# @cache.memoize()
# def preprocess_files_multi_sample_same_condition(list_of_contents):
#     estimated_df = load_data(list_of_contents[0])
#     annotation = load_annotation(list_of_contents[1])
#     if (list_of_contents[2] is not None):
#         true_expression_df = load_data(list_of_contents[2])
#         df = preprocess_multi_sample_df_same_condition(estimated_df,true_expression_df,annotation)
    
#     return df
