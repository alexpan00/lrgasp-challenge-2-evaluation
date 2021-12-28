from jinja2 import Environment, FileSystemLoader,select_autoescape
import os
import shutil
import argparse
from preprocess import *
from report_generator.plot import make_plots
from report_generator.table import generate_table
from plot_func.multi_method_plotter import Multi_method_plotter
import pickle
import pandas as pd
all_sections = [
    {'name':'Gene features','id':'gene_features','plots':[]},
    {'name':'Methods Legend','id':'method_legend','legend':[]},
    {'name':'Estimation error','id':'estimation_error','plots':[],'table':{}},
    {'name':'Resolution Entropy','id':'resolution_entropy','plots':[],'table':{}},
    {'name':'Consistency','id':'consistency','plots':[]},
    {'name':'Reproducibility','id':'reproducibility','plots':[]},
    {'name':'Fold change based evaluation','id':'fold_change','plots':[],'table':{}},
    {'name':'Split Statistics','id':'statistics','plots':[]},
]
def generate_output(output,output_path):
    des_assets_path = os.path.join(output_path,'assets') 
    src_assets_dir = os.path.join(os.path.dirname(__file__),'report_generator/templates/assets')
    for src_dir, dirs, files in os.walk(src_assets_dir):
        dst_dir = src_dir.replace(src_assets_dir, des_assets_path, 1)
        if not os.path.exists(dst_dir):
            os.makedirs(dst_dir)
        for file_ in files:
            src_file = os.path.join(src_dir, file_)
            dst_file = os.path.join(dst_dir, file_)
            if os.path.exists(dst_file):
                os.remove(dst_file)
            shutil.copy(src_file, dst_dir)
    with open(os.path.join(output_path,'Report.html'),'w') as f:
        f.write(output)
    results = []
    with open('{}/plot.pkl'.format(output_path),'rb') as f:
        while 1:
            try:
                res = pickle.load(f)
                results.append(res)
            except:
                break
    with open('{}/plot.txt'.format(output_path),'w') as f:
        for res in results:
            names = ['','','data','mean','error']
            for item,name in zip(res,names):
                if type(item) == pd.DataFrame:
                    f.write('{}\n'.format(name))
                    f.write(item.to_csv())
                    f.write('\n')
                else:
                    f.write('{}\n'.format(item))
            f.write('\n\n')
    os.remove('{}/plot.pkl'.format(output_path))
def preprocess_file(quantif_res_path,annotation_path,truth_path,is_multi_sample,is_multi_method,is_long_read,ground_truth_given,K_value_selection):
    input_paths = [[quantif_res_path],[annotation_path],[truth_path]]
    if (is_multi_method == False):
        if (is_multi_sample == False):
            df,anno_df = preprocess_single_sample(input_paths,1,is_long_read,K_value_selection)
        else:
            df,anno_df = preprocess_multi_sample_diff_condition(input_paths,ground_truth_given,is_long_read,K_value_selection)
        return df,anno_df
    else:
        if (is_multi_sample == False):
            dfs,anno_df,method_names = preprocess_single_sample_multi_method(input_paths,1,is_long_read,K_value_selection)
        else:
            dfs,anno_df,method_names = preprocess_multi_sample_multi_method(input_paths,ground_truth_given,is_long_read,K_value_selection)
        return dfs,anno_df,method_names
    
def render(quantif_res_path,annotation_path,truth_path,output_path,is_multi_sample,is_multi_method,is_long_read,ground_truth_given,K_value_selection):
    env = Environment(loader=FileSystemLoader(os.path.join(os.path.dirname(__file__),'report_generator/templates')),autoescape=select_autoescape(['html']))
    template = env.get_template('base.html')
    sections_num = [0]
    if (is_multi_method):
        sections_num = sections_num + [1]
    if (ground_truth_given):
        sections_num  = sections_num + [2]
    if (is_multi_sample):
        sections_num = sections_num + [3,4,5,7]
        if (ground_truth_given):
            sections_num  = sections_num + [6]
    else:
        sections_num = sections_num + [3,7]
    sections_num.sort()
    sections = [all_sections[i] for i in sections_num]
    if (is_multi_method == False):
        df,anno_df = preprocess_file(quantif_res_path,annotation_path,truth_path,is_multi_sample,is_multi_method,is_long_read,ground_truth_given,K_value_selection)
        args = (df,anno_df)
        sections = make_plots(args,output_path,is_multi_sample,is_multi_method,is_long_read,ground_truth_given,K_value_selection,sections)
        sections = generate_table(args,output_path,is_multi_sample,is_multi_method,is_long_read,ground_truth_given,K_value_selection,sections)

    else:
        dfs,anno_df,method_names = preprocess_file(quantif_res_path,annotation_path,truth_path,is_multi_sample,is_multi_method,is_long_read,ground_truth_given,K_value_selection) 
        with open('{}/df.pkl'.format(output_path),'wb') as f:
            pickle.dump(dfs,f)
        args = (dfs,anno_df,method_names)
        sections = make_plots(args,output_path,is_multi_sample,is_multi_method,is_long_read,ground_truth_given,K_value_selection,sections)
        sections = generate_table(args,output_path,is_multi_sample,is_multi_method,is_long_read,ground_truth_given,K_value_selection,sections)
    output = template.render(sections=sections)
    generate_output(output,output_path)
def parse_arguments():
    """
    Parse the arguments
    """
    parser = argparse.ArgumentParser(description="Quantification evaluation reporter",add_help=True)
    
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-a','--annotation', type=str, help="The path of annotation file [GTF]",required=True)
    requiredNamed.add_argument('-r','--result', type=str, help="The path of quantification result file [TSV\ZIP]",required=True)
    requiredNamed.add_argument('-o','--output', type=str, help="The path of output directory",required=True)
     
    requiredNamed.add_argument('--num_method',  type=str,help="Whether multi method data given ['Single' or 'Multi']")
    requiredNamed.add_argument('--num_samples',   type=str, help="Whether multi sample data given ['Single' or 'Multi']")
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument('-t','--truth', type=str, help="The path of true expression file [TSV]")
    optional.add_argument('--seq',  type=str,help="Whether long read data given ['LongRead' or 'ShortRead'] [default:ShortRead]",default='ShortRead')
    optional.add_argument('--K_value_selection',  type=str,help="Which K value calculation['Condition_number','K_value','Generalized_condition_number'] [default:Generalized_condition_number]",default='Generalized_condition_number')
    # optional = parser.add_argument_group('optional arguments')
    # optional.add_argument('--num_iterations',type=int,default=100, help="Number of iterations for EM algorithm [default:100]")
    
    args = parser.parse_args()
    ground_truth_given = args.truth is not None
    if args.num_method == 'Multi':
        is_multi_method = True
        if not args.result.endswith('.zip'):
            raise Exception('Invalid file format given')
    else:
        is_multi_method = False
        if not args.result.endswith('.tsv'):
            raise Exception('Invalid file format given')
    is_multi_sample = True if args.num_samples == 'Multi' else False
    is_long_read = True if args.seq == 'LongRead' else False
    if (not ground_truth_given) and (not is_multi_sample):
        raise Exception('No evaluation can be done for single sample data without ground truth!')
    render(args.result,args.annotation,args.truth,args.output,is_multi_sample,is_multi_method,is_long_read,ground_truth_given,args.K_value_selection)
if __name__ == "__main__":
    parse_arguments()
