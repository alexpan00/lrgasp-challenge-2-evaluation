import plotly.graph_objects as go
from static_data import *
from preprocess import *
table_comment = {
    'consistency':'Consistency Measure (CM) is calculated for C threshold  = 1',
    # 'consistency':'Resolution Entropy (RE) is calculated by the',
    'reproducibility':'Reproducibility Measure (RM) is calculated by the L2 norm of the standard deviation',
}
def generate_table(args,output_path,is_multi_sample,is_multi_method,is_long_read,K_value_selection,sections):
    data, columns = [],[]
    if (is_multi_method == False):
        df,anno_df = args
        if (is_multi_sample == False):
            data.append(prepare_single_sample_table_metrics(df['true_abund'], df['estimated_abund'],df))
            columns = [m for m in single_sample_table_metrics if m['id'] in data[0]]
        else:
            data.append(prepare_multi_sample_diff_conditon_table_metrics(df,True))
            columns = [m for m in multi_sample_diff_condition_table_metrics if m['id'] in data[0]]
    else:
        dfs,anno_df,method_names = args
        if (is_multi_sample == False):
            for df,method_name in zip(dfs,method_names):
                row_dict = prepare_single_sample_table_metrics(df['true_abund'], df['estimated_abund'],df)
                row_dict['method'] = method_name
                data.append(row_dict)
            columns = [{'name': 'Method', 'id': 'method'}]+[m for m in single_sample_table_metrics if m['id'] in data[0]]
        else:
            for df,method_name in zip(dfs,method_names):
                row_dict = prepare_multi_sample_diff_conditon_table_metrics(df,True)
                row_dict['method'] = method_name
                data.append(row_dict)
            columns = [{'name': 'Method', 'id': 'method'}]+[m for m in multi_sample_diff_condition_table_metrics if m['id'] in data[0]]
    
    for i in range(len(sections)):
        shown_columns = []
        if sections[i]['id'] == 'estimation_error':
            shown_columns = ['Normalized Root Mean Square Error', 'Median Relative Difference', "Spearman's rho",'Mean Abundance Recovery Rate']
        elif sections[i]['id'] == 'resolution_entropy':
            shown_columns = ['Resolution Entropy']
        elif sections[i]['id'] == 'consistency':
            shown_columns = ['Consistency Measure']
        elif sections[i]['id'] == 'reproducibility':
            shown_columns = ['Reproducibility Measure']
        elif sections[i]['id'] == 'fold_change':
            shown_columns = ['Precision', 'Recall', 'Accuracy', 'F1 Score', 'AUC']
        if len(shown_columns) > 0:
            if (is_multi_method==True):
                shown_columns.insert(0,'Method')
            selected_columns = [c for c in columns if c['name'] in shown_columns]
            columns_values = [c['name'] for c in selected_columns]
            columns_format = ['{:'+c['format']['specifier']+'}' if 'format' in c else '{}' for c in selected_columns]
            cells_data = [[columns_format[j].format(row[selected_columns[j]['id']])  for j in range(len(selected_columns))] for row in data]
            
            sections[i]['table'] = {'header':columns_values,'rows':cells_data}
            if (sections[i]['id'] in table_comment):
                sections[i]['table']['comment'] = table_comment[sections[i]['id']]
    return sections