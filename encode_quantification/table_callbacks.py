from dash.dependencies import Input, Output, State, ALL, MATCH
from dash.exceptions import PreventUpdate
from app import app
from preprocess import *
from preprocess_util import *
from static_data import *

@app.callback(
    [Output('statistics_table', 'data'),
     Output('statistics_table', 'columns'), ],
    Input('on_data_load', 'children'),
    [State('data_sample_option', 'value'),
    State('multi_methods_option','value'),
     State('ground_truth_given', 'value'),
     State({'type': "upload_file", 'index': ALL}, 'children'),
     State("replicate_column_selector", 'value')])
def update_table(on_data_load, data_sample_option, multi_methods_option,ground_truth_given_val, list_of_contents, replicate_column):
    if (on_data_load is None):
        raise PreventUpdate
    list_of_contents = [c for c in list_of_contents[0:2] if c is not None] + list_of_contents[2:]
    ground_truth_given = True if ground_truth_given_val == 'Yes' else False
    data, columns = [],[]
    if (multi_methods_option=='single_method'):
        if data_sample_option == 'single_sample':
            df,anno_df = preprocess_single_sample(list_of_contents,replicate_column)
            data.append(prepare_single_sample_table_metrics(df['true_abund'], df['estimated_abund'],df))
            columns = [m for m in single_sample_table_metrics if m['id'] in data[0]]
        elif data_sample_option == 'multi_sample_diff_condition':
            df,anno_df = preprocess_multi_sample_diff_condition(list_of_contents,ground_truth_given)
            data.append(prepare_multi_sample_diff_conditon_table_metrics(df,ground_truth_given))
            columns = [m for m in multi_sample_diff_condition_table_metrics if m['id'] in data[0]]
    elif (multi_methods_option=='multi_method'):
        if data_sample_option == 'single_sample':
            dfs,anno_df,method_names = preprocess_single_sample_multi_method(list_of_contents,replicate_column)
            for df,method_name in zip(dfs,method_names):
                row_dict = prepare_single_sample_table_metrics(df['true_abund'], df['estimated_abund'],df)
                row_dict['method'] = method_name
                data.append(row_dict)
            columns = [{'name': 'Method', 'id': 'method'}]+[m for m in single_sample_table_metrics if m['id'] in data[0]]
        elif data_sample_option == 'multi_sample_diff_condition':
            dfs,anno_df,method_names = preprocess_multi_sample_multi_method(list_of_contents,ground_truth_given)
            for df,method_name in zip(dfs,method_names):
                row_dict = prepare_multi_sample_diff_conditon_table_metrics(df,ground_truth_given)
                row_dict['method'] = method_name
                data.append(row_dict)
            columns = [{'name': 'Method', 'id': 'method'}]+[m for m in multi_sample_diff_condition_table_metrics if m['id'] in data[0]]
    # elif data_sample_option == 'multi_sample_same_condition':
    #     data = calculate_statistics_multi_sample_same_condition(list_of_contents)
    #     columns = [m for m in multi_sample_same_condition_table_metrics if m in data[0]]
    #     # data = [{'alfc': 0, 'f1': 0, 'auc': 0}]
    return data, columns
