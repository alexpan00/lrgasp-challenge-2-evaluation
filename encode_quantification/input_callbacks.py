import uuid

import dash_uploader as du
from dash.dependencies import Input, Output, State, ALL, MATCH
from dash.exceptions import PreventUpdate
import dash_html_components as html
import dash_core_components as dcc
import dash_table

from app import app
from input_layouts import upload_style
from preprocess import *

@du.callback(
    output = Output({'type': "upload_file", 'index': MATCH}, 'children'),
    id = {'type': "du_uploader", 'index': MATCH},
)
def on_file_uploaded(filename):
    if filename is None:
        raise PreventUpdate
    return filename


@app.callback(Output('upload_ground_truth_div', 'children'),
              Input('ground_truth_given', 'value'))
def show_ground_truth_upload(ground_truth_given):
    if (ground_truth_given == 'No'):
        return
    return html.Div(["True expression data",
                        du.Upload(
                            id={'type': "du_uploader", 'index': 3},
                             text = 'Drag and Drop or Click here to select a file',
                            max_file_size=1800,
                            filetypes=['tsv'],
                            upload_id=uuid.uuid1(),
                        ),
                        html.Div(id={'type': "upload_file", 'index': 3},style={'display':'none'}),
                        # dcc.Upload(
                        #     id={'type': "upload_file", 'index': 3},
                        #     children=html.Div([
                        #         'Drag and Drop or ',
                        #         html.A('Select a file')
                        #     ]),
                        #     style=upload_style,
                        # ), 
                        html.Div(id={'type': "file_preview", 'index':3})])
@app.callback(Output('single_method_quantif_uploader','style'),
                Output('multi_method_quantif_uploader','style'),
              Input('multi_methods_option', 'value'))
def show_multi_method_upload(multi_methods_option):
    if (multi_methods_option == 'multi_method'):
        return {'display':'none'}, {'display':'block'}
    else:
        return {'display':'block'}, {'display':'none'}

@app.callback(Output('replicate_column_selector_div', 'children'),
              Input({'type': "upload_file", 'index': ALL}, 'children'),
              Input('data_sample_option', 'value'),
              Input('multi_methods_option','value'))
def show_replicate_column_selector(list_of_contents,data_sample_option,multi_methods_option):
    if (data_sample_option == 'multi_sample_diff_condition'):
        return [dcc.Input(id='replicate_column_selector',style={'display':'none'})]
    if ((list_of_contents[0] is None) & (list_of_contents[1] is None)):
        raise PreventUpdate
    list_of_contents = [c for c in list_of_contents[0:2] if c is not None] + list_of_contents[2:]
    if (multi_methods_option=='multi_method'):
        dfs,_ = load_zipped_data(list_of_contents[0])
        return html.Div(["Select which replicate to analyze (between 1 and {}): ".format(dfs[0].shape[1]-1),dcc.Input(
            id="replicate_column_selector",
            type='number',
            value = 1,
            min = 1,
            max = dfs[0].shape[1]-1,
        )])
    elif (multi_methods_option=='single_method'):
        df = load_data(list_of_contents[0])
        return html.Div(["Select which replicate to analyze (between 1 and {}): ".format(df.shape[1]-1),dcc.Input(
            id="replicate_column_selector",
            type='number',
            value = 1,
            min = 1,
            max = df.shape[1] - 1,
        )])


def parse_table(contents):
    try:
        filename = contents[0].split('/')[-1]
        if 'tsv' in filename:
            df = load_data(contents)
            return html.Div('File name: {}, Number of lines: {}, Number of replicates: {}'.format(filename,df.shape[0],df.shape[1]-1))
        elif 'gtf' in filename:
            df = load_data(contents)
            return html.Div('File name: {}, Number of lines in file: {}'.format(filename,df.shape[0]))
        elif 'zip' in filename:
            dfs,_ = load_zipped_data(contents)
            return html.Div('Zip file name: {}, Number of files in zip file: {}'.format(filename,len(dfs)))
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])


@app.callback(Output({'type': "file_preview", 'index': MATCH}, 'children'),
              Input({'type': "upload_file", 'index': MATCH}, 'children'))
def on_upload_file(contents):
    if contents is not None:
        children = parse_table(contents)
        return children

@app.callback([Output('error_prompt', 'children'),Output('on_data_load', 'children')],Input('confirm_data','n_clicks'),
[State({'type': "upload_file", 'index': ALL}, 'children'),
State('annotation_option','value'),
State('data_sample_option','value'),
State('ground_truth_given','value'),
State('multi_methods_option','value'),
State("replicate_column_selector",'value')])
def on_data_load(n_clicks, list_of_contents,annotation_option,data_sample_option,ground_truth_given_val,multi_methods_option,replicate_column):
    if (n_clicks == 0):
        raise PreventUpdate
    list_of_contents.insert(2, [annotation_option])
    list_of_contents = [c for c in list_of_contents[0:2] if c is not None] + list_of_contents[2:]
    if (list_of_contents[0] is None):
        return 'No quantification result is uploaded!',None
    if (annotation_option is None):
        return 'No annotation is selected!',None
    if (((ground_truth_given_val == 'No') or (list_of_contents[2] is None)) and (data_sample_option == 'single_sample')):
        return 'Single sample data selected and no true expression data is uploaded. No available analysis!',None
    if ((ground_truth_given_val == 'Yes') and (list_of_contents[2] is None)):
        return 'True expression data is not uploaded.',None
    ground_truth_given = True if ground_truth_given_val == 'Yes' else False
    if (multi_methods_option == 'single_method'):
        if (data_sample_option == 'single_sample'):
            preprocessed_df,anno_df = preprocess_single_sample(list_of_contents,replicate_column)
        else:
            preprocessed_df,anno_df = preprocess_multi_sample_diff_condition(list_of_contents,ground_truth_given)
    elif (multi_methods_option == 'multi_method'):
        if (data_sample_option == 'single_sample'):
            preprocessed_df,anno_df,_ = preprocess_single_sample_multi_method(list_of_contents,replicate_column)
        else:
            preprocessed_df,anno_df,_ = preprocess_multi_sample_multi_method(list_of_contents,ground_truth_given)
    return None,''
@app.callback(Output('plot_tabs_div', 'children'),Input('confirm_data','n_clicks'))
def clear_tabs(n_clicks):
    return dcc.Tabs(id='plot_tabs', children=[])
