import uuid

import dash_uploader as du
import dash_core_components as dcc
import dash_html_components as html
from static_data import filter_options, data_sample_options,multi_methods_options
upload_style = {
    'width': '100%',
    'height': '40px',
    'lineHeight': '60px',
    'borderWidth': '1px',
    'borderStyle': 'dashed',
    'borderRadius': '5px',
    'textAlign': 'center',
    'margin': '10px',
}
input_layout = html.Div(children=[
    html.H2(children='Data'),
    html.Label(["Data types",
            dcc.RadioItems(
                options=data_sample_options, id='data_sample_option', value='single_sample', labelStyle={'display': 'inline-block'})]),
    html.Label(["Multi methods",
        dcc.RadioItems(
            options=multi_methods_options, id='multi_methods_option', value='single_method', labelStyle={'display': 'inline-block'})]),
    html.Div(["Quantification result",
                html.Div(id='single_method_quantif_uploader',children=[du.Upload(
                    id={'type': "du_uploader", 'index': 1},
                    text = 'Drag and Drop or Click here to select a file',
                    max_files = 1,
                    max_file_size=4096,
                    filetypes=['tsv'],
                    upload_id=uuid.uuid1(),
                ),html.Div(id={'type': "upload_file", 'index':1},style={'display':'none'}),html.Div(id={'type': "file_preview", 'index':1})]),
                html.Div(id='multi_method_quantif_uploader',children=[du.Upload(
                    id={'type': "du_uploader", 'index': 2},
                    text = 'Drag and Drop or Click here to select a zip file',
                    max_files = 1,
                    max_file_size=4096,
                    filetypes=['zip'],
                    upload_id=uuid.uuid1(),
                ),html.Div(id={'type': "upload_file", 'index':2},style={'display':'none'}),html.Div(id={'type': "file_preview", 'index':2})]),
                ]),
                # dcc.Upload(
                #     id={'type': "upload_file", 'index': 1},
                #     children=html.Div([
                #         'Drag and Drop or ',
                #         html.A('Select a file')
                #     ]),
                #     style=upload_style,
                # ),
    html.Div(id='replicate_column_selector_div',children=[dcc.Input(id='replicate_column_selector',style={'display':'none'})]),
    html.Div(["Annotation",
                du.Upload(
                    id={'type': "du_uploader", 'index': 4},
                     text = 'Drag and Drop or Click here to select a file',
                    max_file_size=1800,
                    filetypes=['gtf'],
                    upload_id=uuid.uuid1(),
                ),
                html.Div(id={'type': "upload_file", 'index':4},style={'display':'none'}),
                # dcc.Upload(
                #     id={'type': "upload_file", 'index': 2},
                #     children=html.Div([
                #         'Drag and Drop or ',
                #         html.A('Select a file')
                #     ]),
                #     style=upload_style,
                # ),
                html.Div(id={'type': "file_preview", 'index': 4})]),

    html.Label(["Ground truth is given",
                dcc.RadioItems(
                    options=[{'label': 'Yes', 'value': 'Yes'},{'label': 'No', 'value': 'No'}], id='ground_truth_given', value='Yes', labelStyle={'display': 'inline-block'})]),
    html.Div(id='upload_ground_truth_div'),
    # html.Label(["Multiple methods results are given",
    #             dcc.RadioItems(
    #                 options=[{'label': i, 'value': i} for i in ['Yes', 'No']], id='multi_methods_given', value='No', labelStyle={'display': 'inline-block'})]),
    html.Button('Submit', id='confirm_data', n_clicks=0),
    html.Div(id='error_prompt',style={'color': 'red','font-weight': 'bold'}),
    html.Div(id='on_data_load',style={'display':'none'})])
