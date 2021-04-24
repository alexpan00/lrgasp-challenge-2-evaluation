import dash_core_components as dcc
import dash_html_components as html

from dash_extensions import Download
from static_data import filter_options, transform_options, data_sample_options

plot_layout = html.Div(children=[
    html.Div(children=[
        html.H2(children='Plots'),
        html.Div(children=[
            html.Div(children=[dcc.Dropdown(id='plot_figure_name')]),
            html.Div(id='transform_option_groups',children=[html.Label(["X axis transform", dcc.RadioItems(
            options=transform_options, id='x_axis_transform', value='log')], style={'display': 'none'}),
        html.Label(["Y axis transform", dcc.RadioItems(
            options=transform_options, id='y_axis_transform', value='log')],style={'display': 'none'})]),
            html.Label(["Scale", dcc.Dropdown(
                options=filter_options, id='scale', value='all')]),
            html.Button('Add plot', id='add_plot', n_clicks=0),
        ]),
        html.Div(id='plot_tabs_div',children=dcc.Tabs(id='plot_tabs', children=[])),
        dcc.Loading(
        children=[html.Button(
            "Download all figures", id="download-all-figures-btn"), Download(id='download-all-figures-action')],
        type="circle")
    ])
])
