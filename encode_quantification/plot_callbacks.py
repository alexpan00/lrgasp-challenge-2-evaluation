import io
import zipfile
import base64
import dash_html_components as html
import dash_core_components as dcc
from plotly.io import write_image
import plotly.io as pio
import plotly.graph_objects as go

from dash.dependencies import Input, Output, State, ALL, MATCH
from dash.exceptions import PreventUpdate
from dash_extensions import Download
from dash_extensions.snippets import send_bytes

from app import app
from preprocess import *
from static_data import *
from plot_func.single_sample_plotter import Single_sample_plotter
from plot_func.multi_sample_plotter import Multi_sample_plotter
from plot_func.single_sample_multi_method_plotter import Single_sample_multi_method_plotter
from plot_func.multi_sample_multi_method_plotter import Multi_sample_multi_method_plotter
from plot_func.plot_util import define_theme,define_write_to_file_theme

define_theme()
pio.kaleido.scope.default_width = None
pio.kaleido.scope.default_height = None
@app.callback(Output('plot_figure_name', 'options'),
              Input('confirm_data', 'n_clicks'),
              [State('data_sample_option', 'value'),
               State('ground_truth_given', 'value') ])
def update_options(n_clicks, data_sample_option, ground_truth_given):
    if (n_clicks == 0):
        raise PreventUpdate
    options = []
    if data_sample_option == 'single_sample':
        options = [{'label': i, 'value': i} for i in single_sample_plot_figures]
    elif ((data_sample_option == 'multi_sample_diff_condition') & (ground_truth_given=='Yes')):
        options = [{'label': i, 'value': i}
            for i in multi_sample_diff_condition_with_ground_truth_plot_figures]
    elif ((data_sample_option == 'multi_sample_diff_condition') & (ground_truth_given=='No')):
        options = [{'label': i, 'value': i}
            for i in multi_sample_diff_condition_without_ground_truth_plot_figures]
    return options


@app.callback(
    Output('plot_tabs', 'children'),
    Input('add_plot', 'n_clicks'),
    [State('on_data_load', 'children'),
     State('ground_truth_given','value'),
     State('plot_figure_name', 'value'),
     State('scale', 'value'),
     State('plot_tabs', 'children'),
     State({'type': "upload_file", 'index': ALL}, 'children'),
     State('annotation_option','value'),
     State('data_sample_option', 'value'),
     State('multi_methods_option','value'),
     State("replicate_column_selector", 'value')])
def update_graph(n_clicks, on_data_load, ground_truth_given_val,plot_figure_name, scale,plot_tabs_children, list_of_contents,annotation_option, data_sample_option,multi_methods_option, replicate_column):
    if ((n_clicks == 0) | (on_data_load is None)):
        raise PreventUpdate
    list_of_contents.insert(2, [annotation_option])
    list_of_contents = [c for c in list_of_contents[0:2] if c is not None] + list_of_contents[2:]
    ground_truth_given = True if ground_truth_given_val == 'Yes' else False
    if (multi_methods_option=='single_method'):
        if (data_sample_option == 'single_sample'):
            plot_df,anno_df = preprocess_single_sample(
                list_of_contents, replicate_column)
            plotter = Single_sample_plotter(plot_df,anno_df)
            fig = plotter.plot(plot_figure_name,scale)
        elif (data_sample_option == 'multi_sample_diff_condition'):
            plot_df,anno_df = preprocess_multi_sample_diff_condition(list_of_contents,ground_truth_given)
            plotter = Multi_sample_plotter(plot_df,anno_df)
            fig = plotter.plot(plot_figure_name, scale,ground_truth_given)
    elif (multi_methods_option=='multi_method'):
        if (data_sample_option == 'single_sample'):
            plot_dfs,anno_df,method_names = preprocess_single_sample_multi_method(
                list_of_contents, replicate_column)
            plotter = Single_sample_multi_method_plotter(plot_dfs,anno_df,method_names)
            fig = plotter.plot(plot_figure_name,scale)
        elif (data_sample_option == 'multi_sample_diff_condition'):
            plot_dfs,anno_df,method_names = preprocess_multi_sample_multi_method(
                list_of_contents, ground_truth_given)
            plotter = Multi_sample_multi_method_plotter(plot_dfs,anno_df,method_names)
            fig = plotter.plot(plot_figure_name,scale,ground_truth_given)
    fig.update_xaxes(exponentformat='e',automargin=True)
    fig.update_yaxes(exponentformat='e',automargin=True)
    if (scale != 'all'):
        shown_plot_figure_name = '{} for {}'.format(plot_figure_name, scale)
    else:
        shown_plot_figure_name = plot_figure_name
    fig.update_layout(title=shown_plot_figure_name, title_x=0.5)
    tab = dcc.Tab(label=shown_plot_figure_name, id={'type': 'plot_tab', 'index': n_clicks}, children=[dcc.Graph(id={'type': 'plot_graph', 'index': n_clicks}, figure=fig),
                                                                                                                           dcc.Loading(
        id={'type': "figure-creating-loading", 'index': n_clicks},
        children=[html.Button(
            "Download figure", id={'type': "download-figure-btn", 'index': n_clicks}), Download(id={'type': "download-figure-action", 'index': n_clicks})],
        type="circle",
    )])
    plot_tabs_children.insert(0, tab)
    return plot_tabs_children


@app.callback(Output({'type': 'download-figure-action', 'index': MATCH}, 'data'),
              Input({'type': "download-figure-btn", 'index': MATCH}, 'n_clicks'),
               [State({'type': 'plot_graph', 'index': MATCH}, 'figure'),
               State({'type': 'plot_tab', 'index': MATCH}, 'label')])
def download_figure(n_clicks, figure, figure_title):
    if (n_clicks is None):
        raise PreventUpdate
    with io.BytesIO() as zip_buffer:
        with zipfile.ZipFile(zip_buffer, mode="a") as zf:
            fig = go.Figure(figure)
            for fmt in ['png','pdf']:
                buf = fig.to_image(format=fmt)
                zf.writestr("{}.{}".format(figure_title,fmt), buf)
        content = base64.b64encode(zip_buffer.getvalue()).decode()
    return dict(content=content, filename='{}.zip'.format(figure_title),  mime_type='application/zip',base64=True)

@app.callback(Output('download-all-figures-action', 'data'),
              Input('download-all-figures-btn', 'n_clicks'),
               [State({'type': 'plot_graph', 'index': ALL}, 'figure'),
               State({'type': 'plot_tab', 'index': ALL}, 'label')])
def download_all_figures(n_clicks, figures, figure_titles):
    if (n_clicks is None):
        raise PreventUpdate
    with io.BytesIO() as zip_buffer:
        with zipfile.ZipFile(zip_buffer, mode="a") as zf:
            for figure,figure_title in zip(figures,figure_titles):
                fig = go.Figure(figure)
                for fmt in ['png','pdf']:
                    buf = fig.to_image(format=fmt)
                    zf.writestr("{}.{}".format(figure_title,fmt), buf)
        content = base64.b64encode(zip_buffer.getvalue()).decode()
    return dict(content=content, filename='figures.zip', mime_type='application/zip',base64=True)

# def update_multi_sample_same_condition_graph(plot_df, plot_figure_name, x_axis_column_name, y_axis_column_name, x_axis_transform, y_axis_transform, scale,smoothing):
#     if plot_figure_name == 'coefficient of variation vs estimated abundance scatter':
#         fig = plot_scatter(x_axis_column_name, y_axis_column_name,
#                         x_axis_transform, y_axis_transform, scale, plot_figure_name, plot_df)
#         return fig
