import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import flask
import os

from app import app
from plot_layouts import plot_layout
from table_layouts import table_layout
from input_layouts import input_layout
import plot_callbacks,table_callbacks,input_callbacks

app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    html.H1(children='Quantification'),
    html.Div(id='page-content'),
])
server = app.server


@app.callback(Output('page-content', 'children'),
              Input('url', 'pathname'))
def routing(pathname):
    # if pathname == '/':
    #     return upload_layout
    if pathname == '/':
        return [input_layout,table_layout,plot_layout]
    else:
        return '404'
# @app.server.route('/assets/plots/<path:path>')
# def serve_static(path):
#     root_dir = os.getcwd()
#     return flask.send_from_directory(
#         os.path.join(root_dir, 'assets/plots'), path
#     )


if __name__ == '__main__':
    app.run_server(debug=True)
