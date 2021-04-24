import dash_core_components as dcc
import dash_html_components as html
import dash_table

table_layout = html.Div(children=[
    html.Div(children=[
        html.H2(children='Statistics'),
        dash_table.DataTable(
            id='statistics_table',
            style_cell={'fontSize':22, 'font-family':'Helvetica'},
            style_data={
            'height': 'auto',
            'width': 'auto'
        },
        )])
])