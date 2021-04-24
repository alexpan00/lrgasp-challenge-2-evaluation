import plotly.express as px
import plotly.graph_objects as go

from plot_func.plot_util import *

class Multi_method_plotter:
    def __init__(self,plot_dfs,anno_df,method_names):
        self.plot_dfs = plot_dfs
        self.anno_df = anno_df
        self.method_names = method_names
    def plot_dist(self,x_axis_column_name, scale):
        anno_df = filter_by_scale(scale, self.anno_df)
        if (x_axis_column_name == 'K_value'):
            anno_df, custom_sort = get_k_val_dist(anno_df, x_axis_column_name)
        else:
            anno_df, custom_sort = get_group_range(anno_df, x_axis_column_name)
        anno_df = anno_df.groupby(by='group_range').count().reset_index()
        anno_df = anno_df.rename(columns={'isoform': 'Count'}).sort_values(
                by=['group_range'], key=lambda col: custom_sort(col))
        fig = px.bar(anno_df, x='group_range', y='Count')
        # fig = px.histogram(plot_df, x=x_axis_column_name, log_x=True if x_axis_transform ==
        #                    'log' else False, log_y=True if y_axis_transform == 'log' else False,width=1000,height=1000)
        fig.update_layout(
            xaxis_title=on_plot_shown_label[x_axis_column_name],
            yaxis_title="Count",
            width=fig_size['rec']['width'],
            height=fig_size['rec']['height'],
            template= themes['small_single']
        )
        return fig