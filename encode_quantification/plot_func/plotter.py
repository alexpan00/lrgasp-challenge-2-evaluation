import plotly.express as px
import plotly.graph_objects as go

from plot_func.plot_util import *
import static_data
class Plotter:
    def __init__(self,plot_df,anno_df):
        self.plot_df = plot_df
        self.anno_df = anno_df
    def plot_dist(self,x_axis_column_name, scale):
        anno_df = self.anno_df
        pd.Series({'num_isoforms':anno_df.shape[0],'num_genes':anno_df['gene'].unique().shape[0]}).to_csv(f'{static_data.output_dir}/num_genes_isoforms.tsv',sep='\t')
        anno_df = filter_by_scale(scale, self.anno_df)
        if (x_axis_column_name in ['K_value','num_isoforms']):
            anno_df, custom_sort = get_k_val_dist(anno_df, x_axis_column_name)
        else:
            anno_df, custom_sort = get_group_range(anno_df, x_axis_column_name)
        mean_x_axis_column = anno_df[['gene',x_axis_column_name]].groupby(by='gene').mean().mean()
        out_df = pd.DataFrame({x_axis_column_name:mean_x_axis_column})
        # out_df.to_csv(f'{static_data.output_dir}/{x_axis_column_name}_mean.tsv',sep='\t')
        anno_df = anno_df.groupby(by='group_range').count().reset_index()
        if (x_axis_column_name in ['num_isoforms']):
            anno_df = anno_df.rename(columns={'gene': 'Count'}).sort_values(
                    by=['group_range'], key=lambda col: custom_sort(col))
        else:
            anno_df = anno_df.rename(columns={'isoform': 'Count'}).sort_values(
                    by=['group_range'], key=lambda col: custom_sort(col))
        fig = px.bar(anno_df, x='group_range', y='Count')
        out_df = pd.DataFrame({'group_range':anno_df['group_range'],'Count':anno_df['Count']})
        # out_df.to_csv(f'{static_data.output_dir}/{x_axis_column_name}_distribution.tsv',sep='\t')
        # fig = px.histogram(plot_df, x=x_axis_column_name, log_x=True if x_axis_transform ==
        #                    'log' else False, log_y=True if y_axis_transform == 'log' else False,width=1000,height=1000)
        fig.update_layout(
            xaxis_title=on_plot_shown_label[x_axis_column_name],
            yaxis_title="Count",
            width=fig_size['rec']['width'],
            height=fig_size['rec']['height'],
            template= themes['large_single']
        )
        return fig