import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd
import math
import pickle
from static_data import ARR_ranges, on_plot_shown_label,fig_size,color_schemes,themes
from preprocess_util import *
from plot_func.plot_util_multi_method import *
from plot_func.multi_method_plotter import Multi_method_plotter
class Single_sample_multi_method_plotter(Multi_method_plotter):
    def __init__(self,plot_dfs,anno_df,method_names,output_path):
        Multi_method_plotter.__init__(self,plot_dfs,anno_df,method_names,output_path)
    def plot_dist(self,x_axis_column_name, scale):
        return Multi_method_plotter.plot_dist(self,x_axis_column_name, scale)
    def plot_arr(self,x_axis_column_name,scale):
        fig = go.Figure()
        for plot_df,method_name in zip(self.plot_dfs,self.method_names):
            plot_df = filter_by_scale(scale, plot_df)
            plot_df, custom_sort = get_group_range(plot_df, x_axis_column_name,None,None)
            plot_df = plot_df.groupby(by='group_range').count().reset_index()
            plot_df['Frequency'] = plot_df['isoform']/plot_df['isoform'].sum()
            fill_na_df = pd.DataFrame({'group_range':['{:.0%}-{:.0%}'.format(i/10,(i+1)/10) for i in range(1,10)]+['<=10%','>100%'],'Frequency':[0 for i in range(11)]}).set_index('group_range')
            plot_df = plot_df.append(fill_na_df.loc[fill_na_df.index.difference(pd.Index(plot_df['group_range']))].reset_index())
            plot_df = plot_df.sort_values(
                by=['group_range'], key=lambda col: custom_sort(col, ARR_ranges))
            fig.add_trace(go.Bar(x=plot_df['group_range'],
                        y=plot_df['Frequency'],name=method_name))
        fig.update_layout(
            xaxis_title='Abundance Recovery Rate',
            yaxis_title='Frequency',
            width=fig_size['rec']['width'],height=fig_size['rec']['height'],template=themes['medium_single']
        )
        return fig
    def plot_resolution_entropy(self,scale):
        fig = go.Figure()
        for i,plot_df,method_name in zip(range(len(self.method_names)),self.plot_dfs,self.method_names):
            plot_df = filter_by_scale(scale, plot_df)
            RE = [get_resolution_entropy(plot_df['estimated_abund'],100)]
            fig.add_trace(go.Bar(x=['Resolution Entropy'],
                        y=RE,name=method_name,marker_color=color_schemes[i]))
        fig.update_layout(
            width=fig_size['small_rec']['width'],height=fig_size['small_rec']['height'],template=themes['medium_single']
        )
        return fig
    def plot_corr_scatter(self,x_axis_column_name, y_axis_column_name,scale):
        fig = make_subplots(rows=1, cols=len(self.plot_dfs), vertical_spacing=0.2, horizontal_spacing=0.1,subplot_titles=self.method_names)
        x_maxs,y_maxs = [],[]
        for plot_df,method_name,i in zip(self.plot_dfs,self.method_names,range(len(self.plot_dfs))):
            plot_df = filter_by_scale(scale, plot_df)
            plot_df = plot_df[(np.log2(plot_df[x_axis_column_name]+1) >=1) &  (np.log2(plot_df[y_axis_column_name]+1) >= 1)]
            x = plot_df[x_axis_column_name]
            y = plot_df[y_axis_column_name]
            x = np.log2(x + 1)
            y = np.log2(y + 1)
            x_maxs.append(x.max())
            y_maxs.append(y.max())
            x,y,density = get_density(x,y)
            fig.add_trace(go.Scattergl(x=x, y=y, mode='markers', name='Value',marker=dict(size=5,color=density,colorscale='viridis'),showlegend=False), col=i+1,row=1)
            fig.add_trace(go.Histogram2dContour(x=x, y=y, name='Density',contours={'coloring':'none','showlabels':True},showlegend=False),  col=i+1,row=1)
        x_title = 'Log2(True TPM+1)'
        y_title = 'Log2(Estimated TPM+1)'
        fig.update_layout(autosize=False,width=fig_size['square']['width']*len(self.plot_dfs),height=fig_size['square']['height'],template=themes['medium_multi'])
        fig.update_xaxes(title_text=x_title,range=[1,max(x_maxs+y_maxs)])
        fig.update_yaxes(title_text=y_title,range=[1, max(x_maxs+y_maxs)])
        return fig
    def plot_std_scatter(self,x_axis_column_name, y_axis_column_name,scale):
        fig = make_subplots(rows=1, cols=len(self.plot_dfs), vertical_spacing=0.2, horizontal_spacing=0.1,subplot_titles=self.method_names)
        x_maxs = []
        for plot_df,method_name,i in zip(self.plot_dfs,self.method_names,range(len(self.plot_dfs))):
            plot_df = filter_by_scale(scale, plot_df)
            plot_df = plot_df[(np.log2(plot_df[x_axis_column_name]+1) >=1)]
            x = plot_df[x_axis_column_name]
            y = plot_df[y_axis_column_name]
            x = np.log2(x + 1)
            x_maxs.append(x.max())
            x,y,density = get_density(x,y)
            fig.add_trace(go.Scattergl(x=x, y=y, mode='markers', name='Value',marker=dict(size=5,color=density,colorscale='viridis'),showlegend=False), col=i+1,row=1)
            fig.add_trace(go.Histogram2dContour(x=x, y=y, name='Density',contours={'coloring':'none','showlabels':True},showlegend=False),  col=i+1,row=1)
        x_title = 'Log2(Estimated TPM+1)'
        y_title = 'COV'
        fig.update_layout(autosize=False,width=fig_size['small_square']['width']*len(self.plot_dfs),height=fig_size['small_square']['height'],template=themes['large_single'])
        fig.update_xaxes(title_text=x_title,range=[1,max(x_maxs)])
        fig.update_yaxes(title_text=y_title)
        return fig
    def plot_grouped_curve(self,x_axis_column_name, y_axis_column_names,scale):
        figure_cols = math.ceil(math.sqrt(len(y_axis_column_names)))
        figure_rows = math.ceil(len(y_axis_column_names)/ figure_cols)
        # figure_rows = math.ceil(math.sqrt(len(y_axis_column_names)))
        # figure_cols = math.ceil(len(y_axis_column_names)/ figure_rows)
        fig = make_subplots(rows=figure_rows, cols=figure_cols, vertical_spacing=0.25, horizontal_spacing=0.1)
        ranges,max_threshold = prepare_ranges(self.plot_dfs[0],x_axis_column_name)
        # f = open('{}/plot.pkl'.format(self.output_path),'ab')
        for plot_df,method_name,j in zip(self.plot_dfs,self.method_names,range(len(self.plot_dfs))):
            plot_df = filter_by_scale(scale, plot_df)
            plot_df, custom_sort = get_group_range(plot_df, x_axis_column_name,ranges,max_threshold)
            # fig = go.Figure()
            # figure_rows = math.ceil(math.sqrt(len(y_axis_column_names)))
            # figure_cols = math.ceil(len(y_axis_column_names)/ figure_rows)
            for i in range(len(y_axis_column_names)):
                row_num = math.ceil((i+1)/figure_cols)
                col_num = i % figure_cols+1
                y_axis_column_name = y_axis_column_names[i]

                if ((y_axis_column_name in ['mrd']) & (x_axis_column_name=='K_value')):
                    group_series = plot_df.groupby(by='group_range').apply(lambda df: get_single_sample_metric(
                        y_axis_column_name, df['true_abund'], df['estimated_abund'],df)).to_frame().reset_index()
                else:
                    group_series = plot_df.groupby(by='group_range').apply(lambda df: get_single_sample_metric(
                        y_axis_column_name, df['true_abund'], df['estimated_abund'],df)).to_frame().reset_index()
                group_series = group_series.rename(columns={0: y_axis_column_name}).sort_values(
                    by=['group_range'], key=lambda col: custom_sort(col))
                if (y_axis_column_name in 'nrmse'):
                    group_series[y_axis_column_name] = np.log2(group_series[y_axis_column_name]+1)
                if (y_axis_column_name in ['nrmse','mrd','mean_arr','spearmanr','RE']):
                    fig.add_trace(go.Bar(x=group_series['group_range'], y=group_series[y_axis_column_name]
                    , name='{}'.format(method_name),marker_color=color_schemes[j],showlegend=False), row=row_num, col=col_num)
                else:
                    fig.add_trace(go.Scatter(x=group_series['group_range'], y=group_series[y_axis_column_name],
                                        mode='lines+markers', name='{}'.format(method_name),marker_color=color_schemes[j],showlegend=False), row=row_num, col=col_num)
                pickle.dump([method_name, x_axis_column_name,group_series],f)
                fig.update_xaxes(
                    title_text=on_plot_shown_label[x_axis_column_name],tickangle = 45, row=row_num, col=col_num)
                fig.update_yaxes(
                    title_text=on_plot_shown_label[y_axis_column_name], row=row_num, col=col_num)
                if (y_axis_column_name=='nrmse'):
                    fig.update_yaxes(title_text='Log2(NRMSE+1)', row=row_num, col=col_num)
        f.close()
        fig.update_traces(showlegend=True,col=1,row=1)
        fig.update_layout(
            autosize=False,
            width=fig_size['rec']['width']*figure_cols,height=fig_size['rec']['height']*figure_rows)
        return fig
    # def plot_corr_box_plot(self,x_axis_column_name,y_axis_column_name,scale):
    #     fig = make_subplots(rows=1, cols=len(self.plot_dfs),subplot_titles=self.method_names)
    #     shared_bins_cond = None
    #     for plot_df,method_name,j in zip(self.plot_dfs,self.method_names,range(len(self.plot_dfs))):
    #         plot_df = filter_by_scale(scale, plot_df)
    #         plot_df = plot_df[(np.log2(plot_df[x_axis_column_name]+1) >=1) &  (np.log2(plot_df[y_axis_column_name]+1) >= 1)]
    #         df,shared_bins_cond = prepare_corr_box_plot_data(np.log2(plot_df[y_axis_column_name]+1),np.log2(plot_df[x_axis_column_name]+1),shared_bins_cond)
    #         fig.add_trace(go.Box(x=df['true_abund'],y=df['estimated_abund'],name=method_name),col=j+1,row=1)
    #     fig.update_xaxes(title_text='Log2(True abundance+1)')
    #     fig.update_yaxes(title_text='Log2(Estimated abundance+1)')
    #     fig.update_layout(
    #         autosize=False,showlegend=True,width=fig_size['small_square']['width']*len(self.plot_dfs),height=fig_size['small_square']['height'],template=themes['small_multi'])
    #     return fig
    def plot_corr_box_plot(self,x_axis_column_name,y_axis_column_name,scale):
        fig = make_subplots(rows=1, cols=1)
        shared_bins_cond = None
        for plot_df,method_name,j in zip(self.plot_dfs,self.method_names,range(len(self.plot_dfs))):
            plot_df = filter_by_scale(scale, plot_df)
            # plot_df = plot_df[(np.log2(plot_df[x_axis_column_name]+1) >=1) &  (np.log2(plot_df[y_axis_column_name]+1) >= 1)]
            df,shared_bins_cond = prepare_corr_box_plot_data(np.log2(plot_df[y_axis_column_name]+1),np.log2(plot_df[x_axis_column_name]+1),shared_bins_cond)
            fig.add_trace(go.Box(x=df['true_abund'],y=df['estimated_abund'],name=method_name),col=1,row=1)
        fig.update_xaxes(title_text='Log2(True TPM+1)')
        fig.update_yaxes(title_text='Log2(Estimated TPM+1)')
        fig.update_layout(boxmode = "group",
            autosize=False,showlegend=True,width=fig_size['square']['width'],height=fig_size['square']['height'],template=themes['large_single'])
        return fig
    def plot(self,plot_figure_name, scale):
        x_axis_column_name = single_sample_plot_figures[plot_figure_name]['x']
        y_axis_column_name = single_sample_plot_figures[plot_figure_name]['y']
        if y_axis_column_name == 'dist':
            fig = self.plot_dist(x_axis_column_name, scale)
        elif plot_figure_name == 'Histogram of Abundance Recovery Rate':
            fig = self.plot_arr(x_axis_column_name, scale)
        elif plot_figure_name in ["Statistics with different K values",'Statistics with different isoform lengths','Statistics with different numbers of exons','Statistics with different numbers of isoforms','Statistics with different expression level']:
            fig = self.plot_grouped_curve(x_axis_column_name,y_axis_column_name,scale)
        elif plot_figure_name in ['Correlation of estimated abundance and ground truth']:
            fig = self.plot_corr_scatter(x_axis_column_name, y_axis_column_name, scale)
        elif plot_figure_name in ['coefficient of variation vs estimated abundance scatter']:
            fig = self.plot_std_scatter(x_axis_column_name, y_axis_column_name, scale)
        elif plot_figure_name == 'Correlation Boxplot of estimated abundance and ground truth':
            fig = self.plot_corr_box_plot(x_axis_column_name,y_axis_column_name,scale)
        elif plot_figure_name == 'Resolution Entropy':
            fig = self.plot_resolution_entropy(scale)
        try:
            fig.update_layout(title=plot_figure_name, title_x=0.5)
            fig.update_xaxes(exponentformat='e',automargin=True)
            fig.update_yaxes(exponentformat='e',automargin=True)
        except:
            print(plot_figure_name)
        return fig