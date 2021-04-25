import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd
import math

from static_data import ARR_ranges, on_plot_shown_label,fig_size,color_schemes,themes
from preprocess_util import *
from plot_func.plot_util import *
from plot_func.plotter import Plotter

class Single_sample_plotter(Plotter):
    def __init__(self,plot_df,anno_df):
        Plotter.__init__(self,plot_df,anno_df)
    def plot_dist(self,x_axis_column_name, scale):
        return Plotter.plot_dist(self,x_axis_column_name, scale)
    def plot_arr(self,x_axis_column_name, scale):
        plot_df = filter_by_scale(scale, self.plot_df)
        plot_df, custom_sort = get_group_range(plot_df, x_axis_column_name)
        plot_df = plot_df.groupby(by='group_range').count().reset_index()
        plot_df['Frequency'] = plot_df['isoform']/plot_df['isoform'].sum()
        fill_na_df = pd.DataFrame({'group_range':['{:.0%}-{:.0%}'.format(i/10,(i+1)/10) for i in range(1,10)]+['<=10%','>100%'],'Frequency':[0 for i in range(11)]}).set_index('group_range')
        plot_df = plot_df.append(fill_na_df.loc[fill_na_df.index.difference(pd.Index(plot_df['group_range']))].reset_index())
        plot_df = plot_df.sort_values(
            by=['group_range'], key=lambda col: custom_sort(col, ARR_ranges))
        fig = px.bar(plot_df, x='group_range',
                        y='Frequency')
        fig.update_layout(
            xaxis_title='Abundance Recovery Rate',
            yaxis_title='Frequency',
            width=fig_size['rec']['width'],height=fig_size['rec']['height'],template= themes['small_single']
        )
        return fig
    def plot_grouped_curve(self,x_axis_column_name, y_axis_column_names, scale):
        plot_df = filter_by_scale(scale, self.plot_df)
        plot_df, custom_sort = get_group_range(plot_df, x_axis_column_name)
        # fig = go.Figure()
        # figure_rows = math.ceil(math.sqrt(len(y_axis_column_names)))
        # figure_cols = math.ceil(len(y_axis_column_names)/ figure_rows)
        figure_rows = math.ceil(math.sqrt(len(y_axis_column_names)))
        figure_cols = math.ceil(len(y_axis_column_names)/ figure_rows)
        fig = make_subplots(rows=figure_rows, cols=figure_cols)
        for i in range(len(y_axis_column_names)):
            row_num = math.ceil((i+1)/figure_cols)
            col_num = i % figure_cols+1
            y_axis_column_name = y_axis_column_names[i]
            if ((y_axis_column_name in ['mrd']) & (x_axis_column_name == 'K_value')):
                group_series = plot_df.groupby(by='group_range').apply(lambda df: get_single_sample_metric(
                    y_axis_column_name, df['true_abund'], df['estimated_abund'],plot_df, True )).to_frame().reset_index()
            else:
                group_series = plot_df.groupby(by='group_range').apply(lambda df: get_single_sample_metric(
                    y_axis_column_name, df['true_abund'], df['estimated_abund'],plot_df)).to_frame().reset_index()
            group_series = group_series.rename(columns={0: y_axis_column_name}).sort_values(
                by=['group_range'], key=lambda col: custom_sort(col))
            fig.add_trace(go.Scatter(x=group_series['group_range'], y=group_series[y_axis_column_name],
                                        mode='lines+markers', name=on_plot_shown_label[y_axis_column_name]), row=row_num, col=col_num)
            fig.update_xaxes(
                title_text=on_plot_shown_label[x_axis_column_name],tickangle = 20, row=row_num, col=col_num)
            fig.update_yaxes(
                title_text=on_plot_shown_label[y_axis_column_name], row=row_num, col=col_num)
        fig.update_layout(
            autosize=False,
            width=fig_size['rec']['width']*figure_cols,height=fig_size['rec']['height']*figure_rows,template= themes['small_multi'])
        return fig
    def plot_corr_scatter(self,x_axis_column_name, y_axis_column_name,  scale):
        plot_df = filter_by_scale(scale, self.plot_df)
        plot_df = plot_df[(np.log2(plot_df[x_axis_column_name]+1) >=1) & (np.log2(plot_df[y_axis_column_name]+1) >= 1)]
        x = np.log2(plot_df[x_axis_column_name] + 1)
        y = np.log2(plot_df[y_axis_column_name] + 1)
        x_max = x.max()
        y_max = y.max()
        fig = go.Figure()
        x,y,density = get_density(x,y)
        fig.add_trace(go.Scattergl(x=x, y=y, mode='markers', name='Value',marker=dict(size=5,color=density,colorscale='viridis')))
        fig.add_trace(go.Histogram2dContour(x=x, y=y, name='Density',contours={'coloring':'none','showlabels':True}))
        fig.update_layout(showlegend=False, autosize=False,width=fig_size['square']['width'],height=fig_size['square']['height'],template= themes['large_single'])
        fig.update_xaxes(title_text='Log2(True abundance+1)',range=[1,max(x_max,y_max)])
        fig.update_yaxes(title_text='Log2(Estimated abundance+1)',range=[1,max(x_max,y_max)])
        return fig
    def plot_std_scatter(self,x_axis_column_name, y_axis_column_name, scale):
        plot_df = filter_by_scale(scale, self.plot_df)
        x = np.log2(plot_df[x_axis_column_name] + 1)
        y = plot_df[y_axis_column_name]
        x_max = x.max()
        fig = go.Figure()
        x,y,density = get_density(x,y)
        fig.add_trace(go.Scattergl(x=x, y=y, mode='markers', name='Value',marker=dict(size=5,color=density,colorscale='viridis')))
        fig.add_trace(go.Histogram2dContour(x=x, y=y, name='Density',contours={'coloring':'none','showlabels':True}))
        fig.update_layout(showlegend=False, autosize=False,width=fig_size['square']['width'],height=fig_size['square']['height'],template= themes['large_single'])
        fig.update_xaxes(title_text='Log2(Estimated abundance+1)',range=[1,x_max])
        fig.update_yaxes(title_text='std')
        return fig
    def plot_std_curve(self,x_axis_column_name,y_axis_column_name,scale):
        plot_df = filter_by_scale(scale, self.plot_df)
        n_bins = 1000
        degree = 5
        plot_df['range'] = pd.cut(plot_df[x_axis_column_name],bins=n_bins)
        grouped = plot_df[['range',y_axis_column_name]].groupby('range').mean().reset_index()
        grouped[x_axis_column_name] = [np.mean([interval.left,interval.right]) for interval in grouped['range']]
        grouped = grouped.dropna()
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=grouped[x_axis_column_name],y=grouped[y_axis_column_name],marker_color=color_schemes[0]))
        c = np.polynomial.polynomial.polyfit(grouped[x_axis_column_name],grouped[y_axis_column_name],degree)
        grouped['smoothed_{}'.format(y_axis_column_name)] = np.polynomial.polynomial.polyval(grouped[x_axis_column_name],c)
        fig.add_trace(go.Scatter(x=grouped[x_axis_column_name],y=grouped['smoothed_{}'.format(y_axis_column_name)],mode='lines'))
        auc = np.trapz(grouped['smoothed_{}'.format(y_axis_column_name)],grouped[x_axis_column_name])
        fig.add_annotation(x=grouped[x_axis_column_name].max()*0.95, y=grouped[y_axis_column_name].max()*0.95,
            text="{:.3f}".format(auc),showarrow=False)
        fig.update_layout(xaxis_title= 'Log2(Estimated abundance+1)',
            yaxis_title= 'std',autosize=False,width=fig_size['square']['width'],height=fig_size['square']['height'],template= themes['small_single'],showlegend=False)
        return fig
    def plot_corr_box_plot(self,x_axis_column_name,y_axis_column_name,scale):
        plot_df = filter_by_scale(scale, self.plot_df)
        plot_df = plot_df[(np.log2(plot_df[x_axis_column_name]+1) >=1) &  (np.log2(plot_df[y_axis_column_name]+1) >= 1)]
        df,_ = prepare_corr_box_plot_data(np.log2(plot_df[y_axis_column_name]+1),np.log2(plot_df[x_axis_column_name]+1))
        fig = go.Figure()
        # fig.add_trace(go.Box(x=df['true_abund'],y=df['estimated_abund'],boxpoints='all',jitter=0.3))
        fig.add_trace(go.Box(x=df['true_abund'],y=df['estimated_abund']))
        fig.update_layout(
            xaxis_title= 'Log2(True abundance+1)',
            yaxis_title= 'Log2(Estimated abundance+1)',
            autosize=False,showlegend=False,width=fig_size['square']['width'],height=fig_size['square']['height'],template= themes['small_single'])
        return fig

    def plot_resolution_entropy(self,scale):
        fig = go.Figure()
        plot_df = filter_by_scale(scale, self.plot_df)
        RE = [get_resolution_entropy(plot_df['estimated_abund'],100)]
        fig.add_trace(go.Bar(x=['Resolution Entropy'],
                    y=RE,marker_color=color_schemes[0]))
        fig.update_layout(
            width=fig_size['small_rec']['width'],height=fig_size['small_rec']['height'],template=themes['medium_single']
        )
        return fig
    def plot(self,plot_figure_name, scale):
        x_axis_column_name = single_sample_plot_figures[plot_figure_name]['x']
        y_axis_column_name = single_sample_plot_figures[plot_figure_name]['y']
        if y_axis_column_name == 'dist':
            fig = self.plot_dist(x_axis_column_name, scale)
        elif plot_figure_name == 'Histogram of Abundance Recovery Rate':
            fig = self.plot_arr(x_axis_column_name, scale)
        elif plot_figure_name in ["Statistics with different K values",'Statistics with different isoform lengths','Statistics with different numbers of exons','Statistics with different expression level']:
            fig = self.plot_grouped_curve(x_axis_column_name,y_axis_column_name,scale)
        elif plot_figure_name in ['Correlation of estimated abundance and ground truth']:
            fig = self.plot_corr_scatter(x_axis_column_name, y_axis_column_name, scale)
        elif plot_figure_name in ['Standard deviation vs estimated abundance scatter']:
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
