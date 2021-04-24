import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd
import math

from static_data import ARR_ranges, on_plot_shown_label,fig_size,color_schemes,themes
from sklearn.metrics import roc_curve,precision_recall_curve,average_precision_score,roc_auc_score
from preprocess_util import *
from plot_func.plot_util import *
from plot_func.plotter import Plotter

class Multi_sample_plotter(Plotter):
    def __init__(self,plot_df,anno_df):
        Plotter.__init__(self,plot_df,anno_df)
    def plot_dist(self,x_axis_column_name, scale):
        return Plotter.plot_dist(self,x_axis_column_name, scale)
    def plot_arr(self,x_axis_column_name, scale):
        fig = make_subplots(rows=2, cols=1,horizontal_spacing=0.1,vertical_spacing=0.1,row_titles=['Condition 1','Condition 2'])
        plot_df = filter_by_scale(scale, self.plot_df)
        arr_columns = [x for x in list(plot_df.columns) if 'arr_' in x]
        arr_dfs = []
        for i in range(len(arr_columns)):
            arr_df = plot_df[[arr_columns[i]]].rename(
                columns={arr_columns[i]: 'arr'})
            arr_df, custom_sort = get_group_range(arr_df, x_axis_column_name)
            arr_df = arr_df.groupby(by='group_range').count()
            arr_df['Frequency_{}'.format(i+1)] = arr_df['arr']/arr_df['arr'].sum()
            arr_dfs.append(arr_df['Frequency_{}'.format(i+1)])
        # plot_df = pd.concat(arr_dfs).sort_values(
        #     by=['group_range'], key=lambda col: custom_sort(col, ARR_ranges))
        for row,start_col,end_col in zip([1,2],[0,len(arr_dfs)//2],[len(arr_dfs)//2,len(arr_dfs)]):
            plot_df = pd.concat(arr_dfs[start_col:end_col],axis=1)
            plot_df = plot_df.reset_index()
            fill_na_df = pd.DataFrame({'group_range':['{:.0%}-{:.0%}'.format(i/10,(i+1)/10) for i in range(1,10)]+['<=10%','>100%'],'Frequency':[0 for i in range(11)],'Error':[0 for i in range(11)]}).set_index('group_range')
            plot_df = plot_df.append(fill_na_df.loc[fill_na_df.index.difference(pd.Index(plot_df['group_range']))].reset_index())
            plot_df['Frequency'] = plot_df.mean(axis=1)
            plot_df['Error'] = plot_df.std(axis=1)
            plot_df = plot_df.sort_values(by=['group_range'], key=lambda col: custom_sort(col, ARR_ranges))
            fig.add_trace(go.Bar(x=plot_df['group_range'],y=plot_df['Frequency'],error_y=dict(type='data',array=plot_df['Error'])),row=row,col=1)
            # fig.update_traces(box_visible=True, meanline_visible=True)
        # fig.update_traces(box_visible=True, meanline_visible=True)
        fig.update_xaxes(title_text='Abundance Recovery Rate')
        fig.update_yaxes(title_text='Frequency')
        fig.update_layout(
            width=fig_size['rec']['width'],height=fig_size['rec']['height']*2,template= themes['small_multi']
        )
        return fig
    def plot_stats_box(self,scale):
        plot_df = filter_by_scale(scale, self.plot_df)
        [cond1_metric_dicts,cond2_metric_dicts] = prepare_stats_box_plot_data(plot_df)
        col_num = len(cond1_metric_dicts)
        fig = make_subplots(rows=1, cols=col_num,horizontal_spacing=0.1,vertical_spacing=0.1,)
        for i,cond1_metric_dict,cond2_metric_dict in zip(range(col_num),cond1_metric_dicts,cond2_metric_dicts):
            M = on_plot_shown_label[cond1_metric_dict['Metric']]
            fig.add_trace(go.Bar(x=[M],y=[cond1_metric_dict['Mean']],error_y=dict(type='data',array=[cond1_metric_dict['Error']]),name='Condition 1',showlegend=False,marker_color=color_schemes[0]),row=1,col=i+1)
            fig.add_trace(go.Bar(x=[M],y=[cond2_metric_dict['Mean']],error_y=dict(type='data',array=[cond2_metric_dict['Error']]),name='Condition 2',showlegend=False,marker_color=color_schemes[1]),row=1,col=i+1)
        fig.update_traces(showlegend=True,col=1,row=1)
        fig.update_layout(
            width=fig_size['small_rec']['width']*col_num,
            height=fig_size['small_rec']['height'],
            boxmode='group', 
            template= themes['small_multi'])
        return fig
    def plot_roc(self,x_axis_column_name,y_axis_column_name,scale):
        plot_df = filter_by_scale(scale, self.plot_df)
        fpr,tpr,_ = roc_curve((plot_df['true_alfc'] > plot_df['true_alfc'].median()),plot_df['alfc'])
        auc = roc_auc_score((plot_df['true_alfc'] > plot_df['true_alfc'].median()),plot_df['alfc'])
        df = pd.DataFrame({'fpr':fpr,'tpr':tpr})
        fig = go.Figure(go.Scatter(x=df[x_axis_column_name],y=df[y_axis_column_name],mode='lines'))
        fig.update_layout(
            xaxis_title=on_plot_shown_label[x_axis_column_name],
            yaxis_title=on_plot_shown_label[y_axis_column_name],
            width=fig_size['square']['width'],
            height=fig_size['square']['height'],
            template= themes['large_single']
        )
        fig.update_xaxes(range=[0,1])
        fig.update_yaxes(range=[0,1])
        fig.add_annotation(x=0.9, y=0.05,
            text="AUC: {:.3f}".format(auc),showarrow=False)
        return fig
    def plot_pr(self,x_axis_column_name,y_axis_column_name,scale):
        plot_df = filter_by_scale(scale, self.plot_df)
        precision,recall,_ = precision_recall_curve((plot_df['true_alfc'] > plot_df['true_alfc'].median()),plot_df['alfc'])
        aps = average_precision_score((plot_df['true_alfc'] > plot_df['true_alfc'].median()),plot_df['alfc'])
        df = pd.DataFrame({'precision':precision,'recall':recall})
        fig = go.Figure(go.Scatter(x=df[x_axis_column_name],y=df[y_axis_column_name],mode='lines'))
        fig.update_layout(
            xaxis_title=on_plot_shown_label[x_axis_column_name],
            yaxis_title=on_plot_shown_label[y_axis_column_name],
            width=fig_size['square']['width'],
            height=fig_size['square']['height'],
            template= themes['large_single']
        )
        fig.update_xaxes(range=[0,1])
        fig.update_yaxes(range=[0,1])
        fig.add_annotation(x=0.9, y=0.05,
            text="Average precision: {:.3f}".format(aps),showarrow=False)
        return fig
    def plot_grouped_violin(self,x_axis_column_name, y_axis_column_names, scale):
        plot_df = filter_by_scale(scale, self.plot_df)
        plot_df, custom_sort = get_group_range(plot_df, x_axis_column_name)
        # fig = go.Figure()
        figure_rows = math.ceil(math.sqrt(len(y_axis_column_names)))
        figure_cols = math.ceil(len(y_axis_column_names)/ figure_rows)

        fig = make_subplots(rows=figure_rows, cols=figure_cols, vertical_spacing=0.2, horizontal_spacing=0.1)
        for i in range(len(y_axis_column_names)):
            row_num = math.ceil((i+1)/figure_cols)
            col_num = i % figure_cols+1
            y_axis_column_name = y_axis_column_names[i]
            if y_axis_column_name in ['precision', 'recall', 'accuracy', 'auc', 'f1', 'CM','RM']:
                group_series = plot_df.groupby(by='group_range').apply(
                    lambda df: prepare_grouped_violin_data(y_axis_column_name, df)).to_frame().reset_index()
                group_series = group_series.rename(columns={0: y_axis_column_name}).sort_values(
                    by=['group_range'], key=lambda col: custom_sort(col))
                fig.add_trace(go.Scatter(x=group_series['group_range'], y=group_series[y_axis_column_name],
                                        name=on_plot_shown_label[y_axis_column_name], mode='lines+markers'), row=row_num, col=col_num)
            else:
                if ((y_axis_column_name in ['mrd']) & (x_axis_column_name=='K_value')):
                    group_series = plot_df.groupby(by='group_range').apply(lambda df: prepare_grouped_violin_data(
                        y_axis_column_name, df,True)).explode().to_frame().reset_index()
                else:
                    group_series = plot_df.groupby(by='group_range').apply(lambda df: prepare_grouped_violin_data(
                                y_axis_column_name, df)).explode().to_frame().reset_index()
                group_series = group_series.rename(columns={0: y_axis_column_name})
                group_series[y_axis_column_name] = group_series[y_axis_column_name].astype(float)
                mean_series = group_series.groupby('group_range').mean()
                error_series = group_series.groupby('group_range').std()
                error_bar_group_series = pd.DataFrame({'Mean':mean_series[y_axis_column_name],'Error':error_series[y_axis_column_name]},index=mean_series.index).reset_index()
                error_bar_group_series = error_bar_group_series.sort_values(
                            by=['group_range'], key=lambda col: custom_sort(col))
                fig.add_trace(go.Bar(x=error_bar_group_series['group_range'],y=error_bar_group_series['Mean'],error_y=dict(type='data',array=error_bar_group_series['Error']),name=on_plot_shown_label[y_axis_column_name],showlegend=False),row=row_num,col=col_num)
                
                # fig.add_trace(go.Violin(x=group_series['group_range'], y=group_series[y_axis_column_name],
                #                         name=on_plot_shown_label[y_axis_column_name], box_visible=True, meanline_visible=True), row=row_num, col=col_num)
            fig.update_xaxes(
                title_text=on_plot_shown_label[x_axis_column_name],tickangle = 45, row=row_num, col=col_num)
            fig.update_yaxes(
                title_text=on_plot_shown_label[y_axis_column_name], row=row_num, col=col_num)
        # fig.update_traces(line=dict(width=3))
        fig.update_layout(
            autosize=False,
            # legend=dict(yanchor="top",y=0.1, xanchor="right",x=0.6),
            width=fig_size['rec']['width']*figure_cols,
            height=fig_size['rec']['height']*figure_rows,template= themes['small_multi'])
        
        return fig
    def plot_corr_scatter(self,x_axis_column_names, y_axis_column_names, scale):
        plot_df = filter_by_scale(scale, self.plot_df)
        subplot_titles = ['Condition {}'.format(i+1) for i in range(len(y_axis_column_names))]
        fig = make_subplots(rows=1, cols=len(y_axis_column_names),subplot_titles=subplot_titles)
        # estimated_columns = [x for x in list(plot_df.columns) if 'estimated_abund_' in x]
        # true_columns = [x for x in list(plot_df.columns) if 'true_abund_' in x]
        x_maxs,y_maxs = [],[]
        col_num = 0
        for x_axis_column_name,y_axis_column_name in zip(x_axis_column_names,y_axis_column_names):
            col_num += 1
            df = plot_df[(plot_df[x_axis_column_name] >=1) &  (plot_df[y_axis_column_name] >= 1)]
            x,y,density = get_density(df[x_axis_column_name],df[y_axis_column_name])
            fig.add_trace(go.Scatter(x=x, y=y, mode='markers', name='Value',marker=dict(size=5,color=density,colorscale='viridis')),row=1, col=col_num)
            fig.add_trace(go.Histogram2dContour(x=x, y=y, name='Density',contours={'coloring':'none','showlabels':True}),row=1, col=col_num)
            x_maxs.append(x.max())
            y_maxs.append(y.max())
        x_title = 'Log2(true abundance+1)'
        y_title = 'Log2(Estimated abundance+1)'
        fig.update_xaxes(title_text=x_title,range=[1,max(x_maxs+y_maxs)])
        fig.update_yaxes(title_text=y_title,range=[1,max(x_maxs+y_maxs)])
        fig.update_layout(showlegend=False,autosize=False,width=fig_size['square']['width']*2,height=fig_size['square']['height'],template= themes['large_single'])
        return fig
    def plot_std_scatter(self,x_axis_column_names, y_axis_column_names,  scale):
        plot_df = filter_by_scale(scale,self.plot_df)
        subplot_titles = ['Condition {}'.format(i+1) for i in range(len(y_axis_column_names))]
        fig = make_subplots(rows=1, cols=len(y_axis_column_names),subplot_titles=subplot_titles)
        colorbars = [dict(len=1.05, x=0.45,y=0.49),dict(len=1.05,  x=1.0 , y=0.49)]
        col_num = 0
        x_maxs,y_maxs = [],[]
        for x_axis_column_name,y_axis_column_name in zip(x_axis_column_names,y_axis_column_names):
            col_num += 1
            x = plot_df[x_axis_column_name]
            y = plot_df[y_axis_column_name]
            x,y,density = get_density(x,y)
            fig.add_trace(go.Scattergl(x=x, y=y, mode='markers', name='Value',marker=dict(size=5,color=density,colorscale='viridis')),row=1, col=col_num)
            fig.add_trace(go.Histogram2dContour(x=x, y=y, name='Density',contours={'coloring':'none','showlabels':True}),row=1, col=col_num)
            x_maxs.append(max(x))
            y_maxs.append(max(y))
        x_title = 'Log2(Estimated abundance+1)'
        y_title = 'std'
        fig.update_xaxes(title_text=x_title,range=[1,max(x_maxs)])
        fig.update_yaxes(title_text=y_title,range=[0,max(y_maxs)])
        fig.update_layout(showlegend=False,autosize=False,width=fig_size['square']['width']*2,height=fig_size['square']['height'],template= themes['small_multi'])
        return fig
    def plot_multi_sample_std_curve(self,x_axis_column_names,y_axis_column_names,scale):
        plot_df = filter_by_scale(scale, self.plot_df)
        n_bins = 1000
        degree = 5
        subplot_titles = ['Condition {}'.format(i+1) for i in range(len(y_axis_column_names))]
        # fig = make_subplots(cols=len(x_axis_column_names),rows=1,subplot_titles=subplot_titles)
        fig = go.Figure()
        x_mins,y_mins,x_maxs,y_maxs,aucs = [],[],[],[],[]
        for i,x_axis_column_name,y_axis_column_name in zip(range(len(x_axis_column_names)),x_axis_column_names,y_axis_column_names):
            df = plot_df.copy()
            df['range'] = pd.cut(df[x_axis_column_name],bins=n_bins)
            grouped = df[['range',y_axis_column_name]].groupby('range').mean().reset_index()
            grouped[x_axis_column_name] = [np.mean([interval.left,interval.right]) for interval in grouped['range']]
            grouped = grouped.dropna()
            # fig.add_trace(go.Scatter(x=grouped[x_axis_column_name],y=grouped[y_axis_column_name],marker_color=color_schemes[0]),col=i+1,row=1)
            c = np.polynomial.polynomial.polyfit(grouped[x_axis_column_name],grouped[y_axis_column_name],degree)
            grouped['smoothed_{}'.format(y_axis_column_name)] = np.polynomial.polynomial.polyval(grouped[x_axis_column_name],c)
            fig.add_trace(go.Scattergl(x=grouped[x_axis_column_name],y=grouped['smoothed_{}'.format(y_axis_column_name)],mode='lines',marker_color=color_schemes[i],name=subplot_titles[i]))
            auc = np.trapz(grouped['smoothed_{}'.format(y_axis_column_name)],grouped[x_axis_column_name])
            aucs.append(auc)
            x_mins.append(grouped[x_axis_column_name].min())
            x_maxs.append(grouped[x_axis_column_name].max())
            y_mins.append(grouped[y_axis_column_name].min())
            y_maxs.append(grouped[y_axis_column_name].max())
        fig.add_annotation(x=max(x_maxs)*0.8, y=max(y_maxs)*0.85,text='ASDC',showarrow=False)
        for i in range(len(x_axis_column_names)):
            fig.add_annotation(x=max(x_maxs)*0.8, y=max(y_maxs)*(0.8-i*0.1),
                text="Condition {}: {:.3f}".format(i+1,aucs[i]),showarrow=False)
        fig.update_xaxes(range=[min(x_mins),max(x_maxs)])
        fig.update_yaxes(range=[min(y_mins),max(y_maxs)])
        fig.update_layout(
            xaxis_title= 'Log2(Estimated abundance+1)',
            yaxis_title= 'std',
            autosize=False,showlegend=True,width=fig_size['square']['width']*2,height=fig_size['square']['height'],template= themes['large_multi'])
        return fig
    def plot_consistency_measure_curve(self,x_axis_column_names,y_axis_column_names,scale):
        plot_df = filter_by_scale(scale, self.plot_df)
        CM_list,C_ranges = prepare_consistency_measure_plot_data(plot_df)
        auc = np.trapz(CM_list,C_ranges)
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=C_ranges,y=CM_list,mode='lines',name='Consistency Measure'))
        fig.add_annotation(x=np.max(C_ranges)*0.8, y=np.max(CM_list)*0.85,text='ACMC',showarrow=False)
        fig.add_annotation(x=np.max(C_ranges)*0.8, y=np.max(CM_list)*0.8,
                text="{:.3f}".format(auc),showarrow=False)
        fig.update_layout(
            xaxis_title= 'C threshold',
            yaxis_title= 'Consistency Measure',
            autosize=False,showlegend=True,width=fig_size['square']['width'],height=fig_size['square']['height'],template= themes['small_single'])
        return fig
    def plot_resolution_entropy(self,scale):
        fig = make_subplots(rows=2, cols=1,horizontal_spacing=0.1,vertical_spacing=0.1,row_titles=['Condition 1','Condition 2'])
        plot_df = filter_by_scale(scale, plot_df)
        [cond1_metric_dicts,cond2_metric_dicts] = prepare_stats_box_plot_data(plot_df)
        cond1_metric_dicts = [i for i in cond1_metric_dicts if i['Metric'] == 'RE']
        cond2_metric_dicts = [i for i in cond1_metric_dicts if i['Metric'] == 'RE']
        col_num = len(cond1_metric_dicts)
        for i,cond1_metric_dict,cond2_metric_dict in zip(range(col_num),cond1_metric_dicts,cond2_metric_dicts):
            M = on_plot_shown_label[cond1_metric_dict['Metric']]
            fig.add_trace(go.Bar(x=[M],y=[cond1_metric_dict['Mean']],error_y=dict(type='data',array=[cond1_metric_dict['Error']]),showlegend=False,marker_color=color_schemes[j]),row=1,col=i+1)
            fig.add_trace(go.Bar(x=[M],y=[cond2_metric_dict['Mean']],error_y=dict(type='data',array=[cond2_metric_dict['Error']]),showlegend=False,marker_color=color_schemes[j]),row=2,col=i+1)
        fig.update_traces(showlegend=True,col=1,row=1)
        fig.update_layout(
            width=540,
            height=540,
            boxmode='group', 
            template= themes['small_multi'])
        return fig
    def plot_multi_corr_box_plot(self,x_axis_column_names,y_axis_column_names,scale):
        plot_df = filter_by_scale(scale, self.plot_df)
        x_mins,x_maxs,y_mins,y_maxs = [],[],[],[]
        subplot_titles = ['Condition {}'.format(i+1) for i in range(len(y_axis_column_names))]
        fig = make_subplots(rows=1, cols=len(y_axis_column_names),subplot_titles=subplot_titles)

        estimated_columns = [x for x in list(plot_df.columns) if 'estimated_abund_' in x]
        true_columns = [x for x in list(plot_df.columns) if 'true_abund_' in x]
        for col,start_col,end_col in zip([1,2],[0,len(estimated_columns)//2],[len(estimated_columns)//2,len(estimated_columns)]):
            df = pd.DataFrame()
            df['estimated_abund'] = np.log2(plot_df[estimated_columns[start_col:end_col]].values.flatten()+1)
            df['true_abund'] = np.log2(plot_df[true_columns[start_col:end_col]].values.flatten()+1)
            df = df[(df['true_abund'] >=1) &  (df['estimated_abund'] >= 1)]
            df,_ = prepare_corr_box_plot_data(df['estimated_abund'], df['true_abund'])
            fig.add_trace(go.Box(x=df['true_abund'],y=df['estimated_abund']),col=col,row=1)
            # fig.add_trace(go.Box(x=df['true_abund'],y=df['estimated_abund'],boxpoints='all',jitter=0.3),col=col,row=1)
        x_title = 'Log2(True abundance+1)'
        y_title = 'Log2(Estimated abundance+1)'
        fig.update_xaxes(title_text=x_title)
        fig.update_yaxes(title_text=y_title)
        fig.update_layout(
            autosize=False,showlegend=False,width=fig_size['square']['width']*2,height=fig_size['square']['height'],template= themes['small_single'])
        return fig
    def plot(self,plot_figure_name,scale,ground_truth_given):
        if (ground_truth_given):
            x_axis_column_name = multi_sample_diff_condition_with_ground_truth_plot_figures[plot_figure_name]['x']
            y_axis_column_name = multi_sample_diff_condition_with_ground_truth_plot_figures[plot_figure_name]['y']
        else:
            x_axis_column_name = multi_sample_diff_condition_without_ground_truth_plot_figures[plot_figure_name]['x']
            y_axis_column_name = multi_sample_diff_condition_without_ground_truth_plot_figures[plot_figure_name]['y']
        if plot_figure_name == 'ROC curves for performance of quantification':
            fig = self.plot_roc(x_axis_column_name,y_axis_column_name,scale)
        elif plot_figure_name == 'PR curves for performance of quantification':
            fig = self.plot_pr(x_axis_column_name,y_axis_column_name,scale)
        elif plot_figure_name in ["Statistics with different K values",'Statistics with different isoform lengths','Statistics with different numbers of exons','Statistics with different expression level']:
            fig = self.plot_grouped_violin(x_axis_column_name,y_axis_column_name,scale)
        elif plot_figure_name in ['Correlation of estimated abundance and ground truth']:
            fig = self.plot_corr_scatter(x_axis_column_name, y_axis_column_name, scale)
        elif plot_figure_name in ['Standard deviation vs estimated abundance scatter']:
            fig = self.plot_std_scatter(x_axis_column_name, y_axis_column_name, scale)
        elif y_axis_column_name == 'dist':
            fig = self.plot_dist(x_axis_column_name, scale)
        elif plot_figure_name == 'Histogram of Abundance Recovery Rate':
            fig = self.plot_arr(x_axis_column_name, scale)
        elif plot_figure_name == 'Estimation Error for different conditions':
            fig = self.plot_stats_box(scale)
        elif plot_figure_name == 'Standard deviation vs estimated abundance curve':
            fig = self.plot_multi_sample_std_curve(x_axis_column_name,y_axis_column_name,scale)
        elif plot_figure_name == 'Consistency Measure curve':
            fig = self.plot_consistency_measure_curve(x_axis_column_name,y_axis_column_name,scale)
        elif plot_figure_name == 'Correlation Boxplot of estimated abundance and ground truth':
            fig = self.plot_multi_corr_box_plot(x_axis_column_name,y_axis_column_name,scale)
        elif plot_figure_name == 'Resolution Entropy':
            fig = self.plot_resolution_entropy(scale)
        try:
            fig.update_layout(title=plot_figure_name, title_x=0.5)
            fig.update_xaxes(exponentformat='e',automargin=True)
            fig.update_yaxes(exponentformat='e',automargin=True)
        except:
            print(plot_figure_name)
        return fig