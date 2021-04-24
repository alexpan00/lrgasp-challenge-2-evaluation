import plotly.graph_objects as go
import plotly.io as pio
import numpy as np
import pandas as pd
import math

from static_data import ARR_ranges, on_plot_shown_label,fig_size,color_schemes,themes
from preprocess_util import *

def filter_by_scale(scale, plot_df):
    if scale == 'k>=1':
        plot_df = plot_df[plot_df['K_value'] >= 1]
    elif scale == 'k<1':
        plot_df = plot_df[plot_df['K_value'] < 1]
    elif scale == 'All':
        plot_df = plot_df
    return plot_df
def prepare_ranges(plot_df,groupby):
    if groupby == 'K_value':
        # ranges = K_value_ranges
        # plot_df.loc[:, 'group_range'] = pd.cut(
        #     plot_df[groupby], ranges).astype(str)
        # plot_df.loc[plot_df[groupby] > ranges[-1],
        #             'group_range'] = '>{}'.format(ranges[-1])
        # plot_df.loc[plot_df[groupby] == ranges[0],
        #             'group_range'] = ' {}'.format(ranges[0])
        # plot_df.loc[plot_df[groupby] < ranges[0],
        #             'group_range'] = '<{}'.format(ranges[0])
        # qcutted = pd.qcut(plot_df[plot_df[groupby]<1][groupby], 9,duplicates='drop')
        # categories = qcutted.cat.categories
        # qcutted_str = qcutted.astype(str)
        # qcutted_str[qcutted_str == str(categories[0])] = '(0, {}]'.format(categories[0].right)
        # qcutted_str[qcutted_str == str(categories[-1])] = '({}, 1)'.format(categories[-1].left)
        # plot_df.loc[plot_df[groupby]<1, 'group_range'] = qcutted_str
        # plot_df.loc[plot_df[groupby]>=1, 'group_range'] = '>= 1'
        if (len(plot_df[groupby].unique())<9):
            n_bins = len(plot_df[groupby].unique())
        else:
            n_bins = 9
        qcutted,categories = pd.qcut(plot_df[groupby],n_bins ,duplicates='drop',retbins=True)
        def custom_sort(col):
            vals = []
            for val in col.tolist():
                if ',' in val:
                    vals.append(float(val.split(',')[1][1:-1]))
                else:
                    # vals.append(float(val[2:]))
                    vals.append(float('inf'))
            return pd.Series(vals)
        return categories,None
    elif groupby in ['isoform_length']:
        def custom_sort(col):
            vals = []
            for val in col.tolist():
                if ',' in str(val):
                    vals.append(float(val.split(',')[1][1:-1]))
                else:
                    # vals.append(float(val[1:]))
                    vals.append(float('inf'))
            return pd.Series(vals)
        plot_df[groupby] = plot_df[groupby].astype(int)
       
        if plot_df[groupby].max() > 3000:
            max_threshold = 4000
            lower, higher = int(plot_df[groupby].min()), 4000
            step_size = 400
        else:
            max_threshold = 2100
            lower, higher = int(plot_df[groupby].min()), 2100
            step_size = 200
        # # max_threshold = np.ceil(np.percentile(plot_df[groupby], 80))
        # # lower, higher = int(plot_df.min()), int(plot_df.max())
        # # step_size = int(math.ceil((higher - lower)/n_bins))
        n_bins = 10
        
        edges = [lower] + list(
            range(step_size, higher+1, step_size))
        cutted,categories = pd.cut(
            plot_df.loc[plot_df[groupby] <= max_threshold, groupby], bins=edges,include_lowest=True,retbins=True)
        return categories,max_threshold
    else:
        plot_df[groupby] = plot_df[groupby].astype(int)
        max_threshold = np.ceil(np.percentile(plot_df[groupby], 90))
        if (len(plot_df.loc[plot_df[groupby] <= max_threshold, groupby].unique())<10):
            n_bins = len(plot_df.loc[plot_df[groupby] <= max_threshold, groupby].unique())
        else:
            n_bins = 10
        qcutted,categories = pd.qcut(plot_df.loc[plot_df[groupby] <= max_threshold, groupby], n_bins,labels=False,duplicates='drop',retbins=True)
        # lower, higher = temp_df.min(), temp_df.max()
        # if (len(plot_df.loc[plot_df[groupby] <= max_threshold, groupby].unique())<10):
        #     n_bins = len(plot_df.loc[plot_df[groupby] <= max_threshold, groupby].unique())
        # else:
        #     n_bins = 10
        # edges = list(
        #     range(int(lower-1), int(higher), int(math.ceil((higher - lower)/n_bins))))
        # edges.append(higher)
        # plot_df.loc[plot_df[groupby] <= max_threshold, 'group_range'] = pd.cut(
        #     temp_df, bins=edges).astype('str')
        # plot_df.loc[plot_df[groupby] > max_threshold,
        #             'group_range'] = '>{}'.format(max_threshold)
        return categories,max_threshold
def get_group_range(plot_df, groupby,ranges,max_threshold):
    if groupby == 'K_value':
        # ranges = K_value_ranges
        # plot_df.loc[:, 'group_range'] = pd.cut(
        #     plot_df[groupby], ranges).astype(str)
        # plot_df.loc[plot_df[groupby] > ranges[-1],
        #             'group_range'] = '>{}'.format(ranges[-1])
        # plot_df.loc[plot_df[groupby] == ranges[0],
        #             'group_range'] = ' {}'.format(ranges[0])
        # plot_df.loc[plot_df[groupby] < ranges[0],
        #             'group_range'] = '<{}'.format(ranges[0])
        # qcutted = pd.cut(plot_df[plot_df[groupby]<1][groupby], bins=ranges)
        # categories = qcutted.cat.categories
        # qcutted_str = qcutted.astype(str)
        # qcutted_str[qcutted_str == str(categories[0])] = '(0, {}]'.format(categories[0].right)
        # qcutted_str[(qcutted_str == str(categories[-1]))|(qcutted_str=='nan')] = '({}, 1)'.format(categories[-1].left)
        # plot_df.loc[plot_df[groupby]<1, 'group_range'] = qcutted_str
        # plot_df.loc[plot_df[groupby]>=1, 'group_range'] = '>= 1'

        qcutted = pd.cut(plot_df[groupby], bins=ranges)
        categories = qcutted.cat.categories
        qcutted_str = qcutted.astype(str)
        # qcutted_str[qcutted_str == str(categories[0])] = '(1, {}]'.format(categories[0].right)
        qcutted_str[qcutted_str == str(categories[-1])] = '>{}'.format(categories[-1].left)
        plot_df['group_range'] = qcutted_str
        plot_df = plot_df[plot_df['group_range']!='nan']
        def custom_sort(col):
            vals = []
            for val in col.tolist():
                if ',' in val:
                    vals.append(float(val.split(',')[1][1:-1]))
                else:
                    # vals.append(float(val[2:]))
                    vals.append(float('inf'))
            return pd.Series(vals)
    elif groupby == 'arr':
        def custom_sort(col, ranges):
            vals = []
            for val in col.tolist():
                if val == '<={:.0%}'.format(ranges[1]):
                    vals.append('0')
                else:
                    vals.append(val)
            return pd.Series(vals)
        ranges = ARR_ranges
        plot_df.loc[:, 'group_range'] = pd.cut(plot_df[groupby], ranges).apply(
            lambda x: '{:.0%}-{:.0%}'.format(x.left, x.right) if x is not None else 'nan').astype(str)
        plot_df.loc[plot_df[groupby] > ranges[-1],
                    'group_range'] = '>{:.0%}'.format(ranges[-1])
        plot_df.loc[plot_df[groupby] <= ranges[1],
                    'group_range'] = '<={:.0%}'.format(ranges[1])
        plot_df = plot_df.dropna()
    else:
        def custom_sort(col):
            vals = []
            for val in col.tolist():
                if ',' in val:
                    vals.append(float(val.split(',')[1][1:-1]))
                else:
                    # vals.append(float(val[1:]))
                    vals.append(float('inf'))
            return pd.Series(vals)
        plot_df[groupby] = plot_df[groupby].astype(int)
        ranges = ranges.tolist()
        if ((len(ranges)==2) & (ranges[0]==ranges[1])):
            ranges[0] = ranges[0] - 1
        plot_df.loc[plot_df[groupby] <= max_threshold, 'group_range'] = pd.cut(plot_df.loc[plot_df[groupby] <= max_threshold, groupby], bins=ranges)
        group_range_series = plot_df.loc[plot_df[groupby] <= max_threshold, 'group_range'].apply(lambda x: str(pd.Interval(left=int(round(x.left)), right=int(round(x.right)))) if type(x)==pd.Interval else x)
        plot_df.loc[plot_df[groupby] <= max_threshold, 'group_range'] = group_range_series
        if (plot_df.loc[plot_df[groupby] > max_threshold].shape[0]>0):
            plot_df.loc[plot_df[groupby] > max_threshold,
                            'group_range'] = '>{}'.format(max_threshold)
        plot_df['group_range'] = plot_df['group_range'].astype(str)
        plot_df = plot_df[plot_df['group_range']!='nan']
        if ((len(ranges)==2) & (ranges[0]==ranges[1]-1)):
            plot_df['group_range'] = str(ranges[1])
    return plot_df, custom_sort
def get_density(x,y):
    from scipy.stats import gaussian_kde
    try:
        xy = np.vstack([x,y])
        z = gaussian_kde(xy)(xy)
    except:
        z = y.value_counts()[y]
    # idx = z.argsort()
    # x, y, z = x[idx], y[idx], z[idx]
    return x,y,z