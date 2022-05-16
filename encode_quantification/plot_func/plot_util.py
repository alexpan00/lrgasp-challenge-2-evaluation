import plotly.graph_objects as go
import plotly.io as pio
import numpy as np
import pandas as pd
import math

from static_data import *
from preprocess_util import *


def define_theme():
    # naming a layout theme for future reference
    pio.templates["encode"] = go.layout.Template(
        layout_colorway=color_schemes,
        data_scatter=[dict(line=dict(width=5))]
    )
    pio.templates["large"] = go.layout.Template(
        layout_font=dict(family="Helvetica", size=16),
        layout_title_font=dict(family="Helvetica", size=19),
    )
    pio.templates["ultralarge"] = go.layout.Template(
        layout_font=dict(family="Helvetica", size=16),
        layout_title_font=dict(family="Helvetica", size=19),
    )
    pio.templates["medium"] = go.layout.Template(
        layout_font=dict(family="Helvetica", size=16),
        layout_title_font=dict(family="Helvetica", size=19),
    )

    # pio.templates.default = "encode"
    pio.templates.default = "presentation+encode"


def define_write_to_file_theme():
    # naming a layout theme for future reference
    pio.templates["encode"] = go.layout.Template(
        layout_colorway=color_schemes,
        layout_font=dict(family="Arial Black", size=22),
        layout_title_font=dict(family="Arial Black", size=27),
        data_scatter=[dict(line=dict(width=5))]
    )
    pio.templates["large"] = go.layout.Template(
        layout_font=dict(family="Arial Black", size=35),
        layout_title_font=dict(family="Arial Black", size=40),

    )
    pio.templates["ultralarge"] = go.layout.Template(
        layout_font=dict(family="Arial Black", size=35),
        layout_title_font=dict(family="Arial Black", size=40),
    )
    pio.templates["medium"] = go.layout.Template(
        layout_font=dict(family="Arial Black", size=25),
        layout_title_font=dict(family="Arial Black", size=30),
    )

    # pio.templates.default = "encode"
    pio.templates.default = "presentation+encode"


def filter_by_scale(scale, plot_df):
    if scale == 'k>=1':
        plot_df = plot_df[plot_df['K_value'] >= 1]
    elif scale == 'k<1':
        plot_df = plot_df[plot_df['K_value'] < 1]
    elif scale == 'All':
        plot_df = plot_df
    return plot_df


def get_k_val_dist(plot_df, groupby):
    if groupby == 'K_value':
        ranges = condition_number_ranges
        cutted = pd.cut(
            plot_df[groupby], ranges, include_lowest=True)
        categories = cutted.cat.categories
        plot_df.loc[:, 'group_range'] = cutted.astype(str)
        plot_df.loc[plot_df[groupby] > ranges[-1],
                    'group_range'] = '>{}'.format(ranges[-1])
        plot_df.loc[plot_df['group_range'] == str(
            categories[0]), 'group_range'] = '[{}, {}]'.format(ranges[0], categories[0].right)

        def custom_sort(col):
                vals = []
                for val in col.tolist():
                    if ',' in val:
                        vals.append(float(val.split(',')[1][1:-1]))
                    else:
                        vals.append(float(val[1:]))
                return pd.Series(vals)
        return plot_df, custom_sort
    elif groupby == 'num_isoforms':
        plot_df = plot_df[['gene', groupby]].drop_duplicates()
        def custom_sort(col):
            vals = []
            for val in col.tolist():
                if ',' in val:
                    vals.append(float(val.split(',')[1][1:-1]))
                else:
                    # vals.append(float(val[1:]))
                    vals.append(float('inf'))
            return pd.Series(vals)
        ranges = num_isoforms_range
        cutted = pd.cut(plot_df[groupby], ranges,include_lowest=True)
        categories = cutted.cat.categories
        plot_df.loc[:, 'group_range'] = cutted.apply(lambda x:str(x)).astype(str)
        plot_df.loc[plot_df[groupby] > ranges[-1],
                    'group_range'] = '>{}'.format(ranges[-1])
        plot_df.loc[plot_df['group_range'] == str(
            categories[0]), 'group_range'] = '[{}, {})'.format(int(ranges[0]), int(categories[0].right))
        # Only cut 80% of data and leave 20% together
        # temp_df = plot_df
        # temp_df[groupby] = plot_df[groupby].astype(int)
        # max_threshold = np.ceil(np.percentile(plot_df[groupby], 80))
        # temp_df = plot_df.loc[plot_df[groupby] <= max_threshold, groupby]
        # lower, higher = temp_df.min(), temp_df.max()
        # n_bins = 10
        # edges = list(
        #     range(int(lower-1), int(higher), int(math.ceil((higher - lower)/n_bins))))
        # edges.append(higher)
        # plot_df.loc[plot_df[groupby] <= max_threshold, 'group_range'] = pd.cut(
        #     temp_df, bins=edges).astype('str')
        # plot_df.loc[plot_df[groupby] > max_threshold,
        #             'group_range'] = '>{}'.format(max_threshold)
        return plot_df, custom_sort


def get_group_range(plot_df, groupby):
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
        ranges = condition_number_ranges
        cutted = pd.cut(plot_df[groupby], ranges,include_lowest=True)
        categories = cutted.cat.categories
        plot_df.loc[:, 'group_range'] = cutted.astype(str)
        plot_df.loc[plot_df[groupby] > ranges[-1],
                    'group_range'] = '>{}'.format(ranges[-1])
        plot_df.loc[plot_df['group_range'] == str(categories[0]),'group_range'] = '[{},{}]'.format(ranges[0],categories[0].right)

        # qcutted = pd.qcut(plot_df[groupby], 9,duplicates='drop')
        # categories = qcutted.cat.categories
        # qcutted_str = qcutted.astype(str)
        # # qcutted_str[qcutted_str == str(categories[0])] = '(1, {}]'.format(categories[0].right)
        # plot_df['group_range'] = qcutted_str
        plot_df = plot_df[plot_df['group_range'] != 'nan']
        # plot_df.loc[plot_df[groupby]>=1, 'group_range'] = '>= 1'

        def custom_sort(col):
            vals = []
            for val in col.tolist():
                if ',' in val:
                    vals.append(float(val.split(',')[1][1:-1]))
                else:
                    # vals.append(float(val[2:]))
                    vals.append(float('inf'))
            return pd.Series(vals)
        return plot_df, custom_sort
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
        return plot_df, custom_sort
    elif groupby in ['true_abund','ave_true_abund', 'log2_true_abund','ave_estimated_abund#1','ave_estimated_abund#2']:
        def custom_sort(col):
            vals = []
            for val in col.tolist():
                if ',' in val:
                    vals.append(float(val.split(',')[1][1:-1]))
                else:
                    # vals.append(float(val[1:]))
                    vals.append(float('inf'))
            return pd.Series(vals)
        ranges = abund_range
        cutted = pd.cut(plot_df[groupby], ranges,include_lowest=True)
        categories = cutted.cat.categories
        plot_df.loc[:, 'group_range'] = cutted.astype(str)
        plot_df.loc[plot_df[groupby] > ranges[-1],
                    'group_range'] = '>{}'.format(ranges[-1])
        plot_df.loc[plot_df['group_range'] == str(categories[0]),'group_range'] = '[{},{}]'.format(ranges[0],categories[0].right)

        # qcutted = pd.qcut(plot_df[groupby], 9,duplicates='drop')
        # categories = qcutted.cat.categories
        # qcutted_str = qcutted.astype(str)
        # # qcutted_str[qcutted_str == str(categories[0])] = '(1, {}]'.format(categories[0].right)
        # plot_df['group_range'] = qcutted_str
        plot_df = plot_df[plot_df['group_range'] != 'nan']
        # n_bins = 6
        # plot_df = plot_df[plot_df[groupby] >= 0]
        # max_threshold = np.percentile(plot_df[groupby], 99)
        # cutted = pd.cut(plot_df[groupby][plot_df[groupby]
        #                 < max_threshold], bins=n_bins-1)
        # categories = cutted.cat.categories
        # plot_df.loc[plot_df[groupby] < max_threshold,
        #     'group_range'] = cutted.astype('str')
        # plot_df.loc[plot_df[groupby] <= categories[0].right+0.001,
        #     'group_range'] = '[0, {}]'.format(categories[0].right)

        # plot_df.loc[plot_df[groupby] > categories[-1].right - 0.001,
        #     'group_range'] = '>{:.3f}'.format(categories[-1].right)
        return plot_df, custom_sort
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
        ranges = isoform_length_ranges
        cutted = pd.cut(plot_df[groupby], ranges,include_lowest=True)
        categories = cutted.cat.categories
        plot_df.loc[:, 'group_range'] = cutted.astype(str)
        plot_df.loc[plot_df[groupby] > ranges[-1],
                    'group_range'] = '>{}'.format(ranges[-1])
        plot_df.loc[plot_df['group_range'] == str(categories[0]),'group_range'] = '[{},{}]'.format(ranges[0],categories[0].right)
        plot_df = plot_df[plot_df['group_range'] != 'nan']
        # plot_df[groupby] = plot_df[groupby].astype(int)
        # plot_df = plot_df.dropna()
        # if plot_df[groupby].max() > 3000:
        #     max_threshold = 4000
        #     lower, higher = int(plot_df[groupby].min()), 4000
        #     step_size = 400
        # else:
        #     max_threshold = 2100
        #     lower, higher = int(plot_df[groupby].min()), 2100
        #     step_size = 200
        # # # max_threshold = np.ceil(np.percentile(plot_df[groupby], 80))
        # # # lower, higher = int(plot_df.min()), int(plot_df.max())
        # # # step_size = int(math.ceil((higher - lower)/n_bins))
        # n_bins = 10
        # edges = [lower] + list(
        #     range(step_size, higher+1, step_size))
        # plot_df.loc[plot_df[groupby] <= max_threshold, 'group_range'] = pd.cut(
        #     plot_df.loc[plot_df[groupby] <= max_threshold, groupby], bins=edges).apply(lambda x: str(pd.Interval(left=int(round(x.left)), right=int(round(x.right)))) if x is not None else 'nan').astype(str)
        # plot_df.loc[plot_df[groupby] > max_threshold,
        #             'group_range'] = '>{}'.format(max_threshold)
        # plot_df = plot_df[plot_df['group_range'] != 'nan']
        return plot_df, custom_sort
    elif groupby in ['num_exons','num_isoforms']:
        def custom_sort(col):
            vals = []
            for val in col.tolist():
                if ',' in val:
                    vals.append(float(val.split(',')[1][1:-1]))
                else:
                    # vals.append(float(val[1:]))
                    vals.append(float('inf'))
            return pd.Series(vals)
        if groupby == 'num_exons':
            ranges = num_exons_range
        else:
            ranges = num_isoforms_range
        cutted = pd.cut(plot_df[groupby], ranges,include_lowest=True)
        categories = cutted.cat.categories
        plot_df.loc[:, 'group_range'] = cutted.apply(lambda x:str(x)).astype(str)
        plot_df.loc[plot_df[groupby] > ranges[-1],
                    'group_range'] = '>{}'.format(ranges[-1])
        plot_df.loc[plot_df['group_range'] == str(
            categories[0]), 'group_range'] = '[{}, {}]'.format(int(ranges[0]), int(categories[0].right))
        plot_df = plot_df[plot_df['group_range'] != 'nan']
        return plot_df, custom_sort
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
        # Only cut 80% of data and leave 20% together
        temp_df = plot_df
        temp_df[groupby] = plot_df[groupby].astype(int)
        max_threshold = np.ceil(np.percentile(plot_df[groupby], 80))
        temp_df = plot_df.loc[plot_df[groupby] <= max_threshold, groupby]
        lower, higher = temp_df.min(), temp_df.max()
        n_bins = 10
        edges = list(
            range(int(lower-1), int(higher), int(math.ceil((higher - lower)/n_bins))))
        edges.append(higher)
        plot_df.loc[plot_df[groupby] <= max_threshold, 'group_range'] = pd.cut(
            temp_df, bins=edges).astype('str')
        plot_df.loc[plot_df[groupby] > max_threshold,
                    'group_range'] = '>{}'.format(max_threshold)
        return plot_df, custom_sort
def get_density(x,y):
    from scipy.stats import gaussian_kde
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    # idx = z.argsort()
    # x, y, z = x[idx], y[idx], z[idx]
    return x,y,z


