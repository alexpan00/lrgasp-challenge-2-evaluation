from preprocess import *
import os
import plotly
import plotly.graph_objects as go
import plotly.io as pio
from plot_func.single_sample_plotter import Single_sample_plotter
from plot_func.multi_sample_plotter import Multi_sample_plotter
from plot_func.single_sample_multi_method_plotter import Single_sample_multi_method_plotter
from plot_func.multi_sample_multi_method_plotter import Multi_sample_multi_method_plotter
from static_data import color_schemes
def define_write_to_file_theme():
    # naming a layout theme for future reference
    pio.templates["encode"] = go.layout.Template(
        layout_colorway=color_schemes,
        data_scatter=[dict(line=dict(width=5))]
    )
    pio.templates["small"] = go.layout.Template(
        layout_font=dict(family="Arial Black", size=22),
        layout_title_font = dict(family="Arial Black", size=27),
    )
    pio.templates["large"] = go.layout.Template(
        layout_font=dict(family="Arial Black", size=35),
        layout_title_font = dict(family="Arial Black", size=40),
        
    )
    pio.templates["ultralarge"] = go.layout.Template(
        layout_font=dict(family="Arial Black", size=40),
        layout_title_font = dict(family="Arial Black", size=45),
    )
    pio.templates["medium"] = go.layout.Template(
        layout_font=dict(family="Arial Black", size=25),
        layout_title_font = dict(family="Arial Black", size=30),
    )

    # pio.templates.default = "encode"
    pio.templates.default = "presentation+encode+small"
def define_presentation_theme():
    # naming a layout theme for future reference
    pio.templates["encode"] = go.layout.Template(
        layout_font=dict(family="Helvetica"),
        layout_title_font = dict(family="Helvetica"),
        layout_colorway=color_schemes,
        data_scatter=[dict(line=dict(width=5))]
    )
    pio.templates["large"] = go.layout.Template(
        layout_font=dict(family="Helvetica", size=16),
        layout_title_font = dict(family="Helvetica", size=19),
    )
    pio.templates["ultralarge"] = go.layout.Template(
        layout_font=dict(family="Helvetica", size=16),
        layout_title_font = dict(family="Helvetica", size=19),
    )
    pio.templates["medium"] = go.layout.Template(
        layout_font=dict(family="Helvetica", size=16),
        layout_title_font = dict(family="Helvetica", size=19),
    )
    pio.templates.default = "presentation+encode"
def write_plot_files(fig,plot_figure_name,output_path):
    plot_path = os.path.join(output_path,'plots')
    if not os.path.exists(plot_path):
        os.makedirs(plot_path)
    pdf_path = os.path.join(plot_path,plot_figure_name+'.pdf')
    png_path = os.path.join(plot_path,plot_figure_name+'.png')
    for i in fig['layout']['annotations']:
        i['font'] = dict(size=35)
    fig['layout']['legend_title']['font']['size'] = 50
    fig.write_image(pdf_path,format='pdf')
    fig.write_image(png_path,format='png')
def get_html_str(fig,plot_figure_name,output_path):
    plot_path = os.path.join(output_path,'plots')
    if not os.path.exists(plot_path):
        os.makedirs(plot_path)
    html_path = os.path.join(plot_path,plot_figure_name+'.html')
    if (fig['layout']['width'] is not None):

        if fig['layout']['width'] <=270:
            fig['layout']['height'] = 270/fig['layout']['width'] * fig['layout']['height'] * 1.25
            fig['layout']['width'] = 270 * 1.25
        elif fig['layout']['width'] <= 960:
            fig['layout']['height'] = 960/fig['layout']['width'] * fig['layout']['height'] * 1.25
            fig['layout']['width'] = 960 * 1.25
        else:
            fig['layout']['height'] = 1920/fig['layout']['width'] * fig['layout']['height']
            fig['layout']['width'] = 1920
        
    if plot_figure_name in ['Resolution Entropy']:
        fig['layout']['height'] = 540
        fig['layout']['width'] = 540
    if plot_figure_name in ['Consistency Measure curve']:
        fig['layout']['height'] = 960/fig['layout']['width'] * fig['layout']['height'] * 1.25
        fig['layout']['width'] = 960 * 1.25
    fig['layout']['title'] = None
    html_str = fig.to_html(html_path,full_html=False,include_plotlyjs=False)
    return html_str
def generate_method_legend(method_names):
    header = ['Method name','Color']
    rows = [{'name':name,'color':color_schemes[i]} for i,name in zip(range(len(method_names)),method_names)]
    return {'header':header,'rows':rows}
def make_plots(args,output_path,is_multi_sample,is_multi_method,is_long_read,ground_truth_given,K_value_selection,sections):
    section_indices =  {sections[index]['id']:index for index in range(len(sections))}
    if (is_multi_method == False):
        df,anno_df = args
        if (is_multi_sample == False):
            plotter = Single_sample_plotter(df,anno_df)
            for plot_figure_name in single_sample_plot_figures:
                define_write_to_file_theme()
                fig = plotter.plot(plot_figure_name, 'all')
                write_plot_files(fig,plot_figure_name,output_path)
                define_presentation_theme()
                fig = plotter.plot(plot_figure_name, 'all')
                html_str = get_html_str(fig,plot_figure_name,output_path)
                section_id = single_sample_plot_figures[plot_figure_name]['type']
                sections[section_indices[section_id]]['plots'].append({'title':plot_figure_name,'html_str':html_str})
        else:
            plotter = Multi_sample_plotter(df,anno_df)
            if (ground_truth_given):
                figures = multi_sample_diff_condition_with_ground_truth_plot_figures
            else:
                figures = multi_sample_diff_condition_without_ground_truth_plot_figures
            for plot_figure_name in figures:
                define_write_to_file_theme()
                fig = plotter.plot(plot_figure_name, 'all',ground_truth_given)
                write_plot_files(fig,plot_figure_name,output_path)
                define_presentation_theme()
                fig = plotter.plot(plot_figure_name, 'all',ground_truth_given)
                html_str = get_html_str(fig,plot_figure_name,output_path)
                section_id = figures[plot_figure_name]['type']
                sections[section_indices[section_id]]['plots'].append({'title':plot_figure_name,'html_str':html_str})
    else:
        dfs,anno_df,method_names = args
        if (is_multi_sample == False):
            sections[1]['legend'] = generate_method_legend(method_names)
            plotter = Single_sample_multi_method_plotter(dfs,anno_df,method_names)
            for plot_figure_name in single_sample_plot_figures:
                define_write_to_file_theme()
                fig = plotter.plot(plot_figure_name, 'all')
                write_plot_files(fig,plot_figure_name,output_path)
                define_presentation_theme()
                fig = plotter.plot(plot_figure_name, 'all')
                html_str = get_html_str(fig,plot_figure_name,output_path)
                section_id = single_sample_plot_figures[plot_figure_name]['type']
                sections[section_indices[section_id]]['plots'].append({'title':plot_figure_name,'html_str':html_str})
        else:
            sections[1]['legend'] = generate_method_legend(method_names)
            plotter = Multi_sample_multi_method_plotter(dfs,anno_df,method_names)
            if (ground_truth_given):
                figures = multi_sample_diff_condition_with_ground_truth_plot_figures
            else:
                figures = multi_sample_diff_condition_without_ground_truth_plot_figures
            for plot_figure_name in figures:
                define_write_to_file_theme()
                fig = plotter.plot(plot_figure_name, 'all',ground_truth_given)
                write_plot_files(fig,plot_figure_name,output_path)
                define_presentation_theme()
                fig = plotter.plot(plot_figure_name, 'all',ground_truth_given)
                html_str = get_html_str(fig,plot_figure_name,output_path)
                section_id = figures[plot_figure_name]['type']
                sections[section_indices[section_id]]['plots'].append({'title':plot_figure_name,'html_str':html_str})
    return sections
    
    