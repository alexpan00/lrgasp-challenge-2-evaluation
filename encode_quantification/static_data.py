# X_options = [
#     {'label': 'Normalized Root Mean Square Error', 'value': 'NRMSE'},
#     {'label': 'Median Relative Difference', 'value': 'MRD'},
#     {'label': "Spearman's rho", 'value': 'rho'}
# ]
# Y_options = [
#     {'label': 'Distribution','value': 'dist'}
# ]
color_schemes =["#C0504D",
                "#4F81BD",
            
			"#FF8C00",
            "#808000",
            "#0000FF",
			"#339966",
			"#E88471",	
            "#009392",
            "#333333",
            "#045275",]
fig_size = {
    'small_rec':{'width':270,'height':540},
    'rec':{'width':960,'height':540},
    'square':{'width':960,'height':960},
    'small_square':{'width':540,'height':540},
    'large_square':{'width':1080,'height':1080},
    'large_rec':{'width':1920,'height':1080},
}
themes = {
    'small_single':'presentation+encode',
    'medium_single':'medium+encode',
    'large_single':'large+encode',
    'large_multi':'ultralarge+encode',
    'medium_multi':'large+encode',
    'small_multi':'medium+encode'
}
on_plot_shown_label = {
    'CM':'CM',
    'nrmse': 'NRMSE',
    'mrd': 'MRD',
    'spearmanr':"Spearman",
    'precision': 'Precision',
    'recall': 'Recall',
    'accuracy':'Accuracy',
    'f1': 'F1 Score',
    'auc': 'AUC',
    'fpr': 'False Positive Rate',
    'tpr': 'True Positive Rate',
    'isoform_length': 'Isoform length',
    'true_abund': 'True abundance',
    'log2_true_abund': 'Log2(True abundance+1)',
    'std': 'std',
    'estimated_abund': 'Estimated abundance',
    'RM': 'RM',
    'RE': 'RE',
    'ave_estimated_abund':'Log2(Esimated abundance+1)',
    'ave_true_abund':'Log2(True abundance+1)',
    'num_exons': '# of exons',
    'K_value': 'K value',
    'mean_arr':'Mean ARR'
}
# on_plot_shown_label = {
#     'CM':'Consistency Measure',
#     'nrmse': 'Normalized Root Mean Square Error',
#     'mrd': 'Median Relative Difference',
#     'spearmanr':"Spearman's rho",
#     'precision': 'Precision',
#     'recall': 'Recall',
#     'accuracy':'Accuracy',
#     'f1': 'F1 Score',
#     'auc': 'AUC',
#     'fpr': 'False Positive Rate',
#     'tpr': 'True Positive Rate',
#     'isoform_length': 'Isoform length',
#     'true_abund': 'True abundance',
#     'log2_true_abund': 'Log2(True abundance+1)',
#     'std': 'std',
#     'estimated_abund': 'Estimated abundance',
#     'RM': 'Reproducibility Measure',
#     'RE': 'Resolution Entropy',
#     'ave_estimated_abund':'Log2(Estimated abundance+1)',
#     'ave_true_abund':'Log2(true_abund+1)',
#     'num_exons': 'Number of exons',
#     'K_value': 'K value'
# }
data_sample_options = [
    {'label': 'Single sample', 'value': 'single_sample'},
    {'label': 'Multiple sample under different condition', 'value': 'multi_sample_diff_condition'}]
        # {'label': 'Multiple sample under same condition', 'value': 'multi_sample_same_condition'}
multi_methods_options = [
    {'label': 'Single method', 'value': 'single_method'},
    {'label': 'Multiple method', 'value': 'multi_method'}]
single_sample_table_metrics = [
    {'name': 'Normalized Root Mean Square Error', 'id': 'nrmse','type':'numeric','format':{'specifier':'.3f'},},
    {'name': 'Median Relative Difference', 'id': 'mrd','type':'numeric','format':{'specifier':'.3f'}},
    {'name':'Mean Abundance Recovery Rate','id':'mean_arr','type':'numeric','format':{'specifier':'.1%'}},
    {'name': "Spearman's rho", 'id': 'spearmanr','type':'numeric','format':{'specifier':'.3f'}},
    {'name': "Resolution Entropy", 'id': 'RE','type':'numeric','format':{'specifier':'.3f'}},
]
multi_sample_diff_condition_table_metrics = [
    {'name': 'Normalized Root Mean Square Error', 'id': 'nrmse','type':'numeric','format':{'specifier':'.3f'}},
    {'name': 'Median Relative Difference', 'id': 'mrd','type':'numeric','format':{'specifier':'.3f'}},
    {'name':'Mean Abundance Recovery Rate','id':'mean_arr','type':'numeric','format':{'specifier':'.1%'}},
    {'name': "Spearman's rho", 'id': 'spearmanr','type':'numeric','format':{'specifier':'.3f'}},

    {'name': "Resolution Entropy", 'id': 'RE','type':'numeric','format':{'specifier':'.3f'}},

    {'name': 'Consistency Measure', 'id': 'CM','type':'numeric','format':{'specifier':'.3f'}},
    
    {'name': 'Reproducibility Measure', 'id': 'RM','type':'numeric','format':{'specifier':'.3f'}},

    {'name': 'Precision','id': 'precision','type':'numeric','format':{'specifier':'.3f'}},
    {'name': 'Recall','id': 'recall','type':'numeric','format':{'specifier':'.3f'}},
    {'name': 'Accuracy','id': 'accuracy','type':'numeric','format':{'specifier':'.3f'}},
    {'name': 'F1 Score','id': 'f1','type':'numeric','format':{'specifier':'.3f'}},
    {'name': 'AUC','id': 'auc','type':'numeric','format':{'specifier':'.3f'}},
]
multi_sample_diff_condition_without_ground_truth_table_metrics = [
    {'name': 'Consistency Measure', 'id': 'CM','type':'numeric','format':{'specifier':'.3f'}},
    {'name': 'Reproducibility Measure', 'id': 'RM','type':'numeric','format':{'specifier':'.3f'}},
    {'name': "Resolution Entropy", 'id': 'RE','type':'numeric','format':{'specifier':'.3f'}},
]

single_sample_plot_figures = {
    'Distribution of K values':{'x':'K_value','y':'dist','type':'gene_features'},
    'Distribution of isoform lengths':{'x':'isoform_length','y':'dist','type':'gene_features'},
    'Distribution of numbers of exons':{'x':'num_exons','y':'dist','type':'gene_features'},
    'Histogram of Abundance Recovery Rate':{'x':'arr','y':'freq','type':'estimation_error'},
    'Correlation of estimated abundance and ground truth':{'x':'true_abund','y':'estimated_abund','type':'estimation_error'},
    'Correlation Boxplot of estimated abundance and ground truth':{'x':'true_abund','y':'estimated_abund','type':'estimation_error'},
    'Resolution Entropy':{'x':'','y':'RE','type':'resolution_entropy'},
    'Statistics with different K values':{'x':'K_value','y':[m['id'] for m in single_sample_table_metrics],'type':'statistics'},
    'Statistics with different isoform lengths':{'x':'isoform_length','y':[m['id'] for m in single_sample_table_metrics],'type':'statistics'},
    'Statistics with different numbers of exons':{'x':'num_exons','y':[m['id'] for m in single_sample_table_metrics],'type':'statistics'},
    'Statistics with different expression level':{'x':'log2_true_abund','y':[m['id'] for m in single_sample_table_metrics],'type':'statistics'},
    # 'Standard deviation vs estimated abundance scatter':{'x':'estimated_abund','y':'std','type':'estimation_error'},
}
multi_sample_diff_condition_with_ground_truth_plot_figures = {
    'Distribution of K values':{'x':'K_value','y':'dist','type':'gene_features'},
    'Distribution of isoform lengths':{'x':'isoform_length','y':'dist','type':'gene_features'},
    'Distribution of numbers of exons':{'x':'num_exons','y':'dist','type':'gene_features'},

    'Estimation Error for different conditions':{'x':'metrics','y':['nrmse','mrd','spearmanr'],'type':'estimation_error'},
    'Histogram of Abundance Recovery Rate':{'x':'arr','y':'freq','type':'estimation_error'},
    'Correlation of estimated abundance and ground truth':{'x':['ave_true_abund#1','ave_true_abund#2'],'y':['ave_estimated_abund#1','ave_estimated_abund#2'],'type':'estimation_error'},
    'Correlation Boxplot of estimated abundance and ground truth':{'x':['ave_true_abund#1','ave_true_abund#2'],'y':['ave_estimated_abund#1','ave_estimated_abund#2'],'type':'estimation_error'},
    
    'Resolution Entropy for different conditions':{'x':'','y':'RE','type':'resolution_entropy'},
    
    'Consistency Measure curve':{'x':'C threshold','y':'CM','type':'consistency'},
    # 'Standard deviation vs estimated abundance scatter':{'x':['ave_estimated_abund#1','ave_estimated_abund#2'],'y':['std#1','std#2'],'type':'reproducibility'},
    'Standard deviation vs estimated abundance curve':{'x':['ave_estimated_abund#1','ave_estimated_abund#2'],'y':['std#1','std#2'],'type':'reproducibility'},
    'ROC curves for performance of quantification':{'x':'fpr','y':'tpr','type':'fold_change'},
    'PR curves for performance of quantification':{'x':'recall','y':'precision','type':'fold_change'},

    'Statistics with different K values':{'x':'K_value','y':[m['id'] for m in multi_sample_diff_condition_table_metrics],'type':'statistics'},
    'Statistics with different isoform lengths':{'x':'isoform_length','y':[m['id'] for m in multi_sample_diff_condition_table_metrics],'type':'statistics'},
    'Statistics with different numbers of exons':{'x':'num_exons','y':[m['id'] for m in multi_sample_diff_condition_table_metrics],'type':'statistics'},
    'Statistics with different expression level':{'x':'ave_true_abund','y':[m['id'] for m in multi_sample_diff_condition_table_metrics],'type':'statistics'},
    
    
    
}
multi_sample_diff_condition_without_ground_truth_plot_figures = {
    'Distribution of K values':{'x':'K_value','y':'dist','type':'gene_features'},
    'Distribution of isoform lengths':{'x':'isoform_length','y':'dist','type':'gene_features'},
    'Distribution of numbers of exons':{'x':'num_exons','y':'dist','type':'gene_features'},
    'Standard deviation vs estimated abundance curve':{'x':['ave_estimated_abund#1','ave_estimated_abund#2'],'y':['std#1','std#2'],'type':'reproducibility'},
    'Resolution Entropy for different conditions':{'x':'','y':'RE','type':'resolution_entropy'},
    'Statistics with different K values':{'x':'K_value','y':[m['id'] for m in multi_sample_diff_condition_without_ground_truth_table_metrics],'type':'statistics'},
    'Statistics with different isoform lengths':{'x':'isoform_length','y':[m['id'] for m in multi_sample_diff_condition_without_ground_truth_table_metrics],'type':'statistics'},
    'Statistics with different numbers of exons':{'x':'num_exons','y':[m['id'] for m in multi_sample_diff_condition_without_ground_truth_table_metrics],'type':'statistics'},
}
transform_options = [
    {'label': 'Linear', 'value': 'linear'},
    {'label': 'Log', 'value': 'log'}]

filter_options = [
    {'label': 'No filter', 'value': 'all'},
    {'label': 'K value >= 1', 'value': 'k>=1'},
    {'label': 'K value < 1', 'value': 'k<1'},
]

annot_options = [
    {'label':'lrgasp_gencode_v38_sirvs(human)','value':'human'},
    {'label':'lrgasp_gencode_vM27_sirvs(mouse)','value':'mouse'},
    {'label':'Ensembl_Homo_sapiens.GRCh38.104.chr(human)','value':'ensembl_human'},
]
# K_value_ranges = [i/10 for i in range(11)] + [i*2.5 for i in range(1,11)]
K_value_ranges = [i/10 for i in range(11)]
condition_number_ranges = [i*2.5 for i in range(10)] + [i*25 for i in range(1,21,2)]
generalized_condition_number_ranges = [i for i in range(11)] + [i*25 for i in range(1,11)]
ARR_ranges = [i/10 for i in range(11)]