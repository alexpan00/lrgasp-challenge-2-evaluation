from preprocess import load_annotation
import pickle
def save_k_dict(annot_path,save_path):
    with open(save_path,'wb') as f:
        k_val_dict = load_annotation([annot_path],True,'Condition_number')
        pickle.dump(k_val_dict,f)
path = '/fs/project/PCON0009/Au-scratch2/haoran/quantification_evaluation/example_scripts/gtf'
# for annot in ['LRGASP/lrgasp_gencode_v38_sirvs.gtf','LRGASP/lrgasp_gencode_vM27_sirvs.gtf']:
#     annot_path = '{}/{}'.format(path,annot)
#     save_file_name = '{}.pkl'.format(annot.split('/')[1].split('.')[0])
#     savepath = '{}/{}'.format('/fs/project/PCON0009/Au-scratch2/haoran/_projects/LRGASP_visualization/encode_quantification/library/k_value_dicts',save_file_name)
#     save_k_dict(annot_path,savepath)
for annot in ['ensembl/Homo_sapiens.GRCh38.104.chr.gtf']:
    annot_path = '{}/{}'.format(path,annot)
    save_file_name = '{}.pkl'.format(annot.split('/')[1].split('.')[0])
    savepath = '{}/{}'.format('/fs/project/PCON0009/Au-scratch2/haoran/_projects/LRGASP_visualization/encode_quantification/library/k_value_dicts',save_file_name)
    save_k_dict(annot_path,savepath)
