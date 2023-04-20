from .parse_annotation import parse_annotation
from .construct_feature_matrix import calculate_all_condition_number
from collections import defaultdict
import time
import numpy as np
from numpy import linalg as LA

import time
import os
def check_region_type(region_name):
    if ((region_name.count(':') == 1) and ('-' not in region_name)):
        return 'one_exon'
    elif ((region_name.count(':') == 2) and ('-' not in region_name)):
        return 'two_exons'
    elif (region_name.count('-') == 1):
        return 'one_junction'
    elif ((region_name.count(':') > 2) and ('-' not in region_name)):
        return 'exons'
    else:
        return 'junctions'
def calculate_exon_min_read_mapped_length(exon_region_name,point_dict,exon_position):
    assert '-' not in exon_region_name
    points = exon_region_name.split(':')
    if len(points) == 2:
        return 1
    if exon_position == 'left':
        return point_dict[points[-1]] - point_dict[points[1]] + 1 + 1
    elif exon_position == 'right':
        return point_dict[points[-2]] - point_dict[points[0]] + 1 + 1
    elif exon_position == 'center':
        return point_dict[points[-2]] - point_dict[points[1]] + 2 + 1
def filter_regions_num_exons(gene_regions_dict,genes_regions_len_dict):
    new_gene_regions_dict = defaultdict(lambda:defaultdict(dict))
    new_genes_regions_len_dict = defaultdict(lambda:defaultdict(dict))
    for rname in gene_regions_dict:
        for gname in gene_regions_dict[rname]:
            regions_set = set()
            for region in gene_regions_dict[rname][gname]:
                if check_region_type(region) in ['one_exon','one_junction']:
                    regions_set.add(region)
            for region_name in regions_set:
                new_gene_regions_dict[rname][gname][region_name] = gene_regions_dict[rname][gname][region_name]
                new_genes_regions_len_dict[rname][gname][region_name] = genes_regions_len_dict[rname][gname][region_name]
    return new_gene_regions_dict,new_genes_regions_len_dict
def check_valid_region(chr_name,gene_name,region_name,genes_regions_len_dict,max_read_len,min_read_len,point_dict,READ_JUNC_MIN_MAP_LEN):
    is_region_valid = False
    if (genes_regions_len_dict[chr_name][gene_name][region_name] > min_read_len):              
        if (max_read_len is not None):
            if check_region_type(region_name) in ['one_exon','two_exons']:
                is_region_valid = True
            elif check_region_type(region_name) == 'one_junction':
                exons = region_name.split('-')
                for i,exon in zip(range(len(exons)),exons):
                    points = exon.split(':')
                    if (i == 0):
                        first_exon_length = point_dict[points[-1]] - point_dict[points[0]] + 1
                    elif i == len(exons) - 1:
                        last_exon_length = point_dict[points[-1]] - point_dict[points[0]] + 1
                if (first_exon_length > READ_JUNC_MIN_MAP_LEN) and (last_exon_length > READ_JUNC_MIN_MAP_LEN):
                    if region_name.count(':') == 2:
                        is_region_valid = True
                    else:
                        [exon_1,exon_2] = region_name.split('-')
                        if 2 * READ_JUNC_MIN_MAP_LEN - 2 + calculate_exon_min_read_mapped_length(exon_1,point_dict,exon_position='left') + calculate_exon_min_read_mapped_length(exon_2,point_dict,exon_position='right') <= max_read_len:
                            is_region_valid = True
            elif check_region_type(region_name) == 'exons':
                if calculate_exon_min_read_mapped_length(region_name,point_dict,exon_position='center') <= max_read_len:
                    is_region_valid = True
            elif check_region_type(region_name) == 'junctions':
                exons = region_name.split('-')
                exon_length = 0
                for i,exon in zip(range(len(exons)),exons):
                    points = exon.split(':')
                    if (i == 0):
                        first_exon_length = point_dict[points[-1]] - point_dict[points[0]] + 1
                    elif i == len(exons) - 1:
                        last_exon_length = point_dict[points[-1]] - point_dict[points[0]] + 1
                    else:
                        exon_length += point_dict[points[-1]] - point_dict[points[0]] + 1
                if (first_exon_length > READ_JUNC_MIN_MAP_LEN) and (last_exon_length > READ_JUNC_MIN_MAP_LEN):
                    if exon_length + 2 * READ_JUNC_MIN_MAP_LEN - 2 + calculate_exon_min_read_mapped_length(exons[0],point_dict,exon_position='left') + calculate_exon_min_read_mapped_length(exons[-1],point_dict,exon_position='right') <= max_read_len:
                        is_region_valid = True
        else:
            is_region_valid = True
    return is_region_valid

def filter_regions_read_length(gene_regions_dict,gene_points_dict,genes_regions_len_dict,READ_JUNC_MIN_MAP_LEN,min_read_len,max_read_len=None):
    new_gene_regions_dict = defaultdict(lambda:defaultdict(dict))
    new_genes_regions_len_dict = defaultdict(lambda:defaultdict(dict))
    for chr_name in gene_regions_dict:
        for gene_name in gene_regions_dict[chr_name]:
            point_dict = {}
            for coord in gene_points_dict[chr_name][gene_name]:
                point_dict['P{}'.format(gene_points_dict[chr_name][gene_name][coord])] = int(coord)
            for region_name in gene_regions_dict[chr_name][gene_name]:
                is_region_valid = check_valid_region(chr_name,gene_name,region_name,genes_regions_len_dict,max_read_len,min_read_len,point_dict,READ_JUNC_MIN_MAP_LEN)
                if (is_region_valid):
                    new_gene_regions_dict[chr_name][gene_name][region_name] = gene_regions_dict[chr_name][gene_name][region_name]
                    new_genes_regions_len_dict[chr_name][gene_name][region_name] = genes_regions_len_dict[chr_name][gene_name][region_name]
    return new_gene_regions_dict,new_genes_regions_len_dict
def filter_long_read_regions(gene_regions_dict,genes_regions_len_dict):
    new_gene_regions_dict = defaultdict(lambda:defaultdict(dict))
    new_genes_regions_len_dict = defaultdict(lambda:defaultdict(dict))
    for rname in gene_regions_dict:
        for gname in gene_regions_dict[rname]:
            regions_set = set()
            isoform_region_dict = defaultdict(lambda:set())
            for region in gene_regions_dict[rname][gname]:
                for isoform in gene_regions_dict[rname][gname][region]:
                    isoform_region_dict[isoform].add(region)
            for isoform in isoform_region_dict:
                max_region_exon_num = 0
                longest_region = ''
                for region in isoform_region_dict[isoform]:
                    region_exon_num = region.count(':')
                    if max_region_exon_num < region_exon_num:
                        max_region_exon_num = region_exon_num
                        longest_region = region
                if max_region_exon_num <= 1:
                    if longest_region != '':
                        regions_set.add(longest_region)
                else:
                    for region in isoform_region_dict[isoform]:
                        region_exon_num = region.count(':')
                        if (region_exon_num >= max_region_exon_num - 1):
                            regions_set.add(region)
            for region_name in regions_set:
                new_gene_regions_dict[rname][gname][region_name] = gene_regions_dict[rname][gname][region_name]
                new_genes_regions_len_dict[rname][gname][region_name] = genes_regions_len_dict[rname][gname][region_name]
    return new_gene_regions_dict,new_genes_regions_len_dict
                                
                        


def parse_reference_annotation(ref_file_path,threads,READ_LEN,READ_JUNC_MIN_MAP_LEN,is_long_read):
    [gene_exons_dict, gene_points_dict, gene_isoforms_dict, genes_regions_len_dict,
        _, gene_regions_dict, gene_isoforms_length_dict,raw_isoform_exons_dict] = parse_annotation(ref_file_path, threads,READ_LEN, READ_JUNC_MIN_MAP_LEN)
    if (is_long_read):
        gene_regions_dict,genes_regions_len_dict = gene_regions_dict,genes_regions_len_dict
    else:
        # gene_regions_dict,genes_regions_len_dict = filter_regions_num_exons(gene_regions_dict,genes_regions_len_dict)
        gene_regions_dict,genes_regions_len_dict = filter_regions_num_exons(gene_regions_dict,genes_regions_len_dict)
    # for chr_name in gene_regions_dict:
    #     for gene_name in gene_regions_dict[chr_name].copy():
    #         if (len(gene_regions_dict[chr_name][gene_name]) == 0):
    #             for dic in [gene_points_dict,gene_isoforms_dict,gene_regions_dict,genes_regions_len_dict,gene_isoforms_length_dict,raw_isoform_exons_dict]:
    #                 if chr_name in dic and gene_name in dic[chr_name]:
    #                     del dic[chr_name][gene_name]
    num_isoforms = 0
    for rname in gene_isoforms_dict:
        for gname in gene_isoforms_dict[rname]:
            num_isoforms += len(gene_isoforms_dict[rname][gname])
    print(num_isoforms)
    return gene_points_dict,gene_isoforms_dict,gene_regions_dict,genes_regions_len_dict,gene_isoforms_length_dict,raw_isoform_exons_dict
def get_kvalues_dict(gtf_io,is_long_read,K_value_selection,READ_LEN=150,READ_JUNC_MIN_MAP_LEN=0):
    # threads = len(os.sched_getaffinity(0))
    threads = 1
    print('Start calculating K values...')
    start = time.time()
    _,gene_isoforms_dict,gene_regions_dict,_,gene_isoforms_length_dict,raw_isoform_exons_dict = parse_reference_annotation(gtf_io,threads,READ_LEN,READ_JUNC_MIN_MAP_LEN,is_long_read)
    print('Parsing annoation done in {}s !'.format(time.time()-start))
    start = time.time()
    gene_matrix_dict = calculate_all_condition_number(gene_isoforms_dict,gene_regions_dict,allow_multi_exons=is_long_read)
    kvalues_dict,num_exon_dict,isoform_length_dict,isoform_gene_dict,num_isoforms_dict = {},{},{},{},{}
    for chr_name in gene_matrix_dict:
        for gene_name in gene_matrix_dict[chr_name]:
            # if (len(raw_isoform_exons_dict[chr_name][gene_name])>1):
            for isoform_name in raw_isoform_exons_dict[chr_name][gene_name]:
                K_value_choice = 0
                if (K_value_selection=='Condition_number'):
                    K_value_choice = 1
                elif (K_value_selection=='K_value'):
                    K_value_choice = 0
                elif (K_value_selection=='Generalized_condition_number'):
                    K_value_choice = 2
                kvalues_dict[isoform_name] = gene_matrix_dict[chr_name][gene_name]['condition_number'][K_value_choice]
                num_exon_dict[isoform_name] = len(raw_isoform_exons_dict[chr_name][gene_name][isoform_name]['start_pos'])
                isoform_length_dict[isoform_name] = gene_isoforms_length_dict[chr_name][gene_name][isoform_name]
                num_isoforms_dict[isoform_name] =  len(raw_isoform_exons_dict[chr_name][gene_name])
                isoform_gene_dict[isoform_name] = gene_name

    print('Calculating K values done in {}s !'.format(time.time()-start))
    return kvalues_dict,num_exon_dict,isoform_length_dict,isoform_gene_dict,num_isoforms_dict

