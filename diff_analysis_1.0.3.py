#!/usr/bin/env python3
# -*- coding:utf-8 -*-


"""A short description of the module -- called a docstring."""


# =====================================================================================================================
# imports
# =====================================================================================================================
import sys
import os

main_description = '''
diff_analysis.py
Differential analysis
Author: Vincent Hu 

This wrapper will prepare R scripts to perform differential analysis GLM test on each factor
listed in meta table, as well as on all/selected interaction terms. Traditional pair-wise tests between samples groups
for post-hoc, and run the tests in R. FDR correction for each set of p-value(s) will be calculated.

For normal models, standard linear models in R will be used in analysis

Python 3, R, Rscript, edgeR (R library) are required (in your PATH).

You can generate the R script only without running (for debug) by setting the -nx flag.

Example: diff_analysis.py -c counts.txt -m meta.txt -o CRI_123 -nx -pca -pcal 

Suggestion: Please load anaconda environment to make sure you have Python 3 in your path
\"export $PATH:/group/bioinformatics/software/anaconda3/bin\" # If you have not set up Anaconda environment

'''

if os.environ['CONDA_DEFAULT_ENV'] != 'python3.7':
    raise Exception('Please load Anaconda environment by \"source activate python3.7\"')

import os.path
import argparse
# import re
import itertools
# import xlsxwriter
import numpy as np
import pandas as pd
# import pandas.io.formats.excel # TURN THIS ON if use pandas 0.20.1 or later is used !!!!!!!!!!!!!!!!
# import logging
# import collections
import subprocess as sp
import modules
import txt2xlsx


# =====================================================================================================================
# Global variables
# =====================================================================================================================
anno_table_file = '~/REF/annotation/RNAseq/compiled_annotation_list.txt'
anno_dic = {}
# =====================================================================================================================
# Class definitions
# =====================================================================================================================


# =====================================================================================================================
# Function definitions
# =====================================================================================================================

# ===============================================
#
# ===============================================
def is_valid_table(arg):
    # print(arg)
    abs_path = os.path.abspath(arg)
    if not os.path.isfile(abs_path):
        msg = 'Error: the input table ' + arg + ' can not be found!'
        raise argparse.ArgumentTypeError(msg)
        # sys.exit(1)
    return abs_path


# ===============================================
#
# ===============================================
def is_valid_meta_table(arg):
    print(arg)
    abs_path = os.path.abspath(arg)
    if os.path.isfile(abs_path):
        return abs_path
    else:
        raise argparse.ArgumentError('Error: the input meta table ' + arg + ' can not be found!')
        # sys.exit(1)


# ===============================================
#
# ===============================================
def check_non_neg_int(arg):
    print(arg)
    value = int(arg)
    if value < 0:
        raise argparse.ArgumentTypeError('Error: ' + arg + ' is an invalid negative value!')
        # sys.exit(1)
    return value


# ===============================================
#
# ===============================================
def check_pos_int(arg):
    print(arg)
    value = int(arg)
    if value <= 0:
        raise argparse.ArgumentTypeError('Error: ' + arg + ' is an invalid non-positive value!')
        # sys.exit(1)
    return value


# ===============================================
#
# ===============================================
def check_float(arg):
    print(arg)
    try:
        value = float(arg)
    except ValueError:
        raise argparse.ArgumentTypeError('Error: ' + arg + ' is NOT a valid float number!')
        # sys.exit(1)
    return value


# ===============================================
#
# ===============================================
def check_valid_single_factors_list(arg):
    print(arg)
    if not arg.isdigit():
        raise ValueError('Error: the factors should be specified as an the position indices (start from 1)!')
        # sys.exit(1)
    elif int(arg) < 1:
        raise ValueError('Error: invalid factors specified! Should start from 1.')
    return int(arg)


# ===============================================
#
# ===============================================
def check_valid_interaction_list(arg):  # Expected format for arg: ['1:2:3', '2:3:4']
    if len(arg) == 0:
        raise ValueError('Error: empty interaction list')
        # sys.exit(1)
    print(arg)
    arg_split = arg.strip().split(':')
    s_sort = []
    if len(arg_split) < 2:
        raise ValueError('Error: invalid interaction term(s) in command line')
    for t in arg_split:
        if not t.isdigit():  # Make sure each factor is a number
            raise ValueError('Error: invalid interaction term(s) in command line')
        s_sort.append(int(t))
    if len(s_sort) == 0:
        raise ValueError('Error: invalid interaction term(s) in command line!')
    if len(np.unique(s_sort)) < len(s_sort):
        raise ValueError('Error: duplicated factors in interaction term(s) in command line!')
    s_sort.sort()
    return ':'.join(str(s) for s in s_sort)


# ===============================================
# Process the command line input
# ===============================================
def getargs():
    # global variables
    global anno_dic
    global main_description

    epi_log = '''
    This program requires the python packages:
    1. argparse
    2. subprocess
    2. pandas
    3. numpy
    '''

    help_help = '''
This wrapper will prepare edgeR commands for a omnibus differential analysis test for each factor
listed in -m, including all/selected interaction terms, plus all pair-wise tests between groups
for post-hoc, and run them in R. FDR correction for each set of p-value(s) will be calculated.

For normal models, standard linear models in R will be used in analysis

Python 3, R, Rscript, edgeR (R library) are required (in your PATH).    
    '''

    version_help = '''
    v1.0.3 
    '''

    # command_file_help = '''\
    # pipeline command file (e.g. pipeline_commands.sh).
    # Steps are separated by an empty line.
    # Command lines within each step will be run simultaneously.
    # '''

    source_dir_help = '''
    The source dictionary that contains the .
    '''

    recursive_help = '''
    recursive_help
    '''

    file_type_help = '''
    file_type_help
    '''

    out_dir_help = '''\
    out_dir_help
    default is None.
    '''

    hard_copy_help = '''\
    hard_copy_help
    '''

    mapping_file_help = '''\
    mapping_file_help
    default is 1.
    '''

    regular_expression_help = '''\
    regular_expression_help
    default is "KEY".
    '''

    anno_help = '''The gene annotation compiled from GTF files. 
This file will be used to annotate the normalized expression and differential analysis tabs in the output excel file
Choices include:
'''
    for s in list(anno_dic.keys()):
        anno_help += '\t' + s + '\t' + anno_dic[s] + '\n'

    usage_help = '''diff_analysis.py -c INPUT_COUNT_TABLE
                        -m INPUT_META_TABLE
                        -o OUTPUT_PREFIX
                        [-mo {nb,normal}] [-nx]
                        [-f FILTER_BY_GENE_SUM] [-F FILTER_BY_SAMPLE_SUM]
                        [-s INPUT_SPIKEIN] [-nne] [-nnd] [-l] [-ng] [-np]
                        [-ni] [-bf [BATCH_FACTORS [BATCH_FACTORS ...]]]
                        [-cf [CONT_FACTORS [CONT_FACTORS ...]]]
                        [-ef [EXCLUDE_FACTORS [EXCLUDE_FACTORS ...]]]
                        [-ei [EXCLUDE_INTERACTIONS [EXCLUDE_INTERACTIONS ...]]]
                        [-ii [INCLUDE_INTERACTIONS [INCLUDE_INTERACTIONS ...]]]
                        [-rd] [-sb SET_BCV] [-dm {single,multiple}] [-ow]
                        [-pca] [-pcal] [-pbcv] [-pbox] [-h] [-v]
    '''
    arg_parser = argparse.ArgumentParser(description=main_description, epilog=epi_log, formatter_class=argparse.RawTextHelpFormatter, add_help=False, usage=usage_help)

    ##############################
    #    1. positional arguments #
    ##############################
    # arg_parser.add_argument(dest='input_file_list', default=[], nargs='*',
    #                         help='Input tab-delimited text files or path containing tab-delimited text files to be merged into the output xlsx file. \n'
    #                              'Multiple files/paths can be inputted at the same time. \n'
    #                              'Use *(asterisk) as wildcard.',
    #                         type=lambda _x_: is_valid_input_file(_x_))

    ###############################
    #    2. required arguments    #
    ###############################
    required_group = arg_parser.add_argument_group("Required arguments")

    # required_group.add_argument("-s", dest="source_dir", action="store", required=True, default=None, help=source_dir_help)

    # required_group.add_argument('-i', dest='input_file_list', action='append', required=True, default=[],
    #                            help='Input txt files or path containing txt files to be merged into the output xlsx file',
    #                            type=lambda _x_: is_valid_input_file(_x_))
    required_group.add_argument('-c', dest='input_count_table', action='store', required=True, default=None,
                                help='Input count table, must be tab-delimited text file',
                                type=lambda _x_: is_valid_table(_x_))
    required_group.add_argument('-m', dest='input_meta_table', action='store', required=True, default=None,
                                help='Input meta table, must be tab-delimited text file',
                                type=lambda _x_: is_valid_table(_x_))
    required_group.add_argument('-o', dest='output_prefix', action='store', required=True, default=None,
                                help='''Output prefix. Output files will written to the folder where the INPUT_COUNT_TABLE is, including:
        _count_cleaned.txt\t\tFiltered count table.
                          \t\tThe genes/transcripts with sum of counts lower than -f are removed.
                          \t\tThe samples with sum of counts lower than -F were removed.
        _meta_cleaned.txt\t\tFiltered meta table. The factors/terms with one levels were removed, after sample filtering in count table cleaning.
        _spikein_lib_cleaned.txt\tThe spike-in library size table with filtered samples removed, as in count table cleaning.
        _edger_script.r\t\t\tedgeR script for running differential analysis
        _norm.txt\t\t\tNormalized expression table
        _diff.txt\t\t\tDifferential analysis table''')
    ###############################
    #    3. optional arguments    #
    ###############################
    optional_group = arg_parser.add_argument_group("Optional arguments")
    optional_group.add_argument('-mo', '--model', dest='model', action='store',
                                default='nb', choices=['nb', 'normal'],
                                help='''The underlying distribution of read counts. Default negative binomial: \"-mo nb\"
If use normal distribution, specify \"-mo normal\"''')
    optional_group.add_argument('-nx', '--NoeXecution', dest='do_exec', action='store_false',
                                help='If flag -nx is set, this script will generate the R script only, without any other results')
    optional_group.add_argument('-f', '--filterGenes', dest='filter_by_gene_sum', action='store',
                                default=1, type=check_non_neg_int,
                                help='Filtering genes if sum of raw count across all samples is less than this threshold.'
                                     'Default: \"-f 1\"')
    optional_group.add_argument('-F', '--filterSamples', dest='filter_by_sample_sum', action="store",
                                default=1, type=check_non_neg_int,
                                help='Filtering  the samples if sum of raw counts across all genes is less than this threshold.'
                                     'Default: \"-F 1\"')
    optional_group.add_argument('-s', '--spikein', dest='input_spikein', action='store', default=None,
                                help='Spikein lib size table',
                                type=lambda _x_: is_valid_table(_x_))
    optional_group.add_argument('-nne', '--NoNormExprs', dest='do_norm_exprs', action='store_false',
                                help='''Exclude normalization factor in expression.
This flag only works when negative binomial model is used (\"-mo nb\")
This flag will be neglected in normal model (\"-m normal\")''')
    optional_group.add_argument('-nnd', '--NoNormDiff', dest='do_norm_diff', action='store_false',
                                help='''Exclude normalization factor in differential analysis.
This flag will only work when negative binomial model is used (\"-mo nb\")and will be neglected in normam model \"-mo normal\"''')
    optional_group.add_argument('-l', '--linear', dest='do_logcpm', action='store_false',
                                help='''Output the linear normalized expression in _norm.txt.
By default log 2 transformed normalized expression will be written to _norm.txt.
A default value of 0.05 is added before logarithm.''')
    optional_group.add_argument('-ng', '--NoGlm', dest='do_glm', action='store_false',
                                help='Flag for excluding omnibus test(s) from analysis')
    optional_group.add_argument('-np', '--NoPosthoc', dest='do_post', action='store_false',
                                help='Flag for excluding post-hoc tests from analysis')
    optional_group.add_argument('-ni', '--NoInteractions', dest='do_inter', action='store_false',
                                help='Flag for excluding all interaction terms from modeling')
    optional_group.add_argument('-bf', '--BatchFactors', dest='batch_factors', action='store', default=[],
                                type=check_valid_single_factors_list, nargs='*',
                                help='''List of batch effect factors, starting from 1, as the 2nd column in meta table.
Example: -bf 2 3 5''')
    optional_group.add_argument('-cf', '--ContinuousFactors', dest='cont_factors', action='store', default=[],
                                type=check_valid_single_factors_list, nargs='*',
                                help='''List of continuous factors, starting from 1, as the 2nd column in meta table.
Example: -cf 2 3 5''')
    optional_group.add_argument('-ef', '--ExcludeFactors', dest='exclude_factors', default=[],
                                type=check_valid_single_factors_list, nargs='*',
                                help='''List of factors to be excluded from single factor test. Use the same format as \"-b\"
Example: -ef 2 3''')
    optional_group.add_argument('-ei', '--ExcludeInteractions', dest='exclude_interactions', default=[], nargs='*',
                                type=check_valid_interaction_list,
                                help='''List of interaction terms to be excluded in the full model.
The factors in each interaction term are separated by \":\"
Note: This flag will be override by the -ii flag.
Example: -ei 1:2:3 2:3:4''')
    optional_group.add_argument('-ii', '--IncludeInteractions', dest='include_interactions', default=[], nargs='*',
                                type=check_valid_interaction_list,
                                help='''List of interaction terms to be included in the model.
The factors in each interaction term are separated by \":\"
Note: This flag will be override the default full-rank interaction terms among all factors given in input meta file
Example: -ii 1:2:3 4:5:6''')
    optional_group.add_argument('-rd', '--RecalcDisp', dest='re_disp', action='store_true',
                                help='Re-estimate dispersion for each pairwise post-hoc test')
    optional_group.add_argument('-sb', '--SetBCV', dest='set_bcv', action='store', type=check_float,
                                help='''Arbitrarily set the Biological Coefficient of Variation(BCV), in un-replicated design. BCV is the square root of dispersion!!!
Typical Value:\tHuman:                                 0.4
              \tMouse:                                 0.2
              \tGenetically identical model organisms: 0.1
              \tTechnical replicates:                  0.01''')
    optional_group.add_argument('-dm', '--DispMethod', dest='disp_method', action='store', default='single',
                                choices=['single', 'multiple'],
                                help='Method for calculate dispersion in edgeR. Default "\"-dm single\"')
    optional_group.add_argument('-ow', '--OverWrite', dest='overwrite', action='store_true',
                                help='Overwrite the existing output files.')
    optional_group.add_argument('-pca', '--PCA', dest='do_pca', action='store_true',
                                help='Making PCA plot for samples')
    optional_group.add_argument('-pcal', '--PCALabel', dest='do_pca_label', action='store_true',
                                help='''Label the samples in PCA plot
Note: if -pca is NOT set, -pcal flag will be neglected''')
    optional_group.add_argument('-pcasub', '--PCASubtitle', dest='pca_subtitle', action='store', default='',
                                help='''Subtitle for the PCA plot''')
    optional_group.add_argument('-pbcv', '--PlotBCV', dest='do_bcv', action='store_true',
                                help='''Make dispersion plot
Note: for negative binomial model only!''')
    optional_group.add_argument('-pbox', '--PlotBOX', dest='do_box', action='store_true',
                                help='''Make box plot
Note: for negative binomial model only!''')
    optional_group.add_argument('-npp', '--NoPostProcess', dest='do_postprocess', action='store_false',
                                help='''Skip post processing. If -npp flag set, annotation of following output tables will be skipped:
1. Normalized expression. \t\t\"_norm_annotated.txt\" will be saved.
2. Differential analysis result. \t\"_diff_annotated.txt\" will be saved. 
Note: If -nx (--NoExecution) is set, post processing will be skipped
''')
    optional_group.add_argument('-an', '--AnNotation', dest='do_anno', action='store',
                                choices=list(anno_dic.keys()),
                                help=anno_help)
    optional_group.add_argument('-nxlsx', '--NoXLSX', dest='do_xlsx', action='store_false',
                                help='''Skip generate excel summary file, which contains the following tabs:
1. SAMPLES\t\t\tSample meta info
2. COUNTS\t\t\tCleaned count table, from _count_cleaned.txt
3. log2 Normalized Expression\tAnnotated/Unannotated log2 transformed normalized expression, from _norm_annotated.txt
4. Differential Analysis\tAnnotated/Unannotated differential analysis table, from _diff_annotated.txt
''')
    ###############################
    #    4. other arguments    #
    ###############################
    other_group = arg_parser.add_argument_group("Other arguments")
    other_group.add_argument('-h', '--help', action="help", help=help_help)
    other_group.add_argument('-v', '--version', action="version", version='%(prog)s: version 1.0',
                             help=version_help)

    # ==============================================
    # parse args
    # ==============================================
    args = arg_parser.parse_args()
    return args


# ===============================================
#  function
# ===============================================
def pre_process(args,
                meta_table,
                count_table,
                spikein_table,
                factor_dic,
                batch_factor_dic,
                factor_dic_wo_batch,
                cont_factor_dic,
                interaction_dic,
                factor_list,
                batch_factor_list,
                cont_factor_list,
                factor_list_wo_ex,
                factor_list_wo_batch,
                interaction_list,
                count_dirname):

    re_disp = args.re_disp
    do_glm = args.do_glm
    do_post = args.do_post
    do_pca = args.do_pca

    # if args.set_bcv:
    #     do_glm = False

    # =====================================================
    # =====================================================
    # Check intersection of samples between count table vs meta table
    # =====================================================
    removed_genes = np.array([])
    removed_samples = np.array([])
    # count_table_filtered = count_table

    # filter count table to remove empty genes and samples
    if args.model == 'nb':
        # Filter the count table according to given thresholds in -f filter_by_gene_sum and -F filter_by_sample_sum
        non_empty_row_indices = np.array(count_table.sum(axis=1)) >= args.filter_by_gene_sum
        non_empty_col_indices = np.array(count_table.sum(axis=0)) >= args.filter_by_sample_sum

        removed_samples = np.array(count_table.columns.values)[np.logical_not(non_empty_col_indices)]
        removed_genes = np.array(count_table.index.values)[np.logical_not(non_empty_row_indices)]

        count_table = count_table.iloc[non_empty_row_indices, non_empty_col_indices]

        # To-do: write list of filtered samples and genes to summary
    # print the removed samples and genes
    print('')
    if len(removed_genes) > 0:
        print('After filtering count table:')
        print('\t', len(removed_genes), '\tgenes/transcripts removed', sep='')
    # =====================================================
    # check the usable samples as the intersection of column names of meta_table vs. header of count_table_filtered
    # =====================================================
    common_samples = np.intersect1d(meta_table.index.values, count_table.columns.values)
    if args.do_pca and len(common_samples) < 4:
        do_pca = False
        print('Warning: number of samples is less than threshold 4, disable PCA plotting!')

    # Remove samples from count table if not in common samples
    if len(common_samples) < count_table.shape[1]:
        removed_samples = np.append(removed_samples, np.setdiff1d(np.array(count_table.columns.values), common_samples))

    print('\t', str(len(removed_samples)), '\tsample removed: ', ', '.join(s for s in removed_samples), sep='')

    # common_sample_index = np.searchsorted(count_table.columns.values, common_samples)
    # Write the filtered count_table to file
    count_table = count_table.loc[:, common_samples]
    count_table.to_csv(count_dirname + args.output_prefix + '_count_cleaned.txt', sep='\t', encoding='utf-8')
    print('The cleaned count table is saved in: ', count_dirname + args.output_prefix + '_count_cleaned.txt', sep='')
    print('\t', count_table.shape[0], '\tgenes/transcripts/rows', sep='')
    print('\t', count_table.shape[1], '\tsamples/columns', sep='')
    print('')

    # check spike-in table with correct number of samples
    if not spikein_table.empty:
        if len(np.intersect1d(np.array(spikein_table.index.values), common_samples)) < len(common_samples):
            raise ValueError('Error: incorrect spike-in library size table!!!')
        spikein_table = spikein_table.loc[common_samples]  # remove the extra samples in spike-in lib size table
        spikein_table.to_csv(count_dirname + args.output_prefix + '_spikein_lib_cleaned.txt', sep='\t', encoding='utf-8')

    # filter meta tables for empty samples and invalid factors

    if len(common_samples) < meta_table.shape[0]:
        removed_samples = np.setdiff1d(np.array(meta_table.index.values), common_samples)
        print('After filtering meta table')
        print('\t', str(len(removed_samples)), '\t samples removed: ', ', '.join(s for s in removed_samples))

    meta_table = meta_table.loc[common_samples]

    # The indicator for factors with more than one value levels, after filtering
    factor_valid_indicator = [len(pd.unique(meta_table.iloc[:, i])) > 1 for i in range(0, meta_table.shape[1])]
    # factor_invalid_indicator = np.logical_not(factor_valid_indicator)

    # factor_invalid_pos_array = np.array([i+1 for i in range(0, len(factor_invalid_indicator)) if factor_invalid_indicator[i]])

    num_factors = meta_table.shape[1]
    tmp_factor_list = list(range(1, num_factors+1))

    # Setup the mapping dictionary btw the original vs. filtered factor indices
    # In factor_mapping_dic, keys are the original indices, values are the new indices
    # Note: the range of factor indices start from 1 !!!!!!!
    factor_mapping_dic = {}
    new_index = 1
    for i in range(meta_table.shape[1]):
        if factor_valid_indicator[i]:
            factor_mapping_dic[i+1] = new_index
            new_index += 1

    # map the total list of factors
    factor_list_mapped = []
    if tmp_factor_list:
        for i in tmp_factor_list:
            if factor_valid_indicator[i-1]:
                factor_list_mapped.append(factor_mapping_dic[i])

    # map the exclude factors
    exclude_factors = []
    if args.exclude_factors:
        for e_factor in args.exclude_factors:
            if len(args.exclude_factors) > num_factors or max(args.exclude_factors) > num_factors or min(args.exclude_factors) < 1:
                raise ValueError("Invalid list of factors to be excluded!!!")
            if factor_valid_indicator[e_factor-1]:
                exclude_factors.append(factor_mapping_dic[e_factor])

    # map the batch effect factors
    batch_factors = []
    if args.batch_factors:
        for b_factor in args.batch_factors:
            if len(args.batch_factors) >= num_factors or max(args.batch_factors) > num_factors or min(args.batch_factors) < 1:
                raise ValueError("Error: invalid list of batch-effect factor(s)")
            if factor_valid_indicator[b_factor-1]:
                batch_factors.append(factor_mapping_dic[b_factor])

    # map the continuous factors
    cont_factors = []
    if args.cont_factors:
        if len(args.cont_factors) >= num_factors or max(args.cont_factors) > num_factors or min(args.cont_factors) < 1:
            raise ValueError("Error: invalid list of continuous factor(s)")
        for c_factor in args.cont_factors:
            if factor_valid_indicator[c_factor-1]:
                cont_factors.append(factor_mapping_dic[c_factor])

    # map the interaction terms
    include_interactions = []
    exclude_interactions = []
    if args.include_interactions:  # Include interaction terms (-ii) are set, discard the exclude interactions terms ( -ei)
        if len(np.unique(args.include_interactions)) < len(args.include_interactions):
            raise ValueError('Error: duplication in include interaction terms')
        for inter_term in args.include_interactions:
            inter_factors = inter_term.split(':')    # Each interactions are joined by ':' (Example   1:2:3)
            inter_factors_integer = [int(s) for s in inter_factors]
            if len(inter_factors_integer) > num_factors:
                raise ValueError('Error: over-length include interaction term!')
                # sys.exit(1)
            if max(inter_factors_integer) > num_factors or min(inter_factors_integer) < 1:
                raise ValueError('Error: incorrect factor indices in include interaction terms!')
                # sys.exit(1)
            inter_factors_integer_mapped = []
            for i_factor in inter_factors_integer:
                if factor_valid_indicator[i_factor-1]:
                    inter_factors_integer_mapped.append(factor_mapping_dic[i_factor])
            inter_term_mapped = ':'.join((str(i_factor) for i_factor in inter_factors_integer_mapped))
            include_interactions.append(inter_term_mapped)
    elif args.exclude_interactions:  # if -ii is not used, use -ei flag to exclude the specified interaction terms
        if len(np.unique(args.exclude_interactions)) < len(args.exclude_interactions):
            raise ValueError('Error: duplication in exclude interaction terms')
        for inter_term in args.exclude_interactions:
            inter_factors = inter_term.split(':')
            inter_factors_integer = [int(s) for s in inter_factors]
            if len(inter_factors_integer) > num_factors:
                raise ValueError('Error: over-length exclude interaction term!')
                # sys.exit(1)
            if max(inter_factors_integer) > num_factors or min(inter_factors_integer) < 1:
                raise ValueError('Error: incorrect factor indices in exclude interaction terms!')
                # sys.exit(1)
            inter_factors_integer_mapped = []
            for i_factor in inter_factors_integer:
                if factor_valid_indicator[i_factor-1]:
                    inter_factors_integer_mapped.append(factor_mapping_dic[i_factor])
            inter_term_mapped = ':'.join((str(i_factor) for i_factor in inter_factors_integer_mapped))
            exclude_interactions.append(inter_term_mapped)

    # ================================================================
    # ================================================================
    # Remove the invalid factors with less than one levels
    meta_table = meta_table.iloc[:, factor_valid_indicator]

    meta_table.to_csv(count_dirname + args.output_prefix + '_meta_cleaned.txt', sep='\t', encoding='utf-8')

    # Initialize the dic and list for single factors
    num_factors = meta_table.shape[1]
    for factor in factor_list_mapped:
        factor_dic[factor] = True
        factor_dic_wo_batch[factor] = True
        factor_list.append(factor)

    # ===========================================
    # Remove exclude factors from factor_dic
    tmp_factor_list_wo_ex = factor_list_mapped.copy()
    if exclude_factors:
        if len(exclude_factors) > num_factors or max(exclude_factors) > num_factors or min(exclude_factors) < 1:
            raise ValueError("Error: Invalid list of factor(s) to be excluded!!!")
            # sys.exit(1)
        for e_factor in exclude_factors:
            del factor_dic[e_factor]
            tmp_factor_list_wo_ex.remove(e_factor)
    for i in tmp_factor_list_wo_ex:
        factor_list_wo_ex.append(i)

    # ===========================================
    # Process batch effect factors
    tmp_factor_list_wo_batch = factor_list_mapped.copy()
    # Collect list of batch effect factors into a dictionary, positions of factors as keys
    if batch_factors:
        if len(batch_factors) >= num_factors or max(batch_factors) > num_factors or min(batch_factors) < 1:
            raise ValueError("Error: invalid list of batch-effect factor(s)")
            # sys.exit(1)
        for b_factor in batch_factors:
            batch_factor_dic[b_factor] = True
            del factor_dic_wo_batch[b_factor]
            tmp_factor_list_wo_batch.remove(b_factor)
            batch_factor_list.append(b_factor)
    for i in tmp_factor_list_wo_batch:
        factor_list_wo_batch.append(i)

    # ==========================================================
    # Process interactions
    # ===========================================
    # generate a full rank of interactions, using the indices of all factor list, !!!without batch-effect factors!!!
    if len(factor_list_wo_batch) > 1:
        for length in range(2, len(factor_list_wo_batch)+1):
            for subset in itertools.combinations(factor_list_wo_batch, length):
                inter_term_sorted = ':'.join(str(factor) for factor in subset)
                interaction_dic[inter_term_sorted] = True
                interaction_list.append(inter_term_sorted)

        if include_interactions:  # Include interaction terms (-ii) are set, discard the exclude interactions terms ( -ei)
            tmp_interaction_dic = {}
            tmp_interaction_list = []
            for inter_term in include_interactions:
                inter_factors = inter_term.split(':')    # Each interactions are joined by ':' (Example   1:2:3
                inter_factors_integer = [int(s) for s in inter_factors]
                if len(inter_factors_integer) > num_factors:
                    raise ValueError('Error: over-length interaction term!')
                    # sys.exit(1)
                if max(inter_factors_integer) > num_factors or min(inter_factors_integer) < 1:
                    raise ValueError('Error: incorrect factor indices in interaction terms!')
                    # sys.exit(1)
                inter_term_sorted = ':'.join((str(factor) for factor in inter_factors_integer))
                if inter_term_sorted not in interaction_dic:    # if this inter term contains batch factor
                    raise ValueError('Error: the interaction terms contain batch effect factors')
                    # sys.exit(1)
                tmp_interaction_dic[inter_term_sorted] = True
            tmp_inter_indices = []
            for key in interaction_dic:
                if key in interaction_list:
                    tmp_inter_indices.append(interaction_list.index(key))
                else:
                    raise ValueError('Error: invalid interaction terms in include interaction list')
                    # sys.exit(1)
                tmp_inter_indices.sort()
            tmp_interaction_list = [interaction_list[tmp] for tmp in tmp_interaction_list]
            interaction_dic = tmp_interaction_dic
            interaction_list = tmp_interaction_list
            print('\tPair-wise comparisons will be made using GLM/LM models')
            # force re-estimating dispersion for each pair-wise comparison
            print('\tDispersion for each pair-wise comparisons will be re-estimated!')
            re_disp = True  # **********************
        elif exclude_interactions:  # if -ii is not used, use -ei flag to exclude the specified interaction terms
            for ex_inter_term in exclude_interactions:
                ex_inter_factors = ex_inter_term.split(':')
                ex_inter_factors_integer = [int(s) for s in ex_inter_factors]
                if len(ex_inter_factors_integer) > num_factors:
                    raise ValueError('Error: over-length interaction term!')
                    # sys.exit(1)
                if max(ex_inter_factors_integer) > num_factors or min(ex_inter_factors_integer) < 1:
                    raise ValueError('Error: incorrect factor indices in interaction terms!')
                    # sys.exit(1)
                ex_inter_term_sorted = ':'.join((str(factor) for factor in ex_inter_factors_integer))
                if ex_inter_term_sorted in interaction_dic:
                    del interaction_dic[ex_inter_term_sorted]
                    interaction_list.remove(ex_inter_term_sorted)
                else:
                    raise ValueError('Error: ' + ex_inter_term_sorted + ' is not a valid interaction term to be excluded')
                    # sys.exit(1)

        # # put the interactions terms into interaction_list for edgeR preparation
        # if len(interaction_dic) > 0:
        #     for key in interaction_dic:
        #         interaction_list.append(key)

    # ==================
    # continuous factors
    if cont_factors:
        if len(args.cont_factors) >= num_factors or max(args.cont_factors) > num_factors or min(args.cont_factors) < 1:
            raise ValueError("Error: invalid list of continuous factor(s)")
            # sys.exit(1)
        for c_factor in args.cont_factors:
            cont_factor_dic[c_factor] = True
            cont_factor_list.append(c_factor)
        do_post = False  # *****************************

    # Return
    return [num_factors, re_disp, do_glm, do_post, do_pca]

# # ===============================================
# #  function
# # ===============================================
# def post_process(args):
#     if args.do_box:


# ===============================================
#  function prep_edger_script
# ===============================================
def prep_edger_script(args,
                      factor_dic,
                      batch_factor_dic,
                      factor_dic_wo_batch,
                      cont_factor_dic,
                      interaction_dic,
                      factor_list,
                      factor_name_list,
                      batch_factor_list,
                      factor_list_wo_ex,
                      factor_list_wo_batch,
                      interaction_list,
                      re_disp,
                      do_glm,
                      do_post,
                      do_pca,
                      count_dirname):
    # The edgeR script

    spikein = False
    # ===========================================
    # Start pre-processing
    # *******************************************
    script = ['#!/usr/bin/env RScript', 'rm(list=ls())', 'options(stringsAsFactors = FALSE)',
              '', 'suppressMessages(library(\"methods\"))']

    if args.model == 'nb':
        script.append('suppressMessages(library(\"edgeR\"))')
    script.append('')

    script.append('# Read in count table, each row is a gene/transcript, each col is a sample')
    script.append('count <- read.table(\"' + count_dirname + args.output_prefix + '_count_cleaned.txt\", header=T, row.names=1, sep=\"\\t\")')
    script.append('')

    script.append('# Replace missing values (NA) in count table with row means (gene means)')
    script.append('na.indices <- which(is.na(count), arr.ind=T)')
    script.append('if (nrow(na.indices) > 0 )    count[na.indices] <- rowMeans(count, na.rm=T)[na.indices[,1]]')
    script.append('')

    script.append('# Read in meta table, each row is a sample, each col is a sample')
    script.append('meta <- read.table(\"' + count_dirname + args.output_prefix + '_meta_cleaned.txt\", header=T, row.names=1, sep=\"\\t\")')
    script.append('rownames(meta) <- make.names(rownames(meta))')
    script.append('')

    script.append('# Re-order the samples/cols in count table according to the meta table')
    script.append('count <- count[, rownames(meta)]')
    script.append('# Re-order the samples/rows in meta table according to the count table')
    script.append('meta <- as.data.frame(meta[colnames(count), ])')
    script.append('')

    # set up simple version of design with all factors, except the batch effect factors
    if len(factor_list_wo_batch) > 0:
        script.append('simpledesign <- factor(paste(meta[,' + '], meta[,'.join(str(s) for s in factor_list_wo_batch) + '], sep=\".\"))')
        script.append('')

        # make a list of all pair-wise comparisons to run
        if do_post:
            # define no-batch design, for picking which pair-wise comparisons to run
            script.append('# Finding pair-wise comparisons to run')
            script.append('design.nobatch <- data.frame(meta[,c(' + ','.join(str(s) for s in factor_list_wo_batch) + ')])')
            # combine simple design with full no-batch design
            script.append('combined_design <- cbind(simpledesign, design.nobatch)')
            script.append('pairwise.list <- c()')
            # make a list of all pairwise comparisons to run
            script.append('for (i in 1:(length(levels(factor(simpledesign)))-1)) {')
            script.append('    for(j in (i+1):length(levels(factor(simpledesign)))) {')
            # check that we only changed 1 level, otherwise leave from this iteration of the loop
            script.append('        # Check that only 1 factor level changed for this pair-wise comparison')
            script.append('        a <- data.frame(combined_design[combined_design$simpledesign == levels(factor(simpledesign))[i], 2:ncol(combined_design)])[1, ]')
            script.append('        b <- data.frame(combined_design[combined_design$simpledesign == levels(factor(simpledesign))[j], 2:ncol(combined_design)])[1, ]')
            script.append('        diff <- 0')
            script.append('        for( k in 1:length(a) )')
            script.append('            if( a[k] != b[k] )    diff <- diff + 1')
            script.append('        if( diff == 1 )  pairwise.list <- rbind(pairwise.list, c(i,j))')
            script.append('    }')
            script.append('}')
            script.append('')

    # read in spike-in controls if provided
    if args.input_spikein:
        script.append('# Reading in spikein lib size and matching it to count table')
        script.append('spikes <- read.table(\"' + count_dirname + args.output_prefix + '_spikein_lib_cleaned.txt\", header=T, row.names=1, sep=\"\\t\")')
        script.append('rownames(spikes) <- make.names(rownames(spikes))')
        script.append('spikes <- as.data.frame(spikes[colnames(count),])')
        spikein = True
    # *******************************************
    # End pre-processing
    # ===========================================

    # ==============================================================================================
    # ==============================================================================================

    # {{===========================================

    # Start making design matrix
    for i in factor_list:
        if i in cont_factor_dic:
            script.append('F' + str(i) + ' <- as.numeric(meta[, ' + str(i) + '])')
        else:
            script.append('F' + str(i) + ' <- factor(meta[, ' + str(i) + '])')
    factor_line = ' ,F'.join(str(i) for i in factor_list)

    # set up formula; num_terms tracks number of terms in model
    num_tests = len(factor_list_wo_ex)
    formula = ' + '.join('F'+str(i) for i in factor_list_wo_ex)
    header = '\"#ID\"' + ','.join('\"'+factor_name_list[i-1]+':PValue\"' for i in factor_list_wo_ex)

    # include interaction terms in formula if requested
    if args.do_inter:
        for key in interaction_list:
            indices_list = key.split(':')
            header = header + ',\"' + '::'.join(factor_name_list[int(i)-1] for i in indices_list) + ':PValue\"'
            formula += ' + ' + ':'.join('F'+i for i in indices_list)
            num_tests = num_tests + 1

    # reset header to just #ID if no glm testing
    if not do_glm:
        header = '\"#ID\"'

    if do_glm or do_post:
        script.append('design <- model.matrix(~' + formula + ' )')
        script.append('')

    # set dispersion if requested
    disp_val = ''
    # disp_comment = ''
    if args.set_bcv:
        script.append('bcv <- ' + str(args.set_bcv))
        disp_val = ', dispersion=bcv^2'
        # disp_comment = '#'

    # End marking design matrix
    # ===========================================}}

    # ======================================================================================
    # ======================================================================================

    # ===========================================
    # edgeR differential analysis
    if args.model == 'nb':
        script.append('diff <- DGEList(counts=count)')
        # change library size to spike-ins if provided
        if spikein:
            script.append('diff$samples$lib.size <- spikes[,1]')
        elif args.do_norm_diff:
            script.append('diff <- calcNormFactors(diff)')
        # get normalized values in transcripts per million (cpm)
        script.append('')
        script.append('# Writing normalized expression')

        # printf OUT ("%s\n", & write_stats("diff", "$out_lib.txt", $writetolib, "Normalization" ) );
        script.append('norm <- cpm(diff, log=' + str(args.do_logcpm).upper() + ', normalized.lib.size=' + str(args.do_norm_exprs).upper() + ')')
        script.append('')

        # {{=========================================
        #  Start batch effect correction if requested
        if batch_factor_list:
            script.append('# Perform batch effect correction')
            script.append('batch_design <- factor(paste(' + ','.join('meta[,'+str(i)+']' for i in batch_factor_list) + '))')
            script.append('norm <- removeBatchEffect(norm, batch_design)')
            script.append('')

        # write out normalized values
        script.append('n = c(\"ID\", colnames(norm))')
        script.append('write.table(t(n), \"' + count_dirname + args.output_prefix + '_norm.txt\", sep=\"\\t\", row.names=F, col.names=F, quote=F)')
        script.append('write.table(norm, \"' + count_dirname + args.output_prefix + '_norm.txt\", sep=\"\\t\", col.names=F, quote=F, append=T)')
        script.append('')
        # End batch effect correction
        # ==========================================}}

        # compute averages per sample group and write those as well
        # & average_data(factors)
        script.append('datalist <- unlist(as.data.frame(norm))')
        script.append('genes <- rep(rownames(norm),length(simpledesign))')
        script.append('factors <- rep(as.vector(simpledesign),each=nrow(norm))')
        script.append('avg <- aggregate(datalist ~ genes + factors, FUN=mean)')
        script.append('avg2 <- matrix(avg$datalist, nrow=nrow(norm))')
        script.append('rownames(avg2) <- avg$genes[1:nrow(avg2)]')
        script.append('colnames(avg2) <- avg[avg$genes==avg$genes[1],]$factors')
        script.append('avg2 <- avg2[rownames(norm),]')
        script.append('avg_n <- c(\"ID\",colnames(avg2))')
        script.append('write.table(t(avg_n),\"' + count_dirname + args.output_prefix + '_avg.txt\",sep=\"\\t\",row.names=F,col.names=F,quote=F)')
        script.append('write.table(avg2,\"' + count_dirname + args.output_prefix + '_avg.txt\",sep=\"\\t\",col.names=F,quote=F,append=T)')

        # start on differential stats
        # k tracks number of tests done
        count_test = 0

        # test each factor/factor interaction
        if do_glm:
            # prepare GLM model
            script.append('# Computing normalization and dispersion parameters for GLM fitting')
            if not spikein and args.do_norm_diff:
                script.append('diff <- calcNormFactors(diff)')

            if args.disp_method == 'single':
                script.append('diff <- suppressWarnings(estimateDisp(diff, design))')
            else:
                script.append('diff <- suppressWarnings(estimateGLMCommonDisp(diff, design))')
                script.append('diff <- suppressWarnings(estimateGLMTrendedDisp(diff, design))')
                script.append('diff <- suppressWarnings(estimateGLMTagwiseDisp(diff, design))')

            if args.do_bcv:
                script.append('pdf(\"' + count_dirname + args.output_prefix + '_multigroupBCV.pdf\")')
                script.append('plotBCV(diff)')
                script.append('dev.off()')

            # printf OUT ("%s\n", & write_stats( "diff", "$out_lib.txt", $writetolib, "multi-group test" ) );

            script.append('')

            # {{==========================================================
            # do differential with GLMs
            script.append('# perform differential analysis with GLM')
            if args.set_bcv: # If BCV is set, use LRT modeling
                script.append('fit <- glmFit(diff, design ' + disp_val + ')')
            else:
                script.append('fit <- glmQLFit(diff, design ' + disp_val + ')')
            # ============================================================}}

            # write out fitted coefficients
            script.append('n <- c(\"ID\", colnames(fit$coefficients))')
            script.append('write.table(t(n), file = \"' + count_dirname + args.output_prefix + '_coefficients.txt\", sep=\"\\t\", row.names=F, col.names=F, quote=F)')
            script.append('write.table(fit$coefficients, file = \"' + count_dirname + args.output_prefix + '_coefficients.txt\", sep=\"\\t\", col.names=F, quote=F, append=T)')

            # {{==============================================
            # start tests for each factor and interaction term
            # ************************************************

            script.append('assign <- attr(design, "assign")')
            script.append('coef.total <- length(assign)')
            script.append('coef.num <- unique(assign)')
            script.append('coef.count <- 1')

            # test single factors
            for factor in factor_list:
                if factor in factor_dic:
                    script.append('')
                    script.append('# Testing factor ' + str(factor) + ' of ' + str(len(factor_dic)) + ': ' + factor_name_list[factor-1])
                    script.append('coef.range <- which(assign == coef.count)')
                    script.append('coef.start <- min(coef.range)')
                    script.append('coef.end <- max(coef.range)')
                    if args.set_bcv:
                        script.append('lrt <- glmLRT(fit, coef=coef.start:coef.end)')
                        script.append('FDR <- p.adjust(lrt$table$PValue,method=\"BH\")')
                        script.append('lrt$table <- cbind(lrt$table, FDR)')
                        script.append('current_result <- cbind(lrt$table$PValue, lrt$table$FDR)')
                        script.append('colnames(current_result) <- c(\"' + factor_name_list[factor-1] + ': PValue\",\"' + factor_name_list[factor-1] + ': FDR\")')
                        script.append('rownames(current_result) <- rownames(lrt$table)')
                    else:
                        script.append('qlf <- glmQLFTest(fit, coef=coef.start:coef.end)')
                        script.append('FDR <- p.adjust(qlf$table$PValue,method=\"BH\")')
                        script.append('qlf$table <- cbind(qlf$table, FDR)')
                        script.append('current_result <- cbind(qlf$table$PValue, qlf$table$FDR)')
                        script.append('colnames(current_result) <- c(\"' + factor_name_list[factor-1] + ': PValue\",\"' + factor_name_list[factor-1] + ': FDR\")')
                        script.append('rownames(current_result) <- rownames(qlf$table)')
                    if count_test == 0:
                        script.append('result <- current_result')
                    else:
                        script.append('result <- cbind(result, current_result)')
                    count_test = count_test + 1
                    script.append('coef.count <- coef.count + 1')

            script.append('')
            script.append('# ====================')

            # interaction tests
            if args.do_inter:
                # track if we've used a factor in an interaction term, so we only adjust the DF once
                for key in interaction_list:
                    key_items = key.split(':')
                    key_name = '*'.join(factor_name_list[int(s)-1] for s in key_items)
                    script.append('')
                    script.append('# Testing for interaction ' + key + ': ' + key_name)
                    script.append('coef.range <- which(assign == coef.count)')
                    script.append('coef.start <- min(coef.range)')
                    script.append('coef.end <- max(coef.range)')
                    if args.set_bcv:
                        script.append('lrt <- glmLRT(fit, coef=coef.start:coef.end)')
                        script.append('FDR <- p.adjust(lrt$table$PValue, method=\"BH\")')
                        script.append('lrt$table <- cbind(lrt$table, FDR)')
                        script.append('current_result <- cbind(lrt$table$PValue, lrt$table$FDR)')
                        script.append('colnames(current_result) <- c(\"' + key_name + ': PValue\",\"' + key_name + ' :FDR\")')
                        script.append('rownames(current_result) <- rownames(lrt$table)')
                    else:
                        script.append('qlf <- glmQLFTest(fit, coef=coef.start:coef.end)')
                        script.append('FDR <- p.adjust(qlf$table$PValue, method=\"BH\")')
                        script.append('qlf$table <- cbind(qlf$table, FDR)')
                        script.append('current_result <- cbind(qlf$table$PValue, qlf$table$FDR)')
                        script.append('colnames(current_result) <- c(\"' + key_name + ': PValue\",\"' + key_name + ' :FDR\")')
                        script.append('rownames(current_result) <- rownames(qlf$table)')
                    if count_test == 0:
                        script.append('result <- current_result')
                    else:
                        script.append('result <- cbind(result, current_result)')
                    count_test = count_test + 1
                    script.append('coef.count <- coef.count + 1')

        if do_post:
            # do pair-wise post-hoc tests
            script.append('')
            script.append('# Perform pairwise post-hoc tests, between all pairs of non-batch effect factor combinations')
            if batch_factor_list:
                script.append('simpledesign_2 <- model.matrix(~simpledesign+F' + '+F'.join(str(i) for i in batch_factor_list) + ')')
            else:
                script.append('simpledesign_2 <- model.matrix(~simpledesign)')
            script.append('simplediff <- DGEList(counts=count, group=simpledesign)')
            # change library size to spike-ins if provided
            if spikein:
                script.append('simplediff$samples$lib.size <- spikes[,1]')
            if not spikein and args.do_norm_diff:
                script.append('simplediff <- calcNormFactors(simplediff)')
            if args.disp_method == 'single':
                script.append('simplediff <- suppressWarnings(estimateDisp(simplediff, simpledesign_2))')
            else:
                script.append('simplediff = suppressWarnings(estimateCommonDisp(simplediff, simpledesign2))')
                script.append('simplediff = suppressWarnings(estimateTagwiseDisp(simplediff, simpledesign2))')

            if args.do_bcv:
                script.append('pdf(\"' + count_dirname + args.output_prefix + '_pairwiseBCV.pdf\")')
                script.append('plotBCV(simplediff)')
                # script.append('$dev.off()')

            # printf OUT ("${disp_comment}%s\n", & write_stats( "simplediff", "$out_lib.txt", $writetolib, "all factor combinations (post-hoc)" ) );
            script.append('for(i in 1:nrow(pairwise.list)) {')
            # if redo_disp, then make subset data sets and designs for each pair-wise test
            if re_disp:
                script.append('    # Estimate dispersion for current pair-wise comparison')
                script.append('    a <- (simpledesign==levels(factor(simpledesign))[pairwise.list[i,1]])')
                script.append('    b <- (simpledesign==levels(factor(simpledesign))[pairwise.list[i,2]])')
                script.append('    simpledesign_2 <- factor(simpledesign[as.logical(a+b)])')

                batch_count = 1
                for b_factor in batch_factor_list:
                    script.append('    batch_factor.' + str(batch_count) + ' <- factor(F' + str(b_factor) + '[as.logical(a+b)])')
                    batch_count += 1
                script.append('    count_2 <- count[, as.logical(a+b)]')
                # if no batch effect specified, include the groups in the DGEList object; otherwise we'll include it as a model later
                if batch_factor_list:
                    script.append('    simplediff <- DGEList(counts=count_2, group=simpledesign_2)')
                else:
                    script.append('    simplediff <- DGEList(counts=count_2)')

                # add batch factors if present
                if batch_factor_list:
                    script.append('    simpledesign_2 <- model.matrix(~ simpledesign_2+batch_factor.' + '+'.join(str(b_factor) for b_factor in batch_factor_list)+')')
                else:
                    script.append('    simpledesign_2 <- model.matrix(~ simpledesign_2)')
                # change library size to spike-ins if provided
                if spikein:
                    script.append('    spikes_2 <- spikes[as.logical(a+b),]')
                    script.append('    simplediff$samples$lib.size <- spikes_2')

                if not spikein and args.do_norm_diff:
                    script.append('    simplediff <- calcNormFactors(simplediff)')

                #  pick which dispersion estimation method we're using
                if args.disp_method == 'single':
                    script.append('    simplediff <- suppressWarnings(estimateDisp(simplediff,simpledesign_2))')
                else:
                    script.append('    simplediff <- suppressWarnings(estimateCommonDisp(simplediff,simpledesign_2))')
                    script.append('    simplediff <- suppressWarnings(estimateTagwiseDisp(simplediff,simpledesign_2))')

                # run exactTest or GLM model, depending on batch effect correction
                if batch_factor_list:
                    script.append('    test <- exactTest(simplediff, pair=c(1,2) ' + disp_val + ')')
                elif args.set_bcv:
                    script.append('    fit <- glmFit(simplediff, simpledesign_2)')
                    script.append('    test <- glmLRT(fit, coef=2)')
                else:
                    script.append('    fit <- glmQLFit(simplediff, simpledesign_2)')
                    script.append('    test <- glmQLTest(fit, coef=2)')

                # plot dispersion if requested
                if args.do_bcv:
                    script.append('    pdfname <- sprintf(\"' + count_dirname + args.output_prefix + '_pairwiseBCV.%s-%s.pdf\", test$comparison[2], test$comparison[1])')
                    script.append('    pdf(pdfname)')
                    script.append('    plotBCV(simplediff)')
                    # print OUT "${disp_comment}    dev.off()\n";

                # printf OUT ("${disp_comment}    %s\n", & write_stats( "simplediff", "$out_lib.txt", $writetolib, "pairwise:",
                # "test\$comparison[2]", "\"vs\"", "test\$comparison[1]" ) );

            else:
                script.append('    test <- exactTest(simplediff, pair=c(pairwise.list[i,1],pairwise.list[i,2]) ' + disp_val + ')')

            # names of comparisons
            if batch_factor_list:
                script.append('    thisname <- paste(test$comparison[2], test$comparison[1], sep=\"/\")')
            else:
                script.append('    thisname <- paste(levels(factor(simpledesign))[pairwise.list[i,2]], '
                              'levels(factor(simpledesign))[pairwise.list[i,1]], sep=\"/\")')

            script.append('    FDR <- p.adjust(test$table$PValue,method=\"BH\")')
            script.append('    test$table <- cbind(test$table, FDR)')

            # re-calculate logFC from averaged values
            if batch_factor_list:
                script.append('    col1 <- (colnames(avg2) == levels(factor(simpledesign))[pairwise.list[i,1]])')
                script.append('    col2 <- (colnames(avg2) == levels(factor(simpledesign))[pairwise.list[i,2]])')
                script.append('    test$table$logFC = avg2[,col1] - avg2[,col2]')

            script.append('    colnames(test$table) <- paste(thisname, \":\", colnames(test$table))')
            script.append('    if ( exists(\"result\") ) {')
            script.append('        result <- cbind(result, test$table[,1:2], test$table[,(ncol(test$table)-1):ncol(test$table)] )')
            script.append('    }')
            script.append('    else {')
            script.append('        result <- cbind(test$table[,1:2], test$table[,(ncol(test$table)-1):ncol(test$table)] )')
            script.append('    }')
            script.append('}')
            script.append('')

        if do_post or do_glm:
            # write differential to output
            script.append('# Write differential analysis results')
            script.append('n <- c(\"ID\",colnames(result))')
            script.append('write.table(t(n),\"' + count_dirname + args.output_prefix + '_diff.txt\",sep=\"\\t\",row.names=F,col.names=F,quote=F)')
            script.append('write.table(result,\"' + count_dirname + args.output_prefix + '_diff.txt\",sep=\"\\t\",col.names=F,quote=F,append=T)')
            script.append('')
            script.append('')

    else:  # linear model !!!!!!!!!!!!
        # pre-define header
        script.append('header <- c(' + header + ')')
        if do_post:
            script.append('for(i in 1:nrow(pairwise.list)) {')
            script.append('    thisname <- paste(levels(factor(simpledesign))[pairwise.list[i,2]], levels(factor(simpledesign))[pairwise.list[i,1]], sep=\"/\")')
            script.append('    header <- c(header, paste(thisname,\":logFC\"), paste(thisname,\": PValue\"))')
            script.append('}')

        # duplicate data for batch effect correction
        if batch_factor_list:
            script.append('data.corrected <- count')

        # set up table for stats
        script.append('diff.res <- data.frame(row.names(count))')
        script.append('for ( r in 1:nrow(count) ) {')
        script.append('    for( j in 2:length(header) ) {')
        script.append('        diff.res[r,j] <- NA')
        script.append('    }')
        script.append('}')
        script.append('colnames(diff.res) <- header')
        # loop over all rows for tests
        script.append('for ( r in 1:nrow(count) ) {')
        script.append('    p <- c(rownames(count)[r])')
        # linear model tests for all factors
        if do_glm:
            script.append('    x <- data.frame(t(count[r,]), ' + factor_line + ')')
            script.append('    colnames(x)[1] <- \"Value\"')
            script.append('    m <- lm( Value ~ $formula, data=x)')
            script.append('    anova = summary(aov(m))')
            script.append('    this.pvalue <- anova[[1]][[\"Pr(>F)\"]]')
            script.append('    p <- c(rownames(data)[r], this.pvalue[1:length(this.pvalue)-1])')

        # re-do group averages with a batch effect correction if needed
        # replace this row in data.corrected with intercept (coefficients[1]) + residuals
        if batch_factor_list:
            script.append('    x <- data.frame(t(count[r,]), ' + factor_line + ')')
            script.append('    colnames(x)[1] <- \"Value\"')
            script.append('    m <- lm( Value ~ F' + '+F'.join(str(b_factor) for b_factor in batch_factor_list) + ', data=x)')
            script.append('    data.corrected[r,] <- m$coefficients[1] + m$residuals')

        # run pair-wise tests
        if do_post:
            script.append('    for(i in 1:nrow(pairwise.list)) {')
            script.append('    a <- (simpledesign==levels(factor(simpledesign))[pairwise.list[i,1]])')
            script.append('    b <- (simpledesign==levels(factor(simpledesign))[pairwise.list[i,2]])')

            # if no batch effect, run t-test
            if not batch_factor_list:
                # use "try" to avoid errors with essentially constant data
                script.append('    test <- try(t.test( count[r,as.logical(a)], count[r,as.logical(b)] ), silent=T)')
                script.append('    if( is(test,\"try-error\") ) {')
                script.append('        this.pvalue <- 1')
                script.append('    } else {')
                script.append('        this.pvalue <- test$p.value')
                script.append('    }')
                script.append('    mean1 <- mean(t(data[r,as.logical(a)]))')
                script.append('    mean2 <- mean(t(data[r,as.logical(b)]))')
                if args.do_logcpm:
                    script.append('    logfc <- mean1 - mean2')
                else:
                    script.append('    logfc <- log2(mean1/mean2)')
                script.append('    p <- c(p, logfc, this.pvalue)')

            else:
                # otherwise run linear model for this pair of values + batch factor
                script.append('    simpledesign_2 <- factor(simpledesign[as.logical(a+b)])')
                for i in range(len(batch_factor_list)):
                    script.append('    batch_factor.' + str(i) + ' <- factor(F' + str(batch_factor_list[i]) + '[as.logical(a+b)])')

                script.append('    x <- data.frame(t(count[r,as.logical(a+b)]), simpledesign_2' +
                              ''.join('batch_factor.'+str(i) for i in range(len(batch_factor_list))) + ')')

                script.append('    colnames(x)[1] <- \"Value\"')
                script.append('    m <- lm( Value ~ simpledesign_2' +
                              ''.join('+batch_factor.'+str(i) for i in range(len(batch_factor_list))) + ', data=x)')
                script.append('    anova <- summary(aov(m))')
                script.append('    this.pvalue <- anova[[1]][[\"Pr(>F)\"]][1]')
                script.append('    mean1 <- mean(t(data.corrected[r,as.logical(a)]))')
                script.append('    mean2 <- mean(t(data.corrected[r,as.logical(b)]))')
                if args.do_logcpm:
                    script.append('    logfc <- mean1 - mean2')
                else:
                    script.append('    logfc <- log2(mean1/mean2)')

                script.append('    p <- c(p, logfc, this.pvalue)')

            script.append('    }')

        script.append('    diff.res[r,] <- t(p)')
        script.append('}')

        # compute group averages; replace data with data.corrected if we did a batch effect correction
        if batch_factor_list:
            script.append('data <- data.corrected')
            script.append('data.corrected_n <- c(\"ID\",colnames(data.corrected))')
            script.append('write.table(t(data.corrected_n),\"' + count_dirname + args.output_prefix+'_norm.txt\",sep=\"\\t\",quote=F,row.names=F,col.names=F)')
            script.append('write.table(data.corrected,\"' + count_dirname + args.output_prefix+'_norm.txt\",sep=\"\\t\",quote=F,append=T,col.names=F)')

        # & average_data($factors);
        # do FDR corrections: for all PValue columns
        if do_glm or do_post:
            script.append('plist <- seq(1:ncol(diff.res))[grepl(\"PValue\",colnames(diff.res),fixed=T)]')
            script.append('for( i in seq(length(plist),1,-1) ) {')
            script.append('    pvals <- diff.res[,plist[i]]')
            script.append('    qname <- gsub(\"PValue\",\"Qvalue\",colnames(diff.res[plist[i]]))')
            script.append('    FDR <- p.adjust(pvals, method=\"BH\")')
            script.append('    diff.res1 <- diff.res[,1:plist[i]]')
            script.append('    if( plist[i] < ncol(diff.res) ) {')
            script.append('        diff.res2 <- data.frame(FDR, diff.res[,(plist[i]+1):ncol(diff.res)],check.names=F)')
            script.append('    } else {')
            script.append('        diff.res2 <- data.frame(FDR)')
            script.append('    }')
            script.append('    colnames(diff.res2)[1] <- qname')
            script.append('    diff.res = data.frame(diff.res1, diff.res2, check.names=F)')
            script.append('}')
            script.append('write.table(diff.res,\"' + count_dirname + args.output_prefix+'_diff.txt\",sep=\"\\t\",col.names=T,row.names=F,quote=F)')
            script.append('')
            script.append('')

    # ===========================================
    # ===========================================
    if args.do_box:
        script.append('data <- read.table(\"' + count_dirname + args.output_prefix + '_norm.txt\", header=T, row.names=1, sep=\"\\t\"')
        script.append('pdf(\"' + count_dirname + args.output_prefix + 'box.pdf\", height=11, width=8.5')
        script.append('boxplot(data, ylab=\"Log2 CPM\"')
        script.append('dev.off()')
        script.append('')

    if do_pca:
        script.append('# ==================================================================')
        script.append('# PCA analysis')
        script.append('# ==================================================================')
        script.append('data <- read.table(\"' + count_dirname + args.output_prefix + '_norm.txt\", header=T, row.names=1, sep=\"\\t\")')
        script.append('meta <- read.table(\"' + count_dirname + args.output_prefix + '_meta_cleaned.txt\", header=T, row.names=1, sep=\"\\t\")')
        script.append('suppressMessages(library(\"ggplot2\"))')
        script.append('suppressMessages(library(\"grid\"))')
        script.append('suppressMessages(library(\"gridExtra\"))')
        script.append('suppressMessages(library(\"FactoMineR\"))')
        script.append('suppressMessages(library(\"factoextra\"))')
        script.append('suppressMessages(library(\"corrplot\"))')
        script.append('suppressMessages(library(\"ggpubr\"))')
        script.append('')

        script.append('data.pca <- PCA(t(data), graph = FALSE)')
        script.append('')
        script.append('pdf(\"' + count_dirname + args.output_prefix + '_PCA.pdf\", width=11, height=8.5)')
        script.append('')
        script.append('p <- fviz_pca_ind(data.pca, axes = c(1, 2), habillage = simpledesign, mean.point = FALSE, pointsize = 5, labelsize=4, repel = T)')
        script.append('ggpubr::ggpar(p,')
        script.append('              title = \"Principal Component Analysis\",')
        script.append('              font.main = c(20, \"bold\", \"black\"),')
        script.append('              subtitle = \"' + args.pca_subtitle + '\",')
        script.append('              xlab = paste0(\"PC1 (\", format(round(get_eig(data.pca)[1, 2], 1), nsmall=1), \"%)\"),')
        script.append('              ylab = paste0(\"PC2 (\", format(round(get_eig(data.pca)[2, 2], 1), nsmall=1), \"%)\"),')
        script.append('              font.x = 16, font.y = 16,')
        script.append('              legend.title = \"Group\", legend = "right", font.legend = c(16, \"bold\", \"black\"),')
        script.append('              ggtheme = theme_minimal(), palette = \"jco\"')
        script.append('              )')
        script.append('')
        script.append('p <- fviz_pca_ind(data.pca, axes = c(1, 3), habillage = simpledesign, mean.point = FALSE, pointsize = 5, labelsize=4, repel = T)')
        script.append('ggpubr::ggpar(p,')
        script.append('              title = \"Principal Component Analysis\",')
        script.append('              font.main = c(20, \"bold\", \"black\"),')
        script.append('              subtitle = \"' + args.pca_subtitle + '\",')
        script.append('              xlab = paste0(\"PC1 (\", format(round(get_eig(data.pca)[1, 2], 1), nsmall=1), \"%)\"),')
        script.append('              ylab = paste0(\"PC3 (\", format(round(get_eig(data.pca)[3, 2], 1), nsmall=1), \"%)\"),')
        script.append('              font.x = 16, font.y = 16,')
        script.append('              legend.title = \"Group\", legend = "right", font.legend = c(16, \"bold\", \"black\"),')
        script.append('              ggtheme = theme_minimal(), palette = \"jco\"')
        script.append('              )')
        script.append('')
        script.append('p <- fviz_pca_ind(data.pca, axes = c(2, 3), habillage = simpledesign, mean.point = FALSE, pointsize = 5, labelsize=4, repel = T)')
        script.append('ggpubr::ggpar(p,')
        script.append('              title = \"Principal Component Analysis\",')
        script.append('              font.main = c(20, \"bold\", \"black\"),')
        script.append('              subtitle = \"' + args.pca_subtitle + '\",')
        script.append('              xlab = paste0(\"PC2 (\", format(round(get_eig(data.pca)[2, 2], 1), nsmall=1), \"%)\"),')
        script.append('              ylab = paste0(\"PC3 (\", format(round(get_eig(data.pca)[3, 2], 1), nsmall=1), \"%)\"),')
        script.append('              font.x = 16, font.y = 16,')
        script.append('              legend.title = \"Group\", legend = "right", font.legend = c(16, \"bold\", \"black\"),')
        script.append('              ggtheme = theme_minimal(), palette = \"jco\"')
        script.append('              )')

        script.append('')
        script.append('dev.off()')

        # script.append('data.pca <- prcomp(t(data))')
        # script.append('data.pca.x <- as.data.frame(data.pca$x)')
        # script.append('#data.pca.x$Group <- do.call(rbind,')
        # script.append('#                            lapply(1:nrow(meta), function(x) {  data.frame(id = df[x,1], val = paste(unlist(df[x,1:ncol(df)]), collapse = \"\")) } )')
        # script.append('#                            )')
        # script.append('data.pca.x$Group <- simpledesign')
        # script.append('')
        # script.append('theme <- theme(panel.background = element_blank(),')
        # script.append('               panel.border=element_rect(fill=NA),')
        # script.append('               panel.grid.major = element_blank(),')
        # script.append('               panel.grid.minor = element_blank(),')
        # script.append('               strip.background=element_blank(),')
        # script.append('               axis.text.x=element_text(colour=\"black\"),')
        # script.append('               axis.text.y=element_text(colour=\"black\"),')
        # script.append('               axis.ticks=element_line(colour=\"black\"),')
        # script.append('               plot.margin=unit(c(1,1,1,1),\"line\"))')
        # script.append('')
        # script.append('percentage <- round(data.pca$sdev / sum(data.pca$sdev) * 100, 2)')
        # script.append('percentage <- paste(colnames(data.pca.x), \"(\", paste(as.character(percentage), \"%\", \")\", sep=\"\"))')
        # script.append('')
        # script.append('pdf(\"' + count_dirname + args.output_prefix + '_PCA.pdf\", width=11, height=8.5)')
        # script.append('')
        # if args.do_pca_label:
        #     script.append('p <- ggplot(data.pca.x, aes(x=PC1, y=PC2, color=Group, label=row.names(data.pca.x)))')
        #     script.append('p <- p + geom_point() + theme + geom_text(size=3, nudge_y = -3) + xlab(percentage[1]) + ylab(percentage[2])')
        #     script.append('print(p)')
        #     script.append('')
        #     script.append('p <- ggplot(data.pca.x, aes(x=PC1, y=PC3, color=Group, label=row.names(data.pca.x)))')
        #     script.append('p <- p + geom_point() + theme + geom_text(size=3, nudge_y = -3) + xlab(percentage[1]) + ylab(percentage[3])')
        #     script.append('print(p)')
        #     script.append('')
        #     script.append('p <- ggplot(data.pca.x, aes(x=PC2, y=PC3, color=Group, label=row.names(data.pca.x)))')
        #     script.append('p <- p + geom_point() + theme + geom_text(size=3, nudge_y = -3) + xlab(percentage[2]) + ylab(percentage[3])')
        #     script.append('print(p)')
        # else:
        #     script.append('p <- ggplot(df.pca.x, aes(x=PC1, y=PC2, color=Group))')
        #     script.append('p <- p + geom_point() + theme + xlab(percentage[1]) + ylab(percentage[2])')
        #     script.append('print(p)')
        #     script.append('')
        #     script.append('p <- ggplot(df.pca.x, aes(x=PC1, y=PC3, color=Group))')
        #     script.append('p <- p + geom_point() + theme + xlab(percentage[1]) + ylab(percentage[3])')
        #     script.append('print(p)')
        #     script.append('')
        #     script.append('p <- ggplot(df.pca.x, aes(x=PC2, y=PC3, color=Group))')
        #     script.append('p <- p + geom_point() + theme + xlab(percentage[2]) + ylab(percentage[3])')
        #     script.append('print(p)')
        # script.append('')
        # script.append('dev.off()')

    ###################################################
    out_file = open(count_dirname + args.output_prefix + '_rscript.r', 'w')
    for line in script:
        out_file.write(line)
        out_file.write('\n')
    out_file.close()

# =======================================<<<<<<<<
# End
# =======================================<<<<<<<<


# ===============================================
#  function prep_deseq2_script
# ===============================================
def prep_deseq2_script(args,
                      factor_dic,
                      batch_factor_dic,
                      factor_dic_wo_batch,
                      cont_factor_dic,
                      interaction_dic,
                      factor_list,
                      factor_name_list,
                      batch_factor_list,
                      factor_list_wo_ex,
                      factor_list_wo_batch,
                      interaction_list,
                      re_disp,
                      do_glm,
                      do_post,
                      do_pca,
                      count_dirname):
    # The edgeR script

    spikein = False
    # ===========================================
    # Start pre-processing
    # *******************************************
    script = ['#!/usr/bin/env RScript', 'rm(list=ls())', 'options(stringsAsFactors = FALSE)',
              '', 'suppressMessages(library(\"methods\"))']

    if args.model == 'nb':
        script.append('suppressMessages(library(\"DESeq2\"))')
    script.append('')

    script.append('# Read in count table, each row is a gene/transcript, each col is a sample')
    script.append('count <- read.table(\"' + count_dirname + args.output_prefix + '_count_cleaned.txt\", header=T, row.names=1, sep=\"\\t\")')
    script.append('')

    script.append('# Replace missing values (NA) in count table with row means (gene means)')
    script.append('na.indices <- which(is.na(count), arr.ind=T)')
    script.append('if (nrow(na.indices) > 0 )    count[na.indices] <- rowMeans(count, na.rm=T)[na.indices[,1]]')
    script.append('')

    script.append('# Read in meta table, each row is a sample, each col is a sample')
    script.append('meta <- read.table(\"' + count_dirname + args.output_prefix + '_meta_cleaned.txt\", header=T, row.names=1, sep=\"\\t\")')
    script.append('rownames(meta) <- make.names(rownames(meta))')
    script.append('')

    script.append('# Re-order the samples/cols in count table according to the meta table')
    script.append('count <- count[, rownames(meta)]')
    script.append('# Re-order the samples/rows in meta table according to the count table')
    script.append('meta <- as.data.frame(meta[colnames(count), ])')
    script.append('')

    # set up simple version of design with all factors, except the batch effect factors
    if len(factor_list_wo_batch) > 0:
        script.append('simpledesign <- factor(paste(meta[,' + '], meta[,'.join(str(s) for s in factor_list_wo_batch) + '], sep=\".\"))')
        script.append('')

        # make a list of all pair-wise comparisons to run
        if do_post:
            # define no-batch design, for picking which pair-wise comparisons to run
            script.append('# Finding pair-wise comparisons to run')
            script.append('design.nobatch <- data.frame(meta[,c(' + ','.join(str(s) for s in factor_list_wo_batch) + ')])')
            # combine simple design with full no-batch design
            script.append('combined_design <- cbind(simpledesign, design.nobatch)')
            script.append('pairwise.list <- c()')
            # make a list of all pairwise comparisons to run
            script.append('for (i in 1:(length(levels(factor(simpledesign)))-1)) {')
            script.append('    for(j in (i+1):length(levels(factor(simpledesign)))) {')
            # check that we only changed 1 level, otherwise leave from this iteration of the loop
            script.append('        # Check that only 1 factor level changed for this pair-wise comparison')
            script.append('        a <- data.frame(combined_design[combined_design$simpledesign == levels(factor(simpledesign))[i], 2:ncol(combined_design)])[1, ]')
            script.append('        b <- data.frame(combined_design[combined_design$simpledesign == levels(factor(simpledesign))[j], 2:ncol(combined_design)])[1, ]')
            script.append('        diff <- 0')
            script.append('        for( k in 1:length(a) )')
            script.append('            if( a[k] != b[k] )    diff <- diff + 1')
            script.append('        if( diff == 1 )  pairwise.list <- rbind(pairwise.list, c(i,j))')
            script.append('    }')
            script.append('}')
            script.append('')

    # read in spike-in controls if provided
    if args.input_spikein:
        script.append('# Reading in spikein lib size and matching it to count table')
        script.append('spikes <- read.table(\"' + count_dirname + args.output_prefix + '_spikein_lib_cleaned.txt\", header=T, row.names=1, sep=\"\\t\")')
        script.append('rownames(spikes) <- make.names(rownames(spikes))')
        script.append('spikes <- as.data.frame(spikes[colnames(count),])')
        spikein = True
    # *******************************************
    # End pre-processing
    # ===========================================

    # ==============================================================================================
    # ==============================================================================================

    # {{===========================================

    # Start making design matrix
    for i in factor_list:
        if i in cont_factor_dic:
            script.append('F' + str(i) + ' <- as.numeric(meta[, ' + str(i) + '])')
        else:
            script.append('F' + str(i) + ' <- factor(meta[, ' + str(i) + '])')
    factor_line = ' ,F'.join(str(i) for i in factor_list)

    # set up formula; num_terms tracks number of terms in model
    num_tests = len(factor_list_wo_ex)
    formula = ' + '.join('F'+str(i) for i in factor_list_wo_ex)
    header = '\"#ID\"' + ','.join('\"'+factor_name_list[i-1]+':PValue\"' for i in factor_list_wo_ex)

    # include interaction terms in formula if requested
    if args.do_inter:
        for key in interaction_list:
            indices_list = key.split(':')
            header = header + ',\"' + '::'.join(factor_name_list[int(i)-1] for i in indices_list) + ':PValue\"'
            formula += ' + ' + ':'.join('F'+i for i in indices_list)
            num_tests = num_tests + 1

    # reset header to just #ID if no glm testing
    if not do_glm:
        header = '\"#ID\"'

    if do_glm or do_post:
        script.append('design <- model.matrix(~' + formula + ' )')
        script.append('')

    # set dispersion if requested
    disp_val = ''
    # disp_comment = ''
    if args.set_bcv:
        script.append('bcv <- ' + str(args.set_bcv))
        disp_val = ', dispersion=bcv^2'
        # disp_comment = '#'

    # End marking design matrix
    # ===========================================}}

    # ======================================================================================
    # ======================================================================================

    # ===========================================
    # edgeR differential analysis
    if args.model == 'nb':
        script.append('diff <- DGEList(counts=count)')
        # change library size to spike-ins if provided
        if spikein:
            script.append('diff$samples$lib.size <- spikes[,1]')
        elif args.do_norm_diff:
            script.append('diff <- calcNormFactors(diff)')
        # get normalized values in transcripts per million (cpm)
        script.append('')
        script.append('# Writing normalized expression')

        # printf OUT ("%s\n", & write_stats("diff", "$out_lib.txt", $writetolib, "Normalization" ) );
        script.append('norm <- cpm(diff, log=' + str(args.do_logcpm).upper() + ', normalized.lib.size=' + str(args.do_norm_exprs).upper() + ')')
        script.append('')

        # {{=========================================
        #  Start batch effect correction if requested
        if batch_factor_list:
            script.append('# Perform batch effect correction')
            script.append('batch_design <- factor(paste("' + ','.join('meta[,'+str(i)+']' for i in batch_factor_list) + '))')
            script.append('norm <- removeBatchEffect(norm, batch_design)')
            script.append('')

        # write out normalized values
        script.append('n = c(\"ID\", colnames(norm))')
        script.append('write.table(t(n), \"' + count_dirname + args.output_prefix + '_norm.txt\", sep=\"\\t\", row.names=F, col.names=F, quote=F)')
        script.append('write.table(norm, \"' + count_dirname + args.output_prefix + '_norm.txt\", sep=\"\\t\", col.names=F, quote=F, append=T)')
        script.append('')
        # End batch effect correction
        # ==========================================}}

        # compute averages per sample group and write those as well
        # & average_data(factors)

        # start on differential stats
        # k tracks number of tests done
        count_test = 0

        # test each factor/factor interaction
        if do_glm:
            # prepare GLM model
            script.append('# Computing normalization and dispersion parameters for GLM fitting')
            if not spikein and args.do_norm_diff:
                script.append('diff <- calcNormFactors(diff)')

            if args.disp_method == 'single':
                script.append('diff <- suppressWarnings(estimateDisp(diff, design))')
            else:
                script.append('diff <- suppressWarnings(estimateGLMCommonDisp(diff, design))')
                script.append('diff <- suppressWarnings(estimateGLMTrendedDisp(diff, design))')
                script.append('diff <- suppressWarnings(estimateGLMTagwiseDisp(diff, design))')

            if args.do_bcv:
                script.append('pdf(\"' + count_dirname + args.output_prefix + '_multigroupBCV.pdf\")')
                script.append('plotBCV(diff)')
                script.append('dev.off()')

            # printf OUT ("%s\n", & write_stats( "diff", "$out_lib.txt", $writetolib, "multi-group test" ) );

            script.append('')

            # {{==========================================================
            # do differential with GLMs
            script.append('# perform differential analysis with GLM')
            if args.set_bcv: # If BCV is set, use LRT modeling
                script.append('fit <- glmFit(diff, design ' + disp_val + ')')
            else:
                script.append('fit <- glmQLFit(diff, design ' + disp_val + ')')
            # ============================================================}}

            # write out fitted coefficients
            script.append('n <- c(\"ID\", colnames(fit$coefficients))')
            script.append('write.table(t(n), file = \"' + count_dirname + args.output_prefix + '_coefficients.txt\", sep=\"\\t\", row.names=F, col.names=F, quote=F)')
            script.append('write.table(fit$coefficients, file = \"' + count_dirname + args.output_prefix + '_coefficients.txt\", sep=\"\\t\", col.names=F, quote=F, append=T)')

            # {{==============================================
            # start tests for each factor and interaction term
            # ************************************************

            script.append('assign <- attr(design, "assign")')
            script.append('coef.total <- length(assign)')
            script.append('coef.num <- unique(assign)')
            script.append('coef.count <- 1')

            # test single factors
            for factor in factor_list:
                if factor in factor_dic:
                    script.append('')
                    script.append('# Testing factor ' + str(factor) + ' of ' + str(len(factor_dic)) + ': ' + factor_name_list[factor-1])
                    script.append('coef.range <- which(assign == coef.count)')
                    script.append('coef.start <- min(coef.range)')
                    script.append('coef.end <- max(coef.range)')
                    if args.set_bcv:
                        script.append('lrt <- glmLRT(fit, coef=coef.start:coef.end)')
                        script.append('FDR <- p.adjust(lrt$table$PValue,method=\"BH\")')
                        script.append('lrt$table <- cbind(lrt$table, FDR)')
                        script.append('current_result <- cbind(lrt$table$PValue, lrt$table$FDR)')
                        script.append('colnames(current_result) <- c(\"' + factor_name_list[factor-1] + ': PValue\",\"' + factor_name_list[factor-1] + ': FDR\")')
                        script.append('rownames(current_result) <- rownames(lrt$table)')
                    else:
                        script.append('qlf <- glmQLFTest(fit, coef=coef.start:coef.end)')
                        script.append('FDR <- p.adjust(qlf$table$PValue,method=\"BH\")')
                        script.append('qlf$table <- cbind(qlf$table, FDR)')
                        script.append('current_result <- cbind(qlf$table$PValue, qlf$table$FDR)')
                        script.append('colnames(current_result) <- c(\"' + factor_name_list[factor-1] + ': PValue\",\"' + factor_name_list[factor-1] + ': FDR\")')
                        script.append('rownames(current_result) <- rownames(qlf$table)')
                    if count_test == 0:
                        script.append('result <- current_result')
                    else:
                        script.append('result <- cbind(result, current_result)')
                    count_test = count_test + 1
                    script.append('coef.count <- coef.count + 1')

            script.append('')
            script.append('# ====================')

            # interaction tests
            if args.do_inter:
                # track if we've used a factor in an interaction term, so we only adjust the DF once
                for key in interaction_list:
                    key_items = key.split(':')
                    key_name = '*'.join(factor_name_list[int(s)-1] for s in key_items)
                    script.append('')
                    script.append('# Testing for interaction ' + key + ': ' + key_name)
                    script.append('coef.range <- which(assign == coef.count)')
                    script.append('coef.start <- min(coef.range)')
                    script.append('coef.end <- max(coef.range)')
                    if args.set_bcv:
                        script.append('lrt <- glmLRT(fit, coef=coef.start:coef.end)')
                        script.append('FDR <- p.adjust(lrt$table$PValue, method=\"BH\")')
                        script.append('lrt$table <- cbind(lrt$table, FDR)')
                        script.append('current_result <- cbind(lrt$table$PValue, lrt$table$FDR)')
                        script.append('colnames(current_result) <- c(\"' + key_name + ': PValue\",\"' + key_name + ' :FDR\")')
                        script.append('rownames(current_result) <- rownames(lrt$table)')
                    else:
                        script.append('qlf <- glmQLFTest(fit, coef=coef.start:coef.end)')
                        script.append('FDR <- p.adjust(qlf$table$PValue, method=\"BH\")')
                        script.append('qlf$table <- cbind(qlf$table, FDR)')
                        script.append('current_result <- cbind(qlf$table$PValue, qlf$table$FDR)')
                        script.append('colnames(current_result) <- c(\"' + key_name + ': PValue\",\"' + key_name + ' :FDR\")')
                        script.append('rownames(current_result) <- rownames(qlf$table)')
                    if count_test == 0:
                        script.append('result <- current_result')
                    else:
                        script.append('result <- cbind(result, current_result)')
                    count_test = count_test + 1
                    script.append('coef.count <- coef.count + 1')

        if do_post:
            # do pair-wise post-hoc tests
            script.append('')
            script.append('# Perform pairwise post-hoc tests, between all pairs of non-batch effect factor combinations')
            if batch_factor_list:
                script.append('simpledesign_2 <- model.matrix(~simpledesign+F' + '+F'.join(str(i) for i in batch_factor_list) + ')')
            else:
                script.append('simpledesign_2 <- model.matrix(~simpledesign)')
            script.append('simplediff <- DGEList(counts=count, group=simpledesign)')
            # change library size to spike-ins if provided
            if spikein:
                script.append('simplediff$samples$lib.size <- spikes[,1]')
            if not spikein and args.do_norm_diff:
                script.append('simplediff <- calcNormFactors(simplediff)')
            if args.disp_method == 'single':
                script.append('simplediff <- suppressWarnings(estimateDisp(simplediff, simpledesign_2))')
            else:
                script.append('simplediff = suppressWarnings(estimateCommonDisp(simplediff, simpledesign2))')
                script.append('simplediff = suppressWarnings(estimateTagwiseDisp(simplediff, simpledesign2))')

            if args.do_bcv:
                script.append('pdf(\"' + count_dirname + args.output_prefix + '_pairwiseBCV.pdf\")')
                script.append('plotBCV(simplediff)')
                # script.append('$dev.off()')

            # printf OUT ("${disp_comment}%s\n", & write_stats( "simplediff", "$out_lib.txt", $writetolib, "all factor combinations (post-hoc)" ) );
            script.append('for(i in 1:nrow(pairwise.list)) {')
            # if redo_disp, then make subset data sets and designs for each pair-wise test
            if re_disp:
                script.append('    # Estimate dispersion for current pair-wise comparison')
                script.append('    a <- (simpledesign==levels(factor(simpledesign))[pairwise.list[i,1]])')
                script.append('    b <- (simpledesign==levels(factor(simpledesign))[pairwise.list[i,2]])')
                script.append('    simpledesign_2 <- factor(simpledesign[as.logical(a+b)])')

                batch_count = 1
                for b_factor in batch_factor_list:
                    script.append('    batch_factor.' + str(batch_count) + ' <- factor(F' + str(b_factor) + '[as.logical(a+b)])')
                    batch_count += 1
                script.append('    count_2 <- count[, as.logical(a+b)]')
                # if no batch effect specified, include the groups in the DGEList object; otherwise we'll include it as a model later
                if batch_factor_list:
                    script.append('    simplediff <- DGEList(counts=count_2, group=simpledesign_2)')
                else:
                    script.append('    simplediff <- DGEList(counts=count_2)')

                # add batch factors if present
                if batch_factor_list:
                    script.append('    simpledesign_2 <- model.matrix(~ simpledesign_2+batch_factor.' + '+'.join(str(b_factor) for b_factor in batch_factor_list)+')')
                else:
                    script.append('    simpledesign_2 <- model.matrix(~ simpledesign_2)')
                # change library size to spike-ins if provided
                if spikein:
                    script.append('    spikes_2 <- spikes[as.logical(a+b),]')
                    script.append('    simplediff$samples$lib.size <- spikes_2')

                if not spikein and args.do_norm_diff:
                    script.append('    simplediff <- calcNormFactors(simplediff)')

                #  pick which dispersion estimation method we're using
                if args.disp_method == 'single':
                    script.append('    simplediff <- suppressWarnings(estimateDisp(simplediff,simpledesign_2))')
                else:
                    script.append('    simplediff <- suppressWarnings(estimateCommonDisp(simplediff,simpledesign_2))')
                    script.append('    simplediff <- suppressWarnings(estimateTagwiseDisp(simplediff,simpledesign_2))')

                # run exactTest or GLM model, depending on batch effect correction
                if batch_factor_list:
                    script.append('    test <- exactTest(simplediff, pair=c(1,2) ' + disp_val + ')')
                elif args.set_bcv:
                    script.append('    fit <- glmFit(simplediff, simpledesign_2)')
                    script.append('    test <- glmLRT(fit, coef=2)')
                else:
                    script.append('    fit <- glmQLFit(simplediff, simpledesign_2)')
                    script.append('    test <- glmQLTest(fit, coef=2)')

                # plot dispersion if requested
                if args.do_bcv:
                    script.append('    pdfname <- sprintf(\"' + count_dirname + args.output_prefix + '_pairwiseBCV.%s-%s.pdf\", test$comparison[2], test$comparison[1])')
                    script.append('    pdf(pdfname)')
                    script.append('    plotBCV(simplediff)')
                    # print OUT "${disp_comment}    dev.off()\n";

                # printf OUT ("${disp_comment}    %s\n", & write_stats( "simplediff", "$out_lib.txt", $writetolib, "pairwise:",
                # "test\$comparison[2]", "\"vs\"", "test\$comparison[1]" ) );

            else:
                script.append('    test <- exactTest(simplediff, pair=c(pairwise.list[i,1],pairwise.list[i,2]) ' + disp_val + ')')

            # names of comparisons
            if batch_factor_list:
                script.append('    thisname <- paste(test$comparison[2], test$comparison[1], sep=\"/\")')
            else:
                script.append('    thisname <- paste(levels(factor(simpledesign))[pairwise.list[i,2]], '
                              'levels(factor(simpledesign))[pairwise.list[i,1]], sep=\"/\")')

            script.append('    FDR <- p.adjust(test$table$PValue,method=\"BH\")')
            script.append('    test$table <- cbind(test$table, FDR)')

            # re-calculate logFC from averaged values
            if batch_factor_list:
                script.append('    col1 <- (colnames(avg2) == levels(factor(simpledesign))[pairwise.list[i,1]])')
                script.append('    col2 <- (colnames(avg2) == levels(factor(simpledesign))[pairwise.list[i,2]])')
                script.append('    test$table$logFC = avg2[,col1] - avg2[,col2]')

            script.append('    colnames(test$table) <- paste(thisname, \":\", colnames(test$table))')
            script.append('    if ( exists(\"result\") ) {')
            script.append('        result <- cbind(result, test$table[,1:2], test$table[,(ncol(test$table)-1):ncol(test$table)] )')
            script.append('    }')
            script.append('    else {')
            script.append('        result <- cbind(test$table[,1:2], test$table[,(ncol(test$table)-1):ncol(test$table)] )')
            script.append('    }')
            script.append('}')
            script.append('')

        if do_post or do_glm:
            # write differential to output
            script.append('# Write differential analysis results')
            script.append('n <- c(\"ID\",colnames(result))')
            script.append('write.table(t(n),\"' + count_dirname + args.output_prefix + '_diff.txt\",sep=\"\\t\",row.names=F,col.names=F,quote=F)')
            script.append('write.table(result,\"' + count_dirname + args.output_prefix + '_diff.txt\",sep=\"\\t\",col.names=F,quote=F,append=T)')
            script.append('')
            script.append('')

    else:  # linear model !!!!!!!!!!!!
        # pre-define header
        script.append('header <- c(' + header + ')')
        if do_post:
            script.append('for(i in 1:nrow(pairwise.list)) {')
            script.append('    thisname <- paste(levels(factor(simpledesign))[pairwise.list[i,1]], levels(factor(simpledesign))[pairwise.list[i,2]], sep=\"/\")')
            script.append('    header <- c(header, paste(thisname,\":logFC\"), paste(thisname,\": PValue\"))')
            script.append('}')

        # duplicate data for batch effect correction
        if batch_factor_list:
            script.append('data.corrected <- count')

        # set up table for stats
        script.append('diff.res <- data.frame(row.names(count))')
        script.append('for ( r in 1:nrow(count) ) {')
        script.append('    for( j in 2:length(header) ) {')
        script.append('        diff.res[r,j] <- NA')
        script.append('    }')
        script.append('}')
        script.append('colnames(diff.res) <- header')
        # loop over all rows for tests
        script.append('for ( r in 1:nrow(count) ) {')
        script.append('    p <- c(rownames(count)[r])')
        # linear model tests for all factors
        if do_glm:
            script.append('    x <- data.frame(t(count[r,]), ' + factor_line + ')')
            script.append('    colnames(x)[1] <- \"Value\"')
            script.append('    m <- lm( Value ~ $formula, data=x)')
            script.append('    anova = summary(aov(m))')
            script.append('    this.pvalue <- anova[[1]][[\"Pr(>F)\"]]')
            script.append('    p <- c(rownames(data)[r], this.pvalue[1:length(this.pvalue)-1])')

        # re-do group averages with a batch effect correction if needed
        # replace this row in data.corrected with intercept (coefficients[1]) + residuals
        if batch_factor_list:
            script.append('    x <- data.frame(t(count[r,]), ' + factor_line + ')')
            script.append('    colnames(x)[1] <- \"Value\"')
            script.append('    m <- lm( Value ~ F' + '+F'.join(str(b_factor) for b_factor in batch_factor_list) + ', data=x)')
            script.append('    data.corrected[r,] <- m$coefficients[1] + m$residuals')

        # run pair-wise tests
        if do_post:
            script.append('    for(i in 1:nrow(pairwise.list)) {')
            script.append('    a <- (simpledesign==levels(factor(simpledesign))[pairwise.list[i,1]])')
            script.append('    b <- (simpledesign==levels(factor(simpledesign))[pairwise.list[i,2]])')

            # if no batch effect, run t-test
            if not batch_factor_list:
                # use "try" to avoid errors with essentially constant data
                script.append('    test <- try(t.test( count[r,as.logical(a)], count[r,as.logical(b)] ), silent=T)')
                script.append('    if( is(test,\"try-error\") ) {')
                script.append('        this.pvalue <- 1')
                script.append('    } else {')
                script.append('        this.pvalue <- test$p.value')
                script.append('    }')
                script.append('    mean1 <- mean(t(data[r,as.logical(a)]))')
                script.append('    mean2 <- mean(t(data[r,as.logical(b)]))')
                if args.do_logcpm:
                    script.append('    logfc <- mean1 - mean2')
                else:
                    script.append('    logfc <- log2(mean1/mean2)')
                script.append('    p <- c(p, logfc, this.pvalue)')

            else:
                # otherwise run linear model for this pair of values + batch factor
                script.append('    simpledesign_2 <- factor(simpledesign[as.logical(a+b)])')
                for i in range(len(batch_factor_list)):
                    script.append('    batch_factor.' + str(i) + ' <- factor(F' + str(batch_factor_list[i]) + '[as.logical(a+b)])')

                script.append('    x <- data.frame(t(count[r,as.logical(a+b)]), simpledesign_2' +
                              ''.join('batch_factor.'+str(i) for i in range(len(batch_factor_list))) + ')')

                script.append('    colnames(x)[1] <- \"Value\"')
                script.append('    m <- lm( Value ~ simpledesign_2' +
                              ''.join('+batch_factor.'+str(i) for i in range(len(batch_factor_list))) + ', data=x)')
                script.append('    anova <- summary(aov(m))')
                script.append('    this.pvalue <- anova[[1]][[\"Pr(>F)\"]][1]')
                script.append('    mean1 <- mean(t(data.corrected[r,as.logical(a)]))')
                script.append('    mean2 <- mean(t(data.corrected[r,as.logical(b)]))')
                if args.do_logcpm:
                    script.append('    logfc <- mean1 - mean2')
                else:
                    script.append('    logfc <- log2(mean1/mean2)')

                script.append('    p <- c(p, logfc, this.pvalue)')

            script.append('    }')

        script.append('    diff.res[r,] <- t(p)')
        script.append('}')

        # compute group averages; replace data with data.corrected if we did a batch effect correction
        if batch_factor_list:
            script.append('data <- data.corrected')
            script.append('data.corrected_n <- c(\"ID\",colnames(data.corrected))')
            script.append('write.table(t(data.corrected_n),\"' + count_dirname + args.output_prefix+'_norm.txt\",sep=\"\\t\",quote=F,row.names=F,col.names=F)')
            script.append('write.table(data.corrected,\"' + count_dirname + args.output_prefix+'_norm.txt\",sep=\"\\t\",quote=F,append=T,col.names=F)')

        # & average_data($factors);
        # do FDR corrections: for all PValue columns
        if do_glm or do_post:
            script.append('plist <- seq(1:ncol(diff.res))[grepl(\"PValue\",colnames(diff.res),fixed=T)]')
            script.append('for( i in seq(length(plist),1,-1) ) {')
            script.append('    pvals <- diff.res[,plist[i]]')
            script.append('    qname <- gsub(\"PValue\",\"Qvalue\",colnames(diff.res[plist[i]]))')
            script.append('    FDR <- p.adjust(pvals, method=\"BH\")')
            script.append('    diff.res1 <- diff.res[,1:plist[i]]')
            script.append('    if( plist[i] < ncol(diff.res) ) {')
            script.append('        diff.res2 <- data.frame(FDR, diff.res[,(plist[i]+1):ncol(diff.res)],check.names=F)')
            script.append('    } else {')
            script.append('        diff.res2 <- data.frame(FDR)')
            script.append('    }')
            script.append('    colnames(diff.res2)[1] <- qname')
            script.append('    diff.res = data.frame(diff.res1, diff.res2, check.names=F)')
            script.append('}')
            script.append('write.table(diff.res,\"' + count_dirname + args.output_prefix+'_diff.txt\",sep=\"\\t\",col.names=T,row.names=F,quote=F)')
            script.append('')
            script.append('')

    # ===========================================
    # ===========================================
    if args.do_box:
        script.append('data <- read.table(\"' + count_dirname + args.output_prefix + '_norm.txt\", header=T, row.names=1, sep=\"\\t\"')
        script.append('pdf(\"' + count_dirname + args.output_prefix + 'box.pdf\", height=11, width=8.5')
        script.append('boxplot(data, ylab=\"Log2 CPM\"')
        script.append('dev.off()')
        script.append('')

    if do_pca:
        script.append('# ==================================================================')
        script.append('# PCA analysis')
        script.append('# ==================================================================')
        script.append('data <- read.table(\"' + count_dirname + args.output_prefix + '_norm.txt\", header=T, row.names=1, sep=\"\\t\")')
        script.append('meta <- read.table(\"' + count_dirname + args.output_prefix + '_meta_cleaned.txt\", header=T, row.names=1, sep=\"\\t\")')
        script.append('suppressMessages(library(\"ggplot2\"))')
        script.append('suppressMessages(library(\"grid\"))')
        script.append('suppressMessages(library(\"gridExtra\"))')
        script.append('')
        script.append('data.pca <- prcomp(t(data))')
        script.append('data.pca.x <- as.data.frame(data.pca$x)')
        script.append('#data.pca.x$Group <- do.call(rbind,')
        script.append('#                            lapply(1:nrow(meta), function(x) {  data.frame(id = df[x,1], val = paste(unlist(df[x,1:ncol(df)]), collapse = \"\")) } )')
        script.append('#                            )')
        script.append('data.pca.x$Group <- simpledesign')
        script.append('')
        script.append('theme <- theme(panel.background = element_blank(),')
        script.append('               panel.border=element_rect(fill=NA),')
        script.append('               panel.grid.major = element_blank(),')
        script.append('               panel.grid.minor = element_blank(),')
        script.append('               strip.background=element_blank(),')
        script.append('               axis.text.x=element_text(colour=\"black\"),')
        script.append('               axis.text.y=element_text(colour=\"black\"),')
        script.append('               axis.ticks=element_line(colour=\"black\"),')
        script.append('               plot.margin=unit(c(1,1,1,1),\"line\"))')
        script.append('')
        script.append('percentage <- round(data.pca$sdev / sum(data.pca$sdev) * 100, 2)')
        script.append('percentage <- paste(colnames(data.pca.x), \"(\", paste(as.character(percentage), \"%\", \")\", sep=\"\"))')
        script.append('')
        script.append('pdf(\"' + count_dirname + args.output_prefix + '_PCA.pdf\", width=11, height=8.5)')
        script.append('')
        if args.do_pca_label:
            script.append('p <- ggplot(data.pca.x, aes(x=PC1, y=PC2, color=Group, label=row.names(data.pca.x)))')
            script.append('p <- p + geom_point() + theme + geom_text(size=3, nudge_y = -3) + xlab(percentage[1]) + ylab(percentage[2])')
            script.append('print(p)')
            script.append('')
            script.append('p <- ggplot(data.pca.x, aes(x=PC1, y=PC3, color=Group, label=row.names(data.pca.x)))')
            script.append('p <- p + geom_point() + theme + geom_text(size=3, nudge_y = -3) + xlab(percentage[1]) + ylab(percentage[3])')
            script.append('print(p)')
            script.append('')
            script.append('p <- ggplot(data.pca.x, aes(x=PC2, y=PC3, color=Group, label=row.names(data.pca.x)))')
            script.append('p <- p + geom_point() + theme + geom_text(size=3, nudge_y = -3) + xlab(percentage[2]) + ylab(percentage[3])')
            script.append('print(p)')
        else:
            script.append('p <- ggplot(df.pca.x, aes(x=PC1, y=PC2, color=Group))')
            script.append('p <- p + geom_point() + theme + xlab(percentage[1]) + ylab(percentage[2])')
            script.append('print(p)')
            script.append('')
            script.append('p <- ggplot(df.pca.x, aes(x=PC1, y=PC3, color=Group))')
            script.append('p <- p + geom_point() + theme + xlab(percentage[1]) + ylab(percentage[3])')
            script.append('print(p)')
            script.append('')
            script.append('p <- ggplot(df.pca.x, aes(x=PC2, y=PC3, color=Group))')
            script.append('p <- p + geom_point() + theme + xlab(percentage[2]) + ylab(percentage[3])')
            script.append('print(p)')
        script.append('')
        script.append('dev.off()')

    out_file = open(count_dirname + args.output_prefix + '_rscript.r', 'w')
    for line in script:
        out_file.write(line)
        out_file.write('\n')
    out_file.close()

# =======================================<<<<<<<<
# End
# =======================================<<<<<<<<

# ===============================================
#  function prep_deseq2_script
# ===============================================
def prep_voom_script(args,
                      factor_dic,
                      batch_factor_dic,
                      factor_dic_wo_batch,
                      cont_factor_dic,
                      interaction_dic,
                      factor_list,
                      factor_name_list,
                      batch_factor_list,
                      factor_list_wo_ex,
                      factor_list_wo_batch,
                      interaction_list,
                      re_disp,
                      do_glm,
                      do_post,
                      do_pca,
                      count_dirname):
    # The edgeR script

    spikein = False
    # ===========================================
    # Start pre-processing
    # *******************************************
    script = ['#!/usr/bin/env RScript', 'rm(list=ls())', 'options(stringsAsFactors = FALSE)',
              '', 'suppressMessages(library(\"methods\"))']

    if args.model == 'nb':
        script.append('suppressMessages(library(\"limma\"))')
        script.append('suppressMessages(library(\"tidyverse\"))')
    script.append('')

    script.append('# Read in count table, each row is a gene/transcript, each col is a sample')
    script.append('count <- read.table(\"' + count_dirname + args.output_prefix + '_count_cleaned.txt\", header=T, row.names=1, sep=\"\\t\")')
    script.append('')

    script.append('# Replace missing values (NA) in count table with row means (gene means)')
    script.append('na.indices <- which(is.na(count), arr.ind=T)')
    script.append('if (nrow(na.indices) > 0 )    count[na.indices] <- rowMeans(count, na.rm=T)[na.indices[,1]]')
    script.append('')

    script.append('# Read in meta table, each row is a sample, each col is a sample')
    script.append('meta <- read.table(\"' + count_dirname + args.output_prefix + '_meta_cleaned.txt\", header=T, row.names=1, sep=\"\\t\")')
    script.append('rownames(meta) <- make.names(rownames(meta))')
    script.append('')

    script.append('# Re-order the samples/cols in count table according to the meta table')
    script.append('count <- count[, rownames(meta)]')
    script.append('# Re-order the samples/rows in meta table according to the count table')
    script.append('meta <- as.data.frame(meta[colnames(count), ])')
    script.append('')

    # set up simple version of design with all factors, except the batch effect factors
    if len(factor_list_wo_batch) > 0:
        script.append('simpledesign <- factor(paste(meta[,' + '], meta[,'.join(str(s) for s in factor_list_wo_batch) + '], sep=\".\"))')
        script.append('')

        # make a list of all pair-wise comparisons to run
        if do_post:
            # define no-batch design, for picking which pair-wise comparisons to run
            script.append('# Finding pair-wise comparisons to run')
            script.append('design.nobatch <- data.frame(meta[,c(' + ','.join(str(s) for s in factor_list_wo_batch) + ')])')
            # combine simple design with full no-batch design
            script.append('combined_design <- cbind(simpledesign, design.nobatch)')
            script.append('pairwise.list <- c()')
            # make a list of all pairwise comparisons to run
            script.append('for (i in 1:(length(levels(factor(simpledesign)))-1)) {')
            script.append('    for(j in (i+1):length(levels(factor(simpledesign)))) {')
            # check that we only changed 1 level, otherwise leave from this iteration of the loop
            script.append('        # Check that only 1 factor level changed for this pair-wise comparison')
            script.append('        a <- data.frame(combined_design[combined_design$simpledesign == levels(factor(simpledesign))[i], 2:ncol(combined_design)])[1, ]')
            script.append('        b <- data.frame(combined_design[combined_design$simpledesign == levels(factor(simpledesign))[j], 2:ncol(combined_design)])[1, ]')
            script.append('        diff <- 0')
            script.append('        for( k in 1:length(a) )')
            script.append('            if( a[k] != b[k] )    diff <- diff + 1')
            script.append('        if( diff == 1 )  pairwise.list <- rbind(pairwise.list, c(i,j))')
            script.append('    }')
            script.append('}')
            script.append('')

    # read in spike-in controls if provided
    if args.input_spikein:
        script.append('# Reading in spikein lib size and matching it to count table')
        script.append('spikes <- read.table(\"' + count_dirname + args.output_prefix + '_spikein_lib_cleaned.txt\", header=T, row.names=1, sep=\"\\t\")')
        script.append('rownames(spikes) <- make.names(rownames(spikes))')
        script.append('spikes <- as.data.frame(spikes[colnames(count),])')
        spikein = True
    # *******************************************
    # End pre-processing
    # ===========================================

    # ==============================================================================================
    # ==============================================================================================

    # {{===========================================

    # Start making design matrix
    for i in factor_list:
        if i in cont_factor_dic:
            script.append('F' + str(i) + ' <- as.numeric(meta[, ' + str(i) + '])')
        else:
            script.append('F' + str(i) + ' <- factor(meta[, ' + str(i) + '])')
    factor_line = ' ,F'.join(str(i) for i in factor_list)

    # set up formula; num_terms tracks number of terms in model
    num_tests = len(factor_list_wo_ex)
    formula = ' + '.join('F'+str(i) for i in factor_list_wo_ex)
    header = '\"#ID\"' + ','.join('\"'+factor_name_list[i-1]+':PValue\"' for i in factor_list_wo_ex)

    # include interaction terms in formula if requested
    if args.do_inter:
        for key in interaction_list:
            indices_list = key.split(':')
            header = header + ',\"' + '::'.join(factor_name_list[int(i)-1] for i in indices_list) + ':PValue\"'
            formula += ' + ' + ':'.join('F'+i for i in indices_list)
            num_tests = num_tests + 1

    # reset header to just #ID if no glm testing
    if not do_glm:
        header = '\"#ID\"'

    if do_glm or do_post:
        script.append('design <- model.matrix(~' + formula + ' )')
        script.append('')

    # set dispersion if requested
    disp_val = ''
    # disp_comment = ''
    if args.set_bcv:
        script.append('bcv <- ' + str(args.set_bcv))
        disp_val = ', dispersion=bcv^2'
        # disp_comment = '#'

    # End marking design matrix
    # ===========================================}}

    # ======================================================================================
    # ======================================================================================

    # ===========================================
    # edgeR differential analysis
    if args.model == 'nb':
        script.append('diff <- DGEList(counts=count)')
        # change library size to spike-ins if provided
        if spikein:
            script.append('diff$samples$lib.size <- spikes[,1]')
        elif args.do_norm_diff:
            script.append('diff <- calcNormFactors(diff)')
        # get normalized values in transcripts per million (cpm)
        script.append('')
        script.append('# Writing normalized expression')

        # printf OUT ("%s\n", & write_stats("diff", "$out_lib.txt", $writetolib, "Normalization" ) );
        script.append('norm <- cpm(diff, log=' + str(args.do_logcpm).upper() + ', normalized.lib.size=' + str(args.do_norm_exprs).upper() + ')')
        script.append('')

        # {{=========================================
        #  Start batch effect correction if requested
        if batch_factor_list:
            script.append('# Perform batch effect correction')
            script.append('batch_design <- factor(paste("' + ','.join('meta[,'+str(i)+']' for i in batch_factor_list) + '))')
            script.append('norm <- removeBatchEffect(norm, batch_design)')
            script.append('')

        # write out normalized values
        script.append('n = c(\"ID\", colnames(norm))')
        script.append('write.table(t(n), \"' + count_dirname + args.output_prefix + '_norm.txt\", sep=\"\\t\", row.names=F, col.names=F, quote=F)')
        script.append('write.table(norm, \"' + count_dirname + args.output_prefix + '_norm.txt\", sep=\"\\t\", col.names=F, quote=F, append=T)')
        script.append('')
        # End batch effect correction
        # ==========================================}}

        # compute averages per sample group and write those as well
        # & average_data(factors)

        # start on differential stats
        # k tracks number of tests done
        count_test = 0

        # test each factor/factor interaction
        if do_glm:
            # prepare GLM model
            script.append('# Computing normalization and dispersion parameters for GLM fitting')
            if not spikein and args.do_norm_diff:
                script.append('diff <- calcNormFactors(diff)')

            if args.disp_method == 'single':
                script.append('diff <- suppressWarnings(estimateDisp(diff, design))')
            else:
                script.append('diff <- suppressWarnings(estimateGLMCommonDisp(diff, design))')
                script.append('diff <- suppressWarnings(estimateGLMTrendedDisp(diff, design))')
                script.append('diff <- suppressWarnings(estimateGLMTagwiseDisp(diff, design))')

            if args.do_bcv:
                script.append('pdf(\"' + count_dirname + args.output_prefix + '_multigroupBCV.pdf\")')
                script.append('plotBCV(diff)')
                script.append('dev.off()')

            # printf OUT ("%s\n", & write_stats( "diff", "$out_lib.txt", $writetolib, "multi-group test" ) );

            script.append('')

            # {{==========================================================
            # do differential with GLMs
            script.append('# perform differential analysis with GLM')
            if args.set_bcv: # If BCV is set, use LRT modeling
                script.append('fit <- glmFit(diff, design ' + disp_val + ')')
            else:
                script.append('fit <- glmQLFit(diff, design ' + disp_val + ')')
            # ============================================================}}

            # write out fitted coefficients
            script.append('n <- c(\"ID\", colnames(fit$coefficients))')
            script.append('write.table(t(n), file = \"' + count_dirname + args.output_prefix + '_coefficients.txt\", sep=\"\\t\", row.names=F, col.names=F, quote=F)')
            script.append('write.table(fit$coefficients, file = \"' + count_dirname + args.output_prefix + '_coefficients.txt\", sep=\"\\t\", col.names=F, quote=F, append=T)')

            # {{==============================================
            # start tests for each factor and interaction term
            # ************************************************

            script.append('assign <- attr(design, "assign")')
            script.append('coef.total <- length(assign)')
            script.append('coef.num <- unique(assign)')
            script.append('coef.count <- 1')

            # test single factors
            for factor in factor_list:
                if factor in factor_dic:
                    script.append('')
                    script.append('# Testing factor ' + str(factor) + ' of ' + str(len(factor_dic)) + ': ' + factor_name_list[factor-1])
                    script.append('coef.range <- which(assign == coef.count)')
                    script.append('coef.start <- min(coef.range)')
                    script.append('coef.end <- max(coef.range)')
                    if args.set_bcv:
                        script.append('lrt <- glmLRT(fit, coef=coef.start:coef.end)')
                        script.append('FDR <- p.adjust(lrt$table$PValue,method=\"BH\")')
                        script.append('lrt$table <- cbind(lrt$table, FDR)')
                        script.append('current_result <- cbind(lrt$table$PValue, lrt$table$FDR)')
                        script.append('colnames(current_result) <- c(\"' + factor_name_list[factor-1] + ': PValue\",\"' + factor_name_list[factor-1] + ': FDR\")')
                        script.append('rownames(current_result) <- rownames(lrt$table)')
                    else:
                        script.append('qlf <- glmQLFTest(fit, coef=coef.start:coef.end)')
                        script.append('FDR <- p.adjust(qlf$table$PValue,method=\"BH\")')
                        script.append('qlf$table <- cbind(qlf$table, FDR)')
                        script.append('current_result <- cbind(qlf$table$PValue, qlf$table$FDR)')
                        script.append('colnames(current_result) <- c(\"' + factor_name_list[factor-1] + ': PValue\",\"' + factor_name_list[factor-1] + ': FDR\")')
                        script.append('rownames(current_result) <- rownames(qlf$table)')
                    if count_test == 0:
                        script.append('result <- current_result')
                    else:
                        script.append('result <- cbind(result, current_result)')
                    count_test = count_test + 1
                    script.append('coef.count <- coef.count + 1')

            script.append('')
            script.append('# ====================')

            # interaction tests
            if args.do_inter:
                # track if we've used a factor in an interaction term, so we only adjust the DF once
                for key in interaction_list:
                    key_items = key.split(':')
                    key_name = '*'.join(factor_name_list[int(s)-1] for s in key_items)
                    script.append('')
                    script.append('# Testing for interaction ' + key + ': ' + key_name)
                    script.append('coef.range <- which(assign == coef.count)')
                    script.append('coef.start <- min(coef.range)')
                    script.append('coef.end <- max(coef.range)')
                    if args.set_bcv:
                        script.append('lrt <- glmLRT(fit, coef=coef.start:coef.end)')
                        script.append('FDR <- p.adjust(lrt$table$PValue, method=\"BH\")')
                        script.append('lrt$table <- cbind(lrt$table, FDR)')
                        script.append('current_result <- cbind(lrt$table$PValue, lrt$table$FDR)')
                        script.append('colnames(current_result) <- c(\"' + key_name + ': PValue\",\"' + key_name + ' :FDR\")')
                        script.append('rownames(current_result) <- rownames(lrt$table)')
                    else:
                        script.append('qlf <- glmQLFTest(fit, coef=coef.start:coef.end)')
                        script.append('FDR <- p.adjust(qlf$table$PValue, method=\"BH\")')
                        script.append('qlf$table <- cbind(qlf$table, FDR)')
                        script.append('current_result <- cbind(qlf$table$PValue, qlf$table$FDR)')
                        script.append('colnames(current_result) <- c(\"' + key_name + ': PValue\",\"' + key_name + ' :FDR\")')
                        script.append('rownames(current_result) <- rownames(qlf$table)')
                    if count_test == 0:
                        script.append('result <- current_result')
                    else:
                        script.append('result <- cbind(result, current_result)')
                    count_test = count_test + 1
                    script.append('coef.count <- coef.count + 1')

        if do_post:
            # do pair-wise post-hoc tests
            script.append('')
            script.append('# Perform pairwise post-hoc tests, between all pairs of non-batch effect factor combinations')
            if batch_factor_list:
                script.append('simpledesign_2 <- model.matrix(~simpledesign+F' + '+F'.join(str(i) for i in batch_factor_list) + ')')
            else:
                script.append('simpledesign_2 <- model.matrix(~simpledesign)')
            script.append('simplediff <- DGEList(counts=count, group=simpledesign)')
            # change library size to spike-ins if provided
            if spikein:
                script.append('simplediff$samples$lib.size <- spikes[,1]')
            if not spikein and args.do_norm_diff:
                script.append('simplediff <- calcNormFactors(simplediff)')
            if args.disp_method == 'single':
                script.append('simplediff <- suppressWarnings(estimateDisp(simplediff, simpledesign_2))')
            else:
                script.append('simplediff = suppressWarnings(estimateCommonDisp(simplediff, simpledesign2))')
                script.append('simplediff = suppressWarnings(estimateTagwiseDisp(simplediff, simpledesign2))')

            if args.do_bcv:
                script.append('pdf(\"' + count_dirname + args.output_prefix + '_pairwiseBCV.pdf\")')
                script.append('plotBCV(simplediff)')
                # script.append('$dev.off()')

            # printf OUT ("${disp_comment}%s\n", & write_stats( "simplediff", "$out_lib.txt", $writetolib, "all factor combinations (post-hoc)" ) );
            script.append('for(i in 1:nrow(pairwise.list)) {')
            # if redo_disp, then make subset data sets and designs for each pair-wise test
            if re_disp:
                script.append('    # Estimate dispersion for current pair-wise comparison')
                script.append('    a <- (simpledesign==levels(factor(simpledesign))[pairwise.list[i,1]])')
                script.append('    b <- (simpledesign==levels(factor(simpledesign))[pairwise.list[i,2]])')
                script.append('    simpledesign_2 <- factor(simpledesign[as.logical(a+b)])')

                batch_count = 1
                for b_factor in batch_factor_list:
                    script.append('    batch_factor.' + str(batch_count) + ' <- factor(F' + str(b_factor) + '[as.logical(a+b)])')
                    batch_count += 1
                script.append('    count_2 <- count[, as.logical(a+b)]')
                # if no batch effect specified, include the groups in the DGEList object; otherwise we'll include it as a model later
                if batch_factor_list:
                    script.append('    simplediff <- DGEList(counts=count_2, group=simpledesign_2)')
                else:
                    script.append('    simplediff <- DGEList(counts=count_2)')

                # add batch factors if present
                if batch_factor_list:
                    script.append('    simpledesign_2 <- model.matrix(~ simpledesign_2+batch_factor.' + '+'.join(str(b_factor) for b_factor in batch_factor_list)+')')
                else:
                    script.append('    simpledesign_2 <- model.matrix(~ simpledesign_2)')
                # change library size to spike-ins if provided
                if spikein:
                    script.append('    spikes_2 <- spikes[as.logical(a+b),]')
                    script.append('    simplediff$samples$lib.size <- spikes_2')

                if not spikein and args.do_norm_diff:
                    script.append('    simplediff <- calcNormFactors(simplediff)')

                #  pick which dispersion estimation method we're using
                if args.disp_method == 'single':
                    script.append('    simplediff <- suppressWarnings(estimateDisp(simplediff,simpledesign_2))')
                else:
                    script.append('    simplediff <- suppressWarnings(estimateCommonDisp(simplediff,simpledesign_2))')
                    script.append('    simplediff <- suppressWarnings(estimateTagwiseDisp(simplediff,simpledesign_2))')

                # run exactTest or GLM model, depending on batch effect correction
                if batch_factor_list:
                    script.append('    test <- exactTest(simplediff, pair=c(1,2) ' + disp_val + ')')
                elif args.set_bcv:
                    script.append('    fit <- glmFit(simplediff, simpledesign_2)')
                    script.append('    test <- glmLRT(fit, coef=2)')
                else:
                    script.append('    fit <- glmQLFit(simplediff, simpledesign_2)')
                    script.append('    test <- glmQLTest(fit, coef=2)')

                # plot dispersion if requested
                if args.do_bcv:
                    script.append('    pdfname <- sprintf(\"' + count_dirname + args.output_prefix + '_pairwiseBCV.%s-%s.pdf\", test$comparison[2], test$comparison[1])')
                    script.append('    pdf(pdfname)')
                    script.append('    plotBCV(simplediff)')
                    # print OUT "${disp_comment}    dev.off()\n";

                # printf OUT ("${disp_comment}    %s\n", & write_stats( "simplediff", "$out_lib.txt", $writetolib, "pairwise:",
                # "test\$comparison[2]", "\"vs\"", "test\$comparison[1]" ) );

            else:
                script.append('    test <- exactTest(simplediff, pair=c(pairwise.list[i,1],pairwise.list[i,2]) ' + disp_val + ')')

            # names of comparisons
            if batch_factor_list:
                script.append('    thisname <- paste(test$comparison[2], test$comparison[1], sep=\"/\")')
            else:
                script.append('    thisname <- paste(levels(factor(simpledesign))[pairwise.list[i,2]], '
                              'levels(factor(simpledesign))[pairwise.list[i,1]], sep=\"/\")')

            script.append('    FDR <- p.adjust(test$table$PValue,method=\"BH\")')
            script.append('    test$table <- cbind(test$table, FDR)')

            # re-calculate logFC from averaged values
            if batch_factor_list:
                script.append('    col1 <- (colnames(avg2) == levels(factor(simpledesign))[pairwise.list[i,1]])')
                script.append('    col2 <- (colnames(avg2) == levels(factor(simpledesign))[pairwise.list[i,2]])')
                script.append('    test$table$logFC = avg2[,col1] - avg2[,col2]')

            script.append('    colnames(test$table) <- paste(thisname, \":\", colnames(test$table))')
            script.append('    if ( exists(\"result\") ) {')
            script.append('        result <- cbind(result, test$table[,1:2], test$table[,(ncol(test$table)-1):ncol(test$table)] )')
            script.append('    }')
            script.append('    else {')
            script.append('        result <- cbind(test$table[,1:2], test$table[,(ncol(test$table)-1):ncol(test$table)] )')
            script.append('    }')
            script.append('}')
            script.append('')

        if do_post or do_glm:
            # write differential to output
            script.append('# Write differential analysis results')
            script.append('n <- c(\"ID\",colnames(result))')
            script.append('write.table(t(n),\"' + count_dirname + args.output_prefix + '_diff.txt\",sep=\"\\t\",row.names=F,col.names=F,quote=F)')
            script.append('write.table(result,\"' + count_dirname + args.output_prefix + '_diff.txt\",sep=\"\\t\",col.names=F,quote=F,append=T)')
            script.append('')
            script.append('')

    else:  # linear model !!!!!!!!!!!!
        # pre-define header
        script.append('header <- c(' + header + ')')
        if do_post:
            script.append('for(i in 1:nrow(pairwise.list)) {')
            script.append('    thisname <- paste(levels(factor(simpledesign))[pairwise.list[i,1]], levels(factor(simpledesign))[pairwise.list[i,2]], sep=\"/\")')
            script.append('    header <- c(header, paste(thisname,\":logFC\"), paste(thisname,\": PValue\"))')
            script.append('}')

        # duplicate data for batch effect correction
        if batch_factor_list:
            script.append('data.corrected <- count')

        # set up table for stats
        script.append('diff.res <- data.frame(row.names(count))')
        script.append('for ( r in 1:nrow(count) ) {')
        script.append('    for( j in 2:length(header) ) {')
        script.append('        diff.res[r,j] <- NA')
        script.append('    }')
        script.append('}')
        script.append('colnames(diff.res) <- header')
        # loop over all rows for tests
        script.append('for ( r in 1:nrow(count) ) {')
        script.append('    p <- c(rownames(count)[r])')
        # linear model tests for all factors
        if do_glm:
            script.append('    x <- data.frame(t(count[r,]), ' + factor_line + ')')
            script.append('    colnames(x)[1] <- \"Value\"')
            script.append('    m <- lm( Value ~ $formula, data=x)')
            script.append('    anova = summary(aov(m))')
            script.append('    this.pvalue <- anova[[1]][[\"Pr(>F)\"]]')
            script.append('    p <- c(rownames(data)[r], this.pvalue[1:length(this.pvalue)-1])')

        # re-do group averages with a batch effect correction if needed
        # replace this row in data.corrected with intercept (coefficients[1]) + residuals
        if batch_factor_list:
            script.append('    x <- data.frame(t(count[r,]), ' + factor_line + ')')
            script.append('    colnames(x)[1] <- \"Value\"')
            script.append('    m <- lm( Value ~ F' + '+F'.join(str(b_factor) for b_factor in batch_factor_list) + ', data=x)')
            script.append('    data.corrected[r,] <- m$coefficients[1] + m$residuals')

        # run pair-wise tests
        if do_post:
            script.append('    for(i in 1:nrow(pairwise.list)) {')
            script.append('    a <- (simpledesign==levels(factor(simpledesign))[pairwise.list[i,1]])')
            script.append('    b <- (simpledesign==levels(factor(simpledesign))[pairwise.list[i,2]])')

            # if no batch effect, run t-test
            if not batch_factor_list:
                # use "try" to avoid errors with essentially constant data
                script.append('    test <- try(t.test( count[r,as.logical(a)], count[r,as.logical(b)] ), silent=T)')
                script.append('    if( is(test,\"try-error\") ) {')
                script.append('        this.pvalue <- 1')
                script.append('    } else {')
                script.append('        this.pvalue <- test$p.value')
                script.append('    }')
                script.append('    mean1 <- mean(t(data[r,as.logical(a)]))')
                script.append('    mean2 <- mean(t(data[r,as.logical(b)]))')
                if args.do_logcpm:
                    script.append('    logfc <- mean1 - mean2')
                else:
                    script.append('    logfc <- log2(mean1/mean2)')
                script.append('    p <- c(p, logfc, this.pvalue)')

            else:
                # otherwise run linear model for this pair of values + batch factor
                script.append('    simpledesign_2 <- factor(simpledesign[as.logical(a+b)])')
                for i in range(len(batch_factor_list)):
                    script.append('    batch_factor.' + str(i) + ' <- factor(F' + str(batch_factor_list[i]) + '[as.logical(a+b)])')

                script.append('    x <- data.frame(t(count[r,as.logical(a+b)]), simpledesign_2' +
                              ''.join('batch_factor.'+str(i) for i in range(len(batch_factor_list))) + ')')

                script.append('    colnames(x)[1] <- \"Value\"')
                script.append('    m <- lm( Value ~ simpledesign_2' +
                              ''.join('+batch_factor.'+str(i) for i in range(len(batch_factor_list))) + ', data=x)')
                script.append('    anova <- summary(aov(m))')
                script.append('    this.pvalue <- anova[[1]][[\"Pr(>F)\"]][1]')
                script.append('    mean1 <- mean(t(data.corrected[r,as.logical(a)]))')
                script.append('    mean2 <- mean(t(data.corrected[r,as.logical(b)]))')
                if args.do_logcpm:
                    script.append('    logfc <- mean1 - mean2')
                else:
                    script.append('    logfc <- log2(mean1/mean2)')

                script.append('    p <- c(p, logfc, this.pvalue)')

            script.append('    }')

        script.append('    diff.res[r,] <- t(p)')
        script.append('}')

        # compute group averages; replace data with data.corrected if we did a batch effect correction
        if batch_factor_list:
            script.append('data <- data.corrected')
            script.append('data.corrected_n <- c(\"ID\",colnames(data.corrected))')
            script.append('write.table(t(data.corrected_n),\"' + count_dirname + args.output_prefix+'_norm.txt\",sep=\"\\t\",quote=F,row.names=F,col.names=F)')
            script.append('write.table(data.corrected,\"' + count_dirname + args.output_prefix+'_norm.txt\",sep=\"\\t\",quote=F,append=T,col.names=F)')

        # & average_data($factors);
        # do FDR corrections: for all PValue columns
        if do_glm or do_post:
            script.append('plist <- seq(1:ncol(diff.res))[grepl(\"PValue\",colnames(diff.res),fixed=T)]')
            script.append('for( i in seq(length(plist),1,-1) ) {')
            script.append('    pvals <- diff.res[,plist[i]]')
            script.append('    qname <- gsub(\"PValue\",\"Qvalue\",colnames(diff.res[plist[i]]))')
            script.append('    FDR <- p.adjust(pvals, method=\"BH\")')
            script.append('    diff.res1 <- diff.res[,1:plist[i]]')
            script.append('    if( plist[i] < ncol(diff.res) ) {')
            script.append('        diff.res2 <- data.frame(FDR, diff.res[,(plist[i]+1):ncol(diff.res)],check.names=F)')
            script.append('    } else {')
            script.append('        diff.res2 <- data.frame(FDR)')
            script.append('    }')
            script.append('    colnames(diff.res2)[1] <- qname')
            script.append('    diff.res = data.frame(diff.res1, diff.res2, check.names=F)')
            script.append('}')
            script.append('write.table(diff.res,\"' + count_dirname + args.output_prefix+'_diff.txt\",sep=\"\\t\",col.names=T,row.names=F,quote=F)')
            script.append('')
            script.append('')

    # ===========================================
    # ===========================================
    if args.do_box:
        script.append('data <- read.table(\"' + count_dirname + args.output_prefix + '_norm.txt\", header=T, row.names=1, sep=\"\\t\"')
        script.append('pdf(\"' + count_dirname + args.output_prefix + 'box.pdf\", height=11, width=8.5')
        script.append('boxplot(data, ylab=\"Log2 CPM\"')
        script.append('dev.off()')
        script.append('')

    if do_pca:
        script.append('# ==================================================================')
        script.append('# PCA analysis')
        script.append('# ==================================================================')
        script.append('data <- read.table(\"' + count_dirname + args.output_prefix + '_norm.txt\", header=T, row.names=1, sep=\"\\t\")')
        script.append('meta <- read.table(\"' + count_dirname + args.output_prefix + '_meta_cleaned.txt\", header=T, row.names=1, sep=\"\\t\")')
        script.append('suppressMessages(library(\"ggplot2\"))')
        script.append('suppressMessages(library(\"grid\"))')
        script.append('suppressMessages(library(\"gridExtra\"))')
        script.append('')
        script.append('data.pca <- prcomp(t(data))')
        script.append('data.pca.x <- as.data.frame(data.pca$x)')
        script.append('#data.pca.x$Group <- do.call(rbind,')
        script.append('#                            lapply(1:nrow(meta), function(x) {  data.frame(id = df[x,1], val = paste(unlist(df[x,1:ncol(df)]), collapse = \"\")) } )')
        script.append('#                            )')
        script.append('data.pca.x$Group <- simpledesign')
        script.append('')
        script.append('theme <- theme(panel.background = element_blank(),')
        script.append('               panel.border=element_rect(fill=NA),')
        script.append('               panel.grid.major = element_blank(),')
        script.append('               panel.grid.minor = element_blank(),')
        script.append('               strip.background=element_blank(),')
        script.append('               axis.text.x=element_text(colour=\"black\"),')
        script.append('               axis.text.y=element_text(colour=\"black\"),')
        script.append('               axis.ticks=element_line(colour=\"black\"),')
        script.append('               plot.margin=unit(c(1,1,1,1),\"line\"))')
        script.append('')
        script.append('percentage <- round(data.pca$sdev / sum(data.pca$sdev) * 100, 2)')
        script.append('percentage <- paste(colnames(data.pca.x), \"(\", paste(as.character(percentage), \"%\", \")\", sep=\"\"))')
        script.append('')
        script.append('pdf(\"' + count_dirname + args.output_prefix + '_PCA.pdf\", width=11, height=8.5)')
        script.append('')
        if args.do_pca_label:
            script.append('p <- ggplot(data.pca.x, aes(x=PC1, y=PC2, color=Group, label=row.names(data.pca.x)))')
            script.append('p <- p + geom_point() + theme + geom_text(size=3, nudge_y = -3) + xlab(percentage[1]) + ylab(percentage[2])')
            script.append('print(p)')
            script.append('')
            script.append('p <- ggplot(data.pca.x, aes(x=PC1, y=PC3, color=Group, label=row.names(data.pca.x)))')
            script.append('p <- p + geom_point() + theme + geom_text(size=3, nudge_y = -3) + xlab(percentage[1]) + ylab(percentage[3])')
            script.append('print(p)')
            script.append('')
            script.append('p <- ggplot(data.pca.x, aes(x=PC2, y=PC3, color=Group, label=row.names(data.pca.x)))')
            script.append('p <- p + geom_point() + theme + geom_text(size=3, nudge_y = -3) + xlab(percentage[2]) + ylab(percentage[3])')
            script.append('print(p)')
        else:
            script.append('p <- ggplot(df.pca.x, aes(x=PC1, y=PC2, color=Group))')
            script.append('p <- p + geom_point() + theme + xlab(percentage[1]) + ylab(percentage[2])')
            script.append('print(p)')
            script.append('')
            script.append('p <- ggplot(df.pca.x, aes(x=PC1, y=PC3, color=Group))')
            script.append('p <- p + geom_point() + theme + xlab(percentage[1]) + ylab(percentage[3])')
            script.append('print(p)')
            script.append('')
            script.append('p <- ggplot(df.pca.x, aes(x=PC2, y=PC3, color=Group))')
            script.append('p <- p + geom_point() + theme + xlab(percentage[2]) + ylab(percentage[3])')
            script.append('print(p)')
        script.append('')
        script.append('dev.off()')

    out_file = open(count_dirname + args.output_prefix + '_rscript.r', 'w')
    for line in script:
        out_file.write(line)
        out_file.write('\n')
    out_file.close()

# =======================================<<<<<<<<
# End
# =======================================<<<<<<<<


# =======================================<<<<<<<<
# annotate_table
# =======================================<<<<<<<<

def annotate_table(args, count_dirname):
    anno_file = anno_dic[args.do_anno]
    anno_table = pd.DataFrame()
    if os.path.isfile(anno_file):
        try:
            anno_table = pd.read_csv(anno_file, sep='\t', header=0, index_col=0, comment='#')
        except IOError:
            print('Error: unable to open annotation table' + anno_file)
            exit(1)

    norm_table = pd.DataFrame()
    if os.path.isfile(count_dirname + args.output_prefix + '_norm.txt'):
        try:
            norm_table = pd.read_csv(count_dirname + args.output_prefix + '_norm.txt', sep='\t', header=0, index_col=0, comment='#')
        except IOError:
            print('Error: norm table ', count_dirname, args.output_prefix, '_norm.txt can not be opened correctly!', sep='')
            exit(1)
    frame = [anno_table, norm_table]
    combined = pd.concat(frame, axis=1, join='inner')
    combined.to_csv(count_dirname + args.output_prefix + '_norm_annotated.txt', sep='\t', index=True, index_label='GeneID', encoding='utf-8')

    diff_table = pd.DataFrame()
    if os.path.isfile(count_dirname + args.output_prefix + '_diff.txt'):
        try:
            diff_table = pd.read_csv(count_dirname + args.output_prefix + '_diff.txt', sep='\t', header=0, index_col=0, comment='#')
        except IOError:
            print('Error: diff table ', count_dirname, args.output_prefix, '_diff.txt can not be opened correctly!', sep='')
            exit(1)
    frame = [anno_table, diff_table]
    combined = pd.concat(frame, axis=1, join='inner')
    combined.to_csv(count_dirname + args.output_prefix + '_diff_annotated.txt', sep='\t', index=True, index_label='GeneID', encoding='utf-8')

# =======================================<<<<<<<<
# End
# =======================================<<<<<<<<


# =======================================<<<<<<<<
# post_process
# =======================================<<<<<<<<
def post_process(args, count_dirname):
    if args.do_anno:
        annotate_table(args, count_dirname)

    if args.do_xlsx:
        if args.do_anno:
            input_file_list = [args.input_meta_table,
                               count_dirname + args.output_prefix + '_count_cleaned.txt',
                               count_dirname + args.output_prefix + '_norm_annotated.txt',
                               count_dirname + args.output_prefix + '_diff_annotated.txt']
        else:
            input_file_list = [args.input_meta_table,
                               count_dirname + args.output_prefix + '_count_cleaned.txt',
                               count_dirname + args.output_prefix + '_norm.txt',
                               count_dirname + args.output_prefix + '_diff.txt']
        sheet_name_list = ['SAMPLES', 'COUNTS', 'log2 Normalized Expression', 'Differential Analysis']
        txt2xlsx.write_xlsx(count_dirname + args.output_prefix + '.xlsx',
                            input_file_list,
                            sheet_name_list)


# ===============================================
# Main function
# ===============================================
def main():

    # Variables
    # count_table_input_file_name = ''            # file name for count table, sample names are col names
    # meta_table_input_file_name = ''             # file name for meta table. The factor names are the col names
    # spike_in_table_input_file_name = ''         # file name for spike-in count table
    # out_put_prefix = ''                         # prefix for output files

    factor_dic = {}                             # Dictionary for included factors in analysis, all_factors - exclude_factors
    batch_factor_dic = {}                       # Dictionary for batch effect factor(s)
    factor_dic_wo_batch = {}                    # Dictionary for factors exclude batch factors, all_factors - batch_factors
    cont_factor_dic = {}                        # Dictionary for continuous factors
    interaction_dic = {}                        # Dictionary for factor interactions

    factor_list = []
    cont_factor_list = []
    batch_factor_list = []
    factor_list_wo_ex = []
    factor_list_wo_batch = []
    interaction_list = []

    # num_factors = 0
    # re_disp = False
    # do_post = False


    # ===========================
    anno_list_table = pd.DataFrame()
    if os.path.isfile(anno_table_file):
        try:
            anno_list_table = pd.read_csv(anno_table_file, sep='\t', header=0, index_col=0, comment='#')
        except IOError:
            print('Error: Can not open list of compiled gene annotations ', anno_table_file, '!', sep='')
            sys.exit(1)
   
    if len(anno_list_table.index.values) == 0 or len(anno_list_table.columns.values) != 1:
        raise ValueError('Error: incorrect list of compiled gene annotations!')
   
    if len(np.unique(anno_list_table.index.values)) < anno_list_table.shape[0]:
        raise ValueError('Error: duplicated names in list of compiled annotations!')
   
    for i in range(anno_list_table.shape[0]):
        anno_dic[anno_list_table.index.values[i]] = anno_list_table.iloc[i, 0].strip()


    # ===========================
    #  get command line args

    args = getargs()

    # Read in the input files and preProcess
    # Note: in meta table, each row is a sample, each column is an variable!!!!!!
    meta_table = pd.DataFrame()
    if os.path.isfile(args.input_meta_table):
        try:
            meta_table = pd.read_csv(args.input_meta_table, sep='\t', header=0, index_col=0, comment='#')
        except IOError:
            print('Error: unable to open meta table ', args.input_meta_table, '!', sep='')
            sys.exit(1)

    # meta_table.set_index(meta_table.columns.values[0])
    if len(np.unique(meta_table.index.values)) < meta_table.shape[0]:
        raise ValueError('Error: duplicated sample names in meta table!')

    if len(np.unique(meta_table.columns.values)) < meta_table.shape[1]:
        raise ValueError('Error: duplicated factor names in meta table!')

    print('Successfully opened ORIGINAL meta table: ', args.input_meta_table, sep='')
    print('\t', meta_table.shape[0], '\tsamples/rows', sep='')
    print('\t', meta_table.shape[1], '\tfactors/columns: ', ', '.join(i for i in meta_table.columns.values), sep='')

    # Get the base path of meta_table
    meta_dirname = os.path.dirname(args.input_meta_table)

    # Note: in count table, each row is a gene, each column is a sample!!!!
    count_table = pd.DataFrame()
    if os.path.isfile(args.input_count_table):
        try:
            count_table = pd.read_csv(args.input_count_table, sep='\t', header=0, index_col=0, comment='#')
        except IOError:
            print('Error: unable to open count table: ', args.input_count_table, '!', sep='')
            sys.exit(1)

    if len(np.unique(count_table.index.values)) < count_table.shape[0]:
        raise ValueError('Error: duplicated gene/transcript names in count table!')

    if len(np.unique(count_table.columns.values)) < count_table.shape[1]:
        raise ValueError('Error: duplicated sample names in count table!')

    print('Successfully opened ORIGINAL count table: ', args.input_count_table, sep='')
    print('\t', count_table.shape[0], '\tgenes/transcripts/rows', sep='')
    print('\t', count_table.shape[1], '\tsamples/columns', sep='')

    count_dirname = os.path.dirname(args.input_count_table)
    if not os.path.samefile(meta_dirname, count_dirname):
        print('Warning: count table and meta table are in different folders')
        print('\tAnalysis result will be saved in ' + count_dirname)
    count_dirname += '/'

    # ===========================================

    spikein_table = pd.DataFrame()
    if args.input_spikein:
        if os.path.isfile(args.input_spikein):
            try:
                spikein_table = pd.read_csv(args.input_spikein, sep='\t', header=0, index_col=0, comment='#')
            except IOError:
                print('Error, unable to open spike-in library size table ', args.input_spikein, '!', sep='')
                sys.exit(1)

        if len(np.unique(spikein_table.index.values)) < spikein_table.shape[0]:
            raise ValueError('Error: duplicated sample names in spikein table')

        if spikein_table.shape[1] < 1:
            raise ValueError('Error: invalid number of columns in spikein table')

        print('Successfully opened spike-in library size table with', args.input_spikein, ': ',
              spikein_table.shape[0], ' samples ORIGINALLY!!!', sep='')

        # spikein_dirname = os.path.dirname(args.spikein_table)
    # ========================================================================================

    [num_factors, re_disp, do_glm, do_post, do_pca] = pre_process(args,
                                                                  meta_table,
                                                                  count_table,
                                                                  spikein_table,
                                                                  factor_dic,
                                                                  batch_factor_dic,
                                                                  factor_dic_wo_batch,
                                                                  cont_factor_dic,
                                                                  interaction_dic,
                                                                  factor_list,
                                                                  batch_factor_list,
                                                                  cont_factor_list,
                                                                  factor_list_wo_ex,
                                                                  factor_list_wo_batch,
                                                                  interaction_list,
                                                                  count_dirname)

    # ==========================================================================================

    factor_name_list = meta_table.columns.values

    # ==========================================================================================
    # Prepare the R script for differential analysis
    prep_edger_script(args,
                      factor_dic,
                      batch_factor_dic,
                      factor_dic_wo_batch,
                      cont_factor_dic,
                      interaction_dic,
                      factor_list,
                      factor_name_list,
                      batch_factor_list,
                      factor_list_wo_ex,
                      factor_list_wo_batch,
                      interaction_list,
                      re_disp,
                      do_glm,
                      do_post,
                      do_pca,
                      count_dirname)

    # ===========================================
    # Prepare run.sh for qsub
    # sh_script = []
    # sh_script.append('#!/bin/bash -l')
    # sh_script.append('#PBS -j oe')
    # sh_script.append('#PBS -l walltime=100:00:00')
    # sh_script.append('#PBS -l mem=12gb')
    # sh_script.append('#PBS -l nodes=1:ppn=4')
    # sh_script.append('#PBS -N ' + args.output_prefix)
    # sh_script.append('#PBS -o ${PBS_JOBNAME}.o${PBS_JOBID}.log')
    # sh_script.append('cd $PBS_O_WORKDIR')
    #
    # sh_script.append('PYTHONPATH=""')
    # sh_script.append('PERL5LIB="":')
    # sh_script.append('R_LIBS=""')
    # sh_script.append('GEM_HOME=""')
    # sh_script.append('PATH=$(getconf PATH)')
    # sh_script.append('source /etc/bashrc')
    # sh_script.append('source /etc/profile')
    # sh_script.append('module purge')
    # sh_script.append('PATH="$PATH:/usr/local/bin"')
    # sh_script.append('module load gcc java-jdk')
    # sh_script.append('module load R/3.5.0')
    # sh_script.append('Rscript ' + args.output_prefix + '_rscript.r > ' + args.output_prefix + '_rscript.log')
    #
    # out_file = open(count_dirname + args.output_prefix + '_rscript_run.sh', 'w')
    # for line in sh_script:
    #     out_file.write(line)
    #     out_file.write('\n')
    # out_file.close()

    # qsub_command =  args.output_prefix + '_rscript_run.sh'

    if not args.do_exec:
        print('NO execution flag set, R script written without running!')
        exit(0)

    # os.environ['PATH'] = sp.run(['getconf', 'PATH']).stdout.decode('ascii')
    # modules.module('purge')
    modules.module('load', 'gcc')
    modules.module('load', 'java-jdk')
    # modules.module('load', 'R/3.5.0')

    process = sp.run('Rscript ' + args.output_prefix + '_rscript.r > ' + args.output_prefix + '_rscript.log', shell=True, capture_output=True)

    exit_status = process.returncode
    if exit_status == 1:
        print('Error: R script for differential analysis failed to run!')
        sys.exit(1)
    print('R script for differential analysis successfully performed!')

    if args.do_postprocess:
        post_process(args, count_dirname)

# =========================================================
# End of main function
# =========================================================


# ====================================================================================================================================
# RUN MAIN function
# ====================================================================================================================================
if __name__ == '__main__':
    main()
