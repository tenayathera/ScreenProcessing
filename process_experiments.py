# merge counts files into a data table, combine reads from multiple sequencing runs,
#  filter by read counts, generate phenotype scores, average replicates

import pandas as pd
import os
import sys
import numpy as np
import scipy as sp
from scipy import stats
import fnmatch
import argparse

from expt_config_parser import parse_expt_config, parse_library_config
from fastqgz_to_counts import make_directory, print_now
from screen_analysis import PlottingObject

defaultLibConfigName = 'library_config.txt'


# a screen processing pipeline that requires just a config file and a directory of supported libraries
# error checking in config parser is fairly robust, so not checking for input errors here
def process_experiments_from_config(config_file, library_directory, generate_plots='png'):
    """load in the supported libraries and sublibraries"""
    try:
        libraries_to_sublibraries, libraries_to_tables = parse_library_config(
            os.path.join(library_directory, defaultLibConfigName))
    except ValueError as err:
        print(' '.join(err.args))
        return

    expt_parameters, parse_status, parse_string = parse_expt_config(config_file, libraries_to_sublibraries)

    print_now(parse_string)

    if parse_status > 0:  # Critical errors in parsing
        print('Exiting due to experiment config file errors\n')
        return

    make_directory(expt_parameters['output_folder'])
    outbase = os.path.join(expt_parameters['output_folder'], expt_parameters['experiment_name'])
    plotting = None

    if generate_plots != 'off':
        plot_directory = os.path.join(expt_parameters['output_folder'], expt_parameters['experiment_name'] + '_plots')
        make_directory(plot_directory)

        plotting = PlottingObject(new_directory=plot_directory, new_image_extension=generate_plots,
                                  new_plot_with_pylab=False)

    # load in library table and filter to requested sublibraries
    print_now('Accessing library information')

    library_table = pd.read_csv(os.path.join(library_directory, libraries_to_tables[expt_parameters['library']]),
                                sep='\t', header=0, index_col=0).sort_index()
    sublib_column = library_table.apply(lambda row: row['sublibrary'].lower() in expt_parameters['sublibraries'],
                                        axis=1)

    if sum(sublib_column) == 0:
        print('After limiting analysis to specified sublibraries, no elements are left')
        return

    library_table[sublib_column].to_csv(outbase + '_librarytable.txt', sep='\t')

    # load in counts, create table of total counts in each and each file as a column
    print_now('Loading counts data')

    column_dict = dict()
    for tup in sorted(expt_parameters['counts_file_list']):
        if tup in column_dict:
            print('Asserting that tuples of condition, replicate, and count file should be unique; '
                  'are the cases where this should not be enforced?')
            raise Exception('condition, replicate, and count file combination already assigned')

        count_series = read_counts_file(tup[2]).reset_index().drop_duplicates('id').set_index(
            'id')  # for now also dropping duplicate ids in counts for overlapping linc sublibraries
        count_series = library_table[sublib_column].align(count_series, axis=0, join='left', fill_value=0)[
            1]  # expand series to fill 0 for every missing entry

        column_dict[tup] = count_series['counts']  # [sublib_column] #then shrink series to only desired sublibraries

    # print column_dict
    counts_table = pd.DataFrame(column_dict)  # , index=library_table[sublib_column].index)
    counts_table.to_csv(outbase + '_rawcountstable.txt', sep='\t')
    counts_table.sum().to_csv(outbase + '_rawcountstable_summary.txt', sep='\t')

    # merge counts for same conditions/replicates, and create summary table
    # save scatter plot before each merger, and histogram of counts post mergers
    print_now('Merging experiment counts split across lanes/indexes')

    expt_groups = counts_table.groupby(level=[0, 1], axis=1)
    merged_counts_table = expt_groups.aggregate(np.sum)
    merged_counts_table.to_csv(outbase + '_mergedcountstable.txt', sep='\t')
    merged_counts_table.sum().to_csv(outbase + '_mergedcountstable_summary.txt', sep='\t')

    if generate_plots != 'off' and max(expt_groups.count().iloc[0]) > 1:
        print_now('-generating scatter plots of counts pre-merger')

        temp_data_dict = {'library': library_table[sublib_column],
                          'premerged counts': counts_table,
                          'counts': merged_counts_table}

        for (phenotype, replicate), countsCols in expt_groups:
            if len(countsCols.columns) == 1:
                continue

            else:
                plotting.premerged_counts_scatter_matrix(temp_data_dict, phenotype, replicate)

    if generate_plots != 'off':
        print_now('-generating sgRNA read count histograms')

        temp_data_dict = {'library': library_table[sublib_column],
                          'counts': merged_counts_table}

        for (phenotype, replicate), countsCol in merged_counts_table.items():
            plotting.counts_histogram(temp_data_dict, phenotype, replicate)

    # create pairs of columns for each comparison, filter to na, then generate sgRNA phenotype score
    print_now('Computing sgRNA phenotype scores')

    growth_value_dict = {(tup[0], tup[1]): tup[2] for tup in expt_parameters['growth_value_tuples']}

    phenotype_list = list(set([item[0] for item in expt_parameters['condition_tuples']]))
    replicate_list = sorted(list(set([item[1] for item in expt_parameters['counts_file_list']])))

    phenotype_score_dict = dict()
    for (phenotype, condition1, condition2) in expt_parameters['condition_tuples']:
        for replicate in replicate_list:
            column1 = merged_counts_table[(condition1, replicate)]
            column2 = merged_counts_table[(condition2, replicate)]
            filt_cols = filter_low_counts(pd.concat((column1, column2), axis=1), expt_parameters['filter_type'],
                                          expt_parameters['minimum_reads'])

            score = compute_phenotype_score(filt_cols[(condition1, replicate)], filt_cols[(condition2, replicate)],
                                            library_table[sublib_column], growth_value_dict[(phenotype, replicate)],
                                            expt_parameters['pseudocount_behavior'], expt_parameters['pseudocount'])

            phenotype_score_dict[(phenotype, replicate)] = score

    if generate_plots != 'off':
        temp_data_dict = {'library': library_table[sublib_column],
                          'counts': merged_counts_table,
                          'phenotypes': pd.DataFrame(phenotype_score_dict)}

        print_now('-generating phenotype histograms and scatter plots')

        for (phenotype, condition1, condition2) in expt_parameters['condition_tuples']:
            for replicate in replicate_list:
                plotting.counts_scatter(temp_data_dict, condition1, replicate, condition2, replicate,
                                        color_by_phenotype_condition=phenotype,
                                        color_by_phenotype_replicate=replicate)

                plotting.phenotype_histogram(temp_data_dict, phenotype, replicate)
                plotting.sg_r_n_as_passing_filter_hist(temp_data_dict, phenotype, replicate)

    # scatterplot sgRNAs for all replicates, then average together and add columns to phenotype score table
    if len(replicate_list) > 1:
        print_now('Averaging replicates')

        for phenotype in phenotype_list:
            rep_cols = pd.DataFrame(
                {(phen, rep): col for (phen, rep), col in phenotype_score_dict.items() if phen == phenotype})
            phenotype_score_dict[(phenotype, 'ave_' + '_'.join(replicate_list))] = rep_cols.mean(axis=1,
                                                                                                 skipna=False)
            # average nan and real to nan; otherwise this could lead to data points with just one rep informing results

    phenotype_table = pd.DataFrame(phenotype_score_dict)
    phenotype_table.to_csv(outbase + '_phenotypetable.txt', sep='\t')

    if len(replicate_list) > 1 and generate_plots != 'off':
        temp_data_dict = {'library': library_table[sublib_column],
                          'phenotypes': phenotype_table}

        print_now('-generating replicate phenotype histograms and scatter plots')

        for phenotype, phengroup in phenotype_table.groupby(level=0, axis=1):
            for i, ((p, rep1), col1) in enumerate(phengroup.items()):
                if rep1[:4] == 'ave_':
                    plotting.phenotype_histogram(temp_data_dict, phenotype, rep1)

                for j, ((p2, rep2), col2) in enumerate(phengroup.items()):
                    if rep2[:4] == 'ave_' or j <= i:
                        continue

                    else:
                        plotting.phenotype_scatter(temp_data_dict, phenotype, rep1, phenotype, rep2)

                        # generate pseudogenes
    neg_table = phenotype_table.loc[library_table[sublib_column].loc[:, 'gene'] == 'negative_control', :]

    if expt_parameters['generate_pseudogene_dist'] != 'off' and len(expt_parameters['analyses']) > 0:
        print('Generating a pseudogene distribution from negative controls')
        sys.stdout.flush()

        pseudo_table_list = []
        pseudo_lib_tables = []
        neg_values = neg_table.values
        neg_columns = neg_table.columns

        if expt_parameters['generate_pseudogene_dist'].lower() == 'manual':
            for pseudogene in range(expt_parameters['num_pseudogenes']):
                rand_indices = np.random.randint(0, len(neg_table), expt_parameters['pseudogene_size'])
                pseudo_table = neg_values[rand_indices, :]
                pseudo_index = [f'pseudo_{pseudogene}_{i}' for i in range(expt_parameters['pseudogene_size'])]
                # so pseudogenes aren't treated as duplicates
                pseudo_seqs = [f'seq_{pseudogene}_{i}' for i in range(expt_parameters['pseudogene_size'])]
                pseudo_table_list.append(pd.DataFrame(pseudo_table, index=pseudo_index, columns=neg_columns))
                pseudo_lib = pd.DataFrame({'gene': [f'pseudo_{pseudogene}'] * expt_parameters['pseudogene_size'],
                                           'transcripts': ['na'] * expt_parameters['pseudogene_size'],
                                           'sequence': pseudo_seqs}, index=pseudo_index)
                pseudo_lib_tables.append(pseudo_lib)

        elif expt_parameters['generate_pseudogene_dist'].lower() == 'auto':
            for pseudogene, (gene, group) in enumerate(
                    library_table[sublib_column].drop_duplicates(['gene', 'sequence']).groupby('gene')):
                if gene == 'negative_control':
                    continue
                for transcript, (transcriptName, transcriptGroup) in enumerate(group.groupby('transcripts')):
                    rand_indices = np.random.randint(0, len(neg_table), len(transcriptGroup))
                    pseudo_table = neg_values[rand_indices, :]
                    pseudo_index = [f'pseudo_{pseudogene}_{transcript}_{i}' for i in range(len(transcriptGroup))]
                    pseudo_seqs = [f'seq_{pseudogene}_{transcript}_{i}' for i in range(len(transcriptGroup))]
                    pseudo_table_list.append(pd.DataFrame(pseudo_table, index=pseudo_index, columns=neg_columns))
                    pseudo_lib = pd.DataFrame({'gene': [f'pseudo_{pseudogene}'] * len(transcriptGroup),
                                               'transcripts': [f'pseudo_transcript_{transcript}'] * len(
                                                   transcriptGroup),
                                               'sequence': pseudo_seqs}, index=pseudo_index)
                    pseudo_lib_tables.append(pseudo_lib)

        else:
            print('generate_pseudogene_dist parameter not recognized, defaulting to off')

        phenotype_table = phenotype_table.append(pd.concat(pseudo_table_list))
        library_table_gene_analysis = library_table[sublib_column].append(pd.concat(pseudo_lib_tables))
    else:
        library_table_gene_analysis = library_table[sublib_column]

    # compute gene scores for replicates, averaged reps, and pseudogenes
    if len(expt_parameters['analyses']) > 0:
        print('Computing gene scores')
        sys.stdout.flush()

        phenotype_table_deduplicated = phenotype_table.loc[
            library_table_gene_analysis.drop_duplicates(['gene', 'sequence']).index]
        if expt_parameters['collapse_to_transcripts']:
            gene_groups = phenotype_table_deduplicated.loc[
                          library_table_gene_analysis.loc[:, 'gene'] != 'negative_control',
                          :].groupby([library_table_gene_analysis['gene'], library_table_gene_analysis['transcripts']])
        else:
            gene_groups = phenotype_table_deduplicated.loc[
                          library_table_gene_analysis.loc[:, 'gene'] != 'negative_control',
                          :].groupby(library_table_gene_analysis['gene'])

        analysis_tables = []
        for analysis in expt_parameters['analyses']:
            print('--' + analysis)
            sys.stdout.flush()

            analysis_tables.append(
                apply_gene_score_function(gene_groups, neg_table, analysis, expt_parameters['analyses'][analysis]))

        gene_table = pd.concat(analysis_tables, axis=1).reorder_levels([1, 2, 0], axis=1).sort_index(axis=1)
        gene_table.to_csv(outbase + '_genetable.txt', sep='\t')

        # collapse the gene-transcript indices into a single score for a gene by best MW p-value, where applicable
        if expt_parameters['collapse_to_transcripts'] and 'calculate_mw' in expt_parameters['analyses']:
            print('Collapsing transcript scores to gene scores')
            sys.stdout.flush()

            gene_table_collapsed = score_gene_by_best_transcript(gene_table)
            gene_table_collapsed.to_csv(outbase + '_genetable_collapsed.txt', sep='\t')

    if generate_plots != 'off':
        if 'calculate_ave' in expt_parameters['analyses'] and 'calculate_mw' in expt_parameters['analyses']:
            temp_data_dict = {'library': library_table[sublib_column],
                              'gene scores': gene_table_collapsed if expt_parameters[
                                  'collapse_to_transcripts'] else gene_table}

            for (phenotype, replicate), gtable in gene_table_collapsed.groupby(level=[0, 1], axis=1):
                if len(replicate_list) == 1 or replicate[:4] == 'ave_':  # just plot averaged reps where available
                    plotting.volcano_plot(temp_data_dict, phenotype, replicate, label_hits=True)

    print('Done!')



def score_gene_by_best_transcript(gene_table):
    """ given a gene table indexed by both gene and transcript, score genes by the best m-w p-value
    per phenotype/replicate"""
    gene_table_trans_groups = gene_table.reorder_levels([2, 0, 1], axis=1)['Mann-Whitney p-value']\
        .reset_index().groupby('gene')

    best_transcript_frame = gene_table_trans_groups.apply(get_best_transcript)

    tup_list = []
    best_trans_list = []
    for tup, group in gene_table.groupby(level=list(range(2)), axis=1):
        tup_list.append(tup)
        cur_frame = gene_table.loc[list(zip(best_transcript_frame.index, best_transcript_frame[tup])), tup]
        best_trans_list.append(cur_frame.reset_index().set_index('gene'))

    return pd.concat(best_trans_list, axis=1, keys=tup_list)


def get_best_transcript(group):
    # set the index to be transcripts and then get the index with the lowest p-value for each cell
    return group.set_index('transcripts').drop(('gene', ''), axis=1).idxmin()


# return Series of counts from a counts file indexed by element id
def read_counts_file(counts_file_name):
    counts_table = pd.read_csv(counts_file_name, header=None, delimiter='\t', names=['id', 'counts'])
    counts_table.index = counts_table['id']
    return counts_table['counts']


# return DataFrame of library features indexed by element id
def read_library_file(library_fasta_file_name, element_type_func, gene_name_func, misc_func_list=None):
    element_list = []
    with open(library_fasta_file_name) as infile:
        id_line = infile.readline()
        while id_line != '':
            seq_line = infile.readline()
            if id_line[0] != '>' or seq_line is None:
                raise ValueError('Error parsing fasta file')

            element_list.append((id_line[1:].strip(), seq_line.strip()))

            id_line = infile.readline()

    element_ids, element_seqs = list(zip(*element_list))
    library_table = pd.DataFrame(np.array(element_seqs), index=np.array(element_ids), columns=['aligned_seq'],
                                 dtype='object')

    library_table['element_type'] = element_type_func(library_table)
    library_table['gene_name'] = gene_name_func(library_table)

    if misc_func_list is not None:
        col_list = [library_table]
        for miscFunc in misc_func_list:
            col_list.append(miscFunc(library_table))
        if len(col_list) != 1:
            library_table = pd.concat(col_list, axis=1)

    return library_table


# print all counts file paths, to assist with making an experiment table
def print_counts_file_paths(base_directory_path_list):
    print('Make a tab-delimited file with the following columns:')
    print('counts_file\texperiment\tcondition\treplicate_id')
    print('and the following list in the counts_file column:')
    for basePath in base_directory_path_list:
        for root, dirs, filenames in os.walk(basePath):
            for filename in fnmatch.filter(filenames, '*.counts'):
                print(os.path.join(root, filename))


def merge_counts_for_experiments(experiment_file_name, library_table):
    expt_table = pd.read_csv(experiment_file_name, delimiter='\t')
    print(expt_table)

    # load in all counts independently
    counts_cols = []
    for countsFile in expt_table['counts_file']:
        counts_cols.append(read_counts_file(countsFile))

    counts_table = pd.concat(counts_cols, axis=1, keys=expt_table['counts_file']).align(library_table, axis=0)[0]

    counts_table = counts_table.fillna(value=0)  # nan values are 0 values, will use nan to filter out elements later

    # print counts_table.head()

    # convert counts columns to experiments, summing when reads across multiple lanes
    expt_tuples = [
        (expt_table.loc[row, 'experiment'], expt_table.loc[row, 'condition'], expt_table.loc[row, 'replicate_id']) for
        row
        in expt_table.index]
    expt_tuples_to_runs = dict()
    for i, tup in enumerate(expt_tuples):
        if tup not in expt_tuples_to_runs:
            expt_tuples_to_runs[tup] = []
        expt_tuples_to_runs[tup].append(expt_table.loc[i, 'counts_file'])

    # print expt_tuples_to_runs

    expt_columns = []
    for tup in sorted(expt_tuples_to_runs.keys()):
        if len(expt_tuples_to_runs[tup]) == 1:
            expt_columns.append(counts_table[expt_tuples_to_runs[tup][0]])
        else:
            column = counts_table[expt_tuples_to_runs[tup][0]]
            for i in range(1, len(expt_tuples_to_runs[tup])):
                column += counts_table[expt_tuples_to_runs[tup][i]]

            expt_columns.append(column)

    # print len(expt_columns), expt_columns[-1]

    expts_table = pd.concat(expt_columns, axis=1, keys=sorted(expt_tuples_to_runs.keys()))
    expts_table.columns = pd.MultiIndex.from_tuples(sorted(expt_tuples_to_runs.keys()))
    # print expts_table

    # mergedTable = pd.concat([libraryTable,counts_table,expts_table],axis=1,
    #                          keys = ['library_properties','raw_counts', 'merged_experiments'])

    return counts_table, expts_table


# filter out reads if /all/ reads for an expt accross replicates/conditions < min_reads
def filter_counts_per_experiment(min_reads, expts_table, library_table):
    experiment_groups = []

    expt_tuples = expts_table.columns

    expt_set = set([tup[0] for tup in expt_tuples])
    for expt in expt_set:
        expt_df = expts_table[[tup for tup in expt_tuples if tup[0] == expt]]
        expt_df_under_min = (expt_df < min_reads).all(axis=1)
        expt_df_filtered = expt_df.align(expt_df_under_min[not expt_df_under_min], axis=0, join='right')[0]
        experiment_groups.append(expt_df_filtered)

        print(expt, len(expt_df_under_min[expt_df_under_min]))

    result_table = pd.concat(experiment_groups, axis=1).align(library_table, axis=0)[0]

    return result_table


# more flexible read filtering
# keep row if either both/all columns are above threshold, or if either/any column is
# in other words, mask if any column is below threshold or only if all columns are below
def filter_low_counts(counts_columns, filter_type, filter_threshold):
    if filter_type == 'both' or filter_type == 'all':
        fail_filter_column = counts_columns.apply(lambda row: min(row) < filter_threshold, axis=1)
    elif filter_type == 'either' or filter_type == 'any':
        fail_filter_column = counts_columns.apply(lambda row: max(row) < filter_threshold, axis=1)
    else:
        raise ValueError('filter type not recognized or not implemented')

    result_table = counts_columns.copy()
    result_table.loc[fail_filter_column, :] = np.nan

    return result_table


def compute_phenotype_score(counts1, counts2, library_table, growth_value, pseudocount_behavior, pseudocount_value,
                            norm_to_negs=True):
    """compute phenotype scores for any given comparison of two conditions"""
    combined_counts = pd.concat([counts1, counts2], axis=1)

    # pseudocount
    if pseudocount_behavior == 'default' or pseudocount_behavior == 'zeros only':
        default_behavior = lambda row: row if min(row) != 0 else row + pseudocount_value
        combined_counts_pseudo = combined_counts.apply(default_behavior, axis=1)
    elif pseudocount_behavior == 'all values':
        combined_counts_pseudo = combined_counts.apply(lambda row: row + pseudocount_value, axis=1)
    elif pseudocount_behavior == 'filter out':
        combined_counts_pseudo = combined_counts.copy()
        zero_rows = combined_counts.apply(lambda row: min(row) <= 0, axis=1)
        combined_counts_pseudo.loc[zero_rows, :] = np.nan
    else:
        raise ValueError('Pseudocount behavior not recognized or not implemented')

    total_counts = combined_counts_pseudo.sum()
    counts_ratio = float(total_counts[0]) / total_counts[1]

    # compute neg control log2 enrichment
    if norm_to_negs:
        neg_counts = combined_counts_pseudo.align(library_table[library_table['gene'] == 'negative_control'],
                                                  axis=0, join='inner')[0]
        # print negCounts
    else:
        neg_counts = combined_counts_pseudo
    neglog2e = neg_counts.apply(calc_log2e, counts_ratio=counts_ratio, growth_value=1, wt_log2_e=0, axis=1).median()
    # print neglog2e

    # compute phenotype scores
    scores = combined_counts_pseudo.apply(calc_log2e, counts_ratio=counts_ratio, growth_value=growth_value,
                                          wt_log2_e=neglog2e, axis=1)

    return scores


def calc_log2e(row, counts_ratio, growth_value, wt_log2_e):
    return (np.log2(counts_ratio * row[1] / row[0]) - wt_log2_e) / growth_value


# average replicate phenotype scores
def average_phenotype_scores(score_table):
    expt_tuples = score_table.columns
    expts_to_replicates = dict()
    for tup in expt_tuples:
        if (tup[0], tup[1]) not in expts_to_replicates:
            expts_to_replicates[(tup[0], tup[1])] = set()
        expts_to_replicates[(tup[0], tup[1])].add(tup[2])

    averaged_columns = []
    labels = []
    for expt in expts_to_replicates:
        expt_df = score_table[[(expt[0], expt[1], rep_id) for rep_id in expts_to_replicates[expt]]]
        averaged_columns.append(expt_df.mean(axis=1))
        labels.append((expt[0], expt[1], 'ave_' + '_'.join(expts_to_replicates[expt])))

    result_table = pd.concat(averaged_columns, axis=1, keys=labels).align(score_table, axis=0)[0]
    result_table.columns = pd.MultiIndex.from_tuples(labels)

    return result_table


def compute_gene_scores(library_table, score_table, norm_to_negs=True):
    gene_groups = score_table.groupby(library_table['gene_name'])

    scored_columns = []
    for expt in score_table.columns:
        if norm_to_negs:
            neg_array = np.ma.array(data=score_table[expt].loc[gene_groups.groups['negative_control']].dropna(),
                                    mask=False)
        else:
            neg_array = np.ma.array(data=score_table[expt].dropna(), mask=False)

        col_list = []
        group_list = []
        for name, group in gene_groups:
            if name == 'negative_control':
                continue
            col_list.append(
                geneStats(group[expt], neg_array))  # group[expt].apply(geneStats, axis = 0, neg_array = neg_array))
            group_list.append(name)

        scored_columns.append(pd.DataFrame(np.array(col_list), index=group_list, columns=['KS', 'KS_sign', 'MW']))

    # return scoredColumns
    return pd.concat(scored_columns, axis=1, keys=score_table.columns)


# apply gene scoring functions to pre-grouped tables of phenotypes
def apply_gene_score_function(grouped_phenotype_table, negative_table, analysis, analysis_param_list):
    if analysis == 'calculate_ave':
        num_to_average = analysis_param_list[0]
        if num_to_average <= 0:
            means = grouped_phenotype_table.aggregate(np.mean)
            counts = grouped_phenotype_table.count()
            result = pd.concat([means, counts], axis=1,
                               keys=['average of all phenotypes', 'average of all phenotypes_sgRNAcount'])
        else:
            means = grouped_phenotype_table.apply(lambda x: average_best_n(x, num_to_average))
            counts = grouped_phenotype_table.count()
            result = pd.concat([means, counts], axis=1,
                               keys=[f'average phenotype of strongest {num_to_average}', 'sgRNA count_avg'])
    elif analysis == 'calculate_mw':
        pvals = grouped_phenotype_table.apply(lambda x: apply_m_w(x, negative_table))
        counts = grouped_phenotype_table.count()
        result = pd.concat([pvals, counts], axis=1, keys=['Mann-Whitney p-value', 'sgRNA count_MW'])
    elif analysis == 'calculate_nth':
        nth = analysis_param_list[0]
        pvals = grouped_phenotype_table.aggregate(
            lambda x: sorted(x, key=abs, reverse=True)[nth - 1] if nth <= len(x) else np.nan)
        counts = grouped_phenotype_table.count()
        result = pd.concat([pvals, counts], axis=1, keys=[f'{nth}th best score', 'sgRNA count_nth best'])
    else:
        raise ValueError('Analysis {analysis} not recognized or not implemented')

    return result


def average_best_n(group, num_to_average):
    return group.apply(lambda column: np.mean(sorted(column.dropna(), key=abs, reverse=True)[:num_to_average]) if len(
        column.dropna()) > 0 else np.nan)


def apply_m_w(group, negative_table):
    if int(sp.__version__.split('.')[1]) >= 17:  # implementation of the "alternative flag":
        return group.apply(lambda column:
                           stats.mannwhitneyu(column.dropna().values, negative_table[column.name].dropna().values,
                                              alternative='two-sided')[1] if len(column.dropna()) > 0 else np.nan)
    else:
        return group.apply(
            lambda column: stats.mannwhitneyu(column.dropna().values, negative_table[column.name].dropna().values)[
                               1] * 2 if len(
                column.dropna()) > 0 else np.nan)  # pre v0.17 stats.mannwhitneyu is one-tailed!!


def parse_gkfile(gk_file_name):
    """parse a tab-delimited file with column headers: experiment, replicate_id, G_value, K_value
       (calculated with martin's parse_growthdata.py)"""
    gkdict = dict()

    with open(gk_file_name, 'rU') as infile:
        for line in infile:
            if line.split('\t')[0] == 'experiment':
                continue
            else:
                linesplit = line.strip().split('\t')
                gkdict[(linesplit[0], linesplit[1])] = (float(linesplit[2]), float(linesplit[3]))

    return gkdict


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate sgRNA- and gene-level phenotypes based on '
                                                 'sequencing read counts, as specified by the experiment config file.')
    parser.add_argument('Config_File', help='Experiment config file specifying screen analysis settings '
                                            '(see accomapnying BLANK and DEMO files).')
    parser.add_argument('Library_File_Directory',
                        help='Directory containing reference library tables and the library_config.txt file.')

    parser.add_argument('--plot_extension', default='png',
                        help='Image extension for plot files, or \"off\". Default is png.')

    args = parser.parse_args()
    # print args

    process_experiments_from_config(args.Config_File, args.Library_File_Directory, args.plot_extension.lower())
