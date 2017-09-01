"""
Functions for interpreting the results in assay_results.tsv

assay_results.tsv:
    rows = metabolites
    cols = cell lines (6 sensitive, 6 resistant)
    final col is the p-value of a t-test between the abundance of the given metabolite
        in sensitive and resistant cell lines


"""

import pandas as pd
import numpy as np
import scipy.stats as st
import math
import itertools
import shutil
import gzip
import os
import sys

base_dir = os.path.dirname(os.path.realpath(__file__))
input_data_dir = base_dir + '/../data/input/'
metadata_dir = base_dir + '/../metadata/'
output_data_dir = base_dir + '/../data/output/'
sys.path += [base_dir + '/../../Rockefeller_Metabolomics']


def select_sif_edges(edge, gzipped_sif, input_data_dir, output_data_dir):
    """
    Filters full pathway commons network file for only the relationships
    characterized by an edge of interest (i.e. 'used-to-produce')
    """
    # i.e. 'PathwayCommons9.All.hgnc.sif.gz'
    sif_path = input_data_dir + gzipped_sif

    # i.e. 'used-to-produce'
    desired_edge = edge
    outf_name = output_data_dir + desired_edge + '.sif'

    # open gzipped sif file
    with gzip.open(sif_path, 'rt') as inf:
        sif = inf.readlines()

    outfh = open(outf_name, "w")

    for line in sif:
        if line == '\n':
            break
        line = line.strip('\n')
        current_edge = line.split('\t')[1]
        if current_edge == desired_edge:
            outfh.write(line + '\n')

    outfh.close()

    # compress your new file
    with open(outf_name, 'rb') as not_zipped, gzip.open(outf_name + '.gz', 'wb') as zipped:
        shutil.copyfileobj(not_zipped, zipped)


# TODO: Verify the log base used to generate these values
# TODO: What to do with the geometric mean?
# TODO: include in create_network?
def add_geometric_cols(assay_file, input_data_dir):
    """
    Calculates 6 vectors:
        geometric means of sensitive cell line metabolite abundance changes
        geometric means of resistant cell line metabolite abundance changes
        geometric means of all cell line metabolite abundance changes
        geometric variance of sensitive cell line metabolite abundance changes
        geometric variance of resistant cell line metabolite abundance changes
        geometric variance of all cell line metabolite abundance changes

    Geometric mean: the anti-log of the arithmetic mean of log-transformed values.
    Geometric variance: the anti-log of the arithmetic variance of log-transformed values.

    Generates assay_results_extended.tsv
    Returns assay results in the form of pd.DataFrame (including geom. means)
    """

    # NOTE: see https://github.com/pandas-dev/pandas/issues/16452
    # In assay_results.tsv, the floats are rounded to 14 digits after the decimal.
    # However this is not how they are stored in memory, so pandas loads the "true"
    # value. I could restrict the number of decimals that are loaded here but that
    # doesn't change how many are actually used in calculations.
    # I've chosen to restrict the decimals only at the point of writing it to file.
    assay_results_path = input_data_dir + assay_file
    assay_results = pd.read_csv(assay_results_path,
                                sep='\t',
                                index_col=0,
                                header=0,
                                na_values = 'nd')
    assay_results.rename(columns={'Unnamed: 13': 'ttest_p'},
                         inplace=True)

    # calculate vector of geometric means of each metabolite sens & resis cell lines
    # the values in assay_results are log-transformed (I'm assuming base is e...)
    for metab in assay_results.index:
        resistant = assay_results.ix[metab,:6]
        sensitive = assay_results.ix[metab,6:12]
        overall = assay_results.ix[metab,:12]

        for count, group in enumerate([resistant, sensitive, overall]):
            arith_mean = np.mean(group)
            arith_var = np.var(group)
            geom_mean = math.e**arith_mean
            geom_var = math.e**arith_var
            if count == 0:
                assay_results.ix[metab, 'resistant_gmean'] = geom_mean
                assay_results.ix[metab, 'resistant_gvar'] = geom_var
            if count == 1:
                assay_results.ix[metab, 'sensitive_gmean'] = geom_mean
                assay_results.ix[metab, 'sensitive_gvar'] = geom_var
            if count == 2:
                assay_results.ix[metab, 'overall_gmean'] = geom_mean
                assay_results.ix[metab, 'overall_gvar'] = geom_var

    # TODO: Determine how to not have it replace digits after the 14th decimal point with zeros
    # note that setting float_format='%.14f' causes sum(), np.sum(), etc. to output "inf" :|
    assay_results.to_csv(assay_results_path.replace('.tsv','_extended.tsv'),
                         sep='\t',
                         #float_format='%.14f',
                         na_rep='NaN')
    return assay_results


def add_arith_mean_cols(assay_results_df, input_dir):
    """
    Adds an arithmetic mean column for each of the sensitive and resistant subsets.
    as well as the overall.

    """
    for metab in assay_results_df.index:
        resistant = assay_results_df.ix[metab, :6]
        sensitive = assay_results_df.ix[metab, 6:12]
        overall = assay_results_df.ix[metab, :12]

        for count, group in enumerate([resistant, sensitive, overall]):
            arith_mean = np.mean(group)
            arith_var = np.var(group)
            if count == 0:
                assay_results_df.ix[metab, 'resistant_amean'] = arith_mean
                assay_results_df.ix[metab, 'resistant_avar'] = arith_var
            if count == 1:
                assay_results_df.ix[metab, 'sensitive_amean'] = arith_mean
                assay_results_df.ix[metab, 'sensitive_avar'] = arith_var
            if count == 2:
                assay_results_df.ix[metab, 'overall_amean'] = arith_mean
                assay_results_df.ix[metab, 'overall_avar'] = arith_var

    assay_results_df.to_csv(input_dir + 'assay_results_extended.tsv',
                         sep='\t',
                         na_rep='NaN')

    return assay_results_df


# TODO: kill this. it doesn't make a lick of sense.
def add_spearmans_corr_col(assay_results_df):
    """
    Takes assay_results pd.DataFrame and calculates the Spearman's rank correlation
    between the sensitive and resistant cell lines for each metabolite.
    Inserts a column with the spearman's rho.
    Inserts a column with the p-value of this correlation.
    Returns the amended pd.DataFrame
    """
    # prepare to collect the spearman's correlation and p value
    sp_corr = []
    sp_pval = []

    # calculate spearman's correlation
    for metab in assay_results_df.index:
        resistant = assay_results_df.ix[metab][:6]
        sensitive = assay_results_df.ix[metab][6:12]
        # calculate spearman's correlation
        # calculate as though NaNs are not present
        spearmanr = st.spearmanr(resistant,
                                 sensitive,
                                 nan_policy='omit')
        sp_corr.append(spearmanr.correlation)
        sp_pval.append(spearmanr.pvalue)

    assay_results_df['spearman_corr'] = sp_corr
    assay_results_df['spearman_pval'] = sp_pval

    assay_results_df.to_csv(input_data_dir + 'assay_results_extended.tsv',
                         sep='\t',
                         #float_format='%.14f',
                         na_rep='NaN')

    return assay_results_df


def add_fc_between_geom_means(assay_results_df):
    """
    Calculates the fold change between the geometric means of the sensitive
    and resistant cell lines for each metabolite. Adds this column to the df.
    Writes out and returns the new dataframe.

    If sens gmean = 0.5 and res gmean = 1.0, FC = 1.0
    if sens gmean = 1.0 and res gmean = 0.5, FC = -0.5
    """
    assay_results_df['fc_gmeans'] = assay_results_df['resistant_gmean']/assay_results_df['sensitive_gmean'] - 1
    assay_results_df.to_csv(input_data_dir + 'assay_results_extended.tsv',
                         sep='\t',
                         na_rep='NaN')

    return assay_results_df


def add_fc_between_arith_means(assay_results_df):
    """
    Calculates the fold change between the geometric means of the sensitive
    and resistant cell lines for each metabolite. Adds this column to the df.
    Writes out and returns the new dataframe.

    If sens gmean = 0.5 and res gmean = 1.0, FC = 1.0
    if sens gmean = 1.0 and res gmean = 0.5, FC = -0.5
    """
    assay_results_df['fc_ameans'] = assay_results_df['resistant_amean']/assay_results_df['sensitive_amean'] - 1
    assay_results_df.to_csv(input_data_dir + 'assay_results_extended.tsv',
                         sep='\t',
                         na_rep='NaN')

    return assay_results_df


def add_ind_ttest_col(assay_results_df):
    """
    Seeing if I can replicate their ttest results.
    (Yes, this is how they did it).
    """
    ttest_statistic = []
    ttest_pval = []
    for metab in assay_results_df.index:
        resistant = assay_results_df.ix[metab][:6]
        sensitive = assay_results_df.ix[metab][6:12]
        ttest = st.ttest_ind(resistant,
                             sensitive,
                             nan_policy='omit')
        # omitting the nans breaks it...
        ttest_statistic.append(ttest.statistic)
        ttest_pval.append(ttest.pvalue)

    assay_results_df['T_stat'] = ttest_statistic
    assay_results_df['ttest_pval'] = ttest_pval

    return assay_results_df


def add_mann_whitney_u_col(assay_results_df):
    """
    Add a Mann Whitney U column
    """

    mwu_statistic = []
    mwu_pval = []
    for metab in assay_results_df.index:
        resistant = assay_results_df.ix[metab][:6]
        sensitive = assay_results_df.ix[metab][6:12]
        mwu = st.mannwhitneyu(resistant,
                             sensitive)
        #                     nan_policy='omit')
        # omitting the nans breaks it...
        mwu_statistic.append(mwu.statistic)
        mwu_pval.append(mwu.pvalue)

    assay_results_df['MWU_stat'] = mwu_statistic
    assay_results_df['MWU_pval'] = mwu_pval

    return assay_results_df


# i.e. all 12 malate FCs vs. all 12 aspartate FCs
# can i do something where it's [6 sens malate vs 6 sens aspartate] vs [6 res malate vs 6 res aspartate]?
def make_pairwise_metab_spearman_matrix(assay_results_df, approach='all'):
    """
    If approach == all: computes spearman correlations between all pairs of profiled metabs
    Correlates two groups of n=12. Writes out 2 matrices: spearman's rho and pval

    Such that the x axis and y axis share the same labels.
        i.e.
                    malate  UDP-hexose  aspartate   ... arginine
        malate        X         X           X               X
        UDP-hexose              X           X               X
        aspartate                           X               X
        ...
        arginine                                            X

    If approach = 'separate': generates 4 matrices. That is, rho and pval matrices are computed
    separately for sensitive and resistant cell lines.
    """
    # isolate the fold change values
    assay_results_df = assay_results_df.ix[:,:12]

    # get an itertools.combinations object of metabolite pairs
    pairs = itertools.combinations(assay_results_df.index, 2)

    if approach=='all':

        # initialize empty dataframes
        corr_df = pd.DataFrame(index=assay_results_df.index,
                               columns=assay_results_df.index)

        pval_df = pd.DataFrame(index=assay_results_df.index,
                               columns=assay_results_df.index)

        for pair in pairs:
            metab1 = assay_results_df.ix[pair[0]][:12]
            metab2 = assay_results_df.ix[pair[1]][:12]

            # calculate spearman's correlation
            # calculate as though NaNs are not present
            spearmanr = st.spearmanr(metab1,
                                     metab2,
                                     nan_policy='omit')

            corr_df.ix[pair[0], pair[1]] = spearmanr.correlation
            pval_df.ix[pair[0], pair[1]] = spearmanr.pvalue

            corr_df.to_csv(input_data_dir + 'pairwise_spearman_corr_coef.tsv',
                           sep='\t',
                           na_rep='NaN')
            pval_df.to_csv(input_data_dir + 'pairwise_spearman_pval.tsv',
                           sep='\t',
                           na_rep='NaN')

    if approach=='separate':

        # initialize empty dataframes
        res_corr_df = pd.DataFrame(index=assay_results_df.index,
                               columns=assay_results_df.index)

        res_pval_df = pd.DataFrame(index=assay_results_df.index,
                               columns=assay_results_df.index)

        sens_corr_df = pd.DataFrame(index=assay_results_df.index,
                               columns=assay_results_df.index)

        sens_pval_df = pd.DataFrame(index=assay_results_df.index,
                               columns=assay_results_df.index)

        for pair in pairs:

            res_metab1 = assay_results_df.ix[pair[0]][:6]
            res_metab2 = assay_results_df.ix[pair[1]][:6]

            sens_metab1 = assay_results_df.ix[pair[0]][6:12]
            sens_metab2 = assay_results_df.ix[pair[1]][6:12]

            res_spearmanr = st.spearmanr(res_metab1,
                                          res_metab2,
                                          nan_policy='omit')

            sens_spearmanr = st.spearmanr(sens_metab1,
                                         sens_metab2,
                                         nan_policy='omit')

            res_corr_df.ix[pair[0], pair[1]] = res_spearmanr.correlation
            res_pval_df.ix[pair[0], pair[1]] = res_spearmanr.pvalue

            sens_corr_df.ix[pair[0], pair[1]] = sens_spearmanr.correlation
            sens_pval_df.ix[pair[0], pair[1]] = sens_spearmanr.pvalue

            res_corr_df.to_csv(input_data_dir + 'resistant_pairwise_spearman_corr_coef.tsv',
                           sep='\t',
                           na_rep='NaN')
            res_pval_df.to_csv(input_data_dir + 'resistant_pairwise_spearman_pval.tsv',
                           sep='\t',
                           na_rep='NaN')

            sens_corr_df.to_csv(input_data_dir + 'sensitive_pairwise_spearman_corr_coef.tsv',
                           sep='\t',
                           na_rep='NaN')
            sens_pval_df.to_csv(input_data_dir + 'sensitive_pairwise_spearman_pval.tsv',
                           sep='\t',
                           na_rep='NaN')


def ln_assay_results(assay_file, input_dir):
    """
    Produces an output file which contains natural log values of all original assay results
    Returns the log'd data frame
    """
    assay_results_path = input_dir + assay_file
    assay_results = pd.read_csv(assay_results_path,
                                sep='\t',
                                index_col=0,
                                header=0,
                                na_values = 'nd')
    log_fcs = np.log(assay_results.ix[:,:12])
    log_fcs.to_csv(input_dir + 'assay_results_log.tsv',
                   sep='\t',
                   na_rep='NaN')

    return log_fcs
