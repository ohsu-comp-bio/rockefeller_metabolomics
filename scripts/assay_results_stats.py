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
import math

# TODO: Verify the log base used to generate these values
def get_geometric_mean(assay_file, input_data_dir):
    """
    Calculates a 2 vectors:
        1 vector of geometric means of sensitive cell line metabolite abundance changes
        1 vector of geometric means of resistant cell line metabolite abundance changes

    Geometric mean: the anti-log of the arithmetic mean of log-transformed values.

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
        sensitive = assay_results.ix[metab,:6]
        resistant = assay_results.ix[metab,6:12]
        s_arith_mean = np.mean(sensitive)
        r_arith_mean = np.mean(resistant)
        s_geom_mean = math.e**s_arith_mean
        r_geom_mean = math.e**r_arith_mean
        assay_results.ix[metab, 'sensitive_gmean'] = s_geom_mean
        assay_results.ix[metab, 'resistant_gmean'] = r_geom_mean

    # TODO: Determine how to not have it replace digits after the 14th decimal point with zeros
    assay_results.to_csv(assay_results_path.replace('.tsv','_extended.tsv'),
                         sep='\t',
                         float_format='%.14f',
                         na_rep='NaN')

    return assay_results