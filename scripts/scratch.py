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
        sensitive = assay_results.ix[metab,:6]
        resistant = assay_results.ix[metab,6:12]
        overall = assay_results.ix[metab,:12]

        for count, group in enumerate([sensitive, resistant, overall]):
            arith_mean = np.mean(group)
            arith_var = np.var(group)
            geom_mean = math.e**arith_mean
            geom_var = math.e**arith_var
            if count == 0:
                assay_results.ix[metab, 'sensitive_gmean'] = geom_mean
                assay_results.ix[metab, 'sensitive_gvar'] = geom_var
            if count == 1:
                assay_results.ix[metab, 'resistant_gmean'] = geom_mean
                assay_results.ix[metab, 'resistant_gvar'] = geom_var
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