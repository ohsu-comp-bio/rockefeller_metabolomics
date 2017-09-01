def generate_drug_mut_assoc_bp(auc_df, anova_df, clf_type, plots_dir, id, run):
    """
    Generates a figure with 1 boxplot per high-quality classifier.
    Data points for each boxplot are drug-mutation associations from auc_df.
    They are colored according to whether their are expected to be
    positively correlated (blue), negatively correlated (green),
    or not correlated (black) as per the ANOVA scores in anova_df.
    auc_df and anova_df should bear the same shape:
        cols = mutations
        rows = drugs
    """

    # get the sign ("direction") of the ANOVA score in iorio_assoc
    # this will be used to determine the color of datapoints in boxplots
    assoc_directions = {drug: {'pos_anova': [], 'neg_anova': []} for drug in anova_df.columns}
    for drug in anova_df.columns:
        row_indexer = 0
        for anova_score in anova_df[drug]:
            if anova_score < 0:
                mut = anova_df.index[row_indexer]
                assoc_directions[drug]['neg_anova'].append(mut)
            if anova_score > 0:
                mut = anova_df.index[row_indexer]
                assoc_directions[drug]['pos_anova'].append(mut)
            row_indexer += 1

    # generate boxplots of tcga_auc
    fig = plt.figure(figsize=(14, 8))
    bp = auc_df.boxplot(showfliers=False)
    ax = fig.add_subplot(111)
    ax.grid(False)
    axes = plt.gca()
    plt.title(clf_type + "classifier behavior")
    plt.ylabel("AUC")
    plt.xlabel("Drug Classifier")
    colcount = 0
    for drug in auc_df.columns:
        colcount += 1

        # lists muts
        pos_anova_muts = assoc_directions[drug]['pos_anova']  # color these blue
        neg_anova_muts = assoc_directions[drug]['neg_anova']  # color these red
        zero_anova_muts = auc_df[drug].index.drop(pos_anova_muts).drop(neg_anova_muts)

        # get the drug column (type is pd.Series)
        aucs = auc_df[drug]

        # generate series for different data colors (auc values)
        blue_points = aucs[pos_anova_muts]
        red_points = aucs[neg_anova_muts]
        black_points = aucs[zero_anova_muts]
        colors = itertools.cycle(["r.", "k.", "b."])
        point_size = itertools.cycle([9.0, 3.0, 9.0])

        for colorgroup in [red_points, black_points, blue_points]:
            # prepare to add jitter
            x = np.random.normal(colcount, 0.08, len(colorgroup))
            plt.plot(x, colorgroup, next(colors), alpha=0.6, markersize=next(point_size))
        plt.xticks(rotation='45')

    plt.axhline(y=0.50, c="0.75")
    axes.set_ylim([0.0, 1.0])
    axes.set_yticks(np.arange(0, 1.1, 0.1))
    plt.tight_layout()
    # plt.show(block=False)
    # plt.waitforbuttonpress(0)
    # print("Press any button to close the plots and proceed.")
    plt.savefig(plots_dir + clf_type + '/' + run + '/'
                + id + '_' + clf_type + '_' + run + '_behavior_boxplots.png')
    plt.close(fig)


def generate_ccle_resp_bp(pred_ccle_resp, patient_resp, clf_type, plots_dir, id, run):
    """
    Generates boxplots of (1) predicted and (2) actual cell line responses
    to each drug. Draws a line representative of the predicted patient/sample
    response for comparison.
    pred_ccle_resp (dict):
        key (str): drug name (only those for well behaved classifiers)
        value (pd.Series): vector of cell line predicted responses to drug-key
    patient_resp (pd.Series):
        key (str): drug name
        value (float): predicted patient response to drug-key
    clf_type (str):
        i.e. 'ElasticNet', 'rForest', or 'SVRrbf'
    prefix (str):
        user-specified prefix for output plot file names
    """

    # load CCLE response data relevant to the well-behaving drug classifiers
    actual_ccle_resp = get_drug_ioria(pred_ccle_resp.keys())

    # for each drug make a separate boxplot
    for drug in pred_ccle_resp.keys():
        # TODO: make a line for patient response

        # make a little pandas dataframe
        # TODO: refactor this. i rushed and i'm sure there's a better way.
        single_drug_df = pd.DataFrame({'Measured': actual_ccle_resp[drug].values})
        single_drug_df['Predicted'] = pd.Series(pred_ccle_resp[drug].values)

        fig = plt.figure(figsize=(8, 8))
        bp = single_drug_df.boxplot(showfliers=False, vert=False)
        ax = fig.add_subplot(111)
        ax.grid(False)
        axes = plt.gca()
        plt.title(drug + " responses: actual and predicted by the " + clf_type + " classifier")
        plt.xlabel("AUC")

        # add the points with jitter
        colcount = 0
        for colname in single_drug_df:
            colcount += 1
            y = np.random.normal(colcount, 0.08, single_drug_df.shape[0])
            plt.plot(single_drug_df[colname], y, "c.", alpha=0.2)

        plt.axvline(x=0.50, c="0.75")
        axes.set_xlim([0.0, 1.1])
        axes.set_xticks(np.arange(0, 1.1, 0.1))
        plt.tight_layout()

        plt.savefig(plots_dir + clf_type + '/' + run + '/'
            + id + '_' + clf_type + '_' + drug.replace('/', '__').replace(' ', '_')
            + '_' + run + '_ccle_response_bp.png'
            )
        plt.close(fig)
