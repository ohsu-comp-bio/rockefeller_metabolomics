"""
Functions for generating boxplots.

Author: Hannah Manning
Date: August 30, 2017
"""

import itertools
import numpy as np
import pandas as pd

def generate_metab_boxplots_all(assay_df, input_dir, restrict_y=False, log_e=False):
    """
    Generates 112 boxplots side by side and saves them in data/plots.
    Can restrict the ylim by setting restrict_y to desired tuple (i.e. (0,10))
    """

    assay_df_profiled = assay_df.ix[:,:12]
    assay_df_profiled = assay_df_profiled.T

    base_outf_name = input_dir.replace('input', 'plots') + 'rawfc_boxplots_'
    ylabel = 'Fold change'

    if log_e:
        base_outf_name = base_outf_name.replace('raw', 'log')
        assay_df_profiled = np.log(assay_df_profiled)
        ylabel = 'Log(fold change)'

    if restrict_y:
        base_outf_name = base_outf_name.replace('boxplots',
                                                'boxplots_' +
                                                str(restrict_y[0]) +
                                                '_' + str(restrict_y[1]) + 'ylim')

    res_lines = assay_df.columns[:6]
    sens_lines = assay_df.columns[6:12]

    chunks = []
    prev_i = 0
    # this isn't doing it right
    for i in range(16,128,16):
        chunks.append(assay_df_profiled.ix[:,prev_i:i])
        prev_i = i

    for chunk_count, chunk in enumerate(chunks):
        fig = plt.figure(figsize=(14, 8))
        chunk.boxplot(showfliers=False)
        ax = fig.add_subplot(111)
        ax.grid(False)
        axes = plt.gca()
        plt.xticks(rotation='45')
        plt.ylabel(ylabel)
        plt.title('All measurements for metabs: ' + chunk.columns[0] + ' through ' + chunk.columns[-1])
        # important that this goes at the end

        # add points
        colcount = 1
        for metab in chunk.columns:
            red_points = chunk[metab][res_lines]
            cyan_points = chunk[metab][sens_lines]
            colors = itertools.cycle(["r.", "c."])

            for colorgroup in [red_points, cyan_points]:
                # prepare to add jitter
                x = np.random.normal(colcount, 0.08, len(colorgroup))
                plt.plot(x, colorgroup, next(colors), alpha=0.6)
            colcount += 1

        if restrict_y:
            plt.ylim(i for i in restrict_y)

        plt.tight_layout()
        plt.savefig(base_outf_name + str(chunk_count) +'.png')
        plt.close(fig)


def generate_cell_line_boxplot(assay_df, input_dir, restrict_y = False, log_e=False):
    """
    Generates 1 plot containing 12 boxes and saves it to data/plots/
    Each dot is 1 metabolite's fc in that cell line

    ylim can be set with restrict_y
        restrict_y = (0,5)
    """

    assay_df_profiled = assay_df.ix[:,:12]
    outfile_name = input_dir.replace('input','plots') + 'rawfc_boxplot_by_cell_line.png'

    if log_e:
        assay_df_profiled = np.log(assay_df_profiled)
        outfile_name = outfile_name.replace('raw', 'log')

    res_df = assay_df.ix[:6]
    sens_df = assay_df.ix[6:12]

    fig = plt.figure(figsize=(14, 8))
    bp = assay_df_profiled.boxplot(showfliers=False)
    ax = fig.add_subplot(111)
    ax.grid(False)
    axes = plt.gca()
    if not log_e:
        plt.title('Raw FCs by cell line')
        plt.ylabel('Fold change')
    if log_e:
        plt.title('Natural log FCs by cell line')
        plt.ylabel('Log(fold change)')

    colcount = 1
    for cell_line in assay_df_profiled.columns:
        points = assay_df_profiled[cell_line]
        x = np.random.normal(colcount, 0.08, len(points))
        plt.plot(x, points, "b.", alpha=0.6)
        colcount += 1

    if restrict_y:
        plt.ylim(i for i in restrict_y)
        plt.tight_layout()
        outfile_name = outfile_name.replace('.png', '_ylim_' + str(restrict_y[0]) + '_' + str(restrict_y[1]) + '.png')
        plt.savefig(outfile_name)
        plt.close(fig)

    if not restrict_y:
        plt.tight_layout()
        plt.savefig(outfile_name)
        plt.close(fig)
