"""
Calls on other functions to generate two networks of small molecule
interactions relevant to a user-specified matrix of log-transformed,
fold changes in the abundance of small molecules in (1) sensitive and
(2) resistant cell lines following pathway perturbation.

Example command:
    python create_sens_res_networks.py -s used-to-produce.sif -a assay_results_extended.tsv

Author: Hannah Manning <manningh@ohsu.edu>
Date: August 24, 2017
"""

from filter_sif import *
from metabolite_mapping import *
from sif_chebi_ids import *
from user_input_chebis import *
from assay_results_stats import *
from format_sif import *
from pull_from_metadata import *
from utils import *

import pandas as pd
import numpy as np
import argparse
import sys
import os

base_dir = os.path.dirname(os.path.realpath(__file__))
input_dir = base_dir + '/../data/input/'
metadata_dir = base_dir + '/../metadata/'
output_dir = base_dir + '/../data/output/'
sys.path += [base_dir + '/../../Rockefeller_Metabolomics']


def main():
    # take user-specified names of files
    parser = argparse.ArgumentParser()

    parser.add_argument("-s", "--sif",
                        type=str,
                        help="Name of used-to-produce network file in data/input/\n"
                             "(i.e. 'used-to-produce.sif')")
    parser.add_argument("-a", "--assay_file",
                        type=str,
                        help="Name of assay results file\n"
                             "(i.e. 'assay_results_extended.tsv')")
    parser.add_argument("-c", "--minimum_correlation",
                        type=str,
                        default=0.5,
                        help="Min absolute value of spearman correlation\n"
                             "between two nodes in order for them\n"
                             "to be kept in resulting networks\n"
                             "default = 0.5")
    parser.add_argument("-l", "--linker_lenience",
                        default=4,
                        type=str,
                        help="Max number of .sif relatonships in which an unprofiled\n"
                             "'linker' node may be present in order to be included\n"
                             "in the distance-of-2 network\n"
                             "default = 4")

    args = parser.parse_args()
    sif = args.sif
    assay_file = args.assay_file
    min_corr = args.minimum_correlation
    ll = args.linker_lenience

    # read in the assay results dataframe
    # note that it already contains the following added columns:
    # sensitive_gmean, sensitive_gvar, resistant_gmean, resistant_gvar,
    # overall_gmean, overall_gvar, spearman_corr, spearman_pval, fc_gmeans
    assay_results_path = input_dir + assay_file
    assay_df = pd.read_csv(assay_results_path,
                                sep='\t',
                                index_col=0,
                                header=0,
                                na_values = 'nd')

    # get the list of "significantly altered" metabs as per the t-test
    signif_metab_names = assay_df[assay_df['ttest_p'] <= 0.05].index.tolist()

    # build dictionary {'metabolite': ['CHEBI:1234', 'CHEBI:4567', ... ]
    name_to_chebis_map = make_name_chebi_map(assay_df,
                                             input_dir,
                                             metadata_dir)

    # get list of all chebis either exactly representing to or derived from profiled metabolites:
    chebis = [chebi for derived_group in list(name_to_chebis_map.values()) for chebi in derived_group]


    # filter used-to-produce sif into BOTH, EITHER, and DISTANCE-OF-2 sifs for all, sensitive, and resistant
    # data
    [all_output_locations, sens_output_locations, res_output_locations] = filter_sif_by_chebi(chebis,
                                                                        name_to_chebis_map,
                                                                        min_corr,
                                                                        ll,
                                                                        input_dir,
                                                                        output_dir)
    # TODO: REPEAT THIS FOR BOTH AND EITHER SIFS. BEING LAZY AND ONLY DOING DIST2 SIFS FOR NOW
    all_dist2_chebi_location = all_output_locations[2]
    sens_dist2_chebi_location = sens_output_locations[2]
    res_dist2_chebi_location = res_output_locations[2]

    # convert the all chebi sif to the named version
    all_dist2_named_loc = convert_chebi_sif_to_named_sif(all_dist2_chebi_location, chebis, name_to_chebis_map, input_dir)
    sens_dist2_named_loc = convert_chebi_sif_to_named_sif(sens_dist2_chebi_location, chebis, name_to_chebis_map, input_dir)
    res_dist2_named_loc = convert_chebi_sif_to_named_sif(res_dist2_chebi_location, chebis, name_to_chebis_map, input_dir)

    # todo: make sure you repeat this for the non-dist2 ones as well
    # remove weakly correlated from the 'all' one only (because that data's not available in separate approach)
    all_dist2_named_loc = remove_weakly_correlated(all_dist2_named_loc, input_dir)

    # specify formatting for the named all sif
    all_dist2_format_path = specify_named_chibe_formatting_overall(signif_metab_names, all_dist2_named_loc)

    # a different approach is needed to format the sensitive and resistant sifs
    sens_dist2_format_path = specify_named_chibe_formatting_subset(sens_dist2_named_loc, 'sensitive')
    res_dist2_format_path = specify_named_chibe_formatting_subset(res_dist2_named_loc, 'resistant')



if __name__ == '__main__':
    main()