"""
Calls on other functions to generate two networks of small molecule
interactions relevant to a user-specified matrix of log-transformed,
fold changes in the abundance of small molecules in (1) sensitive and
(2) resistant cell lines following pathway perturbation.

Example commands:
    python create_networks.py
    python create_networks.py -c1 0.2 -c2 0.5 -p .10


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
import argparse
import glob
import sys
import os

base_dir = os.path.dirname(os.path.realpath(__file__))
input_dir = base_dir + '/../data/input/'
metadata_dir = base_dir + '/../metadata/'
sys.path += [base_dir + '/../../Rockefeller_Metabolomics']


def main():
    # take user-specified names of files
    parser = argparse.ArgumentParser()

    parser.add_argument("-s", "--sif",
                        type=str,
                        default='used-to-produce.sif',
                        help="Name of used-to-produce network file in data/input/\n"
                             "(default = 'used-to-produce.sif')")
    parser.add_argument("-a", "--assay_file",
                        type=str,
                        default='assay_results_extended.tsv',
                        help="Name of assay results file\n"
                             "(default = 'assay_results_extended.tsv')")
    parser.add_argument("-c1", "--minimum_corr_dist1",
                        type=float,
                        default=0.5,
                        help="Min absolute value of spearman correlation\n"
                             "between two nodes in order for them\n"
                             "to be kept in 'both' relationships\n"
                             "default = 0.5")
    parser.add_argument("-c2", "--minimum_corr_dist2",
                        type=float,
                        default=0.75,
                        help="Min absolute value of spearman correlation\n"
                             "between two nodes in order for them\n"
                             "to be kept in 'distance of 2' relationships\n"
                             "default = 0.5")
    parser.add_argument("-ll", "--linker_lenience",
                        type=int,
                        default=4,
                        help="Max number of .sif relationships in which an unprofiled\n"
                             "'linker' node may be present in order to be included\n"
                             "in the distance-of-2 network\n"
                             "default = 4")
    parser.add_argument("-p", "--maximum_mwu_p",
                        type=float,
                        default=0.10,
                        help="The maximum p value allowed for each metabolite\n"
                             "(mann whitney u test between sensitive and resistant)\n"
                             "(default = 0.10)")

    args = parser.parse_args()
    sif_name = args.sif
    assay_file = args.assay_file
    min_corr_both = args.minimum_corr_dist1
    min_corr_dist2 = args.minimum_corr_dist2
    ll = args.linker_lenience
    max_p = args.maximum_mwu_p

    # generate the run's informative title
    output_dir = base_dir + '/../data/output/'
    (output_dir, run_tag) = name_outfile_path(output_dir, min_corr_both, min_corr_dist2, max_p)

    # make the new output dir if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # read in the extended assay results dataframe
    assay_results_path = input_dir + assay_file
    assay_df = pd.read_csv(assay_results_path,
                                sep='\t',
                                index_col=0,
                                header=0,
                                na_values = 'nd')

    # get the list of sufficiently altered metabs as per Mann Whitney
    signif_metab_names = assay_df[assay_df['MWU_pval'] <= max_p].index.tolist()

    # build dictionary {'metabolite': ['CHEBI:1234', 'CHEBI:4567', ... ]
    name_to_chebis_map = make_name_chebi_map(signif_metab_names,
                                             input_dir,
                                             metadata_dir)

    # get list of all chebis either (1) exactly representing or (2) derived from profiled metabolites:
    chebis = [chebi for derived_group in list(name_to_chebis_map.values()) for chebi in derived_group]


    # filter used-to-produce sif into BOTH, EITHER, and DISTANCE-OF-2 sifs for all, sensitive, and resistant
    # data (ChEBI IDs only)
    # TODO: in "all" filtration, does it makes sense to use the fc of ameans of fcs? MWU?
    ch_paths = filter_sif_by_chebi(sif_name,
                                       chebis,
                                       name_to_chebis_map,
                                       min_corr_both,
                                       min_corr_dist2,
                                       ll,
                                       input_dir,
                                       output_dir,
                                       run_tag)

    # convert ChEBI sifs to named sifs
    named_paths = []
    for ch_path in ch_paths:
        named_paths.append(convert_chebi_sif_to_named_sif(ch_path, chebis, name_to_chebis_map, input_dir))

    # specify formatting for the named sifs
    for named_path in named_paths:
        if 'all' in named_path:
            specify_named_chibe_formatting_overall(signif_metab_names, named_path)
        if 'res' in named_path:
            specify_named_chibe_formatting_subset(named_path, 'resistant')
        if 'sens' in named_path:
            specify_named_chibe_formatting_subset(named_path, 'sensitive')

    # remove the chebi sifs
    for filename in glob.glob(output_dir + "*chebi*"):
        os.remove(filename)

if __name__ == '__main__':
    main()