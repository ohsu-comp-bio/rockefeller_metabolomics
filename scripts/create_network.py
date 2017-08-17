"""
Calls on other functions to generate a network of small molecule
interactions relevant to a user-specified matrix of log-transformed,
fold changes in the abundance of small molecules in sensitive and
resistant cell lines following pathway perturbation.

Example command:
    python create_network.py -e used-to-produce
        -s PathwayCommons9.All.hgnc.sif.gz
        -a assay_results.tsv

Author: Hannah Manning <manningh@ohsu.edu>
Date: August 8, 2017
"""

from filter_sif_edge import *
from metabolite_mapping import *
from sif_chebi_ids import *
from user_input_chebis import *
from generate_filtered_sif import *
from assay_results_stats import *

import argparse
import sys
import os

base_dir = os.path.dirname(os.path.realpath(__file__))
input_data_dir = base_dir + '/../data/input/'
metadata_dir = base_dir + '/../metadata/'
output_data_dir = base_dir + '/../data/output/'
sys.path += [base_dir + '/../../Rockefeller_Metabolomics']

def main():
    # take user-specified names of files
    parser = argparse.ArgumentParser()

    parser.add_argument("-e", "--edge",
                        type=str,
                        help="PathwayCommons edge (i.e. 'used-to-produce')")
    parser.add_argument("-s", "--sif",
                        type=str,
                        help="Name of gzipped PathwayCommons sif file in data/\n"
                        "(i.e. 'PathwayCommons9.All.hgnc.sif.gz'")
    parser.add_argument("-a", "--assay_file",
                        type=str,
                        help="Name of assay results file (i.e. 'assay_results.tsv')")

    args = parser.parse_args()
    edge = args.edge
    sif = args.sif
    assay_data = args.assay_file

    print("Paring greater .sif file by user-specified edge: " + edge)
    select_sif_edges(edge, sif, input_data_dir, output_data_dir)
    metabs = get_metabolites(assay_data, input_data_dir, output_data_dir)

    [signif_chebis, signif_no_chebi_match] = get_significant_metabs(assay_data, input_data_dir)

    if len(signif_no_chebi_match) > 0:
        print("Exact ChEBI ID matches could not be found for the following " +
              str(len(signif_no_chebi_match)) + " metabolites:\n")
        for name in signif_no_chebi_match:
            print(name)

    print("Identifying significant IDs with no exact ChEBI ID match...")
    signif_no_chebi_match = map_all_metabs_to_chebi_ids(metabs, signif_no_chebi_match, output_data_dir)

    print("Verifying presence of exact ChEBI ID matches in " + edge + ".sif...")
    misfit_loc = verify_exact_chebis('used-to-produce.sif', output_data_dir)

    # ask user to supply an file with additional ChEBI IDs to include in network
    path_to_addl_chebis = ask_user_for_chebis(misfit_loc, signif_no_chebi_match, input_data_dir)

    print("Filtering " + edge + ".sif by ChEBI IDs of interest...")
    [both_sif_path, either_sif_path, dist_2_path] = filter_sif_by_chebi(edge,
                                                           path_to_addl_chebis,
                                                           output_data_dir)

    print("Specifying .sif format for ChiBE visualization...")
    specify_chibe_formatting(signif_chebis, both_sif_path)
    specify_chibe_formatting(signif_chebis, either_sif_path)
    specify_chibe_formatting(signif_chebis, dist_2_path)

    print("Converting ChEBI IDs in BOTH .sif...")
    convert_chebi_sif_to_named_sif(both_sif_path, input_data_dir, output_data_dir)

    print("Converting ChEBI IDs in EITHER .sif...")
    convert_chebi_sif_to_named_sif(either_sif_path, input_data_dir, output_data_dir)


if __name__ == '__main__':
    main()