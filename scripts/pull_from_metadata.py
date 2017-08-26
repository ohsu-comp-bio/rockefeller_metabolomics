"""
Functions that extract data from the user-provided metadata file.

Author: Hannah Manning
Date: August 21, 2017
"""

from metabolite_mapping import *
import sys

base_dir = os.path.dirname(os.path.realpath(__file__))
input_data_dir = base_dir + '/../data/input/'
metadata_dir = base_dir + '/../metadata/'
output_data_dir = base_dir + '/../data/output/'
sys.path += [base_dir + '/../../Rockefeller_Metabolomics']

# TODO: not currently used in create_network.py
def collect_all_chebis_being_used(metadata_dir):
    """
    Another slapped-together function that will serve a purpose right now
    but maybe not later.
    Gets CHEBI IDs from exact_matches.tsv and from user_added_chebis_ALL.txt
    Writes them out to all_chebis_of_interest.txt.

    AND MAPS THEM TO THEIR NON-CHEBI-ID NAMES.
    """
    meta_path = metadata_dir + 'meta_data_user_specified_chebis.txt'
    meta_file = open(meta_path, "r")
    meta = meta_file.readlines()

    exact_matches_path = output_data_dir + 'exact_matches.tsv'
    exact_matches_file = open(exact_matches_path, "r")
    exact_matches = exact_matches_file.readlines()

    all_chebis_of_interest = open(output_data_dir + 'all_chebis_of_interest.txt', "w")

    for line in exact_matches:
        line = line.strip('\n')
        parts = line.split('\t')
        name = parts[0]
        chebi_id = parts[1]
        all_chebis_of_interest.write(chebi_id + '\t' + name + '\n')

    exact_matches_file.close()

    for line in meta:
        if line.startswith('CHEBI:'):
            parts = line.split('\t')
            chebi_id = parts[0]
            name = parts[2]
            all_chebis_of_interest.write(chebi_id + '\t' + name + '\n')

    user_added_file.close()
    meta_file.close()
    all_chebis_of_interest.close()

def group_chebis_of_same_parent():
    """
    Maps all ChEBI IDs being used to the original metabolite of interest
    provided in assay_results.tsv.
    Returns a dictionary where:
        key = named metabolite
                (i.e. 'aspartate')
        value = list of derived ChEBI IDs
                (i.e. ['CHEBI:29990', 'CHEBI:29991', 'CHEBI:29995']

    meta_data_path = meta_data_dir + meta_data_user_specified_chebis.txt
    """
    # get all metabolites (by name) from assay_file
    metabs = get_metabolites('assay_results_extended.tsv', input_data_dir)

    # initialize the dictionary
    metab_dict = {}
    for i in metabs:
        metab_dict[i] = []

    meta_data = pd.read_csv(metadata_dir + 'meta_data_user_specified_chebis.txt',
                            sep='\t',
                            comment='#',
                            index_col=False
                            )

    # parse metadata to add ChEBI IDs to dictionary as values
    for metab in metab_dict:
        metab_rows = meta_data[meta_data['METAB_OF_INTEREST'] == metab]
        for i in metab_rows.index:
            metab_dict[metab].append(metab_rows['ADDED_CHEBI_ID'][i])

    # collect false alarm exact matches
    failed_exact_fh = open(output_data_dir + 'exact_matches_not_in_sif.tsv', 'r')
    failed_exact = failed_exact_fh.readlines()

    failed = []
    for line in failed_exact:
        line = line.strip('\n')
        line = line.split('\t')
        failed.append(line[1])

    failed_exact_fh.close()

    exact_matches_fh = open(output_data_dir + 'exact_matches.tsv', 'r')
    exact_matches = exact_matches_fh.readlines()

    for line in exact_matches:
        line = line.strip('\n')
        line = line.split('\t')

        name = line[0]
        cheb = line[1]

        if cheb not in failed:
            metab_dict[name].append(cheb)

    exact_matches_fh.close()

    return metab_dict


def append_user_added_chebis_to_map(name_to_chebis_dict, path_to_addl_chebis):
    """
    Takes a dictionary with profiled metabolite names as keys (i.e. 'malate', 'aspartate')
    and lists of derived ChEBI IDs as values (i.e. malate: [CHEBI:1234, CHEBI:4567])

    Adds user-curated ChEBI IDs to the value lists. These are drawn path_to_addl_chebis.

    Returns enhanced dictionary.
    """
    # load metadata
    addl_chebis = path_to_addl_chebis
    metadata = pd.read_csv(path_to_addl_chebis,
                           sep='\t',
                           comment='#',
                           usecols=[0, 1],
                           index_col=0)

    # get ChEBI IDs relevant to metab of interest:
    for name in name_to_chebis_dict:
        [name_to_chebis_dict[name].append(chebi) for chebi in metadata[metadata['METAB_OF_INTEREST'] == name].index]

    return name_to_chebis_dict