"""
Functions for
    (1) retrieving metabolites from a matrix of assay results,
    (2) selecting significantly altered metabolites
    (3) mapping metabolite names to their ChEBI IDs


Author: Hannah Manning <manningh@ohsu.edu>
Date: August 8, 2017
"""


import libchebipy as lc
import pandas as pd
import os
from utils import *
import pull_from_metadata

def make_name_chebi_map(assay_df,
                        input_dir,
                        metadata_dir):
    """
    Parses assay results dataframe to generate a dictionary with profiled metabolite names
    as keys. The value for each metabolite name key is a list of ChEBI IDs that were derived
    from that name (i.e. the exact ChEBI match and any isomers, conjugate acids/bases, etc.)

    Returns this dictionary.
    """
    # Initialize a dictionary to contain compound names as keys and all derived
    # ChEBI IDs as values (in a list)
    name_to_chebis_map = {metab: [] for metab in assay_df.index.tolist()}

    # read in the compounds dataframe for mapping names to ChEBI IDs
    compounds_df = pd.read_csv(input_dir + 'compounds.tsv',
                               sep='\t',
                               usecols=[2, 5],
                               header=0,
                               index_col=0)

    # search for exact matches in compounds.tsv
    for metab in name_to_chebis_map:
        exact_chebi_list = compounds_df[compounds_df['NAME'] == metab].index.tolist()
        if len(exact_chebi_list) > 0:
            name_to_chebis_map[metab].append(exact_chebi_list[0])

    # set location of data:
    path_to_addl_chebis = metadata_dir + 'meta_data_user_specified_chebis.txt'

    # add metadata chebis to this map:
    name_to_chebis_map = pull_from_metadata.append_user_added_chebis_to_map(name_to_chebis_map, path_to_addl_chebis)

    return name_to_chebis_map


def get_significant_metab_chebis(assay_file, input_data_dir):
    """
    Returns a 2 lists of ChEBI IDs for metabolites whose fold change in abundance achieved
    significance at the p < 0.05 level according to a t-test between sensitive and
    resistant strains. If exact matches are available in the ChEBI database, those ChEBI IDs
    are entered into list 1. If no matches are found, the names are entered into list 2.
    """
    # i.e. 'assay_results.tsv'
    assay_results_path = input_data_dir + assay_file
    assay_results = pd.read_csv(assay_results_path, sep='\t', index_col=0, header=0)
    assay_results.rename(columns={'Unnamed: 13': 'ttest_p'}, inplace=True)
    signif = assay_results[assay_results['ttest_p'] < 0.05]
    signif_metabs = signif.index.tolist()

    manually_curated_signif = open(input_data_dir + 'user_added_chebis_SIGNIF_ONLY.txt', 'r')
    man_cur_sig_fh = manually_curated_signif.readlines()

    # see chebi_id_metadata for explanations of the pre-specified signif chebi ids
    signif_chebi_ids = [ch_id.strip('\n') for ch_id in man_cur_sig_fh]
    signif_no_chebi_ids = []
    for name in signif_metabs:
        chebi_ID_obj = lc.search(name, exact=True)
        if len(chebi_ID_obj) > 0:
            signif_chebi_ids.append('CHEBI:' + str(chebi_ID_obj[0]._ChebiEntity__chebi_id))
        elif len(chebi_ID_obj) == 0:
            signif_no_chebi_ids.append(name)

    return signif_chebi_ids, signif_no_chebi_ids


def map_all_metabs_to_chebi_ids(metabs, signif_no_chebi_ids, output_data_dir):
    """
    Maps a list of metabolites to their ChEBI IDs, if available.
    See: https://github.com/libChEBI/libChEBIpy/blob/master/libchebipy/_chebi_entity.py


    Returns:
        exact_matches.tsv:
            2 col .tsv with provided metabolite name and matched CHEBI ID
                i.e.
                    malate	CHEBI:25115

        all_close_matches.tsv:
            2 col .tsv with provided metabolite name and dictionaries of
            possible matching CHEBI names and IDs for all metabolites without exact matches
                i.e.
                    GSH	[{'Gsh-prostaglandin A1': '5548'}, {'S-Decyl GSH': '8955'}]

        signif_close_matches.tsv: 2 col .tsv
            same format as all_close_matches.tsv but includes only significant metabolites
            with no exact CHEBI id matches

    """
    names_map = {}
    no_exact_match = {}
    signif_no_exact_match = {}
    for name in metabs:
        chebi_ID_obj = lc.search(name, exact=True)
        # works for 72 of the 112
        if len(chebi_ID_obj) > 0:
            names_map[name] = 'CHEBI:' + str(chebi_ID_obj[0]._ChebiEntity__chebi_id)
            # TODO: search used-to-produce sif for chebi_ID_obj[0]._ChebiEntity__chebi_id
        elif len(chebi_ID_obj) == 0:
            # if an exact match is not possible, don't use exact match
            chebi_ID_obj = lc.search(name)
            no_exact_match[name] = []
            for i in range(0,len(chebi_ID_obj)):
                # fill a dictionary with desired_name: {alternative_name, alternative_id}
                alt_id = str(chebi_ID_obj[i]._ChebiEntity__chebi_id)
                chebi_entity = lc.ChebiEntity(alt_id)
                no_exact_match[name].append({chebi_entity.get_name(): 'CHEBI:' + alt_id})
                if name in signif_no_chebi_ids:
                    signif_no_exact_match[name] = []
                    for i in range(0,len(chebi_ID_obj)):
                        alt_id = str(chebi_ID_obj[i]._ChebiEntity__chebi_id)
                        chebi_entity = lc.ChebiEntity(alt_id)
                        signif_no_exact_match[name].append(
                            {chebi_entity.get_name(): 'CHEBI:' + alt_id})

    exact_outf = open(output_data_dir + "exact_matches.tsv", "w")
    for k, v in names_map.items():
        exact_outf.write(str(k) + '\t' + str(v) + '\n')
    exact_outf.close()

    all_close_matches_outf = open(output_data_dir + "all_close_matches.tsv", "w")
    for k, v in no_exact_match.items():
        all_close_matches_outf.write(str(k) + '\t' + str(v) + '\n')
    all_close_matches_outf.close()

    signif_close_matches_path = output_data_dir + "signif_close_matches.tsv"
    signif_close_matches_outf = open(signif_close_matches_path, "w")
    for k, v in signif_no_exact_match.items():
        signif_close_matches_outf.write(str(k) + '\t' + str(v) + '\n')
    signif_close_matches_outf.close()

    return signif_close_matches_path


def convert_chebi_sif_to_named_sif(chebi_sif_path, chebis, name_to_id_map, input_dir):
    """
    Use compounds.tsv to convert ChEBI IDs in a .sif to biochemist-readable names.
    Writes a new file with same name as chebi sif, but replaces "chebi" with "named".
    In the future: make this bidirectional.
    """

    # load the sif that is going to be converted
    ch_sif_fh = open(chebi_sif_path, 'r')
    ch_sif = ch_sif_fh.readlines()

    name_sif_path = chebi_sif_path.replace('chebi', 'named')
    name_sif_fh = open(name_sif_path, 'w')

    # load compounds.tsv to search for names of unprofiled linkers that won't appear in name_to_id_map
    compounds_df = pd.read_csv(input_dir + 'compounds.tsv',
                               sep='\t',
                               usecols=[2, 5],
                               header=0,
                               index_col=0)

    # make a list of tuples (A,B) where A and B are names of metabolites
    named_node_pairs = []

    for relationship in ch_sif:
        # separate the .sif line into pieces
        # TODO: for some reason the normal approach isn't working. debug
        [chebi_a, edge, chebi_b] = get_parts_of_sif_line(relationship)

        # preset the names to chebi ids in case there is no name match
        name_a = chebi_a
        name_b = chebi_b

        # search for matches
        # if chebi_a represents or is derived from a profiled metabolite...
        if chebi_a in chebis:
            # find its parent name, if it has one
            for name, ch_list in name_to_id_map.items():
                if chebi_a in name_to_id_map[name]:
                    name_a = name

        # elif it's unprofiled
        elif chebi_a not in chebis:
            # it's probably a linker, try to find its name in compounds df
            if chebi_a in compounds_df.index:
                name_a = compounds_df.ix[chebi_a].tolist()[0]

        # repeat for chebi_b
        if chebi_b in chebis:
            for name, ch_list in name_to_id_map.items():
                if chebi_b in name_to_id_map[name]:
                    name_b = name
        elif chebi_b not in chebis:
            if chebi_b in compounds_df.index:
                name_b = compounds_df.ix[chebi_b].tolist()[0]

        named_node_pairs.append((name_a, name_b))
        # todo: look for duplicates in a second iteration

    # remove duplicates from named_node_pairs
    deduped_named_node_pairs = set(named_node_pairs)

    for pair in deduped_named_node_pairs:
        name_a = pair[0]
        name_b = pair[1]
        new_line = "{}\t{}\t{}\n".format(name_a, 'used-to-produce' , name_b)
        name_sif_fh.write(new_line)

    ch_sif_fh.close()
    name_sif_fh.close()

    return name_sif_path