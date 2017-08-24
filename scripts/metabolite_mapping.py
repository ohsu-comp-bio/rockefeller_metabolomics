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


def get_metabolites(assay_file, input_data_dir):
    """
    Gets list of metabolites from assay_results.tsv.

    assay_results.tsv:
        -matrix with cell lines as cols and metabolites as row names
        -final column is the p-value resulting from a t-test between
            sensitive and resistant cell lines
    """

    # i.e. 'assay_results.tsv'
    assay_results_path = input_data_dir + assay_file
    assay_results = pd.read_csv(assay_results_path, sep='\t', index_col=0, header=0)
    metabolites = assay_results.index.tolist()

    return metabolites


def get_significant_metabs(assay_file, input_data_dir):
    """
    Returns a list of metabolites whose fold change in abundance achieved significance
    at the p < 0.05 level according to a t-test between sensitive and resistant strains.
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


def convert_chebi_sif_to_named_sif(chebi_sif_path, input_data_dir, output_data_dir):
    """
    Use compounds.tsv to convert ChEBI IDs in a .sif to biochemist-readable names.
    Writes a new file with same name as chebi sif, but replaces "chebi" with "named".
    In the future: make this bidirectional.
    """

    # load compounds.tsv
    compounds_df = pd.read_csv(input_data_dir + 'compounds.tsv',
                               sep = '\t',
                               header=0,
                               index_col=2)

    # generate a ChEBI-ID-to-name map in the form of a 1-col pd.DataFrame
    # index = ChEBI IDs
    id_name_map = compounds_df[['NAME']]

    # load the sif whose members are labeled with ChEBI IDs
    chebi_sif_fh = open(chebi_sif_path, "r")
    chebi_sif = chebi_sif_fh.readlines()

    name_sif_path = chebi_sif_path.replace('chebi','named')
    name_sif_fh = open(name_sif_path, "w")

    for line in chebi_sif:
        line = line.strip('\n')
        parts = line.split('\t')
        chebi_a = parts[0]
        chebi_b = parts[2]
        edge = parts[1]

        # set default "name" value as chebi ID in case there is no match
        name_a = chebi_a
        name_b = chebi_b
        if chebi_a in id_name_map.index:
            name_a = id_name_map.ix[chebi_a,0]
        elif chebi_a not in id_name_map.index:
            print("Failed to map " + chebi_a + " to name.")
        if chebi_b in id_name_map.index:
            name_b = id_name_map.ix[chebi_b,0]
        elif chebi_b not in id_name_map.index:
            print("Failed to map " + chebi_b + " to name.")

        new_line = "{}\t{}\t{}\n".format(name_a, edge, name_b)
        name_sif_fh.write(new_line)

    chebi_sif_fh.close()
    name_sif_fh.close()

    # change format file to reflect IDs

    # if there is a format file
    if os.path.isfile(chebi_sif_path.replace('.sif','.format')):
        # set the path to that format file
        format_filepath = chebi_sif_path.replace('.sif','.format')
        # load the format file's contents into a dataframe for ease of manipulation
        format_df = pd.read_csv(format_filepath, sep='\t', header=None)
        # collect the chebis column
        format_chebis = pd.Series(format_df.ix[:,1])

        for ch_id in format_chebis:
            # preset name to the ChEBI ID in case it doesn't map
            name = ch_id
            if ch_id in id_name_map.index:
                name = id_name_map.get_value(ch_id, 'NAME')
            format_chebis.replace(ch_id, name, inplace=True)
        format_df[1] = format_chebis
        format_df.to_csv(name_sif_path.replace('.sif', '.format'),
                         sep='\t', header=False, index=False)