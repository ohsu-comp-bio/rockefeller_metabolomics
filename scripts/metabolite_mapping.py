import libchebipy as lc
import pandas as pd


def get_metabolites(assay_file):
    """
    Gets list of metabolites from assay_results.tsv.
    """

    # i.e. 'assay_results.tsv'
    assay_results_path = input_data_dir + assay_file
    assay_results = pd.read_csv(assay_results_path, sep='\t', index_col=0, header=0)
    metabolites = assay_results.index.tolist()

    return metabolites


def get_significant_metabs(assay_file):
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
    signif_chebi_ids = [id.strip('\n') for id in man_cur_sig_fh]
    signif_no_chebi_ids = []
    for name in signif_metabs:
        chebi_ID_obj = lc.search(name, exact=True)
        if len(chebi_ID_obj) > 0:
            signif_chebi_ids.append('CHEBI:' + str(chebi_ID_obj[0]._ChebiEntity__chebi_id))
        elif len(chebi_ID_obj) == 0:
            signif_no_chebi_ids.append(name)

    return signif_chebi_ids, signif_no_chebi_ids


def map_all_metabs_to_chebi_ids(metabs, signif_no_chebi_ids):
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
                        signif_no_exact_match[name].append({chebi_entity.get_name(): 'CHEBI:' + alt_id})

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