"""
Functions for
    (1) verifying that exact matches (provided name to a CheBI ID) map back
        to the unfiltered, single edge-type .sif.
    (2) extracting all members of a .sif

Author: Hannah Manning <manningh@ohsu.edu>
Date: August 8, 2017
"""

def verify_exact_chebis(sif_file, output_data_dir):
    """
    Checks whether exact CHEBI ID matches appear in the used-to-produce sif.
    If they do not, they are written to exact_matches_not_in_sif.tsv
    :return:
    """
    # exact_df = pd.read_csv(output_data_dir + "exact_matches.tsv", sep='\t')
    # chebis = exact_df.ix[:, 1].tolist()

    exact_path = output_data_dir + "exact_matches.tsv"
    exact_fh = open(exact_path, "r")
    exact = exact_fh.readlines()

    sif_chebi_list = make_list_of_sif_chebis(sif_file, output_data_dir)

    misfits = {}
    for exact_match in exact:
        exact_match = exact_match.strip('\n')
        entities = exact_match.split('\t')
        given_name = entities[0]
        matched_chebi = entities[1]
        if matched_chebi not in sif_chebi_list:
            misfits[given_name] = matched_chebi

    misfit_loc = output_data_dir + 'exact_matches_not_in_sif.tsv'
    misfit_outf = open(misfit_loc, 'w')
    for k, v in misfits.items():
        misfit_outf.write(str(k) + '\t' + str(v) + '\n')
    misfit_outf.close()

    return misfit_loc

def make_list_of_sif_chebis(sif_file, output_data_dir):
    """
    Parses a sif to generate a list of all ChEBI ids that it contains.
    """
    sif_path = output_data_dir + sif_file
    sif_fh = open(sif_path, "r")
    sif = sif_fh.readlines()

    chebis = []
    for relationship in sif:
        relationship = relationship.strip('\n')
        parts = relationship.split('\t')
        chebis.append(parts[0])
        chebis.append(parts[2])

    return chebis