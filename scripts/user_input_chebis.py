"""
Function for integrating user input in the generation of .sifs.
"""

import os

def ask_user_for_chebis(misfit_loc, signif_no_match_loc, input_data_dir):
    """
    Provides location of file with suggested matches for non-matching, significant chebis.
    Asks user to select which to keep and store them in a separate file.
    Verifies the existence of this file and returns path to user-specified file or None.
    """
    signif_close_matches_file = signif_no_match_loc.replace('/scripts/..','')
    misfit_loc = misfit_loc.replace('/scripts/..','')

    print("{}{}".format("\nPlease review ", signif_close_matches_file))
    print("Select CHEBI IDs from that file that you would like\n"
          "to include in the network by generating a file in input/\n"
          "Format: 1 CHEBI ID per line i.e.\n"
          "CHEBI:12345\n"
          "CHEBI:23456\n")
    print("{}{}".format("You may also wish to review\n", misfit_loc))
    print("Exact CHEBI IDs in that file do not appear in the used-to-produce.sif.\n"
          "You may search https://www.ebi.ac.uk/chebi/ for alternatives to add to\n"
          "the file described above.")


    # ask user where they stored their manually curated chebis
    addl_chebis_file = input("Enter the name of your additional CHEBI IDs file\n"
                             "(or 'n' to skip): ")
    if addl_chebis_file != 'n':
        while not os.path.exists(input_data_dir + addl_chebis_file):
            print(addl_chebis_file + " not found.\n "
                                     "Ensure that file name is correct and that it is located in\n" +
                                     input_data_dir)
            addl_chebis_file = input("Enter the name of the file (or 'n' to skip): ")
            # if the user changes her mind and decides not to use a file, break loop
            if addl_chebis_file == 'n':
                break

    if addl_chebis_file == 'n':
        return None

    else:
        return input_data_dir + addl_chebis_file


def make_addl_chebis_file_from_metadata(metadata_dir):
    """
    This is a bad idea but it's happening because I'm lazy.
    Extract additional chebis from metadata/meta_data_user_specified_chebis.txt
    and build user_added_chebis_ALL.txt
    """
    notes_loc = metadata_dir + 'meta_data_user_specified_chebis.txt'
    notes_fh = open(notes_loc, "r")
    notes = notes_fh.readlines()

    addl_chebis = []
    for line in notes:
        if line.startswith('CHEBI:'):
            line = line.split('\t')
            addl_chebis.append(line[0])

    addl_chebis_file = open(input_data_dir + 'user_added_chebis_ALL.txt', 'w')
    for i in addl_chebis:
        addl_chebis_file.write(i + '\n')
    addl_chebis_file.close()


def collect_all_chebis_being_used(metadata_dir):
    """
    Another slapped-together function that will serve a purpose right now
    but maybe not later.
    Gets CHEBI IDs from exact_matches.tsv and from user_added_chebis_ALL.txt
    Writes them out to all_chebis_of_interest.txt.

    AND MAPS THEM TO THEIR NON-CHEBI-ID NAMES.
    """
    user_added_path = input_data_dir + 'user_added_chebis_ALL.txt'
    user_added_file = open(user_added_path, "r")
    user_added = user_added_file.readlines()

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