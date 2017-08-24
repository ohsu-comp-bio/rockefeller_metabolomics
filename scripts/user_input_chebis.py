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

    # print("{}{}".format("\nPlease review ", signif_close_matches_file))
    print("You may specify additional ChEBI IDs to include in network generation\n"
          "by generating a file in input/\n\n"
          "Format: 1 CHEBI ID per line i.e.\n"
          "CHEBI:12345\n"
          "CHEBI:23456\n")
    print("Create 1 such file containing ALL additional ChEBIS and 1 file containing\n"
          "ONLY the significantly altered ChEBIs (for network formatting purposes)")
    # print("{}{}".format("You may also wish to review\n", misfit_loc))
    # print("Exact CHEBI IDs in that file do not appear in the used-to-produce.sif.\n"
    print("You may search https://www.ebi.ac.uk/chebi/ for alternatives to add to\n"
          "the files described above.")


    # ask user where they stored their manually curated chebis
    addl_chebis_file = input("Enter the name of your 'ALL' additional ChEBI IDs file\n"
                             "(or 'n' to skip): ")

    # TODO: CONTINUE
    signif_addl_chebis_file = input("Enter the name of your 'SIGNIFICANT' additional ChEBI IDs file\n"
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

    if signif_addl_chebis_file != 'n':
        while not os.path.exists(input_data_dir + signif_addl_chebis_file):
            print(signif_addl_chebis_file + " not found.\n "
                                     "Ensure that file name is correct and that it is located in\n" +
                                     input_data_dir)
            signif_addl_chebis_file = input("Enter the name of the file (or 'n' to skip): ")
            # if the user changes her mind and decides not to use a file, break loop
            if signif_addl_chebis_file == 'n':
                break

    all_addl = None
    signif_addl = None

    if addl_chebis_file != 'n':
        all_addl = input_data_dir + addl_chebis_file

    if signif_addl_chebis_file != 'n':
        signif_addl = input_data_dir + signif_addl_chebis_file

    return all_addl, signif_addl


def make_addl_chebis_files_from_metadata(metadata_dir):
    """
    This might be a bad idea but it's happening because I'm lazy.
    Extract additional chebis from metadata/meta_data_user_specified_chebis.txt
    and build user_added_chebis_ALL.txt
    """
    notes_loc = metadata_dir + 'meta_data_user_specified_chebis.txt'
    notes_fh = open(notes_loc, "r")
    notes = notes_fh.readlines()

    addl_chebis = []
    signif_addl_chebis = []
    for line in notes:
        if line.startswith('CHEBI:'):
            line = line.split('\t')
            chebi_id = line[0]
            signif = line[3]
            addl_chebis.append(chebi_id)
            if signif == 'YES':
                signif_addl_chebis.append(chebi_id)

    addl_chebis_file = open(input_data_dir + 'user_added_chebis_ALL.txt', 'w')
    signif_addl_chebis_file = open(input_data_dir + 'user_added_chebis_SIGNIF_ONLY.txt', 'w')

    for i in addl_chebis:
        addl_chebis_file.write(i + '\n')
    addl_chebis_file.close()

    for i in signif_addl_chebis:
        signif_addl_chebis_file.write(i + '\n')
    signif_addl_chebis_file.close()
