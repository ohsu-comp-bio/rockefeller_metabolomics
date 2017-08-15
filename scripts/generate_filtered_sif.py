"""
Functions for constructing a .sif network file containing desired
relationships between provided ChEBI IDs. Resulting network file
will contain either ChEBI IDs or corresponding human-readable names,
depending on user input.

Specified relationships may be "BOTH", "EITHER", or "DISTANCE-OF-2".
    BOTH: both nodes must be present in list of ChEBI IDs of interest.
    EITHER: presence of either node in list of ChEBI IDs of interest is
        sufficient for inclusion in the resulting .sif.
    DISTANCE-OF-2: builds a .sif as in "BOTH" but also includes pairs where,
        for example, nodeA and nodeC are of interest in:
            nodeA--->nodeB
            nodeB--->nodeC

Author: Hannah Manning <manningh@ohsu.edu>
Date: August 8, 2017
"""

import pandas as pd

def filter_sif_by_chebi(edge, path_to_addl_chebis):
    """
    Generates a list of chebi IDs from a file with the format of exact_matches.tsv
     col1: name provided by researchers
     col2: chebi ID
     (i.e. 'malate	CHEBI:25115')
    """
    sif_loc = output_data_dir + edge + ".sif"
    sif_fh = open(sif_loc, "r")
    exact_df = pd.read_csv(output_data_dir + "exact_matches.tsv", sep='\t')

    sif = sif_fh.readlines()
    chebis = exact_df.ix[:,1].tolist()
    if path_to_addl_chebis is not None:
        addl_chebis_file = open(path_to_addl_chebis, "r")
        addl_chebis_fh = addl_chebis_file.readlines()
        for chb in addl_chebis_fh:
            chb = chb.strip('\n')
            chebis.append(chb)
    both_filtered_sif = "chebi-filtered-" + edge + "-BOTH-MEMBERS.sif"
    either_filtered_sif = "chebi-filtered-" + edge + "-EITHER-MEMBER.sif"

    both_matched_file = open(output_data_dir + both_filtered_sif, "w")
    either_matched_file = open(output_data_dir + either_filtered_sif, "w")
    for relationship in sif:
        relationship = relationship.strip('\n')
        relationship_parts = relationship.split('\t')
        members = [relationship_parts[0], relationship_parts[2]]

        # to generate a sif where 1 member match is sufficient:
        if any(smallmol in chebis for smallmol in members):
            either_matched_file.write(relationship + '\n')

        # require that both members are present in the matched chebi IDs
        if all(smallmol in chebis for smallmol in members):
            both_matched_file.write(relationship + '\n')

    sif_fh.close()
    both_matched_file.close()
    either_matched_file.close()

    return output_data_dir + both_filtered_sif, output_data_dir + either_filtered_sif


def specify_chibe_formatting(sig_chebi_ids, sif_path):
    """
    Specify formatting for visualization with Chisio BioPAX Editor (ChiBE)
    For the sake of clarity: ChiBE != ChEBI
    Generates a .format file in the same directory as .sif used by ChiBE


    Notes from Ozgun:
    If ChiBE finds a .format file with the same name at the same directory of the SIF file,
    it will use that format file to decorate the SIF graph. Colors are given in RGB (red green blue).

    color: The background color of the node
    textcolor: Color of the node text (name)
    bordercolor: Border color
    borderwidth: Width of border, integer unfortunately, default is 1 and 2 is wide enough to stress.
    tooltip: Sets the text that will show up when the mouse is over the node.

    The .format file is tab delimited text file. Example:
    node CXCR4 color 23 65 13
    node RTN4 textcolor 0 0 0
    node IL6ST bordercolor 180 23 14
    node CARD9 borderwidth 2
    node IRS1 tooltip

    This is tooltip text
    """

    format_path = sif_path.replace('sif','format')
    format_fh = open(format_path, "w")

    sif_fh = open(sif_path, "r")
    sif = sif_fh.readlines()
    sif_chebis = []
    for relationship in sif:
        relationship = relationship.strip('\n')
        relationship_members = relationship.split('\t')
        if relationship_members[0] not in sif_chebis:
            sif_chebis.append(relationship_members[0])
        if relationship_members[2] not in sif_chebis:
            sif_chebis.append(relationship_members[2])

    present_sig_chebi_ids = []
    for id in sig_chebi_ids:
        if id in sif_chebis:
            present_sig_chebi_ids.append(id)

    for id in present_sig_chebi_ids:
        format_fh.write("node\t" + id + "\tcolor\t23 222 176\n")
    format_fh.close()


def convert_sif_chebis_to_names(chid_to_name_map, chebi_sif, new_sif_name):
    """
    Uses a 2 col tsv (all_chebis_of_interest.txt) to convert chebi ID's in a
    filtered sif back to their original names.

    chid_to_name_map = name of 2 col tsv where col 1 = "CHEBI:12345" and col 2 is name
    chebi_sif = filtered sif i.e. chebi-filtered-used-to-produce-BOTH-MEMBERS.sif
    new_sif_name = i.e. "named-filtered-used-to-produce-BOTH-MEMBERS.sif"

    Returns new sif file with name specified
    """
    chebi_sif_path = output_data_dir + chebi_sif
    chebi_sif_file = open(chebi_sif_path, "r")
    chebi_sif = chebi_sif_file.readlines()

    name_map = pd.read_table(output_data_dir + chid_to_name_map, sep = '\t', header=None, index_col=0)

    sif_outf = open(output_data_dir + new_sif_name, "w")

    for relationship in chebi_sif:
        relationship = relationship.strip('\n')
        parts = relationship.split('\t')
        chebi_a = parts[0]
        edge = parts[1]
        chebi_b = parts[2]

        # traverse the name map file
        if chebi_a in name_map.index:
            name_a = name_map.get_value(chebi_a, 1)
        elif chebi_a not in name_map.index:
            print("Failed ChEBI ID to name map: " + chebi_a)
            name_a = chebi_a
        if chebi_b in name_map.index:
            name_b = name_map.get_value(chebi_b, 1)
        elif chebi_b not in name_map.index:
            print("Failed ChEBI ID to name map: " + chebi_b)
            name_b = chebi_b
        sif_outf.write(name_a + '\t' + edge + '\t' + name_b + '\n')

    chebi_sif_file.close()
    sif_outf.close()
    # change format file to reflect IDs
    if os.path.isfile(chebi_sif_path.replace('.sif','.format')):
        format_filepath = chebi_sif_path.replace('.sif','.format')
        format_df = pd.read_csv(format_filepath, sep='\t', header=None)
        chebis = pd.Series(format_df.ix[:,1])

        for id in chebis:
            if id in name_map.index:
                name = name_map.get_value(id, 1)
                chebis.replace(id, name, inplace=True)
        format_df[1] = chebis
        format_df.to_csv(output_data_dir + new_sif_name.replace('.sif', '.format'),
                         sep='\t', header=False, index=False)