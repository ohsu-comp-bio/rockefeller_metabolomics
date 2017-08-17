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

def filter_sif_by_chebi(edge, path_to_addl_chebis, output_data_dir):
    """
    Uses ChEBI IDs (from exact-matches and from a user-added file containing
    additional ChEBI IDs of interest) to filter a .sif with a single edge type
    (i.e. used-to-produce).

    Builds 3 .sif files:
        both_filtered_sif: contains relationships wherein both ChEBI IDs were of interest
        either_filtered_sif: contains relationships wherein either ChEBI ID was of interest
        dist_2_filtered_sif: contains relationships as in "both_filtered_sif" AND relationships
            wherein A and C are of interest but some intermediate B is not.
            i.e., both of the following would be kept, despite B not being of initial interest:
                    A used-to-produce B
                    B used-to-produce C

    """

    # load single-edged .sif network file (i.e. used-to-produce.sif)
    sif_loc = output_data_dir + edge + ".sif"
    sif_fh = open(sif_loc, "r")
    sif = sif_fh.readlines()

    # load the ChEBI IDs which matched exactly to provided names
    exact_df = pd.read_csv(output_data_dir + "exact_matches.tsv", sep='\t')

    # make a list of all the ChEBI IDs in exact and append all the user-added
    # ChEBI IDs to it
    chebis = exact_df.ix[:,1].tolist()
    if path_to_addl_chebis is not None:
        addl_chebis_file = open(path_to_addl_chebis, "r")
        addl_chebis_fh = addl_chebis_file.readlines()
        for chb in addl_chebis_fh:
            chb = chb.strip('\n')
            chebis.append(chb)

    # generate names of outfiles
    both_filtered_sif = "chebi-filtered-" + edge + "-BOTH-MEMBERS.sif"
    either_filtered_sif = "chebi-filtered-" + edge + "-EITHER-MEMBER.sif"
    dist_2_filtered_sif = "chebi-filtered-" + edge + "-DIST-OF-2.sif"

    # open outfiles for writing
    both_matched_file = open(output_data_dir + both_filtered_sif, "w")
    either_matched_file = open(output_data_dir + either_filtered_sif, "w")
    dist_2_matched_file = open(output_data_dir + dist_2_filtered_sif, "w")

    # prepare to collect all "either" relationships for distance-of-2 assessment
    eithers = []

    # build BOTH and EITHER networks in one iteration
    for relationship in sif:
        relationship = relationship.strip('\n')
        relationship_parts = relationship.split('\t')
        members = [relationship_parts[0], relationship_parts[2]]

        # to generate a sif where 1 member match is sufficient:
        # TODO: store these for use in determining distance-of-2 relationships
        if any(smallmol in chebis for smallmol in members):
            either_matched_file.write(relationship + '\n')
            eithers.append(relationship)

        # require that both members are present in the matched chebi IDs
        # write it also to the distance-of-2 network file
        if all(smallmol in chebis for smallmol in members):
            both_matched_file.write(relationship + '\n')
            dist_2_matched_file.write(relationship + '\n')

    # close single-edge sif, both, and either files
    sif_fh.close()
    both_matched_file.close()
    either_matched_file.close()

    # build remainder of distance-of-2 network in a second iteration.
    # prepare to store relationships members based on whether they
    # are of interest ("profiled")

    profiled_producer = []
    unprofiled_producer = []

    profiled_produced = []
    unprofiled_produced = []

    for rel in eithers:
        parts = rel.split('\t')
        producer = parts[0]
        edge = parts[1]
        produced = parts[2]
        if producer in chebis:
            profiled_producer.append(producer)
            unprofiled_produced.append(produced)
        if produced in chebis:
            profiled_produced.append(produced)
            unprofiled_producer.append(producer)

    # prepare to capture unprofiled entities which are
    # produced by a profiled entity. Then evaluate whether
    # they also produce a profiled entity.
    middlemen = []
    for u_pd in unprofiled_produced:
        if u_pd in unprofiled_producer:
            middlemen.append(u_pd)

    # iterate through eithers a final time
    for rel in eithers:
        parts = rel.split('\t')
        producer = parts[0]
        edge = parts[1]
        produced = parts[2]

        if producer in profiled_producer:
            if produced in middlemen:
                dist_2_matched_file.write(rel + '\n')
        if produced in profiled_produced:
            if producer in middlemen:
                dist_2_matched_file.write(rel + '\n')

    dist_2_matched_file.close()

    return output_data_dir + both_filtered_sif, \
           output_data_dir + either_filtered_sif, \
           output_data_dir + dist_2_filtered_sif



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
    for ch_id in sig_chebi_ids:
        if ch_id in sif_chebis:
            present_sig_chebi_ids.append(ch_id)

    for ch_id in present_sig_chebi_ids:
        format_fh.write("node\t" + ch_id + "\tcolor\t23 222 176\n")
    format_fh.close()

    return format_path
