"""
Functions to generate .format files for .sif networks to be visualized in ChiBE.

A useful tool for picking complementary colors: https://www.sessions.edu/color-calculator/

Author: Hannah Manning
Date: August 15, 2017
"""

from pull_from_metadata import *
import numpy as np
import sys

base_dir = os.path.dirname(os.path.realpath(__file__))
input_data_dir = base_dir + '/../data/input/'
metadata_dir = base_dir + '/../metadata/'
output_data_dir = base_dir + '/../data/output/'
sys.path += [base_dir + '/../../Rockefeller_Metabolomics']

def specify_chibe_formatting(sig_chebi_ids, sif_path, assay_file):
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

    """

    # TODO: instead of feeding it signif_chebi_ids, just use assay_results_df

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

    # color the borders of significantly altered metabs red
    for ch_id in present_sig_chebi_ids:
        format_fh.write("node\t" + ch_id + "\tbordercolor\t255 0 0\n")
        format_fh.write("node\t" + ch_id + "\tborderwidth\t3\n")

    # read in the extended results to get spearman correlation coefficients
    assay_results = pd.read_csv(input_data_dir + assay_file,
                                sep='\t',
                                index_col=0,
                                header=0,
                                na_values = 'nd')

    # get the name to all-IDs dictionary
    name_chebis_map = group_chebis_of_same_parent()

    # TODO: Change to using fold change instead of corr
    for name, chebi_list in name_chebis_map.items():
        fc = assay_results['fc_gmeans'][name]
        rgb = calculate_node_color(fc)
        rgb_str = "{} {} {}\n".format(rgb[0], rgb[1], rgb[2])
        for ch_id in chebi_list:
            format_fh.write("node\t" + ch_id + "\tcolor\t" + rgb_str)

    format_fh.close()

    return format_path

def calculate_node_color(fold_change):
    """
    Calculates node color based on fold change.
    Returns node color in RGB format (i.e. np.array([100, 150, 100]))
    Assumes color is RGB between [0, 0, 0] and [255, 255, 255]
    """

    # TODO: Deal with fold changes not being out of 1 or -1...
    gray = [220, 220, 220]
    red = [255, 89, 0]
    darkblue = [0, 111, 255]

    white = np.array([255, 255, 255])
    blue = np.array([71, 151, 255])
    orange = np.array([255, 150, 94])


    blue_vector = white - blue
    orange_vector = white - orange

    if fold_change < -1.0:
        rgb = darkblue

    # if the spearman correlation is negative, set the node to a gradation of orange
    if -1.0 <= fold_change < -0.1:
        rgb = blue + blue_vector * (1 + fold_change)
        rgb = [int(i) for i in rgb.tolist()]

    if -0.1 <= fold_change <= 0.1:
        rgb = gray

    # if the spearman correlation is positive, set the node to a gradation of blue
    if 0.1 < fold_change <= 1.0:
        rgb = orange + orange_vector * (1 - fold_change)
        rgb = [int(i) for i in rgb.tolist()]

    if fold_change > 1.0:
        rgb = red

    return rgb
