"""
Function(s) for filtering a .sif network file by a given edge
(i.e. used-to-produce).

Author: Hannah Manning <manningh@ohsu.edu>
Date: August 8, 2017
"""

import shutil
import gzip


def select_sif_edges(edge, gzipped_sif, input_data_dir, output_data_dir):
    """
    Filters full pathway commons network file for only the relationships
    characterized by an edge of interest (i.e. 'used-to-produce')
    """
    # i.e. 'PathwayCommons9.All.hgnc.sif.gz'
    sif_path = input_data_dir + gzipped_sif

    # i.e. 'used-to-produce'
    desired_edge = edge
    outf_name = output_data_dir + desired_edge + '.sif'

    # open gzipped sif file
    with gzip.open(sif_path, 'rt') as inf:
        sif = inf.readlines()

    outfh = open(outf_name, "w")

    for line in sif:
        if line == '\n':
            break
        line = line.strip('\n')
        current_edge = line.split('\t')[1]
        if current_edge == desired_edge:
            outfh.write(line + '\n')

    outfh.close()

    # compress your new file
    with open(outf_name, 'rb') as not_zipped, gzip.open(outf_name + '.gz', 'wb') as zipped:
        shutil.copyfileobj(not_zipped, zipped)