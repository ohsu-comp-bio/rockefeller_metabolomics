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
import numpy as np
import operator
from utils import *
from pull_from_metadata import *

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

    # create a dict where keys = names and vals = list of derived chebi ids
    name_ids_dict = group_chebis_of_same_parent()

    # load pairwise spearman matrix
    pairwise_spearman_matrix = pd.read_csv(input_data_dir +
                                           'pairwise_spearman_corr_coef.tsv',
                                           sep='\t',
                                           index_col=0)

    # set minimum correlation for inclusion of edges
    min_corr = 0.5

    # reminder: relationship members are ChEBI IDs
    for relationship in sif:

        relationship = relationship.strip('\n')
        relationship_parts = relationship.split('\t')
        members = [relationship_parts[0], relationship_parts[2]]

        # drop the relationship if both members are profiled and
        # they have a low correlation
        drop = filter_by_pairwise_correlation(members,
                                              name_ids_dict,
                                              pairwise_spearman_matrix,
                                              min_corr)
        # if the pair has a decent correlation...
        if not drop:
            # build BOTH and EITHER networks in one iteration
            # to generate a sif where 1 member match is sufficient:
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
        [producer, edge, produced] = get_parts_of_sif_line(rel)

        # if the producer is in the list of profiled chebis...
        if producer in chebis:
            profiled_producer.append(producer)
            unprofiled_produced.append(produced)

        # if the produced is in teh list of profiled chebis...
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

    # remove duplicates
    # note that almost every unprofiled produced is also an unprofiled producer
    # means almost all are middlemen (that's ok)
    middlemen = list(set(middlemen))

    # initialize dictionaries to keep track of who produces whom in middlemen relationships
    # produces_mm = {middleman1: [list of metabs that produce it]}
    # produced_by_mm = {middleman1: [list of metabs it produces]}
    produces_mm = {key: [] for key in middlemen}
    produced_by_mm = {key: [] for key in middlemen}
    for rel in eithers:
        [producer, edge, produced] = get_parts_of_sif_line(rel)
        if producer in middlemen:
            produced_by_mm[producer].append(produced)
        if produced in middlemen:
            produces_mm[produced].append(producer)

    # make dictionary to describe Profiled_a --> Unprofiled_b --> Profiled_c
    # but only fill it with the ones with sufficient correlation
    pup = {}
    for middleman, its_producers in produces_mm.items():
        pup[middleman] = []
        for a in its_producers:
            for c in produced_by_mm[middleman]:
                drop = filter_by_pairwise_correlation([a,c],
                                                      name_ids_dict,
                                                      pairwise_spearman_matrix,
                                                      min_corr)
                if not drop:
                    pup[middleman].append((a,c))

    pre_dist2 = []

    # store all the relationships which meet profiledA --> unprofiledB --> profiledA
    # for later destruction
    cycles = [] # change type?

    # generate dictionaries of profiled producer: produced and produced: producer
    # cryptic...
    # dict1 k = profiled producer
    # dict1 v = unprofiled produced
    # dict2 k = unprofiled producer
    # dict2 v = profiled produced
    dict1 = {}
    dict2 = {}

    # only add those to pre_dist2 if the unprofiled entity is also a middleman
    for rel in eithers:
        [producer, edge, produced] = get_parts_of_sif_line(rel)

        # if producer is profiled and produced is a middleman, keep it
        if producer in profiled_producer:
            if produced in middlemen:
                # make all profiled producers keys in dict1
                # append all their unprofiled produced entities as vals
                if producer in dict1:
                    dict1[producer].append(produced)
                if producer not in dict1:
                    dict1[producer] = [produced]
                # keeping the whole line
                pre_dist2.append(rel)

        # if produced is profiled and producer is a middleman, keep it
        if produced in profiled_produced:
            if producer in middlemen:
                # make all unprofiled producers keys in dict2
                # append all their profiled produced entities as vals
                if producer in dict2:
                    dict2[producer].append(produced)
                if producer not in dict2:
                    dict2[producer] = [produced]
                # keeping the whole line
                pre_dist2.append(rel)

    # for profiled producer, unprofiled produced list in dict1
    for k, vals in dict1.items():
        # iterate over the unprofiled produced entities
        for v in vals:
            # if unprofiled produced is an unprofiled producer
            if v in dict2:
                # and if it is an unprofiled producer of our dict1 profiled producer
                if k in dict2[v]:
                    # then it is a cyclic A --> B --> A relationship and it needs to go
                    cycles.append((k,v))

    acyclic_predist2 = []
    for rel in pre_dist2:
        [producer, edge, produced] = get_parts_of_sif_line(rel)
        acyclic = True
        for cyclic_rel in cycles:
            # if the producer --> produced relationship is in the list of cyclic ones...
            if all(entity in [producer, produced] for entity in cyclic_rel):
                # flag it as cyclic
                acyclic = False
        # if the relationship is acyclic (at distance of 2), then keep it
        if acyclic:
            acyclic_predist2.append(rel)

    # determine which edges to keep in distance of 2 sif based on correlation of A and C
    # in example profiledA --> unprofiledB --> profiledC
    # pup[B] = [(A,C), ... ]
    print("{}{}".format("the length of acyclic_predist was: ", str(len(acyclic_predist2))))
    acyc_corr_predist2 = []
    for rel in acyclic_predist2:
        [producer, edge, produced] = get_parts_of_sif_line(rel)
        # start by assuming the edge will be dropped
        drop = True

        # if the producer is a middleman (unprofiled),
        if producer in pup:
            B = producer    # for the sake of clarity
            for pair in pup[B]:
                # if the second member of any of these tuples (A, C) matches produced
                # that means the correlation between those two profiled compounds (A and C)
                # was deemed sufficient to maintain an edge between B (producer) and C (produced)
                if pair[1] == produced:
                    drop = False

        # if the produced entity is a middleman (unprofiled),
        if produced in pup:
            B = produced
            for pair in pup[B]:
                # if the first member of any of these tuples (A, C) matches producer
                # ''
                # was deemed sufficient to maintain an edge between A (producer) and B (produced)
                if pair[0] == producer:
                    drop = False

        if not drop:
            acyc_corr_predist2.append(rel)
    print("{}{}".format("the length of acyc_corr_predist2 is now: ", str(len(acyc_corr_predist2))))

    # go through the pre_dist relationships and only keep the ones whose middlemen
    # are "non-ubiquitous" (i.e. appears less than 2 times as a producer and as produced)
    # return a dataframe containing unprofiled producers
    # and how many edges they yield in pre_dist2
    [sorted_producers, sorted_produced] = sort_sif_list_by_biggest_producer(acyc_corr_predist2, middlemen)

    # view histogram of these values
    # the vast majority of counts are below 6.
    # f = plt.hist(testing456[testing456[1]>5], 50, facecolor = 'blue', alpha = 0.6)
    # plt.close()

    # print("original shape of sorted_producers: " + str(shape(sorted_producers)))
    # print("original shape of sorted_produced: " + str(shape(sorted_produced)))

    # i've decided to exclude sif relationships in which the middleman
    # takes part as either a producer or produced entity in more than 4 relationships
    sorted_producers = sorted_producers[sorted_producers[1] < 4]
    sorted_produced = sorted_produced[sorted_produced[1] < 4]

    # print("shape after restricting sorted_producers: " + str(shape(sorted_producers)))
    # print("shape after restricting sorted_produced: " + str(shape(sorted_produced)))

    ok_unprofiled_producers = sorted_producers.index
    ok_unprofiled_produced = sorted_produced.index

    # get the union between the ok lists, use only those
    ok_unprofiled = ok_unprofiled_producers & ok_unprofiled_produced

    # only include relationships from either where the unprofiled entity isn't too promiscuous
    for rel in acyc_corr_predist2:
        [producer, edge, produced] = get_parts_of_sif_line(rel)

        if producer in chebis:
            if produced in ok_unprofiled:
                dist_2_matched_file.write(rel + '\n')
        if produced in chebis:
            if producer in ok_unprofiled:
                dist_2_matched_file.write(rel + '\n')

    dist_2_matched_file.close()

    return output_data_dir + both_filtered_sif, \
           output_data_dir + either_filtered_sif, \
           output_data_dir + dist_2_filtered_sif


def sort_sif_list_by_biggest_producer(sif_list, middlemen):
    """
    Take a list of node-edge-node relationships (as in either_sif_fh.readlines())
    and a list of "middlemen" in this .sif.

    middlemen are unprofiled entities which appear as both producers and produced
    entities in the sif.

    Returns a sorted pd.DataFrame of unprofiled producers and the
    number of node-edge-node relationships for which they are the first node.
    index = ChEBI IDs, single col = number f edges for which it is a producer
    """
    producer_counts = {}
    produced_counts = {}
    for line in sif_list:
        line = line.strip('\n')
        line = line.split('\t')
        producer = line[0]
        produced = line[2]
        if producer in middlemen:
            if producer not in producer_counts:
                producer_counts[producer] = 1
            if producer in producer_counts:
                producer_counts[producer] += 1
        if produced in middlemen:
            if produced not in produced_counts:
                produced_counts[produced] = 1
            if produced in produced_counts:
                produced_counts[produced] += 1

    # generate a 1col dataframe with index = ChEBI IDs and
    # col1 = number of edges for which it's a producer
    producer_counts_sorted = pd.DataFrame(sorted(producer_counts.items(),
                                                 key=operator.itemgetter(1),
                                                 reverse=True))
    producer_counts_sorted = producer_counts_sorted.set_index(0)

    # do same for produced
    produced_counts_sorted = pd.DataFrame(sorted(produced_counts.items(),
                                                 key=operator.itemgetter(1),
                                                 reverse=True))
    produced_counts_sorted = produced_counts_sorted.set_index(0)

    return producer_counts_sorted, produced_counts_sorted


# TODO: filter sif by correlation between the two metabs
def filter_by_pairwise_correlation(pair,
                                   name_ids_dict,
                                   pairwise_spearman_matrix,
                                   minimum_corr):
    """
    Assesses whether an edge connecting the given pair has sufficient correlation
    to justify its inclusion in the downstream .sif.

    Draws correlations from a matrix of spearman correlations
    (i.e. data/input/pairwise_spearman_corr_coef.tsv)
    (which was generated by assay_results_stats.py).

    Important example:
        relationship1:  profiledA --> unprofiledB --> profiledC
        relationship2:  profiledD --> unprofiledB --> profiledC

         Let's say that A and C are significantly correlated, but D and C are not.
         In this case, the B node will be maintained in the graph.
         The following edges will be maintained:
            profiledA --> unprofiledB
            unprofiledB --> profiledC
         And the following will be excluded:
            profiledD --> unprofiledB

    pair = [ChEBI_ID_1, ChEBI_ID_2]
    name_ids_dict is provided by group_chebis_of_same_parent()
    minimum_corr (float) is a value between 0 and 1 that designates the minimum correlation
        (or anti-correlation) for the pair to be considered sufficient

    Returns True if correlation is sufficient to include the relationship in the .sif.
    """
    # make appropriate null values
    drop = False
    parent1 = None
    parent2 = None

    chebi1 = pair[0]
    chebi2 = pair[1]

    for name, chebi_ids in name_ids_dict.items():
        if chebi1 in name_ids_dict[name]:
            parent1 = name
        if chebi2 in name_ids_dict[name]:
            parent2 = name

    # if both parent names are assigned, we know they're both profiled
    # we therefore have correlation data to assess whether the edge should be dropped
    # if the parent name was not assigned, then this ID has the potential to be an
    # unprofiled linker node and the edge should be maintained (drop = False) for future use
    # distance of 2 filtering
    if all([parent1, parent2]):
        # ensure that you're not trying to assess the correlation of 2 chebi ids derived
        # from the same parent
        if parent1 != parent2:

            # half of the matrix is NaNs to avoid redundancy
            # if your indexing brings you to this half of the matrix, swap your row & col vals
            if np.isnan(pairwise_spearman_matrix.ix[parent1, parent2]):
                if abs(pairwise_spearman_matrix.ix[parent2, parent1]) < minimum_corr:
                    drop = True

            elif abs(pairwise_spearman_matrix.ix[parent1, parent2]) < minimum_corr:
                drop = True

    return drop