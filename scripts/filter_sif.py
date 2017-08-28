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
import os


def filter_sif_by_chebi(chebis, name_ids_dict, min_corr,
                        linker_lenience, input_dir, output_dir):
    """
    Takes a list of all profiled or derived chebis.

    Repeats the following process 3 times:
        once for sensitive cell lines only
        once for resistant cell lines only
        once for all cell lines together

            Filters input/used-to-produce.sif in order to produce 3 filtered networks:
                chebi-BOTH.sif
                    both nodes are profiled
                        profiled: exactly match or be derived from a metabolite in the
                        original dataset
                chebi-EITHER.sif
                    at least one node in the relationship must be profiled
                chebi-DIST2.sif
                    either both nodes are profiled, OR
                    the unprofiled node in an "either" relationship serves as a linker
                    to another profiled node

    Note: these networks contain ChEBI IDs and not human-readable names.

    With the exception of non-both "EITHER.sif" relationships, the correlation between
    the profiled nodes must exceed that specified in min_corr in order to be present
        i.e. 0.7 indicates that A and C in A-->C or in A-->B-->C must be at least
        70% correlated or anticorrelated if A and C are profiled and B is not.

    Returns: output_locations_all, output_locations_sens, output_locations_res
        where each list of output_locations contains the paths to the
            both, either, and distance of 2 sif networks, respectively
    """
    # read in used-to-produce.sif to be filtered down
    sif_fh = open(input_dir + 'used-to-produce.sif', 'r')
    sif = sif_fh.readlines()

    # generate names of outfiles
    all_both_sif = "all_chebi_BOTH_members.sif"
    all_either_sif = "all_chebi_EITHER_member.sif"
    all_dist2_sif = "all_chebi_DIST2.sif"

    sens_both_sif = "sens_chebi_BOTH_members.sif"
    sens_either_sif = "sens_chebi_EITHER_member.sif"
    sens_dist2_sif = "sens_chebi_DIST2.sif"

    res_both_sif = "res_chebi_BOTH_members.sif"
    res_either_sif = "res_chebi_EITHER_member.sif"
    res_dist2_sif = "res_chebi_DIST2.sif"

    # open outfiles for writing
    all_both_file = open(output_dir + all_both_sif, "w")
    all_either_file = open(output_dir + all_either_sif, "w")
    all_dist2_file = open(output_dir + all_dist2_sif, "w")

    sens_both_file = open(output_dir + sens_both_sif, "w")
    sens_either_file = open(output_dir + sens_either_sif, "w")
    sens_dist2_file = open(output_dir + sens_dist2_sif, "w")

    res_both_file = open(output_dir + res_both_sif, "w")
    res_either_file = open(output_dir + res_either_sif, "w")
    res_dist2_file = open(output_dir + res_dist2_sif, "w")

    # prepare to collect .sif relationships in which either member is profiled
    all_eithers = []
    sens_eithers = []
    res_eithers = []

    # read in the metabolite-by-metabolite spearman correlation matrices
    all_corr_matrix = pd.read_csv(input_dir + 'all_pairwise_spearman_corr_coef.tsv',
                                           sep='\t',
                                           index_col=0)

    sens_corr_matrix = pd.read_csv(input_dir + 'sensitive_pairwise_spearman_corr_coef.tsv',
                                            sep='\t',
                                            index_col=0)

    res_corr_matrix = pd.read_csv(input_dir + 'resistant_pairwise_spearman_corr_coef.tsv',
                                            sep='\t',
                                            index_col=0)

    # begin filtration
    for relationship in sif:
        [producer, edge, produced] = get_parts_of_sif_line(relationship)

        # drop a relationship if both members are profiled and the correlation is low
        all_drop = filter_by_pairwise_corr([producer, produced],
                                           name_ids_dict,
                                           all_corr_matrix,
                                           min_corr)

        sens_drop = filter_by_pairwise_corr([producer, produced],
                                            name_ids_dict,
                                            all_corr_matrix,
                                            min_corr)

        res_drop = filter_by_pairwise_corr([producer, produced],
                                           name_ids_dict,
                                           all_corr_matrix,
                                           min_corr)

        # if the overall pair has a decent correlation
        # build BOTH and EITHER networks in one iteration
        if not all_drop:

            # generate a sif where 1 member match is sufficient:
            if any(smallmol in chebis for smallmol in [producer, produced]):
                all_either_file.write(relationship)
                all_eithers.append(relationship)

            # require that both members are present in the matched chebi IDs
            # write it also to the distance-of-2 network file
            if all(smallmol in chebis for smallmol in [producer, produced]):
                all_both_file.write(relationship)
                all_dist2_file.write(relationship)

        # repeat for sensitive
        if not sens_drop:
            # generate a sif where 1 member match is sufficient:
            if any(smallmol in chebis for smallmol in [producer, produced]):
                sens_either_file.write(relationship)
                sens_eithers.append(relationship)

            # require that both members are present in the matched chebi IDs
            # write it also to the distance-of-2 network file
            if all(smallmol in chebis for smallmol in [producer, produced]):
                sens_both_file.write(relationship)
                sens_dist2_file.write(relationship)

        # repeat for resistant
        if not res_drop:
            # generate a sif where 1 member match is sufficient:
            if any(smallmol in chebis for smallmol in [producer, produced]):
                res_either_file.write(relationship)
                res_eithers.append(relationship)

            # require that both members are present in the matched chebi IDs
            # write it also to the distance-of-2 network file
            if all(smallmol in chebis for smallmol in [producer, produced]):
                res_both_file.write(relationship)
                res_dist2_file.write(relationship)

    # The BOTH and EITHER sifs are now complete (with ChEBI IDs, not names yet)

    # close used-to-produce.sif
    sif_fh.close()

    # close the "BOTH.sif" outfiles
    all_both_file.close()
    sens_both_file.close()
    res_both_file.close()

    # close the "EITHER.sif" outfiles
    all_either_file.close()
    sens_either_file.close()
    res_either_file.close()

    # Now build the remainder of the distance-of-2 networks through further iteration

    # store relationship members based on whether their role and profiled status
    (all_profiled_producer,
     all_profiled_produced,
     all_unprofiled_producer,
     all_unprofiled_produced) = \
        categorize_eithers_metabs(all_eithers, chebis)

    (sens_profiled_producer,
     sens_profiled_produced,
     sens_unprofiled_producer,
     sens_unprofiled_produced) = \
        categorize_eithers_metabs(sens_eithers, chebis)

    (res_profiled_producer,
     res_profiled_produced,
     res_unprofiled_producer,
     res_unprofiled_produced) = \
        categorize_eithers_metabs(res_eithers, chebis)

    # identify unprofiled "linkers":
    # capture unprofiled entities which are produced by a profiled entity. Then evaluate
    # whether they also produce a profiled entity.
    all_linkers = find_linkers(all_unprofiled_producer, all_unprofiled_produced)
    sens_linkers = find_linkers(sens_unprofiled_producer, sens_unprofiled_produced)
    res_linkers = find_linkers(res_unprofiled_producer, res_unprofiled_produced)

    # keep track of who produces whom in linker-involved relationships
    (all_produces_linker, all_produced_by_linker) = who_produces_whom(all_eithers, all_linkers)
    (sens_produces_linker, sens_produced_by_linker) = who_produces_whom(sens_eithers, sens_linkers)
    (res_produces_linker, res_produced_by_linker) = who_produces_whom(res_eithers, res_linkers)

    # Make a dictionary where keys are linkers and values are lists of tuples
    # so the relationship ProfiledA --> UnprofiledB --> ProfiledC is represented: pup[B] = [(A,C)]
    all_pup = make_profiled_unprofiled_profiled_dict(all_produces_linker,
                                                     all_produced_by_linker,
                                                     name_ids_dict,
                                                     all_corr_matrix,
                                                     min_corr)

    sens_pup = make_profiled_unprofiled_profiled_dict(sens_produces_linker,
                                                      sens_produced_by_linker,
                                                      name_ids_dict,
                                                      all_corr_matrix,
                                                      min_corr)

    res_pup = make_profiled_unprofiled_profiled_dict(res_produces_linker,
                                                     res_produced_by_linker,
                                                     name_ids_dict,
                                                     all_corr_matrix,
                                                     min_corr)

    # TODO: make less cryptic
    # prepare a first pass for collecting distance-of-2 relationships
    # add relationships to pre_dist2 if not a both relationship and if unprofiled entity is a
    # correlation-approved linker
    (all_pre_dist2, all_dict1, all_dict2) = build_pre_distance_of_2(all_eithers,
                                                                    all_linkers,
                                                                    all_profiled_producer,
                                                                    all_profiled_produced)

    (sens_pre_dist2, sens_dict1, sens_dict2) = build_pre_distance_of_2(sens_eithers,
                                                                       sens_linkers,
                                                                       sens_profiled_producer,
                                                                       sens_profiled_produced)

    (res_pre_dist2, res_dict1, res_dict2) = build_pre_distance_of_2(res_eithers,
                                                                    res_linkers,
                                                                    res_profiled_producer,
                                                                    res_profiled_produced)

    # identify cycles wherein one profiled entity is linked to itself by an unprofiled entity
    all_cycles = find_unprofiled_cycles(all_dict1, all_dict2)
    sens_cycles = find_unprofiled_cycles(sens_dict1, sens_dict2)
    res_cycles = find_unprofiled_cycles(res_dict1, res_dict2)

    # remove cycles from pre_dist2 lists
    all_acyclic_pre_dist2 = remove_cycles(all_pre_dist2, all_cycles)
    sens_acyclic_pre_dist2 = remove_cycles(sens_pre_dist2, sens_cycles)
    res_acyclic_pre_dist2 = remove_cycles(res_pre_dist2, res_cycles)

    # a final correlation check is required for distance of two relationships
    all_acyclic_corr_pre_dist2 = pre_dist2_corr_filter(all_acyclic_pre_dist2, all_pup)
    sens_acyclic_corr_pre_dist2 = pre_dist2_corr_filter(sens_acyclic_pre_dist2, sens_pup)
    res_acyclic_corr_pre_dist2 = pre_dist2_corr_filter(res_acyclic_pre_dist2, res_pup)

    # go through the pre_dist relationships and only keep the ones whose linkers
    # are "non-ubiquitous" (i.e. appears less than 2 times as a producer and as produced)
    # return a dataframe containing unprofiled producers
    # and how many edges they yield in pre_dist2

    # start by sorting
    (all_sorted_producers, all_sorted_produced) = sort_sif_list_by_biggest_producer(all_acyclic_corr_pre_dist2,
                                                                                    all_linkers)
    (sens_sorted_producers, sens_sorted_produced) = sort_sif_list_by_biggest_producer(sens_acyclic_corr_pre_dist2,
                                                                                      sens_linkers)
    (res_sorted_producers, res_sorted_produced) = sort_sif_list_by_biggest_producer(res_acyclic_corr_pre_dist2,
                                                                                    res_linkers)

    # i've decided to exclude sif relationships in which the linker
    # takes part as either a producer or produced entity in more than the number of
    # relationships set in linker_lenience (ll)
    # collect the linkers that don't link in excess of linker_lenience
    all_ok_unprofiled_producers = all_sorted_producers[all_sorted_producers[1] < linker_lenience].index
    all_ok_unprofiled_produced = all_sorted_produced[all_sorted_produced[1] < linker_lenience].index

    sens_ok_unprofiled_producers = sens_sorted_producers[sens_sorted_producers[1] < linker_lenience].index
    sens_ok_unprofiled_produced = sens_sorted_produced[sens_sorted_produced[1] < linker_lenience].index

    res_ok_unprofiled_producers = res_sorted_producers[res_sorted_producers[1] < linker_lenience].index
    res_ok_unprofiled_produced = res_sorted_produced[res_sorted_produced[1] < linker_lenience].index

    # get the union of the two lists (for each type in all, sens, and res). Proceed only with those.
    all_ok_unprofiled = all_ok_unprofiled_producers & all_ok_unprofiled_produced
    sens_ok_unprofiled = sens_ok_unprofiled_producers & sens_ok_unprofiled_produced
    res_ok_unprofiled = res_ok_unprofiled_producers & res_ok_unprofiled_produced

    # only include relationships from either where the unprofiled entity isn't too promiscuous
    # write those to file
    write_dist2_out(all_dist2_file, all_acyclic_corr_pre_dist2, all_ok_unprofiled, chebis)
    write_dist2_out(sens_dist2_file, sens_acyclic_corr_pre_dist2, sens_ok_unprofiled, chebis)
    write_dist2_out(res_dist2_file, res_acyclic_corr_pre_dist2, res_ok_unprofiled, chebis)

    all_dist2_file.close()
    sens_dist2_file.close()
    res_dist2_file.close()

    output_locations_all = [output_dir + all_both_sif, output_dir + all_either_sif, output_dir + all_dist2_sif]
    output_locations_sens = [output_dir + sens_both_sif, output_dir + sens_either_sif, output_dir + sens_dist2_sif]
    output_locations_res = [output_dir + res_both_sif, output_dir + res_either_sif, output_dir + res_dist2_sif]

    return output_locations_all, output_locations_sens, output_locations_res


def filter_sif_by_chebi_old(edge, path_to_addl_chebis, output_data_dir):
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
        drop = filter_by_pairwise_corr(members,
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

        # if the produced is in the list of profiled chebis...
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
                drop = filter_by_pairwise_corr([a, c],
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


def filter_by_pairwise_corr(pair,
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
        # drop self-directing edges
        if parent1 == parent2:
            drop = True
        if parent1 != parent2:

            # half of the matrix is NaNs to avoid redundancy
            # if your indexing brings you to this half of the matrix, swap your row & col vals
            if np.isnan(pairwise_spearman_matrix.ix[parent1, parent2]):
                if abs(pairwise_spearman_matrix.ix[parent2, parent1]) < minimum_corr:
                    drop = True

            elif abs(pairwise_spearman_matrix.ix[parent1, parent2]) < minimum_corr:
                drop = True

    return drop


def categorize_eithers_metabs(eithers, chebis):
    """
    Takes a list of .sif relationships and categorizes each node as one of the following:
        profiled_producer
        profiled_produced
        unprofiled_producer
        unprofiled_produced

    Returns [profiled_producer, profiled_produced, unprofiled_producer, unprofiled_produced]
    """

    profiled_producer = []
    profiled_produced = []
    unprofiled_producer = []
    unprofiled_produced = []

    for relationship in eithers:
        [producer, edge, produced] = get_parts_of_sif_line(relationship)
        # if the producer is in the list of profiled chebis...
        if producer in chebis:
            profiled_producer.append(producer)
            unprofiled_produced.append(produced)

        # if the produced is in teh list of profiled chebis...
        if produced in chebis:
            profiled_produced.append(produced)
            unprofiled_producer.append(producer)

    return profiled_producer, profiled_produced, unprofiled_producer, unprofiled_produced


def find_linkers(unprofiled_producers, unprofiled_produced):
    """
    Identifies unprofiled metabolites that serve to connect two profiled metabolites.
    unprofiled_producers is a list of unprofiled entities that produce a profiled metabolite
    unprofiled_produced is a list of unprofiled entities that are produced by a profiled metabolite
    """
    # prepare to capture unprofiled entities which are
    # produced by a profiled entity. Then evaluate whether
    # they also produce a profiled entity.
    linkers = []
    for u_pd in unprofiled_produced:
        if u_pd in unprofiled_producers:
            linkers.append(u_pd)

    # remove duplicates
    # note that almost every unprofiled produced is also an unprofiled producer
    # means almost all are middlemen (that's ok)
    linkers = list(set(linkers))

    return linkers


def who_produces_whom(eithers, linkers):
    """
    Stores information regarding profiledA --> unprofiledB --> profiledC relationships.

    produces_linker = {linker: [list of metabs that produce it] ... }
    produced_by_linker = {linker: [list of metabs it produces] ... }

    Returns produces_linker, produced_by_linker
    """

    # initialize empty dictionaries
    produces_linker = {key: [] for key in linkers}
    produced_by_linker = {key: [] for key in linkers}

    for relationship in eithers:
        [producer, edge, produced] = get_parts_of_sif_line(relationship)
        if producer in linkers:
            produced_by_linker[producer].append(produced)
        if produced in linkers:
            produces_linker[produced].append(producer)

    return produces_linker, produced_by_linker


def make_profiled_unprofiled_profiled_dict(produces_linker,
                                           produced_by_linker,
                                           name_ids_dict,
                                           corr_matrix,
                                           min_corr):
    """
    Generates a dictionary to describe Profiled_a --> Unprofiled_b --> Profiled_c
    Only fills it with the relationships where A and C exhibit sufficient correlation
    """
    pup = {}
    for linker, its_producers in produces_linker.items():
        pup[linker] = []
        for a in its_producers:
            for c in produced_by_linker[linker]:
                drop = filter_by_pairwise_corr([a, c],
                                               name_ids_dict,
                                               corr_matrix,
                                               min_corr)
                if not drop:
                    pup[linker].append((a,c))

    return pup

def build_pre_distance_of_2(eithers, linkers, profiled_producer, profiled_produced):
    """
    Builds a first pass of the distance-of-2 network.
    Does so by generating the following dictionaries
        dict1
            profiled producer: produced
        dict2
            profiled produced: producer
    In other words:
        dict1 k = profiled producer
        dict1 v = unprofiled produced
        dict2 k = unprofiled producer
        dict2 v = profiled produced

    Returns a list of distance-of-two relationships as well as dict1 and dict2

    """
    # TODO: make less cryptic
    dict1 = {}
    dict2 = {}

    # initialize pre_dist2
    pre_dist2 = []

    # only add those to pre_dist2 if the unprofiled entity is also a linker
    for rel in eithers:
        [producer, edge, produced] = get_parts_of_sif_line(rel)

        # if producer is profiled and produced is a linker, keep it
        if producer in profiled_producer:
            if produced in linkers:
                # make all profiled producers keys in dict1
                # append all their unprofiled produced entities as vals
                if producer in dict1:
                    dict1[producer].append(produced)
                if producer not in dict1:
                    dict1[producer] = [produced]
                # keeping the whole line
                pre_dist2.append(rel)

        # if produced is profiled and producer is a linker, keep it
        if produced in profiled_produced:
            if producer in linkers:
                # make all unprofiled producers keys in dict2
                # append all their profiled produced entities as vals
                if producer in dict2:
                    dict2[producer].append(produced)
                if producer not in dict2:
                    dict2[producer] = [produced]
                # keeping the whole line
                pre_dist2.append(rel)

    return pre_dist2, dict1, dict2


def find_unprofiled_cycles(dict1, dict2):
    """
    Identifies cycles wherein one profiled entity is linked to itself by an unprofiled entity
    input dictionaries (see build_pre_distance_of_2):
        dict1 k = profiled producer
        dict1 v = unprofiled produced
        dict2 k = unprofiled producer
        dict2 v = profiled produced
    """
    cycles = []

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

    return cycles

def remove_cycles(pre_dist2, cycles):
    """
    Removes cyclic relationships (as defined in find_unprofiled_cycles) from predist2.
    Returns acycic_predist2
    """
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

    return acyclic_predist2


def pre_dist2_corr_filter(acyclic_predist2, pup):
    """
    Determine which edges to keep in distance of 2 sif based on correlation of A and C
        In example profiledA --> unprofiledB --> profiledC
        Recall that pup[B] = [(A,C), ... ]
    """
    acyc_corr_predist2 = []
    for rel in acyclic_predist2:
        rel = rel.strip('\n')
        [producer, edge, produced] = get_parts_of_sif_line(rel)
        # start by assuming the edge will be dropped
        drop = True

        # if the producer is a linker (unprofiled),
        if producer in pup:
            B = producer    # for the sake of clarity
            for pair in pup[B]:
                # if the second member of any of these tuples (A, C) matches produced
                # that means the correlation between those two profiled compounds (A and C)
                # was deemed sufficient to maintain an edge between B (producer) and C (produced)
                if pair[1] == produced:
                    drop = False

        # if the produced entity is a linker (unprofiled),
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

    return acyc_corr_predist2


def write_dist2_out(dist2_file, acyc_corr_predist2, ok_unprofiled, chebis):
    """
    Writes out distance_of_2 sif including only the linkers which have passed all filtration steps.

    """
    for rel in acyc_corr_predist2:
        [producer, edge, produced] = get_parts_of_sif_line(rel)

        if producer in chebis:
            if produced in ok_unprofiled:
                dist2_file.write(rel + '\n')
        if produced in chebis:
            if producer in ok_unprofiled:
                dist2_file.write(rel + '\n')


def remove_weakly_correlated(all_sif_path, input_dir):
    """
    Parses .sifs generated from complete (non subsetted) data to remove
    nodes which exhibited negligible change between sensitive and resistant
    mean fold change (mean = geometric mean).
    Writes out a file of the same name with "correlated" at the beginning of its name

    """
    # load assay data
    assay_df = pd.read_csv(input_dir + 'assay_results_extended.tsv',
                                sep='\t',
                                index_col=0,
                                header=0,
                                na_values = 'nd')

    # identify metabolites whose fold change from sensitive to resistant geometric means
    # is too weak to include in the "all" graphs.
    weak_metabs = assay_df[abs(assay_df['fc_gmeans']) < 0.1].index.tolist()

    # load sif:
    in_sif_fh = open(all_sif_path, 'r')
    in_sif = in_sif_fh.readlines()

    # open outfile
    out_sif_path = all_sif_path.replace('all','all_correlated')
    out_sif_fh = open(out_sif_path, 'w')

    # iterate through the sif and throw out any edges which contain a weak_metab
    for relationship in in_sif:
        [producer, edge, produced] = get_parts_of_sif_line(relationship)
        if not any(x in weak_metabs for x in [producer, produced]):
            out_sif_fh.write(relationship)

    in_sif_fh.close()
    out_sif_fh.close()

    os.remove(all_sif_path)

    return out_sif_path





