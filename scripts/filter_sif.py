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


def filter_sif_by_chebi(sif_name, chebis, name_ids_dict, min_corr_both, min_corr_dist2,
                        linker_lenience, input_dir, output_dir, run_tag):
    """
    Takes a list of all profiled or derived chebis.

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

    returns filenames
    """
    # read in used-to-produce.sif to be filtered down
    sif_fh = open(input_dir + sif_name, 'r')
    sif = sif_fh.readlines()

    # open the outfiles and store their handles in dict
    # i.e. {'sens_chebi_EITHER_members.sif': <_io.TextIOWrapper
    # name='/Users/manningh/Desktop/safebox/playtime/sens_chebi_EITHER_members.sif' mode='w' encoding='UTF-8'>}
    outf_dict = {}
    for subset in ['all', 'res', 'sens']:
        for filt in ['BOTH', 'EITHER', 'DIST2']:
            sif_name = subset + '_chebi_' + filt + '_members.sif'
            fh = open(output_dir + run_tag + sif_name, "w")
            outf_dict[sif_name] = fh

    # prepare to collect .sif relationships in which either member is profiled
    eithers = []

    # read in the metabolite-by-metabolite spearman correlation matrices
    all_corr_matrix = pd.read_csv(input_dir + 'all_pairwise_spearman_corr_coef.tsv',
                                           sep='\t',
                                           index_col=0)

    # begin filtration
    for relationship in sif:
        [producer, edge, produced] = get_parts_of_sif_line(relationship)

        # drop a relationship if both members are profiled and the correlation is low
        drop = filter_by_pairwise_corr([producer, produced],
                                           name_ids_dict,
                                           all_corr_matrix,
                                           min_corr_both)

        # if the overall pair has a decent correlation
        # build BOTH and EITHER networks in one iteration
        if not drop:

            # generate a sif where 1 member must be profiled to maintain the relationship:
            if any(smallmol in chebis for smallmol in [producer, produced]):
                for file_name, fh in outf_dict.items():
                    if 'EITHER' in file_name:
                        fh.write(relationship)
                        eithers.append(relationship)

            # generate a sif where both members must be profiled to maintain the relationship
            # write it also to the distance-of-2 network file
            if all(smallmol in chebis for smallmol in [producer, produced]):
                for file_name, fh in outf_dict.items():
                    if any(filt_type in file_name for filt_type in ('BOTH','DIST2')):
                        fh.write(relationship)

    # The BOTH and EITHER sifs are now complete (with ChEBI IDs, not names yet)
    # close used-to-produce.sif
    sif_fh.close()

    # close the BOTH and EITHER .sif outfiles
    for file_name, fh in outf_dict.items():
        # if any(filt_type in file_name for filt_type in ('BOTH','EITHER')):
        # close all and reopen DIST2 later
        fh.close()

    # Now build the remainder of the distance-of-2 networks through further iteration
    # store relationship members based on their role and profiled status
    (profiled_producer,
     profiled_produced,
     unprofiled_producer,
     unprofiled_produced) = \
        categorize_eithers_metabs(eithers, chebis)

    # identify unprofiled "linkers":
    # capture unprofiled entities which are produced by a profiled entity. Then evaluate
    # whether they also produce a profiled entity.
    linkers = find_linkers(unprofiled_producer, unprofiled_produced)

    # keep track of who produces whom in linker-involved relationships
    (produces_linker, produced_by_linker) = who_produces_whom(eithers, linkers)

    # Make a dictionary where keys are linkers and values are lists of tuples
    # so the relationship ProfiledA --> UnprofiledB --> ProfiledC is represented: pup[B] = [(A,C)]
    pup = make_profiled_unprofiled_profiled_dict(produces_linker,
                                                 produced_by_linker,
                                                 name_ids_dict,
                                                 all_corr_matrix,
                                                 min_corr_dist2)

    # TODO: make less cryptic
    # prepare a first pass for collecting distance-of-2 relationships
    # add relationships to pre_dist2 if not a both relationship and if unprofiled entity is a
    # correlation-approved linker
    (pre_dist2, dict1, dict2) = build_pre_distance_of_2(eithers,
                                                        linkers,
                                                        profiled_producer,
                                                        profiled_produced)

    # identify cycles wherein one profiled entity is linked to itself by an unprofiled entity
    cycles = find_unprofiled_cycles(dict1, dict2)

    # remove cycles from pre_dist2 lists
    acyclic_pre_dist2 = remove_cycles(pre_dist2, cycles)

    # a final correlation check is required for distance of two relationships
    acyclic_corr_pre_dist2 = pre_dist2_corr_filter(acyclic_pre_dist2, pup)

    # go through the pre_dist relationships and only keep the ones whose linkers
    # are "non-ubiquitous" (i.e. appears less than 2 times as a producer and as produced)
    # return a dataframe containing unprofiled producers
    # and how many edges they yield in pre_dist2

    # start by sorting
    (sorted_producers, sorted_produced) = sort_sif_list_by_biggest_producer(acyclic_corr_pre_dist2,
                                                                                    linkers)
    '''
    # for debugging
    sorted_producers.to_csv(output_dir + "sorted_producers.tsv", sep='\t')
    sorted_produced.to_csv(output_dir + "sorted_produced.tsv", sep='\t')
    '''
    # todo: determine why ok_unprofiled producers and produced aren't getting anything
    # todo: determine why sorted_producers and sorted_produced now don't have any linkers that link less than 4
    # i've decided to exclude sif relationships in which the linker
    # takes part as either a producer or produced entity in more than the number of
    # relationships set in linker_lenience (ll)
    # collect the linkers that don't link in excess of linker_lenience
    ok_unprofiled_producers = sorted_producers[sorted_producers[1] < linker_lenience].index
    ok_unprofiled_produced = sorted_produced[sorted_produced[1] < linker_lenience].index

    # get the union of the two lists (for each type in all, sens, and res). Proceed only with those.
    ok_unprofiled = ok_unprofiled_producers & ok_unprofiled_produced

    # only include relationships from either where the unprofiled entity isn't too promiscuous
    # write those to file
    for file_name, fh in outf_dict.items():
        if 'DIST2' in file_name:
            fh = open(output_dir + run_tag + file_name, "a")
            for rel in acyclic_corr_pre_dist2:
                [producer, edge, produced] = get_parts_of_sif_line(rel)
                if producer in chebis:
                    if produced in ok_unprofiled:
                        fh.write(rel + '\n')
                if produced in chebis:
                    if producer in ok_unprofiled:
                        fh.write(rel + '\n')
            # close the dist2 files
            fh.close()

    return [output_dir + run_tag + filename for filename in outf_dict.keys()]


def sort_sif_list_by_biggest_producer(sif_list, middlemen):
    """
    Take a list of node-edge-node relationships (as in either_sif_fh.readlines())
    and a list of "middlemen" in this .sif.

    middlemen are unprofiled entities which appear as both producers and produced
    entities in the sif.

    Returns a sorted pd.DataFrame of unprofiled producers and the
    number of node-edge-node relationships for which they are the first node, same for produced
    index = ChEBI IDs, single col = number of edges for which it is a producer
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



