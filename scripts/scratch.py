def sort_sif_by_biggest_producer():
    dist_2_sif_fh = open(output_data_dir + 'named-filtered-used-to-produce-DIST-OF-2.sif', 'r')
    dist_2_sif = dist_2_sif_fh.readlines()

    producer_counts = {}
    for line in dist_2_sif:
        line = line.strip('\n')
        line = line.split('\t')
        producer = line[0]
        if producer not in producer_counts:
            producer_counts[producer] = 1
        if producer in producer_counts:
            producer_counts[producer] += 1

    # generate a list of tuples (chemical entity, number of times it appears as a producer in sif)
    producer_counts_sorted = sorted(producer_counts.items(), key=operator.itemgetter(1), reverse=True)



