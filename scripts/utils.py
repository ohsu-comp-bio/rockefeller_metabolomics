"""
Basic utility functions.
"""


def add_lines_to_list(file_path, base_list):
    """
    Appends every line in the file specified by file_path to
    the list provided.

    Returns the list
    """
    fh = open(file_path, "r")
    lines = fh.readlines()
    for line in lines:
        line = line.strip('\n')
        base_list.append(line)

    return base_list

def get_parts_of_sif_line(line):
    """
    Returns entity1, edge, and entity2 from a .sif relationship
    """
    line = line.strip('\n')
    parts = line.split('\t')
    entity1 = parts[0]
    edge = parts[1]
    entity2 = parts[2]

    return entity1, edge, entity2