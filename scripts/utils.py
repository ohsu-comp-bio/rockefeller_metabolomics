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

def name_outfile_path(data_dir, min_corr1, min_corr2, max_pval):
    """
    Generates the name of an outfile.
    """

    def float_to_str(flt):
        if flt == 0.0:
            return "00"
        elif flt == 1.0:
            return "100"
        else:
            newstring = str(flt)
            newstring = newstring.strip("0")
            if "." in newstring:
                parts = newstring.split(".")
                newstring = parts[1]
            if len(newstring) == 1:
                newstring = newstring + "0"
            return newstring

    mc1 = float_to_str(min_corr1)
    mc2 = float_to_str(min_corr2)
    mpv = float_to_str(max_pval)

    out_path = data_dir + "min_corrs_" + str(mc1) + "_" + str(mc2) + "_" + "max_p_" + str(mpv) + '/'

    run_tag = mc1 + "_" + mc2 + "_" + mpv + "_"

    return out_path, run_tag