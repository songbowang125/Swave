import numpy as np

import itertools


def find_continuous_val(val_list):
    """

    find continuous val in a list or string.

    return format: [ [val, [indexes in raw list]], [], ..... ]

    input example:  AAABBBCCAAA
    output example : [['A', [0, 1, 2]], ['B', [3, 4, 5]], ['C', [6, 7]], ['A', [8, 9, 10]]]
    """

    continuous_res = []

    current_pointer = 0

    for val, group in itertools.groupby(val_list):
        group_len = len(list(group))

        group_start = current_pointer               # # open interval
        group_end = current_pointer + group_len     # # closed interval

        continuous_res.append([val, [i for i in range(group_start, group_end)]])

        # # update pointer
        current_pointer = group_end

    return continuous_res


def calculate_seq_similarity_larger_than(thresh, seq_orient, seq1, seq2, extend_orient):
    from src.pack_dotplot.op_dotplot import Dotplot

    # dotplot = Dotplot(seq1, seq2, 10, extend_orient, stride_size=None, given_x_hash_table=None)
    if seq_orient == "forward":
        dotplot = Dotplot(seq1, seq2, 10, str(len(seq1)) + str(len(seq2)) + extend_orient, stride_size=None, given_x_hash_table=None, skip_reverse=True)
    else:
        dotplot = Dotplot(seq1, seq2, 10, str(len(seq1)) + str(len(seq2)) + extend_orient, stride_size=None, given_x_hash_table=None, skip_forward=True)

    # dotplot.to_png(out_img=True)

    detail_stride_size = dotplot.stride_size

    projection_x, _ = dotplot.get_project_x()
    projection_x_len = len(projection_x)

    projection_y, _ = dotplot.get_project_y()
    # projection_y_len = len(projection_y)

    projection_x_summary = find_continuous_val(projection_x)
    # print(projection_x_summary)
    if extend_orient == "left":
        for value, include_indexes in projection_x_summary:
            if value == 0:
                continue

            index_start = include_indexes[0]
            remaining_len = projection_x_len - index_start

            if remaining_len == 0:
                return 0

            # print(index_start, np.count_nonzero(projection_x[index_start: ]) / remaining_len, np.count_nonzero(projection_y[: remaining_len]) / remaining_len )
            if np.count_nonzero(projection_x[index_start: ]) / remaining_len >= thresh:

                # # check the project in y
                if seq_orient == "forward" and np.count_nonzero(projection_y[-remaining_len:]) / remaining_len >= thresh:
                    return remaining_len * detail_stride_size

                if seq_orient == "reverse" and np.count_nonzero(projection_y[: remaining_len]) / remaining_len >= thresh:
                    return remaining_len * detail_stride_size

    if extend_orient == "right":
        projection_x_summary.reverse()
        for value, include_indexes in projection_x_summary:
            if value == 0:
                continue

            index_end = include_indexes[-1]

            remaining_len = index_end

            if remaining_len == 0:
                return 0

            if np.count_nonzero(projection_x[: index_end + 1]) / remaining_len >= thresh:
                # # check the project in y
                if seq_orient == "forward" and np.count_nonzero(projection_y[: remaining_len]) / remaining_len >= thresh:
                    return remaining_len * detail_stride_size
                if seq_orient == "reverse" and np.count_nonzero(projection_y[-remaining_len:]) / remaining_len >= thresh:
                    return remaining_len * detail_stride_size

    return 0


def calculate_seq_similarity(seq1, seq2):

    if len(seq1) == 0 or len(seq2) == 0:
        return 0

    # # using pairwise align for sequence similarity
    # return max(alignment[2] / alignment[4] for alignment in pairwise2.align.globalxx(seq1, seq2))
    # return max(alignment[2] for alignment in pairwise2.align.globalxx(seq1, seq2)) / min(len(seq1), len(seq2))

    # # using dotplot projection for sequence similarity
    from src.pack_dotplot.op_dotplot import Dotplot

    # # use the shorter seq as x_seq, and projection to the shorter seq
    if len(seq1) < len(seq2):
        dotplot = Dotplot(seq1, seq2, 10, "tmp", stride_size=1, given_x_hash_table=None)
    else:
        dotplot = Dotplot(seq2, seq1, 10, "tmp", stride_size=1, given_x_hash_table=None)

    projection_x, _ = dotplot.get_project_x()
    projection_x_rev = dotplot.get_project_x_rev()

    return np.count_nonzero(projection_x + projection_x_rev) / len(projection_x)


def calculate_seq_repeat_ratio(seq1):
    from src.pack_dotplot.op_dotplot import Dotplot

    # # use the shorter seq as x_seq, and projection to the shorter seq

    if len(seq1) == 0:
        return np.inf

    # allowed_seq_len = 10000
    # seq1 = seq1[: allowed_seq_len]

    # dotplot = Dotplot(seq1, seq1, 10, "tmp", stride_size=1, given_x_hash_table=None)
    dotplot = Dotplot(seq1, seq1, 10, "tmp", stride_size=None, given_x_hash_table=None)

    # dotplot.to_png(out_img=True)
    project_y, _ = dotplot.get_project_y()
    return np.average(project_y)


def is_kmer_similar(kmer1, kmer2, mismatch_thresh=0):

    mismatch_cnt = 0

    for i in range(len(kmer1)):

        kmer1_base = kmer1[i]
        kmer2_base = kmer2[i]

        if kmer1_base != kmer2_base:
            mismatch_cnt += 1

    if mismatch_cnt <= mismatch_thresh:
        return True
    else:
        return False


def reverse_complement_seq(seq):
    """
    given a seq, get its reverse complementary seq
    """
    inv_seq = ""
    for i in range(len(seq) - 1, -1, -1):
        bp = seq[i]

        if bp == 'A' or bp == "a":
            inv_bp = 'T'
        elif bp == 'T' or bp == "t":
            inv_bp = 'A'
        elif bp == 'C' or bp == "c":
            inv_bp = 'G'
        elif bp == 'G' or bp == "g":
            inv_bp = 'C'
        else:
            inv_bp = 'N'

        inv_seq += inv_bp

    return inv_seq

