
import numpy as np
from itertools import groupby
from src.pack_graph.op_seq import reverse_complement_seq, is_kmer_similar


class Diag:
    def __init__(self, x_start, x_end, y_start, y_end, orient, offset):

        self.x_start = x_start
        self.x_end = x_end
        self.y_start = y_start
        self.y_end = y_end

        self.orient = orient
        self.offset = offset

        self.true_reverse = False

    def __eq__(self, other):
        if not isinstance(other, Diag):
            return False

        return self.x_start == other.x_start and self.x_end == other.x_end and self.y_start == other.y_start and self.y_end == other.y_end and self.orient == other.orient

    def to_string(self):

        return "{}-{}, {}-{}, {}, {}".format(self.x_start, self.x_end, self.y_start, self.y_end, self.orient, self.offset)


def find_line_diag(matrix, min_line_len, return_format="dict", flip_lr=False, ignore_offset=None):

    rows, cols = matrix.shape

    if return_format == "dict":
        candidate_line_diags = {}
    else:
        candidate_line_diags = []

    for offset in range(-rows + 1, cols):

        if offset == ignore_offset:
            continue

        diag = matrix.diagonal(offset)

        current_pointer = 0

        for val, _ in groupby(diag):
            group_len = len(list(_))
            group_start = current_pointer
            group_end = current_pointer + group_len - 1

            if val == 1 and group_len >= min_line_len:

                if offset <= 0:
                    if flip_lr:
                        diag_obj = Diag(cols - (group_end + 1), cols - group_start, group_start + abs(offset), group_end + 1 + abs(offset), "reverse", offset)
                    else:
                        diag_obj = Diag(group_start, group_end + 1, group_start + abs(offset), group_end + 1 + abs(offset), "forward", offset)
                else:

                    if flip_lr:
                        diag_obj = Diag(cols - (group_end + 1 + offset), cols - (group_start + offset), group_start, group_end + 1, "reverse", offset)
                    else:
                        diag_obj = Diag(group_start + offset, group_end + 1 + offset, group_start, group_end + 1, "forward", offset)

                if return_format == "dict":
                    if offset not in candidate_line_diags:
                        candidate_line_diags[offset] = []

                    candidate_line_diags[offset].append(diag_obj)

                else:
                    candidate_line_diags.append(diag_obj)

            current_pointer = group_end + 1

    # if return_format == "dict":

        # for diag in candidate_line_diags[423]:
        #     print("candi", diag.to_string())

        # neighbor_thresh = 1
        # for major_offset in range(0, cols, neighbor_thresh + 1):
        #     linear_diags_by_shift(major_offset, candidate_line_diags, neighbor_thresh, min_line_len, rows, cols)
        #
        # for major_offset in range(-rows + 1, 0, neighbor_thresh + 1):
        #     linear_diags_by_shift(major_offset, candidate_line_diags, neighbor_thresh, min_line_len, rows, cols)

        # for diag in candidate_line_diags[423]:
        #     print("after candi", diag.to_string())

    return candidate_line_diags


def linear_diags_by_shift(major_offset, candidate_line_diags, neighbor_thresh, min_line_len,  n_rows, n_cols):
    shift1_offset = major_offset + neighbor_thresh
    shift2_offset = major_offset - neighbor_thresh

    # # STEP: do self linear
    for offset in [major_offset, shift1_offset, shift2_offset]:
        if offset not in candidate_line_diags:
            continue

        if offset != major_offset and offset == 0:
            continue
        offset_include_diags = candidate_line_diags[offset]

        if len(offset_include_diags) > 1:
            will_delete_diags = []

            current_diag = offset_include_diags[0]
            next_diag_pointer = 1

            while next_diag_pointer < len(offset_include_diags):
                next_diag = offset_include_diags[next_diag_pointer]

                if next_diag.y_start - current_diag.y_end <= min_line_len:
                    will_delete_diags.append(next_diag)

                    current_diag.x_start = min(current_diag.x_start, next_diag.x_start)
                    current_diag.y_start = min(current_diag.y_start, next_diag.y_start)

                    current_diag.x_end = max(current_diag.x_end, next_diag.x_end)
                    current_diag.y_end = max(current_diag.y_end, next_diag.y_end)

                else:
                    current_diag = next_diag

                next_diag_pointer += 1

            for diag in will_delete_diags:
                offset_include_diags.remove(diag)

    # # STEP: do shift linear
    if major_offset in candidate_line_diags:
        major_offset_include_diags = candidate_line_diags[major_offset]

        for major_diag in major_offset_include_diags:

            for shift_offset in [shift1_offset, shift2_offset]:
                if shift_offset not in candidate_line_diags:
                    continue

                shift_offset_include_diags = candidate_line_diags[shift_offset]

                will_delete_diags = []

                for shift_diag in shift_offset_include_diags:

                    # # sort on reads (y)
                    # # target diag is the latter one
                    if shift_diag.y_start >= major_diag.y_end:
                        distance_on_y = shift_diag.y_start - major_diag.y_end
                        distance_on_x = shift_diag.x_start - major_diag.x_end if major_diag.orient == "forward" else major_diag.x_start - shift_diag.x_end

                    # # base diag in the latter one
                    elif major_diag.y_start >= shift_diag.y_end:
                        distance_on_y = major_diag.y_start - shift_diag.y_end
                        distance_on_x = major_diag.x_start - shift_diag.x_end if major_diag.orient == "forward" else shift_diag.x_start - major_diag.x_end

                    # # base and target are overlapped on y
                    else:
                        continue

                    # linear_thresh = max([major_diag_len, shift_diag_len])
                    linear_thresh = min_line_len

                    if distance_on_y <= linear_thresh and distance_on_x <= linear_thresh:
                        major_diag.x_start = min([shift_diag.x_start, major_diag.x_start])
                        major_diag.y_start = min([shift_diag.y_start, major_diag.y_start])
                        major_diag.x_end = max([shift_diag.x_end, major_diag.x_end])
                        major_diag.y_end = max([shift_diag.y_end, major_diag.y_end])

                        major_diag_y_len = major_diag.y_end - major_diag.y_start
                        major_diag_x_len = major_diag.x_end - major_diag.x_start

                        if major_diag_x_len != major_diag_y_len:

                            if major_diag_x_len > major_diag_y_len:

                                diff_length = major_diag_x_len - major_diag_y_len

                                if major_diag.y_end + diff_length >= n_rows:
                                    major_diag.x_end -= diff_length
                                else:
                                    major_diag.y_end += diff_length

                            else:
                                diff_length = major_diag_y_len - major_diag_x_len

                                if major_diag.orient == "forward":
                                    if major_diag.x_end + diff_length >= n_cols:
                                        major_diag.y_end -= diff_length
                                    else:
                                        major_diag.x_end += diff_length
                                else:
                                    if major_diag.x_start - diff_length < 0:
                                        major_diag.y_start += diff_length
                                    else:
                                        major_diag.x_start -= diff_length

                        will_delete_diags.append(shift_diag)

                for diag in will_delete_diags:
                    shift_offset_include_diags.remove(diag)

        # for offset in candidate_line_diags:
        #     offset_include_diags = candidate_line_diags[offset]
        #
        #     # print(offset, [diag.to_string() for diag in offset_include_diags])
        #
        #     if len(offset_include_diags) > 1:
        #         will_delete_diags = []
        #
        #         current_diag = offset_include_diags[0]
        #         next_diag_pointer = 1
        #
        #         while next_diag_pointer < len(offset_include_diags):
        #             next_diag = offset_include_diags[next_diag_pointer]
        #
        #             if next_diag.y_start - current_diag.y_end <= min_line_len:
        #                 will_delete_diags.append(next_diag)
        #
        #                 current_diag.x_start = min(current_diag.x_start, next_diag.x_start)
        #                 current_diag.y_start = min(current_diag.y_start, next_diag.y_start)
        #
        #                 current_diag.x_end = max(current_diag.x_end, next_diag.x_end)
        #                 current_diag.y_end = max(current_diag.y_end, next_diag.y_end)
        #
        #             else:
        #                 current_diag = next_diag
        #
        #             next_diag_pointer += 1
        #
        #         for diag in will_delete_diags:
        #             offset_include_diags.remove(diag)


def boost_diags(dotplot_type, unique_diags, x2y_matrix_len_x, x2y_matrix_len_y, raw_project_x, raw_project_x_rev, augment_coeff):


    # if dotplot_type == "ref2alt":
    # unique_diags = sorted(unique_diags, key=lambda x: (x.orient, x.y_start))
    # else:
    #     unique_diags = sorted(unique_diags, key=lambda x: (x.orient, x.y_start), reverse=True)
    # unique_diags = sorted(unique_diags, key=lambda x: (x.y_end - x.y_start), reverse=True)[: 8]

    unique_diags = sorted(unique_diags, key=lambda x: (x.orient, x.y_start))

    first_diag = sorted([diag for diag in unique_diags if diag.orient == "forward"], key=lambda x: (x.x_start + x.y_start))[0]
    last_diag = sorted([diag for diag in unique_diags if diag.orient == "forward"], key=lambda x: (x.x_end + x.y_end))[-1]
    other_diags = sorted([diag for diag in unique_diags if diag not in [first_diag, last_diag]], key=lambda x: (x.y_end - x.y_start), reverse=True)


    # print()
    # print(dotplot_type, x2y_matrix_len_x, x2y_matrix_len_y,)
    # for diag in unique_diags:
    #     print("final", diag.to_string())

    # # STEP: update the projection, and ignore the overlap on y
    # x2y_matrix_len_x = x2y_matrix.shape[1]
    # x2y_matrix_len_y = x2y_matrix.shape[0]

    y_boost_flag = np.zeros(x2y_matrix_len_y)

    # for diag in unique_diags:
    for diag in [first_diag, last_diag] + other_diags:
        # print("final", diag.to_string())

        # print(diag.to_string())
        diag_x_positions = np.array(range(diag.x_start, diag.x_end))
        diag_x_positions_rev = np.flip(diag_x_positions)

        diag_y_positions = np.array(range(diag.y_start, diag.y_end))

        allowed_boost_y_positions = np.intersect1d(diag_y_positions,  np.where(y_boost_flag != 1))

        if diag.orient == "forward":
            raw_project_x[diag_x_positions[np.where(np.in1d(diag_y_positions, allowed_boost_y_positions))]] += augment_coeff

            # print(np.where(np.in1d(diag_y_positions, allowed_boost_y_positions)))
            # if len(diag_x_positions[np.where(np.in1d(diag_y_positions, allowed_boost_y_positions))]) > 0:
            #     print(min(diag_x_positions[np.where(np.in1d(diag_y_positions, allowed_boost_y_positions))]), max(diag_x_positions[np.where(np.in1d(diag_y_positions, allowed_boost_y_positions))]))
            # else:
            #     print([])
        else:

            raw_project_x[diag_x_positions_rev[np.where(np.in1d(diag_y_positions, allowed_boost_y_positions))]] += augment_coeff

            # print(diag_x_positions_rev[np.where(np.in1d(diag_y_positions, allowed_boost_y_positions))])
            # if len(diag_x_positions_rev[np.where(np.in1d(diag_y_positions, allowed_boost_y_positions))]) > 0:
            #     print(min(diag_x_positions_rev[np.where(np.in1d(diag_y_positions, allowed_boost_y_positions))]), max(diag_x_positions_rev[np.where(np.in1d(diag_y_positions, allowed_boost_y_positions))]))
            # else:
            #     print([])

            # # this is a reverse diag, if it also exist in the rev matrix, then we consider it
            if diag.true_reverse:
                raw_project_x_rev[diag_x_positions_rev[np.where(np.in1d(diag_y_positions, allowed_boost_y_positions))]] += augment_coeff

        y_boost_flag[allowed_boost_y_positions] = 1

    # # The last position, due to the potentially truncated sequence, is commonly unmatched,
    # # therefore, we maintained the last with the previous one
    forward_diags = sorted([diag for diag in unique_diags if diag.orient == "forward"], key=lambda x: x.x_start)

    if not len(forward_diags) == 0:
        start_diag = forward_diags[0]
        start_trim_len = start_diag.x_start + 2

        raw_project_x[0: start_trim_len] = raw_project_x[start_trim_len]
        raw_project_x_rev[0: start_trim_len] = raw_project_x_rev[start_trim_len]

        forward_diags = sorted(forward_diags, key=lambda x: x.x_end)

        end_diag = forward_diags[-1]
        end_trim_len = x2y_matrix_len_x - end_diag.x_end + 2

        raw_project_x[-end_trim_len: -1] = raw_project_x[-end_trim_len]
        raw_project_x_rev[-end_trim_len: -1] = raw_project_x_rev[-end_trim_len]

        # print(x2y_matrix_len_x, x2y_matrix_len_y, start_trim_len, end_trim_len)
    # for diag in forward_diags:
    #     print("\t", diag.to_string())

    # raw_project_x[-1] = raw_project_x[-2]
    # raw_project_x_rev[-1] = raw_project_x_rev[-2]


def find_and_denoise_diags(x2x_matrix, x2y_matrix, x2y_matrix_rev, stride_size):
    """

    """

    x2y_matrix_len_x = x2y_matrix.shape[1]
    x2y_matrix_len_y = x2y_matrix.shape[0]

    max_shift_len = 50 / stride_size

    x2x_matrix_fliplr = np.fliplr(x2x_matrix)
    x2y_matrix_fliplr = np.fliplr(x2y_matrix)

    # print(stride_size, 30 / stride_size, 50 / stride_size)
    x2x_matrix_line_diags = find_line_diag(x2x_matrix, min_line_len=max(1, 30 / stride_size), return_format="dict", ignore_offset=0)
    x2x_matrix_line_diags_fliplr = find_line_diag(x2x_matrix_fliplr, min_line_len=max(1, 30 / stride_size), flip_lr=True, return_format="dict", ignore_offset=0)

    # for diag in x2x_matrix_line_diags:
    #     print(diag.to_string())

    # # allocate diags into bins to accelerate the compare process
    x2x_bin_size = 100
    x2x_matrix_line_diags_bins = {}

    for offset in x2x_matrix_line_diags:
        for diag in x2x_matrix_line_diags[offset]:
            diag_bin = int(diag.x_start / x2x_bin_size)

            if diag_bin not in x2x_matrix_line_diags_bins:
                x2x_matrix_line_diags_bins[diag_bin] = []
            x2x_matrix_line_diags_bins[diag_bin].append(diag)

    for offset in x2x_matrix_line_diags_fliplr:
        for diag in x2x_matrix_line_diags_fliplr[offset]:
            diag_bin = int(diag.x_start / x2x_bin_size)

            if diag_bin not in x2x_matrix_line_diags_bins:
                x2x_matrix_line_diags_bins[diag_bin] = []
            x2x_matrix_line_diags_bins[diag_bin].append(diag)

    unique_diags = []

    for mode in ["raw", "fliplr"]:

        if mode == "raw":
            x2y_matrix_line_diags = find_line_diag(x2y_matrix, min_line_len=max(5, 50 / stride_size), return_format="dict")
        else:
            x2y_matrix_line_diags = find_line_diag(x2y_matrix_fliplr, min_line_len=max(5, 50 / stride_size), flip_lr=True, return_format="dict")

        for offset in x2y_matrix_line_diags:

            for diag in x2y_matrix_line_diags[offset]:

                # # for the first and last match diag, they are main bone of the image, so we directly keep them
                if ((diag.x_start in [0] and diag.y_start in [0])
                        or (diag.x_end in [x2y_matrix_len_x - 1, x2y_matrix_len_x] and diag.y_end in [x2y_matrix_len_y - 1, x2y_matrix_len_y])):
                    unique_diags.append(diag)
                    continue

                diag_start, diag_end = diag.x_start, diag.x_end

                match_flag = False

                diag_bin = int(diag.x_start / x2x_bin_size)
                if diag_bin in x2x_matrix_line_diags_bins:
                    for target_diag in x2x_matrix_line_diags_bins[diag_bin]:
                        if target_diag.x_start - max_shift_len <= diag_start <= diag_end <= target_diag.x_end + max_shift_len:
                        # if abs(target_diag.x_start - diag_start) <= max_shift_len and abs(diag_end - target_diag.x_end) <= max_shift_len:
                            match_flag = True
                            break

                if match_flag is False:
                    unique_diags.append(diag)

    # # STEP: sort
    unique_diags = sorted(unique_diags, key=lambda x: x.y_start)

    # # STEP: for the reversed diags in x2y_matrix, check if they are really in the x2y_matrix_rev
    x2y_matrix_rev_line_diags = find_line_diag(np.fliplr(x2y_matrix_rev), min_line_len=max(5, 50 / stride_size), flip_lr=True, return_format="dict")

    x2y_matrix_rev_line_diags_bins = {}
    for offset in x2y_matrix_rev_line_diags:
        for diag in x2y_matrix_rev_line_diags[offset]:
            diag_bin = int(diag.x_start / x2x_bin_size)

            if diag_bin not in x2y_matrix_rev_line_diags_bins:
                x2y_matrix_rev_line_diags_bins[diag_bin] = []
            x2y_matrix_rev_line_diags_bins[diag_bin].append(diag)

    for diag in unique_diags:
        if diag.orient == "forward":
            continue

        diag_bin = int(diag.x_start / x2x_bin_size)
        if diag_bin in x2y_matrix_rev_line_diags_bins:
            for target_diag_rev in x2y_matrix_rev_line_diags_bins[diag_bin]:
                if diag == target_diag_rev:
                    diag.true_reverse = True

    if len(unique_diags) == 0:
        return unique_diags

    # # remove overlapped with the left and right-most diags (anchor diags)
    first_diag = sorted([diag for diag in unique_diags if diag.orient == "forward"], key=lambda x: (x.x_start + x.y_start))[0]
    last_diag = sorted([diag for diag in unique_diags if diag.orient == "forward"], key=lambda x: (x.x_end + x.y_end))[-1]

    overlapped_diags = []
    for diag in unique_diags:
        if diag == first_diag or diag == last_diag:
            continue

        if diag.orient == "forward" and (diag.y_start < first_diag.y_end or diag.y_end > last_diag.y_start):
            overlapped_diags.append(diag)

    for diag in overlapped_diags:
        unique_diags.remove(diag)

    # print()
    # for diag in unique_diags:
    #     print('raw', diag.to_string())

    return unique_diags


def segment_projections_into_matrix(x2x_dotplot_project_x, x2y_dotplot_project_x, x2x_dotplot_project_x_rev, x2y_dotplot_project_x_rev, ideal_x2y_project_type, augment_coeff):
    """


    """

    projection_matrix = []  # # format [[segment_start, segment_end, segment_x2x, segment_x2y, segment_x2y_rev, segment_label], [], [], ....]

    dotplot_x2y_project_subtract = x2x_dotplot_project_x - x2y_dotplot_project_x
    dotplot_x2y_project_rev_subtract = x2x_dotplot_project_x_rev - x2y_dotplot_project_x_rev

    # # segment using iteration
    previous_vals = [augment_coeff, dotplot_x2y_project_subtract[0], dotplot_x2y_project_rev_subtract[0]]    # # [x2x_val, x2y_val, x2y_val_rev]
    previous_pointer = 0

    for i in range(1, len(x2x_dotplot_project_x) + 1):

        if i == len(x2x_dotplot_project_x):
            new_vals = [-1, -1, -1]     # # meet the end
        else:
            new_vals = [augment_coeff, dotplot_x2y_project_subtract[i], dotplot_x2y_project_rev_subtract[i]]

        # # find a segment
        if new_vals != previous_vals or i == len(x2x_dotplot_project_x):

            segment_start = previous_pointer
            segment_end = i

            projection_matrix.append([segment_start, segment_end, previous_vals[0], previous_vals[1], previous_vals[2], "None"])

            previous_vals = new_vals
            previous_pointer = i

    # # STEP: fix
    # print("{}\t{}\t{}".format(np.shape(projection_matrix)[0], np.shape(x2x_dotplot_project_x)[0], round(np.shape(projection_matrix)[0] / np.shape(x2x_dotplot_project_x)[0], 3)))
    # print(np.shape(projection_matrix)[0] / np.shape(x2x_dotplot_project_x)[0])

    # # the proportion of projection to the raw dotplot, this is also the compression number
    if np.shape(projection_matrix)[0] / np.shape(x2x_dotplot_project_x)[0] > 0.8:
        projection_matrix = [[0, len(x2x_dotplot_project_x), augment_coeff, 0, 0, 0, "None"]]
        return projection_matrix, True      #  # bad flag as True

    return projection_matrix, False
