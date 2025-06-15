import logging
import numpy as np
import os
import cv2 as cv
import math
from src.pack_graph.op_seq import reverse_complement_seq
from matplotlib import pyplot as plt
from src.pack_dotplot.op_projection import find_and_denoise_diags, segment_projections_into_matrix, Diag, boost_diags
from src.pack_graph.op_seq import calculate_seq_similarity_larger_than, is_kmer_similar


class Dotplot:
    def __init__(self, seq_x, seq_y, kmer_size, out_prefix, against="auto", stride_size=None, given_x_hash_table=None, given_y_hash_table=None, skip_forward=False, skip_reverse=False):

        self.seq_x = seq_x.upper()
        self.seq_y = seq_y.upper()

        self.seq_x_len = len(seq_x)
        self.seq_y_len = len(seq_y)

        self.kmer_size = kmer_size

        self.skip_forward = skip_forward
        self.skip_reverse = skip_reverse

        if stride_size is None:
            self.stride_size = calculate_stride_size(self.seq_x, self.seq_y)
        else:
            self.stride_size = stride_size

        self.out_prefix = out_prefix

        self.matrix = np.zeros((int(self.seq_y_len / self.stride_size) + 1, int(self.seq_x_len / self.stride_size) + 1))
        self.matrix_rev = np.zeros((int(self.seq_y_len / self.stride_size) + 1, int(self.seq_x_len / self.stride_size) + 1))

        # if self.seq_x_len >= 50000 or self.seq_y_len >= 50000:
        #     self.create_matrix_split()
        # else:

        if against == "auto":
            if self.seq_x_len >= self.seq_y_len:
                self.create_matrix_against_x(given_x_hash_table)
            else:
                self.create_matrix_against_y(given_y_hash_table)

        elif against == "x":
            self.create_matrix_against_x(given_x_hash_table)

        elif against == "y":
            self.create_matrix_against_y(given_y_hash_table)

        else:
            logging.error("No such against axis: {}. Choose from [auto, x, y]".format(against))

        self.matrix = self.matrix[:-1, :-1]
        self.matrix_rev = self.matrix_rev[:-1, :-1]

        # self.seq_x_hash_table = None
        # self.seq_y_hash_table = None

    def rotate_to_alt2ref(self):
        self.matrix = np.fliplr(np.rot90(self.matrix, k=-1))
        self.matrix_rev = np.fliplr(np.rot90(self.matrix_rev, k=-1))

        # # switch seq x and y
        tmp_seq_x = self.seq_x
        self.seq_x = self.seq_y
        self.seq_y = tmp_seq_x

        # # switch seq x and y
        tmp_seq_x_len = self.seq_x_len
        self.seq_x_len = self.seq_y_len
        self.seq_y_len = tmp_seq_x_len

    def get_project_x(self, augment=False):

        project_x = np.sum(self.matrix, axis=0)

        augment_coeff = int(100 * np.average(project_x))

        if augment:
            project_x[np.diag(self.matrix) == 1] += augment_coeff

        return project_x, augment_coeff

    def get_project_x_rev(self, baseline=0):

        """
        Why we need base line: not all position of x2y is reversely aligned. But all positions in x2x is augmented.
                                therefore, we need use the x2x augment value as the base line
        """
        project_x_rev = baseline + np.sum(self.matrix_rev, axis=0)

        return project_x_rev

    def get_project_y(self, augment=False):

        project_y = np.sum(self.matrix, axis=1)

        augment_coeff = int(100 * np.average(project_y))

        if augment:
            project_y[np.diag(self.matrix) == 1] += augment_coeff

        return project_y, augment_coeff

    def get_project_y_rev(self, baseline=0):

        project_y_rev = baseline + np.sum(self.matrix_rev, axis=1)

        return project_y_rev

    def create_matrix_against_x(self, given_x_hash_table=None):

        if given_x_hash_table is None:
            self.seq_x_hash_table = HashTable(self.seq_x, kmer_size=self.kmer_size, stride_size=self.stride_size)
        else:
            self.seq_x_hash_table = given_x_hash_table

        # # the much detaier way
        pos_on_y = 0
        while pos_on_y < self.seq_y_len:
            kmer_str = self.seq_y[pos_on_y: pos_on_y + self.kmer_size]

            index_on_y = int(pos_on_y / self.stride_size)

            # # for original kmer
            if not self.skip_forward:
                indexes_on_x = self.seq_x_hash_table.find_all(kmer_str)
                if indexes_on_x is not None:
                    for index_on_x in indexes_on_x:
                        self.matrix[index_on_y, index_on_x] = 1

            # # for reversed kmer
            if not self.skip_reverse:
                kmer_str_reversed = reverse_complement_seq(kmer_str)

                indexes_on_x = self.seq_x_hash_table.find_all(kmer_str_reversed)

                if indexes_on_x is not None:
                    for index_on_x in indexes_on_x:
                        self.matrix[index_on_y, index_on_x] = 1
                        self.matrix_rev[index_on_y, index_on_x] = 1

            pos_on_y += self.stride_size

    def create_matrix_against_y(self, given_y_hash_table=None):

        if given_y_hash_table is None:
            self.seq_y_hash_table = HashTable(self.seq_y, kmer_size=self.kmer_size, stride_size=self.stride_size)
        else:
            self.seq_y_hash_table = given_y_hash_table

        # # the much detaier way
        pos_on_x = 0
        while pos_on_x < self.seq_x_len:

            kmer_str = self.seq_x[pos_on_x: pos_on_x + self.kmer_size]

            index_on_x = int(pos_on_x / self.stride_size)

            # # for original kmer
            if not self.skip_forward:
                indexes_on_y = self.seq_y_hash_table.find_all(kmer_str)

                if indexes_on_y is not None:
                    for index_on_y in indexes_on_y:
                        self.matrix[index_on_y, index_on_x] = 1

            # # for reversed kmer
            if not self.skip_reverse:
                kmer_str_reversed = reverse_complement_seq(kmer_str)

                indexes_on_y = self.seq_y_hash_table.find_all(kmer_str_reversed)

                if indexes_on_y is not None:
                    for index_on_y in indexes_on_y:
                        self.matrix[index_on_y, index_on_x] = 1
                        self.matrix_rev[index_on_y, index_on_x] = 1

            pos_on_x += self.stride_size

    # def create_matrix(self, given_x_hash_table=None):
    #
    #     if given_x_hash_table is None:
    #         self.seq_x_hash_table = HashTable(self.seq_x, kmer_size=self.kmer_size, stride_size=self.stride_size)
    #     else:
    #         self.seq_x_hash_table = given_x_hash_table
    #
    #     # # the much detaier way
    #     # pos_on_y = 0
    #     # while pos_on_y < self.seq_y_len:
    #     #
    #     #     index_on_y = int(pos_on_y / self.stride_size)
    #     #
    #     #     for i in range(self.stride_size):
    #     #         kmer_str = self.seq_y[pos_on_y + i: pos_on_y + i + self.kmer_size]
    #     #
    #     #         # # for original kmer
    #     #         indexes_on_x = self.seq_x_hash_table.find_all(kmer_str)
    #     #
    #     #         if indexes_on_x is not None:
    #     #             for index_on_x in indexes_on_x:
    #     #                 self.matrix[index_on_y, index_on_x] = 1
    #     #
    #     #         # # for reversed kmer
    #     #         kmer_str_reversed = reverse_complement_seq(kmer_str)
    #     #
    #     #         indexes_on_x = self.seq_x_hash_table.find_all(kmer_str_reversed)
    #     #
    #     #         if indexes_on_x is not None:
    #     #             for index_on_x in indexes_on_x:
    #     #                 self.matrix[index_on_y, index_on_x] = 1
    #     #                 self.matrix_rev[index_on_y, index_on_x] = 1
    #     #
    #     #     pos_on_y += self.stride_size
    #
    #
    #     # # the raw  way
    #     pos_on_y = 0
    #     while pos_on_y < self.seq_y_len:
    #
    #         kmer_str = self.seq_y[pos_on_y: pos_on_y + self.kmer_size]
    #
    #         index_on_y = int(pos_on_y / self.stride_size)
    #
    #         # # for original kmer
    #         indexes_on_x = self.seq_x_hash_table.find_all(kmer_str)
    #
    #         if indexes_on_x is not None:
    #             for index_on_x in indexes_on_x:
    #                 self.matrix[index_on_y, index_on_x] = 1
    #
    #         # # for reversed kmer
    #         kmer_str_reversed = reverse_complement_seq(kmer_str)
    #
    #         indexes_on_x = self.seq_x_hash_table.find_all(kmer_str_reversed)
    #
    #         if indexes_on_x is not None:
    #             for index_on_x in indexes_on_x:
    #                 self.matrix[index_on_y, index_on_x] = 1
    #                 self.matrix_rev[index_on_y, index_on_x] = 1
    #
    #         pos_on_y += self.stride_size

    # def create_matrix_split(self, split_window_size=15000):
    #
    #     max_window_num = round(min(self.seq_x_len, self.seq_y_len) / split_window_size)
    #
    #     for window_index in range(max_window_num):
    #
    #         # # get the split seq x and seq y by the window index
    #         window_start_x = window_index * split_window_size
    #         window_end_x = (window_index + 1) * split_window_size
    #
    #         window_start_y = window_index * split_window_size
    #         window_end_y = (window_index + 1) * split_window_size
    #
    #         if window_index == (max_window_num - 1):    # # last window
    #             window_end_x = self.seq_x_len
    #             window_end_y = self.seq_y_len
    #
    #         split_seq_x = self.seq_x[window_start_x: window_end_x]
    #         split_seq_y = self.seq_y[window_start_y: window_end_y]
    #
    #         split_seq_x_len = len(split_seq_x)
    #         split_seq_y_len = len(split_seq_y)
    #
    #         # # begin to generate
    #         split_seq_x_hash_table = HashTable(split_seq_x, kmer_size=self.kmer_size, stride_size=self.stride_size)
    #
    #         pos_on_split_y = 0
    #
    #         while pos_on_split_y < split_seq_y_len:
    #             kmer_str = split_seq_y[pos_on_split_y: pos_on_split_y + self.kmer_size]
    #
    #             index_on_split_y = int(pos_on_split_y / self.stride_size)
    #
    #             # # for original kmer
    #             indexes_on_split_x = split_seq_x_hash_table.find_all(kmer_str)
    #
    #             if indexes_on_split_x is not None:
    #                 for index_on_split_x in indexes_on_split_x:
    #                     self.matrix[index_on_split_y + int(window_start_y / self.stride_size), index_on_split_x + int(window_start_x / self.stride_size)] = 1
    #
    #             # # for reversed kmer
    #             kmer_str_reversed = reverse_complement_seq(kmer_str)
    #
    #             indexes_on_split_x = split_seq_x_hash_table.find_all(kmer_str_reversed)
    #
    #             if indexes_on_split_x is not None:
    #                 for index_on_split_x in indexes_on_split_x:
    #                     self.matrix[index_on_split_y + int(window_start_y / self.stride_size), index_on_split_x + int(window_start_x / self.stride_size)] = 1
    #                     self.matrix_rev[index_on_split_y + int(window_start_y / self.stride_size), index_on_split_x + int(window_start_x / self.stride_size)] = 1
    #
    #             pos_on_split_y += self.stride_size

    def get_seq_x_hash_table(self):
        return self.seq_x_hash_table

    def to_png(self, reverse=False, out_img=False):
        self.dotplot_file = self.out_prefix + ".dotplot.png"

        if reverse:

            if np.max(self.matrix_rev) == 0:
                matrix_resize_norm = 255 * np.ones(np.shape(self.matrix_rev))

            else:
                matrix_resize_norm = 255 * abs(self.matrix_rev - np.max(self.matrix_rev)) / (np.max(self.matrix_rev) - np.min(self.matrix_rev))
        else:

            if np.max(self.matrix) == 0:
                matrix_resize_norm = 255 * np.ones(np.shape(self.matrix))

            else:
                matrix_resize_norm = 255 * abs(self.matrix - np.max(self.matrix)) / (np.max(self.matrix) - np.min(self.matrix))

        if out_img:
            cv.imwrite(self.dotplot_file, matrix_resize_norm)

        return matrix_resize_norm

    def to_csv(self):
        self.matrix_file = self.out_prefix + ".matrix.txt"

        np.savetxt(self.matrix_file, self.matrix[self.seq_x_len:, ], delimiter="\t", fmt='%d')


class HashTable:

    def __init__(self, seq, kmer_size, stride_size):

        self.seq = seq
        # self.seq_reversed = reverse_complement_seq(self.seq)

        self.hash_table = {}

        # # STEP: generate hash table for original seq, a much detailer way by traversing the stride window
        pos_on_seq = 0
        while pos_on_seq < len(seq):

            index_on_seq = int(pos_on_seq / stride_size)

            for i in range(stride_size):

                kmer_str = seq[pos_on_seq + i: pos_on_seq + i + kmer_size]

                kmer_hash = self.calculate_hash_value(kmer_str)

                if self.hash_table.get(kmer_hash) is None:
                    self.hash_table[kmer_hash] = []

                self.hash_table[kmer_hash].append(index_on_seq)

            pos_on_seq += stride_size

        # # # STEP: generate hash table for original seq, a easy way by the whole stride window
        # pos_on_seq = 0
        # while pos_on_seq < len(seq):
        #
        #     kmer_str = seq[pos_on_seq: pos_on_seq + kmer_size]
        #
        #     index_on_seq = int(pos_on_seq / stride_size)
        #
        #     kmer_hash = self.calculate_hash_value(kmer_str)
        #
        #     if self.hash_table.get(kmer_hash) is None:
        #         self.hash_table[kmer_hash] = []
        #
        #     self.hash_table[kmer_hash].append(index_on_seq)
        #
        #     pos_on_seq += stride_size

    def calculate_hash_value(self, seq_str):

        seq_hash = seq_str

        return seq_hash

    def find_all(self, seq_str):

        seq_hash = self.calculate_hash_value(seq_str)

        return self.hash_table.get(seq_hash)


def calculate_stride_size(seq_x, seq_y):
    # # automatedly calculate the stride size
    # # we aim to make the size of dotplot at [150, 15000].
    # # why 150: the min_sv_length is 50, and the expand length is equal to 50, therefore, the minimum seq length is 150
    max_seq_len = max([len(seq_x), len(seq_y)])

    if max_seq_len <= 1500:
        stride_size = 1

    elif max_seq_len <= 15000:
        stride_size = 10
    # else:
    #     stride_size = 30
    elif max_seq_len <= 75000:
        stride_size = 30

    elif max_seq_len <= 100000:
        stride_size = 50

    else:
        stride_size = 100

    return stride_size


def output_synthesized_dotplot(output_prefix, dotplot_matrix):

    dotplot_file = output_prefix + ".dotplot.png"
    # cv.imwrite(dotplot_file, x2y)
    #
    # if x2y_rev is not None:
    #     cv.imwrite(dotplot_file, np.concatenate((x2x, x2y, x2y_rev), axis=0))
    # else:
    #     cv.imwrite(dotplot_file, np.concatenate((x2x, x2y), axis=0))

    cv.imwrite(dotplot_file, dotplot_matrix)


def output_synthesized_dotplot_project(dotplot_output_prefix, x2x_dotplot_project_x, x2y_dotplot_project_x, x2x_dotplot_project_x_rev, x2y_dotplot_project_x_rev):

    plt.figure(figsize=(10, 10))
    ax1 = plt.subplot(611)
    ax1.plot(x2x_dotplot_project_x)
    ax1.set_xticks([])

    ax2 = plt.subplot(612)
    ax2.plot(x2y_dotplot_project_x)
    ax2.set_xticks([])

    ax3 = plt.subplot(613)
    ax3.plot(x2x_dotplot_project_x - x2y_dotplot_project_x)
    ax3.set_xticks([])

    ax4 = plt.subplot(614)
    ax4.plot(x2x_dotplot_project_x_rev)
    ax4.set_xticks([])

    ax5 = plt.subplot(615)
    ax5.plot(x2y_dotplot_project_x_rev)
    ax5.set_xticks([])

    ax6 = plt.subplot(616)
    ax6.plot(x2x_dotplot_project_x_rev - x2y_dotplot_project_x_rev)

    plt.savefig(dotplot_output_prefix + ".proj_x.png")
    plt.close()


# def detail_compare_kmer(kmer_x, kmer_y):
#
#     detail_kmer_size = 10
#     detail_kmer_stride = 1
#
#     kmer_x_hash = HashTable(kmer_x, detail_kmer_size, detail_kmer_stride)
#
#     pos_on_y = 0
#
#     while pos_on_y < len(kmer_y):
#
#         sub_kmer_str_y = kmer_y[pos_on_y: pos_on_y + detail_kmer_size]
#
#         # # for original kmer
#         indexes_on_x = kmer_x_hash.find_all(sub_kmer_str_y)
#
#         if indexes_on_x is not None:
#
#             return True
#
#         # # for reversed kmer
#         kmer_str_reversed = reverse_complement_seq(sub_kmer_str_y)
#
#         indexes_on_x = kmer_x_hash.find_all(kmer_str_reversed)
#
#         if indexes_on_x is not None:
#             return True
#
#         pos_on_y += detail_kmer_stride
#
#     return False


def base_level_dotplot_extension(unique_diags, x_seq, y_seq, x2y_matrix, raw_project_x, raw_project_x_rev, raw_kmer_size, stride_size):
    # print()
    # for diag in unique_diags:
    #     print("before base", diag.to_string())

    # # STEP: set param
    x2y_matrix_shape = x2y_matrix.shape
    x2y_matrix_len_x = x2y_matrix_shape[1]
    x2y_matrix_len_y = x2y_matrix_shape[0]

    re_check_kmer_size = 10
    re_check_length_start = 1 * raw_kmer_size
    re_check_length_end = 1 * raw_kmer_size

    # # STEP: perform base level boost
    for i in range(len(unique_diags)):
        diag = unique_diags[i]
        try:
            if i == 0:
                max_augment_start = np.inf
                max_augment_end = unique_diags[i + 1].x_start - unique_diags[i].x_end

            elif i == len(unique_diags) - 1:
                max_augment_start = unique_diags[i].x_start - unique_diags[i - 1].x_end
                max_augment_end = np.inf
            else:
                max_augment_start = unique_diags[i].x_start - unique_diags[i - 1].x_end
                max_augment_end = unique_diags[i + 1].x_start - unique_diags[i].x_end

            # # < 0 means there is a candidate overlap between current diag and its previous or later one
            if max_augment_start < 0:
                max_augment_start = np.inf
            if max_augment_end < 0:
                max_augment_end = np.inf

        except:
            max_augment_start = np.inf
            max_augment_end = np.inf

        # # STEP: adjust the coordinate by the orientataion of the diag
        if diag.orient == "forward":
            diag_x_start = diag.x_start
            diag_x_end = diag.x_end
            diag_y_start = diag.y_start
            diag_y_end = diag.y_end

        else:
            diag_x_start = diag.x_start
            diag_x_end = diag.x_end
            diag_y_start = x2y_matrix_len_y - diag.y_end - 1
            diag_y_end = x2y_matrix_len_y - diag.y_start - 1

            y_seq = reverse_complement_seq(y_seq)

        # # STEP: deal the start seq of the diag line
        if not (diag_x_start - 1 < 0 and diag_y_start - 1 < 0):

            diag_start_seq_x = x_seq[diag_x_start * stride_size - re_check_length_start: diag_x_start * stride_size]
            diag_start_seq_y = y_seq[diag_y_start * stride_size - re_check_length_start: diag_y_start * stride_size]

            re_check_length_start = min(re_check_length_start, len(diag_start_seq_x), len(diag_start_seq_y))

            diag_start_match_flag = False

            # # traverse on seq_x
            max_match_length, index_on_x, index_on_y = 0, 0, 0
            for index_on_x in range(re_check_length_start, re_check_kmer_size, -1):
                # # traverse on seq_y
                for index_on_y in range(re_check_length_start, re_check_kmer_size, -1):
                    # # why: the sub kmer comparison mush start at the left or end, we cannot compare inside the kmer
                    if not (index_on_x == re_check_length_start or index_on_y == re_check_length_start):
                        continue

                    diag_start_seq_x_sub_seq = diag_start_seq_x[index_on_x - re_check_kmer_size: index_on_x]
                    diag_start_seq_y_sub_seq = diag_start_seq_y[index_on_y - re_check_kmer_size: index_on_y]

                    # # match the sub seq
                    if is_kmer_similar(diag_start_seq_x_sub_seq, diag_start_seq_y_sub_seq, mismatch_thresh=1):
                        diag_start_match_flag = True

                        # # find the max match length
                        max_match_length += re_check_kmer_size
                        max_match_flag = "meet_sub_seq_end"  # # choose from meet_sub_seq_end or meet_mis_match
                        for i in range(min(index_on_x, index_on_y)):

                            if diag_start_seq_x[index_on_x - re_check_kmer_size - i] == diag_start_seq_y[index_on_y - re_check_kmer_size - i]:
                                max_match_length += 1
                            else:
                                max_match_flag = "meet_mis_match"
                                break

                        # # update the diag
                        if max_match_flag == "meet_mis_match":
                            augment_length = min(math.ceil(max_match_length / stride_size), max_augment_start)
                            # augment_length = min(round(max_match_length / stride_size), max_augment_start)

                        else:
                            # # if meet_sub_seq_end, then the whole sub seq is matched although it is not long enough, so we consider this as the whole match of the sub seq
                            augment_length = min(math.ceil(re_check_length_start / stride_size), max_augment_start)
                            # augment_length = min(round(re_check_seq_length / stride_size), max_augment_start)

                        for i in range(augment_length):
                            if diag.x_start - 1 >= 0:

                                # # the start is a open-interval, therefore, we first -=1 then update the project
                                diag.x_start -= 1
                                diag.y_start -= 1

                                raw_project_x[diag.x_start] += 1
                                if diag.orient == "reverse":
                                    raw_project_x_rev[diag.x_start] += 1

                        break

                if diag_start_match_flag is True:
                    break

        if not (diag_x_end + 1 >= x2y_matrix_len_x and diag_y_end + 1 >= x2y_matrix_len_y):

            diag_end_seq_x = x_seq[diag_x_end * stride_size: diag_x_end * stride_size + re_check_length_end]
            diag_end_seq_y = y_seq[diag_y_end * stride_size: diag_y_end * stride_size + re_check_length_end]

            re_check_length_end = min(re_check_length_end, len(diag_end_seq_x), len(diag_end_seq_y))
            # print(diag_end_seq_x, diag_end_seq_y, re_check_length_end)

            diag_end_match_flag = False

            # # traverse on seq_x
            max_match_length, index_on_x, index_on_y = 0, 0, 0
            for index_on_x in range(re_check_length_end - re_check_kmer_size):
                # # traverse on seq_y
                for index_on_y in range(re_check_length_end - re_check_kmer_size):

                    # # why: the sub kmer comparison mush start at the left or end, we cannot compare inside the kmer
                    if not (index_on_x == 0 or index_on_y == 0):
                        continue

                    diag_end_seq_x_sub_seq = diag_end_seq_x[index_on_x: index_on_x + re_check_kmer_size]
                    diag_end_seq_y_sub_seq = diag_end_seq_y[index_on_y: index_on_y + re_check_kmer_size]

                    # # match the sub seq
                    if is_kmer_similar(diag_end_seq_x_sub_seq, diag_end_seq_y_sub_seq, mismatch_thresh=1):
                        diag_end_match_flag = True

                        # # find the max match length
                        max_match_length += re_check_kmer_size
                        max_match_flag = "meet_sub_seq_end"  # # choose from meet_sub_seq_end or meet_mis_match

                        for i in range(re_check_length_end - max(index_on_x, index_on_y) - re_check_kmer_size):
                            if diag_end_seq_x[index_on_x + re_check_kmer_size + i] == diag_end_seq_y[index_on_y + re_check_kmer_size + i]:
                                max_match_length += 1
                            else:
                                max_match_flag = "meet_mis_match"
                                break

                        # # update the diag
                        if max_match_flag == "meet_mis_match":
                            augment_length = min(math.ceil(max_match_length / stride_size), max_augment_end)
                            # augment_length = min(round(max_match_length / stride_size), max_augment_end)

                        else:
                            # # if meet_sub_seq_end, then the whole sub seq is matched although it is not long enough, so we consider this as the whole match of the sub seq
                            augment_length = min(math.ceil(re_check_length_end / stride_size), max_augment_end)
                            # augment_length = min(round(re_check_seq_length / stride_size), max_augment_end)

                        # # update the diag and boost the raw project
                        for i in range(augment_length):
                            if diag.x_end < x2y_matrix_len_x:

                                # # the end is a close-interval, therefore, we first update the project then +=1
                                raw_project_x[diag.x_end] += 1
                                if diag.orient == "reverse":
                                    raw_project_x_rev[diag.x_end] += 1

                                diag.x_end += 1
                                diag.y_end += 1

                        break

                if diag_end_match_flag is True:
                    break

    return unique_diags


def linear_level_dotplot_extension(unique_diags, x2y_matrix):
    x2y_matrix_shape = x2y_matrix.shape
    x2y_matrix_len_x = x2y_matrix_shape[1]
    x2y_matrix_len_y = x2y_matrix_shape[0]

    # # STEP: sort
    unique_diags = sorted(unique_diags, key=lambda x: x.y_start)

    #
    # # STEP: check linear diags, threahold dot number is min_sv_size / stride_size
    # # why do this: SNPs would case break of diags, which would be predicted as del + inv
    linear_diags = []

    for i in range(len(unique_diags) - 1, 0, -1):

        current_diag = unique_diags[i]
        current_diag_len = current_diag.y_end - current_diag.y_start

        previous_diag = unique_diags[i - 1]
        previous_diag_len = previous_diag.y_end - previous_diag.y_start

        if current_diag.orient != previous_diag.orient:
            continue

        # # the first diag, the major anchor one, does not involve in this process
        if (current_diag.x_start == 0 and current_diag.y_start == 0) or (previous_diag.x_start == 0 and previous_diag.y_start == 0):
            continue
        # if (current_diag.x_start == 0 and current_diag.y_start == 0) or (previous_diag.x_start == 0 and previous_diag.y_start == 0):
        #     linear_thresh = max(5, 50 / stride_size)
        # else:

        linear_thresh = max([current_diag_len, previous_diag_len])

        if current_diag.orient == "reverse" and ((current_diag.y_start - previous_diag.y_end) == (previous_diag.x_start - current_diag.x_end)) and (current_diag.y_start - previous_diag.y_end) <= linear_thresh:

            previous_diag.x_start = min([previous_diag.x_start, current_diag.x_start])
            previous_diag.y_start = min([previous_diag.y_start, current_diag.y_start])
            previous_diag.x_end = max([previous_diag.x_end, current_diag.x_end])
            previous_diag.y_end = max([previous_diag.y_end, current_diag.y_end])

            # previous_diag_x_length = previous_diag.x_end - previous_diag.x_start
            # previous_diag_y_length = previous_diag.y_end - previous_diag.y_start
            #
            # if previous_diag_x_length > previous_diag_y_length:
            #     previous_diag_length_diff = previous_diag_x_length - previous_diag_y_length
            #     # previous_diag.y_end += previous_diag_length_diff
            #     previous_diag.x_end -= previous_diag_length_diff
            #
            # elif previous_diag_x_length < previous_diag_y_length:
            #     previous_diag_length_diff = previous_diag_y_length - previous_diag_x_length
            #     # previous_diag.x_end += previous_diag_length_diff
            #     previous_diag.y_end -= previous_diag_length_diff

            linear_diags.append(current_diag)

        if current_diag.orient == "forward" and ((current_diag.y_start - previous_diag.y_end) == (current_diag.x_start - previous_diag.x_end)) and (current_diag.y_start - previous_diag.y_end) <= linear_thresh:

            previous_diag.x_start = min([previous_diag.x_start, current_diag.x_start])
            previous_diag.y_start = min([previous_diag.y_start, current_diag.y_start])
            previous_diag.x_end = max([previous_diag.x_end, current_diag.x_end])
            previous_diag.y_end = max([previous_diag.y_end, current_diag.y_end])

            linear_diags.append(current_diag)

    for diag in linear_diags:
        unique_diags.remove(diag)

    # print()
    # for diag in unique_diags:
    #     print('after', diag.to_string())

    # # remove diags that are fully covered by others
    unique_diags = sorted(unique_diags, key=lambda x: (x.y_end - x.y_start), reverse=True)

    full_covered_diags = []
    for i in range(len(unique_diags) - 1, -1, -1):

        for j in range(i - 1, -1, -1):
            base_diag = unique_diags[i]

            if ((base_diag.x_start in [0] and base_diag.y_start in [0])
                    or (base_diag.x_end in [x2y_matrix_len_x - 1, x2y_matrix_len_x] and base_diag.y_end in [x2y_matrix_len_y - 1, x2y_matrix_len_y])):
                continue

            target_diag = unique_diags[j]

            if target_diag.y_start - 1 <= base_diag.y_start <= base_diag.y_end <= target_diag.y_end + 1:
                full_covered_diags.append(base_diag)
                break

    for diag in full_covered_diags:
        unique_diags.remove(diag)

    return unique_diags


def kmer10_level_dotplot_extension(unique_diags, x_seq, y_seq, x2y_matrix, stride_size, thresh=0.5):
    """
    this is for the diag extension in repetitive regions.

    Although there might be no exact match in repetitive region, we seek for similar match using detail dotplot (with small kmers and strides)
    """
    # unique_diags = sorted(unique_diags, key=lambda x: (x.y_start + x.x_start))
    unique_diags = sorted(unique_diags, key=lambda x: x.y_start)
    # unique_diags = sorted(unique_diags, key=lambda x: (x.y_start + x.y_end) / 2)

    # print()
    # print()
    # print(np.shape(x2y_matrix), len(x_seq), len(y_seq))
    #
    x2y_matrix_shape = np.shape(x2y_matrix)
    x2y_matrix_shape_x = x2y_matrix_shape[1]
    x2y_matrix_shape_y = x2y_matrix_shape[0]

    unique_diags.insert(0, Diag(0, 0, 0, 0, "forward", 0))
    unique_diags.append(Diag(x2y_matrix_shape_x, x2y_matrix_shape_x, x2y_matrix_shape_y, x2y_matrix_shape_y, "forward", 0))

    # print()
    # for diag in unique_diags:
    #     print("before kmer10", diag.to_string())

    for diag_index in range(1, len(unique_diags) - 1):
        current_diag = unique_diags[diag_index]

        if current_diag.x_start == 0 and current_diag.y_start == 0:
            continue

        # # STEP: deal with the previous gap
        previous_diag = unique_diags[diag_index - 1]
        gap_length = current_diag.y_start - previous_diag.y_end

        # print()
        # print("start", current_diag.to_string(), stride_size)
        # print("previous: ", previous_diag.to_string(), "gap_length", gap_length)

        if gap_length > 0:
            gap_start_on_y = previous_diag.y_end
            gap_end_on_y = current_diag.y_start

            if current_diag.orient == "forward":
                gap_start_on_x = max(current_diag.x_start - gap_length, 0)
                gap_end_on_x = current_diag.x_start
            else:
                gap_start_on_x = current_diag.x_end
                gap_end_on_x = min(current_diag.x_end + gap_length, x2y_matrix_shape_x)

            # if gap_start_on_y < 0 or gap_start_on_x < 0 or gap_end_on_x >= x2y_matrix_shape_x or gap_end_on_y >= x2y_matrix_shape_y:
            #     continue

            # gap_seq_on_y = y_seq[gap_start_on_y * stride_size: (gap_end_on_y + 1) * stride_size]
            # gap_seq_on_x = x_seq[gap_start_on_x * stride_size: (gap_end_on_x + 1) * stride_size]
            gap_seq_on_y = y_seq[gap_start_on_y * stride_size: gap_end_on_y * stride_size]
            gap_seq_on_x = x_seq[gap_start_on_x * stride_size: gap_end_on_x * stride_size]

            # print(gap_start_on_x, gap_end_on_x, gap_start_on_y, gap_end_on_y)

            if gap_start_on_y == gap_end_on_y or gap_start_on_x == gap_end_on_x:
                continue

            if current_diag.orient == "forward":
                similar_len = calculate_seq_similarity_larger_than(thresh, current_diag.orient, gap_seq_on_x, gap_seq_on_y, "left")
            else:
                similar_len = calculate_seq_similarity_larger_than(thresh, current_diag.orient, gap_seq_on_x, gap_seq_on_y, "right")

            # print(similar_len)

            if current_diag.orient == "forward":
                current_diag.x_start -= int(similar_len / stride_size)
                current_diag.y_start -= int(similar_len / stride_size)
            else:
                current_diag.x_end += int(similar_len / stride_size)
                current_diag.y_start -= int(similar_len / stride_size)

        # # STEP: deal with the latter gap
        latter_diag = unique_diags[diag_index + 1]
        gap_length = latter_diag.y_start - current_diag.y_end
        # print("latter ", latter_diag.to_string(), "gap_length", gap_length)

        if gap_length > 0:

            gap_start_on_y = current_diag.y_end
            gap_end_on_y = latter_diag.y_start

            if current_diag.orient == "forward":
                gap_start_on_x = current_diag.x_end
                gap_end_on_x = min(current_diag.x_end + gap_length, x2y_matrix_shape_x)
            else:
                gap_start_on_x = max(current_diag.x_start - gap_length, 0)
                gap_end_on_x = current_diag.x_start

            # gap_seq_on_y = y_seq[gap_start_on_y * stride_size: (gap_end_on_y + 1) * stride_size]
            # gap_seq_on_x = x_seq[gap_start_on_x * stride_size: (gap_end_on_x + 1) * stride_size]
            gap_seq_on_y = y_seq[gap_start_on_y * stride_size: gap_end_on_y * stride_size]
            gap_seq_on_x = x_seq[gap_start_on_x * stride_size: gap_end_on_x * stride_size]

            # print(gap_start_on_x, gap_end_on_x, gap_start_on_y, gap_end_on_y)

            if gap_start_on_y == gap_end_on_y or gap_start_on_x == gap_end_on_x:
                continue

            if current_diag.orient == "forward":
                similar_len = calculate_seq_similarity_larger_than(thresh, current_diag.orient, gap_seq_on_x, gap_seq_on_y, "right")
            else:
                similar_len = calculate_seq_similarity_larger_than(thresh, current_diag.orient, gap_seq_on_x, gap_seq_on_y, "left")

            # print(len(gap_seq_on_x), len(gap_seq_on_y), similar_len)

            if current_diag.orient == "forward":
                current_diag.x_end += int(similar_len / stride_size)
                current_diag.y_end += int(similar_len / stride_size)
            else:
                current_diag.x_start -= int(similar_len / stride_size)
                current_diag.y_end += int(similar_len / stride_size)

    #     print("end", current_diag.to_string())
    #     # # detail dotplot on the front side
    #     # if diag
    # print("after")
    # for diag in unique_diags:
    #     print(diag.to_string())

    return unique_diags[1: -1]


def rotate_diags_to_alt2ref(unique_diags):

    rotated_diags = []

    for diag in unique_diags:
        new_diag = Diag(diag.y_start, diag.y_end, diag.x_start, diag.x_end, diag.orient, diag.offset)
        new_diag.true_reverse = diag.true_reverse

        rotated_diags.append(new_diag)

    return rotated_diags


def generate_dotplot_and_projection(ref_seq, alt_seq, snarl_ref_start, dotplot_output_prefix, dotplot_stride_size, options, ideal_csv=None):
    """

    for the ref seq and alt seq, generate dotplots and projections, which are further segmented into matrixes

    """

    # import psutil
    # start_memory = psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024

    dotplot_output_prefix = dotplot_output_prefix.replace(">", "+")
    dotplot_output_prefix = dotplot_output_prefix.replace("<", "-")
    dotplot_output_prefix = dotplot_output_prefix.replace("*", "none")

    dotplot_projection_bg_dict = {}     # # format {"ref2alt": projection_matrix for x2x, "alt2ref": projection_matrix for x2x}
    dotplot_projection_dict = {}    # # format {"ref2alt": projection_matrix, "alt2ref": projection_matrix}

    dotplot_id = -1

    x2x_dotplot = None
    x2y_dotplot = None

    unique_diags = None
    for dotplot_type in ["ref2alt", "alt2ref"]:
        dotplot_id += 1

        # print()
        # print(dotplot_type, len(ref_seq), len(alt_seq), dotplot_stride_size)

        dotplot_projection_bg_dict[dotplot_type] = []
        dotplot_projection_dict[dotplot_type] = []

        if dotplot_type is "ref2alt":
            x2x_dotplot = Dotplot(ref_seq, ref_seq, options.kmer_size, dotplot_output_prefix + ".{}.{}".format(dotplot_id, dotplot_type), stride_size=dotplot_stride_size, given_x_hash_table=None, given_y_hash_table=None)

            x2y_dotplot = Dotplot(ref_seq, alt_seq, options.kmer_size, dotplot_output_prefix + ".{}.{}".format(dotplot_id, dotplot_type), stride_size=dotplot_stride_size, given_x_hash_table=x2x_dotplot.get_seq_x_hash_table(), given_y_hash_table=None)
            # x2y_dotplot = Dotplot(ref_seq, alt_seq, options.kmer_size, dotplot_output_prefix + ".{}.{}".format(dotplot_id, dotplot_type), stride_size=dotplot_stride_size, given_x_hash_table=x2x_dotplot.get_seq_x_hash_table())

            x2x_dotplot_project_x, augment_coeff = x2x_dotplot.get_project_x(augment=True)
            x2x_dotplot_project_x_rev = x2x_dotplot.get_project_x_rev(baseline=augment_coeff)

            x2y_dotplot_project_x, _ = x2y_dotplot.get_project_x(augment=False)
            x2y_dotplot_project_x_rev = x2y_dotplot.get_project_x_rev(baseline=augment_coeff)

            unique_diags = find_and_denoise_diags(x2x_dotplot.matrix, x2y_dotplot.matrix, x2y_dotplot.matrix_rev, x2x_dotplot.stride_size)

            if len(unique_diags) == 0:
                dotplot_projection_dict['ref2alt'] = "Bad"
                dotplot_projection_dict['alt2ref'] = "Bad"

                return dotplot_projection_dict

            # # do two rounds of lienar extension
            unique_diags = linear_level_dotplot_extension(unique_diags, x2y_dotplot.matrix)
            unique_diags = linear_level_dotplot_extension(unique_diags, x2y_dotplot.matrix)

            unique_diags = kmer10_level_dotplot_extension(unique_diags, x2y_dotplot.seq_x, x2y_dotplot.seq_y, x2y_dotplot.matrix, dotplot_stride_size)
            unique_diags = base_level_dotplot_extension(unique_diags, x2y_dotplot.seq_x, x2y_dotplot.seq_y, x2y_dotplot.matrix, x2y_dotplot_project_x, x2y_dotplot_project_x_rev, options.kmer_size, x2x_dotplot.stride_size)

            x2y_matrix_len_x, x2y_matrix_len_y = x2y_dotplot.matrix.shape[1], x2y_dotplot.matrix.shape[0]

        else:

            x2x_dotplot = Dotplot(alt_seq, alt_seq, options.kmer_size, dotplot_output_prefix + ".{}.{}".format(dotplot_id, dotplot_type), stride_size=dotplot_stride_size, given_x_hash_table=None, given_y_hash_table=None)

            # x2y_dotplot = Dotplot(alt_seq, ref_seq, options.kmer_size, dotplot_output_prefix + ".{}.{}".format(dotplot_id, dotplot_type), stride_size=dotplot_stride_size, given_x_hash_table=x2x_dotplot.get_seq_x_hash_table(), given_y_hash_table=None)
            # x2y_dotplot = Dotplot(alt_seq, ref_seq, options.kmer_size, dotplot_output_prefix + ".{}.{}".format(dotplot_id, dotplot_type), stride_size=dotplot_stride_size, given_x_hash_table=x2x_dotplot.get_seq_x_hash_table())
            # x2y_dotplot.rotate_to_alt2ref()

            x2x_dotplot_project_x, augment_coeff = x2x_dotplot.get_project_x(augment=True)
            x2x_dotplot_project_x_rev = x2x_dotplot.get_project_x_rev(baseline=augment_coeff)

            x2y_dotplot_project_x, _ = x2y_dotplot.get_project_y(augment=False)
            x2y_dotplot_project_x_rev = x2y_dotplot.get_project_y_rev(baseline=augment_coeff)

            unique_diags = rotate_diags_to_alt2ref(unique_diags)

            x2y_matrix_len_x, x2y_matrix_len_y = x2y_dotplot.matrix.shape[0], x2y_dotplot.matrix.shape[1]

        # print()

        boost_diags(dotplot_type, unique_diags, x2y_matrix_len_x, x2y_matrix_len_y, x2y_dotplot_project_x, x2y_dotplot_project_x_rev, augment_coeff)

        projection_matrix, bad_flag = segment_projections_into_matrix(x2x_dotplot_project_x, x2y_dotplot_project_x, x2x_dotplot_project_x_rev, x2y_dotplot_project_x_rev, None, augment_coeff)
        dotplot_projection_dict[dotplot_type].extend(projection_matrix)

        # try:
        #     if dotplot_type == "ref2alt":
        #
        #         output_synthesized_dotplot(os.path.join(options.img_output_path, str(snarl_ref_start) + ".{}.{}.{}.{}".format(len(alt_seq), dotplot_id, dotplot_type, random.randint(1, 20))), x2y_dotplot.to_png())
        #         # output_synthesized_dotplot_project(os.path.join(options.img_output_path, str(snarl_ref_start) + ".{}.{}.{}".format(len(alt_seq), dotplot_id, dotplot_type)), x2x_dotplot_project_x, x2y_dotplot_project_x, x2x_dotplot_project_x_rev, x2y_dotplot_project_x_rev, )
        #
        #     # else:
        #     #     x2y_dotplot.rotate_to_alt2ref()
        #     #
        #     #     output_synthesized_dotplot(os.path.join(options.img_output_path, str(snarl_ref_start) + ".{}.{}.{}".format(len(alt_seq), dotplot_id, dotplot_type)), x2y_dotplot.to_png())
        #     #     output_synthesized_dotplot_project(os.path.join(options.img_output_path, str(snarl_ref_start) + ".{}.{}.{}".format(len(alt_seq), dotplot_id, dotplot_type)), x2x_dotplot_project_x, x2y_dotplot_project_x, x2x_dotplot_project_x_rev, x2y_dotplot_project_x_rev, )
        #
        # except OSError:     # # OSError: [Errno 36] File name too long
        #     pass

        if bad_flag is True:
            dotplot_projection_dict['ref2alt'] = "Bad"
            dotplot_projection_dict['alt2ref'] = "Bad"

            return dotplot_projection_dict

        # np.save(dotplot_output_prefix + ".{}.{}.npy".format(dotplot_id, dotplot_type), dotplot_projection)

    # end_memory = psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024
    #
    # print("generate memory cost GB ", end_memory - start_memory)

    return dotplot_projection_dict
