import os
import numpy as np
from src.pack_sv.op_vcf import output_vcf_header, generate_sample_level_vcf, output_vcf_records
from src.pack_model.model import index2label, label2index
from multiprocessing import Pool
import pysam
import logging
import sys
import traceback
from src.pack_graph.op_seq import calculate_seq_similarity, calculate_seq_repeat_ratio, find_continuous_val
import re


class SSV:

    def __init__(self, type, chr, start, end, alt_seq=None):

        self.id = "{}-{}-{}-{}".format(chr, start, end, type)

        self.type = type

        self.chr = chr
        self.start = start
        self.end = end

        self.alt_seq = alt_seq

        if self.type == "INS" or self.type == "insertion":
            self.length = len(self.alt_seq)
        else:
            # self.length = self.end - self.start
            self.length = self.end - self.start + 1

        self.source_chr = None
        self.source_start = None
        self.source_end = None

        self.insert_chr = None
        self.insert_start = None
        self.insert_end = None

        self.raw_dotplot = None
        self.raw_chr = None
        self.raw_start = None
        self.raw_end = None

        self.projection_matrix = None

        self.id = None

        self.within_inv = False
        self.within_repeat = False

    def set_id(self, id):
        self.id = id

    def refine_dup_events(self, inserted_chr, inserted_pos, tandem_thresh=None):
        """
        set the inserted chrom and pos for dup-like event
        """

        # if tandem_thresh is not None:
        #     if inserted_chr == self.chr and (abs(inserted_pos - self.start) < tandem_thresh or abs(inserted_pos - self.end) < tandem_thresh or (self.start <= inserted_pos <= self.end)):
        #         self.type = self.type.replace("DUP", "tDUP")
        #     else:
        #         self.type = self.type.replace("DUP", "dDUP")

        self.source_chr = self.chr
        self.source_start = self.start
        self.source_end = self.end

        self.insert_chr = inserted_chr
        self.insert_start = inserted_pos
        self.insert_end = inserted_pos + 1

    def set_raw_dotplot(self, raw_dotplot):
        self.raw_dotplot = raw_dotplot

    def set_raw_projection(self, projection_matrix):
        self.projection_matrix = projection_matrix

    def set_raw_cords(self, raw_chr, raw_start, raw_end):
        self.raw_chr = raw_chr
        self.raw_start = raw_start
        self.raw_end = raw_end

    def string_format(self):
        return "{},{}-{}-{},{}".format(self.type, self.chr, self.start, self.end, self.length)

    def to_string(self):

        if "DUP" in self.type or "duplication" in self.type:
            return "{}: {}-{}-{}-{}, {}".format(self.type, self.source_start, self.source_end, self.insert_start, self.insert_end, self.length)

        else:
            return "{}: {}-{}, {}".format(self.type, self.start, self.end, self.length)


class CSV:

    def __init__(self, ssv_components):
        self.ssv_components = ssv_components

        self.update_csv_info()

        self.score = "NA"

    def set_score(self, score):
        if score < 0 or score > 1.5:
            self.score = "LowQual"
        else:
            self.score = "PASS"

    def update_csv_info(self):

        if len(self.ssv_components) == 0:
            self.chr = -1
            self.start = -1
            self.end = -1
            self.length = -1
            self.type = "+".join([ssv.type for ssv in self.ssv_components])

        else:
            self.ssv_components = sorted(self.ssv_components, key=lambda x: (x.start, x.end))

            # # generate CSV info
            self.chr = self.ssv_components[0].chr
            self.start = self.ssv_components[0].start
            self.end = self.ssv_components[-1].end
            self.length = np.sum([ssv.length for ssv in self.ssv_components])

            self.type = "+".join([ssv.type for ssv in self.ssv_components])

    def add_ssv(self, new_ssv):

        self.ssv_components.append(new_ssv)

        # # update CSV info
        self.update_csv_info()

    def string_format(self):
        return "+".join(["{},{}-{}-{},{}".format(ssv.type, ssv.chr, ssv.start, ssv.end, ssv.length) for ssv in self.ssv_components])

    def to_string(self):
        return self.chr, self.start, self.end, self.length, self.score, [ssv.to_string() for ssv in self.ssv_components]


def handle_continuous_ssv_types(ssv_comps):
    ssv_comps = sorted(ssv_comps, key=lambda x: x.start)

    will_removed_ssv_comp = []
    for ssv_comp_index in range(len(ssv_comps) - 1, 0, -1):
        ssv_comp = ssv_comps[ssv_comp_index]

        previous_ssv_comp = ssv_comps[ssv_comp_index - 1]

        if ssv_comp.type != previous_ssv_comp.type:
            continue

        # # NOTE: there is no INS in ssv comps
        # if ssv_comp.type == "INS" and ssv_comp.type == previous_ssv_comp.type and abs(ssv_comp.start - (previous_ssv_comp.end + previous_ssv_comp.length)) < 500:
        #     # # update the insertion info into the previous one
        #     previous_ssv_comp.length += ssv_comp.length
        #     previous_ssv_comp.alt_seq += ssv_comp.alt_seq
        #     will_removed_ssv_comp.append(ssv_comp)

        # if ssv_comp.type in ["DEL", "INV", "DUP", "dDUP", "tDUP", ""] and ssv_comp.type == previous_ssv_comp.type and abs(ssv_comp.start - previous_ssv_comp.end) < 500:
        if ssv_comp.type in ["INV", "DUP", "invDUP"] and abs(ssv_comp.start - previous_ssv_comp.end) < 1000:
            # previous_ssv_comp.end += ssv_comp.length
            # previous_ssv_comp.length += ssv_comp.length

            previous_ssv_comp.end = ssv_comp.end
            previous_ssv_comp.length = previous_ssv_comp.end - previous_ssv_comp.start + 1

            will_removed_ssv_comp.append(ssv_comp)

        if ssv_comp.type in ["DEL"] and abs(ssv_comp.start - previous_ssv_comp.end) < 100:
            # previous_ssv_comp.end += ssv_comp.length
            # previous_ssv_comp.length += ssv_comp.length

            previous_ssv_comp.end = ssv_comp.end
            previous_ssv_comp.length = previous_ssv_comp.end - previous_ssv_comp.start + 1

            will_removed_ssv_comp.append(ssv_comp)

    new_ssv_comps = []
    for ssv_comp in ssv_comps:
        if ssv_comp not in will_removed_ssv_comp:
            new_ssv_comps.append(ssv_comp)

    return new_ssv_comps


def fetec_seq_from_fa_open(seq_fa_path, ref_name, ref_start=None, ref_end=None):
    seq_fa = open(seq_fa_path)

    seq_fa_lines = seq_fa.readlines()

    # # find the line_index of this reference name
    ref_name_index = None
    for line_index in range(len(seq_fa_lines)):
        line = seq_fa_lines[line_index].strip().replace(">", "")
        if ref_name == line:
            ref_name_index = line_index

            break

    # # return seq
    if ref_name_index is not None:
        if ref_start is None and ref_end is None:
            ref_name_seq = seq_fa_lines[ref_name_index + 1].strip()
        else:
            ref_name_seq = seq_fa_lines[ref_name_index + 1].strip()[ref_start: ref_end]

        seq_fa.close()
        return ref_name_seq
    else:
        seq_fa.close()
        return ""


def fetch_ssv_comp_seq(ssv_comp, dotplot_ref_start, seq_fa_path, ref_path_reformat, alt_path_reformat, max_seq_len=None, min_seq_len=None):
    if ssv_comp.type == "INS":
        seq_start = ssv_comp.raw_start - dotplot_ref_start
        if max_seq_len is not None:
            seq_end = min(seq_start + max_seq_len, ssv_comp.raw_start + ssv_comp.length - dotplot_ref_start)
        else:
            seq_end = ssv_comp.raw_start + ssv_comp.length - dotplot_ref_start

        seq_len = seq_end - seq_start
        if min_seq_len is not None and seq_len < min_seq_len:
            extend_len = 1 * seq_len
            seq = fetec_seq_from_fa_open(seq_fa_path, alt_path_reformat, max(0, seq_start - extend_len), seq_end + extend_len)

        else:
            seq = fetec_seq_from_fa_open(seq_fa_path, alt_path_reformat, seq_start, seq_end)

    else:
        seq_start = ssv_comp.start - dotplot_ref_start

        if max_seq_len is not None:
            seq_end = min(seq_start + max_seq_len, ssv_comp.end - dotplot_ref_start)
        else:
            seq_end = ssv_comp.end - dotplot_ref_start

        seq_len = seq_end - seq_start
        if min_seq_len is not None and seq_len < min_seq_len:
            extend_len = 1 * seq_len
            seq = fetec_seq_from_fa_open(seq_fa_path, "{}-{}".format(ref_path_reformat, alt_path_reformat), max(0, seq_start - extend_len), seq_end + extend_len)

        else:
            seq = fetec_seq_from_fa_open(seq_fa_path, "{}-{}".format(ref_path_reformat, alt_path_reformat), seq_start, seq_end)

    return seq


def merge_ref2alt_alt2ref_csv(ref_path, alt_path, ref2alt_comps, alt2ref_comps, seq_fa_path, dotplot_ref_start, dotplot_ref_end, dotplot_stride_size, options):
    """

    """

    ref_path_reformat = ref_path.replace(">", "+").replace("<", "-")
    alt_path_reformat = alt_path.replace(">", "+").replace("<", "-")

    ref2alt_comps = [comp for comp in handle_continuous_ssv_types(ref2alt_comps) if comp.length >= options.min_sv_size]

    ref2alt_comps = deal_with_dup(ref2alt_comps, dotplot_ref_start, seq_fa_path, ref_path_reformat, alt_path_reformat, options)

    ref2alt_comps = deal_with_ins_from_alt2ref(ref2alt_comps, alt2ref_comps, dotplot_ref_start, dotplot_ref_end, seq_fa_path, ref_path_reformat, alt_path_reformat, options)

    ref2alt_comps = find_dup_insert_pos(ref2alt_comps, dotplot_stride_size, options)

    # # deal with mutiple del and ins
    ref2alt_comps = deal_with_multi_ins_del(ref2alt_comps, seq_fa_path, ref_path_reformat, alt_path_reformat, dotplot_ref_start, options)

    # # generate a csv obj
    csv = CSV(handle_continuous_ssv_types(ref2alt_comps))
    csv.set_score(eval_csv_quality(csv, seq_fa_path, ref_path_reformat, alt_path_reformat))

    if len(csv.ssv_components) > options.max_sv_comps and is_scarred_inv(csv.type) is False:

        # ssv_stats = {}
        #
        # for ssv in csv.ssv_components:
        #     if ssv.type not in ssv_stats:
        #         ssv_stats[ssv.type] = 0
        #     ssv_stats[ssv.type] += ssv.length
        #
        # ssv_types = sorted(ssv_stats.keys(), reverse=True)
        #
        # csv.type = "hyperCPX_{}".format("+".join(ssv_types))

        ssv_type_list = [ssv.type for ssv in csv.ssv_components]

        if "INV" in ssv_type_list:
            csv.type = "hyperCPX_INV"
        elif "invDUP" in ssv_type_list:
            csv.type = "hyperCPX_invDUP"
        elif "DUP" in ssv_type_list:
            csv.type = "hyperCPX_DUP"
        elif "DEL" in ssv_type_list and "INS" in ssv_type_list:
            csv.type = "hyperCPX_INS+DEL"
        elif "DEL" in ssv_type_list:
            csv.type = "hyperCPX_DEL"
        elif "INS" in ssv_type_list:
            csv.type = "hyperCPX_INS"
        else:
            csv.type = "hyperCPX"

        csv.ssv_components = [SSV(csv.type, csv.chr, csv.start, csv.end)]

    return csv


def is_scarred_inv(sv_type):

    # # using the much ideal way for keeping scarred inv
    scarred_inv = True
    sv_type_split = sv_type.split("+")
    if not len(sv_type_split) >= 3:
        return False

    if sv_type_split[0] == "INV":
        for i in range(len(sv_type_split)):

            if i % 2 == 0 and sv_type_split[i] != "INV":
                scarred_inv = False
                break

            if i % 2 != 0 and sv_type_split[i] not in ["INS", "DEL"]:
                scarred_inv = False
                break
    else:
        for i in range(len(sv_type_split)):

            if i % 2 != 0 and sv_type_split[i] != "INV":
                scarred_inv = False
                break

            if i % 2 == 0 and sv_type_split[i] not in ["INS", "DEL"]:
                scarred_inv = False
                break

    return scarred_inv


def deal_with_ins_from_alt2ref(ref2alt_comps, alt2ref_comps, dotplot_ref_start, dotplot_ref_end, seq_fa_path, ref_path_reformat, alt_path_reformat, options):
    # #
    dup_exist = False
    for ref2alt_ssv in ref2alt_comps:
        if ref2alt_ssv.type in ["DUP"]:
            dup_exist = True
            break

    # #
    ins_within_inv = []

    for i in range(len(alt2ref_comps)):

        # # ins is del in alt2ref
        if alt2ref_comps[i].type == "DEL":
            if i == 0 or i == len(alt2ref_comps) - 1:
                continue
            try:
                if alt2ref_comps[i - 1].type in ["INV", "invDUP"] and alt2ref_comps[i + 1].type in ["INV", "invDUP"]:
                    # if alt2ref_comps[i - 1].type in ["INV"] and alt2ref_comps[i + 1].type in ["INV"]:
                    ins_within_inv.append(alt2ref_comps[i])
            except IndexError:
                pass

    # #
    alt2ref_comps = [comp for comp in handle_continuous_ssv_types(alt2ref_comps) if comp.type == "DEL" and comp.length >= options.min_sv_size]

    # # traverse those remained ssv in alt2ref, for capturing the INSs that are not reflected in ref2alt
    for alt2ref_ssv in alt2ref_comps:

        new_start = alt2ref_ssv.start

        raw_start = alt2ref_ssv.start
        for tmp_ssv in alt2ref_comps:
            if tmp_ssv.type == "DEL" and raw_start >= tmp_ssv.end:
                new_start -= tmp_ssv.length

            # if tmp_ssv.type == "INV":
            #     new_start -= tmp_ssv.length
        raw_start = new_start
        for tmp_ssv in ref2alt_comps:
            if tmp_ssv.type == "DEL" and raw_start >= tmp_ssv.end:
                new_start += tmp_ssv.length

            # if tmp_ssv.type == "DUP" or tmp_ssv.type == "invDUP": # # NOTE: the DUP in ref2alt is DEL in alt2ref, so we do not need to del the DUP length since we have del all DEL length from alt2ref
            #     new_start -= tmp_ssv.length

        # # add into ref2alt_comps
        if options.min_sv_size <= alt2ref_ssv.length <= options.max_sv_size:

            new_ssv = SSV("INS", alt2ref_ssv.chr, new_start, new_start, alt_seq="N" * alt2ref_ssv.length)

            # print('before', new_ssv.to_string())

            # # check if this insertion is covered by a inversion
            if alt2ref_ssv in ins_within_inv:
                # # reverse the cord of this insertion
                new_ssv = SSV("INS", alt2ref_ssv.chr, dotplot_ref_start + (dotplot_ref_end - new_start), dotplot_ref_start + (dotplot_ref_end - new_start), alt_seq="N" * alt2ref_ssv.length)
                # print('after', new_ssv.to_string())

                # # split the whole inv that is cutted by the ins
                for target_ssv in ref2alt_comps:
                    if target_ssv.type == "INV" and target_ssv.start + 1 <= new_ssv.start <= new_ssv.end <= target_ssv.end - 1:
                        # new_ssv = SSV("INS", alt2ref_ssv.chr, dotplot_ref_start + (dotplot_ref_end - new_start), dotplot_ref_start + (dotplot_ref_end - new_start), alt_seq="N" * alt2ref_ssv.length)

                        split_inv_left = SSV("INV", target_ssv.chr, target_ssv.start, new_ssv.start - 1)
                        split_inv_right = SSV("INV", target_ssv.chr, new_ssv.end + 1, target_ssv.end)
                        if options.min_sv_size <= split_inv_left.length <= options.max_sv_size:
                            ref2alt_comps.append(split_inv_left)
                        if options.min_sv_size <= split_inv_right.length <= options.max_sv_size:
                            ref2alt_comps.append(split_inv_right)
                        ref2alt_comps.remove(target_ssv)

                        break

                new_ssv.set_raw_cords(alt2ref_ssv.chr, alt2ref_ssv.start, alt2ref_ssv.end)
                new_ssv.within_inv = True

                if dup_exist:
                    new_ssv_seq = fetch_ssv_comp_seq(new_ssv, dotplot_ref_start, seq_fa_path, ref_path_reformat, alt_path_reformat, max_seq_len=5000)
                    repeat_ratio = calculate_seq_repeat_ratio(new_ssv_seq)
                    if repeat_ratio > 1.2:
                        new_ssv.within_repeat = True
                    # print(repeat_ratio, len(new_ssv_seq), new_ssv.to_string())

                ref2alt_comps.append(new_ssv)

            else:
                covered_flag = False
                for target_ssv in ref2alt_comps:
                    if target_ssv.type == "INV" and target_ssv.start + 1 <= new_ssv.start <= new_ssv.end <= target_ssv.end - 1:
                        covered_flag = True
                        break

                if not covered_flag:
                    new_ssv.set_raw_cords(alt2ref_ssv.chr, alt2ref_ssv.start, alt2ref_ssv.end)

                    if dup_exist:
                        new_ssv_seq = fetch_ssv_comp_seq(new_ssv, dotplot_ref_start, seq_fa_path, ref_path_reformat, alt_path_reformat, max_seq_len=5000)
                        repeat_ratio = calculate_seq_repeat_ratio(new_ssv_seq)
                        if repeat_ratio > 1.2:
                            new_ssv.within_repeat = True

                        # print(repeat_ratio, len(new_ssv_seq), new_ssv.to_string())

                    ref2alt_comps.append(new_ssv)

    return ref2alt_comps


def find_dup_insert_pos(ref2alt_comps, dotplot_stride_size, options):
    dup_comp_no_insert = []
    consumed_ins = []
    for dup_comp in ref2alt_comps:

        if "DUP" not in dup_comp.type:
            continue

        # # STEP: determine the inserted pos of this dup by finding the most likely insertion in alt2ref
        best_match_ins_index = -1
        best_match_ins_score = -1

        for index in range(len(ref2alt_comps)):

            if ref2alt_comps[index].type != "INS":
                continue

            # # this alt2ref_ssv has been consumed
            if ref2alt_comps[index].length < max(options.min_sv_size, 2 * dotplot_stride_size):
                continue

            if ref2alt_comps[index].within_repeat is True:
                continue

            # #
            # ins_comp_seq = seq_fa.fetch(alt_path_reformat, alt2ref_ssv.start - dotplot_ref_start, alt2ref_ssv.start + alt2ref_ssv.length - dotplot_ref_start)
            # seq_similarity = calculate_seq_similarity(dup_comp_seq, ins_comp_seq)
            # print(seq_similarity, dup_comp.length, alt2ref_ssv.length)

            # #
            # dup_comp_seq = fetch_ssv_comp_seq(dup_comp, dotplot_ref_start, seq_fa_path, ref_path_reformat, alt_path_reformat, allowed_seq_len=None)
            # new_ssv_seq = fetch_ssv_comp_seq(ref2alt_comps[index], dotplot_ref_start, seq_fa_path, ref_path_reformat, alt_path_reformat, allowed_seq_len=None)
            #
            # if dup_comp_seq[: -dotplot_stride_size] in new_ssv_seq or dup_comp_seq[dotplot_stride_size: ] in new_ssv_seq:
            #     cur_match_score = 1 - abs(1 - dup_comp.length / ref2alt_comps[index].length)
            #
            #     if cur_match_score > best_match_ins_score:
            #         best_match_ins_index = index
            #         best_match_ins_score = cur_match_score

            # #
            ins_comp = ref2alt_comps[index]
            length_thresh = max(options.min_sv_size, 2 * dotplot_stride_size)
            if dup_comp.type == "DUP":
                # # exact length similarity
                if abs(dup_comp.length - ins_comp.length) < length_thresh:
                    cur_match_score = 1 - abs(1 - dup_comp.length / ins_comp.length)

                    if cur_match_score > best_match_ins_score:
                        best_match_ins_index = index
                        best_match_ins_score = cur_match_score

                else:
                    # # is a tandem dup
                    if abs(dup_comp.start - ins_comp.start) < length_thresh or abs(dup_comp.end - ins_comp.start) < length_thresh:
                        cur_match_score = 1 - abs(1 - dup_comp.length / ins_comp.length)

                        if cur_match_score > best_match_ins_score:
                            best_match_ins_index = index
                            best_match_ins_score = cur_match_score

            else:
                #
                cur_match_score = 1 - abs(1 - dup_comp.length / ins_comp.length)

                if cur_match_score > best_match_ins_score:
                    best_match_ins_index = index
                    best_match_ins_score = cur_match_score

            # #
            # cur_match_score = 1 - abs(1 - dup_comp.length / ins_comp.length)
            #
            # if cur_match_score > best_match_ins_score:
            #     best_match_ins_index = index
            #     best_match_ins_score = cur_match_score

        if best_match_ins_index != -1:
            dup_comp.refine_dup_events(ref2alt_comps[best_match_ins_index].chr, ref2alt_comps[best_match_ins_index].start, tandem_thresh=None)
            # identical_comp.refine_dup_events(alt2ref_comps[best_match_ins_index].chr, alt2ref_comps[best_match_ins_index].start - identical_comp.length, tandem_thresh=tandem_thresh)

            ref2alt_comps[best_match_ins_index].length -= dup_comp.length

            if ref2alt_comps[best_match_ins_index].length < options.min_sv_size:
                consumed_ins.append(ref2alt_comps[best_match_ins_index])
        else:
            dup_comp_no_insert.append(dup_comp)

    for dup_comp in dup_comp_no_insert:
        ref2alt_comps.remove(dup_comp)

    for ins_comp in consumed_ins:
        ref2alt_comps.remove(ins_comp)

    return ref2alt_comps


def deal_with_dup(ref2alt_comps, dotplot_ref_start, seq_fa_path, ref_path_reformat, alt_path_reformat, options):
    # # check invdup between inversions
    for comp_index in range(len(ref2alt_comps)):
        if ref2alt_comps[comp_index].type == "invDUP":
            if comp_index == 0 or comp_index == len(ref2alt_comps) - 1:
                continue
            try:
                # if ref2alt_comps[comp_index - 1].type == "INV" and ref2alt_comps[comp_index + 1].type == "INV":
                if ((ref2alt_comps[comp_index - 1].type == "INV" and ref2alt_comps[comp_index + 1].type == "INV") or
                        (ref2alt_comps[comp_index - 1].type == "invDUP" and ref2alt_comps[comp_index + 1].type == "INV") or
                        (ref2alt_comps[comp_index - 1].type == "INV" and ref2alt_comps[comp_index + 1].type == "invDUP")):
                    ref2alt_comps[comp_index].type = "INV"

            except IndexError:
                pass

    # # STEP: collect and merge continuous variants
    if options.dup_to_ins:
        ref2alt_comps = [comp for comp in handle_continuous_ssv_types(ref2alt_comps) if comp.type not in ["DUP", "invDUP"] and comp.length >= options.min_sv_size]
    else:
        ref2alt_comps = [comp for comp in handle_continuous_ssv_types(ref2alt_comps) if comp.length >= options.min_sv_size]

        # STEP: for dup-like SSVs, we need the alt2ref results to refine its inserted position

        dup_cnt = len([ssv.type for ssv in ref2alt_comps if 'DUP' in ssv.type])

        # # massive dups
        if dup_cnt >= 3:
            # print("massive dups")
            dup_comp_at_mass = []
            for dup_comp in ref2alt_comps:

                if dup_comp.type in ["DUP", "invDUP"]:
                    dup_comp_at_mass.append(dup_comp)

            for dup_comp in dup_comp_at_mass:
                ref2alt_comps.remove(dup_comp)
        else:

            dup_comp_at_repeats = []
            for dup_comp in ref2alt_comps:

                if dup_comp.type != "DUP":
                    continue

                # # STEP: determine whether the dup is located at repeat regions
                dup_comp_seq = fetch_ssv_comp_seq(dup_comp, dotplot_ref_start, seq_fa_path, ref_path_reformat, alt_path_reformat, max_seq_len=None, min_seq_len=200)

                repeat_ratio = calculate_seq_repeat_ratio(dup_comp_seq)

                # print(repeat_ratio, len(dup_comp_seq), dup_comp.to_string())
                # print(dup_comp_seq)

                if repeat_ratio > 1.2:  # # the dup locates at a repeat region
                    dup_comp_at_repeats.append(dup_comp)
                    continue

            for dup_comp in dup_comp_at_repeats:
                ref2alt_comps.remove(dup_comp)

    return ref2alt_comps


def deal_with_multi_ins_del(ref2alt_comps, seq_fa_path, ref_path_reformat, alt_path_reformat, dotplot_ref_start, options):
    ref2alt_comps = sorted(ref2alt_comps, key=lambda x: (x.start, x.end))

    # for i in range(len(ref2alt_comps)):
    #
    #     # # ins is del in alt2ref
    #     if ref2alt_comps[i].type == "DEL":
    #         if i == 0 or i == len(ref2alt_comps) - 1:
    #             continue
    #         try:
    #             if ref2alt_comps[i - 1].type in ["INV", "invDUP"] and ref2alt_comps[i + 1].type in ["INV", "invDUP"]:
    #             # if alt2ref_comps[i - 1].type in ["INV"] and alt2ref_comps[i + 1].type in ["INV"]:
    #                 ref2alt_comps[i].within_inv = True
    #         except IndexError:
    #             pass

    ref2alt_comps_types = [ssv.type for ssv in ref2alt_comps]

    if ("DEL" in ref2alt_comps_types and "INS" in ref2alt_comps_types) or ref2alt_comps_types.count("INS") > 1 or ref2alt_comps_types.count("DEL") > 1:
        allowed_seq_len = 5000

        # # first round, deal with DEL INS from unlinear mapping due to snps
        identical_comps = []
        ins_comps = sorted([ssv for ssv in ref2alt_comps if ssv.type in ["INS"]], key=lambda x: x.length, reverse=True)
        del_comps = sorted([ssv for ssv in ref2alt_comps if ssv.type in ["DEL"]], key=lambda x: x.length, reverse=True)

        for ins_comp in ins_comps:
            # if ins_comp.within_inv is True:
            #     continue

            if ins_comp in identical_comps:
                continue

            ins_seq = fetch_ssv_comp_seq(ins_comp, dotplot_ref_start, seq_fa_path, ref_path_reformat, alt_path_reformat, max_seq_len=allowed_seq_len)

            # if ins_seq is None:
            #     continue

            for del_comp in del_comps:

                if del_comp in identical_comps:
                    continue

                if abs(ins_comp.length - del_comp.length) <= options.min_sv_size:

                    del_seq = fetch_ssv_comp_seq(del_comp, dotplot_ref_start, seq_fa_path, ref_path_reformat, alt_path_reformat, max_seq_len=allowed_seq_len)

                    # if del_seq in None:
                    #     continue

                    seq_similarity = calculate_seq_similarity(ins_seq, del_seq)

                    if seq_similarity >= 0.5:
                        identical_comps.append(del_comp)
                        identical_comps.append(ins_comp)

                        break

        for identical_comp in identical_comps:
            ref2alt_comps.remove(identical_comp)

        # # second round, choose the largest comp as a anchor, compare the others with it
        del_ins_comps = sorted([ssv for ssv in ref2alt_comps if ssv.type in ["DEL", "INS"]], key=lambda x: x.length, reverse=True)

        if len(del_ins_comps) != 0:
            largest_del_ins_comp = del_ins_comps[0]

            identical_comps = [largest_del_ins_comp]

            # # compare with the largest del or ins comp, and calculate the sequence identity
            # # fetch the sequence for largest del_ins_comp
            largest_seq = fetch_ssv_comp_seq(largest_del_ins_comp, dotplot_ref_start, seq_fa_path, ref_path_reformat, alt_path_reformat, max_seq_len=allowed_seq_len)

            for index in range(1, len(del_ins_comps)):

                cur_del_ins_comp = del_ins_comps[index]

                # if cur_del_ins_comp.within_inv is True:
                #     continue

                cur_seq = fetch_ssv_comp_seq(cur_del_ins_comp, dotplot_ref_start, seq_fa_path, ref_path_reformat, alt_path_reformat, max_seq_len=allowed_seq_len)

                seq_similarity = calculate_seq_similarity(largest_seq, cur_seq)

                if seq_similarity >= 0.5:
                    identical_comps.append(cur_del_ins_comp)

            # # find the final type and length
            final_del_ins_chr = identical_comps[0].chr
            final_del_ins_pos = identical_comps[0].start  # # use the leftmost cord
            final_del_ins_length = 0

            for identical_comp in identical_comps:
                if identical_comp.type == "DEL":
                    final_del_ins_length -= identical_comp.length
                else:
                    final_del_ins_length += identical_comp.length

                # if identical_comp.start < final_del_ins_pos:
                #     final_del_ins_pos = identical_comp.start

                ref2alt_comps.remove(identical_comp)

            if final_del_ins_length >= options.min_sv_size:  # # a insertion

                ref2alt_comps.append(SSV("INS", final_del_ins_chr, final_del_ins_pos, final_del_ins_pos, alt_seq="N" * final_del_ins_length))

            elif final_del_ins_length <= -options.min_sv_size:  # # a deletion

                final_del_ins_length = abs(final_del_ins_length)

                ref2alt_comps.append(SSV("DEL", final_del_ins_chr, final_del_ins_pos, final_del_ins_pos + final_del_ins_length))

            else:
                pass

    return ref2alt_comps


def eval_csv_quality(csv, seq_fa_path, ref_path, alt_path):
    e = 0.0001
    try:

        ref_seq = fetec_seq_from_fa_open(seq_fa_path, "{}-{}".format(ref_path, alt_path))
        alt_seq = fetec_seq_from_fa_open(seq_fa_path, alt_path)

    except:
        return None

    seq_effective_length = len(alt_seq) - len(ref_seq)

    csv_effective_length = 0
    for ssv in csv.ssv_components:
        if ssv.type in ["INS", "DUP", "invDUP"]:
            csv_effective_length += ssv.length
        elif ssv.type in ["DEL"]:
            csv_effective_length -= ssv.length
        else:
            pass

    if -50 <= abs(csv_effective_length) - abs(seq_effective_length) <= 50:
        return 1
    else:
        return round((csv_effective_length + e) / (seq_effective_length + e), 3)


def manual_interpret_inv(segment, raw_prediction):
    # # we calculate a manual type
    x2x_value = float(segment[2])
    x2y_value = float(segment[3])
    x2y_rev_value = float(segment[4])
    allowed_thresh = x2x_value * 0.1

    if x2x_value - allowed_thresh <= abs(x2y_rev_value) <= x2x_value + allowed_thresh:
        manual_prediction = label2index["INV"]
        return manual_prediction

    else:
        return raw_prediction


def prediction_to_ssv_comps(matrix, prediction, dotplot_type, sub_dotplot_ref_chr, sub_dotplot_ref_start, dotplot_stride_size, options):
    ssv_components = []
    for pred_label_index, include_indexes in find_continuous_val(prediction):
        # print(pred_label_index, include_indexes, np.shape(matrix))
        # # 0: None, 1: REF, they are not SVs, so we skip
        if pred_label_index == 0 or pred_label_index == 1:
            continue

        ssv_type = index2label[pred_label_index]

        ssv_start = int(matrix[include_indexes[0]][0])  # # first seg's start pos
        ssv_end = int(matrix[include_indexes[-1]][1]) - 1  # # last seg's end pos. why -1: the matrix is 'Open left and close right'
        # ssv_end = int(matrix[include_indexes[-1]][1])

        ssv_obj = SSV(ssv_type, sub_dotplot_ref_chr, sub_dotplot_ref_start + (ssv_start * dotplot_stride_size), sub_dotplot_ref_start + (ssv_end * dotplot_stride_size), )

        if options.max_sv_size >= ssv_obj.length >= max(options.min_sv_size, 2 * dotplot_stride_size):
            # if options.max_sv_size >= ssv_obj.length >=options.min_sv_size:
            #
            # print(dotplot_type, ssv_type, ssv_start, ssv_end, sub_dotplot_ref_start + (ssv_start * dotplot_stride_size), sub_dotplot_ref_start + (ssv_end * dotplot_stride_size), (ssv_end - ssv_start) * dotplot_stride_size)
            ssv_obj.set_raw_projection([matrix[i] for i in include_indexes])
            ssv_components.append(ssv_obj)

    return ssv_components


def interpret_one_bin(bin_index, snarls, options):
    try:
        logging.info("Interpreting parallel bin index {}".format(bin_index))

        predict_npz_path = os.path.join(options.output_path, "tmp.matrix.{}.res.npz".format(bin_index))

        # raw_matrix_dict = np.load(matrix_npz_path)
        predict_res_dict = np.load(predict_npz_path)

        interpret_res_dict = {}

        for matrix_id in predict_res_dict:
            matrix_id_split = matrix_id.split(".")

            top_snarl_id, sub_ref_start, sub_ref_end, sub_dotplot_ref_chr, sub_dotplot_ref_start, sub_dotplot_ref_end, sub_ref_path, sub_alt_path, dotplot_type, dotplot_stride_size \
                = matrix_id_split[0], int(matrix_id_split[1]), int(matrix_id_split[2]), matrix_id_split[3], int(matrix_id_split[4]), int(matrix_id_split[5]), matrix_id_split[6], matrix_id_split[7], matrix_id_split[8], int(matrix_id_split[9])

            sub_snarl_id = "{}-{}".format(sub_ref_path, sub_alt_path)

            matrix_snarl_npz_path = os.path.join(options.npz_output_path, "{}.npz".format(top_snarl_id.replace(">", "+").replace("<", "-")))
            matrix_snarl_seq_path = os.path.join(options.npz_output_path, "{}.seq.fa".format(top_snarl_id.replace(">", "+").replace("<", "-")))

            raw_matrix_dict = np.load(matrix_snarl_npz_path)

            if matrix_id not in raw_matrix_dict:
                continue

            if top_snarl_id not in interpret_res_dict:
                interpret_res_dict[top_snarl_id] = {}

            if sub_snarl_id not in interpret_res_dict[top_snarl_id]:
                interpret_res_dict[top_snarl_id][sub_snarl_id] = {"ref2alt": None, "alt2ref": None}

            matrix = raw_matrix_dict[matrix_id]
            prediction = predict_res_dict[matrix_id]

            # # fix the wrong prediction of 1 (REF), assign the wrong prediction to its neighbor right prediction
            for index in range(len(prediction)):
                if index == 0 or index == len(prediction) - 1:
                    continue

                if prediction[index] == 1:
                    cur_segment = matrix[index]
                    allowed_thresh = float(cur_segment[2]) * 0.1

                    if prediction[index - 1] != 1:
                        previous_segment = matrix[index - 1]

                        if abs(float(cur_segment[3]) - float(previous_segment[3])) <= allowed_thresh and abs(float(cur_segment[4]) - float(previous_segment[4])) <= allowed_thresh:
                            prediction[index] = prediction[index - 1]

                    if prediction[index + 1] != 1:
                        latter_segment = matrix[index + 1]

                        if abs(float(cur_segment[3]) - float(latter_segment[3])) <= allowed_thresh and abs(float(cur_segment[4]) - float(latter_segment[4])) <= allowed_thresh:
                            prediction[index] = prediction[index + 1]

            ssv_components = prediction_to_ssv_comps(matrix, prediction, dotplot_type, sub_dotplot_ref_chr, sub_dotplot_ref_start, dotplot_stride_size, options)

            if len(ssv_components) == 0:
                # # failed prediction, then apply a manual prediction for inversion since it has a unique wave
                for i in range(len(prediction)):
                    prediction[i] = manual_interpret_inv(matrix[i], prediction[i])

                ssv_components = prediction_to_ssv_comps(matrix, prediction, dotplot_type, sub_dotplot_ref_chr, sub_dotplot_ref_start, dotplot_stride_size, options)

            interpret_res_dict[top_snarl_id][sub_snarl_id][dotplot_type] = ssv_components

            # # STEP: merge the ref2alt and alt2ref CSV if both are already called
            ref2alt_comps = interpret_res_dict[top_snarl_id][sub_snarl_id]["ref2alt"]
            alt2ref_comps = interpret_res_dict[top_snarl_id][sub_snarl_id]["alt2ref"]

            if ref2alt_comps is not None and alt2ref_comps is not None:
                """
                   if not both ref2alt and alt2ref are precessed, 
                   the interpret_res_dict[top_snarl_id][alt_path] will be like {"ref2alt": None, "alt2ref": [ssv_comps]}, is a dict
                """

                # # if ref2alt failed the prediction, then we retrieve ssvs from alt2ref
                if len(ref2alt_comps) == 0 and len(alt2ref_comps) != 0:

                    for comp in alt2ref_comps:
                        # # why inv, inv is shared between the ref2alt and alt2ref since inv generates same dotplot and projections between the ref2alt and alt2ref
                        if comp.type == "INV":
                            retrieved_ssv = SSV("INV", comp.chr, sub_dotplot_ref_start + (sub_dotplot_ref_end - comp.end), sub_dotplot_ref_start + (sub_dotplot_ref_end - comp.start))
                            ref2alt_comps = [retrieved_ssv]

                # # update the store dict
                interpret_res_dict[top_snarl_id][sub_snarl_id] = merge_ref2alt_alt2ref_csv(sub_ref_path, sub_alt_path, ref2alt_comps, alt2ref_comps, matrix_snarl_seq_path, sub_dotplot_ref_start, sub_dotplot_ref_end, dotplot_stride_size, options)

                # print()
                # print("Final CSV:", interpret_res_dict[top_snarl_id][sub_snarl_id].to_string())
                # print("Final CSV:", "{}, {}-{}-{}-{}-{}".format(top_snarl_id, interpret_res_dict[top_snarl_id][sub_snarl_id].chr, interpret_res_dict[top_snarl_id][sub_snarl_id].start, interpret_res_dict[top_snarl_id][sub_snarl_id].end, interpret_res_dict[top_snarl_id][sub_snarl_id].type, interpret_res_dict[top_snarl_id][sub_snarl_id].length))

                """
                    after here, the interpret_res_dict[top_snarl_id][alt_path] will be like [ssv_comps]
                """

        logging.info("Recording parallel bin index {}".format(bin_index))
        vcf_records, vcf_records_split = generate_vcf_records(interpret_res_dict, snarls, options)

        hap_level_vcf_path = os.path.join(options.output_path, "tmp.matrix.{}.res.vcf".format(bin_index))
        hap_level_vcf_fout = open(hap_level_vcf_path, "w")

        output_vcf_header(hap_level_vcf_fout, options)
        output_vcf_records(hap_level_vcf_fout, vcf_records)

        hap_level_vcf_fout.close()

        # # output to split/single vcf file
        hap_level_split_vcf_path = os.path.join(options.output_path, "tmp.matrix.{}.res.split.vcf".format(bin_index))
        hap_level_split_vcf_fout = open(hap_level_split_vcf_path, "w")

        output_vcf_header(hap_level_split_vcf_fout, options)
        output_vcf_records(hap_level_split_vcf_fout, vcf_records_split)

        hap_level_split_vcf_fout.close()

        return ["finish", hap_level_vcf_path, hap_level_split_vcf_path]

    except:
        error_type, error_value, error_trace = sys.exc_info()
        error_log = "Error log: " + str(error_type) + ": " + str(error_value) + '. ' + 'Locate At: ' + str(traceback.extract_tb(error_trace))
        error_info = [bin_index, error_log]

        return error_info


def output_at_top_and_sub_snarl_level(vcf_records, vcf_records_split, top_snarl_id, snarl, snarl_ref_start, snarl_ref_end, snarl_path_csv_dict, snarl_subsnarl_csv_dict, path_subsnarl_dict, options):
    # # STEP: find the ref asm path and these asms that share the same path with ref
    try:
        ref_asm_path = snarl.ref_asm_path
    except AttributeError:
        ref_asm_path = "."

    try:
        ref_path_include_asm_names = snarl.path_asm_dict[ref_asm_path]
    except KeyError:
        ref_path_include_asm_names = []

    # # STEP
    if snarl.reversed_mapping:
        top_snarl_gt_list = ["0:.:.:." for i in range(len(options.alt_asm_names))]  # # hap_gt_format = "GT:Type:Length:BKP"
    else:
        top_snarl_gt_list = [".:.:.:." for i in range(len(options.alt_asm_names))]  # # hap_gt_format = "GT:Type:Length:BKP"

    # # for asms that share the same path with ref_asm, we set them as 0
    for asm_name in ref_path_include_asm_names:
        if asm_name == options.ref_asm_name:
            continue
        top_snarl_gt_list[options.alt_asm_names.index(asm_name)] = "0:.:.:."

    snarl_include_alt_paths = []

    alt_path_cnt = 0
    for alt_path in snarl_path_csv_dict:
        if alt_path == ref_asm_path:
            continue

        alt_path_csv = snarl_path_csv_dict[alt_path]
        alt_path_include_asm_names = snarl.path_asm_dict[alt_path]

        # # this path is predicted, but there is no SV, we set the GT as 0
        if len(alt_path_csv.ssv_components) == 0:
            for asm_name in alt_path_include_asm_names:
                top_snarl_gt_list[options.alt_asm_names.index(asm_name)] = "0:.:.:."
            continue

        alt_path_cnt += 1

        # # update the gt information
        alt_path_csv_type, alt_path_csv_length = alt_path_csv.type, alt_path_csv.length

        # alt_path_csv_type, alt_path_csv_length, alt_path_csv_qual = alt_path_csv.type, alt_path_csv.length, alt_path_csv.score
        #
        # if snarl.reversed_mapping:
        #     if alt_path_csv.score != "LowQual":
        #         alt_path_csv.score = "MediumQual"

        alt_path_csv_bkps = ["{}_{}_{}_{}_{}_{}_{}_{}".format(ssv.type, ssv.length, ssv.source_chr, ssv.source_start, ssv.source_end, ssv.insert_chr, ssv.insert_start, ssv.insert_end) if ssv.type in ["DUP", "invDUP"] else "{}_{}_{}_{}_{}".format(ssv.type, ssv.length, ssv.chr, ssv.start, ssv.end) for ssv in alt_path_csv.ssv_components]

        for asm_name in alt_path_include_asm_names:
            top_snarl_gt_list[options.alt_asm_names.index(asm_name)] = "{}:{}:{}:{}".format(alt_path_cnt, alt_path_csv_type, alt_path_csv_length, ",".join(alt_path_csv_bkps))

        snarl_include_alt_paths.append(alt_path)

    top_snarl_AC, top_snarl_AF, top_snarl_AN, top_snarl_NS = calculate_ac_af_an_ns(top_snarl_gt_list, options)

    if snarl.reversed_mapping:
        top_snarl_qual = "MediumQual"
    else:
        top_snarl_qual = "PASS"

    vcf_records.append([str(snarl.ref_chrom), str(snarl_ref_start), "{}".format(top_snarl_id), ref_asm_path, ",".join(snarl_include_alt_paths), ".", top_snarl_qual, "END={};AC={};AF={};AN={};NS={};LV=0".format(snarl_ref_end, top_snarl_AC, top_snarl_AF, top_snarl_AN, top_snarl_NS), "GT:TYPE:LENGTH:BKPS", "\t".join(top_snarl_gt_list)])

    # # output at subsnarl level

    for sub_snarl_id in snarl_subsnarl_csv_dict:

        sub_snarl_csv = snarl_subsnarl_csv_dict[sub_snarl_id]

        if not len(sub_snarl_csv.ssv_components) > 0:
            continue

        sub_snarl_id_split = sub_snarl_id.split("-")
        sub_snarl_ref_path = sub_snarl_id_split[0]
        sub_snarl_alt_path = sub_snarl_id_split[1]

        # # re-generate the path and snarl info
        sub_snarl_ref_path_include_nodes = re.findall(r'[><]([a-zA-Z0-9]+)', sub_snarl_ref_path)
        sub_snarl_ref_path_include_nodes_orients = re.findall(r'[><]', sub_snarl_ref_path)
        sub_snarl_alt_path_include_nodes = re.findall(r'[><]([a-zA-Z0-9]+)', sub_snarl_alt_path)
        sub_snarl_alt_path_include_nodes_orients = re.findall(r'[><]', sub_snarl_alt_path)

        try:
            sub_snarl_id_only_start_end = "{}{}{}{}".format(sub_snarl_ref_path_include_nodes_orients[0], sub_snarl_ref_path_include_nodes[0], sub_snarl_ref_path_include_nodes_orients[-1], sub_snarl_ref_path_include_nodes[-1])

        # # for reversed mapping snarls, there is no ref path
        except IndexError:
            sub_snarl_id_only_start_end = top_snarl_id

        sub_snarl_ref_path_rm_start_end = "".join("{}{}".format(sub_snarl_ref_path_include_nodes_orients[i], sub_snarl_ref_path_include_nodes[i]) for i in range(1, len(sub_snarl_ref_path_include_nodes) - 1))
        sub_snarl_alt_path_rm_start_end = "".join("{}{}".format(sub_snarl_alt_path_include_nodes_orients[i], sub_snarl_alt_path_include_nodes[i]) for i in range(1, len(sub_snarl_alt_path_include_nodes) - 1))

        if sub_snarl_ref_path_rm_start_end == "":
            sub_snarl_ref_path_rm_start_end = "*"
        if sub_snarl_alt_path_rm_start_end == "":
            sub_snarl_alt_path_rm_start_end = "*"

        # # STEP: generate gt
        sub_snarl_gt_list = []

        for gt_info in top_snarl_gt_list:
            gt_flag = gt_info.split(":")[0]

            # # we now set all gt from partition to 0, and reset the csv gt flags in the below
            # # here is to maintain the '.' gt flag that due to missing paths
            if gt_flag == ".":
                sub_snarl_gt_list.append(".")
            else:
                sub_snarl_gt_list.append("0")

        sub_snarl_include_paths = []
        for path in path_subsnarl_dict:
            if sub_snarl_id in path_subsnarl_dict[path]:
                sub_snarl_include_paths.append(path)

        sub_snarl_include_asm_names = []
        for path in sub_snarl_include_paths:
            sub_snarl_include_asm_names.extend(snarl.path_asm_dict[path])

        for asm_name in sub_snarl_include_asm_names:
            sub_snarl_gt_list[options.alt_asm_names.index(asm_name)] = "1"

        # # output
        csv_bkps = ["{}_{}_{}_{}_{}_{}_{}_{}".format(ssv.type, ssv.length, ssv.source_chr, ssv.source_start, ssv.source_end, ssv.insert_chr, ssv.insert_start, ssv.insert_end) if ssv.type in ["DUP", "invDUP"] else "{}_{}_{}_{}_{}".format(ssv.type, ssv.length, ssv.chr, ssv.start, ssv.end) for ssv in sub_snarl_csv.ssv_components]

        # alt_path_population_allele = round(len(csv_include_asm_names) / len(options.alt_asm_names), 3)
        csv_AC, csv_AF, csv_AN, csv_NS = calculate_ac_af_an_ns(sub_snarl_gt_list, options)

        if snarl.reversed_mapping:
            if sub_snarl_csv.score != "LowQual":
                sub_snarl_csv.score = "MediumQual"

        vcf_records_split.append([str(sub_snarl_csv.chr), str(sub_snarl_csv.start), "{}_{}".format(top_snarl_id, sub_snarl_id_only_start_end), sub_snarl_ref_path_rm_start_end, sub_snarl_alt_path_rm_start_end, ".", sub_snarl_csv.score, "END={};SVLEN={};SVTYPE={};BKPS={};AC={};AF={};AN={};NS={};LV=0".format(sub_snarl_csv.end, sub_snarl_csv.length, sub_snarl_csv.type, ",".join(csv_bkps), csv_AC, csv_AF, csv_AN, csv_NS), "GT", "\t".join(sub_snarl_gt_list)])


def calculate_ac_af_an_ns(gt_list, options):
    AC = 0
    non_called_cnt = 0
    for gt_info in gt_list:
        gt_flag = gt_info.split(":")[0]
        if gt_flag == ".":
            non_called_cnt += 1
        else:
            if gt_flag != "0":
                AC += 1

    AN = len(options.alt_asm_names) - non_called_cnt
    AF = round(AC / AN, 5)
    NS = len(options.spec_samples)

    return AC, AF, AN, NS


def generate_vcf_records(interpret_res_dict, snarls, options):

    vcf_records = []
    vcf_records_split = []

    for top_snarl_id in interpret_res_dict:

        snarl = snarls[top_snarl_id]
        snarl_subsnarl_csv_dict = interpret_res_dict[top_snarl_id]

        # # STEP: check if there are alt_paths that fail the prediction proess
        failed_prediction_alt_paths = []
        for sub_snarl_id in snarl_subsnarl_csv_dict:
            alt_path_csv = snarl_subsnarl_csv_dict[sub_snarl_id]

            if isinstance(alt_path_csv, dict):
                failed_prediction_alt_paths.append(sub_snarl_id)

        for sub_snarl_id in failed_prediction_alt_paths:
            snarl_subsnarl_csv_dict.pop(sub_snarl_id)

        # # generate for the whole snarl
        # with open(os.path.join(options.npz_output_path, "{}.subsnarl.json".format(top_snarl_id.replace(">", "+").replace("<", "-"))), 'r') as file:
        #     path_subsnarl_dict = json.load(file)

        path_subsnarl_dict = {}
        with open(os.path.join(options.npz_output_path, "{}.subsnarl.txt".format(top_snarl_id.replace(">", "+").replace("<", "-"))), 'r') as fin:
            for line in fin:
                line_split = line.strip().split("\t")
                path_subsnarl_dict[line_split[0]] = line_split[1:]

        snarl_path_csv_dict = {}
        for path in path_subsnarl_dict:
            path_include_ssv_comps = []
            for sub_snarl_id in path_subsnarl_dict[path]:
                if sub_snarl_id in snarl_subsnarl_csv_dict:
                    path_include_ssv_comps.extend(snarl_subsnarl_csv_dict[sub_snarl_id].ssv_components)

            snarl_path_csv_dict[path] = CSV(path_include_ssv_comps)

        # # STEP: output the snarl
        output_at_top_and_sub_snarl_level(vcf_records, vcf_records_split, top_snarl_id, snarl, snarl.ref_start, snarl.ref_end, snarl_path_csv_dict, snarl_subsnarl_csv_dict, path_subsnarl_dict, options)

    return vcf_records, vcf_records_split


def interpret_prediction_to_variant_parallel(gfa_obj, parallel_bin_num, options):
    parallel_pool = Pool(options.thread_num)
    parallel_pool_res = []

    snarls = gfa_obj.snarls
    for bin_index in range(parallel_bin_num):
        # interpret_one_bin(bin_index, gfa_minigraph, options)
        # parallel_pool_res.append(parallel_pool.apply_async(interpret_one_bin, (bin_index, snarls, options)))

        # # load snarl infos for each bin
        predict_npz_path = os.path.join(options.output_path, "tmp.matrix.{}.res.npz".format(bin_index))
        predict_res_dict = np.load(predict_npz_path)

        bin_snarls = {}
        for matrix_id in predict_res_dict:
            matrix_id_split = matrix_id.split(".")
            top_snarl_id = matrix_id_split[0]

            bin_snarls[top_snarl_id] = snarls[top_snarl_id]

        parallel_pool_res.append(parallel_pool.apply_async(interpret_one_bin, (bin_index, bin_snarls, options)))

    parallel_pool.close()
    parallel_pool.join()

    # # STEP: generate vcf records
    bin_vcf_paths = []
    bin_split_vcf_paths = []

    errored_bin = []

    for res in parallel_pool_res:
        res = res.get()

        if res[0] != "finish":
            logging.error("Fail interpreting for bin {}. {}".format(res[0], res[1]))
            errored_bin.append(res[0])
            continue

        bin_hap_level_vcf_path, bin_hap_level_split_vcf_path = res[1], res[2]

        bin_vcf_paths.append(bin_hap_level_vcf_path)
        bin_split_vcf_paths.append(bin_hap_level_split_vcf_path)

    if len(errored_bin) != 0:
        logging.error("Exit due to interpreting error")
        exit(-1)

    logging.info("Outputing into final VCFs")

    # # output to population vcf file
    hap_level_vcf_path = os.path.join(options.output_path, "tmp.swave.hap_level.vcf")
    hap_level_vcf_path_fout = open(hap_level_vcf_path, "w")

    vcf_records = []

    for file in bin_vcf_paths:
        with pysam.VariantFile(file) as fin:
            for record in fin:
                vcf_records.append(str(record).strip().split("\t"))

    output_vcf_header(hap_level_vcf_path_fout, options)
    output_vcf_records(hap_level_vcf_path_fout, vcf_records)

    hap_level_vcf_path_fout.close()

    # # output to split/single vcf file
    hap_level_split_vcf_path = os.path.join(options.output_path, "tmp.swave.hap_level.split.vcf")
    hap_level_split_vcf_path_fout = open(hap_level_split_vcf_path, "w")

    vcf_records_split = []
    for file in bin_split_vcf_paths:
        with pysam.VariantFile(file) as fin:
            for record in fin:
                vcf_records_split.append(str(record).strip().split("\t"))

    output_vcf_header(hap_level_split_vcf_path_fout, options)
    output_vcf_records(hap_level_split_vcf_path_fout, vcf_records_split)

    hap_level_split_vcf_path_fout.close()

    # # generate sample-level vcf
    generate_sample_level_vcf(hap_level_vcf_path, options)
    generate_sample_level_vcf(hap_level_split_vcf_path, options, mode="split")


def interpret_add_to_variant(raw_snarls, merged_snarl_records, options):
    """
    interpret the variants from the add mode, including reset the gt tags for each record
    """

    final_records = []

    for snarl_id in merged_snarl_records:

        snarl = merged_snarl_records[snarl_id]

        snarl_include_asm_names = [asm_name for asm_list in raw_snarls[snarl_id].path_asm_dict.values() for asm_name in asm_list]

        for path in snarl:

            # # regenerate GT tags
            path_gt = ["0" if asm_name in snarl_include_asm_names else "." for asm_name in options.alt_asm_names]

            path_include_names = raw_snarls[snarl_id].path_asm_dict[path]
            for asm_name in path_include_names:
                if asm_name == options.ref_asm_name:
                    continue
                path_gt[options.alt_asm_names.index(asm_name)] = "1"

            # # output records
            path_record = snarl[path]
            if path_record == "reference":
                continue

            path_record = path_record.strip().split("\t")[: 9]
            path_record.extend(path_gt)

            final_records.append(path_record)

    final_output_path = os.path.join(os.path.join(options.output_path, "swave.hap_level.split.merged.vcf"))
    final_output_fout = open(final_output_path, "w")

    output_vcf_header(final_output_fout, options)
    output_vcf_records(final_output_fout, final_records)

    final_output_fout.close()

    generate_sample_level_vcf(final_output_path, options, mode="split")

    if os.path.exists(os.path.join(options.output_path, "swave.sample_level.split.vcf")):
        os.remove(os.path.join(options.output_path, "swave.sample_level.split.vcf"))

    if os.path.exists(os.path.join(options.output_path, "swave.sample_level.vcf")):
        os.remove(os.path.join(options.output_path, "swave.sample_level.vcf"))
