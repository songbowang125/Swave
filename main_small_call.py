import datetime
import multiprocessing
import os
import shutil

import numpy as np
from intervaltree import IntervalTree
import pysam
import re
import argparse
import sys


class BED_Row:

    def __init__(self, row_string):
        row_string_split = row_string.strip().split("\t")

        self.ref_chrom, self.ref_start, self.ref_end = row_string_split[0], int(row_string_split[1]), int(row_string_split[2])

        self.start_node, self.end_node = row_string_split[3], row_string_split[4]

        contig_info = row_string_split[5].split(":")

        self.contig_path = contig_info[0]

        if self.contig_path == ".":  # # why alt_path is '.': because the asm has no contig covering this snarl
            self.contig_orient, self.contig_name, self.contig_start, self.contig_end = None, None, None, None
        else:
            self.contig_orient, self.contig_name, self.contig_start, self.contig_end = contig_info[2], contig_info[3], int(contig_info[4]), int(contig_info[5])

    def to_string(self):

        return "{}\t{}".format(self.start_node, self.end_node)




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

def load_from_sample_tsv(input_tsv_path):

    input_dict = {}
    input_list = []
    input_hap_numbers = []

    with open(input_tsv_path) as fin:
        for line in fin:
            if "NAME" in line:
                continue

            line_split = line.strip().split()

            sample_name = line_split[0]
            hap_paths = line_split[1: ]

            input_dict[sample_name] = [os.path.basename(path) for path in hap_paths]
            input_list.extend(hap_paths)

            input_hap_numbers.append(len(hap_paths))

    return input_dict, input_list


def search_snarl_same_path(ref_path, asm_path_list, bed_rows_dict):
    """
    in one snarl, if the sample's path is same with ref's path, then we detect small vars
    """

    samepath_candiate_dict = {}

    # # load ref bed info
    ref_bed_rows = bed_rows_dict[ref_path]

    ref_dict = {}
    for i in range(len(ref_bed_rows) - 1):
        cur_row, next_row = ref_bed_rows[i], ref_bed_rows[i + 1]
        if cur_row.contig_name is None:
            continue
        if next_row.contig_name is not None and max(cur_row.contig_start, next_row.contig_start) <= min(cur_row.contig_end, next_row.contig_end):
            continue
        ref_dict[cur_row.end_node] = cur_row

    # # load sample bed info
    for asm_path in asm_path_list:

        asm_bed_rows = bed_rows_dict[asm_path]

        for i in range(len(asm_bed_rows) - 1):
            cur_row, next_row = asm_bed_rows[i], asm_bed_rows[i + 1]

            node = cur_row.end_node
            if node not in ref_dict:
                continue

            if cur_row.contig_name is None:
                continue
            if next_row.contig_name is not None and max(cur_row.contig_start, next_row.contig_start) <= min(cur_row.contig_end, next_row.contig_end):
                continue

            # # same path found
            if cur_row.contig_path == ref_dict[node].contig_path:
                if node + "_snarl" not in samepath_candiate_dict:
                    samepath_candiate_dict[node + "_snarl"] = {ref_path: [ref_dict[node].contig_name, max(0, ref_dict[node].contig_start - 500), ref_dict[node].contig_end + 500, ref_dict[node].contig_orient]}

                samepath_candiate_dict[node + "_snarl"][asm_path] = [cur_row.contig_name, max(0, cur_row.contig_start - 500), cur_row.contig_end + 500, cur_row.contig_orient]

    return samepath_candiate_dict



def search_non_snarl(ref_path, asm_path_list, bed_rows_dict):

    nonsnarl_candidate_dict = {}

    # # STEP: traverse each sample

    for asm_path in [ref_path] + asm_path_list:

        # # STEP: load sample information
        bed_rows = bed_rows_dict[asm_path]

        # # STEP: get nodes
        for i in range(len(bed_rows) - 1):

            cur_row, next_row = bed_rows[i], bed_rows[i + 1]

            if cur_row.contig_name is None or next_row.contig_name is None:
                continue
            if max(cur_row.contig_start, next_row.contig_start) <= min(cur_row.contig_end, next_row.contig_end):
                continue

            node, contig_orient = cur_row.end_node, cur_row.contig_orient

            if node not in nonsnarl_candidate_dict:
                nonsnarl_candidate_dict[node] = {}

            # # meet the first record for the contig
            # if i == 0 or (i != 0 and bed_rows[i - 1].contig_name != None and bed_rows[i - 1].contig_name != cur_row.contig_name):
            #     if cur_row.contig_orient == "+":
            #         node_contig_range = [cur_row.contig_name, 0, cur_row.contig_start, contig_orient]
            #     else:
            #         with pysam.FastaFile(asm_path) as fin:
            #             contig_length = fin.get_reference_length(cur_row.contig_name)
            #         node_contig_range = [cur_row.contig_name, cur_row.contig_end, contig_length, contig_orient]

            # # meet the last record for the contig
            if cur_row.contig_name != next_row.contig_name:
                # # also different with previous one, indicating this is  alone
                pre_row = bed_rows[i - 1]
                if cur_row.contig_name != pre_row.contig_name:
                    continue


                if cur_row.contig_orient == "+":
                    with pysam.FastaFile(asm_path) as fin:
                        contig_length = fin.get_reference_length(cur_row.contig_name)
                    node_contig_range = [cur_row.contig_name, cur_row.contig_end, contig_length, contig_orient]
                else:
                    node_contig_range = [cur_row.contig_name, 0, cur_row.contig_start, contig_orient]

            else:
                if next_row.contig_start > cur_row.contig_end:
                    node_contig_range = [cur_row.contig_name, cur_row.contig_end, next_row.contig_start, contig_orient]
                elif cur_row.contig_start > next_row.contig_end:
                    node_contig_range = [cur_row.contig_name, next_row.contig_end, cur_row.contig_start, contig_orient]
                else:
                    print("Warning: unable to fetch contig cords at ",  cur_row.to_string())
                    continue

            nonsnarl_candidate_dict[node][asm_path] = node_contig_range

    return nonsnarl_candidate_dict


def find_intervals(nums):
    if not nums:
        return []

    nums = sorted(set(nums))

    res = []
    start = nums[0]

    for i in range(1, len(nums)):
        if nums[i] != nums[i - 1] + 1:
            res.append((start, nums[i - 1]))
            start = nums[i]

    res.append((start, nums[-1]))
    return res

def search_unmapped_snarl(ref_path, bed_rows_dict, unmapped_bed_rows_dict):

    # # STEP: collect
    unmapped_candidate_dict = {}

    for asm_path in unmapped_bed_rows_dict:
        for unmapped_start_line, unmapped_end_line in find_intervals(unmapped_bed_rows_dict[asm_path]):
            try:
                ref_start_bed_row, ref_end_bed_row = bed_rows_dict[ref_path][unmapped_start_line], bed_rows_dict[ref_path][unmapped_end_line + 1]
                ref_start_bed_row_contig, ref_start_bed_row_start, ref_start_bed_row_end, ref_start_orient = ref_start_bed_row.contig_name, ref_start_bed_row.contig_start, ref_start_bed_row.contig_end, ref_start_bed_row.contig_orient
                ref_end_bed_row_contig, ref_end_bed_row_start, ref_end_bed_row_end, ref_end_bed_row_orient = ref_end_bed_row.contig_name, ref_end_bed_row.contig_start, ref_end_bed_row.contig_end, ref_end_bed_row.contig_orient

                asm_start_bed_row, asm_end_bed_row = bed_rows_dict[asm_path][unmapped_start_line], bed_rows_dict[asm_path][unmapped_end_line + 1]
                asm_start_bed_row_contig, asm_start_bed_row_start, asm_start_bed_row_end, asm_start_orient = asm_start_bed_row.contig_name, asm_start_bed_row.contig_start, asm_start_bed_row.contig_end, asm_start_bed_row.contig_orient
                asm_end_bed_row_contig, asm_end_bed_row_start, asm_end_bed_row_end, asm_end_bed_row_orient = asm_end_bed_row.contig_name, asm_end_bed_row.contig_start, asm_end_bed_row.contig_end, asm_end_bed_row.contig_orient

                if asm_start_bed_row_contig != asm_end_bed_row_contig or ref_start_bed_row_contig != ref_end_bed_row_contig:
                    continue

                if asm_start_orient != asm_end_bed_row_orient or ref_start_orient != ref_end_bed_row_orient:
                    continue

                candidate_id = "{}-{}".format(unmapped_start_line, unmapped_end_line)
                if candidate_id not in unmapped_candidate_dict:
                    unmapped_candidate_dict[candidate_id] = {}

                if ref_start_orient == "+":
                    unmapped_candidate_dict[candidate_id][ref_path] = [ref_start_bed_row_contig, ref_start_bed_row_end, ref_end_bed_row_start, ref_start_orient]
                else:
                    unmapped_candidate_dict[candidate_id][ref_path] = [ref_start_bed_row_contig, ref_end_bed_row_end, ref_start_bed_row_start, ref_start_orient]

                if asm_start_orient == "+":
                    unmapped_candidate_dict[candidate_id][asm_path] = [asm_start_bed_row_contig, asm_start_bed_row_end, asm_end_bed_row_start, asm_start_orient]
                else:
                    unmapped_candidate_dict[candidate_id][asm_path] = [asm_start_bed_row_contig, asm_end_bed_row_end, asm_start_bed_row_start, asm_start_orient]

            except IndexError:
                continue

    return unmapped_candidate_dict


def rename_contig(raw_contig_name):
    INVALID_CHARS = r'[<>:"/\\|?*]'

    new_contig_name = re.sub(INVALID_CHARS, "-", raw_contig_name)

    new_contig_name = new_contig_name.rstrip(' .')

    return new_contig_name


def generate_seq_for_candidate(ref_path, id, candidate_var, output_path):

    ref_cord = candidate_var[ref_path]

    output_prefix = os.path.join(output_path, "tmp.{}+++{}+++{}+++{}".format(rename_contig(ref_cord[0]), ref_cord[1], ref_cord[2], id[1:]))

    ref_seq_path = output_prefix + ".ref.fa"
    alt_seq_path = output_prefix + ".alt.fa"

    ref_seq_length = ref_cord[2] - ref_cord[1]

    if ref_seq_length <= 0:
        return 0

    with pysam.FastaFile(ref_path) as fin, open(ref_seq_path, "w") as fout:
        fout.write(">ref\n".format(os.path.basename(ref_path)))

        if ref_cord[3] == "+":
            fout.write("{}\n".format(fin.fetch(ref_cord[0], ref_cord[1], ref_cord[2]).upper()))
        else:
            fout.write("{}\n".format(reverse_complement_seq(fin.fetch(ref_cord[0], ref_cord[1], ref_cord[2]).upper())))

    for alt_path in candidate_var:
        if alt_path == ref_path:
            continue

        alt_cord = candidate_var[alt_path]

        alt_seq_length = alt_cord[2] - alt_cord[1]
        if alt_seq_length <= 0:
            continue
        if alt_seq_length / ref_seq_length > 3:
            continue

        with pysam.FastaFile(alt_path) as fin, open(alt_seq_path, "a") as fout:
            fout.write(">{}\n".format(os.path.basename(alt_path)))
            if alt_cord[3] == "+":
                fout.write("{}\n".format(fin.fetch(alt_cord[0], alt_cord[1], alt_cord[2]).upper()))
            else:
                fout.write("{}\n".format(reverse_complement_seq(fin.fetch(alt_cord[0], alt_cord[1], alt_cord[2]).upper())))

def cigar_to_list(cigar):
    """
    convert cigar string to list
    :return:
    """

    opVals = re.findall(r'(\d+)([\w=])', cigar)
    lengths = [int(opVals[i][0]) for i in range(0, len(opVals))]
    ops = [opVals[i][1] for i in range(0, len(opVals))]

    return ops, lengths


def detect_small_variants(file, output_path, minimap2_path):

    ref_file = os.path.join(output_path, file)
    alt_file = ref_file.replace(".ref.fa", ".alt.fa")

    if not os.path.exists(alt_file):
        return 0

    # # run alignment
    ref_cord_split = ref_file.replace(".ref.fa", "").split("+++")
    ref_contig, ref_start, ref_end = ref_cord_split[0], int(ref_cord_split[1]), int(ref_cord_split[2])

    out_paf_file = ref_file.replace(".ref.fa", ".paf")

    os.system("{} -t 4 -Y -c --eqx -x asm20 {} {} > {} 2> /dev/null".format(minimap2_path, ref_file, alt_file, out_paf_file))
    # os.system("minimap2 -t 4 -Y -c --eqx -x asm20 --secondary=no -s 25000 -K 8G  {} {} > {} 2> /dev/null".format(ref_file, alt_file, out_paf_file))

    # # parse alignment cigar string
    out_call_file = ref_file.replace(".ref.fa", ".call.bed")

    with open(out_paf_file) as paf_fin, pysam.FastaFile(ref_file) as ref_fin, pysam.FastaFile(alt_file) as alt_fin, open(out_call_file, "w") as fout:
        for line in paf_fin:
            if line.startswith("@"):
                continue

            line_split = line.strip().split()

            alt_contig, alt_shift, ref_shift = line_split[0], int(line_split[2]), int(line_split[7])

            cigar_string = line_split[21]
            if "cg:Z" not in cigar_string:
                cigar_string = line_split[22]

            cigar_ops, cigar_ops_length = cigar_to_list(cigar_string.replace("cg:Z:", ""))

            ref_pointer = ref_shift
            alt_pointer = alt_shift

            for i in range(len(cigar_ops)):
                op = cigar_ops[i]
                op_len = cigar_ops_length[i]

                if op == "N" or op == "S":
                    alt_pointer += op_len

                elif op == "I":

                    if op_len < 50:
                        ref_bases = ref_fin.fetch("ref", ref_pointer - 1, ref_pointer)
                        alt_bases = alt_fin.fetch(alt_contig, alt_pointer - 1, alt_pointer + op_len)

                        fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(alt_contig, ref_start + ref_pointer, ref_start + ref_pointer, "INS", op_len, ref_bases, alt_bases))

                    alt_pointer += op_len

                elif op == "D":
                    if op_len < 50:
                        ref_bases = ref_fin.fetch("ref", ref_pointer - 1, ref_pointer + op_len)
                        alt_bases = alt_fin.fetch(alt_contig, alt_pointer - 1, alt_pointer)

                        fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(alt_contig, ref_start + ref_pointer, ref_start + ref_pointer + op_len, "DEL", op_len, ref_bases, alt_bases))

                    ref_pointer += op_len

                elif op in ["M", '=']:
                    ref_pointer += op_len
                    alt_pointer += op_len

                elif op in ["X", "E"]:
                    if op_len < 50:

                        for snp_index in range(op_len):

                            ref_bases = ref_fin.fetch("ref", ref_pointer, ref_pointer + 1)
                            alt_bases = alt_fin.fetch(alt_contig, alt_pointer, alt_pointer + 1)

                            fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(alt_contig, ref_start + ref_pointer + 1, ref_start + ref_pointer + 1, "SNP", 1, ref_bases, alt_bases))

                            ref_pointer += 1
                            alt_pointer += 1
                    else:
                        ref_pointer += op_len
                        alt_pointer += op_len

                else:
                    continue



def hap_gt_to_sample_gt(sample_gt, mode="raw"):

    if mode == "split":
        gt_flags = []
        for hap_gt in sample_gt:
            gt_flags.append(hap_gt)

        return "{}".format("|".join(gt_flags))

    else:
        gt_flags = []


        for hap_gt in sample_gt:
            hap_gt_split = hap_gt.split(":")

            gt_flags.append(hap_gt_split[0])

        return "{}".format("|".join(gt_flags))

def generate_sample_level_vcf(hap_level_path, sample_list, input_dict, mode="raw"):

    sample_level_vcf_path = hap_level_path.replace("hap_level", "sample_level").replace("tmp.", "")

    sample_names = list(input_dict.keys())

    sample_hap_files = []
    for sample in input_dict:
        sample_hap_files.append(input_dict[sample])

    with pysam.VariantFile(hap_level_path) as fin, open(sample_level_vcf_path, "w") as fout:
        # # write header
        fout.write("\n".join(str(fin.header).split("\n")[: -2]) + "\n")
        last_header_line = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
        last_header_line.extend(sample_names)
        fout.write("\t".join(last_header_line) + "\n")

        # # parse the header of the hap level vcf and find the corresponding index
        hap_names = str(fin.header).strip().split("\n")[-1].split("\t")[9: ]

        hap_names_index_on_sample = []    # # store

        for hap_name in hap_names:

            found_index = None
            for index in range(len(sample_hap_files)):

                hap_files = sample_hap_files[index]

                for hap_file in hap_files:
                    if hap_name in hap_file:
                        found_index = index
                        break

                if found_index != None:
                    break

            hap_names_index_on_sample.append(found_index)

            if found_index == None:
                print("Not found this hap file")
                exit()

        # # parse each record
        for record in fin:

            sample_level_gts = [[] for i in range(len(sample_names))]

            record_gts = str(record).strip().split("\t")[9: ]

            for index in range(len(record_gts)):

                sample_level_gts[hap_names_index_on_sample[index]].append(record_gts[index])

            sample_level_gts = [hap_gt_to_sample_gt(hap_gt_list, mode=mode) for hap_gt_list in sample_level_gts]

            fout.write("\t".join(str(record).strip().split("\t")[ :9]) + "\t" + "\t".join(sample_level_gts) + "\n")




def output_vcf_header(vcf_fout, ref_path, sample_list):
    """
    generate vcf header
    """

    print("##fileformat=VCFv4.3", file=vcf_fout)

    # # STEP: add chromosome info
    ref_file = pysam.FastaFile(ref_path)
    chroms = ref_file.references
    for chr in chroms:
        chr_length = ref_file.get_reference_length(chr)
        print('##contig=<ID={0},length={1}>'.format(chr, chr_length), file=vcf_fout)
    ref_file.close()

    # # STEP: add info field
    # print("##QUAL=<QUAL=XXX,Description=\"The SV quality of the SV described in this region\">", file=vcf_fout)

    print("##ID=<ID=Path,Description=\"ID\">", file=vcf_fout)
    print("##REF=<ID=Path,Description=\"Reference bases\">", file=vcf_fout)
    print("##ALT=<ID=Path,Description=\"Alternative bases\">", file=vcf_fout)

    print("##FILTER=<ID=PASS,Description=\"Passed\">", file=vcf_fout)

    print("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the SV/CSV\">", file=vcf_fout)
    print("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV/CSV\">", file=vcf_fout)
    print("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the SV/CSV\">", file=vcf_fout)

    print("##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Total number of alternate alleles in called genotypes\">", file=vcf_fout)
    print("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Estimated allele frequency in the range (0,1]\">", file=vcf_fout)
    print("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">", file=vcf_fout)

    # # STEP: add gt info
    print("##FORMAT=<ID=GT,Number=.,Type=String,Description=\"SV genotype for each haplotype\">", file=vcf_fout)


    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}".format( "\t".join(sample_list)), file=vcf_fout)


def output_vcf(input_path, ref_path, output_path):
    input_dict, input_list = load_from_sample_tsv(input_path)

    # # convert file path to file name
    sample_list = [os.path.basename(input) for input in input_list]

    # # begin
    vcf_path = os.path.join(output_path, "swave_small.hap_level.vcf")
    vcf_split_path = os.path.join(output_path, "swave_small.hap_level.split.vcf")

    with open(vcf_path, "w") as fout, open(vcf_split_path, "w") as fout_split:

        output_vcf_header(fout, ref_path, sample_list)
        output_vcf_header(fout_split, ref_path, sample_list)

        var_cnt = 0
        var_split_cnt = 0

        for call_bed in os.listdir(output_path):

            if ".call.bed" not in call_bed:
                continue

            chrom = call_bed.split(".")[1].split("+++")[0]

            var_collect = {}
            var_collect_split = {}

            with open(os.path.join(output_path, call_bed)) as fin:

                for line in fin:
                    line_split = line.strip().split("\t")
                    sample, start, end, type, length, ref_bases, alt_base = line_split[ : ]

                    allele_id = "{}-{}".format(ref_bases, alt_base)
                    pos_id = "{}-{}-{}-{}".format(chrom, start, end, type)
                    var_id = "{}-{}".format(pos_id, allele_id)

                    # # STEP: update at position-level (also the snarl level for SV)
                    if pos_id not in var_collect:
                        var_collect[pos_id] = {}
                    if allele_id not in var_collect[pos_id]:
                        var_collect[pos_id][allele_id] = []
                    var_collect[pos_id][allele_id].append(sample)

                    # # STEP: update at the split/allele level (also the split or allele level for SV)
                    if var_id not in var_collect_split:
                        var_collect_split[var_id] = []
                    var_collect_split[var_id].append(sample)

            # # STEP: output at position-level (also the snarl level for SV)
            for pos_id in sorted(var_collect.keys()):
                chrom, start, end, type = pos_id.split("-")[:]
                include_alleles = list(var_collect[pos_id].keys())

                # # generate ref and alt bases
                ref_bases = include_alleles[0].split("-")[0]
                alt_bases_list = []
                ac = 0

                # # generate gt flag
                sample_gts = ["0" for i in range(len(sample_list))]
                for allele_cnt in range(len(include_alleles)):
                    allele = include_alleles[allele_cnt]
                    alt_bases_list.append(allele.split("-")[1])

                    for sample in var_collect[pos_id][allele]:
                        sample_gts[sample_list.index(sample)] = str(allele_cnt + 1)     # # must + 1, because the ref allele in 0

                        ac += 1

                # # output
                ns = len(sample_list)
                af = round(ac / ns, 5)
                fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\tAC={};AF={};NS={}\tGT\t{}\n".format(chrom, start, var_cnt, ref_bases, ",".join(alt_bases_list), ".", ".", ac, af, ns, "\t".join(sample_gts)))

                var_cnt += 1

            # # STEP: output the split/allele level (also the split or allele level for SV)
            for var_id in sorted(var_collect_split.keys()):
                chrom, start, end, type, ref_bases, alt_base = var_id.split("-")[:]
                if len(ref_bases) == len(alt_base):
                    svlen = len(alt_base)
                else:
                    svlen = abs(len(alt_base) - len(ref_bases))

                include_samples = var_collect_split[var_id]

                # # generate gt flag
                sample_gts = ["0" for i in range(len(sample_list))]
                for sample in include_samples:
                    sample_gts[sample_list.index(sample)] = "1"

                # # output
                ac = len(include_samples)
                ns = len(sample_list)
                af = round(ac / ns, 5)
                fout_split.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\tSVLEN={};SVTYPE={};AC={};AF={};NS={}\tGT\t{}\n".format(chrom, start, var_split_cnt, ref_bases, alt_base, ".", ".", svlen, type, ac, af, ns, "\t".join(sample_gts)))

                var_split_cnt += 1

    try:
        generate_sample_level_vcf(vcf_path, sample_list, input_dict, )
        generate_sample_level_vcf(vcf_split_path, sample_list, input_dict, mode="split")
    except Exception:

        print("[Warning]: failed to generate sample-level VCF, but you can still use the hap-level VCF for analysis. Thanks for returning the error to the author")
        exit(-1)


def compare(base_vcf, target_vcf, hc_bed):
    def build_tree(hc_bed):
        trees = {}

        with open(hc_bed) as fin:

            for line in fin:
                line_split = line.strip().split("\t")

                chrom, start, end = line_split[0], int(line_split[1]), int(line_split[2])
                if chrom not in trees:
                    trees[chrom] = IntervalTree()
                trees[chrom][start: end + 1] = True
        return trees

    def in_hc_tree(chr, pos, trees):
        return chr in trees and len(trees[chr][pos]) > 0

    if hc_bed is not None:
        hc_trees = build_tree(hc_bed)
    else:
        hc_trees = {}

    # # begin
    base_records = []
    with pysam.VariantFile(base_vcf) as fin:
        for record in fin:
            chrom, start, ref_bases, alt_bases = record.contig, record.start + 1, record.ref.upper(), record.alts[0].upper()
            if len(ref_bases) >= 50 or len(alt_bases) >= 50:
                continue

            if hc_bed is not None and not in_hc_tree(chrom, start, hc_trees):
                continue

            gt_flag = str(record).strip().split("\t")[9].split(":")[0]
            if "." in gt_flag:
                continue

            base_records.append("{}-{}-{}-{}".format(chrom, start, ref_bases, alt_bases))
            # base_records.append("{}-{}".format(chrom, start, ref_bases, alt_bases))

    target_records = []
    with pysam.VariantFile(target_vcf) as fin:
        for record in fin:
            chrom, start, ref_bases, alt_bases = record.contig, record.start + 1, record.ref.upper(), record.alts[0].upper()
            if len(ref_bases) >= 50 or len(alt_bases) >= 50:
                continue

            if hc_bed is not None and not in_hc_tree(chrom, start, hc_trees):
                continue

            target_records.append("{}-{}-{}-{}".format(chrom, start, ref_bases, alt_bases))
            # target_records.append("{}-{}".format(chrom, start, ref_bases, alt_bases))


    base_records = set(base_records)
    target_records = set(target_records)

    tp = len(base_records & target_records)
    fp = len(target_records - base_records)
    fn = len(base_records - target_records)

    p = round(tp / (tp + fp), 3) if tp + fp else 0
    r = round(tp / (tp + fn), 3) if tp + fn else 0
    f1 = round(2 * p * r / (p + r), 3) if p + r else 0

    print(f1, p, r, tp, fp, fn, len(base_records), len(target_records))
    #
    # for record in base_records - target_records:
    #     print(record.replace("-", "\t"))

    # for record in target_records - base_records:
    #     print(record.replace("-", "\t"))

def load_minigraph_beds(ref_path, asm_path_list, bed_folder):

    # # STEP: traverse each sample

    bed_rows_dict = {}
    unmapped_bed_rows_dict = {}  # # format: {asm1: [line1, line2], asm2: [line3, line4]}

    for asm_path in [ref_path] + asm_path_list:
        asm_name = os.path.basename(asm_path)
        asm_bed_path = os.path.join(bed_folder, "{}.bed".format(asm_name))

        if not os.path.exists(asm_bed_path):
            print("Not find 'Minigrah --call' BED for {}. Please run minigraph --call, and store the BED in {}".format(asm_name, bed_folder))
            exit(-1)

        if asm_path not in bed_rows_dict:
            bed_rows_dict[asm_path] = []

        with open(asm_bed_path) as fin:
            line_cnt = 0
            for line in fin:
                bed_row = BED_Row(line)

                # # find unmapped
                if bed_row.contig_path == ".":
                    if asm_path not in unmapped_bed_rows_dict:
                        unmapped_bed_rows_dict[asm_path] = []
                    unmapped_bed_rows_dict[asm_path].append(line_cnt)
                    unmapped_bed_rows_dict[asm_path].append(line_cnt - 1)    # # why: the unmapped row will affect the above row in the downstream codes

                bed_rows_dict[asm_path].append(bed_row)

                line_cnt += 1


    return bed_rows_dict, unmapped_bed_rows_dict


def parse_arguments(arguments=sys.argv[1:]):

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    inout_params = parser.add_argument_group("Input/Output parameters")

    inout_params.add_argument('--input_path', dest="input_path", type=os.path.abspath, required=True, help='Absolute path to input .tsv file')
    inout_params.add_argument('--ref_path', dest="ref_path", type=os.path.abspath, required=True, help='Absolute path to reference')
    inout_params.add_argument('--output_path', dest="output_path", type=os.path.abspath, required=True, help='Absolute path to output')

    inout_params.add_argument('--minigraph_callbed_folder', dest="minigraph_callbed_folder", type=os.path.abspath, required=True, help='Absolute path to the folder storing minigraph --call BED files')
    inout_params.add_argument('--minimap2_path', dest="minimap2_path", required=True, help='Absolute path to minimap2')

    inout_params.add_argument('--thread_num', dest="thread_num", type=int, required=True, help='Thread number')

    return parser.parse_args(arguments)

if __name__ == '__main__':

    options = parse_arguments()
    ref_path = options.ref_path
    input_path = options.input_path
    output_path = options.output_path
    minigraph_callbed_folder = options.minigraph_callbed_folder
    minimap2_path = options.minimap2_path
    thread_num = options.thread_num

    if shutil.which(minimap2_path) is None:
        print("ERROR: not found minimap2 at {}. Please check your path (format should be /path/to/minimap2)".format(minimap2_path))
        exit(-1)

    start_time = datetime.datetime.now()

    if os.path.exists(output_path):
        shutil.rmtree(output_path)

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    # # STEP: load info from bed file
    print("===========STEP1: load info===========")
    input_dict, input_list = load_from_sample_tsv(input_path)
    bed_rows_dict, unmapped_bed_rows_dict = load_minigraph_beds(ref_path, input_list, minigraph_callbed_folder)

    # # STEP: search for candidates
    print("===========STEP2: search candidates===========")
    unmapped_candidate_dict = search_unmapped_snarl(ref_path, bed_rows_dict, unmapped_bed_rows_dict)
    nonsnarl_candidate_dict = search_non_snarl(ref_path, input_list, bed_rows_dict)
    samepath_candiate_dict = search_snarl_same_path(ref_path, input_list, bed_rows_dict)

    # # STEP: generate seq
    print("===========STEP3: prepare sequence===========")
    process_pool = multiprocessing.Pool(thread_num)
    for candidate_dict in [samepath_candiate_dict, nonsnarl_candidate_dict, unmapped_candidate_dict]:
        for id in candidate_dict:
            candidate_var = candidate_dict[id]
            if ref_path not in candidate_var:
                continue
            process_pool.apply_async(generate_seq_for_candidate, (ref_path, id, candidate_var, output_path))

    process_pool.close()
    process_pool.join()

    # # STEP: detect variant
    print("===========STEP4: detect variant===========")
    process_pool = multiprocessing.Pool(int(thread_num / 4))
    for file in os.listdir(output_path):
        if ".ref.fa" not in file or ".fa.fai" in file or "+++" not in file:
            continue
        process_pool.apply_async(detect_small_variants, (file, output_path, minimap2_path))

    process_pool.close()
    process_pool.join()

    print("===========STEP5: output and clean===========")
    output_vcf(input_path, ref_path, output_path)

    for file in os.listdir(output_path):
        if file.startswith("tmp."):
            os.remove(os.path.join(output_path, file))

    end_time = datetime.datetime.now()
    print("Swave small variant calling finished. Cost {} min".format(int((end_time - start_time).seconds / 60)))
