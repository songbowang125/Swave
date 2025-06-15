import logging
import os.path
import re
import sys
import multiprocessing
import numpy as np
import pysam
from src.pack_dotplot.op_dotplot import generate_dotplot_and_projection, calculate_stride_size
from src.pack_graph.op_seq import reverse_complement_seq
import traceback


def get_snarl_path_seq(node_seq_dict, snarl_path, start_node_id=None, end_node_id=None, spec_trim_length=1000):
    """
    Given a snarl path (formatted as >s692>s694>s695), output its seq

    Note if the left and right anchor node is provided, we trim their seq length to get a much short seq
    """

    if snarl_path == "*":
        return ""

    snarl_path_seq = ""

    include_nodes = re.findall(r'[><]([a-zA-Z0-9]+)', snarl_path)
    include_nodes_orients = re.findall(r'[><]', snarl_path)

    for i in range(len(include_nodes)):

        node_id = include_nodes[i]
        node_orient = include_nodes_orients[i]

        # # forward node, direct add sequence
        if node_orient == ">":
            node_seq = node_seq_dict[node_id]
        else:
            node_seq = reverse_complement_seq(node_seq_dict[node_id])

        # # trim node if it is the left or right anchor node
        if node_id == start_node_id:
            start_node_seq_trimmed = node_seq[-spec_trim_length:]

            snarl_path_seq += start_node_seq_trimmed

        elif node_id == end_node_id:
            end_node_seq_trimmed = node_seq[: spec_trim_length]

            snarl_path_seq += end_node_seq_trimmed

        else:
            snarl_path_seq += node_seq

    return snarl_path_seq


def resolve_large_snarl_path(node_seq_dict, ref_path, alt_path, start_node_with_orient, end_node_with_orient, snarl_ref_start, snarl_ref_end):

    # # STEP: check nodes in the ref path, determine if they exist in path
    ref_path_include_nodes = re.findall(r'[><]([a-zA-Z0-9]+)', ref_path)
    ref_path_include_nodes_orients = re.findall(r'[><]', ref_path)
    ref_path_split = [start_node_with_orient] + ["{}{}".format(ref_path_include_nodes_orients[i], ref_path_include_nodes[i]) for i in range(len(ref_path_include_nodes))] + [end_node_with_orient]

    alt_path_include_nodes = re.findall(r'[><]([a-zA-Z0-9]+)', alt_path)
    alt_path_include_nodes_orients = re.findall(r'[><]', alt_path)
    alt_path_split = [start_node_with_orient] + ["{}{}".format(alt_path_include_nodes_orients[i], alt_path_include_nodes[i]) for i in range(len(alt_path_include_nodes))] + [end_node_with_orient]

    # # STEP: find the start pos for each ref nodes
    ref_pos_pointer = snarl_ref_start
    ref_path_include_nodes_pos = []
    for i in range(len(ref_path_include_nodes)):
        ref_path_include_nodes_pos.append(ref_pos_pointer)
        ref_pos_pointer += len(node_seq_dict[ref_path_include_nodes[i]])

    ref_path_split_pos = [snarl_ref_start - len(node_seq_dict[start_node_with_orient[1: ]])] + ref_path_include_nodes_pos + [snarl_ref_end]

    # # STEP: find shared nodes
    shared_nodes_with_orients = []
    for node_with_orient in ref_path_split:
        if node_with_orient in alt_path_split:
            if not (node_with_orient in [start_node_with_orient, end_node_with_orient]) and len(node_seq_dict[node_with_orient[1: ]]) < 1000:
                continue
            shared_nodes_with_orients.append(node_with_orient)

    # # STEP: generate sub snarls
    sub_snarls = []
    for i in range(len(shared_nodes_with_orients) - 1):

        sub_start_node_with_orient = shared_nodes_with_orients[i]
        sub_end_node_with_orient = shared_nodes_with_orients[i + 1]

        sub_start_pos = ref_path_split_pos[ref_path_split.index(sub_start_node_with_orient) + 1]    # # why +1: the start pos is the end pos of the start node, this is also the start pos of the next node
        sub_end_pos = ref_path_split_pos[ref_path_split.index(sub_end_node_with_orient)]

        sub_ref_path_split = ref_path_split[ref_path_split.index(sub_start_node_with_orient): ref_path_split.index(sub_end_node_with_orient) + 1]
        sub_alt_path_split = alt_path_split[alt_path_split.index(sub_start_node_with_orient): alt_path_split.index(sub_end_node_with_orient) + 1]

        if sub_ref_path_split == sub_alt_path_split:
            continue

        sub_snarls.append([sub_start_node_with_orient, sub_end_node_with_orient, sub_start_pos, sub_end_pos, "".join(sub_ref_path_split[1: -1]), "".join(sub_alt_path_split[1: -1])]) # # why [1: -1]: remove the start and end node from path

    return sub_snarls


def plot_subsnarl_by_gepard(sub_snarl_index, path_cnt, sub_ref_path_seq, sub_alt_path_seq, options):


    gepard_seq1_path = os.path.join(options.img_output_path, "gepard.path_{}.sub_{}.seq1.fa".format(path_cnt, sub_snarl_index))
    gepard_seq2_path = os.path.join(options.img_output_path, "gepard.path_{}.sub_{}.seq2.fa".format(path_cnt, sub_snarl_index))

    with open(gepard_seq1_path, "w") as seq1_fout, open(gepard_seq2_path, "w") as seq2_fout:
        seq1_fout.write(">path_{}.sub_{}.seq1\n".format(path_cnt, sub_snarl_index))
        seq1_fout.write(sub_ref_path_seq + "\n")

        seq2_fout.write(">path_{}.sub_{}.seq2\n".format(path_cnt, sub_snarl_index))
        seq2_fout.write(sub_alt_path_seq + "\n")

    gepard_ref2ref_img_path = os.path.join(options.img_output_path, "gepard.path_{}.sub_{}.ref2ref.png".format(path_cnt, sub_snarl_index))
    gepard_ref2alt_img_path = os.path.join(options.img_output_path, "gepard.path_{}.sub_{}.ref2alt.png".format(path_cnt, sub_snarl_index))

    cmd_str = "java -cp ./src/postprocess/pack_gepard/Gepard_altered_noheader.jar org.gepard.client.cmdline.CommandLine  -seq {} {} -matrix ./src/postprocess/pack_gepard/edna.mat -outfile {} -word 15 -silent".format(gepard_seq1_path, gepard_seq1_path, gepard_ref2ref_img_path)
    os.system(cmd_str)

    cmd_str = "java -cp ./src/postprocess/pack_gepard/Gepard_altered_noheader.jar org.gepard.client.cmdline.CommandLine  -seq {} {} -matrix ./src/postprocess/pack_gepard/edna.mat -outfile {} -word 15 -silent".format(gepard_seq1_path, gepard_seq2_path, gepard_ref2alt_img_path)
    os.system(cmd_str)


def plot_snarl_by_gepard(snarl, ref_asm_path, path_asm_dict, node_seq_dict, options):

    # # STEP: fetch sequence for each path
    path_seq_dict = {}
    for path in path_asm_dict:
        # # fetch sequence for this path
        if path == "*":
            path_seq_dict[path] = ""
        else:
            path_seq_dict[path] = get_snarl_path_seq(node_seq_dict, path)

    # extend_len = min(10000, 2 * max([len(path_seq_dict[path]) for path in path_seq_dict] + [snarl.ref_end - snarl.ref_start]))
    extend_len = int(1 * max([len(path_seq_dict[path]) for path in path_seq_dict] + [snarl.ref_end - snarl.ref_start]))

    # # fetch sequence for start and end node
    start_node_id, end_node_id = snarl.start_node_id, snarl.end_node_id
    start_node_orient, end_node_orient = snarl.start_node_orient, snarl.end_node_orient

    # print(">>", snarl_id, "Extend length: ", extend_len, "snarl start:", snarl.ref_start, "snarl end: ", snarl.ref_end)

    start_node_seq = get_snarl_path_seq(node_seq_dict, "{}{}".format(start_node_orient, start_node_id), start_node_id, end_node_id, spec_trim_length=extend_len)
    end_node_seq = get_snarl_path_seq(node_seq_dict, "{}{}".format(end_node_orient, end_node_id), start_node_id, end_node_id, spec_trim_length=extend_len)

    if len(start_node_seq) < extend_len:
        try:
            with pysam.FastaFile(options.ref_asm_path) as ref_file:
                start_node_seq = ref_file.fetch(snarl.ref_chrom, max(snarl.ref_start - int( 1 * extend_len), 0), snarl.ref_start).replace("N", "")
        # # for snarl.ref_chrom not in ref file
        except KeyError:
            pass

    if len(end_node_seq) < extend_len:
        try:
            with pysam.FastaFile(options.ref_asm_path) as ref_file:
                end_node_seq = ref_file.fetch(snarl.ref_chrom, snarl.ref_end, snarl.ref_end + int(1 * extend_len)).replace("N", "")
        # # for snarl.ref_chrom not in ref file
        except KeyError:
            pass

    ref_path_seq = path_seq_dict[ref_asm_path]
    gepard_ref_seq_path = os.path.join(options.img_output_path, "gepard.ref.fa")
    with open(gepard_ref_seq_path, "w") as fout:
        fout.write(">ref\n")
        fout.write(start_node_seq + ref_path_seq + end_node_seq + "\n")

    ref2ref_img_out = os.path.join(options.img_output_path, "gepard.ref2ref.png")
    cmd_str = "java -cp ./src/postprocess/pack_gepard/Gepard_altered_noheader.jar org.gepard.client.cmdline.CommandLine  -seq {} {} -matrix ./src/postprocess/pack_gepard/edna.mat -outfile {} -word 30 -silent".format(gepard_ref_seq_path, gepard_ref_seq_path, ref2ref_img_out)
    os.system(cmd_str)

    alt_path_index = 0
    for alt_path in path_seq_dict:

        if options.spec_path is not None and options.spec_path not in alt_path:
            continue

        alt_path_index += 1
        alt_path_seq = path_seq_dict[alt_path]
        gepard_alt_seq_path = os.path.join(options.img_output_path, "gepard.path_{}.fa".format(alt_path_index))
        with open(gepard_alt_seq_path, "w") as fout:
            fout.write(">path_{}\n".format(alt_path_index))
            fout.write(start_node_seq + alt_path_seq + end_node_seq + "\n")

        ref2alt_img_out = os.path.join(options.img_output_path, "gepard.ref2alt.path_{}.png".format(alt_path_index))
        cmd_str = "java -cp ./src/postprocess/pack_gepard/Gepard_altered_noheader.jar org.gepard.client.cmdline.CommandLine  -seq {} {} -matrix ./src/postprocess/pack_gepard/edna.mat -outfile {} -word 30 -silent".format(gepard_ref_seq_path, gepard_alt_seq_path, ref2alt_img_out)
        os.system(cmd_str)

        os.remove(gepard_alt_seq_path)

    os.remove(gepard_ref_seq_path)


def generate_one_snarl(top_snarl_id, snarl, node_seq_dict, snarl_cnt, total_snarl_cnt, options):
    """
    Process one snarl, generate dotplot and projections
    """

    try:

        if snarl_cnt % 1000 == 0:
            logging.info("Generating {}th snarl. In total {} snarls".format(snarl_cnt, total_snarl_cnt))

        # # STEP: find the path for ref asm
        path_asm_dict = snarl.path_asm_dict
        path_subsnarl_dict = {}

        if len(path_asm_dict) == 0:
            return "finish"

        ref_asm_path = ""
        for path in path_asm_dict:
            # # find the path for ref asm
            if options.ref_asm_name in path_asm_dict[path]:
                ref_asm_path = path
                break

        # # STEP: there is only ref asm path, also meaning that alt asm paths are same with ref asm
        if ref_asm_path != "" and len(path_asm_dict) == 1:
            return "finish"

        # # plot gepard images for refpath and altpath, just for showcase
        # # plot gepard images for refpath and altpath, just for showcase
        # plot_snarl_by_gepard(snarl, ref_asm_path, path_asm_dict, node_seq_dict, options)
        # # end
        # # end

        # print("Processing {}, {} paths covering {} haplotypes".format(snarl.to_string(), len(path_asm_dict), sum(len(value) for value in path_asm_dict.values()) - 1))

        seq_output_path = os.path.join(options.npz_output_path, "{}.seq.fa".format(top_snarl_id.replace(">", "+").replace("<", "-")))
        seq_output_fout = open(seq_output_path, "w")

        snarl_dotplot_dict = {}

        snarl_start_node_id, snarl_end_node_id = snarl.start_node_id, snarl.end_node_id
        snarl_start_node_orient, snarl_end_node_orient = snarl.start_node_orient, snarl.end_node_orient

        snarl_start_node_with_orient = "{}{}".format(snarl_start_node_orient, snarl_start_node_id)
        snarl_end_node_with_orient = "{}{}".format(snarl_end_node_orient, snarl_end_node_id)

        # # for reversed mapping snarl, since we dont no the flanking nodes, the ref start is the start point of first node in the alt path
        if snarl.reversed_mapping:
            snarl.ref_start -= len(node_seq_dict[snarl_start_node_id])
            snarl.ref_end += len(node_seq_dict[snarl_end_node_id])

            if snarl.ref_end - snarl.ref_start > options.max_sv_size:
                return "finish"

        # # traverse each alt path
        bad_cnt = 0
        bad_flag = False
        path_cnt = 0

        processed_sub_snarls = []
        processed_seqs = []
        for alt_asm_path in path_asm_dict:
            path_cnt += 1

            path_subsnarl_dict[alt_asm_path] = []
            if options.spec_path is not None and options.spec_path not in alt_asm_path:
                continue

            # # skip ref asm path
            if alt_asm_path == ref_asm_path:
                continue

            if not snarl.reversed_mapping:
                sub_snarls = resolve_large_snarl_path(node_seq_dict, ref_asm_path, alt_asm_path, snarl_start_node_with_orient, snarl_end_node_with_orient, snarl.ref_start, snarl.ref_end)
            else:
                # # for reversed mapping snarls, there is no ref path so we cannot split it into subsnarl but keep it as one
                sub_snarls = [[snarl_start_node_orient, snarl_end_node_orient, snarl.ref_start, snarl.ref_end, ref_asm_path, alt_asm_path]]

            for sub_snarl_index in range(len(sub_snarls)):

                sub_start_node_with_orient, sub_end_node_with_orient, sub_ref_start, sub_ref_end, sub_ref_path, sub_alt_path = sub_snarls[sub_snarl_index]
                sub_snarl_id, sub_ref_chrom = "{}-{}".format(sub_ref_path, sub_alt_path), snarl.ref_chrom

                if options.spec_path is not None and options.spec_path not in sub_snarl_id:
                    continue

                if sub_snarl_id in processed_sub_snarls:
                    continue

                sub_start_node_id = sub_start_node_with_orient[1: ]
                sub_end_node_id = sub_end_node_with_orient[1: ]

                # # calculate extend len by path seqs
                sub_ref_path_seq, sub_alt_path_seq = get_snarl_path_seq(node_seq_dict, sub_ref_path), get_snarl_path_seq(node_seq_dict, sub_alt_path)

                # print(ref_asm_path, alt_asm_path)
                # print(sub_ref_path, sub_alt_path)

                if sub_ref_path == "":
                    try:
                        with pysam.FastaFile(options.ref_asm_path) as ref_file:
                            sub_ref_path_seq = ref_file.fetch(sub_ref_chrom, sub_ref_start, sub_ref_end).upper()

                    # # for snarl.ref_chrom not in ref file
                    except KeyError:
                        return "finish"

                if "N" in sub_ref_path_seq:
                    return "finish"
                # print(len(sub_ref_path_seq), len(sub_alt_path_seq))

                # print("sub_snarl_id", top_snarl_id, snarl.snarl_id, sub_snarl_id, snarl.ref_chrom, snarl.ref_start, snarl.ref_end)
                # print(snarl.ref_start, snarl.ref_end, sub_ref_start, sub_ref_end)
                extend_len = min(10000, 2 * max([len(sub_ref_path_seq), len(sub_alt_path_seq), (sub_ref_end - sub_ref_start)]))

                # # generate start node seq and end node seq
                sub_start_node_seq = get_snarl_path_seq(node_seq_dict, sub_start_node_with_orient, sub_start_node_id, sub_end_node_id, spec_trim_length=extend_len)
                if len(sub_start_node_seq) < extend_len:
                    try:
                        with pysam.FastaFile(options.ref_asm_path) as ref_file:
                            sub_start_node_seq = ref_file.fetch(sub_ref_chrom, max(sub_ref_start - extend_len, 0), sub_ref_start).replace("N", "").upper()
                    # # for sub_ref_chrom not in ref file
                    except KeyError:
                        pass

                sub_end_node_seq = get_snarl_path_seq(node_seq_dict, sub_end_node_with_orient, sub_start_node_id, sub_end_node_id, spec_trim_length=extend_len)
                if len(sub_end_node_seq) < extend_len:
                    try:
                        with pysam.FastaFile(options.ref_asm_path) as ref_file:
                            sub_end_node_seq = ref_file.fetch(sub_ref_chrom, sub_ref_end, sub_ref_end + extend_len).replace("N", "").upper()
                    # # for sub_ref_chrom not in ref file
                    except KeyError:
                        pass

                # # add start and node seq to path seqs
                sub_ref_path_seq = sub_start_node_seq + sub_ref_path_seq + sub_end_node_seq
                sub_alt_path_seq = sub_start_node_seq + sub_alt_path_seq + sub_end_node_seq
                # print(len(sub_ref_path_seq), len(sub_alt_path_seq))

                # print(len(sub_ref_path_seq), len(sub_alt_path_seq), len(sub_start_node_seq), len(sub_end_node_seq))
                sub_ref_path = sub_start_node_with_orient + sub_ref_path + sub_end_node_with_orient
                sub_alt_path = sub_start_node_with_orient + sub_alt_path + sub_end_node_with_orient
                sub_snarl_id = "{}-{}".format(sub_ref_path, sub_alt_path)

                if sub_snarl_id not in processed_seqs:
                    seq_output_fout.write(">{}\n".format(sub_snarl_id.replace(">", "+").replace("<", "-")))
                    seq_output_fout.write(sub_ref_path_seq + "\n")
                    processed_seqs.append(sub_snarl_id)

                if sub_alt_path not in processed_seqs:
                    seq_output_fout.write(">{}\n".format(sub_alt_path.replace(">", "+").replace("<", "-")))
                    seq_output_fout.write(sub_alt_path_seq + "\n")
                    processed_seqs.append(sub_alt_path)

                # # generate dotplot
                sub_dotplot_ref_start, sub_dotplot_ref_end = sub_ref_start - len(sub_start_node_seq), sub_ref_end + len(sub_end_node_seq)
                sub_dotplot_stride_size = calculate_stride_size(sub_ref_path_seq, sub_alt_path_seq)

                sub_dotplot_id = "{}.{}.{}.{}.{}.{}.{}.{}".format(top_snarl_id, sub_ref_start, sub_ref_end, sub_ref_chrom, sub_dotplot_ref_start, sub_dotplot_ref_end, sub_ref_path, sub_alt_path)
                sub_dotplot_output_prefix = os.path.join(options.img_output_path, sub_dotplot_id)

                # # the following codes are for generate dotplot from gepard
                # # the following codes are for generate dotplot from gepard
                # plot_subsnarl_by_gepard(sub_snarl_index, path_cnt, sub_ref_path_seq, sub_alt_path_seq, options)
                # #
                # #

                dotplot_projection_dict = generate_dotplot_and_projection(sub_ref_path_seq, sub_alt_path_seq, sub_ref_start, sub_dotplot_output_prefix, sub_dotplot_stride_size, options)

                processed_sub_snarls.append(sub_snarl_id)
                path_subsnarl_dict[alt_asm_path].append(sub_snarl_id)

                if not (dotplot_projection_dict['ref2alt'] == "Bad" and dotplot_projection_dict['alt2ref'] == "Bad"):
                    for dotplot_type in dotplot_projection_dict:
                        snarl_dotplot_dict[sub_dotplot_id + ".{}.{}".format(dotplot_type, sub_dotplot_stride_size)] = dotplot_projection_dict[dotplot_type]
                else:
                    bad_cnt += 1

                # break

        if not bad_flag:
            save_npz_path = os.path.join(options.npz_output_path, "{}.npz".format(top_snarl_id.replace(">", "+").replace("<", "-")))

            save_subsnarl_path = os.path.join(options.npz_output_path, "{}.subsnarl.txt".format(top_snarl_id.replace(">", "+").replace("<", "-")))

            try:
                np.savez(save_npz_path, **snarl_dotplot_dict)

                with open(save_subsnarl_path, "w") as fout:
                    for alt_asm_path in path_subsnarl_dict:
                        fout.write("{}\t{}\n".format(alt_asm_path, "\t".join(path_subsnarl_dict[alt_asm_path])))

                # with open(save_subsnarl_path, 'w') as f:
                #     json.dump(path_subsnarl_dict, f, indent=4)

            except:
                if os.path.exists(seq_output_path):
                    os.remove(seq_output_path)

                if os.path.exists(save_npz_path):
                    os.remove(save_npz_path)

                if os.path.exists(save_subsnarl_path):
                    os.remove(save_subsnarl_path)
        else:
            if os.path.exists(seq_output_path):
                os.remove(seq_output_path)

        # # large memory collection
        dotplot_projection_dict = None
        snarl_dotplot_dict = None

        return "finish"

    except:
        error_type, error_value, error_trace = sys.exc_info()
        error_log = "Error log: " + str(error_type) + ": " + str(error_value) + '. ' + 'Locate At: ' + str(traceback.extract_tb(error_trace))
        error_info = [top_snarl_id , error_log]

        # # large memory collection
        dotplot_projection_dict = None
        snarl_dotplot_dict = None

        return error_info


def generate_projection_matrix_for_snarls_parallel(gfa, options):

    # multiprocessing.set_start_method('spawn')
    parallel_pool = multiprocessing.Pool(options.thread_num)
    parallel_pool_res = []

    node_seq_file = pysam.FastaFile(gfa.gfa2fa_path)

    # # STEP: traverse each snarl
    total_snarl_cnt = len(gfa.snarls)

    cnt = 0

    for snarl_id, snarl in sorted(gfa.snarls.items(), key=lambda x: (x[1].ref_end - x[1].ref_start), reverse=True):

        if options.spec_snarl is not None and snarl_id != options.spec_snarl:
            continue

        if snarl.ref_end - snarl.ref_start > options.max_sv_size:
            continue

        # if snarl.ref_chrom not in options.available_chroms:
        #     continue
        try:
            node_seq_dict = {}
            for node in snarl.get_available_nodes():
                if node in node_seq_dict:
                    continue
                node_seq_dict[node] = node_seq_file.fetch(node).upper()
        except:
            node_seq_dict = {}
            continue

        # generate_one_snarl(snarl_id, snarl, node_seq_dict, cnt, total_snarl_cnt, options)
        parallel_pool_res.append(parallel_pool.apply_async(generate_one_snarl, (snarl_id, snarl, node_seq_dict, cnt, total_snarl_cnt, options)))

        cnt += 1

    parallel_pool.close()
    parallel_pool.join()

    node_seq_file.close()

    # # STEP: deal with errors from process pool
    for res in parallel_pool_res:
        res = res.get()

        if res != "finish":
            logging.error("Skip generating for snarl {}. {}".format(res[0], res[1]))
