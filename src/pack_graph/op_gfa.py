import os
import re
import logging
from src.pack_graph.op_snarl import retrieve_snarl_from_minigraph_call_population, retrieve_snarl_from_minigraph_call_single, retrieve_snarl_from_vg_decomposed, retrieve_snarl_from_seqwave_call


class Node:
    def __init__(self, id, sequence):
        self.id = id
        self.sequence = sequence

        self.source = None   # # fomrat as {"id1": [start, end], "id2: [start, end]}

        self.length = len(sequence)

    def set_node_source(self, new_source_ref, new_source_start, new_source_end, overwrite=False):

        if self.source is None:
            self.source = [new_source_ref, new_source_start, new_source_end]
        else:
            if overwrite:
                self.source = [new_source_ref, new_source_start, new_source_end]

    def reset_length(self, new_length):
        self.length = new_length

    def reset_seq(self, seq):
        self.sequence = seq


class GFA:

    def __init__(self):

        # self.gfa_path = gfa_path

        self.nodes = {}
        self.small_nodes = {}
        # self.edges = {}
        # self.paths = {}

        self.snarls = {}

        self.gfa2fa_path = ""

    def set_gfa2fa_path(self, gfa2fa_path):
        self.gfa2fa_path = gfa2fa_path

    def release_memo(self):
        self.small_nodes = {}

    def parse_gfa_file_small_nodes(self, gfa_path, min_sv_size, spec_snarls=None):

        if spec_snarls is not None:
            available_nodes = []

            for snarl_id in spec_snarls:
                if snarl_id not in self.snarls:
                    continue

                available_nodes.extend(re.findall(r'[><]([a-zA-Z0-9]+)', snarl_id))

                for path in self.snarls[snarl_id].path_asm_dict.keys():
                    available_nodes.extend(re.findall(r'[><]([a-zA-Z0-9]+)', path))
        else:
            available_nodes = "all"

        # traverse each line in gfa file
        with open(gfa_path) as fin:
            for line in fin:
                line_split = line.strip().split("\t")

                marker = line_split[0]

                # # is a node
                if marker == "S":
                    # node_id, node_seq = line_split[1], line_split[2]
                    node_id = line_split[1]

                    if available_nodes != "all" and node_id not in available_nodes:
                        continue

                    node_seq = line_split[2]
                    node_length = len(node_seq)

                    if node_length < min_sv_size:
                        node = Node(node_id, "")
                        node.reset_length(node_length)

                        self.small_nodes[node_id] = node

                else:
                    continue

    def parse_gfa_file(self, gfa_path, spec_snarls=None):
        """
        parse the gfa file to get nodes, edges and paths
        """

        if spec_snarls is not None:
            available_nodes = []

            for snarl_id in spec_snarls:
                if snarl_id not in self.snarls:
                    continue

                available_nodes.extend(re.findall(r'[><]([a-zA-Z0-9]+)', snarl_id))

                for path in self.snarls[snarl_id].path_asm_dict.keys():
                    available_nodes.extend(re.findall(r'[><]([a-zA-Z0-9]+)', path))
        else:
            available_nodes = "all"

        # traverse each line in gfa file
        with open(gfa_path) as fin:
            for line in fin:
                line_split = line.strip().split("\t")

                marker = line_split[0]

                # # is a node
                if marker == "S":
                    # node_id, node_seq = line_split[1], line_split[2]
                    node_id = line_split[1]

                    if available_nodes != "all" and node_id not in available_nodes:
                        continue

                    node_seq = line_split[2]
                    node_length = len(node_seq)
                    node = Node(node_id, node_seq)
                    node.reset_length(node_length)

                    # # STEP: get node source
                    node_source_ref = None
                    node_source_start = None
                    node_source_end = None

                    for item in line_split[3: ]:
                        if "SN:Z" in item:
                            node_source_ref = item.split(":")[-1]
                        if "SO:i" in item:
                            node_source_start = int(item.split(":")[-1])
                            node_source_end = node_source_start + node_length - 1
                    if node_source_ref is None or node_source_start is None or node_source_end is None:
                        # node.reset_seq(line_split[2])
                        pass
                    else:
                        node.set_node_source(node_source_ref, node_source_start, node_source_end)

                    self.nodes[node_id] = node

                # # is a edge
                elif marker == "L":
                    continue
                    # break   # # skip edges

                    # start_node, start_node_orient = line_split[1], line_split[2]
                    # end_node, end_node_orient = line_split[3], line_split[4]
                    #
                    # edge_id = "{}-{}".format(start_node, end_node)
                    #
                    # if edge_id not in self.edges:
                    #     self.edges[edge_id] = []
                    #
                    # self.edges[edge_id] = [start_node, start_node_orient, end_node, end_node_orient]

                # # is a path (for pggb)
                elif marker == "P":
                    continue
                    # break
                    # path_id, path = line_split[1], line_split[2]
                    #
                    # if path_id not in self.paths:
                    #     self.paths[path_id] = []
                    #
                    # self.paths[path_id].append(path)

                # # is a path (for cactus)
                elif marker == "W":
                    continue
                    # break
                    # path_id, path = line_split[1], line_split[-1]
                    #
                    # if path_id not in self.paths:
                    #     self.paths[path_id] = []
                    #
                    # self.paths[path_id].append(path)

    # def to_file(self, file_path, refs, overwrite=True):
    #
    #     if os.path.exists(file_path) and overwrite is False:
    #         print("[ERROR] File exits, please specific overwrite")
    #
    #     with open(file_path, "w") as fout:
    #
    #         # # STEP: output nodes
    #         for node_id in self.nodes:
    #
    #             current_node = self.nodes[node_id]
    #
    #             if len(current_node.source.keys()) == 0:
    #                 fout.write("S\t{}\t{}\tLN:i:{}\tSN:Z:{}\tSO:i:{}\n".format(node_id, current_node.sequence, current_node.length, None, None, None))
    #             else:
    #
    #                 # # output by the order of refs
    #                 for node_source_ref in refs:
    #                     if node_source_ref in current_node.source:
    #                         node_source_start = current_node.source[node_source_ref][0]
    #                         node_source_end = current_node.source[node_source_ref][1]
    #                         fout.write("S\t{}\t{}\tLN:i:{}\tSN:Z:{}\tSO:i:{}\tSE:i:{}\n".format(node_id, current_node.sequence, current_node.length, node_source_ref, node_source_start, node_source_end))
    #
    #                         break
    #
    #         # # STEP: output edges
    #         for edge_id in self.edges:
    #             start_node, start_node_orient, end_node, end_node_orient = self.edges[edge_id][0], self.edges[edge_id][1], self.edges[edge_id][2], self.edges[edge_id][3]
    #
    #             fout.write("L\t{}\t{}\t{}\t{}\n".format(start_node, start_node_orient, end_node, end_node_orient))
    #
    #         # # STEP: output paths
    #         for path_id in self.paths:
    #             for sub_path in self.paths[path_id]:
    #                 fout.write("P\t{}\t{}\t*\n".format(path_id, sub_path))


def load_pangraph(options):

    gfa_minigraph = GFA()

    if options.remove_small:
        gfa_minigraph.parse_gfa_file_small_nodes(options.gfa_path, options.min_sv_size)
    # gfa_minigraph.parse_gfa_file(options.gfa_path)

    if options.gfa_source in ["cactus", "pggb"]:
        retrieve_snarl_from_vg_decomposed(gfa_minigraph, options)

    elif options.gfa_source in ["minigraph"]:
        for i in range(len(options.asm_paths)):

            if options.output_mode == "population":
               retrieve_snarl_from_minigraph_call_population(options.asm_names[i], options.asm_paths[i], gfa_minigraph, options)

            elif options.output_mode == "single":
                retrieve_snarl_from_minigraph_call_single(options.asm_names[i], options.asm_paths[i], gfa_minigraph, options)

            else:
                logging.error("No this output mode {}".format(options.output_mode))
                exit()

    else:
        logging.error("No this GFA source {}".format(options.gfa_source))
        exit()

    # gfa_minigraph.parse_gfa_file(options.spec_snarls)
    gfa_minigraph.set_gfa2fa_path(gfatools_gfa2fa(options))

    # resolve_large_snarl_info_sub_snarls(gfa_minigraph, options)
    #
    # simplify_gfa(options.gfa_path)
    #
    # for snarl_id in gfa_minigraph.snarls:
    #     gfa_minigraph.snarls[snarl_id].path_list = list(gfa_minigraph.snarls[snarl_id].path_asm_dict.keys())
    #     gfa_minigraph.snarls[snarl_id].subpath_list = []

    gfa_minigraph.release_memo()

    return gfa_minigraph


def minigraph_add(options):

    gfa_minigraph = GFA()

    # # STEP: run graph increment
    # added_gfa_path = os.path.join(options.output_path, "pg.minigraph.add.gfa")
    #
    # logging.info("Preprocessing: Generate the incremental graph to {}".format(added_gfa_path))
    #
    # if not os.path.exists(added_gfa_path):
    #     cmd_str = "{} -t {} -cxggs {} {} > {} \n".format(options.minigraph, options.thread_num, options.gfa_path, " ".join(options.asm_paths), added_gfa_path)
    #     os.system(cmd_str)
    #
    # options.gfa_path = added_gfa_path

    # # STEP: minigraph call
    for i in range(len(options.asm_paths)):

        asm_call_path = os.path.join(options.output_path, "{}.vcf".format(options.asm_names[i]))     # # then search the output path

        if not os.path.exists(asm_call_path):
            logging.info("Preprocessing: Generate asm VCF file from {}".format(options.asm_paths[i]))

            # cmd_str = "{} -t {} -xasm --call {} {} > {}".format(options.minigraph, options.thread_num, added_gfa_path, options.asm_paths[i], asm_call_path)
            cmd_str = "{} -t {} -xasm --call {} {} > {}".format(options.minigraph, options.thread_num, options.gfa_path, options.asm_paths[i], asm_call_path)

            os.system(cmd_str)

        retrieve_snarl_from_minigraph_call_population(options.asm_names[i], options.asm_paths[i], gfa_minigraph, options)

    gfa_minigraph.set_gfa2fa_path(gfatools_gfa2fa(options))

    existing_snarl_records = retrieve_snarl_from_seqwave_call(options.vcf_path)

    return gfa_minigraph, existing_snarl_records


def minigraph_merge(options, raw_snarls, shared_snarl_records):

    # # load the unique_snarl_records
    unique_snarl_records = retrieve_snarl_from_seqwave_call(os.path.join(options.output_path, "seqwave.hap_level.split.vcf"), reset_population=True)

    merged_snarl_records = unique_snarl_records

    # # merge shared records into the all records
    for snarl_id in shared_snarl_records:
        snarl = shared_snarl_records[snarl_id]

        if snarl_id not in merged_snarl_records:
            merged_snarl_records[snarl_id] = {}

        for path in snarl:
            path_record = snarl[path]

            merged_snarl_records[snarl_id][path] = path_record

    # # output
    from src.pack_sv.op_sv import interpret_add_to_variant
    interpret_add_to_variant(raw_snarls, merged_snarl_records, options)

    # # # for these unique snarls that without a called record, mark they are reference
    # for snarl_id in unique_snarls:
    #
    #     snarl = unique_snarls[snarl_id]
    #
    #     if snarl_id not in unique_snarl_records:
    #         unique_snarl_records[snarl_id] = {}
    #
    #     for alt_path in snarl.path_asm_dict:
    #         if alt_path not in unique_snarl_records[snarl_id]:
    #             unique_snarl_records[snarl_id][alt_path] = "reference"

    # # # merge the unique records into existing records
    # for snarl_id in unique_snarl_records:
    #
    #     snarl = unique_snarl_records[snarl_id]
    #
    #     if snarl_id not in existing_snarl_records:
    #         existing_snarl_records[snarl_id] = {}
    #
    #     # # all the path in unique_snarl_records are new
    #     for alt_path in snarl:
    #         alt_path_record = snarl[alt_path]
    #         existing_snarl_records[snarl_id][alt_path] = alt_path_record


def gfatools_gfa2fa(options):

    gfa2fa_path_1 = options.gfa_path + ".gfa2fa.fa" # # first search the same path as gfa_path
    gfa2fa_path_2 = os.path.join(options.output_path, os.path.basename(options.gfa_path) + ".gfa2fa.fa")   # # then search

    if os.path.exists(gfa2fa_path_1):

        logging.info("Preprocessing: Load gfa2fa file from {}".format(gfa2fa_path_1))
        return gfa2fa_path_1

    elif os.path.exists(gfa2fa_path_2):
        logging.info("Preprocessing: Load gfa2fa file from {}".format(gfa2fa_path_2))

        return gfa2fa_path_2

    else:
        logging.info("Preprocessing: Generate gfa2fa file from {}".format(gfa2fa_path_2))

        cmd_str = "{} gfa2fa {} >{}".format(options.gfatools, options.gfa_path, gfa2fa_path_2)
        os.system(cmd_str)

        return gfa2fa_path_2

