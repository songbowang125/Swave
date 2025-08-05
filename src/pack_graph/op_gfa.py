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

    gfa_obj = GFA()

    if options.remove_small:
        gfa_obj.parse_gfa_file_small_nodes(options.gfa_path, options.min_sv_size)
    # gfa_obj.parse_gfa_file(options.gfa_path)

    if options.gfa_source in ["cactus", "pggb"]:
        retrieve_snarl_from_vg_decomposed(gfa_obj, options)

    elif options.gfa_source in ["minigraph"]:
        for i in range(len(options.asm_paths)):

            if options.output_mode == "population":
               retrieve_snarl_from_minigraph_call_population(options.asm_names[i], options.asm_paths[i], gfa_obj, options)

            elif options.output_mode == "single":
                retrieve_snarl_from_minigraph_call_single(options.asm_names[i], options.asm_paths[i], gfa_obj, options)

            else:
                logging.error("No this output mode {}".format(options.output_mode))
                exit()

    else:
        logging.error("No this GFA source {}".format(options.gfa_source))
        exit()

    # gfa_obj.parse_gfa_file(options.spec_snarls)
    gfa_obj.set_gfa2fa_path(gfatools_gfa2fa(options))

    # resolve_large_snarl_info_sub_snarls(gfa_obj, options)
    #
    # simplify_gfa(options.gfa_path)
    #
    # for snarl_id in gfa_obj.snarls:
    #     gfa_obj.snarls[snarl_id].path_list = list(gfa_obj.snarls[snarl_id].path_asm_dict.keys())
    #     gfa_obj.snarls[snarl_id].subpath_list = []

    gfa_obj.release_memo()

    return gfa_obj


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

