import os
import pysam
import re
import logging


class Snarl:

    def __init__(self, start_node_id, start_node_orient, end_node_id, end_node_orient, ref_chrom, ref_start, ref_end, reversed_mapping=False):

        self.snarl_id = "{}{}{}{}".format(start_node_orient, start_node_id, end_node_orient, end_node_id)

        self.start_node_id = start_node_id
        self.start_node_orient = start_node_orient

        self.end_node_id = end_node_id
        self.end_node_orient = end_node_orient

        self.path_asm_dict = {}

        self.ref_chrom = ref_chrom
        self.ref_start = ref_start
        self.ref_end = ref_end

        self.available_nodes = []

        self.reversed_mapping = reversed_mapping
        # self.sample_gts = sample_gts

        # self.include_CSVs = [None for i in range(len(self.include_paths))]

    def get_available_nodes(self):
        self.available_nodes.extend(re.findall(r'[><]([a-zA-Z0-9]+)', self.snarl_id))
        for path in self.path_asm_dict.keys():
            self.available_nodes.extend(re.findall(r'[><]([a-zA-Z0-9]+)', path))

        return self.available_nodes

    def set_ref_asm_path(self, ref_asm_path):
        self.ref_asm_path = ref_asm_path

    def to_string(self):
        return "{}, {}-{}-{}, {}".format(self.snarl_id, self.ref_chrom, self.ref_start, self.ref_end, self.ref_end - self.ref_start)


def simplify_gfa(input_gfa):
    """
    remove the sequence of nodes in the gfa. This will reduce the size of gfa file
    """

    output_gfa = input_gfa.replace(".gfa", ".simplified.gfa")

    with open(input_gfa) as fin, open(output_gfa, "w") as fout:

        for line in fin:

            if line.startswith("S"):

                line_split = line.strip().split("\t")

                fout.write("{}\tN\t{}\n".format("\t".join(line_split[0: 2]), "\t".join(line_split[3:])))

            else:
                fout.write(line)


def compare_seq_source(base_source, target_source):
    base_source_split = base_source.split(":")
    base_source_chr = base_source_split[0]
    base_source_start = int(base_source_split[1])
    base_source_end = int(base_source_split[2])

    target_source_split = target_source.split(":")
    target_source_chr = target_source_split[0]
    target_source_start = int(target_source_split[1])
    target_source_end = int(target_source_split[2])

    if base_source_chr != target_source_chr:
        return 0.0, None, None, None

    overlap_start = max(base_source_start, target_source_start)
    overlap_end = min(base_source_end, target_source_end)

    if overlap_start < overlap_end:
        new_source_chr = base_source_chr
        new_source_start = min(base_source_start, target_source_start)
        new_source_end = max(base_source_end, target_source_end)

        return (overlap_end - overlap_start) / min(base_source_end - base_source_start + 1, target_source_end - target_source_start + 1), new_source_chr, new_source_start, new_source_end

    if overlap_start > overlap_end and overlap_start - overlap_end < 200:
        new_source_chr = base_source_chr
        new_source_start = min(base_source_start, target_source_start)
        new_source_end = max(base_source_end, target_source_end)

        # return abs(overlap_end - overlap_start) / min(base_source_end - base_source_start + 1, target_source_end - target_source_start + 1), new_source_chr, new_source_start, new_source_end
        return abs(overlap_end - overlap_start) / min(abs(base_source_end - base_source_start) + 1, abs(target_source_end - target_source_start) + 1), new_source_chr, new_source_start, new_source_end

    return 0.0, None, None, None


def remove_snarl_nodes_from_path(raw_path):
    raw_path_include_nodes = re.findall(r'[><]([a-zA-Z0-9]+)', raw_path)
    raw_path_include_nodes_orients = re.findall(r'[><]', raw_path)

    new_path = "".join("{}{}".format(raw_path_include_nodes_orients[i], raw_path_include_nodes[i]) for i in range(1, len(raw_path_include_nodes) - 1))
    if new_path == "":
        new_path = "*"

    return new_path


def retrieve_snarl_from_vg_decomposed(gfa_minigraph, options):

    logging.info("Preprocessing: Load allele info from {}".format(options.decomposed_vcf))

    with pysam.VariantFile(options.decomposed_vcf) as fin:

        sample_order = str(fin.header).strip().split("\n")[-1].split("\t")[9: ]

        for record in fin:

            if not (len(record.ref) >= 50 or max([len(alt_seq) for alt_seq in record.alts]) >= 50):
                continue

            snarl_id = record.id
            snarl_ref_chrom = record.contig
            snarl_ref_start = record.start + 1
            snarl_ref_end = record.stop + 1

            if options.spec_snarl is not None and snarl_id not in options.spec_snarl:
                continue

            snarl_id_include_nodes = re.findall(r'[><]([a-zA-Z0-9]+)', snarl_id)
            snarl_id_include_nodes_orients = re.findall(r'[><]', snarl_id)

            snarl_start_node_id = snarl_id_include_nodes[0]
            snarl_start_node_orient = snarl_id_include_nodes_orients[0]

            snarl_end_node_id = snarl_id_include_nodes[1]
            snarl_end_node_orient = snarl_id_include_nodes_orients[1]

            if not (snarl_start_node_orient == ">" and snarl_end_node_orient == ">"):
                continue

            if snarl_id not in gfa_minigraph.snarls:
                gfa_minigraph.snarls[snarl_id] = Snarl(snarl_start_node_id, snarl_start_node_orient, snarl_end_node_id, snarl_end_node_orient, snarl_ref_chrom, snarl_ref_start, snarl_ref_end)

            snarl_obj = gfa_minigraph.snarls[snarl_id]

            # # obtain all paths form AT, since the path includes the snarl start and end node, we need to remove them
            all_paths = [remove_snarl_nodes_from_path(path) for path in record.info["AT"]]

            # # set ref asm path
            ref_asm_path = all_paths[0]
            if ref_asm_path not in snarl_obj.path_asm_dict:
                snarl_obj.path_asm_dict[ref_asm_path] = []
            snarl_obj.path_asm_dict[ref_asm_path].append(options.ref_asm_name)

            snarl_obj.set_ref_asm_path(ref_asm_path)

            # # set alt asm path by parsing sample gt flags and find the current path
            sample_gts = str(record).strip().split("\t")[9: ]
            for sample_index in range(len(sample_gts)):
                sample_gt = sample_gts[sample_index]
                sample_name = sample_order[sample_index]

                if sample_name not in options.input_dict:
                    continue

                hap_gts = sample_gt.split("|")
                for hap_index in range(len(hap_gts)):
                    hap_gt = hap_gts[hap_index]

                    if hap_gt == ".":
                        continue

                    alt_path = all_paths[int(hap_gt)]
                    hap_name = os.path.basename(options.input_dict[sample_name][hap_index])

                    if alt_path not in snarl_obj.path_asm_dict:
                        snarl_obj.path_asm_dict[alt_path] = []
                    snarl_obj.path_asm_dict[alt_path].append(hap_name)


def get_reverse_path_from_node_list(node_list):

    reversed_path = ""

    for i in range(len(node_list) - 1, -1, -1):

        node_with_orient = node_list[i]

        if ">" in node_with_orient:
            reversed_path += node_with_orient.replace(">", "<")
        else:
            reversed_path += node_with_orient.replace("<", ">")

    return reversed_path


def retrieve_reverse_mapping_snarl(contig_mapping_dict, contig_mapping_stats):

    reverse_mapping_snarl_dict = {}

    for cur_contig in contig_mapping_stats:

        cur_contig_stats = contig_mapping_stats[cur_contig]
        cur_contig_list = contig_mapping_dict[cur_contig]

        cur_contig_strand = None
        if (cur_contig_stats['+'] - cur_contig_stats['-']) > 0.8 * cur_contig_stats['+']:
            cur_contig_strand = "+"

        if (cur_contig_stats['-'] - cur_contig_stats['+']) > 0.8 * cur_contig_stats['+']:
            cur_contig_strand = "-"

        if cur_contig_strand is None:
            continue

        previous_reverse_index = None
        previous_reverse_nodes = []
        previous_reverse_cords = []
        previous_reverse_chrom = None

        for sub_index in range(len(cur_contig_list)):
            sub_mapping = cur_contig_list[sub_index]
            sub_mapping_strand, sub_start_node_with_orient, sub_end_node_orient, sub_alt_path, ref_chrom, ref_start, ref_end = sub_mapping[0], sub_mapping[1], sub_mapping[2], sub_mapping[3], sub_mapping[4], int(sub_mapping[5]), int(sub_mapping[6])

            if sub_mapping_strand == cur_contig_strand or sub_mapping_strand == ".":
                if previous_reverse_index is not None:
                    reverse_ref_start, reverse_ref_end = min(previous_reverse_cords), max(previous_reverse_cords)
                    reverse_start_node, reverse_end_node = previous_reverse_nodes[0], previous_reverse_nodes[-1]
                    if reverse_ref_end - reverse_ref_start > 5000:
                        reverse_mapping_snarl_dict["{}{}".format(reverse_start_node, reverse_end_node)] = [cur_contig, get_reverse_path_from_node_list(previous_reverse_nodes), previous_reverse_chrom, reverse_ref_start, reverse_ref_end]

                previous_reverse_index = None
                previous_reverse_nodes = []
                previous_reverse_cords = []
                previous_reverse_chrom = None


            else:

                if previous_reverse_index is None:
                    previous_reverse_index = sub_index

                previous_reverse_chrom = ref_chrom

                previous_reverse_cords.append(ref_start)
                previous_reverse_cords.append(ref_end)

                # # append nodes from snarl and alt path
                if sub_start_node_with_orient not in previous_reverse_nodes:
                    previous_reverse_nodes.append(sub_start_node_with_orient)

                alt_path_include_nodes = re.findall(r'[><]([a-zA-Z0-9]+)', sub_alt_path)
                alt_path_include_nodes_orients = re.findall(r'[><]', sub_alt_path)
                for i in range(len(alt_path_include_nodes)):
                    previous_reverse_nodes.append("{}{}".format(alt_path_include_nodes_orients[i], alt_path_include_nodes[i]))

                if sub_end_node_orient not in previous_reverse_nodes:
                    previous_reverse_nodes.append(sub_end_node_orient)

    return reverse_mapping_snarl_dict


def retrieve_snarl_from_minigraph_call_population(spec_asm_name, spec_asm_path, gfa_minigraph, options):
    """

    One step,
    traverse the vcf file and fetch the raw snarl information
    """

    # # find or generate the
    spec_asm_vcf_path_1 = os.path.join(os.path.dirname(options.gfa_path), "{}.bed".format(spec_asm_name)) # # fist search the gfa file's path
    spec_asm_vcf_path_2 = os.path.join(options.output_path, "{}.bed".format(spec_asm_name))     # # then search the output path

    if os.path.exists(spec_asm_vcf_path_1):

        logging.info("Preprocessing: Load allele info from {}".format(spec_asm_vcf_path_1))
        spec_asm_vcf_path = spec_asm_vcf_path_1

    elif os.path.exists(spec_asm_vcf_path_2):

        logging.info("Preprocessing: Load allele info from {}".format(spec_asm_vcf_path_2))
        spec_asm_vcf_path = spec_asm_vcf_path_2

    else:

        logging.info("Preprocessing: Generate allele info from {}".format(spec_asm_vcf_path_2))

        os.system("{} -t {} -xasm --call {} {} > {}".format(options.minigraph, options.thread_num, options.gfa_path, spec_asm_path, spec_asm_vcf_path_2))

        spec_asm_vcf_path = spec_asm_vcf_path_2

    contig_mapping_dict = {}
    contig_mapping_stats = {}

    # # STEP: traverse the vcf file and fetch the raw snarl information
    previous_seq_source_contig = None
    with open(spec_asm_vcf_path) as fin:
        for line in fin:
            line_split = line.strip().split("\t")

            snarl_ref_chrom, snarl_ref_start, snarl_ref_end = line_split[0], int(line_split[1]), int(line_split[2])

            snarl_start_node_with_orient = line_split[3]
            snarl_start_node_id = snarl_start_node_with_orient[1:]
            snarl_start_node_orient = snarl_start_node_with_orient[0]

            snarl_end_node_with_orient = line_split[4]
            snarl_end_node_id = snarl_end_node_with_orient[1:]
            snarl_end_node_orient = snarl_end_node_with_orient[0]

            snarl_id = "{}{}{}{}".format(snarl_start_node_orient, snarl_start_node_id, snarl_end_node_orient, snarl_end_node_id, )

            if options.spec_snarl is not None and snarl_id != options.spec_snarl:
                continue

            alt_path_split = line_split[5].split(":")
            alt_path = alt_path_split[0]

            if options.force_reverse:
                # # for saving inverted mapping
                if alt_path == ".":
                    seq_source_contig, seq_source_start, seq_source_end, seq_source_strand = previous_seq_source_contig, 0, 0, "."
                else:
                    seq_source_contig, seq_source_start, seq_source_end, seq_source_strand = alt_path_split[3], int(alt_path_split[4]), int(alt_path_split[5]), alt_path_split[2]

                if seq_source_contig not in contig_mapping_dict:
                    contig_mapping_dict[seq_source_contig] = []
                    contig_mapping_stats[seq_source_contig] = {"+": 0, "-": 0, ".": 0}

                contig_mapping_dict[seq_source_contig].append([seq_source_strand, snarl_start_node_with_orient, snarl_end_node_with_orient, alt_path, snarl_ref_chrom, snarl_ref_start, snarl_ref_end])
                contig_mapping_stats[seq_source_contig][seq_source_strand] += (seq_source_end - seq_source_start)

                previous_seq_source_contig = seq_source_contig

                # #

            if alt_path == ".":  # # why alt_path is '.': because the asm has no contig covering this snarl
                continue

            if options.remove_small:
                # # STEP: remove small nodes with length less than min_sv_size
                alt_path_include_nodes = re.findall(r's\d+', alt_path)
                alt_path_include_nodes_orients = re.findall(r'([^a-zA-Z0-9]+)', alt_path)

                small_node_indexs = []
                for index in range(len(alt_path_include_nodes) - 1, -1, -1):
                    if alt_path_include_nodes[index] in gfa_minigraph.small_nodes:
                    # if gfa_minigraph.small_nodes[alt_path_include_nodes[index]].length < options.min_sv_size:
                        small_node_indexs.append(index)

                for index in small_node_indexs:
                    alt_path_include_nodes.pop(index)
                    alt_path_include_nodes_orients.pop(index)

                alt_path = "".join("{}{}".format(alt_path_include_nodes_orients[i], alt_path_include_nodes[i]) for i in range(len(alt_path_include_nodes)))
                # # done

            if not options.force_reverse:
                if alt_path == "":
                    alt_path = "*"

                if snarl_id not in gfa_minigraph.snarls:
                    gfa_minigraph.snarls[snarl_id] = Snarl(snarl_start_node_id, snarl_start_node_orient, snarl_end_node_id, snarl_end_node_orient, snarl_ref_chrom, snarl_ref_start, snarl_ref_end)

                snarl_obj = gfa_minigraph.snarls[snarl_id]

                # # STEP: add the alt path to gfa's snarls
                if alt_path not in snarl_obj.path_asm_dict:
                    snarl_obj.path_asm_dict[alt_path] = []
                snarl_obj.path_asm_dict[alt_path].append(spec_asm_name)

                if spec_asm_name == options.ref_asm_name:
                    snarl_obj.set_ref_asm_path(alt_path)

    if options.force_reverse:
        # # deal with reversed mapping
        reverse_mapping_snarl_dict = retrieve_reverse_mapping_snarl(contig_mapping_dict, contig_mapping_stats)

        for snarl_id in reverse_mapping_snarl_dict:
            snarl_info = reverse_mapping_snarl_dict[snarl_id]
            # print(snarl_id, snarl_info)

            alt_path, ref_chrom, ref_start, ref_end = snarl_info[1], snarl_info[2], int(snarl_info[3]), int(snarl_info[4])

            if alt_path == "":
                alt_path = "*"

            snarl_id_include_nodes = re.findall(r'[><]([a-zA-Z0-9]+)', snarl_id)
            snarl_id_include_nodes_orients = re.findall(r'[><]', snarl_id)

            snarl_start_node_id, snarl_start_node_orient = snarl_id_include_nodes[0], snarl_id_include_nodes_orients[0]
            snarl_end_node_id, snarl_end_node_orient = snarl_id_include_nodes[1], snarl_id_include_nodes_orients[1]

            if snarl_id not in gfa_minigraph.snarls:
                gfa_minigraph.snarls[snarl_id] = Snarl(snarl_start_node_id, snarl_start_node_orient, snarl_end_node_id, snarl_end_node_orient, ref_chrom, ref_start, ref_end, reversed_mapping=True)

            snarl_obj = gfa_minigraph.snarls[snarl_id]

            # # STEP: add the alt path to gfa's snarls
            if alt_path not in snarl_obj.path_asm_dict:
                snarl_obj.path_asm_dict[alt_path] = []
            snarl_obj.path_asm_dict[alt_path].append(spec_asm_name)

            if spec_asm_name == options.ref_asm_name:
                snarl_obj.set_ref_asm_path(alt_path)


def retrieve_snarl_from_minigraph_call_single(spec_asm_name, spec_asm_path, gfa_minigraph, options):
    """

    Two steps,
    first, traverse the vcf file and fetch the raw snarl information
    second, calculate the final snarl information
    """

    # # find or generate the
    spec_asm_vcf_path_1 = os.path.join(os.path.dirname(options.gfa_path), "{}.bed".format(spec_asm_name)) # # fist search the gfa file's path
    spec_asm_vcf_path_2 = os.path.join(options.output_path, "{}.bed".format(spec_asm_name))     # # then search the output path

    if os.path.exists(spec_asm_vcf_path_1):

        logging.info("Loading allele info from {}".format(spec_asm_vcf_path_1))
        spec_asm_vcf_path = spec_asm_vcf_path_1

    elif os.path.exists(spec_asm_vcf_path_2):

        logging.info("Loading allele info from {}".format(spec_asm_vcf_path_2))
        spec_asm_vcf_path = spec_asm_vcf_path_2

    else:

        logging.info("Generating allele info from {}".format(spec_asm_vcf_path_2))

        os.system("{} -t {} -xasm --call {} {} > {}".format(options.minigraph, options.thread_num, options.gfa_path, spec_asm_path, spec_asm_vcf_path_2))

        spec_asm_vcf_path = spec_asm_vcf_path_2

    # # STEP: traverse the vcf file and fetch the raw snarl information
    source_snarl_dict = {}
    with open(spec_asm_vcf_path) as fin:

        for line in fin:
            line_split = line.strip().split("\t")

            snarl_ref_chrom, snarl_ref_start, snarl_ref_end = line_split[0], int(line_split[1]), int(line_split[2])

            snarl_start_node, snarl_end_node = line_split[3], line_split[4]

            alt_path = line_split[5].split(":")[0]

            seq_source = ":".join(line_split[5].split(":")[3:6])

            if alt_path == ".":  # # why alt_path is '.': because the asm has no contig covering this snarl
                continue

            # # STEP: remove small nodes with length less than min_sv_size
            # alt_path_include_nodes = re.findall(r's\d+', alt_path)
            # alt_path_include_nodes_orients = re.findall(r'([^a-zA-Z0-9]+)', alt_path)
            #
            # small_node_indexs = []
            # for index in range(len(alt_path_include_nodes) - 1, -1, -1):
            #
            #     if gfa_minigraph.nodes[alt_path_include_nodes[index]].length < options.min_sv_size:
            #         small_node_indexs.append(index)
            #
            # for index in small_node_indexs:
            #     alt_path_include_nodes.pop(index)
            #     alt_path_include_nodes_orients.pop(index)
            #
            # alt_path = "".join("{}{}".format(alt_path_include_nodes_orients[i], alt_path_include_nodes[i]) for i in range(len(alt_path_include_nodes)))
            # # done

            whole_path = "{}{}{}".format(snarl_start_node, alt_path, snarl_end_node)

            # # STEP: compare with the previous snarls to check if they are overlaped at the seq_source (the position on contig)
            if len(source_snarl_dict) == 0:
                source_snarl_dict[seq_source] = {"snarl_info": [snarl_ref_chrom, snarl_ref_start, snarl_ref_end], "include_paths": [whole_path]}
            else:

                previous_seq_source = list(source_snarl_dict.keys())[-1]

                overlap_ratio, new_source_chr, new_source_start, new_source_end = compare_seq_source(seq_source, previous_seq_source)

                # # overlapped with the previous one
                if overlap_ratio > 0.0:
                    new_seq_source = "{}:{}:{}".format(new_source_chr, new_source_start, new_source_end)

                    # # update the snarl info
                    new_snarl_ref_chrom = snarl_ref_chrom
                    new_snarl_ref_start = min(source_snarl_dict[previous_seq_source]['snarl_info'][1], snarl_ref_start)
                    new_snarl_ref_end = max(source_snarl_dict[previous_seq_source]['snarl_info'][2], snarl_ref_end)

                    # # update the dict
                    source_snarl_dict.update({new_seq_source: source_snarl_dict.pop(previous_seq_source)})

                    source_snarl_dict[new_seq_source]["snarl_info"] = [new_snarl_ref_chrom, new_snarl_ref_start, new_snarl_ref_end]
                    source_snarl_dict[new_seq_source]["include_paths"].append(whole_path)

                else:
                    source_snarl_dict[seq_source] = {"snarl_info": [snarl_ref_chrom, snarl_ref_start, snarl_ref_end], "include_paths": [whole_path]}

    # # after traverse, add this
    for seq_source in source_snarl_dict:

        include_paths = source_snarl_dict[seq_source]["include_paths"]

        # # for one specific seq_source, if more than one include path, means it's overlapped by multi-snarls
        # # then we need determine the final path
        if len(include_paths) > 1:
            joint_path = "".join(include_paths).replace("*", "")

            # joint_path_include_nodes = re.findall(r's\d+', joint_path)
            # joint_path_include_nodes_orients = re.findall(r'([^a-zA-Z0-9]+)', joint_path)

            joint_path_include_nodes = re.findall(r'[><]([a-zA-Z0-9]+)', joint_path)
            joint_path_include_nodes_orients = re.findall(r'[><]', joint_path)

            # # remove the duplicated node in joint path
            joint_path = ""
            for i in range(len(joint_path_include_nodes)):
                cur_node, cur_node_orient = joint_path_include_nodes[i], joint_path_include_nodes_orients[i]

                if cur_node not in joint_path:
                    joint_path += "{}{}".format(cur_node_orient, cur_node)

            final_path = joint_path
        else:
            final_path = include_paths[0].replace("*", "")

        # # STEP: convert from the final path, to the final snarl information
        # final_path_include_nodes = re.findall(r's\d+', final_path)
        # final_path_include_nodes_orients = re.findall(r'([^a-zA-Z0-9]+)', final_path)

        final_path_include_nodes = re.findall(r'[><]([a-zA-Z0-9]+)', final_path)
        final_path_include_nodes_orients = re.findall(r'[><]', final_path)

        start_node_id, start_node_orient = final_path_include_nodes[0], final_path_include_nodes_orients[0]
        end_node_id, end_node_orient = final_path_include_nodes[-1], final_path_include_nodes_orients[-1]

        alt_path = "".join(["{}{}".format(final_path_include_nodes_orients[i], final_path_include_nodes[i], ) for i in range(1, len(final_path_include_nodes) - 1)])
        if alt_path == "":
            alt_path = "*"

        snarl_id = "{}{}{}{}".format(start_node_orient, start_node_id, end_node_orient, end_node_id)

        if snarl_id not in gfa_minigraph.snarls:
            snarl_info = source_snarl_dict[seq_source]['snarl_info']
            snarl_ref_chrom, snarl_ref_start, snarl_ref_end = snarl_info[0], snarl_info[1], snarl_info[2]

            gfa_minigraph.snarls[snarl_id] = Snarl(start_node_id, start_node_orient, end_node_id, end_node_orient, snarl_ref_chrom, snarl_ref_start, snarl_ref_end)

        snarl_obj = gfa_minigraph.snarls[snarl_id]

        # # STEP: add the alt path to gfa's snarls
        if alt_path not in snarl_obj.path_asm_dict:
            snarl_obj.path_asm_dict[alt_path] = []
        snarl_obj.path_asm_dict[alt_path].append(spec_asm_name)

        if spec_asm_name == options.ref_asm_name:
            snarl_obj.set_ref_asm_path(alt_path)


# def filter_massive_mismatch_snarl(gfa, options):
#
#     # # STEP: traverse each snarl
#     for snarl_id in gfa.snarls:
#         snarl = gfa.snarls[snarl_id]
#
#         # # STEP: collect the ref asm path and alt asm paths
#         ref_asm_path = ""
#         alt_asm_paths = []
#
#         for path in snarl.path_asm_dict:
#             # # find the path for ref asm
#             if options.ref_asm_name in snarl.path_asm_dict[path]:
#                 ref_asm_path = path
#             else:
#                 alt_asm_paths.append(path)
#
#         # # STEP: compare the similarity between ref_asm and each alt_asm
#         for alt_asm_path in alt_asm_paths:


# def match_ref_and_alt_path_for_sub_snarls(entire_ref_path, entire_alt_path):
#
#
#     entire_ref_path = entire_ref_path.replace("*", "")
#     entire_alt_path = entire_alt_path.replace("*", "")
#
#     # if "s9>" in entire_ref_path:
#     #     print()
#     # print()
#     # print()
#     # print("=====Old=====", entire_ref_path)
#     # print("=====Old=====", entire_alt_path)
#
#     sub_snarl_dict = {}  # # format {snarl_id: [start_node_id, start_node_orient, end_node_id, end_node_orient, ref_path, alt_path}
#     #
#     # ref_nodes = re.findall(r's\d+', entire_ref_path)
#     # ref_nodes_orients = re.findall(r'([^a-zA-Z0-9]+)', entire_ref_path)
#     #
#     ref_nodes = re.findall(r'[><]([a-zA-Z0-9]+)', entire_ref_path)
#     ref_nodes_orients = re.findall(r'[><]', entire_ref_path)
#     #
#     # alt_nodes = re.findall(r's\d+', entire_alt_path)
#     # alt_nodes_orients = re.findall(r'([^a-zA-Z0-9]+)', entire_alt_path)
#
#     alt_nodes = re.findall(r'[><]([a-zA-Z0-9]+)', entire_alt_path)
#     alt_nodes_orients = re.findall(r'[><]', entire_alt_path)
#
#     # # STEP: for each ref node, find its index in alt path
#     all_ref_node_indexes_in_alt_path = []
#
#     for ref_node in ref_nodes:
#         try:
#             ref_node_index = alt_nodes.index(ref_node)
#         except ValueError:
#             ref_node_index = -1
#
#         all_ref_node_indexes_in_alt_path.append(ref_node_index)
#
#
#     # # STEP: find candidate snarl
#     snarl_start_node = ref_nodes[0]
#     snarl_start_node_orient = ref_nodes_orients[0]
#     snarl_start_node_index_in_ref_path = 0
#     snarl_start_node_index_in_alt_path = all_ref_node_indexes_in_alt_path[snarl_start_node_index_in_ref_path]
#
#     next_ref_node_index_in_ref_path = 1
#     while next_ref_node_index_in_ref_path < len(ref_nodes):
#
#         next_ref_node_index_in_alt_path = all_ref_node_indexes_in_alt_path[next_ref_node_index_in_ref_path]
#
#         # # cur ref node is missed from alt_path, then move to next
#         if next_ref_node_index_in_alt_path == -1:
#             next_ref_node_index_in_ref_path += 1
#
#         else:
#             snarl_end_node = ref_nodes[next_ref_node_index_in_ref_path]
#             snarl_end_node_orient = ref_nodes_orients[next_ref_node_index_in_ref_path]
#             snarl_end_node_index_in_ref_path = next_ref_node_index_in_ref_path
#             snarl_end_node_index_in_alt_path = all_ref_node_indexes_in_alt_path[snarl_end_node_index_in_ref_path]
#
#             # # compare the path between ref and alt
#             ref_path = "".join(["{}{}".format(ref_nodes_orients[i], ref_nodes[i]) for i in range(snarl_start_node_index_in_ref_path + 1, snarl_end_node_index_in_ref_path)])
#             alt_path = "".join(["{}{}".format(alt_nodes_orients[i], alt_nodes[i]) for i in range(snarl_start_node_index_in_alt_path + 1, snarl_end_node_index_in_alt_path)])
#
#             if ref_path == "":
#                 ref_path = "*"
#             if alt_path == "":
#                 alt_path = "*"
#
#             # # path is different, then we find a sub snarl
#             if ref_path != alt_path:
#                 snarl_id = "{}{}{}{}".format(snarl_start_node_orient, snarl_start_node, snarl_end_node_orient, snarl_end_node)
#
#                 sub_snarl_dict[snarl_id] = [snarl_start_node, snarl_start_node_orient, snarl_end_node, snarl_end_node_orient, ref_path, alt_path]
#
#             # # update info and switch to next
#             snarl_start_node = snarl_end_node
#             snarl_start_node_orient = snarl_end_node_orient
#             snarl_start_node_index_in_ref_path = snarl_end_node_index_in_ref_path
#             snarl_start_node_index_in_alt_path = snarl_end_node_index_in_alt_path
#
#             next_ref_node_index_in_ref_path += 1
#
#     # if len(sub_snarl_dict) == 0:
#     #     print("Unresolved")
#     # else:
#     #     print("Resolved")
#     #
#     # for id in sub_snarl_dict:
#     #     print("=====New=====", id, sub_snarl_dict[id])
#     # print()
#
#     return sub_snarl_dict


# def resolve_large_snarl_info_sub_snarls(gfa, options):
#
#     resolved_sub_snarls = {}
#
#     for large_snarl_id in gfa.snarls:
#
#         large_snarl = gfa.snarls[large_snarl_id]
#
#         # # a small one
#         if large_snarl.ref_end - large_snarl.ref_start < options.max_sv_size:
#             # print(11111, large_snarl.to_string())
#             resolved_sub_snarls[large_snarl_id] = gfa.snarls[large_snarl_id]
#             continue
#
#         # print("Processing\n", large_snarl.to_string())
#
#         # # STEP: find ref_asm path
#         ref_asm_path = ""
#         ref_asm_path_same_asms = [] # # collection of the asms that have the same path with ref asm
#
#         for path in large_snarl.path_asm_dict:
#             # # find the path for ref asm
#             if options.ref_asm_name in large_snarl.path_asm_dict[path]:
#                 ref_asm_path = path
#                 ref_asm_path_same_asms = large_snarl.path_asm_dict[path]
#
#         # # STEP: compare each alt path with ref path
#         for alt_asm_path in large_snarl.path_asm_dict:
#
#             if alt_asm_path == ref_asm_path:
#                 continue
#
#             alt_asm_path_same_asms = large_snarl.path_asm_dict[alt_asm_path]  # # collection of the asms that have the same alt path
#
#             entire_ref_asm_path = "{}{}".format(large_snarl.start_node_orient, large_snarl.start_node_id) + ref_asm_path + "{}{}".format(large_snarl.end_node_orient, large_snarl.end_node_id)
#             entire_alt_asm_path = "{}{}".format(large_snarl.start_node_orient, large_snarl.start_node_id) + alt_asm_path + "{}{}".format(large_snarl.end_node_orient, large_snarl.end_node_id)
#
#             sub_snarl_dict = match_ref_and_alt_path_for_sub_snarls(entire_ref_asm_path, entire_alt_asm_path)    # # format {snarl_id: [sub_snarl_start_node_id, sub_snarl_start_node_orient, sub_snarl_end_node_id, sub_snarl_end_node_orient, ref_path, alt_path}
#
#             for sub_snarl_id in sub_snarl_dict:
#                 sub_snarl_info = sub_snarl_dict[sub_snarl_id]
#
#                 sub_snarl_ref_path = sub_snarl_info[4]
#                 sub_snarl_alt_path = sub_snarl_info[5]
#
#                 first_time_exist = False
#                 if sub_snarl_id not in resolved_sub_snarls:
#                     first_time_exist = True
#
#                     # # STEP: collect the info of this snarl
#                     sub_snarl_start_node_id, sub_snarl_start_node_orient, sub_snarl_end_node_id, sub_snarl_end_node_orient = sub_snarl_info[0], sub_snarl_info[1], sub_snarl_info[2], sub_snarl_info[3]
#
#                     if sub_snarl_start_node_id not in gfa.nodes or sub_snarl_end_node_id not in gfa.nodes:
#                         continue
#
#                     sub_snarl_ref_chrom, sub_snarl_ref_start, sub_snarl_ref_end = large_snarl.ref_chrom, gfa.nodes[sub_snarl_start_node_id].source[2] + 1, gfa.nodes[sub_snarl_end_node_id].source[1]
#
#                     resolved_sub_snarls[sub_snarl_id] = Snarl(sub_snarl_start_node_id, sub_snarl_start_node_orient, sub_snarl_end_node_id, sub_snarl_end_node_orient, sub_snarl_ref_chrom, sub_snarl_ref_start, sub_snarl_ref_end)
#
#                 sub_snarl_obj = resolved_sub_snarls[sub_snarl_id]
#
#                 # # STEP: add the ref and alt path to gfa's snarls
#                 if sub_snarl_ref_path not in sub_snarl_obj.path_asm_dict:
#                     sub_snarl_obj.path_asm_dict[sub_snarl_ref_path] = ref_asm_path_same_asms
#
#                 if sub_snarl_alt_path not in sub_snarl_obj.path_asm_dict:
#                     sub_snarl_obj.path_asm_dict[sub_snarl_alt_path] = []
#                 sub_snarl_obj.path_asm_dict[sub_snarl_alt_path].extend(alt_asm_path_same_asms)
#
#                 sub_snarl_obj.set_ref_asm_path(sub_snarl_ref_path)
#
#                 print("=====Snarl_{}=====".format(first_time_exist), sub_snarl_obj.to_string(), sub_snarl_obj.path_asm_dict)
#
#     # for snarl_id in resolved_sub_snarls:
#     #
#     #     snarl = resolved_sub_snarls[snarl_id]
#         # print("Processing", snarl.to_string(), snarl.path_asm_dict)
#
#     gfa.snarls = resolved_sub_snarls
#
#     # print(1111111111111, gfa.snarls[">s44>s46"].to_string(), gfa.snarls[">s44>s46"].path_asm_dict)


def retrieve_snarl_from_seqwave_call(vcf_path, reset_population=False):
    """
    load the already-called snarls from seqwave output
    """

    snarl_records = {}       # # format {'snarl_id': {path1: record1, path2: record2}}

    with pysam.VariantFile(vcf_path) as vcf_file:

        for record in vcf_file:

            snarl_id, snarl_ref_path, snarl_alt_path = record.id.split("_")[0], record.ref, record.alts[0]

            if snarl_id not in snarl_records:
                snarl_records[snarl_id] = {}

            if snarl_ref_path not in snarl_records[snarl_id]:
                snarl_records[snarl_id][snarl_ref_path] = "reference"

            if snarl_alt_path not in snarl_records[snarl_id]:
                if reset_population:
                    record_split = str(record).strip().split("\t")
                    info_field = record_split[7].split(";")[: -4]   # # remove the AC AF AN NS info
                    info_field.append("AC=0;AF=0;AN=.;NS=.")

                    record_split[7] = ";".join(info_field)

                    snarl_records[snarl_id][snarl_alt_path] = "\t".join(record_split)

                else:
                    snarl_records[snarl_id][snarl_alt_path] = str(record).strip()

    return snarl_records




