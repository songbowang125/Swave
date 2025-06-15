import os
from src.pack_graph.op_gfa import GFA
import pysam
import re
from src.pack_graph.op_seq import reverse_complement_seq


def convert_path_to_seq(path, nodes):

    if path == "*":
        return ""

    path_include_nodes = re.findall(r'[><]([a-zA-Z0-9]+)', path)
    path_include_nodes_orients = re.findall(r'[><]', path)

    path_seq = ""
    for i in range(len(path_include_nodes)):

        node, node_orient = path_include_nodes[i], path_include_nodes_orients[i]

        try:
            if node_orient == ">":
                node_seq = nodes[node].sequence
            else:
                node_seq = reverse_complement_seq(nodes[node].sequence)
        except:
            node_seq = ""

        path_seq += node_seq

    return path_seq


def swave_convert(options):

    if options.output_path is None:
        options.output_path = os.path.dirname(options.vcf_path)

    gfa = GFA()
    gfa.parse_gfa_file(options.gfa_path)

    with pysam.VariantFile(options.vcf_path) as fin, open(os.path.join(options.output_path, os.path.basename(options.vcf_path).replace(".vcf", ".converted.vcf")), "w") as fout, pysam.FastaFile(options.ref_path) as ref_file:

        fout.write(str(fin.header))

        previous_record_chr = None
        previous_record_pos = None

        for record in fin:

            record_split = str(record).strip().split("\t")

            # # convert the path to seq
            ref_path = record_split[3]
            alt_paths = record_split[4].split(",")

            ref_path_seq = convert_path_to_seq(ref_path, gfa.nodes)
            alt_path_seqs = [convert_path_to_seq(alt_path, gfa.nodes) for alt_path in alt_paths]

            snarl_start_node = gfa.nodes[re.findall(r'[><]([a-zA-Z0-9]+)', record.id.split("_")[0])[0]]
            snarl_start_pos = snarl_start_node.source[1] + snarl_start_node.length

            snarl_start_base = ref_file.fetch(record.contig, snarl_start_pos - 1, snarl_start_pos).upper()

            # # deal with empty sequence
            if ref_path_seq == "":
                ref_path_seq = snarl_start_base
            else:
                ref_path_seq = snarl_start_base + ref_path_seq

            for i in range(len(alt_path_seqs)):
                if alt_path_seqs[i] == "":
                    alt_path_seqs[i] = snarl_start_base
                else:
                    alt_path_seqs[i] = snarl_start_base + alt_path_seqs[i]

            # # this fix the error for pangenie-index 'GraphBuilder::GraphBuilder: reference allele given in VCF does not match allele in reference fasta file at that position'
            # # current we don't know why the ref seq provided by graph are not same as by reference file, but we can force it to be same
            if options.force_pangenie:
                ref_allele_seq = ref_file.fetch(record.contig, snarl_start_pos - 1, snarl_start_pos + len(ref_path_seq) - 1).upper()
                if not (ref_path_seq == ref_allele_seq):
                    ref_path_seq = ref_allele_seq

                # # remove overlapped records
                if record.contig == previous_record_chr and snarl_start_pos <= previous_record_pos:
                    continue

                previous_record_chr = record.contig
                previous_record_pos = snarl_start_pos + len(ref_path_seq)

                # # write to file
                record_split[1] = str(snarl_start_pos)
                record_split[3] = ref_path_seq
                record_split[4] = ",".join(alt_path_seqs)

                for i in range(8, len(record_split)):
                    record_split[i] = record_split[i].split(":")[0]

                fout.write("\t".join(record_split) + "\n")
            else:
                record_split[3] = ref_path_seq
                record_split[4] = ",".join(alt_path_seqs)

                fout.write("\t".join(record_split) + "\n")

            # # # fix length by the snarl start pos and record start pos
            # snarl_start_node = gfa.nodes[re.findall(r'[><]([a-zA-Z0-9]+)', record.id.split("_")[0])[0]]
            # snarl_end_node = gfa.nodes[re.findall(r'[><]([a-zA-Z0-9]+)', record.id.split("_")[0])[1]]
            #
            # snarl_start_pos = snarl_start_node.source[1] + snarl_start_node.length
            #
            # record_start_pos = record.start
            # record_start_base = ref_file.fetch(record.contig, record_start_pos, record_start_pos + 1).upper()
            #
            # print(len(ref_path_seq), [len(alt_seq) for alt_seq in alt_path_seqs])
            # print(ref_path_seq)
            # # # trim the seq length
            # if record_start_pos > snarl_start_pos:
            #     diff_length = record_start_pos - snarl_start_pos
            #
            #     print("record_start_pos > snarl_start_pos",  record_start_pos, snarl_start_pos, diff_length)
            #
            #     ref_path_seq = ref_path_seq[diff_length: ]
            #     alt_path_seqs = [alt_seq[diff_length: ] for alt_seq in alt_path_seqs]
            #
            #     # # deal with empty sequence
            #     if ref_path_seq == "":
            #         ref_path_seq = record_start_base
            #
            #     for i in range(len(alt_path_seqs)):
            #         if alt_path_seqs[i] == "":
            #             alt_path_seqs[i] = record_start_base
            #
            # # # append the seq length
            # elif record_start_pos < snarl_start_pos:
            #     diff_length = snarl_start_pos - record_start_pos
            #
            #     append_seq = ref_file.fetch(record.contig, record_start_pos, record_start_pos + diff_length).upper()
            #
            #     ref_path_seq = append_seq + ref_path_seq
            #
            #     for i in range(len(alt_path_seqs)):
            #         if alt_path_seqs[i] == "":
            #             alt_path_seqs[i] = record_start_base + alt_path_seqs[i]
            #
            #     # if diff_length > 50:
            #     print("record_start_pos < snarl_start_pos", record.contig, record_start_pos, snarl_start_pos)
            #
            # else:
            #     pass
            # # # write to file
            # record_split[3] = ref_path_seq
            # record_split[4] = ",".join(alt_path_seqs)
            #
            # fout.write("\t".join(record_split))


def parse_line_nodes(vcf_path_nodes, vcf_path_orients, line):

    line_split = line.strip().split("\t")

    alt_path = line_split[5].split(":")[0]

    if alt_path == ".":  # # why alt_path is '.': because the asm has no contig covering this snarl
        pass
    else:
        # # for start node
        snarl_start_node_id, snarl_start_node_orient = line_split[3][1:], line_split[3][0]

        if not (len(vcf_path_nodes) > 1 and snarl_start_node_id == vcf_path_nodes[-1]):
            vcf_path_nodes.append(snarl_start_node_id)
            vcf_path_orients.append("+" if snarl_start_node_orient == ">" else "-")

        # # for alt path nodes
        alt_path_include_nodes = re.findall(r's\d+', alt_path)
        alt_path_include_nodes_orients = re.findall(r'([^a-zA-Z0-9]+)', alt_path)

        for i in range(len(alt_path_include_nodes)):
            vcf_path_nodes.append(alt_path_include_nodes[i])
            vcf_path_orients.append("+" if alt_path_include_nodes_orients[i] == ">" else "-")

        # # for end node
        snarl_end_node_id, snarl_end_node_orient = line_split[4][1:], line_split[4][0]

        vcf_path_nodes.append(snarl_end_node_id)
        vcf_path_orients.append("+" if snarl_end_node_orient == ">" else "-")

def swave_convert_Plines(options):

    if options.output_path is None:
        options.output_path = os.path.dirname(options.gfa_path)


    ref_vcf_name = os.path.basename(options.ref_vcf_path)

    output_gfa_path = os.path.join(options.output_path, os.path.basename(options.gfa_path).replace(".gfa", ".convert_Plines.gfa"))

    with open(options.gfa_path) as fin, open(output_gfa_path, "w") as fout:
        if options.force_vg:
            for line in fin:
                if line[0] == "S":
                    fout.write("\t".join(line.strip().split("\t")[0: 3]) + "\n")
                else:
                    fout.write(line)
        else:
            fout.writelines(fin.readlines())

        for vcf_file in os.listdir(options.vcf_path):

            if ".vcf" not in vcf_file:
                continue

            vcf_name = os.path.basename(vcf_file)

            print(vcf_name)

            vcf_path_nodes = []
            vcf_path_orients = []

            if vcf_name != ref_vcf_name:
                with open(os.path.join(options.vcf_path, vcf_file)) as fin:
                    cnt = 0
                    for line in fin:
                        cnt += 1

                        parse_line_nodes(vcf_path_nodes, vcf_path_orients, line)

                vcf_path_str = ",".join("{}{}".format(vcf_path_nodes[i], vcf_path_orients[i]) for i in range(len(vcf_path_nodes)))

                vcf_Pline = "P\t{}\t{}\t*\n".format(vcf_name, vcf_path_str)

                fout.write(vcf_Pline)
            else:
                with open(os.path.join(options.vcf_path, vcf_file)) as fin:
                    previous_chrom = None

                    for line in fin:

                        line_chrom = line.strip().split("\t")[0]

                        if line_chrom != previous_chrom:
                            if previous_chrom != None:
                                vcf_path_str = ",".join("{}{}".format(vcf_path_nodes[i], vcf_path_orients[i]) for i in range(len(vcf_path_nodes)))
                                vcf_Pline = "P\t{}\t{}\t*\n".format(previous_chrom, vcf_path_str)

                                fout.write(vcf_Pline)

                            # # reset
                            previous_chrom = line_chrom
                            vcf_path_nodes = []
                            vcf_path_orients = []

                        parse_line_nodes(vcf_path_nodes, vcf_path_orients, line)

                    # # append the last chrom
                    vcf_path_str = ",".join("{}{}".format(vcf_path_nodes[i], vcf_path_orients[i]) for i in range(len(vcf_path_nodes)))
                    # if "chr" not in previous_chrom:
                    #     vcf_Pline = "P\tchr{}\t{}\t*\n".format(previous_chrom, vcf_path_str)
                    # else:
                    vcf_Pline = "P\t{}\t{}\t*\n".format(previous_chrom, vcf_path_str)

                    fout.write(vcf_Pline)