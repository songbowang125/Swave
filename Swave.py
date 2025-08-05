import os.path
import datetime

import logging
import pysam
import sys
import argparse
from src.version import __version__
import pkg_resources


def parse_arguments(arguments=sys.argv[1:]):

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="Swave")

    subparsers = parser.add_subparsers(title='modules', dest='command', required=True)

    # # parameters for call
    parser_call = subparsers.add_parser('call', help='Call structural variants from inputs')

    parser_call.add_argument('--output_path', dest="output_path", type=os.path.abspath, required=True, help='Absolute path to output')
    parser_call.add_argument('--input_path', dest="input_path", type=os.path.abspath, required=True, help='Absolute path to input')
    parser_call.add_argument('--ref_path', dest="ref_path", type=os.path.abspath, required=True, help='Absolute path to reference')
    parser_call.add_argument('--gfa_path', dest="gfa_path", type=os.path.abspath, required=True, help='Absolute path to gfa')

    parser_call.add_argument('--gfa_source', dest="gfa_source", choices=['minigraph', 'cactus', 'pggb'], default='minigraph', help='The source of GFA. Choose from [minigraph, cactus, pggb].')
    parser_call.add_argument('--decomposed_vcf', dest="decomposed_vcf", type=os.path.abspath, required=False, default=None, help='The decomposed VCF from vg (could be processed by vcfbub and vcfwave), required for cactus and pggb graphs')

    parser_call.add_argument('--output_mode', dest="output_mode", choices=['auto', 'population', 'single'], default='auto', help='Mode of output. Choose from [auto, population, single]')

    parser_call.add_argument('--spec_samples', dest='spec_samples', metavar="IN", nargs="+", type=str, default=["all"], help='Candidate samples')
    parser_call.add_argument('--spec_snarl', dest='spec_snarl', type=str, default=None, help="Only call the specific snarl (e.g. '>s1>s2')")
    parser_call.add_argument('--spec_path', dest='spec_path', type=str, default=None, help="Only call the specific path (e.g. '>s3>s4')")

    parser_call.add_argument("--min_sv_size", dest="min_sv_size", type=int, default=50, help='Minimum SV size to detect (default: %(default)s)')
    parser_call.add_argument("--max_sv_size", dest="max_sv_size", type=int, default=1000000, help='Maximum SV size to detect (default: %(default)s)')

    parser_call.add_argument("--max_sv_comps", dest="max_sv_comps", type=int, default=5, help='Maximum number of SV components for detailed output, otherwise, hyperCPX')

    parser_call.add_argument('--dup_to_ins', action="store_true", default=False, help='Report duplications as insertions')
    parser_call.add_argument('--remove_small', action="store_true", default=False, help='Remove nodes smaller than min_sv_size')
    parser_call.add_argument('--force_reverse', action="store_true", default=False, help='Force call reversed mapping snarls')

    parser_call.add_argument('--thread_num', dest="thread_num", type=int, default=1, help='Thread numbers (default: %(default)s)')

    parser_call.add_argument("--minigraph", dest="minigraph", type=str, default="minigraph", help='Path to executable minigraph path')
    parser_call.add_argument("--gfatools", dest="gfatools", type=str, default="gfatools", help='Path to executable gfatools path')

    parser_convert = subparsers.add_parser('convert_seq', help='Convert the graph paths to sequence in the REF and ALT columns of VCF')
    parser_convert.add_argument('--vcf_path', dest="vcf_path", type=os.path.abspath, required=True, help='Absolute path to vcf')
    parser_convert.add_argument('--gfa_path', dest="gfa_path", type=os.path.abspath, required=True, help='Absolute path to gfa')
    parser_convert.add_argument('--ref_path', dest="ref_path", type=os.path.abspath, required=True, help='Absolute path to reference')
    parser_convert.add_argument('--output_path', dest="output_path", type=os.path.abspath, default=None, help='Absolute path to output')
    parser_convert.add_argument('--force_pangenie', action="store_true", default=False, help='Force output sequences to meet the requirement of pangenie index')

    parser_convert_Plines = subparsers.add_parser('convert_Plines', help='Convert the minigraph call VCF into GFA P lines')
    parser_convert_Plines.add_argument('--gfa_path', dest="gfa_path", type=os.path.abspath, required=True, help='Absolute path to gfa')
    parser_convert_Plines.add_argument('--vcf_path', dest="vcf_path", type=os.path.abspath, required=True, help='Absolute path to vcf')
    parser_convert_Plines.add_argument('--ref_vcf_path', dest="ref_vcf_path", type=os.path.abspath, default=None, help='Absolute path to vcf')
    parser_convert_Plines.add_argument('--output_path', dest="output_path", type=os.path.abspath, default=None, help='Absolute path to output')
    parser_convert_Plines.add_argument('--force_vg', action="store_true", default=False, help='Force output sequences to meet the requirement of pangenie index')

    parser_extract_csv = subparsers.add_parser('extract_csv', help='Extract CSV from VCF')
    parser_extract_csv.add_argument('--vcf_path', dest="vcf_path", type=os.path.abspath, required=True, help='Absolute path to vcf')
    parser_extract_csv.add_argument('--spec_csv', dest="spec_csv", choices=['INV', 'DUP', "All"], default=None, help='Specific a relevant CSV type')
    parser_extract_csv.add_argument('--output_path', dest="output_path", type=os.path.abspath, default=None, help='Absolute path to output')

    parser_extract_sample = subparsers.add_parser('extract_sample', help='Extract SVs for specific samples from VCF')
    parser_extract_sample.add_argument('--vcf_path', dest="vcf_path", type=os.path.abspath, required=True, help='Absolute path to vcf')
    parser_extract_sample.add_argument('--spec_samples', dest='spec_samples', metavar="IN", nargs="+", type=str, required=True, help='Candidate samples')
    parser_extract_sample.add_argument('--output_path', dest="output_path", type=os.path.abspath, default=None, help='Absolute path to output')

    return parser.parse_args(arguments)


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

            input_dict[sample_name] = hap_paths
            input_list.extend(hap_paths)

            input_hap_numbers.append(len(hap_paths))

    return input_dict, input_list


def check_arguments(options):

    """
    check auguments and assign defaults
    """

    if options.command in ["add"]:
        options.gfa_source = "minigraph"
        options.spec_snarl = None
        options.spec_path = None
        options.spec_samples = ["all"]
        options.output_mode = "population"

    else:
        if options.gfa_source in ["cactus", "pggb"]:
            if options.decomposed_vcf is None:
                logging.error("For cactus and pggb graphs, the vg decomposed VCF is required")
                exit(-1)

    options.img_output_path = os.path.join(options.output_path, "tmp_imgs")
    options.npz_output_path = os.path.join(options.output_path, "tmp_npzs")

    options.kmer_size = 30

    options.model_path = pkg_resources.resource_filename('src', 'pack_model/LSTM-l1-fc64-bi.pth')

    if not os.path.exists(options.model_path):
        options.model_path = "src/pack_model/LSTM-l1-fc64-bi.pth"

    options.model_device = "cpu"
    options.model_cpu = 8

    # #
    ref_file = pysam.FastaFile(options.ref_path)
    options.available_chroms = ref_file.references
    ref_file.close()

    options.input_dict, options.input_list = load_from_sample_tsv(options.input_path)

    # # STEP: determine sample numbers and output mode
    if options.spec_samples == ["all"]:
        options.spec_samples = list(options.input_dict.keys())

    if options.output_mode == "auto":
        if len(options.spec_samples) > 1:
            options.output_mode = "population"
        else:
            options.output_mode = "single"
    # #
    options.asm_paths = [options.ref_path]
    for sample in options.spec_samples:
        if sample not in options.input_dict:
            logging.error("Input sample name '{}' not shown in the input tsv file".format(sample))
            exit(-1)

        options.asm_paths.extend(options.input_dict[sample])

    # #
    options.ref_asm_path = options.asm_paths[0]
    options.alt_asm_paths = options.asm_paths[1:]

    options.asm_names = [os.path.basename(asm_path) for asm_path in options.asm_paths]
    options.ref_asm_name = options.asm_names[0]
    options.alt_asm_names = options.asm_names[1:]

    # # # STEP: check if all input fasta file exit
    # for asm_path in options.asm_paths:
    #     if not os.path.exists(asm_path):
    #         logging.error("Input file not found {}".format(asm_path))
    #         exit(-1)

    if not os.path.exists(options.output_path):
        os.mkdir(options.output_path)

    if not os.path.exists(options.img_output_path):
        os.mkdir(options.img_output_path)

    if not os.path.exists(options.npz_output_path):
        os.mkdir(options.npz_output_path)


if __name__ == '__main__':
    start_time = datetime.datetime.now()

    # # STEP: load parameters
    options = parse_arguments()

    if options.command in ["call"]:

        check_arguments(options)

        # # STEP: create log file
        log_format = logging.Formatter("%(asctime)s [%(levelname)s]    %(message)s")

        root_logger = logging.getLogger()
        root_logger.setLevel(logging.INFO)

        fileHandler = logging.FileHandler(os.path.join(options.output_path, "swave.log"), mode="w")
        fileHandler.setFormatter(log_format)

        consoleHandler = logging.StreamHandler()
        consoleHandler.setFormatter(log_format)
        root_logger.addHandler(consoleHandler)
        root_logger.addHandler(fileHandler)

        logging.info("Command: {}".format(" ".join(sys.argv)))

        # # STEP: execute commands
        from main_call import swave_call
        swave_call(None, options)

    elif options.command in ["convert_seq", "convert_Plines", "extract_sample"]:

        if options.command == 'convert_seq':
            from main_convert import swave_convert
            swave_convert(options)

        elif options.command == "convert_Plines":
            from main_convert import swave_convert_Plines
            swave_convert_Plines(options)

        elif options.command == "extract_sample":
            from main_extract import swave_extract_sample
            swave_extract_sample(options)

    # # STEP: end
    end_time = datetime.datetime.now()
    cost_time = (end_time - start_time).seconds
    logging.info("Swave {} v{} successfully finished. Time Cost {}s".format(options.command, __version__, cost_time))
