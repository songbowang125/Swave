from src.version import __version__
import pysam


def output_vcf_header(vcf_fout, options):
    """
    generate vcf header
    """

    print("##fileformat=VCFv4.3", file=vcf_fout)
    print("##source=xxx v{0}".format(__version__), file=vcf_fout)

    # # STEP: add chromosome info
    ref_file = pysam.FastaFile(options.ref_asm_path)
    chroms = ref_file.references
    for chr in chroms:
        chr_length = ref_file.get_reference_length(chr)
        print('##contig=<ID={0},length={1}>'.format(chr, chr_length), file=vcf_fout)
    ref_file.close()

    # # STEP: add info field
    # print("##QUAL=<QUAL=XXX,Description=\"The SV quality of the SV described in this region\">", file=vcf_fout)

    print("##ID=<ID=Path,Description=\"Snarl path\">", file=vcf_fout)
    print("##REF=<ID=Path,Description=\"Reference path\">", file=vcf_fout)
    print("##ALT=<ID=Path,Description=\"Alternative paths\">", file=vcf_fout)

    print("##FILTER=<ID=PASS,Description=\"Passed\">", file=vcf_fout)
    print("##FILTER=<ID=LowQual,Description=\"LowQual\">", file=vcf_fout)
    print("##FILTER=<ID=MediumQual,Description=\"MediumQual\">", file=vcf_fout)

    print("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the SV/CSV\">", file=vcf_fout)
    print("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV/CSV\">", file=vcf_fout)
    print("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the SV/CSV\">", file=vcf_fout)
    print("##INFO=<ID=BKPS,Number=.,Type=String,Description=\"Breakpoints of the SV/CSV (length-start-end-insert)\">", file=vcf_fout)

    # print("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Population allele frequency\">", file=vcf_fout)
    # print("##INFO=<ID=Samples,Number=1,Type=Integer,Description=\"Supporting Samples for this SV/CSV\">", file=vcf_fout)

    print("##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Total number of alternate alleles in called genotypes\">", file=vcf_fout)
    print("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Estimated allele frequency in the range (0,1]\">", file=vcf_fout)
    print("##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">", file=vcf_fout)
    print("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">", file=vcf_fout)
    print("##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Level in the snarl tree (0=top level)\">", file=vcf_fout)

    # # STEP: add gt info
    print("##FORMAT=<ID=GT,Number=.,Type=String,Description=\"SV genotype for each haplotype\">", file=vcf_fout)
    print("##FORMAT=<ID=TYPE,Number=.,Type=String,Description=\"SV type for each haplotype\">", file=vcf_fout)
    print("##FORMAT=<ID=LENGTH,Number=.,Type=String,Description=\"SV length for each haplotype\">", file=vcf_fout)
    print("##FORMAT=<ID=QUAL,Number=.,Type=String,Description=\"SV quality for each haplotype\">", file=vcf_fout)
    print("##FORMAT=<ID=BKPS,Number=.,Type=String,Description=\"SV breakpoint for each haplotype\">", file=vcf_fout)

    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}".format( "\t".join(options.alt_asm_names)), file=vcf_fout)


def output_vcf_records(vcf_fout, records):

    # # sort records by chrom and start positions
    records = sorted(records, key=lambda x: (x[0], int(x[1])))

    for record in records:
        print("\t".join(record), file=vcf_fout)

def hap_gt_to_sample_gt(sample_gt, mode="raw"):

    if mode == "split":
        gt_flags = []
        for hap_gt in sample_gt:
            gt_flags.append(hap_gt)

        return "{}".format("|".join(gt_flags))

    else:
        gt_flags = []
        sv_types = []
        sv_lengths = []
        # sv_quals = []
        sv_bkps = []

        for hap_gt in sample_gt:
            hap_gt_split = hap_gt.split(":")

            gt_flags.append(hap_gt_split[0])
            sv_types.append(hap_gt_split[1])
            sv_lengths.append(hap_gt_split[2])
            # sv_quals.append(hap_gt_split[3])
            sv_bkps.append(hap_gt_split[3])

        return "{}:{}:{}:{}".format("|".join(gt_flags), "|".join(sv_types), "|".join(sv_lengths), "|".join(sv_bkps))

        # return "{}:{}:{}:{}:{}".format("|".join(gt_flags), ",".join(sv_types), ",".join(sv_lengths), ",".join(sv_quals), ",".join(sv_bkps))


def generate_sample_level_vcf(hap_level_path, options, mode="raw"):

    sample_level_vcf_path = hap_level_path.replace("hap_level", "sample_level").replace("tmp.", "")

    sample_names = options.spec_samples

    sample_hap_files = []
    for sample in sample_names:
        sample_hap_files.append(options.input_dict[sample])

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

