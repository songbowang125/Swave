import logging
import os
import pysam



def swave_extract_sample(options):
    if options.output_path is None:
        options.output_path = os.path.dirname(options.vcf_path)

    output_vcf_path = os.path.join(options.output_path, os.path.basename(options.vcf_path).replace(".vcf", ".extract_{}.vcf".format("_".join(options.spec_samples))))

    with pysam.VariantFile(options.vcf_path) as fin, open(output_vcf_path, "w") as fout:

        fout.write(str(fin.header))

        sample_list = str(fin.header).strip().split("\n")[-1].split("\t")[9: ]

        # # check existence
        for spec_sample in options.spec_samples:
            if spec_sample not in sample_list:
                logging.error("{} is not included in the VCF".format(spec_sample))
                exit()

        for record in fin:

            sample_gts = str(record).strip().split("\t")[9: ]

            exist_sample_cnt = 0

            for spec_sample in options.spec_samples:
                spec_sample_gt = sample_gts[sample_list.index(spec_sample)].split(":")[0]

                if spec_sample_gt not in ["./.", "0/0", "./0", "0/.", ".|.", "0|0", ".|0", "0|.", ".", "0"]:
                    exist_sample_cnt += 1

            if exist_sample_cnt > 0:
                fout.write(str(record))
