## Swave
Swave detects structural variants (SVs) and complex SVs from a pangenome graph constructed on one or multiple genome assemblies. Swave also supports SNP/Indel calling from 'minigraph' pangenome.  

<div align=left><img width=80% height=80% src="https://github.com/songbowang125/Swave/blob/main/workflow.svg"/></div> 


## License

Swave is free for non-commercial use by academic, government, and non-profit/not-for-profit institutions. A commercial version of the software is available and licensed through Xi’an Jiaotong University.
For more information, please contact with Songbo Wang (songbowang125@163.com) or Kai Ye (kaiye@xjtu.edu.cn).

## Contents

- [Installation](#installation)
- [Usage_1 (Demo data)](#usage_1-demo-data)
- [Usage_2 (User data for minigraph)](#usage_2-user-data-for-minigraph)
- [Usage_3 (User data for minigraph-cactus and pggb)](#usage_3-user-data-for-minigraph-cactus-and-pggb)
- [Usage_4 (SNP and indel calling)](#usage_4-snp-and-indel-calling)
- [Usage_5 (Genotyping short-reads)](#usage_5-genotyping-short-reads)
- [Output of Swave](#output-of-swave)
- [Recommended Settings](#Recommended-Settings)

## Installation

### Operation systems
Linux-based systems required, including but are not limited to: 
* MacOS, Big Sur
* Ubuntu (including Windows Subsystem)
* CentOS

### Install from source

```
## Get the source code
git clone https://github.com/songbowang125/Swave.git
cd Swave

## Create a conda environment for Swave
conda env create -f ./environment.yml 

## Activate the environment
conda activate Swave-env
```
### Additional requirements
For minigraph pangenome graph (**preferred**), ["minigraph"](https://github.com/lh3/minigraph) and ["gfatools"](https://github.com/lh3/gfatools) is requred. 


For MC/PGGB pangenome graph, ["MC"](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md) or ["pggb"](https://github.com/pangenome/pggb), and ["gfatools"](https://github.com/lh3/gfatools) is requred.   


### Run test


You need to enter the source code folder and type 'python Swave.py':
```commandline
cd Swave
python Swave.py -h
```



## Usage_1 (Demo data)

1. Download the all the files from ["Demo link"](https://drive.google.com/drive/folders/1t0L0pBVYk4GhKXg1IcSsVkCw7PtdbaI5?usp=drive_link), and store them in the same folder.
```commandline
The demo data contains a pangenome graph built on T2T reference and three samples (six haplotypes).
```

2. Run the following command:
```commandline
python Swave.py call --input_path assemblies-trio-hg00514.tsv  --ref_path t2t.chr20.fa --output_path ./ --gfa_path pg.minigraph.gfa --gfa_source minigraph
```

3. Check output file:
```commandline
One log file and two result files are generated, including:
    swave.log
    swave.sample_level.vcf (Multi-allelic outputs)
    swave.sample_level.split.vcf (Bi-allelic outputs by splitting the Multi-allelic outputs)
```

## Usage_2 (User data for minigraph)

### 1. Prepare your data:

1.1 Basiclly, you need to prepare the pangenome gfa file and all the used sequence fasta (or fa) files:
  ```commandline
    pangenome.gfa

    referece.fasta
    sample1_hap1.fasta
    sample1_hap2.fasta
    ...
    sampleN_hapN.fasta
  ```

1.2 If you don't have the pangenome gfa, you can construct one with ["Minigraph"](https://github.com/lh3/minigraph):
  ```commandline
    # # construct
    minigraph -t thread_num -cxggs referece.fasta sample1_hap1.fasta ... sampleN_hapN.fasta > pangenome.gfa
    
    # # convert to fasta format (Output format for gfatools: add a suffix .gfa2fa.fa to .gfa. Format cannot be changed)
    gfatools gfa2fa pangenome.gfa > pangenome.gfa.gfa2fa.fa 
  ```

1.3 And then, using a tab-separated assemblies.tsv file to store all the sample files. 

  You can find a demo tsv file in ["Demo tsv link"](https://drive.google.com/file/d/1QgRsWIsFQvLNIgaqXst9LK-cQjROjIgP/view?usp=drive_link).

  Each row in the TSV file represents a sample, with the first column indicating the sample name and the remaining columns listing paths to the haplotype FASTA files.

  **!!!! Use absolute path for each fasta file !!!!**
  ```commandline
    assemblies.tsv: 
        sample1    /path/to/sample1_hap1.fasta    /path/to/sample1_hap2.fasta    /path/to/sample1_hap3.fasta
        ...
        sampleN    /path/to/sampleN_hapN.fasta    ...
  ```

1.4 Call snarls (bubbles) with ["Minigraph --call"](https://github.com/lh3/minigraph) (Optional, but recommended):

For each fasta file (including the reference file), run:
  ```commandline
    minigraph -t thread_num -xasm --call pangenome.gfa referece.fasta > referece.fasta.bed
    minigraph -t thread_num -xasm --call pangenome.gfa sample1_hap1.fasta > sample1_hap1.fasta.bed
    ...
  ```

**!!!! The output BED files must be named by adding a .bed suffix to the Fasta name. !!!!**

**!!!! For example: referece.fasta -> referece.fasta.bed;    sample1_hap1.fasta -> sample1_hap1.fasta.bed !!!!**

**!!!! To save time, it is recommended to run minigraph --call in parallel. !!!!**

### 2. Run Swave
2.1 Check two things:
  ```commandline
    The bed outputs must be named by appending a .bed suffix to the corresponding FASTA filenames. 
  
    All BED files should be placed in the same directory as pangenome.gfa or in the output folder.
  ```
2.2 Run:
  ```commandline
    python Swave.py call --input_path assemblies.tsv  --ref_path referece.fasta --gfa_source minigraph --gfa_path pangenome.gfa --output_path ./
  ```
2.3 Check outputs:
```commandline
  One log file and two result files are generated, including:
    swave.log
    swave.sample_level.vcf (Multi-allelic outputs)
    swave.sample_level.split.vcf (Bi-allelic outputs by splitting the Multi-allelic outputs)
```
## Usage_3 (User data for minigraph-cactus and pggb)

### 1. Prepare your data:
1.1 Basiclly, you need to prepare the pangenome gfa file and the vcf file from vg deconstruct:
  ```commandline
    pangenome.gfa.gz
    
    pangenome.raw.vcf.gz
  ```

And then generate the fasta file for p
```commandline
  # # convert to fasta format (Output format for gfatools: add a suffix .gfa2fa.fa to gfa file. Format cannot be changed)
  gfatools gfa2fa pangenome.gfa.gz > pangenome.gfa.gz.gfa2fa.fa 
```
Swave uses the pangenome.raw.vcf.gz to extract allele information, and doesn't need to run the 'minigraph --call' commands. 

The pangenome.raw.vcf.gz file is an automated output part of the minigraph-cactus and pggb pipeline.

If you need to generate it by your self, run:
  ```commandline
  Single chromosome: 
    vg deconstruct -e -a -p chr1 pangenome.gfa.gz | bgzip > pangenome.raw.vcf.gz && tabix pangenome.raw.vcf.gz
  Multiple chromosome:
    vg deconstruct -e -a -P chr pangenome.gfa.gz | bgzip > pangenome.raw.vcf.gz && tabix pangenome.raw.vcf.gz
  ```
### 2. Run Swave

2.1 Run:
  ```commandline
    python Swave.py call --input_path assemblies.tsv  --ref_path referece.fasta --gfa_source cactus (or pggb) --gfa_path pangenome.gfa.gz --decomposed_vcf pangenome.raw.vcf.gz --output_path ./
  ```

If you meet the error bellow, please add the /path/to/ref.fa to your assemblies.tsv:
  ```
  [ERROR] Fail interpreting for bin 0. Error log: <class 'ValueError'>: 'ref.fa' is not in list. Locate At: [<FrameSummary file Swave/src/pack_sv/op_sv.py, line 902 in interpret_one_bin>, <FrameSummary file Swave/src/pack_sv/op_sv.py, line 1130 in generate_vcf_records>, <FrameSummary file Swave/src/pack_sv/op_sv.py, line 968 in output_at_top_and_sub_snarl_level>]
  ```

2.2 Check outputs:
```commandline
  One log file and two result files are generated, including:
    swave.log
    swave.sample_level.vcf (Multi-allelic outputs)
    swave.sample_level.split.vcf (Bi-allelic outputs by splitting the Multi-allelic outputs)
```
## Usage_4 (SNP and indel calling)
Swave is designed for SV calling from pangenome. However, we have implemented a script (main_small_call.py) for SNP/Indel calling from minigraph pangenome.

The script main_small_call.py requires ["minimap2"](https://github.com/lh3/minimap2) in your env, then run:

```commandline
Note: this SNP/Indel calling script is for minigraph pangenome only!

python main_small_call.py --input_path assemblies.tsv  --ref_path referece.fasta --output_path ./ --minigraph_callbed_folder {folder to minigraph --call bed files} --minimap2_path /path/to/minimap2 --thread_num {}
```

We tested the results using ["HG002 benchmark set"](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.020-20250117/) (Variants: CHM13v2.0_HG2-T2TQ100-V1.1_smvar.vcf.gz; High-confident-regions: CHM13v2.0_HG2-T2TQ100-V1.1_smvar.benchmark.bed).
 And this script achieved F1-score of 0.977, comparable to PAV using assemblies (0.974).


The MC/Pggb pangenome has already included base-level bubbles for SNP/Indels, so you can extract them as SNP/Indel.


## Usage_5 (Genotyping short-reads)

Using the pangenome-level output of Swave and Pangenie, you can genotyping other short-reads samples.

```commandline
 1. Convert the allele paths to allele sequences in the VCF. Please use the *_level.vcf, not the *_level.split.vcf.
    python Swave.py convert_seq --vcf_path swave.sample_level.vcf --gfa_path {} --ref_path referece.fasta --force_pangenie
 
 2. Pangenie index:
    PanGenie-index -v swave.sample_level.converted.vcf -r referece.fasta -t threads -o index-prefix
 
 3. Pangenie genotyping. As Pangenie required, the read1 and read2 should be stored in one fastq file
    PanGenie -f index-prefix -i short-reads.fastq -s sample_name -j threads -t threads -o output-prefix
```

## Output of Swave
### 1.Output files:

One log file and two result files are generated, including:

    swave.log
    swave.sample_level.vcf (Multi-allelic outputs)
    swave.sample_level.split.vcf (Bi-allelic outputs by splitting the Multi-allelic outputs)

The swave.sample_level.split.vcf follows standard VCF format and can be applied for downstream analysis.

### 2. Output filter:
 
Swave offers a three-level confidence tags in the outputted VCF file (column QUAL), including HighQual (renamed to PASS to fulfill the VCF format), MediumQual and LowQual.

The quality reflects the ratio between two lengths:

(1) length difference between reference and alternative sequence, while

(2) length difference between reference and predicted-SV-caused sequence.

If the predicted SV is true-positive, this ratio should be 0. In the contrast, if the SV is less-accurate or false-positive, this ratio would increase and be larger than 0.

SVs with 'LowQual' are recommended to be filtered out for downstream analysis.   

### 3. Get REF and ALT sequence
Both the vcf files contain graph paths in the REF and ALT column, you can convert them to sequence with:

    python Swave.py convert_seq --vcf_path {} --gfa_path {} --ref_path {}

## Recommended Settings
(1) It is best not to use ':' in the  chromosome names. Swave could pass this, but this would cause problems with the VCF or 'minigraph call' BED formats.

(2) If you are benchmarking with HG002 Tier1 callset (HG002_SVs_Tier1_v0.6.vcf), please add parameter '--dup_to_ins' 

(3) For minigraph pangenome, you can move all the 'minigraph --call' BED files to the folder where the .gfa locates. This will save much running time when you re-run Swave.

(4) If you re-run Swave in a previous output folder, please remove all the 'tmp*' files and 'tmp*' folders.

## Contact
If you have any questions, please feel free to contact: songbowang125@163.com or songbowang125@xjtu.edu.cn
