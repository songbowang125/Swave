## Swave
Swave detects structural variants (SVs) and complex SVs from a pangenome graph constructed on one or multiple genome assemblies.  

<div align=left><img width=80% height=80% src="https://github.com/songbowang125/Swave/blob/main/workflow.svg"/></div> 


## License

Swave is free for non-commercial use by academic, government, and non-profit/not-for-profit institutions. A commercial version of the software is available and licensed through Xi’an Jiaotong University.
For more information, please contact with Songbo Wang (songbowang125@163.com) or Kai Ye (kaiye@xjtu.edu.cn).

## Contents

- [Installation](#installation)
- [Usage_1 (Demo data)](#usage_1-demo-data)
- [Usage_2 (User data for minigraph)](#usage_2-user-data-for-minigraph)
- [Usage_3 (User data for minigraph-cactus and pggb)](#usage_3-user-data-for-minigraph-cactus-and-pggb)


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
    sample1_hap3.fasta
    ...
    sampleN_hapN.fasta
  ```

1.2 If you don't have the pangenome gfa, you can construct one with ["Minigraph"](https://github.com/lh3/minigraph):
  ```commandline
    minigraph -t thread_num -cxggs referece.fasta sample1_hap1.fasta ... sampleN_hapN.fasta > pangenome.gfa
  ```

1.3 And then, using a tab-separated assemblies.tsv file to store all the sample files. 

  You can find a demo tsv file in ["Demo tsv link"](https://drive.google.com/file/d/1QgRsWIsFQvLNIgaqXst9LK-cQjROjIgP/view?usp=drive_link).

  Each row in the TSV file represents a sample, with the first column indicating the sample name and the remaining columns listing paths to the haplotype FASTA files.
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

Swave will automatically run the above commands if the bed files are not found in the pangenome.gfa folder.

However, Swave executes them one-by-one. To save time, it is recommended to run them in parallel manually.

### 2. Run Swave
2.1 Check two things:
  ```commandline
    The bed outputs must be named by appending a .bed suffix to the corresponding FASTA filenames. 
    (If you followed the above steps, either manually or automatically, you can ignore this note.) 
  
    All BED files should be placed in the same directory as pangenome.gfa.
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
Swave uses the pangenome.raw.vcf.gz to extract allele information, and doesn't need to run the 'minigraph --call' commands. 

The pangenome.raw.vcf.gz file is an automated output part of the minigraph-cactus and pggb pipeline.

If you need to generate it by your self, run 'vg deconstruct -e -a pangenome.gfa.gz | bgzip > pangenome.raw.vcf.gz && tabix pangenome.raw.vcf.gz'


### 2. Run Swave

2.1 Run:
  ```commandline
    python Swave.py call --input_path assemblies.tsv  --ref_path referece.fasta --gfa_source cactus (or pggb) --gfa_path pangenome.gfa.gz --decomposed_vcf pangenome.raw.vcf.gz --output_path ./
  ```
2.2 Check outputs:
```commandline
  One log file and two result files are generated, including:
    swave.log
    swave.sample_level.vcf (Multi-allelic outputs)
    swave.sample_level.split.vcf (Bi-allelic outputs by splitting the Multi-allelic outputs)
```

## Contact
If you have any questions, please feel free to contact: songbowang125@163.com
