# PepSeq pipeline

The pipeline takes a list of snpeff annoatated vcf files to create a PepSeq Libray that is ready for ordering oligos.

## Software Installation

The PepSeq pipeline is a python package that uses python 3.7 and other packages.

The package provides a yaml file to create the conda/python3.7 environment for the pipeline.

#### Step 1: install Miniconda3
 
```
wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.11.0-Linux-x86_64.sh

bash ~/Miniconda3-py38_4.11.0-Linux-x86_64.sh

Follow the prompts on the installer screens. 
The installer prompts “Do you wish the installer to initialize Anaconda3 by running conda init?” We recommend “yes”.
```


#### Step 2: create a conda environment

```
conda env create -f conda_pepseq_design.yml
```

## Activate the python environment before start running the PepSeq pipeline

```
conda activate pepseq_design
```

make sure you have software tool 'bcftools' install in your shell environment, otherwise load bcftools from TGen HPC node:

```
module load bcftools
```

## Run pepseq_pipeline.py

Program options:

```
    Program arguments:
    
        --mutation_files: a csv file contains a list of vcf files in the following format:
                            sample_id1,ashion_vcf1,phoenix_vcf1
                            sample_id2,ashion_vcf2,phoenix_vcf2
        --library_id:     pepSeq library ID, e.g. TM1
        --caller_count:   minimum caller count for variant consensus calls by Phoenix pipeline
        --peptide_length: peptide length [15]
        --encode_dir:     directory that contains executable and accessory files for oligo encoding
```

## input file format

```
sampleID_1,Ashion_vcf_1,Phoenix_vcf_1
sampleID_2,Ashion_vcf_2,Phoenix_vcf_2
 ... ...
sampleID_n,Ashion_vcf_n,Phoenix_vcf_n
```

## output results

```
Lirary ordering file:      <LibraryID>.best_encodings.orderfile.csv

Lirary oligo encodng file: <LibraryID>.best_encodings.csv

intermediate files:        <scratch> directory
```
