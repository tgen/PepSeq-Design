import os
import sys
import shutil
import argparse
import pandas as pd
import subprocess as proc
import h2o
from utils import *

def main(args):
    # 1. initialize object: tumor mutation peptides
    tm = Neopeptide(args)

    # 2. convert variants to peptides
    tm.make_prot_seq()

    # 3. read peptide file(s)
    tm.read_peptide_file()

    # 4. convert peptides to oligo encodings
    tm.oligo_encoding()

    # 5. run machine learning model to rank best oligo encodings
    tm.oligo_scoring()

if __name__ == "__main__":
    '''Set up before running "pepseq_pipeline.py":
    
        Conda3 or Miniconda3: install anaconda3 or miniconda3
        python3 packages:     run "conda env create -f conda_pepseq_design.yml" to install python3 packages
        bcftools:             "module load bcftools" or install bcftools in your environment
        
    Program arguments:
    
        --mutation_files: a csv file contains a list of vcf files in the following format:
                            sample_id1,ashion_vcf1,phoenix_vcf1
                            sample_id2,ashion_vcf2,phoenix_vcf2
        --library_id:     pepSeq library ID, e.g. TM1
        --caller_count:   minimum caller count for variant consensus calls by Phoenix pipeline
        --peptide_length: peptide length [15]
        --encode_dir:     directory that contains executable and accessory files for oligo encoding
    '''
    parser = argparse.ArgumentParser(prog="pepseq_pipeline.py", description="Run PepSeq pipeline")
    parser.add_argument('-i', '--mutation_files', required=True, help="mutation files in csv format")
    parser.add_argument('-l', '--library_id', help="PepSeq library ID")
    parser.add_argument('-c', '--caller_count', default=4, help="mininum variant caller count for Phoenix pipeline")
    parser.add_argument('-p', '--peptide_length', default=15, help="neopeptide length")
    parser.add_argument('-w', '--working_dir', default="scratch", help="working directory")
    parser.add_argument('-e', '--encode_dir', default="oligo_encoding",
                        help="directory containing executable and accessory files for nucleotide encoding")

    args, remaining_argv = parser.parse_known_args()

    # create scratch dir and "scratch/varcode" dir
    args.varcode_dir = os.path.join(args.working_dir, 'varcode')
    os.makedirs(args.varcode_dir, exist_ok=True)

    main(args)
