import os
import sys
import re
import pandas as pd
import varcode
import subprocess as proc

class Neopeptide:
    def __init__(self, args):
        self.library_id = args.library_id
        self.caller_count = args.caller_count
        self.peptide_length = args.peptide_length
        self.samples = {}
        with open(args.mutation_files) as f:
            for line in f:
                sample_id, ashion_vcf, phoenix_vcf = line.strip().split(',')
                self.samples[sample_id] = [(ashion_vcf, 'Ashion'), (phoenix_vcf, 'Research')]
        self.encode_dir = args.encode_dir
        self.working_dir = args.working_dir
        self.varcode_dir = args.varcode_dir

    def oligo_encoding(self):
        '''This method is adapted from John Altin <jaltin@tgen.org>, Zane Fink <zanef2@illinois.edu>,
           and Erin Kelly <ekelley@tgen.org>

        sources:
        1. https://github.com/LadnerLab/Library-Design/
        2. https://github.com/TGenNorth/TM2_PepSeq/blob/master/analysis/bash/library_design/run_pepseq_design_step1.sh

        :return: None
        '''
        self.named_peptides = os.path.join(self.working_dir, '{}.named_peptides.csv'.format(self.library_id))
        self.encoded_seqs = os.path.join(self.working_dir, '{}.encoded_seqs.csv'.format(self.library_id))
        self.encoded_ratio = os.path.join(self.working_dir, '{}.encoded_ratio.csv'.format(self.library_id))

        cmd_encoding = os.path.join(self.encode_dir, 'oligo_encoding')
        cmd_encoding += ' -i {}'.format(self.named_peptides)
        cmd_encoding += ' -p {}'.format(os.path.join(self.encode_dir, 'codon_weights.csv'))
        cmd_encoding += ' -r {}'.format(self.encoded_ratio)
        cmd_encoding += ' -s {}'.format(self.encoded_seqs)
        cmd_encoding += ' -n 300 -c 2 -t 10000 -g 0.55'

        if not os.path.isfile(self.encoded_seqs):
            proc.call(cmd_encoding, shell=True)

    def oligo_scoring(self):
        '''This method is adapted from John Altin <jaltin@tgen.org>, Zane Fink <zanef2@illinois.edu>,
           and Erin Kelly <ekelley@tgen.org>

        sources:
        1. https://github.com/LadnerLab/Library-Design/
        2. https://github.com/TGenNorth/TM2_PepSeq/blob/master/analysis/bash/library_design/run_pepseq_design_step2.sh

        :return: None
        '''
        self.best_encodings = '{}.best_encodings.csv'.format(self.library_id)
        self.encoding_orderfile = '{}.best_encodings.orderfile.csv'.format(self.library_id)

        cmd_encoding = 'python ' + os.path.join(self.encode_dir, 'encoding_with_nn.py')
        cmd_encoding += ' -r {}'.format(self.encoded_ratio)
        cmd_encoding += ' -s {}'.format(self.encoded_seqs)
        cmd_encoding += ' -m {}'.format(os.path.join(self.encode_dir, 'DeepLearning_model_R_1539970074840_1_20181019'))
        cmd_encoding += ' -o {}'.format(self.best_encodings)
        cmd_encoding += ' --subsample 300 --read_per_loop 10 -n 3'

        if not os.path.isfile(self.best_encodings):
            proc.call(cmd_encoding, shell=True)

        fw = open(self.encoding_orderfile, 'w')
        fw.write('Seq_ID, Nucleotide_encoding_with_adapters\n')
        with open(self.best_encodings) as f:
            for line in f:
                a = line.strip().split(',')
                seq_id, nucleotide = a[0], a[-1]
                fw.write(','.join([seq_id, nucleotide]) + '\n')
        fw.close()

    def read_peptide_file(self):
        # self.samples[sample_id] = [(ashion_vcf, 'Ashion'), (phoenix_vcf, 'Research'), merged_peptide_file]
        w = {}
        for sample_id, sample_files in self.samples.items():
            file = sample_files[2]
            df = pd.read_csv(file, dtype=str, usecols=['variantId', 'effectId', 'gene_name', 'topEffect', 'sourceName',
                                                       'mutPept7', 'mutPept8', 'mutPept9', 'wtPept7', 'wtPept8',
                                                       'wtPept9']).fillna(value='NA')
            df = df.reset_index()  # make sure indexes pair with number of rows
            for index, row in df.iterrows():
                if bool(row['topEffect']) is not True:
                    continue
                if len(row['mutPept7']) == 15:
                    w[row['mutPept7']] = [sample_id, row['variantId'], row['effectId'], row['gene_name'],
                                          row['topEffect'], row['sourceName'], 'mutPept7', file]
                if len(row['mutPept8']) == 15:
                    w[row['mutPept8']] = [sample_id, row['variantId'], row['effectId'], row['gene_name'],
                                          row['topEffect'], row['sourceName'], 'mutPept8', file]
                if len(row['mutPept9']) == 15:
                    w[row['mutPept9']] = [sample_id, row['variantId'], row['effectId'], row['gene_name'],
                                          row['topEffect'], row['sourceName'], 'mutPept9', file]
                if len(row['wtPept7']) == 15:
                    w[row['wtPept7']] = [sample_id, row['variantId'], row['effectId'], row['gene_name'],
                                         row['topEffect'], row['sourceName'], 'wtPept7', file]
                if len(row['wtPept8']) == 15:
                    w[row['wtPept8']] = [sample_id, row['variantId'], row['effectId'], row['gene_name'],
                                         row['topEffect'], row['sourceName'], 'wtPept8', file]
                if len(row['wtPept9']) == 15:
                    w[row['wtPept9']] = [sample_id, row['variantId'], row['effectId'], row['gene_name'],
                                         row['topEffect'], row['sourceName'], 'wtPept9', file]
        # assign peptide_id and write out "named_peptides.csv"
        self.named_peptides = os.path.join(self.working_dir, '{}.named_peptides.csv'.format(self.library_id))
        fh1 = open(self.named_peptides, 'w')
        i = 0
        for k, v in sorted(w.items()):
            i += 1
            peptide_id = '_'.join([v[0], '{:05d}'.format(i)])
            fh1.write('{0},{1}\n'.format(peptide_id, k))
        fh1.close()

        # write out "annotated_peptides.csv"
        self.annotated_peptides = os.path.join(self.working_dir, '{}.annotated_peptides.csv'.format(self.library_id))
        fh2 = open(self.annotated_peptides, 'w')
        csv_header = "sampleId,peptide_source,peptide_seq,variantId,effectId,gene_name,topEffect,pipeline_source,source_file"
        fh2.write(csv_header + '\n')
        for k, v in sorted(w.items()):
            cols = [v[0], v[6], k, v[1], v[2], v[3], v[4], v[5], v[7]]
            fh2.write(','.join(cols) + '\n')
        fh2.close()

    def make_prot_seq(self):
        for sample_id, vcf_list in self.samples.items():
            for vcf_path, vcf_source in vcf_list:
                if vcf_path == '-':
                    continue
                if vcf_source == 'Ashion':
                    peptide_path_ashion = self._make_prot_seq(sample_id, vcf_path, vcf_source, 'GRCh37')
                    print('Ashion peptide file: {}'.format(peptide_path_ashion))
                if vcf_source == 'Research':
                    vcf_file = os.path.basename(vcf_path)
                    if vcf_file.endswith('.vcf.gz'):
                        vcf_root = '.'.join(vcf_file.split('.')[:-2])
                    elif vcf_file.endswith('.vcf'):
                        vcf_root = '.'.join(vcf_file.split('.')[:-1])
                    else:
                        vcf_root = vcf_file
                        print("Warning: {} does not end with .vcf or .vcf.gz".format(vcf_file), file=sys.stderr)
                    vcf_filt = os.path.join(self.varcode_dir, '{}.filt.vcf'.format(vcf_root))
                    cmd_filt = 'bcftools filter -i "INFO/CC>={0}" -Ov -o {1} {2}'.format(self.caller_count, vcf_filt, vcf_path)
                    proc.call(cmd_filt, shell=True)
                    peptide_path_research = self._make_prot_seq(sample_id, vcf_filt, vcf_source, 'hg38')
                    print('TGen peptide file: {}'.format(peptide_path_research))
            # merge peptide.csv files for Ashion and Research/TGen
            if os.path.isfile(peptide_path_ashion) and os.path.isfile(peptide_path_research):
                peptide_path_merged = self._merge_prot_seq(peptide_path_ashion, peptide_path_research)
                self.samples[sample_id].append(peptide_path_merged)
                print('both Ashion and Research mutation files exist, and merged: {}'.format(peptide_path_merged))
            elif os.path.isfile(peptide_path_ashion):
                self.samples[sample_id].append(peptide_path_ashion)
            elif os.path.isfile(peptide_path_research):
                self.samples[sample_id].append(peptide_path_research)
            else:
                os.exit('{} does not produce mutation peptide file'.format(sample_id))

    def _merge_prot_seq(self, peptide_path_ashion, peptide_path_research):
        peptide_research = {}
        with open(peptide_path_research) as f1:
            header = f1.readline()
            for line in f1:
                a = line.strip().split(',')
                if a[1] == "":  # mutPept1 == ""
                    continue
                # dict: key=(mutPept1, ..., mutPept15), value=[gene_name, source_name, line]
                mutPept_list = tuple(a[1:16])
                gene_name = a[17]
                source_name = a[-1]
                if mutPept_list in peptide_research:
                    print('Warning: duplicated mutPept list in Research/TGen peptide set: {}, {}'.format(gene_name, source_name))
                else:
                    peptide_research[mutPept_list] = [gene_name, source_name, line]
        # write out merged '.varCode.csv' file
        peptide_path_merged = peptide_path_ashion.replace('Ashion.varCode.csv', 'Ashion_Research.merged.varCode.csv')
        fw = open(peptide_path_merged, 'w')
        with open(peptide_path_ashion) as f2:
            header = f2.readline()
            fw.write(header)
            for line in f2:
                a = line.strip().split(',')
                if a[1] == "":  # mutPept1 == ""
                    continue
                mutPept_list = tuple(a[1:16])
                gene_name_ashion = a[17]
                source_name_ashion = a[-1]
                # write out record in both Ashion set and Research set
                if mutPept_list in peptide_research:
                    gene_name_research, source_name_research, line_research = peptide_research[mutPept_list]
                    new_line = ','.join(a[:-1]+[source_name_ashion+'|'+source_name_research]) + '\n'
                    fw.write(new_line)
                    # remove mutPept_list from Research/TGen set
                    peptide_research.pop(mutPept_list, None)
                    if gene_name_ashion == gene_name_research:
                        print('gene_name matches: {} and {}'.format(gene_name_ashion, gene_name_research))
                    else:
                        print('mutPept* match BUT gene_name DO NOT match: GRCh37={} and GRCh38={}'.format(gene_name_ashion, gene_name_research))
                else:
                    # write out record only in Research set
                    fw.write(line)
        # write out extra mutation/genes picked up in Research/TGen pipeline
        for item in peptide_research:
            fw.write(item[2])
        fw.close()
        return peptide_path_merged

    def _make_prot_seq(self, sample_id, vcf, vcf_source, genome=None):
        '''This method is adapted from Rebecca Halperin <rhalperin@tgen.org> and Kevin Drenner <kdrenner@tgen.org>

        sources:
        1. NeoNox2 (rhalperin): https://github.com/tgen/NeoNox2/
        2. ${HOME}/NeoNox2/makeProtSeq.py

        :param sample_id:   sample id
        :param vcf:         snpeff annotated vcf file path
        :param vcf_source:  pipeline name
        :param genome:      GRCh37 or hg38
        :return:            peptide file path
        '''
        source_name = '{}_{}'.format(sample_id, vcf_source)
        peptide_file = '{}_{}.varCode.csv'.format(self.library_id, source_name)
        peptide_path = os.path.join(self.varcode_dir, peptide_file)
        vcfVariants = varcode.load_vcf(vcf, genome=genome)
        vcfEffects = vcfVariants.effects()
        nonSilentMutations = vcfEffects.drop_silent_and_noncoding()
        topEffects = nonSilentMutations.top_priority_effect_per_variant()
        effectList = nonSilentMutations.effects
        topEffectId = dict()
        peptideLength = self.peptide_length
        mutPept = "mutPept"
        wtPept = "wtPept"

        for myVariant, myEffect in topEffects.items():
            topEffectId[myVariant.short_description] = '_'.join([myEffect.transcript_id, myEffect.short_description])

        effectDF = pd.DataFrame(index=range(0, len(effectList)),
                                columns=["variantId", "effectId", "gene_name", "effectClass", "topEffect", "mutPept1",
                                         "mutPept2", "mutPept3", "mutPept4", "mutPept5", "mutPept6", "mutPept7",
                                         "mutPept8", "mutPept9", "mutPept10", "mutPept11", "mutPept12", "mutPept13",
                                         "mutPept14", "mutPept15", "wtPept1", "wtPept2", "wtPept3", "wtPept4",
                                         "wtPept5", "wtPept6", "wtPept7", "wtPept8", "wtPept9", "wtPept10",
                                         "wtPept11", "wtPept12", "wtPept13", "wtPept14", "wtPept15"])
        for i in range(0, len(effectList)):
            myEffect = effectList[i]
            effectDF.loc[i, "variantId"] = myEffect.variant.short_description
            effectDF.loc[i, "effectId"] = '_'.join([myEffect.transcript_id, myEffect.short_description])
            effectDF.loc[i, "gene_name"] = myEffect.gene_name
            effectDF.loc[i, "effectClass"] = type(myEffect)
            if myEffect.variant.short_description in topEffectId:
                if effectDF.loc[i, "effectId"] == topEffectId[myEffect.variant.short_description]:
                    effectDF.loc[i, "topEffect"] = True
                else:
                    effectDF.loc[i, "topEffect"] = False
            else:
                effectDF.loc[i, "topEffect"] = False
            for j in range(1, peptideLength + 1):
                if isinstance(myEffect.mutant_protein_sequence, str) and isinstance(myEffect.aa_mutation_start_offset, int):
                    if re.search('del', myEffect.variant.short_description):
                        if re.search('\*', myEffect.short_description):
                            effectDF.loc[i, mutPept + str(j)] = ""
                            effectDF.loc[i, wtPept + str(j)] = ""
                        elif (myEffect.aa_mutation_start_offset + peptideLength) > len(myEffect.mutant_protein_sequence):
                            if j == 1:
                                endLength = len(myEffect.mutant_protein_sequence) - myEffect.aa_mutation_start_offset
                                difference = peptideLength - endLength
                                effectDF.loc[i, mutPept + str(j)] = myEffect.mutant_protein_sequence[
                                    myEffect.aa_mutation_start_offset - difference:myEffect.aa_mutation_start_offset + (endLength)]
                                effectDF.loc[i, wtPept + str(j)] = myEffect.original_protein_sequence[
                                    myEffect.aa_mutation_start_offset - difference:myEffect.aa_mutation_start_offset + (endLength)]
                                if endLength != 0:
                                    endLength = endLength - 1
                                difference = difference + 1
                            else:
                                if endLength != 0:
                                    effectDF.loc[i, mutPept + str(j)] = myEffect.mutant_protein_sequence[
                                        myEffect.aa_mutation_start_offset - difference:myEffect.aa_mutation_start_offset + (endLength)]
                                    effectDF.loc[i, wtPept + str(j)] = myEffect.original_protein_sequence[
                                        myEffect.aa_mutation_start_offset - difference:myEffect.aa_mutation_start_offset + (endLength)]
                                    difference = difference + 1
                                    endLength = endLength - 1
                                else:
                                    effectDF.loc[i, mutPept + str(j)] = ""
                                    effectDF.loc[i, wtPept + str(j)] = ""
                        elif j == 1:
                            effectDF.loc[i, mutPept + str(j)] = myEffect.mutant_protein_sequence[
                                myEffect.aa_mutation_start_offset:myEffect.aa_mutation_start_offset + (peptideLength)]
                            effectDF.loc[i, wtPept + str(j)] = myEffect.original_protein_sequence[
                                myEffect.aa_mutation_start_offset - (j - 1):myEffect.aa_mutation_start_offset + (peptideLength)]
                        else:
                            effectDF.loc[i, mutPept + str(j)] = myEffect.mutant_protein_sequence[
                                myEffect.aa_mutation_start_offset - (j - 1):myEffect.aa_mutation_start_offset + (peptideLength - (j - 1))]
                            effectDF.loc[i, wtPept + str(j)] = myEffect.original_protein_sequence[
                                myEffect.aa_mutation_start_offset - (j - 1):myEffect.aa_mutation_start_offset + (peptideLength - (j - 1))]
                    if re.search('ins', myEffect.variant.short_description):
                        if re.search('\*', myEffect.short_description):
                            effectDF.loc[i, mutPept + str(j)] = ""
                            effectDF.loc[i, wtPept + str(j)] = ""
                        elif (myEffect.aa_mutation_start_offset + peptideLength) > len(myEffect.mutant_protein_sequence):
                            if j == 1:
                                endLength = len(myEffect.mutant_protein_sequence) - myEffect.aa_mutation_start_offset
                                difference = peptideLength - endLength
                                effectDF.loc[i, mutPept + str(j)] = myEffect.mutant_protein_sequence[
                                    myEffect.aa_mutation_start_offset - difference:myEffect.aa_mutation_start_offset + (endLength)]
                                effectDF.loc[i, wtPept + str(j)] = myEffect.original_protein_sequence[
                                    myEffect.aa_mutation_start_offset - difference:myEffect.aa_mutation_start_offset + (endLength)]
                                if endLength != 0:
                                    endLength = endLength - 1
                                difference = difference + 1
                            else:
                                if endLength != 0:
                                    effectDF.loc[i, mutPept + str(j)] = myEffect.mutant_protein_sequence[
                                        myEffect.aa_mutation_start_offset - difference:myEffect.aa_mutation_start_offset + (endLength)]
                                    effectDF.loc[i, wtPept + str(j)] = myEffect.original_protein_sequence[
                                        myEffect.aa_mutation_start_offset - difference:myEffect.aa_mutation_start_offset + (endLength)]
                                    difference = difference + 1
                                    endLength = endLength - 1
                                else:
                                    effectDF.loc[i, mutPept + str(j)] = ""
                                    effectDF.loc[i, wtPept + str(j)] = ""
                        elif j == 1:
                            effectDF.loc[i, mutPept + str(j)] = myEffect.mutant_protein_sequence[
                                myEffect.aa_mutation_start_offset:myEffect.aa_mutation_start_offset + (peptideLength)]
                            effectDF.loc[i, wtPept + str(j)] = myEffect.original_protein_sequence[
                                myEffect.aa_mutation_start_offset - (j - 1):myEffect.aa_mutation_start_offset + (peptideLength)]
                        else:
                            effectDF.loc[i, mutPept + str(j)] = myEffect.mutant_protein_sequence[
                                myEffect.aa_mutation_start_offset - (j - 1):myEffect.aa_mutation_start_offset + (peptideLength - (j - 1))]
                            effectDF.loc[i, wtPept + str(j)] = myEffect.original_protein_sequence[
                                myEffect.aa_mutation_start_offset - (j - 1):myEffect.aa_mutation_start_offset + (peptideLength - (j - 1))]
                    if re.search('>', myEffect.variant.short_description):
                        if re.search('stop-loss', myEffect.short_description):
                            print("HERE YO!!!! WHY YOU NOT PRINTING")
                            print(myEffect.short_description)
                            if j == 1:
                                print(myEffect.mutant_protein_sequence[
                                      myEffect.aa_mutation_start_offset:myEffect.aa_mutation_start_offset + (peptideLength)])
                                print(mutPept + str(j))
                                effectDF.loc[i, mutPept + str(j)] = myEffect.mutant_protein_sequence[
                                    myEffect.aa_mutation_start_offset:myEffect.aa_mutation_start_offset + (peptideLength)]
                                effectDF.loc[i, wtPept + str(j)] = ""
                            else:
                                print(myEffect.mutant_protein_sequence[
                                      myEffect.aa_mutation_start_offset - (j - 1):myEffect.aa_mutation_start_offset + (peptideLength - (j - 1))])
                                print(mutPept + str(j))
                                print(myEffect.original_protein_sequence[
                                      myEffect.aa_mutation_start_offset - (j - 1):myEffect.aa_mutation_start_offset + (peptideLength - (j - 1))])
                                effectDF.loc[i, mutPept + str(j)] = myEffect.mutant_protein_sequence[
                                    myEffect.aa_mutation_start_offset - (j - 1):myEffect.aa_mutation_start_offset + (peptideLength - (j - 1))]
                                effectDF.loc[i, wtPept + str(j)] = ""
                        elif re.search('\*', myEffect.short_description):
                            effectDF.loc[i, mutPept + str(j)] = ""
                            effectDF.loc[i, wtPept + str(j)] = ""
                        elif (myEffect.aa_mutation_start_offset + peptideLength) > len(
                                myEffect.original_protein_sequence):
                            if j == 1:
                                endLength = len(myEffect.mutant_protein_sequence) - myEffect.aa_mutation_start_offset
                                difference = peptideLength - endLength
                                effectDF.loc[i, mutPept + str(j)] = myEffect.mutant_protein_sequence[
                                    myEffect.aa_mutation_start_offset - difference:myEffect.aa_mutation_start_offset + (endLength)]
                                effectDF.loc[i, wtPept + str(j)] = myEffect.original_protein_sequence[
                                    myEffect.aa_mutation_start_offset - difference:myEffect.aa_mutation_start_offset + (endLength)]
                                if endLength != 0:
                                    endLength = endLength - 1
                                difference = difference + 1
                            else:
                                if endLength != 0:
                                    effectDF.loc[i, mutPept + str(j)] = myEffect.mutant_protein_sequence[
                                        myEffect.aa_mutation_start_offset - difference:myEffect.aa_mutation_start_offset + (endLength)]
                                    effectDF.loc[i, wtPept + str(j)] = myEffect.original_protein_sequence[
                                        myEffect.aa_mutation_start_offset - difference:myEffect.aa_mutation_start_offset + (endLength)]
                                    difference = difference + 1
                                    endLength = endLength - 1
                                else:
                                    effectDF.loc[i, mutPept + str(j)] = ""
                                    effectDF.loc[i, wtPept + str(j)] = ""
                        elif (myEffect.aa_mutation_start_offset - peptideLength) < 0:
                            if j == 1:
                                startLength = myEffect.aa_mutation_start_offset
                                difference = peptideLength - startLength
                            if startLength >= 0:
                                effectDF.loc[i, mutPept + str(j)] = myEffect.mutant_protein_sequence[
                                    myEffect.aa_mutation_start_offset - (j - 1):myEffect.aa_mutation_start_offset + (peptideLength - (j - 1))]
                                effectDF.loc[i, wtPept + str(j)] = myEffect.original_protein_sequence[
                                    myEffect.aa_mutation_start_offset - (j - 1):myEffect.aa_mutation_start_offset + (peptideLength - (j - 1))]
                                difference = difference - 1
                                startLength = startLength - 1
                            else:
                                effectDF.loc[i, mutPept + str(j)] = ""
                                effectDF.loc[i, wtPept + str(j)] = ""

                        elif j == 1:
                            effectDF.loc[i, mutPept + str(j)] = myEffect.mutant_protein_sequence[
                                myEffect.aa_mutation_start_offset:myEffect.aa_mutation_start_offset + (peptideLength)]
                            effectDF.loc[i, wtPept + str(j)] = myEffect.original_protein_sequence[
                                myEffect.aa_mutation_start_offset - (j - 1):myEffect.aa_mutation_start_offset + (peptideLength)]
                        else:
                            effectDF.loc[i, mutPept + str(j)] = myEffect.mutant_protein_sequence[
                                myEffect.aa_mutation_start_offset - (j - 1):myEffect.aa_mutation_start_offset + (peptideLength - (j - 1))]
                            effectDF.loc[i, wtPept + str(j)] = myEffect.original_protein_sequence[
                                myEffect.aa_mutation_start_offset - (j - 1):myEffect.aa_mutation_start_offset + (peptideLength - (j - 1))]
        print(myEffect.short_description)
        peptDF = pd.DataFrame(effectDF.groupby(
            ["variantId", "mutPept1", "mutPept2", "mutPept3", "mutPept4", "mutPept5", "mutPept6", "mutPept7", "mutPept8", "mutPept9",
             "mutPept10", "mutPept11", "mutPept12", "mutPept13", "mutPept14", "mutPept15"])["effectId"].unique().apply('|'.join))
        peptDF = peptDF.join(pd.DataFrame(effectDF.groupby(
            ["variantId", "mutPept1", "mutPept2", "mutPept3", "mutPept4", "mutPept5", "mutPept6", "mutPept7", "mutPept8", "mutPept9",
             "mutPept10", "mutPept11", "mutPept12", "mutPept13", "mutPept14", "mutPept15"])["gene_name"].unique().apply('|'.join)))
        peptDF = peptDF.join(pd.DataFrame(effectDF.groupby(
            ["variantId", "mutPept1", "mutPept2", "mutPept3", "mutPept4", "mutPept5", "mutPept6", "mutPept7", "mutPept8", "mutPept9",
             "mutPept10", "mutPept11", "mutPept12", "mutPept13", "mutPept14", "mutPept15"])["topEffect"].any()))
        peptDF = peptDF.join(pd.DataFrame(effectDF.groupby(
            ["variantId", "mutPept1", "mutPept2", "mutPept3", "mutPept4", "mutPept5", "mutPept6", "mutPept7", "mutPept8", "mutPept9",
             "mutPept10", "mutPept11", "mutPept12", "mutPept13", "mutPept14", "mutPept15"])[
             "wtPept1", "wtPept2", "wtPept3", "wtPept4", "wtPept5", "wtPept6", "wtPept7", "wtPept8", "wtPept9",
             "wtPept10", "wtPept11", "wtPept12", "wtPept13", "wtPept14", "wtPept15"].first()))
        peptDF = peptDF.assign(sourceName=source_name)

        peptDF.to_csv(peptide_path)
        return peptide_path

