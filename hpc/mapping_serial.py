from __future__ import division
import pandas as pd
import numpy as np
#from Bio.pairwise2 import format_alignment
import Bio
from Bio.SubsMat.MatrixInfo import blosum62
import re
from Bio import pairwise2
import time
import argparse




def getArgs():
    parser = argparse.ArgumentParser('python')
    parser.add_argument('-folderPath',
                        default='/work/syang29/Monsif_pbs/',
                        required=False,
                        help='the path of folder')
    parser.add_argument('-inputCSV',
                        default='./SNPs_with_crystal_structure',
                        required=True,
                        help='the lst of input')
    parser.add_argument('-out',
                        default='./pdb_aligning_results',
                        required=True,
                        help='the folder of output')

    parser.add_argument('-pdb',
                        default='../Monsif/pdb_only/',
                        required=False,
                        help='the folder of pdb')
    parser.add_argument('-fasta',
                        default='../Monsif/fasta/',
                        required=False,
                        help='the folder of fasta')
    parser.add_argument('-test',
                        default='../Monsif/test/',
                        required=False,
                        help='the folder of test')
    return parser.parse_args()



def read_file(filename):
    f = open(filename, 'r')
    list_ = [i.strip("\n").strip(" ") for i in list(f)]
    f.close()
    return list_


def count_offset_2(seq_x, pos_x,seq_y):

    assert len(seq_x) == len(seq_y)
    seq_x_list = re.findall(r'-+[A-Z]+|^[A-Z]+', seq_x)  # '--------xxxxxxx' or 'xxxxxxxx' list
    seq_y_list = re.findall(r'-+[A-Z]+|^[A-Z]+', seq_y)  # '--------xxxxxxx' or 'xxxxxxxx' list
    len_seq_x = len(seq_x_list)
    len_seq_y = len(seq_y_list)

    length_x = [len(seq_x_list[i]) for i in range(len_seq_x)]
    length_acc_x = [sum(length_x[0:i]) for i in range(0, len_seq_x+1)]
    offset_x = [seq_x_list[i].count('-') for i in range(len_seq_x)]
    pure_length_x = [x - y for x, y in zip(length_x, offset_x)]
    pure_acc_x = [sum(pure_length_x[0:i]) for i in range(0, len_seq_x+1)]

    length_y = [len(seq_y_list[i]) for i in range(len_seq_y)]
    length_acc_y = [sum(length_y[0:i]) for i in range(0, len_seq_y + 1)]
    offset_y = [seq_y_list[i].count('-') for i in range(len_seq_y)]
    pure_length_y = [x - y for x, y in zip(length_y, offset_y)]
    pure_acc_y = [sum(pure_length_y[0:i]) for i in range(0, len_seq_y + 1)]

    pos_x_total = 0
    ###############################################################################################
    for i in range(0, len_seq_x):

        if pure_acc_x[i] < pos_x:
            pos_x_total = pos_x_total + offset_x[i]
        else:
            break
    pos_x_total = pos_x_total + pos_x

    _length = 0
    for i in range(0, len_seq_y):

        if length_acc_y[i] < pos_x_total:
            _length = _length + offset_y[i]

        else:
            break
    # off_y = pos_x_total - _length
    return pos_x_total , _length


def mapping(input_csv, pdb_path, seq_path, total_seq_path, output):
    df = pd.read_csv(input_csv, delim_whitespace=True)
    f_out = open(output, 'w')
    temp_zip = zip(df.wild_type, df.PROTEIN_ID, df.Position, df.Mutant, df.SNP_ID,
                   df.disorder, df.confidence, df.wt_codon, df.mutant_codon)
    temp_zip = sorted(temp_zip, key=lambda x: (x[1], len(x[4])))
    for wt, protein, position, mutant, snp, disorder, confidence, wt_codon, mutant_codon in temp_zip:

        struct_list = read_file(pdb_path + '%s.out' % protein)
        struct_list = np.sort(struct_list)

        seq_list = read_file(seq_path + '%s.seq' % protein)
        x = ''.join(seq_list)
        x_len = len(x)
        position_x = position
        for struct in struct_list:
            y1 = read_file(total_seq_path + '/%s.fasta' % struct)
            y = ''.join(y1[1:])
            y_len = len(y)

            alignment = pairwise2.align.globalds(x, y, blosum62, -12, -1.2)
            for i in range(len(alignment)):
                if i >= 1:
                    break
                pdb_x_seq = alignment[i][0]
                pure_x_seq = ''.join(re.findall(r'[A-Z]+', pdb_x_seq))

                pdb_y_seq = alignment[i][1]
                pure_y_seq = ''.join(re.findall(r'[A-Z]+', pdb_y_seq))

                if x_len != len(pure_x_seq):
                    continue
                if y_len != len(pure_y_seq):
                    continue

                pos_y_t, offset_y = count_offset_2(pdb_x_seq, position_x, pdb_y_seq)
                if pdb_y_seq[pos_y_t - 1] == '-':
                    continue

                position_y = pos_y_t - offset_y
                w_str = '{}\t{:>10}\t{}\t{}\t{}\t{}\t{:<5}\t{:<3}{:>5}\t{:>5}{:>5}{:>3}\n'. \
                    format(protein, snp, wt_codon, mutant_codon, wt,  mutant, str(position_x), disorder,str(confidence), \
                           struct,  str(position_y), pure_y_seq[position_y - 1])
                
                f_out.write(w_str)
    f_out.close()


if __name__ == "__main__":

    args = getArgs()

    input_csv = args.inputCSV

    output = args.out
    pdb_path = args.pdb
    seq_path = args.test
    total_seq_path = args.fasta

    # input_csv = './SNPs_with_crystal_structure'
    #
    # # input_csv = './compare_output/left_struct'
    # output = './pdb_aligning_results'
    # pdb_path = './pdb_only/'
    # seq_path = './test/'
    # total_seq_path = './fasta/'
    start = time.time()

    mapping(input_csv, pdb_path, seq_path, total_seq_path, output)
    end = time.time()
    print'time elapsed :'+str(end - start)
