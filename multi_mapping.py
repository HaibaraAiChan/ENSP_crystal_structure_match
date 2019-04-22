from __future__ import division
import pandas as pd
import numpy as np
import zipfile

from Bio.SubsMat.MatrixInfo import blosum62
import re
from Bio import pairwise2
import time
from multiprocessing import Pool


def read_file(filename):
    f = open(filename, 'r')
    list_ = [i.strip("\n").strip(" ") for i in list(f)]
    f.close()
    return list_


def count_offset_2(seq_x, pos_x, seq_y):

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
    off_y = 0
    _length = 0
    for i in range(0, len_seq_y):

        if length_acc_y[i] < pos_x_total:
            _length = _length + offset_y[i]

        else:
            break
    off_y = pos_x_total - _length
    return off_y


def worker(temp_zip, output):
    f_out = open(output, 'w')
    w_str = ''
    for wt, protein, position, mutant, snp, disorder, confidence, wt_codon, mutant_codon in temp_zip:
        struct_list = read_file(pdb_path + '%s.out' % protein)
        struct_list = np.sort(struct_list)

        seq_list = read_file(seq_path + '%s.seq' % protein)
        x = ''.join(seq_list)
        x_len = len(x)
        position_x = position
        w_str = ''
        for struct in struct_list:
            y1 = read_file(total_seq_path + '/%s.fasta' % struct)
            y = ''.join(y1[1:])
            y_len = len(y)

            alignment = pairwise2.align.globalds(x, y, blosum62, -12, -1.2)
            pdb_x_seq = alignment[0][0]
            pure_x_seq = ''.join(re.findall(r'[A-Z]+', pdb_x_seq))

            pdb_y_seq = alignment[0][1]
            pure_y_seq = ''.join(re.findall(r'[A-Z]+', pdb_y_seq))

            if x_len != len(pure_x_seq):
                continue

            if pdb_y_seq[position - 1] == '-':
                continue
            position_y = count_offset_2(pdb_x_seq, position_x, pdb_y_seq)
            t_w_str = protein + '\t' \
                  + '{:>10}'.format(snp) \
                  + '\t' + wt_codon \
                  + '\t' + mutant_codon \
                  + '\t' + mutant + '\t' + disorder + '\t' \
                  + '{:<5}'.format(str(confidence)) + '\t' \
                  + '{:>5}'.format(struct) + '\t' \
                  + '{:>5}'.format(str(position_x)) \
                  + '{:>3}'.format(wt) \
                  + '{:>5}'.format(str(position_y)) \
                  + '{:>3}'.format(pure_y_seq[position_y - 1])
            w_str = w_str + t_w_str

    f_out.write(w_str)

    f_out.close()


def mapping(input_csv, pdb_path, seq_path, total_seq_path, output):
    df = pd.read_csv(input_csv, delim_whitespace=True)

    temp_zip = zip(df.wild_type, df.PROTEIN_ID, df.Position, df.Mutant, df.SNP_ID,
                   df.disorder, df.confidence, df.wt_codon, df.mutant_codon)
    temp_zip = sorted(temp_zip, key=lambda x: (x[1], len(x[4])))
    size = len(temp_zip)

    P_NUM = 30
    p = Pool(P_NUM)

    for i in range(P_NUM):
        p.apply_async(worker,
                      args=(temp_zip[int(size / P_NUM) * i: int(size / P_NUM) * (i + 1)], output + '_' + str(i)))
    p.close()
    p.join()


if __name__ == "__main__":
    input_csv = './SNPs_with_crystal_structure'
    output = './output/pdb_aligning_results'
    pdb_path = './pdb_only/'
    seq_path = './test/'
    total_seq_path = './fasta/'

    start = time.time()

    mapping(input_csv, pdb_path, seq_path, total_seq_path, output)
    end = time.time()
    print'time elapsed :' + str(end - start)
