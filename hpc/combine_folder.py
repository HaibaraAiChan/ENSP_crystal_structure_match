import os
import numpy as np


def combine(folder, out):
    out_str_list = []
    for f in os.listdir(folder):
        file = open(folder + f, 'r')
        context = file.read()
        out_str_list.append(context)
    out_str_list = np.sort(out_str_list)
    output = ''.join([out_str_list[i] for i in range(len(out_str_list))])
    f_out = open(out, 'w')
    f_out.write(output)
    f_out.close()


if __name__ == '__main__':
    folder = '../output/'
    out = 'map_output_hpc.txt'
    combine(folder, out)
