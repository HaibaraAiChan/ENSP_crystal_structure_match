



import argparse




def getArgs():
    parser = argparse.ArgumentParser('python')
    parser.add_argument('-input',
                        default='../../SNPs_with_crystal_structure',
                        required=False,
                        help='the path of input file')
    parser.add_argument('-output',
                        default='./input/',
                        required=False,
                        help='the result stored in ')
    parser.add_argument('-numLst',
                        default=64,
                        required=False,
                        help='the num of list')
    parser.add_argument('-res',
                        default='input.lst',
                        required=False,
                        help='the input feed to pbs file')
    return parser.parse_args()



def gen_list(file, out, num,res):
    file_1 = open(file, 'r')
    context_1 = file_1.read()
    file_1.close()

    list_1 = context_1.split('\n')[1:]
    title = context_1.split('\n')[0]

    list_1 = filter(None, list_1)
    size_1 = len(list_1)
    div = size_1 / num
    tail = size_1 % num
    if tail != 0:
        div = size_1/(num-1)
        tail = size_1 % (num-1)
    map_1 = [list_1[i:i + div] for i in range(0, len(list_1) - tail, div)]
    if tail != 0:
        map_1.append(list(list_1[-tail:]))

    print len(map_1)

    for i in range(num):
        f = open(out+str(i), 'w')
        str_w = title+'\n'+'\n'.join(map_1[i])
        f.write(str_w)
        f.close()
    f = open(res, 'w')
    lst = [i for i in range(num)]
    str_w = '\n'.join([out+str(i) for i in lst])
    f.write(str_w)
    f.close()


if __name__ == '__main__':
    args =getArgs()
    origin = args.input
    output = args.output
    num_list = args.numLst
    res = args.res

    gen_list(origin, output, num_list, res)
