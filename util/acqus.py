#!/usr/bin/env python


# parse an entry and split array as required
def __parse(dat):
    if dat[0]=='(':
        # array data
        dat = dat.split(')\n')[1].split()
    return dat

# load file into dictionary
def __load_dic(filename):
    f = open(filename,'r')
    dat = f.read()
    f.close()

    fields = [x.strip().strip('$') for x in dat.split('##')]
    dic = {}
    for field in fields:
        if field=='': continue
        x = field.split('= ')
        if len(x)>1:
            dic[x[0]] = __parse(x[1])
    return dic

# get acqus entry
def get_entry(filename, parameter_name, index=None):
    dic = __load_dic(filename)
    if index==None:
        return dic[parameter_name]
    else:
        return dic[parameter_name][index]

if __name__ == "__main__":
    import argparse
    # parse command line arguments
    parser = argparse.ArgumentParser(description='Extract a parameter from bruker acqus format file')
    parser.add_argument('-f','--file',default='acqus',help='input file')
    parser.add_argument('name',help='parameter name')
    parser.add_argument('N',type=int,help='index for arrayed parameters',nargs='?')
    args = parser.parse_args()
    filename = args.file
    parameter_name = args.name
    N = args.N

    print(get_entry(filename, parameter_name, N))
