import numpy as np

def strrep(x):
    return '\t'.join('{:.2f}'.format(y) for y in x.flatten())
#    return '\t'.join(str(y) for y in x.flatten())

def pprint_array(x):
    for row in x:
        print(strrep(row))

