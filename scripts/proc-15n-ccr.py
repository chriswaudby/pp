#!/Users/chris/Library/Enthought/Canopy_64bit/User/bin/python
import sys
import nmrglue as ng
import numpy as np
from copy import deepcopy

propagators = {
    'IaS+ -> a' : {'B':4, 'D1':'sum','D2':'+a'},
    'IaS+ -> b' : {'B':4, 'D1':'sum','D2':'+b'},
    'IaS- -> a' : {'B':3, 'D1':'diff','D2':'-a'},
    'IaS- -> b' : {'B':3, 'D1':'diff','D2':'-b'},
    'IbS+ -> a' : {'B':2, 'D1':'diff','D2':'+a'},
    'IbS+ -> b' : {'B':2, 'D1':'diff','D2':'+b'},
    'IbS- -> a' : {'B':1, 'D1':'sum','D2':'-a'},
    'IbS- -> b' : {'B':1, 'D1':'sum','D2':'-b'}
}

#expt = str(sys.argv[1])

def main():
    process_component('IbS+ -> b','IbS- -> b',output_filename='hbb.fid')
    process_component('IaS+ -> b','IaS- -> b',output_filename='hab.fid')
    process_component('IbS+ -> a','IbS- -> a',output_filename='hba.fid')
    process_component('IaS+ -> a','IaS- -> a',output_filename='haa.fid')

    process_component('IbS+ -> b','IbS- -> b',output_filename='nbb.fid',start='N')
    process_component('IaS+ -> b','IaS- -> b',output_filename='nab.fid',start='N')
    process_component('IbS+ -> a','IbS- -> a',output_filename='nba.fid',start='N')
    process_component('IaS+ -> a','IaS- -> a',output_filename='naa.fid',start='N')



def process_component(prop1, prop2, output_filename, start='H'):
    # import the raw FID data
    #dic,data=ng.pipe.read(expt + '/test.fid')
    dic,data=ng.pipe.read('test.fid')
    udic = ng.pipe.guess_udic(dic, data)    # set the spectral parameters

    dic_echo, data_echo = isolate_propagator(udic, data, start=start, pathway=propagators[prop1])
    dic_antiecho, data_antiecho = isolate_propagator(udic, data, start=start, pathway=propagators[prop2])

    # resolve echo/anti-echo components
    output_dic = deepcopy(dic_echo)
    output_dic[0]['size'] *= 2

    output_data = np.zeros((output_dic[0]['size'],output_dic[1]['size']),dtype='complex64')
    output_data[::2] = data_echo + data_antiecho
    output_data[1::2] = 1j*(data_echo - data_antiecho)

    if start=='N' and 'Ia' in prop1:
        # reverse I(alpha) terms in N-start experiments
        output_data *= -1

    # write processed data back to nmrPipe by directly editing original dictionary
    dic['FDF1APOD'] /= 96
    dic['FDF1TDSIZE'] /= 96
    dic['FDSPECNUM'] /= 96

    #ng.pipe.write(expt + '/' + output_filename, dic, output_data, overwrite=True)
    ng.pipe.write(output_filename, dic, output_data, overwrite=True)

    # # write the processed data back to nmrPipe format (via universal dictionary)
    # C = ng.convert.converter()
    # C.from_universal(output_dic, output_data)
    # ng.pipe.write(expt + '/' + output_filename, *C.to_pipe(), overwrite=True)


def isolate_propagator(input_dic, input_data, start, pathway):
    # start = 'H' or 'N'
    # pathway = dictionary with B, D1 and D2
    # B = 1 to 4
    # D1 = 'sum' or 'diff'
    # D2 = '+a', '+b', '-a' or '-b'
    B = pathway['B']
    D1 = pathway['D1']
    D2 = pathway['D2']

    # make a copy of the input dictionary for the processed data
    dic = deepcopy(input_dic)

    # split into H-start and N-start
    dataH = input_data[::2] + input_data[1::2]
    dataN = input_data[::2] - input_data[1::2]

    # select H-start or N-start to work with...
    if start == 'H':
        data = dataH
    else:
        data = dataN
    dic[0]['size'] /= 2     # update spectrum size

    # make linear combinations of B propagators
    # dataB will become a list with the four combinations
    if B == 1:
        data = -data[0::4,:]+data[2::4,:]+1j*(data[1::4,:]+data[3::4,:])
    elif B == 2:
        data = -data[0::4,:]+data[2::4,:]-1j*(data[1::4,:]+data[3::4,:])
    elif B == 3:
        data = -data[1::4,:]+data[3::4,:]-1j*(data[0::4,:]+data[2::4,:])
    else:
        data = -data[1::4,:]+data[3::4,:]+1j*(data[0::4,:]+data[2::4,:])
    dic[0]['size'] /= 4     # update spectrum size

    # phase cycles for first S3CT propagators
    #psi1 = [0., 240., 120., 0., 240., 120.]
    #psi2 = [0., 240., 120., 180., 60., 300.]
    data1 = np.copy(data)
    data2 = np.copy(data)
    #for i in range(6):
    #    data1[i::6] *= np.exp( 1j * np.pi * psi1[i] / 180. )
    #    data2[i::6] *= np.exp( 1j * np.pi * psi2[i] / 180. )
    data2[1::2] *= -1

    if D1 == 'diff':
        data = data1 + data2
    else:
        data = data1 - data2
    data = data[0::2] + data[1::2] # + data[2::6] + data[3::6] + data[4::6] + data[5::6]
    dic[0]['size'] /= 2     # update spectrum size

    # phase cycles for second S3CT propagators
    psi1 = [0., 120., 240., 0., 120., 240.]
    psi2 = [0., 240., 120., 0., 240., 120.]
    psi3 = [0., 120., 240., 180., 300., 60.]
    psi4 = [0., 240., 120., 180., 60., 300.]
    data1 = np.copy(data)
    data2 = np.copy(data)
    data3 = np.copy(data)
    data4 = np.copy(data)
    for i in range(6):
        data1[i::6] *= np.exp( 1j * np.pi * psi1[i] / 180. )
        data2[i::6] *= np.exp( 1j * np.pi * psi2[i] / 180. )
        data3[i::6] *= np.exp( 1j * np.pi * psi3[i] / 180. )
        data4[i::6] *= np.exp( 1j * np.pi * psi4[i] / 180. )
    if D2 == '+b':
        data = data1 - data3
    elif D2 == '+a':
        data = -1 * (data1 + data3)
    elif D2 == '-b':
        data = -1 * (data2 + data4)
    elif D2 == '-a':
        data = data2 - data4
    data = data[0::6] + data[1::6] + data[2::6] + data[3::6] + data[4::6] + data[5::6]
    dic[0]['size'] /= 6     # update spectrum size

    # finally combine H start phase cycling
    if start == 'H':
        data = data[0::2] - data[1::2]
    else:
        data = data[0::2] + data[1::2]

    dic[0]['size'] /= 2     # update spectrum size


    return dic, data


main()
