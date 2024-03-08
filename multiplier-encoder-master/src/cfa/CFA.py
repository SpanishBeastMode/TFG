import pickle

import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from utils import lp, lp_via_z3

def read_CFA(i, fdir, library_path=lp):
    fp = os.path.join(library_path, fdir, f'CFA{i}')
    with open(fp, 'rb') as file:
        cfa = pickle.load(file)
    return cfa

def controlled_Fulladder(xdict):
    controlled_in1 = (xdict['in1'] and xdict['enable'])
    in2, c_in = xdict['in2'], xdict['c_in']
    return (xdict['c_out'] == ((c_in and (controlled_in1 or in2)) or (controlled_in1 and in2)))\
        and (xdict['out'] == (controlled_in1 ^ in2 ^ c_in))


def CFA_ver2(xdict):
    return controlled_Fulladder(xdict) and (xdict['enable'] == xdict['enable_out'])

def CFA_ver3(xdict):
    return CFA_ver2(xdict) and (xdict['in1'] == xdict['in1_out'])

def CFA_ver3_to_link_B(xdict):
    return controlled_Fulladder(xdict) and (xdict['in1'] == xdict['in1_out'])