from .CFA import *
from .CFAEncoder import CFAEncoder
from .StrictCFAEncoder import StrictCFAEncoder


import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from utils import lp, lp_via_z3

if __name__ == '__main__':
    version = '_ver2'
    from multipliers.multiplier_ver2 import Multiplier_ver2
    multiplier = Multiplier_ver2(3, 3,cfa_index=0, library_path=lp)

    strict_cfa = multiplier.cfa
    cfa_plmt, cfa_ising_model = strict_cfa
    cfa_graph = {
        'nodes': cfa_ising_model['biases'].values(),
        'edges': cfa_ising_model['couplings'].keys()
    }
    f, n_x, n_a = eval(f'CFA{version}'), 7, 4

    # encoder = CFAEncoder(f, n_x, n_a, cfa_plmt, cfa_graph)

    encoder = StrictCFAEncoder(version, f, n_x, n_a, cfa_plmt, cfa_graph)
    encoder.run_solver()



    