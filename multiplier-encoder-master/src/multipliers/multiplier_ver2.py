from .multiplier import Multiplier
from .mixins import CFAMixin, EmbedMixin
from .construct_subgraph import Construct_CFAgraph_for_Multiplier_ver2, Construct_chain_for_adder_carry, construct_vertical_chains

import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from utils import lp

class Multiplier_ver2(CFAMixin, EmbedMixin, Multiplier):
    version = '_ver2'
    propg_direction = {
        'A': (90, 'down'),
        'B': (45, 'down'),
        'out': (60, 'down')
    }
    shared_from = {
        'out': 'in2',
        'c_out': 'c_in',
        'enable_out': 'enable'
    }
    shared_direction = {
        'carry': propg_direction['B'], # also for 'enable'
        'out': propg_direction['out']
    }
    shared_coupling = {
         ('c_out', 'enable_out'): ('c_in', 'enable')
    }

    def __init__(self, *args, cfa_index=0, library_path=lp, **kwargs):
        self.cfa = self.choose_cfa(cfa_index, self.version, library_path=library_path)
        self.cfa_index = cfa_index

        self.Construct_CFAgraph = Construct_CFAgraph_for_Multiplier_ver2
        self.Construct_Adderchain = Construct_chain_for_adder_carry

        super().__init__(*args, **kwargs)


    def plmt_cfas(self):
        return [[self.cfa for _ in range(self.length)] for _ in range(self.width) ]

    def construct_chains(self):
        return {
            'carry': self.construct_chains_for_adder_carries(),
            'A': construct_vertical_chains(self.cfa_graphs, v='in1')
        }

