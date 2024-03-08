from .multiplier import Multiplier
from .mixins import CFAMixin, EmbedMixin
from .construct_subgraph import Construct_CFAgraph_for_Multiplier_ver2, Construct_chain_for_adder_carry, construct_horizontal_chains, construct_vertical_chains

import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from utils import lp

class Multiplier_ver2_to_link_B(CFAMixin, EmbedMixin, Multiplier):
    version = '_ver2_to_link_B'
    propg_direction = {
        'A': (90, 'down'),
        'B': (45, 'down'),
        'out': (60, 'down')
    }
    shared_from = {
        'out': 'in2',
        'c_out': 'c_in'
    }
    shared_direction = {
        'carry':  propg_direction['B'], 
        'out': propg_direction['out']
    }

    def __init__(self, *args, cfa_indexs=[4, 5], library_path=lp):
        self.cfa_pair = self.choose_cfa_pair(cfa_indexs, self.version, library_path=library_path)

        self.Construct_CFAgraph = Construct_CFAgraph_for_Multiplier_ver2
        self.Construct_Adderchain = Construct_chain_for_adder_carry

        super().__init__(*args)

    def plmt_cfas(self):
        return [[self.cfa_pair[j%2] for j in range(self.length)] for _ in range(self.width)]

    def construct_chains(self):
        return {
            'carry': self.construct_chains_for_adder_carries(),
            'B': construct_horizontal_chains(self.cfa_graphs, vs=['enable']*2),
            'A': construct_vertical_chains(self.cfa_graphs)
        }
