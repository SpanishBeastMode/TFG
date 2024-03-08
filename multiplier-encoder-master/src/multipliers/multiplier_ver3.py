from .multiplier import Multiplier
from .mixins import CFAMixin, EmbedMixin
from .construct_subgraph import Construct_CFAgraph_for_Multiplier_ver3, Construct_chain_for_adder_carry
from .cell_utils import origin_4_max_rect

import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from utils import lp


class Multiplier_ver3(CFAMixin, EmbedMixin, Multiplier):
    version='_ver3'
    propg_direction = {
        'A': (60, 'down'),
        'B': (45, 'down'),
        'out': (30, 'down')
    }
    shared_from = {
        'out': 'in2',
        'c_out': 'c_in',
        'enable_out': 'enable',
        'in1_out': 'in1'
    }
    shared_direction = {
        'carry': propg_direction['B'], # also for 'enable'
        'out': propg_direction['out'],
        'A': propg_direction['A']
    }
    shared_coupling = {
         ('c_out', 'enable_out'): ('c_in', 'enable')
    }
    def __init__(self, *args, cfa_index=5, library_path=lp, origin_cellwrpr=origin_4_max_rect()):
        self.cfa = self.choose_cfa(cfa_index, self.version, library_path=library_path)
        self.cfa_index = cfa_index

        self.Construct_CFAgraph = Construct_CFAgraph_for_Multiplier_ver3
        self.Construct_Adderchain = Construct_chain_for_adder_carry

        super().__init__(*args, origin_cellwrpr=origin_cellwrpr)

    def plmt_cfas(self):
        return [[self.cfa for _ in range(self.length)] for _ in range(self.width)]

    def construct_chains(self):
        return {'carry': self.construct_chains_for_adder_carries()}
