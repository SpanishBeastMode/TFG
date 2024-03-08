from unittest import TestCase, main
from parameterized import parameterized_class

import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from src.cfa import *
from src.multipliers import *
from src.utils import lp_via_z3, lp, improved_lp

@parameterized_class(('suffix', 'version', 'library_path', 'cfa_index'), [
    ('CFA', '_ver2', lp, 0),
    ('CFA', '_ver2', improved_lp, 1),
    ('CFA', '_ver3', lp, 0),
    ('StrictCFA', '_ver2', lp_via_z3, 0)
])
class TestCase_StrictCFAEncoder(TestCase):
    suffix = 'StrictCFA'
    version = '_ver2'
    library_path = lp_via_z3
    cfa_index = 0

    def setUp(self):
        cfa_index = self.cfa_index
        strict_cfa = read_CFA(cfa_index, f'{self.suffix}{self.version}', self.library_path)

        if self.library_path != lp_via_z3:
            strict_cfa = CFAMixin()._parse_cfa(strict_cfa)
            cfa_plmt, cfa_ising_model = strict_cfa 
            ordered_biases = {p: cfa_ising_model['biases'][p] for p in cfa_plmt.values()}
            cfa_ising_model['biases'] = ordered_biases
            strict_cfa = (cfa_plmt, cfa_ising_model)


        cfa_plmt, cfa_ising_model = strict_cfa
        cfa_graph = {
            'nodes': cfa_ising_model['biases'].values(),
            'edges': cfa_ising_model['couplings'].keys()
        }
        f, n_x, n_a = eval(f'CFA{self.version}'), 7, 4
        n_x = 7 if self.version == '_ver2' else 8 if self.version == '_ver3' else None
        n_a = 4

        self.cfa_ising_model = cfa_ising_model
        self.encoder = StrictCFAEncoder(self.version, f, n_x, n_a, cfa_plmt, cfa_graph)

    def test_check_cfa_logic(self):
        logic = self.encoder.verify(self.cfa_ising_model)
        self.assertTrue(logic)
        # self.assertGreater(min_gap, 0)


class TestCase_4_Multiplier_ver3_to_link_B_CFA(TestCase):
    version = '_ver3_to_link_B'
    lp = lp
    def setUp(self):
        multiplier = eval(f'Multiplier{self.version}')(3, 3)
        _, strict_cfa = multiplier.cfa_pair

        cfa_plmt, cfa_ising_model = strict_cfa
        cfa_graph = {
            'nodes': cfa_ising_model['biases'].values(),
            'edges': cfa_ising_model['couplings'].keys()
        }
        f, n_x, n_a = eval(f'CFA{self.version}'), 7, 4

        self.cfa_ising_model = cfa_ising_model
        self.encoder = StrictCFAEncoder(self.version, f, n_x, n_a, cfa_plmt, cfa_graph)


    def test_check_cfa_library(self):
        logic = self.encoder.verify(self.cfa_ising_model)
        self.assertTrue(logic)
        # self.assertGreater(min_gap, 0)


if __name__ == '__main__':
    main()
