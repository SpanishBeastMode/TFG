from unittest import TestCase, main
from parameterized import parameterized_class
import itertools

import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from src.multipliers import *

import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from src.multipliers import Multiplier_ver2, Multiplier_ver3

@parameterized_class(('version', 'length', 'width'), [
    ('_ver3_to_link_B', 21, 12),
    ('_ver2_to_link_B', 17, 8)
])
class TestMultiplier_to_link_B(TestCase):
    def setUp(self):
        multiplier = eval(f'Multiplier{self.version}')(self.length, self.width)
        multiplier.encode()
        self.multiplier = multiplier

        expected_bias = {}
        cfa_pair = multiplier.cfa_pair
        for fv, tv in multiplier.shared_from.items():
            ft_cfa_comb = {
                j: (cfa_pair[j], cfa_pair[j]) for j in [0, 1]
            } if fv == 'in1_out' else {
                j: cfa_pair if j == 1 else cfa_pair[::-1] for j in [0, 1]
            }
            expected_bias[(fv, tv)] = {i_or_j: self._expected_shared_qbt_bias([fv, tv], cfa_sequence)
                                                            for i_or_j, cfa_sequence in ft_cfa_comb.items()}
       
        self.expected_bias = expected_bias
        self.expected_chain_strength = -2

        self.shared_graph = self.multiplier.shared_graph


    def _expected_shared_qbt_bias(self, vs, cfas):
        ps = [cfa_plmt[v] for v, (cfa_plmt, _) in zip(vs, cfas)]
        biases = [cfa_ising_model['biases'][p] for p, (_, cfa_ising_model) in zip(ps, cfas)]
        return sum(biases)
       
    def _locate_shared_qbt(self, _in, sq):
        for i, adder_graph in enumerate(self.multiplier.cfa_graphs):
            for j, cfa_graph in enumerate(adder_graph):
                if sq == cfa_graph[_in]:
                    return j
        return None

    def test_encode(self):
        for _out, _in in self.multiplier.shared_from.items():
            for sq in self.shared_graph['nodes'][_in]:
                c = self.multiplier.model['biases'][sq]
                print(f'{_in}: {sq}, {c}')

                j = self._locate_shared_qbt(_in, sq)
                self.assertEqual(c, self.expected_bias[(_out, _in)][j%2])

        for sc in self.shared_graph['edges']['cfa']:
            c = self.multiplier.model['couplings'][sc]
            print(f'{sc}: {c}')
            self.assertEqual(c, self.expected_shared_coupling_strength)

        for k, chains in self.multiplier.chains.items():
            for edg in itertools.chain(*chains):
                edg = tuple(sorted(edg))
                self.assertEqual(self.multiplier.model['couplings'][edg], self.expected_chain_strength)


if __name__ == '__main__':
    main()