from unittest import TestCase, main
from parameterized import parameterized_class

import itertools

import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from src.multipliers import Multiplier_ver2, Multiplier_ver3

@parameterized_class(('version', 'length', 'width'), [
    ('_ver2', 17, 8),
    ('_ver3', 21, 12)
])
class TestMultiplier(TestCase):
    def setUp(self):
        multiplier = eval(f'Multiplier{self.version}')(self.length, self.width)
        multiplier.encode()
        self.multiplier = multiplier

        self.expected_shared_qbt_bias = {tuple(vs): self._expected_shared_qbt_bias(vs) for vs in multiplier.shared_from.items()}
        self.expected_shared_coupling_strength = {tuple(edgs): self._expected_shared_coupling_strength(edgs) for edgs in multiplier.shared_coupling.items()}
        self.expected_chain_strength = -2

        self.shared_graph = self.multiplier.shared_graph

    def _expected_shared_qbt_bias(self, vs):
        cfa_plmt, cfa_ising_model = self.multiplier.cfa
        ps = [cfa_plmt[v] for v in vs]
        return sum([cfa_ising_model['biases'][p] for p in ps])
        
    def _expected_shared_coupling_strength(self, es):
        cfa_plmt, cfa_ising_model = self.multiplier.cfa
        es = [[cfa_plmt[v] for v in vs] for vs in es]
        es = [tuple(sorted(e)) for e in es]
        coulpings = {tuple(sorted(e)): c for e, c in cfa_ising_model['couplings'].items()}
        return sum([coulpings[e] for e in es])

    def test_encode(self):
        for _out, _in in self.multiplier.shared_from.items():
            c = self.expected_shared_qbt_bias[(_out, _in)]
            for sq in self.shared_graph['nodes'][_in]:
                self.assertEqual(self.multiplier.model['biases'][sq], c)

        for edgs in self.multiplier.shared_coupling.items():
            c = self.expected_shared_coupling_strength[tuple(edgs)]
            for sc in self.shared_graph['edges']['cfa']:
                self.assertEqual(self.multiplier.model['couplings'][sc], c)

        for k, chains in self.multiplier.chains.items():
            for edg in itertools.chain(*chains):
                self.assertEqual(self.multiplier.model['couplings'][tuple(sorted(edg))], self.expected_chain_strength)

    def test_compatibility_with_hardware_topology(self):
        if self.version == '_ver2':
            compatiable = self.multiplier.is_compatiable_with_real_graph(self.multiplier.model)
            self.assertTrue(compatiable)


if __name__ == '__main__':
    main()