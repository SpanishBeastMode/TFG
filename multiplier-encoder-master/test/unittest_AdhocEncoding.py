import copy
from unittest import TestCase, main
from parameterized import parameterized_class, parameterized

from itertools import product
from collections import namedtuple

import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from src.generate_inputs import generate_inputs
from src.multipliers import Multiplier_ver2
from src.initialize import Initializer, InputCFALibrary

GapOffset = namedtuple('GapOffset', ['g_min', 'offset'])

@parameterized_class(('length', 'width'), [
    (3, 3),
    (3, 4),
    (5, 5),
    (17, 8)
])
class TestCase_AdhocEncoding(TestCase):
    via = 'adhoc_encoding'
    length, width = 3, 4
    version = '_ver2'
    cfa_index = 0
    def setUp(self):
        self.assertEqual(self.cfa_index, 0)

        length, width = self.length, self.width
        P = generate_inputs(length, width, counter=10)
        self.inputs_list = ({'out': p} for p in P)

    @parameterized.expand([
        # (True),
        (False)
    ])
    def test_adhoc_InputCFA_mapping(self, ancillar_enhanced):
        inputs = list(self.inputs_list)[0]

        initializer = Initializer(inputs, self.length, self.width, via=self.via, using_improved_library=False)

        # initializer = Initializer(inputs, self.multiplier, ancillar_enhanced)
        adhoc_cfa_plmts = initializer._adhoc_InputCFA_mapping(ancillar_enhanced=ancillar_enhanced)

        extracted_GapOffets = [[GapOffset(*[adhoc_cfa_plmts[i][j][1][key] for key in ['g_min', 'offset']]) for j in range(self.multiplier.length)] for i in range(self.multiplier.width)]

        # indirest_test = False
        # if indirest_test:
        #     multiplier = initializer.multiplier
        #     multiplier.embed(adhoc_cfa_plmts=adhoc_cfa_plmts)
        #     model = multiplier.encode(adhoc_cfa_plmts=adhoc_cfa_plmts)
        #     # multiplier.draw(self.multiplier.embedding_graph, model)
        #     compatiable = self.multiplier.is_compatiable_with_real_graph(model)
        #     self.assertTrue(compatiable)

    @parameterized.expand([
        (True),
        (False)
    ])
    def test_coefficient_compatibility(self, ancillar_enhanced):
        for inputs in self.inputs_list:
            initializer = Initializer(inputs, self.length, self.width, via=self.via, using_improved_library=False)
            model = initializer.run()
            
            compatiable = initializer.multiplier.is_compatiable_with_real_graph(model)
            self.assertTrue(compatiable)


if __name__ == '__main__':
    main()