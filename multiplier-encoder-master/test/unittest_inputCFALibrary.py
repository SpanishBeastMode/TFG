from unittest import TestCase, main
from parameterized import parameterized_class

import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from src.cfa import *
from src.multipliers import *
from src.initialize import *


@parameterized_class(('ancillar_enhanced', 'coord'), [
    (False, (3, 1)),
    # (True, (3, 1)) 
])
class TestCase_InputCFAEncoder(TestCase):
    def setUp(self):
        self.inputCFALibrary = InputCFALibrary(ancillar_enhanced=self.ancillar_enhanced)


    def test_InputCFALibrary_logics(self):
        def _verify_logics_in_a_library(alibrary, coord_inputs, strict):
            for cfa, coord_input in zip(alibrary, coord_inputs):
                inputCFAEncoder = InputCFAEncoder(coord_input, *self.inputCFALibrary.args, ancillar_enhanced=self.ancillar_enhanced,
                                                  strict=strict)

                cfa_plmt, cfa_ising_model = cfa
                logic = inputCFAEncoder.verify(cfa_ising_model)
                if not logic: 
                    return logic
            return True

        for coord, coord_inputs in self.inputCFALibrary.interface_cfa_inputs.items():
            clibrary = self.inputCFALibrary.Library[coord]
            strict = self.inputCFALibrary.strictness[coord]
            # if self.ancillar_enhanced and coord == (3, 1):
            #     pass
            if type(clibrary[0]) == tuple:
                logic = _verify_logics_in_a_library(clibrary, coord_inputs, strict)
                self.assertTrue(logic)
            elif type(clibrary[0]) == list:
                for alibrary in clibrary:
                    logic = _verify_logics_in_a_library(alibrary, coord_inputs, strict)
                    self.assertTrue(logic)


if __name__ == '__main__':
    main()
