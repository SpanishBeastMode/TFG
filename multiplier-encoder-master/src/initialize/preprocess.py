from pysmt.shortcuts import *

import itertools

from .CFA_in_pysmt import FormatCFA
from .unit_propagation import unit_propagation

import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from utils import binary_to_spin


class Preprocess():
    def __init__(self, input_qbt_values, multiplier):
        self.input_qbt_values = input_qbt_values
        self.format_cfa_function = FormatCFA(multiplier.version).format_cfa_function

        self.symbols = {}
        self.formated_embeded_multiplier = self.write_embeded_multiplier_in_pysmt(multiplier.cfa_graphs, multiplier.chains)
        self.formated_input_qbt_values = self.write_inputs_in_pysmt(input_qbt_values)

    def _get_symb(self, qbt):
        if qbt in self.symbols.keys():
            symb = self.symbols[qbt]
        else:
            symb = Symbol(f'q_{qbt}', BOOL)
            self.symbols[qbt] = symb
        return symb


    def write_embeded_multiplier_in_pysmt(self, cfa_graphs, chains):

        def _write_embeded_cfa_in_pysmt(cfa_graph):
            vs = {v: self._get_symb(qbt) for v, qbt in cfa_graph.items() if not v.startswith('a')}
            return self.format_cfa_function(vs)

        def _write_embeded_chain_in_pysmt(chain):
            return And([Iff(*[self._get_symb(qbt) for qbt in edg]) for edg in chain])

        formated_multiplier_formula = [[_write_embeded_cfa_in_pysmt(cfa_graph) for cfa_graph in adder_graph] for adder_graph in cfa_graphs]

        formated_chains = {v: [_write_embeded_chain_in_pysmt(ch) for ch in chs] for v, chs in chains.items()}

        return And(itertools.chain(
                        itertools.chain(*formated_multiplier_formula), 
                        itertools.chain(*formated_chains.values())
                        )
                    )

    def write_inputs_in_pysmt(self, inputs):
        return {self._get_symb(qbt): Bool(True) if spin == 1 else Bool(False) for qbt, spin in inputs.items() }

    def preprocess(self):
        results = unit_propagation(self.formated_embeded_multiplier, self.formated_input_qbt_values)
        enhanced_input_qbt_values = {int(str(symb).split('q_')[1]): binary_to_spin(1 if x == Bool(True) else 0) for symb, x in results.items()}
        extras = set(enhanced_input_qbt_values.keys()) - set(self.input_qbt_values)
        print(f'unit propgation determined the values of extra {len(extras)} qbts')
        return enhanced_input_qbt_values