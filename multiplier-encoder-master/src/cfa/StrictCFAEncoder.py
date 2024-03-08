from .CFAEncoder import CFAEncoder
from .cfa_truth_table import *

import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from utils import lp_via_z3

class StrictCFAEncoder(CFAEncoder):
    _qbt_sharing_options = {
        '_ver2': {
            'out': 'in2',
            'c_out': 'c_in',
            'enable_out': 'enable'
        },
        '_ver3': {
            'out': 'in2',
            'c_out': 'c_in',
            'enable_out': 'enable',
            'in1_out': 'in1'
        },
        '_ver3_to_link_B': {
            'out': 'in2',
            'c_out': 'c_in',
            'in1_out': 'in1'
        },
        '_ver2_to_link_B': {
            'out': 'in2',
            'c_out': 'c_in'
        },
        '_ver1_to_link_B': {}
    }
    _coupling_sharing_options = {
        version: {('c_out', 'enable_out'): ('c_in', 'enable')} 
                    for version, shared_from in _qbt_sharing_options.items()
                        if 'c_out' in shared_from.keys() and 'enable_out' in shared_from.keys()
    }

    def __init__(self, qbt_sharing_option, *args, strict=True, **kwargs):
        super().__init__(*args, **kwargs)

        self.strict = strict
        self.qbt_sharing = self.__class__._qbt_sharing_options[qbt_sharing_option]
        self.coupling_sharing = self.__class__._coupling_sharing_options[qbt_sharing_option] if qbt_sharing_option in self.__class__._coupling_sharing_options.keys() else {}

        if hasattr(self, 'inputs'):
            self.qbt_sharing = {f: t for f, t in self.qbt_sharing.items() if f not in self.inputs.keys() and t not in self.inputs.keys() }
            self.coupling_sharing = {f: t for f, t in self.coupling_sharing.items() if not {*f, *t}.intersection(set(self.inputs.keys()))}

        if self.strict:
            self.add_extra_constraints_for_shared_coeffients()

        self.CFA_truth_table = eval(f'get_CFA{qbt_sharing_option}_truth_table')()
        self.qbt_sharing_option = qbt_sharing_option

        self.library_path = os.path.join(lp_via_z3, f'StrictCFA{self.qbt_sharing_option}')

    def _get_bias(self, v, H):
        return H[self.cfa_plmt[v]]

    def _get_coupling(self, vs, J):
        edg = tuple([self.cfa_plmt[v] for v in vs])
        return J[edg] if edg in J.keys() else J[tuple(edg[::-1])]

    def add_extra_constraints_for_shared_coeffients(self):
        # if hasattr(self, 'inputs'):         # without strict constraints
        #     return

        self.constraints += [self.bias_constriant(sum([self._get_bias(v, self.Hdict) for v in vs])) 
                                    for vs in self.qbt_sharing.items()]

        if self.coupling_sharing:
            self.constraints += [self.coupling_constriant(sum([self._get_coupling(edg, self.Jdict) for edg in edgs]))
                                        for edgs in self.coupling_sharing.items()]

        print(f'\t strict constraints: {self.constraints}')


    def verify(self, ising_model):
        # if hasattr(self, 'inputs'):         # without strict constraints
        #     return

        satisfied = super().verify(ising_model)

        if not self.strict:
            return satisfied
        
        shared_qbt_biases = [[self._get_bias(v, ising_model['biases']) for v in vs]
                                    for vs in self.qbt_sharing.items()]
        satisfied = satisfied and all([sum(bs) >= self.h_range[0] and sum(bs) <= self.h_range[1] for bs in shared_qbt_biases])

        if self.coupling_sharing:
            shared_coupling_strengths = [[self._get_coupling(edg, ising_model['couplings']) for edg in edgs]
                                                for edgs in self.coupling_sharing.items()]
            satisfied = satisfied and all([sum(cs) >= self.J_range[0] and sum(cs) <=self.J_range[1] for cs in shared_coupling_strengths])
        return satisfied