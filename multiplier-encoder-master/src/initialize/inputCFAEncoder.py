from z3 import Sum

import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from utils import lp_via_z3
from cfa import StrictCFAEncoder

class InputCFAEncoder(StrictCFAEncoder):

    def __init__(self, inputs, *arg, ancillar_enhanced=True, version='_ver2', **kwargs):
        self.inputs = inputs

        self.old_f, n_x, n_a, cfa_plmt, cfa_graph = arg
        self.old_cfa_plmt = cfa_plmt

        self.n_a, self.cfa_plmt = self.reduced_cfa_plmt(cfa_plmt, inputs, n_a, ancillar_enhanced)
        self.cfa_graph = cfa_graph if ancillar_enhanced else self.reduced_cfa_graph(cfa_plmt, inputs, cfa_graph) 
        self.f = self.reduced_CFA_function

        args = (self.f, n_x - len(inputs.items()), self.n_a, self.cfa_plmt, self.cfa_graph)
        self.ancillar_enhanced = ancillar_enhanced
        
        self.strict = False if inputs else True

        super().__init__(version, *args, **kwargs)
        if ancillar_enhanced:
            # self.CFA_truth_table = [ for in2, in1, enable, carry_in, carry_out, out in self.CFA_truth_table
            #                             if 'in2' in inputs.keys()]
            pass

        self.library_path = os.path.join(lp_via_z3, f'InputCFA{self.qbt_sharing_option}')


    def reduced_cfa_plmt(self, cfa_plmt, inputs, n_a, ancillar_enhanced):
        reduced_cfa_plmt = {v: p for v, p in cfa_plmt.items() if v not in inputs.keys()}
        if not ancillar_enhanced:
            return n_a, reduced_cfa_plmt
        
        new_ancillar_plmt = {f'a{n_a + i}': cfa_plmt[v] for i, v in enumerate(inputs.keys())}
        return n_a + len(inputs.items()), {**reduced_cfa_plmt, **new_ancillar_plmt}


    def reduced_cfa_graph(self, cfa_plmt, inputs, cfa_graph):
        return {'edges': [edg for edg in cfa_graph['edges'] if all([p not in [cfa_plmt[v] for v in inputs.keys()] for p in edg])]}
    
    def reduced_CFA_function(self, reduced_xs):
        xdict = {**reduced_xs, **self.inputs}
        return self.old_f(xdict)


    def _get_coefficient(self, vs):
        return self._get_coupling(vs, self.Jdict) if type(vs) == tuple\
                    else self._get_bias(vs, self.Hdict) if type(vs) == str \
                        else vs


    def _parse_specs(self, specs):
        sums = [Sum([self._get_coefficient(param) for param in item]) for item in specs.items()]
        constraints = [self.coupling_constriant(s) if type(param) == tuple else self.bias_constriant(s) 
                            for (param, _), s in zip(specs.items(), sums)]
        print(f'\t mutural_dependent_specs: \n\t\t\t\t {specs} \n\t\t\t\t {constraints}')
        return constraints

    def impose_extra_dependency_contraints(self, mutural_dependent_specs={}):
        # params = [mutural_dependent_specs] if self.ancillar_enhanced else [self.inputs, mutural_dependent_specs]
        # for specs in params:
        self.constraints += self._parse_specs(mutural_dependent_specs)
