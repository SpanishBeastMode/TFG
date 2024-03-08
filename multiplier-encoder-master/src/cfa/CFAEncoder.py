from collections import defaultdict
from datetime import datetime
from itertools import product
import json
import pickle

from dimod import BinaryQuadraticModel
from z3 import *

from .script import serialize_bqm
from .cfa_truth_table import get_controlledFA_truth_table

import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import utils



class CFAEncoder():
    def __init__(self, f, n_x, n_a, cfa_plmt, cfa_graph, 
                        h_range=[-4, 4], J_range=[-2, 1]):
        self.f = f
        self.n_x, self.n_a = n_x, n_a
        self.cfa_plmt, self.cfa_graph = cfa_plmt, cfa_graph
        self.h_range, self.J_range = h_range, J_range

        self.gap, self.offset, self.Hdict, self.Jdict = self.configure_variables_for_OMT_solving()
        self.constraints = []

        self.CFA_truth_table = get_controlledFA_truth_table

        self.library_path = os.path.join(utils.lp_via_z3, 'PureCFA')

    def bias_constriant(self, h):
        return And(h >= self.h_range[0], h <= self.h_range[1])

    def coupling_constriant(self, j):
        return And(j >= self.J_range[0], j <= self.J_range[1])

    def configure_variables_for_OMT_solving(self):
        gap = Real('gap')
        offset = Real('offset')
        
        Hdict = {p: Real(f'H_{p}') for p in self.cfa_plmt.values()}
        Jdict = {(p0, p1): Real(f'J_{p0}_{p1}') for p0, p1 in self.cfa_graph['edges']}

        return gap, offset, Hdict, Jdict
    
    def impose_constraints_for_CFA_encoding(self):
        def _calcu_energy(pz, offset, Hdict, Jdict):
            linear = [h*pz[p] for p, h in Hdict.items()]
            quadratic = [j*pz[p0]*pz[p1] for (p0, p1), j in Jdict.items()]
            return Sum(linear) + Sum(quadratic) + offset

        energies = {True: defaultdict(list), False: defaultdict(list)}
        for zs in product([-1, 1], repeat=self.n_x+self.n_a):
            pz = dict(zip(self.Hdict.keys(), zs))
            e = _calcu_energy(pz, self.offset, self.Hdict, self.Jdict)
            xdict = {v: utils.spin_to_binary(pz[p]) for v, p in self.cfa_plmt.items() if not v.startswith('a')}
            energies[self.f(xdict)][tuple(xdict.values())].append(e)

        constraints = [self.gap >= 0]
        constraints += [And(And([e >= 0 for e in E]), Or([e == 0 for e in E])) for x, E in energies[True].items()]
        constraints += [And([e >= self.gap for e in E]) for x, E in energies[False].items()]
        
        constraints += [self.bias_constriant(h) for h in self.Hdict.values()]
        constraints += [self.coupling_constriant(j) for j in self.Jdict.values()]

        self.constraints += constraints

    
    def parse_model(self, model):
        return {
            'g_min': model[self.gap].as_fraction(),
            'offset': model[self.offset].as_fraction(),
            'biases': {p: model[h].as_fraction() for p, h in self.Hdict.items()},
            'couplings': {edg: model[j].as_fraction() for edg, j in self.Jdict.items()}
        }

    def save_cfa(self, ising_model, filename=f'CFA0'):
        if not os.path.exists(self.library_path):
            os.mkdir(self.library_path)

        fp = os.path.join(self.library_path, filename)
        with open(fp, 'wb') as file:
            pickle.dump((self.cfa_plmt, ising_model), file)

        with open(f'{fp}.json', 'w') as file:
            json.dump({
                'plmt': self.cfa_plmt,
                'bqm': serialize_bqm(ising_model)
            }, file, indent=4)


    def run_solver(self): 
        slv = Optimize()      
        self.impose_constraints_for_CFA_encoding()
        slv.add(self.constraints)

        slv.maximize(self.gap)
        # print(self.slv.sexpr())

        print(f'start at {datetime.now()}')
        slv.check()
        print(f'end at {datetime.now()}')

        model = slv.model()

        ising_model = None
        if model:
            ising_model = self.parse_model(model)
        return ising_model


    def verify(self, ising_model, close=utils.close):
        bqm = BinaryQuadraticModel.from_ising(*[ising_model[key] for key in ['biases', 'couplings', 'offset']])
        nz = len(ising_model['biases'].items())

        cfa_edict, noncfa_edict, all_elist = defaultdict(list), defaultdict(list), []

        for zs in product([-1, 1], repeat=nz):
            zsA = dict(zip(ising_model['biases'].keys(), zs))
            e = bqm.energy(zsA)

            xas = [utils.spin_to_binary(z) for z in zs]
            xs, _as = tuple(xas[:self.n_x]), tuple(xas[self.n_x:])
            # if xs in self.CFA_truth_table:
            xdict = dict(zip(self.cfa_plmt.keys(), xs))
            if self.f(xdict):
                cfa_edict[xs].append(e)
            else: 
                noncfa_edict[xs].append(e)

            all_elist.append(e)

        def foreach_cfa_e(_asE, ge=0):
            return all([e > ge or close(e, ge) for e in _asE]) and any([close(e, ge) for e in _asE])

        def foreach_noncfa_e(_asE, gap=ising_model['g_min']):
            return all([e > gap or close(e, gap) for e in _asE])

        logic = all([foreach_cfa_e(_asE) for _asE in cfa_edict.values()]) and all([foreach_noncfa_e(_asE) for _asE in noncfa_edict.values()])

        ge, fe = [min([min(_asE) for _asE in o.values()]) for o in [cfa_edict, noncfa_edict]]
        min_gap = fe - ge

        # assert close(min_gap, ising_model['g_min'])

        return logic and close(min_gap, ising_model['g_min'])


