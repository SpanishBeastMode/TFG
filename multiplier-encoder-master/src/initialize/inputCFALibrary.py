from collections import defaultdict
import copy
from itertools import chain, combinations, product

from .inputCFAEncoder import InputCFAEncoder

import os
import pickle
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from multipliers.multiplier_ver2 import Multiplier_ver2
from cfa.CFA import CFA_ver2
from utils import lp_via_z3, lp

class InputCFALibrary():
    def __init__(self, version='_ver2', ancillar_enhanced=True, cfa_index=0):
        assert  version  == '_ver2' and cfa_index == 0

        self.multiplier = eval(f'Multiplier{version}')(3, 4, cfa_index=cfa_index, library_path=lp)
        self.generic_cfa = self.multiplier.cfa
        self.cfa_plmt = self.generic_cfa[0]
        
        cfa_plmt, cfa_ising_model = self.generic_cfa
        cfa_graph = {
            'nodes': cfa_ising_model['biases'].values(),
            'edges': cfa_ising_model['couplings'].keys()
        }
        arg = (eval(f'CFA{version}'), 7, 4)
        self.args = (*arg, cfa_plmt, cfa_graph)

        inputs = {'out': 24}
        self.multiplier.config_inputs(inputs)
        self.multiplier.impose_inputs_on_qbts()

        self.library_path = os.path.join(lp_via_z3, f'{self.__class__.__name__}_for_CFA{version}_{cfa_index}_with_ancillar_enhanced={ancillar_enhanced}')
        if not os.path.exists(self.library_path):
            os.mkdir(self.library_path)

        assert self.multiplier.length >= 3 and self.multiplier.width >= 4
        self.inputCFA_sharing_dependency, self.encode_priority = self.determine_encode_order_and_configure_sharing_constraints()
        self.ancillar_enhanced = ancillar_enhanced

        self.encode_on_the_fly()


    def _get_cfa_sharing_dependency(self, coord):
        assert self.multiplier.__class__.__name__ == 'Multiplier_ver2'

        restricted_by = lambda i, j: {
            'out': ((i+1, j-1), 'in2'), 
            'in2': ((i-1, j+1), 'out'),
            
            'c_out': ((i, j+1), 'c_in'),
            'c_in': ((i, j-1), 'c_out'),
            'enable_out': ((i, j+1), 'enable'),
            'enable': ((i, j-1), 'enable_out'),
            ('c_out', 'enable_out'): ((i, j+1), ('c_in', 'enable')),
            ('c_in', 'enable'): ((i, j-1), ('c_out', 'enable_out')),
        }

        result = restricted_by(*coord)

        valid_coord = lambda i, j: i in range(self.multiplier.width) and j in range(self.multiplier.length)
        result = {f: (coord, t) for f, (coord, t) in result.items() if valid_coord(*coord)}
        return result

    @property
    def _interface_cfas(self):
        return [(i, j) for i in range(self.multiplier.width) for j in range(self.multiplier.length) if i == 0 or j == 0 or i == self.multiplier.width-1]

    @property
    def interface_cfa_sharing_dependency(self):
        return {coord: self._get_cfa_sharing_dependency(coord) for coord in self._interface_cfas}
    
    def _calcu_cfa_freedom(self, i, j):
        return sum([qbt in self.multiplier.qbt_values.keys() for qbt in self.multiplier.cfa_graphs[i][j].values()])

    @property
    def interface_cfa_freedom(self):
        return {(i, j): self._calcu_cfa_freedom(i, j) for i in range(self.multiplier.width) for j in range(self.multiplier.length) }

    def determine_encode_order_and_configure_sharing_constraints(self):
        input_cfa_sharing_dependency = copy.deepcopy(self.interface_cfa_sharing_dependency)
        cfa_freedom = copy.deepcopy(self.interface_cfa_freedom)

        reduced_inputCFA_sharing_dependency = {fcoord: {f: (tcoord, t) for f, (tcoord, t) in dependency.items() if cfa_freedom[fcoord] >= cfa_freedom[tcoord]} 
                                                            for fcoord, dependency in input_cfa_sharing_dependency.items()}

        freedom_group = defaultdict(list)
        for coord, fd in cfa_freedom.items():
            if fd == 0: continue
            freedom_group[fd].append(coord)

        def mutualy_dependent(agroup):
            return [(f, t) for f, t in combinations(agroup, 2) if t in [coord for coord, _ in reduced_inputCFA_sharing_dependency[f].values()]]
            
        conflicts = {fd: mutualy_dependent(agroup) for fd, agroup in freedom_group.items()}

        def _cfa_impacts(coord):
            return sum([coord in [c for c, _ in dependency.values()] for dependency in reduced_inputCFA_sharing_dependency.values()])

        def _break_mutual_dependence(tcoord, fdependency):
            updated_dependency = copy.deepcopy(fdependency)
            for f, (tc, t )in fdependency.items():
                if tc != tcoord: continue
                del updated_dependency[f]
            return updated_dependency


        encode_priority = copy.deepcopy(freedom_group)
        for fd, agroup in freedom_group.items():
            conflicts = mutualy_dependent(agroup)
            if not conflicts: continue
                
            impacts = {c: _cfa_impacts(c) for c in chain(*conflicts)} 
            impacts = dict(sorted(impacts.items(), key=lambda item: -item[1]))

            for conflict in conflicts:
                fcoord, tcoord = sorted(conflict, key=lambda coord: -impacts[coord])
                fdependency = reduced_inputCFA_sharing_dependency[fcoord]
                updated_dependency = _break_mutual_dependence(tcoord, fdependency)
                reduced_inputCFA_sharing_dependency[fcoord] = updated_dependency


            encode_priority[fd] = [c for c in agroup if c not in chain(*conflicts)] + list(impacts.keys())
        
        encode_priority = dict(sorted(encode_priority.items(), key=lambda item: item[0]))
        return reduced_inputCFA_sharing_dependency, encode_priority


    def _cfa_input(self, i, j):
        default_inputs = {}
        if i == 0:
            default_inputs['in2'] = 0
        if j == 0:
            default_inputs['c_in'] = 0
        
        out_inputs = []
        if j == 0 or i == self.multiplier.width-1:
            out_inputs = [{'out': x} for x in [0, 1]]
        if j == self.multiplier.length-1 and i == self.multiplier.width-1:
            c_out_inputs = [{'c_out': x} for x in [0, 1]]
            out_inputs = [{**c_out, **o} for c_out, o in product(c_out_inputs, out_inputs)]
        return [{**out_input, **default_inputs} for out_input in out_inputs] if out_inputs else [default_inputs]

    @property
    def interface_cfa_inputs(self):
        return {coord: self._cfa_input(*coord) for coord in self._interface_cfas}


    def search_by_an_InputCFAEncoder(self, strict, inputs, relied_dependency, fp):
        try:
            with open(fp, 'rb') as file:
                cfa = pickle.load(file)
        except FileNotFoundError:
            encoder = InputCFAEncoder(inputs, *self.args, ancillar_enhanced=self.ancillar_enhanced, strict=strict)
            encoder.impose_extra_dependency_contraints(relied_dependency)
            ising_model = encoder.run_solver()
            # satisfied = encoder.verify(ising_model)
            # assert satisfied
            encoder.save_cfa(ising_model, fp)
            cfa = (encoder.cfa_plmt, ising_model)
        return cfa


    def encode_for_a_combination_of_inputs(self, coord, coord_path, strict, known_sharing_dependency):
        results = []
        for inputs in self._cfa_input(*coord):
            print(f'\tfor inputs: {inputs}')

            fp = os.path.join(coord_path, f'for_inputs={inputs}')
            cfa = self.search_by_an_InputCFAEncoder(strict, inputs, known_sharing_dependency, fp)
            results.append(cfa)
        return results

    def _get_known_coefficient(self, vs, cfa=None):
        cfa_plmt, cfa_model = cfa
        if type(vs) == str:
            c = cfa_model['biases'][cfa_plmt[vs]]
        else:
            couplings = cfa_model['couplings']
            ps = tuple([cfa_plmt[v] for v in vs])
            ps = ps if ps in couplings.keys() else tuple(ps[::-1])
            c = couplings[ps]
        return c


    def _get_bidirectionally_horizontal_selfdependency(self, coord, inputs_index=0, Library=None):
        assert self.multiplier.__class__.__name__ == 'Multiplier_ver2'

        restricted_by = lambda i, j: {
            'c_out': ((i, j), 'c_in'),
            'c_in': ((i, j), 'c_out'),
            'enable_out': ((i, j), 'enable'),
            'enable': ((i, j), 'enable_out'),
            ('c_out', 'enable_out'): ((i, j), ('c_in', 'enable')),
            ('c_in', 'enable'): ((i, j), ('c_out', 'enable_out')),
        }

        result = restricted_by(*coord)

        result = {f: self._get_known_coefficient(t, Library[coord][inputs_index]) for f, (_, t) in result.items()}
        return result
        

    @property
    def out_dependency(self):
        return {fcoord: {tcoord for tcoord, _ in dependency.values() if self.interface_cfa_freedom[tcoord] > 0 and 'out' in self._cfa_input(*tcoord)[0].keys()} 
                        for fcoord, dependency in self.inputCFA_sharing_dependency.items()}


    def encode_on_the_fly(self, reused_coords=[(0, 2), (2, 1)]):
        encode_in_sequence = list(chain(*self.encode_priority.values()))

        Library = defaultdict(list)
        strictness = {}
        for coord in encode_in_sequence:
            print(f'{coord}:')
            coord_path = os.path.join(self.library_path, f'inputCFA_at_{coord}', )
            if not os.path.exists(coord_path):
                os.mkdir(coord_path)

            strict = True if coord in reused_coords else False

            sharing_dependency = self.inputCFA_sharing_dependency[coord]
            print(f'\t sharing_dependency: {sharing_dependency}')

            relied_coords = {tcoord for tcoord, _ in sharing_dependency.values() if self._calcu_cfa_freedom(*tcoord) > 0 and len(self._cfa_input(*tcoord)) > 1}
            print(f'{coord} reply on #{relied_coords} inputCFAs: {sharing_dependency}')
            
            if relied_coords:
                assert len(relied_coords) == 1
                relied_coord, = relied_coords
                for i, resulting_cfa in enumerate(Library[relied_coord]):
                    coord_subpath = os.path.join(coord_path, f'rely_on_inputCFA={relied_coord}-index={i}')
                    if not os.path.exists(coord_subpath):
                        os.mkdir(coord_subpath)

                    known_sharing_dependency = {f: self._get_known_coefficient(t, self.generic_cfa if self.interface_cfa_freedom[tcoord] == 0 else 
                                                                                    resulting_cfa if tcoord == relied_coord else
                                                                                    None) 
                                                    for f, (tcoord, t) in sharing_dependency.items()}
                    results = self.encode_for_a_combination_of_inputs(coord, coord_subpath, strict, known_sharing_dependency)
                    Library[coord].append(results)
            else:
                known_sharing_dependency = {f: self._get_known_coefficient(t, self.generic_cfa if self.interface_cfa_freedom[tcoord] == 0 else
                                                                                Library[tcoord][0]) for f, (tcoord, t) in sharing_dependency.items()}
            
                if coord == (3, 1):
                    inputs_comb = self._cfa_input(*coord)
                    assert len(inputs_comb) == 2
                    for i, inputs in enumerate(inputs_comb):
                        if i == 1:
                            inputs_dependency = self._get_bidirectionally_horizontal_selfdependency(coord, 0, Library=Library)
                            known_sharing_dependency = {**known_sharing_dependency, **inputs_dependency}

                        fp = os.path.join(coord_path, f'for_inputs={inputs}')
                        cfa = self.search_by_an_InputCFAEncoder(strict, inputs, known_sharing_dependency, fp)
                        Library[coord].append(cfa)

                else:
                    Library[coord] = self.encode_for_a_combination_of_inputs(coord, coord_path, strict, known_sharing_dependency)
            strictness[coord] = strict

        self.Library = Library
        self.strictness = strictness
        return Library


    def _check_sharing_compatibility_for_combinations_of_inputCFAs(self, coord, cfas):
        shared_qbts = self.multiplier.shared_from
        shared_couplings = {('c_out', 'enable_out'): ('c_in', 'enable')}

        input_variables = self.interface_cfa_inputs[coord][0].keys()
        reduced_shared_qbts = {f: t for f, t in shared_qbts.items() if f not in input_variables and t not in input_variables}
        reduced_shared_couplings = {f: t for f, t in shared_couplings.items() if not {*f, *t}.intersection(set(input_variables))}

        for vs in reduced_shared_qbts.items():
            s = sum([self._get_known_coefficient(v, cfa) for v, cfa in zip(vs, cfas)])
            assert s >= -4 and s <= 4

        for edgs in reduced_shared_couplings.items():
            s = sum([self._get_known_coefficient(edg, cfa) for edg, cfa in zip(edgs, cfas)])
            assert s >= -2 and s <= 1

        return True

    def verify_horizontal_compatibility_for_reused_interface_cfas(self, reused_coords=[(2, 1), (0, 2)], 
                                                                        depedent_coords={
                                                                            (2, 2): (2, 1), 
                                                                            (2, 0): (2, 1)}):

        for coord in reused_coords:
            for cfas in product(self.Library[coord], repeat=2):
                self._check_sharing_compatibility_for_combinations_of_inputCFAs(coord, cfas)

        for tcoord, fcoord in depedent_coords.items():
            for tcfas, fcfa in zip(self.Library[tcoord], self.Library[fcoord]):
                for tcfa in tcfas:
                    cfas = [fcfa, tcfa] 
                    if tcoord == (2, 0):
                        cfas = cfas[::-1]

                    self._check_sharing_compatibility_for_combinations_of_inputCFAs(fcoord, cfas)
        return True
    

    def choose_library_for_multiplier_inputs_coord(self, coord, inputs, l, w):
        i, j = coord
        Library = copy.deepcopy(self.Library)

        def _maped_library_coord(coord, l, w):
            i, j = coord
            assert i == 0 or j == 0 or i == w-1

            if i == w-1 or j == 0:
                maxL, maxW = max([j for _, j in Library.keys()]), max([i for i, _ in Library.keys()])
                li = maxW if i == w-1 else maxW - 1 if i > 1 else i
                lj = maxL if j == l-1 else j if j == 0 else 1
                lcoord = (li, lj)
            elif i == 0:
                lcoord = (0, 1 if j == 1 else 2)
            else:
                lcoord = None

            print(f'map {coord} to library {lcoord}')
            return lcoord
        
        lcoord = _maped_library_coord(coord, l, w)

        relied_library_index = inputs['out'][i+j+1] if i == w-1 and j == 0 else \
                                inputs['out'][i+j-1] if i == w-1 and j == l-1 else None
        
        clibrary = Library[lcoord] if relied_library_index is None else Library[lcoord][relied_library_index]
        return clibrary
            