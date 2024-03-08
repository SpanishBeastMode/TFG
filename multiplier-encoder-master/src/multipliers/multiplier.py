from collections import defaultdict
from fractions import Fraction
from functools import reduce
import itertools
import os
import pickle
from dimod import BinaryQuadraticModel

import dwave_networkx as dnx

from .cell_utils import multiplier_origin
from .draw_utils import graph_2_dot
from .real_graphs import get_real_graph

import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from utils import close, lp, binary_to_spin, spin_to_binary

advantage_sys_list = {
    'na-west-1':  [4.1, 6.2]
}

class Multiplier():
    __allowed = ('data_path')
    def __init__(self, *args, origin_cellwrpr=multiplier_origin(), m=16, pegasus_ver=advantage_sys_list['na-west-1'][0], **kwargs):
        print('Multiplier')

        self.length,self.width = args
        self.origin_cellwrpr = origin_cellwrpr

        self.m = m
        self.pegasus_ver = pegasus_ver

        for k, v in kwargs.items():
            if k in self.__allowed:
                setattr(self, k, v)

        self.pegasus_to_linear = dnx.pegasus_coordinates(self.m).pegasus_to_linear
        
        self.embed()


    def embed(self, adhoc_cfa_plmts=None):
        self.cfa_plmts = self.plmt_cfas()

        self.cfa_graphs = self.construct_cfa_graphs(adhoc_cfa_plmts=adhoc_cfa_plmts)

        self.chains = self.construct_chains()
        
    @property
    def embedding_graph(self):
        nodes = defaultdict(list)
        for adder in self.cfa_graphs:
            for cfa_graph in adder:
                for v, qbt in cfa_graph.items():
                    nodes[v].append(qbt)
        edges = {'cfa': [list(dnx.pegasus_graph(self.m, node_list=cfa_graph.values()).edges) for adder_graph in self.cfa_graphs for cfa_graph in adder_graph]}
        edges = {**edges, **self.chains}
        edges = {key: list(itertools.chain(*chains)) for key, chains in edges.items()}

        graph = {'nodes': nodes, 'edges': edges}
        return graph
            
    def draw(self, graph, bqm={}, qbt_values={}):
        data_path = '.' if not hasattr(self, 'data_path') else self.data_path
        if not os.path.exists(data_path):
            os.makedirs(data_path)
            
        flag = True if bqm else False
        dot_file, svg_file = [os.path.join(data_path, f'{self.length}bit*{self.width}bit_bqm={flag}{tp}') for tp in ['.dot', '.svg']]

        real_graph = get_real_graph(self.pegasus_ver)
        with open(dot_file, "w") as file:
            graph_2_dot(file, graph, real_graph, bqm, qbt_values, m=self.m)
        os.system("neato -Tsvg {} -o {}".format(dot_file, svg_file))
        os.system("rm {}".format(dot_file))

    def encode(self, adhoc_cfa_plmts=None):
        cfa_plmts = adhoc_cfa_plmts if adhoc_cfa_plmts else self.cfa_plmts
        
        def encode_cfa(cfa_graph, cfa):
            cfa_plmt, cfa_ising_model = cfa
            p2q = {p: cfa_graph[k] for k, p in cfa_plmt.items()}
            mapped = lambda ps: p2q[ps] if type(ps) == int else tuple(sorted([p2q[p] for p in ps]))
            cfa = {k: v if type(v) != dict else {mapped(ps): c for ps, c in v.items()}
                        for k, v in cfa_ising_model.items()}
            return cfa

        encoded_cfas = [encode_cfa(self.cfa_graphs[i][j], cfa_plmts[i][j]) 
                                for i in range(self.width) 
                                    for j in range(self.length)]
    

        def aggriate_coefficients(encoded_cfas):
            model = {'offset': 0, 'biases': defaultdict(Fraction), 'couplings': defaultdict(Fraction)}
            for ecfa in encoded_cfas:
                for k, v in ecfa.items():
                    if type(v) != dict:
                        if k =='offset': model[k] += v
                        continue
                    for qs, c in v.items():
                        qs = qs if type(qs) == int else tuple(sorted(qs))
                        model[k][qs] += c
            return model

        model = aggriate_coefficients(encoded_cfas)

        for k, chains in self.chains.items():
            for edg in itertools.chain(*chains):
                model['offset'] += 2
                model['couplings'][tuple(sorted(edg))] = -2

        self.model = model
        self.cfa_plmts = cfa_plmts
        return model

    def is_compatiable_with_real_graph(self, model):
        real_graph = get_real_graph(self.pegasus_ver)
        all_qbts_valid = all([qbt in real_graph['nodes'] for qbt in model['biases'].keys()])
        all_edges_valid = all([edge in real_graph['edges'] for edge in model['couplings'].keys()])
        return all_qbts_valid and all_edges_valid

    @property
    def shared_graph(self):
        shared_graph = {'nodes': defaultdict(list), 'edges': defaultdict(list)}
        for _out, _in in self.shared_from.items():
            if _in == 'in2':
                shared_qbts = [self.cfa_graphs[i][j][_in] for i in range(1, self.width) for j in range(0, self.length - 1)]
            elif _in in ['c_in', 'enable']:
                shared_qbts = [self.cfa_graphs[i][j][_in] for i in range(0, self.width) for j in range(1, self.length)]
            elif _in == 'in1':
                shared_qbts = [self.cfa_graphs[i][j][_in] for i in range(1, self.width) for j in range(self.length)]
            shared_graph['nodes'][_in] += shared_qbts

        _in_e = ('c_in', 'enable')
        if reduce(lambda x, y: x and y, [v in self.shared_from.values() for v in _in_e]):
            shared_couplings = [tuple(sorted([self.cfa_graphs[i][j][v] for v in _in_e]))
                                        for i in range(0, self.width) for j in range(1, self.length)]

            shared_graph['edges']['cfa'] += shared_couplings
        return shared_graph

    @property
    def interface_qbts(self):
        l, w = self.length, self.width
        initial_zero_value_bits = {
            'in2': [self.cfa_graphs[0][j]['in2'] for j in range(l)],
            'c_in': [self.cfa_graphs[i][0]['c_in'] for i in range(w)]
        }
        input_qbts = {
            'in1': [self.cfa_graphs[0][j]['in1'] for j in range(l)],
            'enable': [self.cfa_graphs[i][0]['enable'] for i in range(w)],
            'out': [self.cfa_graphs[i][0]['out'] for i in range(w)] + \
                        [self.cfa_graphs[w-1][j]['out'] for j in range(1, l)] + \
                            [self.cfa_graphs[w-1][l-1]['c_out']]
        }
        return initial_zero_value_bits, input_qbts
    
    def interface_qbts_for_reduced_graph(self):
        l, w = self.length, self.width
        input_qbts = {
            'in1': [self.cfa_graphs[0][j]['in1'] for j in range(l)],
            'enable': [self.cfa_graphs[i][0]['enable'] for i in range(w)]
        }
        return input_qbts


    def config_inputs(self, inputs):        
        l, w = self.length, self.width

        to_binary = lambda input: [int(bit) for bit in bin(input)[2:][::-1]]

        vs = {}
        for v, input in inputs.items():
            xs = to_binary(input)

            n = l if v == 'A' else w if v == 'B' else l+w if v == 'out' else None
            xs = xs if len(xs) == n else xs + [0 for _ in range(n - len(xs))]

            vs[v]  = xs
            
        self.inputs = vs
        return vs

    def impose_inputs_on_qbts(self, new_input_qbts=None):
        initial_zero_value_bits, input_qbts = self.interface_qbts

        bits = ['in1', 'enable', 'out']
        variables = ['A', 'B', 'out']
        v2b = dict(zip(variables, bits))

        qbt_values = {}
        for v, qbts in initial_zero_value_bits.items():
            zs = [binary_to_spin(0) for _ in range(len(qbts))]
            qbts = new_input_qbts[v] if new_input_qbts else qbts
            qv = dict(zip(qbts, zs))
            qbt_values = {**qbt_values, **qv}

        for v, xs in self.inputs.items():
            zs = [binary_to_spin(x) for x in xs]
            qbts = new_input_qbts[v] if new_input_qbts else input_qbts[v2b[v]]
            qv = dict(zip(qbts, zs))
            qbt_values = {**qbt_values, **qv}

        self.qbt_values = qbt_values
        self.input_variables = variables
        return qbt_values

    def _parse_a_sample_in_terms_of_cfa_behaviors(self, sample):
        cfa_energies = [[] for _ in range(self.width)]
        cfa_excitations = [[] for _ in range(self.width)]
        for i in range(self.width):
            for j in range(self.length):
                (cfa_plmt, cfa_ising_model), cfa_graph = self.cfa_plmts[i][j], self.cfa_graphs[i][j]
                cfa_value = {cfa_plmt[v]: sample[qbt] for v, qbt in cfa_graph.items()}

                cfa_bqm = BinaryQuadraticModel.from_ising(*[cfa_ising_model[key] for key in ['biases', 'couplings', 'offset']])
                energy = cfa_bqm.energy(cfa_value)
                cfa_energies[i].append(energy)

                gap = cfa_ising_model['g_min']
                excitation = (energy >= gap) or close(energy, gap)
                cfa_excitations[i].append(excitation)
        return cfa_energies, cfa_excitations

    def parse_a_sample_of_all_qbt_values(self, sample, inputs_imposed_via='api'):
        cfa_energies, cfa_excitations = self._parse_a_sample_in_terms_of_cfa_behaviors(sample)

        # truth_of_cfas = all([close(cfa_energy, 0) for adder_energy in cfa_energies for cfa_energy in adder_energy])
        truth_of_cfas = all([not cfa_excitation for adder_excitation in cfa_excitations for cfa_excitation in adder_excitation])

        chains = {v: [[(spin_to_binary(sample[fqbt]), spin_to_binary(sample[tqbt])) for fqbt, tqbt in ch] for ch in chs] 
                        for v, chs in self.chains.items()}
        truth_of_chains = {
            k: all([all([x0 == x1 for x0, x1 in ch]) for ch in chs]) for k, chs in chains.items()
        }
        truth_of_chains = all(truth_of_chains.values())

        if inputs_imposed_via in ['flux_bias', 'extra_chains']:
            initial_zero_value_bits, input_qbts = self.interface_qbts
            zero_value_qbt_values = {v: [spin_to_binary(sample[qbt]) for qbt in qbts] for v, qbts in initial_zero_value_bits.items()}
            truth_of_zero_values = all((x==0 for x in itertools.chain(*zero_value_qbt_values.values())))
        elif inputs_imposed_via in ['api', 'adhoc_encoding']:
            truth_of_zero_values = True
            input_qbts = self.interface_qbts_for_reduced_graph()

        input_qbt_values = {v: [spin_to_binary(sample[qbt]) for qbt in qbts] for v, qbts in input_qbts.items()}
        input_numbers = {v: sum([x * 2**i for i, x in enumerate(xs)]) for v, xs in input_qbt_values.items()}

        if 'out' not in input_numbers.keys():
            input_numbers['out'] = sum([x * 2**i for i, x in enumerate(self.inputs['out'])])

        return [input_numbers[v] for v in ['in1', 'enable', 'out']], (truth_of_cfas, truth_of_chains, truth_of_zero_values), cfa_excitations

