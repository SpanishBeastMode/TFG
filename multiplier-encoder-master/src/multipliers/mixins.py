from collections import defaultdict
import pickle

import dwave_networkx as dnx

from .cell_utils import multiplier_origin

import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from cfa.CFA import read_CFA
from utils import lp, lp_via_z3, improved_lp


class CFAMixin():

    def _parse_cfa(self, cfa):
        xplmt, ising_model = cfa
        xplmt = dict(zip(*xplmt))
        zP = ising_model['biases'].keys()
        aP = set(zP) - set(xplmt.values())
        aplmt = {f'a{i}': p for i, p in enumerate(aP)}
        cfa_plmt = {**xplmt, **aplmt}
        return cfa_plmt, ising_model

    def choose_cfa(self, cfa_index=0, version='_ver2', library_path=None):
        cfa = read_CFA(cfa_index, f'CFA{version}', library_path=library_path)
        return self._parse_cfa(cfa) if library_path in [lp, improved_lp] else cfa


    def choose_cfa_pair(self, cfa_indexs=[4, 5], version='_ver3_to_link_B', library_path=None):
        cfa_pair = [self._parse_cfa(read_CFA(i, f'CFA{version}', library_path=library_path)) for i in cfa_indexs]

        def _check_connectivity(cfa_plmts, tn=8):
            ready = True
            for v, p0 in cfa_plmts[0].items():
                if p0 >= tn: continue

                p1 = cfa_plmts[1][v]
                if v in self.shared_from.values():
                    ready = p0 == p1

            return ready

        def _adapt_one_cfa(cfa, exchange={4: 5, 5: 4}): # to check its logic
            cfa_plmt, cfa_ising_model = cfa
            
            adapted_cfa_plmt = {v: exchange[p] if p in exchange.keys() else p for v, p in cfa_plmt.items()}
            
            adapted_cfa_ising_model = defaultdict(dict)
            for k, v in cfa_ising_model.items():
                if type(v) != dict: 
                    adapted_cfa_ising_model[k] = v
                    continue
                for ps, c in v.items():
                    if type(ps) == int:
                        ps = exchange[ps] if ps in exchange.keys() else ps
                    else:
                        ps = tuple([exchange[p] if p in exchange.keys() else p for p in ps])
                    adapted_cfa_ising_model[k][ps] = c
            
            return adapted_cfa_plmt, adapted_cfa_ising_model

        cfa_plmts = [cfa_plmt for cfa_plmt, _ in cfa_pair]
        if not _check_connectivity(cfa_plmts):
            cfa_pair = [cfa_pair[0], _adapt_one_cfa(cfa_pair[1])]        

        return cfa_pair



class EmbedMixin():

    def construct_cfa_graphs(self, adhoc_cfa_plmts=None):
        cfa_graphs = []

        ipointer = self.origin_cellwrpr
        for i in range(self.width):
            cfa_graphs.append([])
            jpointer = ipointer
            for j in range(self.length):
                cfa_plmt, _ = self.cfa_plmts[i][j]
                cfa_graph = self.Construct_CFAgraph(jpointer, self.shared_direction, self.pegasus_to_linear).embed(cfa_plmt, self.shared_from)
                if adhoc_cfa_plmts:
                    p2q = {p: cfa_graph[v] for v, p in cfa_plmt.items()}
                    adhoc_cfa_plmt, _ = adhoc_cfa_plmts[i][j]
                    adhoc_cfa_graph = {new_v: p2q[p] for new_v, p in adhoc_cfa_plmt.items()}
                    cfa_graph = adhoc_cfa_graph

                cfa_graphs[i].append(cfa_graph)
                jpointer = jpointer.next(*self.propg_direction['B'])
            ipointer = ipointer.next(*self.propg_direction['A'])

        return cfa_graphs

    def construct_chains_for_adder_carries(self):
        jpointer = self.origin_cellwrpr
        for j in range(self.length - 1):
            jpointer = jpointer.next(*self.propg_direction['B'])

        ipointer = jpointer
        chains = []
        for i in range(self.width - 1):
            cfa_plmt, _ = self.cfa_plmts[i][0]
            ch = self.Construct_Adderchain(ipointer, self.propg_direction, self.pegasus_to_linear).embed(cfa_plmt, self.shared_from)
            chains.append(ch)
            ipointer = ipointer.next(*self.propg_direction['A'])
        return chains


class EmbedInitializerChain():
    def __init__(self, multiplier, inputs=None):
        assert 'out' in inputs.keys()

        self.pegasus_to_linear = dnx.pegasus_coordinates(multiplier.m).pegasus_to_linear
        self.multiplier = multiplier

    def _stretch_out(self, fv, base_cellwrpr, cfa_plmt, mapping={'out': 'in2', 'c_out': 'c_in'}):
        tv = mapping[fv]
        pos = cfa_plmt[fv]
        direction = fv if fv == 'out' else 'B' if fv == 'c_out' else None
        out_cellwrpr = base_cellwrpr.next(*self.multiplier.propg_direction[direction])
        if pos < 8:
            proxy_pos = cfa_plmt[tv]
            out_qbt = base_cellwrpr.cell[pos]
        else:
            pos = cfa_plmt[tv]
            proxy_pos = 7 if pos in range(4) else 3
            out_qbt = out_cellwrpr.cell[pos]
        out_qbt = self.pegasus_to_linear(out_qbt)
        proxy_qbt = self.pegasus_to_linear(out_cellwrpr.cell[proxy_pos])
        return out_qbt, proxy_qbt

    def _construct_extra_chains_for_out(self):
        ipointer = self.multiplier.origin_cellwrpr()
        boundary_cellwrprs = []
        for i in range(self.multiplier.width):
            ipointer = ipointer if i == 0 else ipointer.next(*self.multiplier.propg_direction['A']) 
            cfa_plmt, _ = self.multiplier.cfa_plmts[i][0]
            boundary_cellwrprs.append((ipointer, cfa_plmt))
        jpointer = ipointer
        for j in range(self.multiplier.length):
            if j == 0: continue
            jpointer = jpointer.next(*self.multiplier.propg_direction['B'])
            cfa_plmt, _ = self.multiplier.cfa_plmts[self.multiplier.width-1][j]
            boundary_cellwrprs.append((jpointer, cfa_plmt))

        proxy_qbts, extra_chains = [], []
        for base_cellwrpr, cfa_plmt in boundary_cellwrprs:
            out_qbt, proxy_qbt = self._stretch_out('out', base_cellwrpr, cfa_plmt)
            proxy_qbts.append(proxy_qbt)
            extra_chains.append((out_qbt, proxy_qbt))

        c_out_qbt, c_out_proxy_qbt = self._stretch_out('c_out', *boundary_cellwrprs[-1])
        proxy_qbts.append(c_out_proxy_qbt)
        extra_chains.append((c_out_qbt, c_out_proxy_qbt))

        self.boundary_cellwrprs = boundary_cellwrprs
        return proxy_qbts, extra_chains

    def _construct_extra_chains_for_default_c_in(self):
        boundary_cellwrprs = self.boundary_cellwrprs[:self.multiplier.width]

        proxy_qbts, extra_chains = [], []
        for base_cellwrpr, cfa_plmt in boundary_cellwrprs:
            pos = cfa_plmt['c_in']
            c_in_qbt = self.pegasus_to_linear(base_cellwrpr.cell[pos])

            assert pos in range(4, 8)
            proxy_cellwrpr, proxy_pos = base_cellwrpr.next(*(0, 'down')), pos
            proxy_qbt = self.pegasus_to_linear(proxy_cellwrpr.cell[proxy_pos])
            
            proxy_qbts.append(proxy_qbt)
            extra_chains.append((c_in_qbt, proxy_qbt))
        return proxy_qbts, extra_chains

    def _construct_extra_chains_for_default_in2(self):
        jpointer = self.multiplier.origin_cellwrpr()
        boundary_cellwrprs = []
        for j in range(self.multiplier.length):
            jpointer = jpointer if j == 0 else jpointer.next(*self.multiplier.propg_direction['B']) 
            cfa_plmt, _ = self.multiplier.cfa_plmts[0][j]
            boundary_cellwrprs.append((jpointer, cfa_plmt))

        proxy_qbts, extra_chains = [], []
        for base_cellwrpr, cfa_plmt in boundary_cellwrprs:
            pos = cfa_plmt['in2']
            in2_qbt = self.pegasus_to_linear(base_cellwrpr.cell[pos])

            assert pos in range(0, 4)
            proxy_cellwrpr, proxy_pos = base_cellwrpr.next(*(self.multiplier.propg_direction['out'][0], 'up')), 7
            proxy_qbt = self.pegasus_to_linear(proxy_cellwrpr.cell[proxy_pos])
            
            proxy_qbts.append(proxy_qbt)
            extra_chains.append((in2_qbt, proxy_qbt))
        return proxy_qbts, extra_chains

    def construct_extra_chains(self):
        return {
            'out': self._construct_extra_chains_for_out(),
            'c_in': self._construct_extra_chains_for_default_c_in(),
            'in2': self._construct_extra_chains_for_default_in2()
        }

