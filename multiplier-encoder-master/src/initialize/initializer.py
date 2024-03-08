from itertools import chain
from dimod import BinaryQuadraticModel

from .inputCFALibrary import InputCFALibrary

import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from multipliers.mixins import EmbedInitializerChain
from cfa.CFA import CFA_ver2
from multipliers.cell_utils import multiplier_origin, multiplier_origin_diagonal_45_shifted_2units
from multipliers.multiplier_ver2 import Multiplier_ver2
from utils import lp, improved_lp, get_data_path

class Initializer():
    def __init__(self, inputs, length, width, multiplier_version='_ver2', via='api', ancillar_enhanced=False, using_improved_library=False, sys=4.1):
        assert 'out' in inputs.keys()

        dp, improved_dp = get_data_path(sys)

        cfa_index, library_path, sys_path = (1, improved_lp, improved_dp) if using_improved_library else \
                                                (0, lp, dp)
        
        data_path = os.path.join(sys_path, f'Multiplier{multiplier_version}_{cfa_index}')

        if sys == 4.1:
            origin_cellwrpr = multiplier_origin_diagonal_45_shifted_2units() if via in ['extra_chains', 'adhoc_encoding'] else multiplier_origin()
        elif sys == 6.2:
            assert via != 'extra_chains'
            pointer = multiplier_origin()
            for _ in range(3):
                pointer = pointer.next(*(0, 'down'))
            pointer = pointer.next(*(45, 'down'))
            origin_cellwrpr = pointer

        kwargs = dict(
            cfa_index=cfa_index,
            library_path=library_path,
            data_path=data_path,
            origin_cellwrpr=origin_cellwrpr,
            pegasus_ver=sys
        )

        multiplier = eval(f'Multiplier{multiplier_version}')(length, width, **kwargs)
        multiplier.encode()
        multiplier.config_inputs(inputs)
        multiplier.impose_inputs_on_qbts()
        self.multiplier = multiplier

        if via == 'adhoc_encoding':
            assert not using_improved_library and kwargs['cfa_index'] == 0
            self.ancillar_enhanced = ancillar_enhanced

        self.via = via
        self.run = eval(f'self.via_{via}')

        self.data_path = os.path.join(data_path, f'inputs_imposed_via_{via}',
                                                    f'ancillar_enhanced={ancillar_enhanced}' if via == 'adhoc_encoding' else '',
                                                        f'{length}bits*{width}bits',
                                                            f'factoring_{inputs["out"]}')


    def via_api(self):
        return self.multiplier.model

    def via_flux_biases(self):
        # self.multiplier.draw(self.multiplier.embedding_graph)
        
        # initial_zero_value_bits, input_qbts = self.multiplier.interface_qbts
        # input_graph = {'nodes': {**initial_zero_value_bits, **input_qbts}, 'edges': {}}
        # self.multiplier.draw(input_graph, qbt_values=self.multiplier.qbt_values)

        assert self.multiplier.is_compatiable_with_real_graph(self.multiplier.model)
        return self.multiplier.model

    def via_extra_chains(self):
        assert self.multiplier.origin_cellwrpr == multiplier_origin_diagonal_45_shifted_2units

        extra_chains = EmbedInitializerChain(self.multiplier, self.multiplier.inputs).construct_extra_chains()
        proxy_qbts_dict = {variable: proxy_qbts for variable, (proxy_qbts, _) in extra_chains.items()}
        proxy_edgs_dict =  {variable: proxy_edgs for variable, (_, proxy_edgs) in extra_chains.items()}
        self.multiplier.impose_inputs_on_qbts(new_input_qbts=proxy_qbts_dict)

        input_graph = {'nodes': proxy_qbts_dict, 'edges': {'cfa': list(chain(*proxy_edgs_dict.values()))}}
        self.multiplier.draw(input_graph)

        model = self.multiplier.model
        for edg in chain(*proxy_edgs_dict.values()):
            edg = tuple(sorted(edg))
            model['couplings'][edg] = -2
            model['offset'] += 2

        return model

    def _choose_inputCFA_for_inputs_coord(self, coord, inputs, inputCFA_library_occupation, l=None, w=None):
        library = inputCFA_library_occupation[coord]
        i, j = coord
        if j == 0 or i == w-1:
            out_index = inputs['out'][i+j]
            if i == w-1 and j == l-1:
                out_index = 2*inputs['out'][i+j+1] + out_index
            inputCFA = library[out_index]
        elif i == 0:
            inputCFA = library[0]
        else:
            inputCFA = None
        return inputCFA

    
    def _adhoc_InputCFA_mapping(self, ancillar_enhanced=None):
        assert self.multiplier.version == '_ver2' and self.multiplier.cfa_index == 0

        encoder = InputCFALibrary(version='_ver2', ancillar_enhanced=ancillar_enhanced)

        inputs = self.multiplier.inputs
        l, w = self.multiplier.length, self.multiplier.width

        inputCFA_library_occupation = {(i, j): encoder.choose_library_for_multiplier_inputs_coord((i, j), inputs, l, w)
                                                for i in range(w) for j in range(l) if i == 0 or j == 0 or i == w-1}
        adhoc_cfa_plmts = [[self._choose_inputCFA_for_inputs_coord((i, j), inputs, inputCFA_library_occupation, l=l, w=w) if i == 0 or j == 0 or i == w-1 else 
                            self.multiplier.cfa for j in range(l)] for i in range(w)]
        return adhoc_cfa_plmts

    def via_adhoc_encoding(self):
        adhoc_cfa_plmts =  self._adhoc_InputCFA_mapping(ancillar_enhanced=self.ancillar_enhanced)

        self.multiplier.embed(adhoc_cfa_plmts=adhoc_cfa_plmts)
        model = self.multiplier.encode(adhoc_cfa_plmts=adhoc_cfa_plmts)

        # self.multiplier.draw(self.multiplier.embedding_graph, model)
        return model