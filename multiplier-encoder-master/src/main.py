from multipliers import *
from initialize import *
from solve import *

import json
import os

# from utils import dp

if __name__ == '__main__':
    # ====================Embedding=======================
    multiplier_versions = [
        Multiplier_ver1_to_link_B,
        Multiplier_ver2_to_link_B,
        Multiplier_ver3_to_link_B,

        Multiplier_ver3,
        Multiplier_ver2
    ]

    length, width, inputs = 8, 8, {'out': 59989}

    for M in multiplier_versions[-1:]:
        multiplier = M(length, width)
        multiplier.draw(multiplier.embedding_graph)

    # model = multiplier.encode()
    # multiplier.draw(multiplier.embedding_graph, model)
    # multiplier.draw(multiplier.shared_graph, model)

    # initial_zero_value_bits, input_qbts = multiplier.interface_qbts
    # input_graph = {'nodes': {**initial_zero_value_bits, **input_qbts}, 'edges': {}}
    # multiplier.draw(input_graph)
    # ====================Embedding=======================

    # assert M == Multiplier_ver2
    # multiplier.config_inputs(inputs)
    # qbt_values = multiplier.impose_inputs_on_qbts()
    # enhanced_qbt_values = Preprocess(qbt_values, multiplier).preprocess()
    # extra_qbt_values = {qbt: z for qbt, z in enhanced_qbt_values.items() if qbt not in qbt_values.keys()}
    # multiplier.draw(multiplier.embedding_graph, qbt_values=extra_qbt_values)

    # ====================Solving=======================
    # Initializer supports the Multiplier_ver2 by default
    if True:
        length, width, inputs = 8, 8, {'out': 59989}
        initializer = Initializer(inputs, length, width, via='flux_biases', using_improved_library=True, sys=4.1)

        solver = Solver(inputs, initializer)

        solver.run_FW_for_diff_Tas(Ta_list=[10])

        solver.run_pFW_for_diff_Tps(Tp_list=[1, 10, 30, 50, 100], Sp_list=[(30+i)/100 for i in range(0, 10, 3)], Ta=10)

    if True:
        length, width, inputs = 13, 8, {'out': 2055941}
        initializer = Initializer(inputs, length, width, via='flux_biases', using_improved_library=True, sys=4.1)
        solver = Solver(inputs, initializer)

        solver.run_pRV_for_diff_Tps(pause_enhanced=False, Tp_list=[1, 10, 30, 50, 100], Sp_list=[(33+i)/100 for i in range(0, 10, 3)], Ta=10)

    if True:
        length, width, inputs = 14, 8, {'out': 4111631}
        initializer = Initializer(inputs, length, width, via='flux_biases', using_improved_library=True, sys=4.1)
        solver = Solver(inputs, initializer)

        solver.run_BFS_based_iterative_reverse_annealing(pause_enhanced=False, Ta=20)

    if True:
        length, width, inputs = 15, 8, {'out': 8219999}
        initializer = Initializer(inputs, length, width, via='flux_biases', using_improved_library=True, sys=4.1)
        solver = Solver(inputs, initializer)

        solver.gready_pRV(pause_enhanced=False, filtered=True)
    # ====================Solving=======================