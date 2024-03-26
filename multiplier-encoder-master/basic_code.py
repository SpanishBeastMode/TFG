from src.multipliers import *
from src.solve import *
from src.generate_inputs import *
from src.initialize import *
from src.generate_inputs import *
from dwave.system.samplers import DWaveSampler
import numpy as np


def solve(inputs, *args, **kwargs):
    initializer = Initializer(inputs, *args, **kwargs)

    solver = Solver(inputs, initializer)

    Ta = 10

    # solutions = solver.forward_annealing(Ta=Ta)

    # ers, solutions = solver.annealing_enhanced_by_pause(params=solver.params, Ta=Ta, anneal_direction='forward') # forward annealing with pause, 10x8 bits, 256271

    # results, solutions = solver.run_pRV_for_diff_Tps(pause_enhanced=True, Tp_list=[100], Sp_list=[(38-i)/100 for i in range(5)], Ta=10) # 11x8 bits, 511789

    results, solutions = solver.gready_pRV(pause_enhanced=False, Tp_list=[100], Sp_list=[(38-i)/100 for i in range(4)], Ta=10)

    # results, solutions = solver.gready_pRV(pause_enhanced=False, filtered=False)

    # solution = solver.run_BFS_based_iterative_reverse_annealing(pause_enhanced=False, Ta=10)



    # For run fw for different tas
    for key, solution in solutions.items():
        for solv, solu in solution.items():
            for sol in solu[0]:
                if sol[0]*sol[1] == sol[2]:
                    print("Solution found: ", sol)

    # For only forward annealing
    # i=0
    # for solution in solutions[1]:
        # print(solution)
        # for sol in solution:
        #     if i < 25:
        #         print(sol)
        #     i += 1
        #     if sol[0]*sol[1] == sol[2]:
        #         print("Solution found: ", sol)

    return


def run_in_batches(*args, counter=1, **kwargs):
    for out in generate_inputs(*args, counter=counter):
        print("first number: ", out)
        # if length == 8 and width == 8:
        #     assert kwargs['using_improved_library']  # which found the ground state for the following inputs
        #     out = 59989
        # if length == 4 and width == 4:
        #     assert kwargs['using_improved_library']  # which found the ground state for the following inputs
        #     out = 143

        inputs = {'out': out}

        solve(inputs, *args, **kwargs)


if __name__ == '__main__':
    vias = ['api', 'flux_biases', 'extra_chains', 'adhoc_encoding']
    size = 8
    length, width = 8, size
    run_in_batches(
        length, width, via=vias[1], using_improved_library=True, sys=4.1)
