from src.multipliers import *
from src.solve import *
from src.generate_inputs import *
from src.initialize import *
from src.generate_inputs import *
from dwave.system.samplers import DWaveSampler


def solve(inputs, *args, **kwargs):
    initializer = Initializer(inputs, *args, **kwargs)

    solver = Solver(inputs, initializer)

    # solution = solver.annealing_enhanced_by_pause(
        # params=solver.params, anneal_direction='forward', Tp=100, Sp_list=[38/100], Ta=10),
    
    solution = solver.run_BFS_based_iterative_reverse_annealing(pause_enhanced=False, Ta=10, Tp_list=[1], Sp_list=[0.4])

    print(solution)

    return


def run_in_batches(*args, counter=1, **kwargs):
    for out in generate_inputs(*args, counter=counter):
        print("first number: ", out)
        # if length == 8 and width == 8:
        #     assert kwargs['using_improved_library']  # which found the ground state for the following inputs
        #     out = 59989

        inputs = {'out': out}

        solve(inputs, *args, **kwargs)


if __name__ == '__main__':
    vias = ['api', 'flux_biases', 'extra_chains', 'adhoc_encoding']
    size = 8
    length, width = 11, size
    run_in_batches(
        length, width, via=vias[1], using_improved_library=True, sys=4.1)
