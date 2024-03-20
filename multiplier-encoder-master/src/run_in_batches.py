from generate_inputs import generate_inputs
from initialize.initializer import Initializer
from solve.solver import Solver


def solve(inputs, *args, **kwargs):
    initializer = Initializer(inputs, *args, **kwargs)

    solver = Solver(inputs, initializer)

    # solver.run_FW_for_diff_Tas(Ta_list=[10, 20])

    # solver.run_pFW_for_diff_Tps(Tp_list=[100], Sp_list=[(38-i)/100 for i in range(4)], Ta=10)

    # solver.run_pRV_for_diff_Tps(pause_enhanced=True, Tp_list=[100], Sp_list=[(38-i)/100 for i in range(4)], Ta=10)

    # solver.gready_pRV(pause_enhanced=False, Tp_list=[100], Sp_list=[(38-i)/100 for i in range(4)], Ta=10)
    # solver.gready_pRV(pause_enhanced=False, filtered=True)

    solver.run_BFS_based_iterative_reverse_annealing(pause_enhanced=False, Ta=20)


def run_in_batches(*args, counter=1, **kwargs):
    for out in generate_inputs(*args, counter=counter):
        # if length == 8 and width == 8:
        #     assert kwargs['using_improved_library']  # which found the ground state for the following inputs
        #     out = 59989     

        inputs = {'out': out}

        solve(inputs, *args, **kwargs)


if __name__ == '__main__':
    vias = ['api', 'flux_biases', 'extra_chains', 'adhoc_encoding']

    size = 8
    length, width = 16, size
    run_in_batches(length, width, via=vias[1], using_improved_library=True, sys=4.1)
