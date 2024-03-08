import math
import os

sp = os.path.dirname(os.path.dirname(__file__))
lp =  os.path.join(sp, 'library')
lp_via_z3 = os.path.join(sp, 'InputCFALibrary_generated_by_z3_solver')
improved_lp = os.path.join(sp, 'improved_library')

epsilon = 1e-6

def get_data_path(sys):
    sysp = os.path.join(sp, f'data_obtained_from_sys{sys}')
    dp = os.path.join(sysp, 'data_generated_by_library')
    improved_dp = os.path.join(sysp, 'data_generated_by_improved_library')
    return dp, improved_dp

def close(energy, e):
    return math.isclose(energy, e, abs_tol=epsilon)

def spin_to_binary(z):
    assert z in [-1, 1]
    return 1 if z == 1 else 0

def binary_to_spin(x):
    assert type(x) == bool or x in [0, 1]
    return 1 if x else -1

def hamming_distance(sample0, sample1):
    xs0, xs1 = [{q: spin_to_binary(i) for q, i in sample.items()} for sample in [sample0, sample1]]
    distance = sum([x ^ xs1[key] for key, x in xs0.items()])
    return distance
    