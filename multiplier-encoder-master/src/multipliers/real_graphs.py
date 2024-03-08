import pickle

from dwave.system import DWaveSampler

from .draw_utils import graph_2_dot

import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from utils import lp

def get_real_graph(ver=4.1):
    fp = os.path.join(lp, f'qpu_graph_{ver}.txt')
    try:
        with open(fp, 'rb') as file:
            real_graph = pickle.load(file)
    except FileNotFoundError:
        qpu = DWaveSampler(solver={'chip_id':'Advantage_system{}'.format(ver)}, region="na-west-1")
        real_graph = {
            'edges': qpu.edgelist,
            'nodes': qpu.nodelist
        }
        with open(fp, 'wb') as file:
            pickle.dump(real_graph, file)
    return real_graph    
    

if __name__ == '__main__':
    ver = 6.2
    ver = 4.1
    real_graph = get_real_graph(ver)

    dot_file, svg_file = [os.path.join(lp, f'qpu_graph_{ver}{suffix}') for suffix in ['.dot', '.svg']]
    with open(dot_file, "w") as file:
        graph_2_dot(file, None, real_graph, m=16)
    os.system("neato -Tsvg {} -o {}".format(dot_file, svg_file))
    os.system("rm {}".format(dot_file))
