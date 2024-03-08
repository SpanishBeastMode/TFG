import pickle
import json
from collections import defaultdict

import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from utils import lp, improved_lp

def read_cfa_in_binary(i, mp):
    fp = os.path.join(mp, f'CFA{i}')
    with open(fp, 'rb') as file:
        cfa = pickle.load(file)
    return cfa

def write_into_binary(i, cfa, mp):
    fp = os.path.join(mp, f'CFA{i}')
    with open(fp, 'wb') as file:
        pickle.dump(cfa, file)

def serialize_bqm(bqm):
    bqm_in_str = defaultdict(dict)
    for k, v in bqm.items():
        if type(v) == dict: 
            for el, c in v.items():
                bqm_in_str[str(el)] = str(c)
        else:
            bqm_in_str[str(k)] = str(v)
    return bqm_in_str

def write_into_json(i, cfa, mp):
    plmt, bqm = cfa
    lp = os.path.join(mp, f'CFA{i}.json')
    cfa = {
        'plmt': dict(zip(*plmt)),
        'bqm': serialize_bqm(bqm)
    }
    with open(lp, 'w') as file:
        json.dump(cfa, file, indent=4)

if __name__ == '__main__':
    # with open('found_cfa_ver2_0.txt', 'rb') as file: 
    #     ss = pickle.load(file) 
    # dd = {
    #     'plmt': ss[0], 
    #     'bqm': serialize_bqm(ss[1])
    #     }
    # with open('found_cfa_ver2_0.json', 'w') as file: 
    #     json.dump(dd, file, indent=4)

    # multiplier_ver = 1

    # # multiplier_name = f'Multiplier_ver{multiplier_ver}'
    # # suffix = ''
    # suffix = '_to_link_B'
    # multiplier_name = f'Multiplier_ver{multiplier_ver}{suffix}'

    # fp =  os.path.join(lp, f'original_cfhs_ver{multiplier_ver}{suffix}_on_sys4_1.txt')
    fp = os.path.join(improved_lp, f'min_e_eq_2_pfs_for_all_cfh_ver2.txt')
    with open(fp, 'rb') as file:
        data = pickle.load(file)

    # if multiplier_ver == 2 and not multiplier_name.endswith('to_link_B'):
    #     data, = data.values()
    #     cfas = data['cfhs']
    # else:
    #     cfas = [cfa for *cfa, _ in data]
    cfas = data

    # mp = os.path.join(lp, f'{multiplier_name}')
    mp = os.path.join(improved_lp, f'CFA_ver2')
    if not os.path.exists(mp):
        os.mkdir(mp)
    for i, cfa in enumerate(cfas):
        write_into_binary(i, cfa, mp)
        write_into_json(i, cfa, mp)


