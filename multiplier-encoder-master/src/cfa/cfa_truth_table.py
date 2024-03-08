from itertools import product

import os

def get_controlledFA_truth_table(filename='FullAdder_truth_table.txt'):
    fp = os.path.join(os.path.dirname(__file__), filename)

    FA_truth_table = []
    with open(fp, "r") as file:
        for line in file.readlines():
            if line.startswith("#"): continue
            values = line.split("\t")
            values = [int(v) for v in values]
            FA_truth_table.append(values)

    controlledFA_truth_table = [(in1, in2, carry_in, ((in1 and enable) ^ in2 ^ carry_in), ((carry_in and ((in1 and enable) or in2)) or ((in1 and enable) and in2)), enable) 
                                for (in1, in2, carry_in, out, carry_out), enable in product(FA_truth_table, [0, 1])]

    def rearrange(controlled_FA_truth_table):
        return [(in2, in1, enable, carry_in, carry_out, out) for (in1, in2, carry_in, out, carry_out, enable) in controlledFA_truth_table]
    
    controlledFA_truth_table = rearrange(controlledFA_truth_table)
    return controlledFA_truth_table

def get_CFA_ver2_truth_table():
    CFA_ver2_truth_table = [(in2, in1, enable, *latter, enable) 
                                for (in2, in1, enable, *latter) in get_controlledFA_truth_table()]
    return CFA_ver2_truth_table

def get_CFA_ver3_truth_table():
    CFA_ver3_truth_table = [(in2, in1,*latter, in1) 
                                for (in2, in1, *latter) in get_CFA_ver2_truth_table()]
    return CFA_ver3_truth_table

def get_CFA_ver3_to_link_B_truth_table():
    CFA_ver3_to_link_B_truth_table = [(in2, in1,*latter, in1) 
                                for (in2, in1, *latter) in get_controlledFA_truth_table()]
    return CFA_ver3_to_link_B_truth_table