
# Multiplier Encoders to Quantum Annealer

## Authors

Jingwen Ding <jingwen.ding@unitn.it>

Giuseppe Spallitta <giuseppe.spallitta@unitn.it>

Roberto Sebastiani <roberto.sebastiani@unitn.it>

## Content

This is the repository containing the code, encodings and data of the experiments reported in the paper:

"Effective Prime Factorization via Quantum Annealing
by Modular Locally-structured Embedding"
by Jingwen Ding, Giuseppe Spallitta & Roberto Sebastiani.

**ArXiv:** arXiv:2310.17574 <https://arxiv.org/abs/2310.17574> 

Under submission for journal publication.

**LICENSE**

The material displayed in this repository (code, encodings and data) refers to the above paper, which is currently under review for journal publication. It is made public only to allow for reproducing the experimental results of the paper. No other usage is allowed without the explicit written permission of the authors.



## Installation

```bash
  pip install -r requirements.txt
```
    

**IMPORTANT WARNING**

The results presented in our paper have been generated using the
D-Wave annealer which was available to us at the time of submission
(i.e. a D-Wave Advantage 4.1), so that the encoded multiplier was
placed in an area of the topology graph in such sistem which was free
of faulty qubits and of faulty couplings.
Thus, if you are using a different system, please make sure that the
area where the multiplier is placed is free of faulty qubits and of
faulty couplings as well, otherwise the process cannot work.

## Usage

# Embedding of the generic multiplier circuit supports:
- three versions constructed with two alternative CFAs, 
- two versions constructed with a single unified CFA
```python
    multiplier = Multiplier_ver1_to_link_B(length, width)
    multiplier = Multiplier_ver2_to_link_B(length, width)
    multiplier = Multiplier_ver3_to_link_B(length, width)

    multiplier = Multiplier_ver2(length, width)
    multiplier = Multiplier_ver3(length, width)
```
Each version has its own library of the CFA(s) to construct the embedding, ./library/.
Moreover, the Multiplier_ver2 has one more library, which are store in ./improved_library/CFA_ver2.
1) ./library/CFA_ver2/: generated based on the definition of the penalty function
2) ./improved_library/CFA_ver2/: generated after the maximal g_min obtained from 1),
with the numbers of the first excited un-satisfying states minimized.

# Four initialization of qubits are supported by the Multiplier_ver2, to solve specific PF problems:
- via D-Wave API:
- via the flux bias of the qubit
- via extra chains
- via ad-hoc encoding of the CFAs involves the qubit to initialize
```python
    vias = ['api', 'flux_biases', 'extra_chains', 'adhoc_encoding']

    length, width, inputs = 8, 8, {'out': 59989}

    initializer = Initializer(inputs, length, width, vias[0], 
                                using_improved_library=True, sys=4.1)
    # using_improved_library = False/True to choose the library for the single CFA model
    # sys=4.1/6.2 supports 17*8bit/13*8bit multiplier respectively
```
# Solving with the D-Wave Advantage system 4.1/6.2 with general strategies
```python
    solver = Solver(inputs, initializer)

    # solving with the standard annealing for different anenaling time Ta
    solver.run_FW_for_diff_Tas(Ta_list=[10, 20])

    # solving with the standard annealing enhanced by thermal relaxation (pause at Sp for Tp time)
    solver.run_pFW_for_diff_Tps(Tp_list=[100], Sp_list=[(38+i)/100 for i in range(4)], Ta=10)

    # solving with the reverse annealign enhanced by thermal relatxition (reverse the annealing to Sp and at Sp pause for Tp, with Ta before/after the pause)
    solver.run_pRV_for_diff_Tps(pause_enhanced=True, Tp_list=[100], Sp_list=[(38+i)/100 for i in range(4)], Ta=10)
    # pause_enhanced=True/False to choose the initial state from the standard annealing or thermal-relaxatio-enhanced annealing
```
# Solving with the systems with iterative strategies
```python
    # iteratively reverse annealing from the obtained low-energy space
    solver.run_BFS_based_iterative_reverse_annealing(pause_enhanced=False, Ta=20)

    # iteratively reverse annealing from a single lower-energy state, with the relatively long pause
    solver.gready_pRV(pause_enhanced=False, filtered=True)
```