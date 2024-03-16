from utils import close, hamming_distance
from collections import defaultdict
import copy
import json
import math
import pickle
import shutil
import numpy as np
from matplotlib import pyplot as plt

from dimod import BinaryQuadraticModel, SampleSet
from dwave.system import DWaveSampler

import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))


class Solver():
    def __init__(self, inputs, initializer, Ta=20, num_reads=1000):
        self.inputs = inputs
        self.initialized_multiplier = initializer.multiplier
        self.data_path = initializer.data_path

        via = initializer.via

        if via == 'adhoc_encoding':
            self.ancillar_enhanced = initializer.ancillar_enhanced

        model = initializer.run()
        bqm = BinaryQuadraticModel.from_ising(
            *[model[k] for k in ['biases', 'couplings', 'offset']])
        if via in ['api', 'extra_chains']:
            bqm.fix_variables(self.initialized_multiplier.qbt_values)
        self.bqm = bqm

        params = {}
        if via == 'flux_biases':
            params['flux_drift_compensation'] = False
            params['flux_biases'] = self._choose_flux_bias_for_qbts(
                self.initialized_multiplier.qbt_values)

        self.via = via

        # num_reads, self.Ta = 1000, 10 if self.initialized_multiplier.length+self.initialized_multiplier.width <= 8 else 20
        params['num_reads'] = num_reads
        # self.Ta = Ta

        self.params = params
        self.schedule_ptn = lambda Sp, Tp, Ta=None: [
            [Ta, Sp], [Tp + Ta, Sp], [Tp + Ta * 2, 1]]

    def _choose_flux_bias_for_qbts(self, qbt_values, strength=1000):
        # phi0 = self.qpu.properties['anneal_offset_step_phi0']
        PHI0 = 1.4303846404537006e-05
        NUM_QUBITS = 5760
        input_qbts_flux_bias = {qbt: spin * PHI0 *
                                strength for qbt, spin in qbt_values.items()}

        # fb_list = [0]*self.qpu.properties['num_qubits']
        fb_list = [0] * NUM_QUBITS
        for qbt, fb in input_qbts_flux_bias.items():
            fb_list[qbt] = fb
        return fb_list

    def parse_sampleset(self, sampleset):
        statistics = defaultdict(lambda: defaultdict(int))

        l, w = self.initialized_multiplier.length, self.initialized_multiplier.width
        total_excitations = [[0 for _ in range(l)] for _ in range(w)]

        for sample, energy, num_occurrences in sampleset.data(['sample', 'energy', 'num_occurrences']):
            sample = {**sample, **self.initialized_multiplier.qbt_values} if self.via == 'api' or (
                self.via == 'adhoc_encoding' and not self.ancillar_enhanced) else sample
            (A, B, out), (truth_of_cfas, truth_of_chains,
                          truth_of_zero_values), cfa_excitations = self.initialized_multiplier.parse_a_sample_of_all_qbt_values(sample)
            statistics[(A, B, out)][(round(energy, 3), (A*B == out),
                                     truth_of_cfas, truth_of_chains)] += num_occurrences

            for i in range(w):
                for j in range(l):
                    total_excitations[i][j] += num_occurrences * \
                        cfa_excitations[i][j]

        return statistics, total_excitations

    def print_statistics(self, statistics, total_excitations):
        sresult = f'inputs: {self.inputs["out"]
                             }/{self.params["num_reads"]} samples\n'

        def fmt(slist): return '{:>30} {:>10} {:>10} {:>15} {:>15} {:>10}\n'.format(
            *slist)

        sresult += fmt(['[output]', '[energy]', '[A*B=out]',
                       '[truth of cfas]', '[truth of chains]', '[num]'])

        num_in_total = 0
        for config, result in statistics.items():
            for keys, num_occurrences in result.items():
                sresult += fmt([str(config), *keys, num_occurrences])
                num_in_total += num_occurrences
        print(num_in_total)

        # ptn = ' '.join(['{:>5}' for _ in range(self.initialized_multiplier.length)])
        # fmt = lambda i, adder_excitation: f'{i}: {ptn}\n'.format(*adder_excitation)

        # sresult += fmt('-', range(self.initialized_multiplier.width) )
        # for i, adder_excitation in enumerate(total_excitations):
        #    sresult += fmt(i, adder_excitation)
        return sresult

    def _extra_representative_samples(self, sampleset):
        representatives = {}

        i = 0
        for data in sampleset.data(['sample', 'energy', 'num_occurrences']):
            if i == 0 or i == self.params['num_reads']/2 or i == self.params['num_reads']-1 or i in [i*self.params['num_reads']/10 for i in range(1, 10, 2)]:
                representatives[i] = data
            i += data[-1]
        return representatives

    def _read_sampleset(self, fp, old_fp=None):
        try:
            with open(f'{fp}.json', 'r') as file:
                sampleset = SampleSet.from_serializable(
                    json.loads(file.read()))
        except:
            if not old_fp:
                return None
            relative_fp = old_fp.split('/')[-1]
            if f'{relative_fp}.txt' in os.listdir(os.path.dirname(old_fp)):
                with open(f'{old_fp}.txt', 'rb') as file:
                    sampleset = pickle.load(file)

                with open(f'{fp}.json', 'w') as file:
                    file.write(json.dumps(
                        sampleset.to_serializable(), indent=4))
                statistics = self.parse_sampleset(sampleset)
                with open(f'{fp}(statistics)', 'w') as file:
                    file.write(self.print_statistics(*statistics))

        return sampleset

    def _run(self, params=None, old_fp=None):

        qpu = DWaveSampler(Solver=f'chip_id: Advantage_system{self.initialized_multiplier.pegasus_ver}', auto_scale=False)
        sampleset = qpu.sample(self.bqm, **params)

        # pickle is not compatible with Sampleset class anymore -> Sampleset().to_serializable()
        # with open(fp, 'wb') as file:
        #     pickle.dump(sampleset, file)

        statistics = self.parse_sampleset(sampleset)
        return statistics

    def forward_annealing(self, Ta, anneal_dir=None, use_old_data=None):
        if not os.path.exists(anneal_dir):
            os.makedirs(anneal_dir)

        sche = f'(Ta={Ta}us)'
        fp = os.path.join(anneal_dir, f'factoring_{self.inputs["out"]}_via_{sche}')

        old_fp = None
        if use_old_data:
            out = self.inputs['out']
            filename = f'FW{Ta}us_factoring_{out}_with_{self.initialized_multiplier.length}*{
                self.initialized_multiplier.width}bit_of_CFA_{self.initialized_multiplier.cfa_index}_symmetric_None'
            old_fp, statistic_fp = [os.path.join(
                anneal_dir, f'{prefix}{filename}') for prefix in ['sampleset_', '']]
            statistic_fp = statistic_fp + '.txt'
            if os.path.exists(statistic_fp):
                os.remove(statistic_fp)

        params = copy.deepcopy(self.params)
        params['annealing_time'] = Ta

        print(f'Ta={Ta}us: \n\tfp = {fp}, \n\tparams["annealing_time"] = {params["annealing_time"]}')
        return self._run(fp, params, old_fp=old_fp)

    def annealing_enhanced_by_pause(self, params, Tp=1, Sp_list=None, Ta=10, anneal_direction=None, anneal_dir=None, use_old_data=None, rounds=None):

        if not Sp_list:
            if anneal_direction == 'forward':
                fSp, step = 30, 3
                counter = 6 if Tp >= 100 else 4
            else:
                fSp, step = 32, 1  # in total region tunning
                counter = 14
            Sp_list = [(fSp+i)/100 for i in range(0, counter*step, step)]

        print(f'with Ta={Ta}us:')
        schedule_from = [0, 0] if anneal_direction == 'forward' else [0, 1]
        rs = {}
        for Sp in Sp_list:
            sche = f'(Ta={Ta}us, Tp={Tp}us, Sp={Sp})'
            # params = copy.deepcopy(self.params)
            params['anneal_schedule'] = [
                schedule_from, *self.schedule_ptn(Sp, Tp, Ta)]

            print(f'\tSp={Sp}|Tp={Tp}us:, \n\t\tparams["anneal_schedule"] = {
                  params["anneal_schedule"]}')
            rs[Sp] = self._run(params)
        return rs

    def run_FW_for_diff_Tas(self, Ta_list=[10, 30, 50, 100], **kwargs):
        anneal_dir = os.path.join(self.data_path, 'forward_annealing')

        results = {}
        for Ta in Ta_list:
            results[Ta] = self.forward_annealing(
                Ta, anneal_dir=anneal_dir, **kwargs)

        # self.plot_forward_annealing_results(results, anneal_dir)
        return results

    def run_pFW_for_diff_Tps(self, Tp_list=[1, 10, 30, 50, 100], Sp_list=None, anneal_direction='forward', Ta=10):
        params = copy.deepcopy(self.params)

        anneal_dir = os.path.join(self.data_path, f'{
                                  anneal_direction}_annealing_enhanced_by_pause', f'with_Ta={Ta}us')

        results = {}
        for Tp in Tp_list:
            results[Tp] = self.annealing_enhanced_by_pause(params, Tp=Tp, Sp_list=Sp_list,
                                                           Ta=Ta, anneal_direction=anneal_direction, anneal_dir=anneal_dir)

        # self.plot_pause_enhanced_forward_annealing_results(results, anneal_dir, Ta)
        return results

    def _get_global_minima_from_pure_FW(self, Ta_list=[10, 30, 50, 100], **kwargs):
        anneal = self.run_FW_for_diff_Tas
        results = anneal(Ta_list=Ta_list, **kwargs)

        grdstate = self._find_groundstate()

        to_optima = None
        for ta, representatives in sorted(results.items(), key=lambda item: item[0]):
            minima = representatives[0]

            optima = hamming_distance(
                minima.sample, grdstate.sample) if grdstate else minima.energy
            this_optima = (ta, minima)
            if to_optima is None:
                to_optima = this_optima
            else:
                Ta, Minima = to_optima
                Optima = hamming_distance(
                    Minima.sample, grdstate.sample) if grdstate else Minima.energy
                if close(optima, Optima):
                    if ta < Ta:
                        to_optima = this_optima
                elif optima < Optima:
                    to_optima = this_optima
        return to_optima

    def _choose_minima_for_next_pRV_guided_by_moves(self, fminima, Sp_results, moves_list, bound=300):
        reduced_Sp_results = {}
        for Sp, representative in Sp_results.items():
            minima = representative[0]
            moves = hamming_distance(minima.sample, fminima.sample)

            te, fe = minima.energy, fminima.energy
            # if (close(te, fe) and int(fe) > 6) or te > fe
            rte, rfe = round(te), round(fe)
            if ((rte == rfe) and rfe > 6) or rte > rfe:
                continue
            reduced_Sp_results[Sp] = (moves, minima)

        if not moves_list:
            if all([moves >= bound for moves, _ in reduced_Sp_results.values()]):
                minimal_moves = None
                ordered_Sp_results = sorted(
                    reduced_Sp_results.items(), key=lambda item: item[1][0])
                (Sp, (min_moves, minima)), *_ = ordered_Sp_results
                return (Sp, minima)

        ordered_Sp_results = sorted(
            reduced_Sp_results.items(), key=lambda item: -item[1][0])
        sm = sum((moves for _, moves in moves_list))
        for Sp, (moves, minima) in ordered_Sp_results:
            this_minima = (Sp, minima)
            if not moves_list:
                return this_minima
            if sm + moves <= bound and Sp > moves_list[-1][0]:
                return this_minima
        return None

    def _get_minTp_to_reach_minima(self, Tp_results, subject='energy'):
        grdstate = self._find_groundstate()

        to_optima = None
        for tp, sp_results in Tp_results.items():
            for sp, representatives in sp_results.items():
                minima = representatives[0]

                optima = hamming_distance(
                    minima.sample, grdstate.sample) if subject == 'distance' and grdstate else minima.energy
                this_optima = ((tp, sp), minima)
                if to_optima is None:
                    to_optima = this_optima
                else:
                    (Tp, Sp), Minima = to_optima
                    Optima = hamming_distance(
                        Minima.sample, grdstate.sample) if subject == 'distance' and grdstate else Minima.energy
                    if close(optima, Optima):
                        if tp < Tp:
                            to_optima = this_optima
                    elif optima < Optima:
                        # elif round(optima) < round(Optima):
                        to_optima = this_optima
        return to_optima

    def _get_global_minima_from_pause_enhanced_FW(self, Tp_list, Sp_list, Ta_list=[10, 30, 50, 100]):
        anneal = self.run_pFW_for_diff_Tps

        to_global_minima = None
        for ta in Ta_list:
            Tp_results = anneal(Tp_list, Sp_list, Ta=ta)

            to_minima = self._get_minTp_to_reach_minima(Tp_results)
            (tp, sp), minima = to_minima
            to_this_minima = ((ta, tp, sp), minima)

            if to_global_minima is None:
                to_global_minima = to_this_minima
            else:
                (Ta, Tp, Sp), Minima = to_global_minima
                if close(minima[1], Minima[1]):
                    if (2*ta + tp) < (2*Ta + Tp):
                        to_global_minima = to_this_minima
                elif minima[1] < Minima[1]:
                    to_global_minima = to_this_minima
        return to_global_minima

    def _choose_initial_state_for_pRV(self, pause_enhanced=False, Tp_list=None, Sp_list=None, Ta_list=[10], **kwargs):
        if pause_enhanced:
            to_minima = self._get_global_minima_from_pause_enhanced_FW(
                Tp_list, Sp_list, Ta_list=Ta_list)
            (ta, tp, sp), minima = to_minima
            to_dir = f'from_pFW_(Ta={ta}us, Tp={tp}us, Sp={
                sp}, minima={round(minima.energy, 3)})'
        else:
            to_minima = self._get_global_minima_from_pure_FW(Ta_list, **kwargs)
            ta, minima = to_minima
            to_dir = f'from_FW_(Ta={ta}us, minima={round(minima.energy, 3)})'
        return to_minima, to_dir

    def run_pRV_for_diff_Tps(self, pause_enhanced=False,
                             Tp_list=[1, 10, 30, 50, 100], Sp_list=None,
                             anneal_direction='reverse', Ta=10, filtered=True):
        from_minima, to_dir = self._choose_initial_state_for_pRV(
            pause_enhanced, Tp_list=Tp_list, Sp_list=Sp_list)
        _, minima = from_minima

        anneal_dir = os.path.join(self.data_path, f'{anneal_direction}_annealing_enhanced_by_pause',
                                  to_dir,
                                  f'with_Ta={Ta}us')

        params = copy.deepcopy(self.params)
        # params['annealing_schedule'] = [[0, 1], *pause_at_Sp_for_Tp(Sp, Tp)]
        params['reinitialize_state'] = False
        params['initial_state'] = minima.sample

        results = {}
        for Tp in Tp_list:
            results[Tp] = self.annealing_enhanced_by_pause(
                params, Tp=Tp, Sp_list=Sp_list, Ta=Ta, anneal_direction=anneal_direction, anneal_dir=anneal_dir)

        # self.plot_pause_enhanced_reverse_annealing_results(results, anneal_dir, Ta, from_minima, abbr='pRV', filtered=filtered)
        return results

    def _choose_a_lower_energy_space_for_initial_states(self, rounds, initial_state, sp_fp, hard_energies={1: 8, 2: 6, 3: 6, 4: 0}):
        fe, minima, _ = initial_state

        reduced_lower_energy_space = defaultdict(list)
        sampleset = self._read_sampleset(sp_fp)
        for sample, energy, num_occurrences in sampleset.data(['sample', 'energy', 'num_occurrences'], sorted_by='energy'):
            if round(energy) >= round(fe):
                break

            dist = hamming_distance(minima, sample)
            # reduced_lower_energy_space[round(energy, 3)].append((sample, dist, num_occurrences))
            reduced_lower_energy_space[round(energy, 3)].append((sample, dist))

        initial_state_list = []
        min_e = min(reduced_lower_energy_space.keys())
        for e, elist in reduced_lower_energy_space.items():
            if int(e) <= hard_energies[rounds]:
                continue

            if initial_state_list and int(e) in [int(energy) for energy, *_ in initial_state_list]:
                continue

            # counter = 0
            # median = int(sum([num_occurrences for *_, num_occurrences in elist])/2)
            # for dist, sample, num_occurrences in sorted(elist, key=lambda item: item[1]):
            #     if counter >= median:
            #         initial_state_list.append((e, sample, dist))
            #         break
            #     counter += num_occurrences

            median = int(len(elist)/2)
            initial_state = sorted(elist, key=lambda item: item[1])[median]
            initial_state_list.append((e, *initial_state))
        return initial_state_list

    def _iterative_pRV(self, minima, Tp, **kwargs):
        params = copy.deepcopy(self.params)
        params['reinitialize_state'] = False
        params['initial_state'] = minima

        Sp_results = self.annealing_enhanced_by_pause(params, Tp=Tp, **kwargs)
        return Sp_results

    def run_BFS_based_iterative_reverse_annealing(self, pause_enhanced=False,
                                                  Tp_list=[1, 10, 30, 50, 100], Sp_list=[(30+i)/100 for i in range(10)],
                                                  Ta=20, anneal_direction='reverse',
                                                  use_old_data=True):
        anneal_dir = os.path.join(
            self.data_path, 'BFS_based_iterated_reverse_annealing')

        l, w = self.initialized_multiplier.length, self.initialized_multiplier.width
        hard_energies = {1: 8, 2: 6, 3: 2, 4: 0} if l == 16 else {
            # tunable for the choice of IRV's initial state
            1: 8, 2: 6, 3: 6, 4: 0} if l == 14 else {}
        tp_list = [i for i in range(1, 10, 2)] if l == 16 else [i for i in range(
            # tunable for each iteration
            1, 10, 2)] + [20, 100] if l == 14 else copy.deepcopy(Tp_list)

        from_minima, to_relative_dir = self._choose_initial_state_for_pRV(
            pause_enhanced, Tp_list, Sp_list, [Ta], use_old_data=use_old_data)
        ta, minima = from_minima
        fe = round(minima.energy, 3)
        initial_state = (fe, minima.sample, None)
        initial_state_list = [initial_state]

        iterative_results_for_initial_states = defaultdict(dict)
        rounds = 1
        iterative_dir = os.path.join(
            anneal_dir, f'IRV{rounds}_' + to_relative_dir)
        kwargs = {
            'Sp_list': Sp_list,
            'Ta': Ta,
            'anneal_direction': anneal_direction,
            'use_old_data': use_old_data
        }
        MINIMAs = {0: ((ta, None, None), minima)}
        while tp_list:
            tp = tp_list.pop(0)
            kwargs['rounds'] = rounds

            colors = plt.colormaps['viridis'](
                np.linspace(0, 1, len(initial_state_list)))
            colors = colors[::-1]

            iterative_results_per_tp = {}
            for initial_state in initial_state_list[::-1]:
                fe, minima, moves = initial_state

                if rounds > 1:
                    to_relative_dir = f'from_a_low_energy_state-(e={fe}, median_moves={
                        moves})'

                to_abs_dir = os.path.join(
                    iterative_dir, '' if rounds == 1 else to_relative_dir)

                kwargs['anneal_dir'] = os.path.join(
                    to_abs_dir, f'with_Ta={Ta}us')

                Sp_results = self._iterative_pRV(minima, tp, **kwargs)
                to_optima = self._get_minTp_to_reach_minima({tp: Sp_results})
                (_, sp), optima = to_optima
                te = round(optima.energy)

                iterative_results_per_tp[(fe, moves)] = (minima, Sp_results)
                # self.plot_BFS_based_iterative_results_per_tp(MINIMAs[rounds-1], rounds, tp, iterative_results_per_tp, to_dir=iterative_dir, colors=colors)

                iterative_results_for_initial_states[(
                    fe, moves)][tp] = (minima, Sp_results)
                # if rounds > 1:
                #     self.plot_BFS_based_iterative_results_for_initial_states(iterative_results_for_initial_states, MINIMAs[rounds-1], rounds, to_dir=iterative_dir)

                if te < 2:
                    return 'found'

                if te < round(MINIMAs[rounds-1][-1].energy):
                    if rounds == 3 and not te < 4:
                        continue
                    MINIMAs[rounds] = ((Ta, tp, sp), optima)

                    rounds += 1
                    iterative_dir = os.path.join(to_abs_dir, f'IRV{rounds}_from_the_IRV{
                                                 rounds-1}(Tp={tp}us, Sp={sp})')
                    sp_fp = os.path.join(kwargs['anneal_dir'], f'pause_of_Tp={tp}us', f'factoring_{
                                         self.inputs["out"]}_via_(Ta={Ta}us, Tp={tp}us, Sp={sp})')
                    initial_state_list = self._choose_a_lower_energy_space_for_initial_states(
                        rounds, initial_state, sp_fp, hard_energies=hard_energies)
                    if rounds == 4:
                        initial_state_list = initial_state_list[0:1]

                    if l == 16:
                        tp_list = [1, 10, 50, 70, 100] if rounds >= 3 else copy.deepcopy(
                            Tp_list)
                    elif l == 14:
                        tp_list = [1, 10, 50] if rounds >= 3 else copy.deepcopy(
                            Tp_list)
                    else:
                        tp_list = copy.deepcopy(Tp_list)

                    iterative_results_for_initial_states = defaultdict(dict)
                    break

    def _recursive_pRV(self, rounds, from_minima, recursive_dir, Tp_list=None, Ta=None, moves_list=None,
                       abbr=None, filtered=None, **kwargs):
        _, minima = from_minima

        params = copy.deepcopy(self.params)
        params['reinitialize_state'] = False
        params['initial_state'] = minima.sample

        anneal_dir = os.path.join(recursive_dir, f'with_Ta={Ta}us')
        results = {}
        for Tp in Tp_list:
            Sp_results = self.annealing_enhanced_by_pause(
                params, Tp=Tp, Ta=Ta, anneal_dir=anneal_dir, **kwargs)
            results[Tp] = Sp_results
            # self.plot_pause_enhanced_reverse_annealing_results(results, anneal_dir, Ta, from_minima,
            #                                                 abbr='pRV', rounds=rounds, filtered=filtered)
            # Sp, minima = self._choose_minima_for_next_pRV_guided_by_moves(minima, Sp_results, moves_list, bound=300)

            # return (Tp, Sp), minima

        # self.plot_pause_enhanced_reverse_annealing_results(results, anneal_dir, Ta, from_minima,
        #                                                     abbr='pRV', rounds=rounds, filtered=filtered)
        return self._get_minTp_to_reach_minima(results)
        # return (Tp, Sp), minima

    def gready_pRV(self, pause_enhanced=False, Tp_list=[1, 10, 30, 50, 100], Sp_list=None, Ta=10, anneal_direction='reverse', filtered=False):
        _, cfa_ising_model = self.initialized_multiplier.cfa
        gap = cfa_ising_model['g_min']

        recursive_dir = os.path.join(
            self.data_path, f'{anneal_direction}_annealing_enhanced_by_pause')

        out = self.inputs['out']
        self.reversed_grdstate_path = os.path.join(
            recursive_dir, f'factoring_{out}_reversed_grdstate_path.json')

        from_minima, to_dir = self._choose_initial_state_for_pRV(
            pause_enhanced, Tp_list, Sp_list)
        abbr = 'pFW' if pause_enhanced else 'FW'

        def _chose_Tps_and_Sps(rounds):
            # how to choose Tp_list, Sp_list, and a state for initializing the next round of reverse annealing is not determined
            # Here they are chosen manually, like the following spcified tp_list(/sp_list) together with self._get_minTp_to_reach_minima()

            l, w = self.initialized_multiplier.length, self.initialized_multiplier.width
            if l == 15 and w == 8:
                if rounds <= 4:
                    search_path = {
                        1: ([1], [0.4]),
                        2: ([100], [0.43]),
                        3: ([200], [0.43]),
                        4: ([200], [0.44])
                    }
                    tp_list, sp_list = search_path[rounds]
                else:
                    tp_list = [100, 200]
                    sp_list = [(40+i)/100 for i in range(6)]
            else:
                if rounds == 1 and l >= 15:
                    tp_list = [1]
                    sp_list = [(36+i)/100 for i in range(10)]
                else:
                    tp_list = [100, 200]
                    sp_list = [(40+i)/100 for i in range(6)]

            return tp_list, sp_list

        kwargs = {
            'Tp_list': Tp_list,
            'Sp_list': Sp_list,
            'Ta': Ta,
            'anneal_direction': anneal_direction,
            'abbr': abbr,
            'filtered': filtered
        }
        kwargs['moves_list'] = []

        counter = 0
        rounds = 1
        recursive_dir = os.path.join(recursive_dir, to_dir)
        while True:
            if rounds > 5:
                break

            kwargs['Tp_list'], kwargs['Sp_list'] = _chose_Tps_and_Sps(rounds)

            (tp, sp), minima = self._recursive_pRV(
                rounds, from_minima, recursive_dir, **kwargs)

            energy = minima.energy
            if not (energy >= gap or close(energy, gap)):
                if not os.path.exists(self.reversed_grdstate_path):
                    grdstate_path_info = {
                        out: os.path.join(recursive_dir, f'with_Ta={Ta}us',
                                          f'pause_of_Tp={tp}us',
                                          f'factoring_{out}_via_(Ta={Ta}us, Tp={tp}us, Sp={sp})')

                    }
                    with open(self.reversed_grdstate_path, 'w') as file:
                        json.dump(grdstate_path_info, file, indent=4)
                break

            # if round(energy) >= round(from_minima[-1].energy):
            #     counter += 1
            #     continue
            # counter += 1
            moves = hamming_distance(from_minima[-1].sample, minima.sample)
            kwargs['moves_list'].append((sp, moves))

            rounds += 1
            from_minima = ((Ta, tp, sp), minima)
            abbr = 'pRV'
            recursive_dir = os.path.join(recursive_dir, f'from_{abbr}_(Ta={Ta}us, Tp={
                                         tp}us, Sp={sp}, minima={round(energy, 3)})')
        return from_minima

    def _find_reversed_grdstate_path(self):
        anneal_approaches = ['forward_annealing'] + [
            f'{anneal_direction}_annealing_enhanced_by_pause' for anneal_direction in ['forward', 'reverse']]

        for anneal in anneal_approaches:
            anneal_dir = os.path.join(self.data_path, anneal)
            for root, dirs, files in os.walk(anneal_dir):
                for filename in files:
                    if filename.endswith(f'{self.inputs["out"]}_reversed_grdstate_path.json'):
                        return os.path.join(root, filename)
        return None

    def _find_groundstate(self):
        if hasattr(self, 'grdstate'):
            return self.grdstate

        reversed_grdstate_path = self._find_reversed_grdstate_path()
        if not reversed_grdstate_path:
            return None

        self.reversed_grdstate_path = reversed_grdstate_path
        with open(reversed_grdstate_path, 'r') as file:
            grdstate_path_info = json.load(file)
        path = grdstate_path_info[str(self.inputs['out'])]
        sampleset = self._read_sampleset(path)
        grdstate = sampleset.first

        self.grdstate = grdstate
        return grdstate

    def _choose_plot_style_for_samples(self, i):
        label, marker = ('minima', '^-') if i == 0 else \
            ('maxima', 'v-') if i == self.params['num_reads']-1 else \
            ('median', 's-') if i == self.params['num_reads']/2 else \
            (f'{i}th', '.--')
        return label, marker

    def _plot_energies_against_Ta(self, ax, results):
        ax.set_xlabel('Ta')
        ax.set_ylabel('Energy')

        num_reads = self.params['num_reads']
        ten_percent = num_reads/10
        ks = [0] + [ten_percent*k for k in range(1, 10, 2)] + [num_reads-1]

        colors = plt.colormaps['viridis'](np.linspace(0, 1, len(ks)))

        Ta_list = list(results.keys())
        for i, k in enumerate(ks):
            elist = [representatives[k][1]
                     for representatives in results.values()]
            label, marker = self._choose_plot_style_for_samples(int(k))
            ax.plot(Ta_list, elist, marker, label=label, color=colors[i])
            if k != 0:
                continue
            for Ta, e in zip(Ta_list, elist):
                ax.annotate(round(e, 3), [Ta, e])

    def _plot_distances_against_Ta(self, ax, results, wrpt_minima, **kwargs):
        wrpt, minima = wrpt_minima
        assert wrpt == 'grdstate'

        ax.set_xlabel('Ta')
        ylabel = 'Distance from the ground state'
        ax.set_ylabel(ylabel)

        num_reads = self.params['num_reads']
        ten_percent = num_reads/10
        ks = [0] + [ten_percent*k for k in range(1, 10, 2)] + [num_reads-1]

        colors = plt.colormaps['viridis'](np.linspace(0, 1, len(ks)))

        Ta_list = list(results.keys())
        for i, k in enumerate(ks):
            distance_list = [hamming_distance(
                representatives[k].sample, minima.sample) for representatives in results.values()]
            label, marker = self._choose_plot_style_for_samples(int(k))
            ax.plot(Ta_list, distance_list, marker,
                    label=label, color=colors[i])
            if k != 0:
                continue
            for Ta, dist in zip(Ta_list, distance_list):
                ax.annotate(dist, [Ta, dist])

    def _plot_energies_against_Sp(self, ax, results, filtered=False):
        ax.set_xlabel('Sp')
        ax.set_ylabel('Energy')

        num_reads = self.params['num_reads']
        ks = [0] if filtered else [0, num_reads/10, num_reads/2, num_reads-1]
        colors = plt.colormaps['viridis'](
            np.linspace(0, 1, len(results.items())))
        colors = colors[::-1]

        i = 0
        for Tp, rs in results.items():
            Sp_list = list(rs.keys())

            for k in ks:
                elist = [representatives[k][1]
                         for representatives in rs.values()]
                label, marker = self._choose_plot_style_for_samples(int(k))
                ax.plot(Sp_list, elist, marker, label=f'(Tp={
                        Tp}us, {label})', color=colors[i])
                if k != 0:
                    continue
                for Sp, e in zip(Sp_list, elist):
                    ax.annotate(round(e, 3), [Sp, e])
            i += 1

    def _plot_distances_against_Sp(self, ax, results, wrpt_minima=None, initialstate_far_from_grdstate=None, filtered=False):
        wrpt, minima = wrpt_minima

        ax.set_xlabel('Sp')
        ylabel = 'Distance from the ground state' if wrpt == 'grdstate' else 'Moves from the initial state'
        ax.set_ylabel(ylabel)

        num_reads = self.params['num_reads']
        ks = [0] if filtered else [0, num_reads/10, num_reads/2, num_reads-1]
        colors = plt.colormaps['viridis'](
            np.linspace(0, 1, len(results.items())))
        colors = colors[::-1]

        i = 0
        for Tp, rs in results.items():
            Sp_list = list(rs.keys())

            for k in ks:
                distance_list = [hamming_distance(
                    representatives[k].sample, minima.sample) for representatives in rs.values()]
                label, marker = self._choose_plot_style_for_samples(int(k))
                ax.plot(Sp_list, distance_list, marker,
                        label=f'(Tp={Tp}us, {label})', color=colors[i])
                if k != 0:
                    continue
                for Sp, dist in zip(Sp_list, distance_list):
                    ax.annotate(dist, [Sp, dist])
            i += 1

        if initialstate_far_from_grdstate:
            ax.axhline(y=initialstate_far_from_grdstate,
                       color='red', linestyle='dotted')

    def _plot_results_wrpt(self, fp, results, title, wrpt_minima=None, against='Sp', **kwargs):
        # if os.path.exists(fp):
        #     return

        if wrpt_minima:
            cols = 2
            fig, axs = plt.subplots(1, cols, figsize=(
                16, 6), constrained_layout=True)
            if against == 'Sp':
                eval(f'self._plot_energies_against_{against}')(
                    axs[0], results, filtered=kwargs['filtered'])
            else:
                eval(f'self._plot_energies_against_{against}')(axs[0], results)
            eval(f'self._plot_distances_against_{against}')(
                axs[1], results, wrpt_minima, **kwargs)
        else:
            fig, ax = plt.subplots(1, 1, figsize=(
                16, 6), constrained_layout=True)
            eval(f'self._plot_energies_against_{against}')(ax, results)

        plt.legend(bbox_to_anchor=(1.05, 1),
                   loc='upper left', borderaxespad=0.)
        plt.tight_layout()
        fig.suptitle(title)
        plt.savefig(fp)
        plt.close()

    def plot_forward_annealing_results(self, results, anneal_dir):
        filename = 'anneal_for_diff_Tas.png'

        title = 'Forward annealing for different times (Ta)'

        grdstate = self._find_groundstate()
        if grdstate:
            wrpt_minima = ('grdstate', grdstate)
            prefix = '(energy-distance)'
        else:
            wrpt_minima = None
            prefix = ''
        fp = os.path.join(anneal_dir, prefix + filename)
        self._plot_results_wrpt(
            fp, results, title, wrpt_minima=wrpt_minima, against='Ta')

    def plot_pause_enhanced_forward_annealing_results(self, results, anneal_dir, Ta):
        filename = f'reprs_pFW_for_diff_Tps|Ta={Ta}us.png'

        title = f'Forward annealing enhanced by different pauses (Tp|Ta={
            Ta}us)'

        grdstate = self._find_groundstate()
        if grdstate:
            fp = os.path.join(anneal_dir, '(energy-distance)' + filename)
            self._plot_results_wrpt(
                fp, results, title, wrpt_minima=('grdstate', grdstate))
        else:
            fp = os.path.join(anneal_dir, filename)
            self._plot_results_wrpt(fp, results, title)

    def plot_pause_enhanced_reverse_annealing_results(self, results, anneal_dir, Ta, from_minima,
                                                      abbr='pRV', rounds=None,
                                                      filtered=False):
        filename = f'reprs_{abbr}_for_diff_Tps|Ta={Ta}us.png'
        filtered_label = '[filtered-largerTp]' if filtered else ''

        sche, minima = from_minima
        min_energy = round(minima.energy, 3)
        if type(sche) == int:
            ta = sche
            how_anneal = f'From a minima (e={min_energy}) ' + \
                f'obtained by FW(Ta={ta}us) \n' + \
                'reverse annealing'
        elif len(sche) == 3:
            ta, tp, sp = sche
            abbr_at_rounds = 'pFW' if abbr == 'pRV' and rounds == 1 else abbr
            how_anneal = f'From a minima (e={min_energy}) ' + \
                f'obtained by {abbr_at_rounds}(Ta={ta}us, Tp={tp}us, Sp={sp}), \n' + \
                'reverse annealing'
        title = f'{how_anneal} enhanced by different pauses (Tp|Ta={Ta}us)'

        kwargs = {'filtered': filtered}
        grdstate = self._find_groundstate()
        if grdstate:
            fp_wrpr_gd = os.path.join(
                anneal_dir, filtered_label + '(energy-distance)' + filename)
            kwargs['initialstate_far_from_grdstate'] = hamming_distance(
                minima.sample, grdstate.sample)
            self._plot_results_wrpt(fp_wrpr_gd, results, title, wrpt_minima=(
                'grdstate', grdstate), **kwargs)

        fp = os.path.join(anneal_dir, filtered_label +
                          '(energy-moves)' + filename)
        self._plot_results_wrpt(fp, results, title, wrpt_minima=(
            'local_minima', minima), **kwargs)

        fp3D = os.path.join(anneal_dir, filtered_label +
                            '(energy-moves-in-3D)' + filename)
        self._plot_results_wrpt_in_3D(
            fp3D, results, title, wrpt_minima=('local_minima', minima))

    def _plot_results_wrpt_in_3D(self, fp, results, title, wrpt_minima=None, against='Tp'):
        if os.path.exists(fp):
            return

        wrpt, minima = wrpt_minima
        assert wrpt == 'local_minima'

        cols = 2
        fig = plt.figure(figsize=(13, 6.3), layout="constrained")
        axs = [fig.add_subplot(1, cols, j+1, projection='3d')
               for j in range(cols)]

        for i, ax in enumerate(axs):
            ax.set_xlabel('Sp')
            ax.set_ylabel('Tp')
            zlabel = 'Energy' if i == 0 else 'Moves from the initial state'
            ax.set_zlabel(zlabel)

        num_reads = self.params['num_reads']
        # ks = [0, num_reads/10, num_reads/2, num_reads-1]
        ks = [0, num_reads/2]

        dual_results = defaultdict(lambda: defaultdict(list))
        Tp_list = list(results.keys())
        for rs in results.values():
            for Sp, representatives in rs.items():
                for k in ks:
                    kth_minima = representatives[k]
                    energy = round(kth_minima.energy, 3)
                    distance = hamming_distance(
                        kth_minima.sample, minima.sample)
                    dual_results[Sp][k].append((energy, distance))

        colors = plt.colormaps['viridis'](
            np.linspace(0, 1, len(dual_results.items())))

        kwargs = {'zdir': 'x'}
        i = 0
        for Sp, k_reprs in dual_results.items():
            kwargs['zs'] = Sp
            kwargs['color'] = colors[i]
            for k, kth_results in k_reprs.items():
                elist = [energy for energy, _ in kth_results]
                distance_list = [dist for _, dist in kth_results]

                label, marker = self._choose_plot_style_for_samples(int(k))
                kwargs['label'] = f'(Sp={Sp}, {label})'

                marker = marker + '-' if label == 'median' else marker
                axs[0].plot(Tp_list, elist, marker, **kwargs)
                axs[1].plot(Tp_list, distance_list, marker, **kwargs)
            i += 1

        plt.legend(bbox_to_anchor=(1.15, 1),
                   loc='upper left', borderaxespad=0.)
        plt.tight_layout()
        fig.suptitle(title)
        plt.savefig(fp)
        plt.close()

    def _plot_energies_against_Sp_for_all_initial_states_with_fixed_Tp(self, ax, iterative_results_per_tp, colors, k=0, prev_optima=None):
        ax.set_xlabel('Sp')
        ax.set_ylabel('Energy')

        i = 0
        for (fe, median_moves), (_, rs) in iterative_results_per_tp.items():
            Sp_list = list(rs.keys())
            min_e_list = [
                representatives[k].energy for representatives in rs.values()]
            color = colors[i]
            ax.plot(Sp_list, min_e_list, '.-',
                    label=f'(moves={median_moves}|e={fe})', color=color)
            for Sp, e in zip(Sp_list, min_e_list):
                ax.annotate(round(e, 3), [Sp, e])

            fe = round(fe, 3)
            ax.axhline(y=fe, color=color, linestyle='dotted')
            ax.annotate(fe, [ax.get_xlim()[-1], fe])
            i += 1

        minima = round(prev_optima[-1].energy, 3)
        if minima not in [round(fe, 3) for fe, _ in iterative_results_per_tp.keys()]:
            ax.axhline(y=minima, color='red', linestyle='dotted')
            ax.annotate(minima, [ax.get_xlim()[-1], minima])

    def _plot_distances_against_Sp_for_all_initial_states_with_fixed_Tp(self, ax, iterative_results_per_tp, colors, wrpt=None, k=0, prev_optima=None):
        ax.set_xlabel('Sp')
        ylabel = 'Distance from the ground state' if wrpt == 'grdstate' else 'Moves from the initial state'
        ax.set_ylabel(ylabel)

        grdstate = self._find_groundstate()

        i = 0
        for (fe, median_moves), (fsample, rs) in iterative_results_per_tp.items():
            Sp_list = list(rs.keys())
            target = grdstate.sample if wrpt == 'grdstate' else fsample
            distance_list = [hamming_distance(
                target, representatives[k].sample) for representatives in rs.values()]
            color = colors[i]
            ax.plot(Sp_list, distance_list, '.-',
                    label=f'(moves={median_moves}|e={fe})', color=color)
            for Sp, dist in zip(Sp_list, distance_list):
                ax.annotate(dist, [Sp, dist])
            if grdstate:
                abs_dist = hamming_distance(fsample, grdstate.sample)
                ax.axhline(y=abs_dist, color=color, linestyle='dotted')
                ax.annotate(abs_dist, [ax.get_xlim()[-1], abs_dist])
            i += 1

        minima = round(prev_optima[-1].energy, 3)
        if grdstate and minima not in [round(fe, 3) for fe, _ in iterative_results_per_tp.keys()]:
            minima_dist = hamming_distance(
                prev_optima[-1].sample, grdstate.sample)
            ax.axhline(y=minima_dist, color='red', linestyle='dotted')
            ax.annotate(minima_dist, [ax.get_xlim()[-1], minima_dist])

    def plot_BFS_based_iterative_results_per_tp(self, prev_optima, rounds, tp, iterative_results_per_tp, to_dir=None, colors=None):
        grdstate = self._find_groundstate()

        prefix, cols = ('(energy-moves-distance)',
                        3) if grdstate else ('(energy-moves)', 2)
        # prefix = '(for all Tps)'
        (Ta, Tp, Sp), optima = prev_optima
        e = round(optima.energy, 3)
        title = f'IRV{rounds}(Tp={tp}us) when initialized in different states obtained by IRV{
            rounds-1}((Tp={Tp}us, Sp={Sp}, minima={e})|Ta={Ta}us)'
        fp = os.path.join(to_dir, f'{prefix}{title}.png')

        if os.path.exists(fp):
            return

        fig, axs = plt.subplots(1, cols, figsize=(
            20, 6), constrained_layout=True)
        fig.suptitle(title)

        *args, kwargs = (iterative_results_per_tp, colors,
                         {'prev_optima': prev_optima})
        for i, (ax, subject) in enumerate(zip(axs, ['energies', 'distances', 'distances'])):
            kwargs = copy.deepcopy(kwargs)
            if i == 2:
                kwargs['wrpt'] = 'grdstate'
            eval(f'self._plot_{subject}_against_Sp_for_all_initial_states_with_fixed_Tp')(
                ax, *args, **kwargs)

        plt.legend(bbox_to_anchor=(1.05, 1),
                   loc='upper left', borderaxespad=0.)
        plt.tight_layout()
        plt.savefig(fp)
        plt.close()

    def _plot_energies_against_Sp_for_all_Tps_with_fixed_initial_state(self, ax, iterative_results_per_initial_state, k=0, prev_optima=None):
        ax.set_xlabel('Sp')
        ax.set_ylabel('Energy')

        colors = plt.colormaps['viridis'](np.linspace(
            0, 1, len(iterative_results_per_initial_state.keys())))
        colors = colors[::-1]

        i = 0
        for tp, (_, rs) in iterative_results_per_initial_state.items():
            Sp_list = list(rs.keys())
            min_e_list = [
                representatives[k].energy for representatives in rs.values()]
            color = colors[i]
            ax.plot(Sp_list, min_e_list, '.-', label=f'Tp={tp}us', color=color)
            for Sp, e in zip(Sp_list, min_e_list):
                ax.annotate(round(e, 3), [Sp, e])
            i += 1

    def _plot_distances_against_Sp_for_all_Tps_with_fixed_initial_state(self, ax, iterative_results_for_Tps, wrpt=None, k=0, prev_optima=None):
        ax.set_xlabel('Sp')
        ylabel = 'Distance from the ground state' if wrpt == 'grdstate' else 'Moves from the initial state'
        ax.set_ylabel(ylabel)

        grdstate = self._find_groundstate()

        colors = plt.colormaps['viridis'](np.linspace(
            0, 1, len(iterative_results_for_Tps.keys())))
        colors = colors[::-1]

        i = 0
        for tp, (fsample, rs) in iterative_results_for_Tps.items():
            Sp_list = list(rs.keys())
            target = grdstate.sample if wrpt == 'grdstate' else fsample
            distance_list = [hamming_distance(
                target, representatives[k].sample) for representatives in rs.values()]
            color = colors[i]
            ax.plot(Sp_list, distance_list, '.-',
                    label=f'Tp={tp}us', color=color)
            for Sp, dist in zip(Sp_list, distance_list):
                ax.annotate(dist, [Sp, dist])
            i += 1

    def _plot_BFS_based_iterative_results_per_initial_state(self, prev_optima, rounds, iterative_results_per_initial_state, to_dir=None):
        prefix = '(for all Tps)'

        (Ta, Tp, Sp), optima = prev_optima
        e = round(optima.energy, 3)

        initial_state, iterative_results_for_Tps = iterative_results_per_initial_state
        fe, moves = initial_state

        title = f'IRV{rounds}{prefix} when initialized in (moves={moves}|e={fe}) obtained by IRV{
            rounds-1}((Tp={Tp}us, Sp={Sp}, minima={e})|Ta={Ta}us)'
        fp = os.path.join(to_dir, f'{prefix}{title}.png')

        # if os.path.exists(fp):
        #     return

        grdstate = self._find_groundstate()
        cols = 3 if grdstate else 2

        fig, axs = plt.subplots(1, cols, figsize=(
            20, 6), constrained_layout=True)
        fig.suptitle(title)

        kwargs = {'prev_optima': prev_optima}
        for i, (ax, subject) in enumerate(zip(axs, ['energies', 'distances', 'distances'])):
            kwargs = copy.deepcopy(kwargs)
            if i == 2:
                kwargs['wrpt'] = 'grdstate'
            eval(f'self._plot_{subject}_against_Sp_for_all_Tps_with_fixed_initial_state')(
                ax, iterative_results_for_Tps, **kwargs)

        fe = round(fe, 3)
        ax = axs[0]
        ax.axhline(y=fe, color='red', linestyle='dotted')
        ax.annotate(fe, [ax.get_xlim()[-1], fe])
        if grdstate:
            for ax in axs[1:]:
                (fsample, _), *_ = iterative_results_for_Tps.values()
                abs_dist = hamming_distance(fsample, grdstate.sample)
                ax.axhline(y=abs_dist, color='red', linestyle='dotted')
                ax.annotate(abs_dist, [ax.get_xlim()[-1], abs_dist])

        plt.legend(bbox_to_anchor=(1.05, 1),
                   loc='upper left', borderaxespad=0.)
        plt.tight_layout()
        plt.savefig(fp)
        plt.close()

    def plot_BFS_based_iterative_results_for_initial_states(self, iterative_results_for_initial_states, *args, **kwargs):
        for iterative_results_per_initial_state in iterative_results_for_initial_states.items():
            self._plot_BFS_based_iterative_results_per_initial_state(
                *args, iterative_results_per_initial_state, **kwargs)
