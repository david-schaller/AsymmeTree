# -*- coding: utf-8 -*-

import os
import numpy as np

from asymmetree.tools import GraphTools
from asymmetree.file_io import ScenarioFileIO
import asymmetree.treeevolve as te

import asymmetree.best_matches as bm


class ComparisonAnalysis:
    
    def __init__(self, scenario, D, epsilon=0.000_000_01):
        
        self.scenario = scenario
        self.D = D
        self.epsilon=epsilon
        

    def compute_statistics(self):
        
        matrix_filename = 'temp.phylip'
        species_filename = 'temp_species.txt'
        tree_filename = 'temp_tree.txt'
        ScenarioFileIO.matrix_to_phylip(matrix_filename, self.scenario.genes,
                                        self.D)
        ScenarioFileIO.species_to_genes(species_filename, self.scenario)
        ScenarioFileIO.write_newick(tree_filename, self.scenario.S)
        
        # ---- true BMG, RBMG (based on l.c.a.) ----
        bmg, rbmg = self.scenario.bmg_subtrees, self.scenario.rbmg_subtrees
        self.possible_edges_bmg = self.scenario.possible_edges_bmg()
        
        # ---- epsilon method ----
        bmg_eps, rbmg_eps, self.time_eps = bm.ebh_qinfer(
                                              self.scenario, matrix_filename,
                                              species_filename, self.epsilon)
        bmg_eps, rbmg_eps = self.scenario.reduce_to_subtrees(bmg_eps, rbmg_eps)
        self.bmg_eps_stats = GraphTools.performance(bmg, bmg_eps)
        self.rbmg_eps_stats = GraphTools.performance(rbmg, rbmg_eps)
        
        # ---- neighbor-joining + midpoint rooting ----
        bmg_nj, rbmg_nj, self.time_nj = bm.bmg_by_tree_reconstruction(
                                            self.scenario, matrix_filename,
                                            supply_rbmg=True,
                                            return_calltime=True)
        bmg_nj, rbmg_nj = self.scenario.reduce_to_subtrees(bmg_nj, rbmg_nj)
        self.bmg_nj_stats = GraphTools.performance(bmg, bmg_nj)
        self.rbmg_nj_stats = GraphTools.performance(rbmg, rbmg_nj)
        
        # ---- quartet method (root subtree outgroups)----
        bmg_qd, rbmg_qd, self.time_qd = bm.quartet_qinfer(
                                            self.scenario, matrix_filename,
                                            species_filename, tree_filename)
        self.bmg_qd_stats = GraphTools.performance(bmg, bmg_qd)
        self.rbmg_qd_stats = GraphTools.performance(rbmg, rbmg_qd)
        
        # ---- quartet method (all outgroups)----
        bmg_qd2, rbmg_qd2, self.time_qd2 = bm.quartet_qinfer(
                                               self.scenario, matrix_filename,
                                               species_filename, tree_filename,
                                               closest_outgroups=True,
                                               incongruence_threshold=0.2)
        self.bmg_qd2_stats = GraphTools.performance(bmg, bmg_qd2)
        self.rbmg_qd2_stats = GraphTools.performance(rbmg, rbmg_qd2)
        
        os.remove(matrix_filename)
        os.remove(species_filename)
        os.remove(tree_filename)

    
    def summary_line(self):
        data = [*self.scenario.rates_and_counts(),
                self.epsilon,
                self.scenario.bmg_subtrees.size(),
                self.scenario.rbmg_subtrees.size(),
                self.possible_edges_bmg,
                *self.bmg_eps_stats, *self.rbmg_eps_stats,
                *self.bmg_nj_stats, *self.rbmg_nj_stats,
                *self.bmg_qd_stats, *self.rbmg_qd_stats,
                *self.bmg_qd2_stats, *self.rbmg_qd2_stats,
                self.time_eps, self.time_nj, self.time_qd, self.time_qd2]
        
        line = '{}' + (len(data)-1) * '\t{}'
        line = line.format(*data)
        return line
    
    
    def write_summary_line(self, filename, *first_columns):
        line = '\n' + (len(first_columns)) * '{}\t'
        line = line.format(*first_columns) + self.summary_line()
        with open(filename, 'a') as f:
            f.write(line)
        
    
    @staticmethod
    def summary_header():
        data = ['D_rate', 'L_rate', 'H_rate',
                'S_count', 'D_count', 'L_count', 'H_count', 'ancient_dupl',
                'epsilon',
                'BMG_size', 'RBMG_size',
                'possible_edges_BMG']
        for mode in ('eps', 'nj', 'qd', 'qd2'):
            for graph_type in ('BMG', 'RBMG'):
                for stat in ('order', 'size', 'tp', 'tn', 'fp', 'fn',
                             'acc', 'prec', 'recall'):
                    data.append(graph_type + "_" + mode + "_" + stat)
        data.extend(['time_eps', 'time_nj', 'time_qd', 'time_qd2'])
        
        header = '{}' + (len(data)-1) * '\t{}'
        header = header.format(*data)
        return header
    
    
    @staticmethod
    def write_header(filename, *first_columns):
        line = (len(first_columns)) * '{}\t'
        line = line.format(*first_columns) + ComparisonAnalysis.summary_header()
        with open(filename, 'w') as f:
            f.write(line)


# ---- simulation parameters ----

S_min, S_max = 3, 20            # range number of species
D_min, D_max = 0.5, 1.0         # range of dupl. rate
L_min, L_max = 0.5, 1.0         # range of loss rate
H_min, H_max = 0.0, 0.0         # range of HGT rate

repeats = 10
epsilon = 0.5

noise_sds = [0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5]
noise_alphas = []#[0.0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]

sim_ID = 0

stats_file = 'results/method_comparison.csv'

# write header of the main statistics file
ComparisonAnalysis.write_header(stats_file, 'SCEN_ID',
                                'noise_type', 'sd', 'true_sd', 'alpha')


for rep in range(repeats):
    
    D = np.random.uniform(low=D_min, high=D_max) if D_max > D_min else 0.0
    L = np.random.uniform(low=L_min, high=L_max) if L_max > L_min else 0.0
    H = np.random.uniform(low=H_min, high=H_max) if H_max > H_min else 0.0
    
    S = te.simulate_species_tree(np.random.randint(S_min, S_max+1), planted=True)
    
    TGT_simulator = te.GeneTreeSimulator(S)
    TGT1 = TGT_simulator.simulate(dupl_rate=D, loss_rate=L, hgt_rate=H)
    TGT1 = te.assign_rates(TGT1, S, base_rate=1,
                           autocorr_variance=0.2,
                           gamma_param=(0.5, 1.0, 2.2),
                           CSN_weights=(1/3, 1/3, 1/3))
    
    TGT2 = TGT_simulator.simulate(dupl_rate=D, loss_rate=L, hgt_rate=H)
    TGT2 = te.assign_rates(TGT2, S, base_rate=1,
                           autocorr_variance=0.2,
                           gamma_param=(0.5, 1.0, 2.2),
                           CSN_weights=(1/3, 1/3, 1/3))
    
    scenario1 = te.Scenario(S, TGT1, D, L, H)
    scenario2 = te.Scenario(S, TGT2, D, L, H)
    
    D1 = scenario1.get_distance_matrix()
    D2 = scenario2.get_distance_matrix()
    
    for noise_sd in noise_sds:
        D1_noisy = te.noisy_matrix(D1, noise_sd)
        D2_noisy = te.noisy_matrix(D2, noise_sd)
        true_sd1 = np.nanstd(np.divide(D1_noisy, D1)) if noise_sd > 0.0 else 0.0
        true_sd2 = np.nanstd(np.divide(D2_noisy, D2)) if noise_sd > 0.0 else 0.0
        
        summary1 = ComparisonAnalysis(scenario1, D1_noisy, epsilon=epsilon)
        summary1.compute_statistics()
        summary1.write_summary_line(stats_file, sim_ID,
                                    'no_bias', noise_sd, true_sd1, 0.0)
        
        summary2 = ComparisonAnalysis(scenario2, D2_noisy, epsilon=epsilon)
        summary2.compute_statistics()
        summary2.write_summary_line(stats_file, sim_ID+1,
                                    'no_bias', noise_sd, true_sd2, 0.0)
        
        print(sim_ID/ (2*repeats), 'Noise:', noise_sd)
    
    for noise_alpha in noise_alphas:
        D1_noisy, D2_noisy = te.convex_linear_comb(D1, D2, alpha=noise_alpha)
        true_sd1 = np.nanstd(np.divide(D1_noisy, D1)) if noise_alpha > 0.0 else 0.0
        true_sd2 = np.nanstd(np.divide(D2_noisy, D2)) if noise_alpha > 0.0 else 0.0
        
        summary1 = ComparisonAnalysis(scenario1, D1_noisy, epsilon=epsilon)
        summary1.compute_statistics()
        summary1.write_summary_line(stats_file, sim_ID,
                                    'bias', 0.0, true_sd1, noise_alpha)
        
        summary2 = ComparisonAnalysis(scenario2, D2_noisy, epsilon=epsilon)
        summary2.compute_statistics()
        summary2.write_summary_line(stats_file, sim_ID+1,
                                    'bias', 0.0, true_sd2, noise_alpha)
        
        print(sim_ID/ (2*repeats), 'alpha:', noise_alpha)
    
    sim_ID += 2
    
    print(sim_ID / (2*repeats))