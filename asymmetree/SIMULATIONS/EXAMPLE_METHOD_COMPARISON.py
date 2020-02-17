# -*- coding: utf-8 -*-

import time, os
import numpy as np

from asymmetree.tools import FileIO, GraphTools
import asymmetree.simulator.TreeSimulator as ts
import asymmetree.simulator.TreeImbalancer as tm
from asymmetree.simulator.Scenario import Scenario
from asymmetree.simulator import NoisyMatrix

from asymmetree.best_matches import TrueBMG
import asymmetree.best_matches.ExtBestHits as ebh
import asymmetree.best_matches.TreeReconstruction as tr
import asymmetree.best_matches.Quartets as qd


class CustomAnalysis:
    
    def __init__(self, scenario, D, epsilon=0.000_000_01):
        
        self.scenario = scenario
        self.D = D
        self.epsilon=epsilon
        

    def compute_statistics(self):
        
        matrix_filename = "temp.phylip"
        species_filename = "temp_species.txt"
        tree_filename = "temp_tree.txt"
        FileIO.matrix_to_phylip(matrix_filename, self.scenario.genes, self.D)
        FileIO.species_to_genes(species_filename, self.scenario)
        FileIO.write_newick(tree_filename, self.scenario.S)
        
        # ---- true BMG, RBMG (based on l.c.a.) ----
        BMG, RBMG = self.scenario.BMG_subtrees, self.scenario.RBMG_subtrees
        self.possible_edges_BMG = self.scenario.possible_edges_BMG()
        
        # ---- epsilon method ----
        start = time.time()
        BMG_eps, RBMG_eps, self.time_eps = ebh.ebh_qinfer(self.scenario, matrix_filename,
                                                          species_filename, self.epsilon)
        BMG_eps, RBMG_eps = self.scenario.reduce_to_subtrees(BMG_eps, RBMG_eps)
        self.BMG_eps_stats = GraphTools.performance(BMG, BMG_eps)
        self.RBMG_eps_stats = GraphTools.performance(RBMG, RBMG_eps)
        
        # ---- neighbor-joining + midpoint rooting ----
        start = time.time()
        nj_tree = tr.neighbor_joining(self.scenario.genes, self.scenario.gene_index, matrix_filename)
        tr.midpoint_rooting(nj_tree)
        BMG_nj, RBMG_nj = TrueBMG.best_match_graphs(nj_tree)
        self.time_nj = time.time() - start
        BMG_nj, RBMG_nj = self.scenario.reduce_to_subtrees(BMG_nj, RBMG_nj)
        self.BMG_nj_stats = GraphTools.performance(BMG, BMG_nj)
        self.RBMG_nj_stats = GraphTools.performance(RBMG, RBMG_nj)
        
        # ---- quadruple method (root subtree outgroups)----
        BMG_qd, RBMG_qd, self.time_qd = qd.quartet_qinfer(self.scenario, matrix_filename,
                                                          species_filename, tree_filename)
        self.BMG_qd_stats = GraphTools.performance(BMG, BMG_qd)
        self.RBMG_qd_stats = GraphTools.performance(RBMG, RBMG_qd)
        
        # ---- quadruple method (all outgroups)----
        BMG_qd2, RBMG_qd2, self.time_qd2 = qd.quartet_qinfer(self.scenario, matrix_filename,
                                                             species_filename, tree_filename,
                                                             closest_outgroups=True,
                                                             incongruence_threshold=0.2)
        self.BMG_qd2_stats = GraphTools.performance(BMG, BMG_qd2)
        self.RBMG_qd2_stats = GraphTools.performance(RBMG, RBMG_qd2)
        
        os.remove(matrix_filename)
        os.remove(species_filename)
        os.remove(tree_filename)

    
    def summary_line(self):
        data = [*self.scenario.rates_and_counts(),
                self.epsilon,
                self.scenario.BMG_subtrees.size(), self.scenario.RBMG_subtrees.size(),
                self.possible_edges_BMG,
                *self.BMG_eps_stats, *self.RBMG_eps_stats,
                *self.BMG_nj_stats, *self.RBMG_nj_stats,
                *self.BMG_qd_stats, *self.RBMG_qd_stats,
                *self.BMG_qd2_stats, *self.RBMG_qd2_stats,
                self.time_eps, self.time_nj, self.time_qd, self.time_qd2]
        
        line = "{}" + (len(data)-1) * "\t{}"
        line = line.format(*data)
        return line
    
    
    def write_summary_line(self, filename, *first_columns):
        line = "\n" + (len(first_columns)) * "{}\t"
        line = line.format(*first_columns) + self.summary_line()
        with open(filename, "a") as f:
            f.write(line)
        
    
    @staticmethod
    def summary_header():
        data = ["D_rate", "L_rate", "H_rate",
                "S_count", "D_count", "L_count", "H_count", "ancient_dupl",
                "epsilon",
                "BMG_size", "RBMG_size",
                "possible_edges_BMG"]
        for mode in ("eps", "nj", "qd", "qd2"):
            for graph_type in ("BMG", "RBMG"):
                for stat in ("order", "size", "tp", "tn", "fp", "fn", "acc", "prec", "recall"):
                    data.append(graph_type + "_" + mode + "_" + stat)
        data.extend(["time_eps", "time_nj", "time_qd", "time_qd2"])
        
        header = "{}" + (len(data)-1) * "\t{}"
        header = header.format(*data)
        return header
    
    
    @staticmethod
    def write_header(filename, *first_columns):
        line = (len(first_columns)) * "{}\t"
        line = line.format(*first_columns) + CustomAnalysis.summary_header()
        with open(filename, "w") as f:
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

stats_file = "results/method_comparison.csv"

# write header of the main statistics file
CustomAnalysis.write_header(stats_file, "SCEN_ID", "noise_type", "sd", "true_sd", "alpha")


for rep in range(repeats):
    
    D = np.random.uniform(low=D_min, high=D_max) if D_max > D_min else 0.0
    L = np.random.uniform(low=L_min, high=L_max) if L_max > L_min else 0.0
    H = np.random.uniform(low=H_min, high=H_max) if H_max > H_min else 0.0
    
    S = ts.build_species_tree(np.random.randint(S_min, S_max+1), planted=True)
    
    TGT1 = ts.build_gene_tree(S, (D,L,H))
    TGT1 = tm.imbalance_tree(TGT1, S, baseline_rate=1,
                            lognormal_v=0.2,
                            gamma_param=(0.5, 1.0, 2.2),
                            weights=(1/3, 1/3, 1/3),
                            copy_tree=False)
    
    TGT2 = ts.build_gene_tree(S, (D,L,H))
    TGT2 = tm.imbalance_tree(TGT2, S, baseline_rate=1,
                            lognormal_v=0.2,
                            gamma_param=(0.5, 1.0, 2.2),
                            weights=(1/3, 1/3, 1/3),
                            copy_tree=False)
    
    scenario1 = Scenario(S, TGT1, (D,L,H))
    scenario2 = Scenario(S, TGT2, (D,L,H))
    
    D1 = scenario1.get_distance_matrix()
    D2 = scenario2.get_distance_matrix()
    
    for noise_sd in noise_sds:
        D1_noisy = NoisyMatrix.noisy_matrix(D1, noise_sd)
        D2_noisy = NoisyMatrix.noisy_matrix(D2, noise_sd)
        true_sd1 = np.nanstd(np.divide(D1_noisy, D1)) if noise_sd > 0.0 else 0.0
        true_sd2 = np.nanstd(np.divide(D2_noisy, D2)) if noise_sd > 0.0 else 0.0
        
        summary1 = CustomAnalysis(scenario1, D1_noisy, epsilon=epsilon)
        summary1.compute_statistics()
        summary1.write_summary_line(stats_file, sim_ID, "no_bias", noise_sd, true_sd1, 0.0)
        
        summary2 = CustomAnalysis(scenario2, D2_noisy, epsilon=epsilon)
        summary2.compute_statistics()
        summary2.write_summary_line(stats_file, sim_ID+1, "no_bias", noise_sd, true_sd2, 0.0)
        
        print(sim_ID/ (2*repeats), "Noise:", noise_sd)
    
    for noise_alpha in noise_alphas:
        D1_noisy, D2_noisy = NoisyMatrix.convex_linear_comb(D1, D2, alpha=noise_alpha)
        true_sd1 = np.nanstd(np.divide(D1_noisy, D1)) if noise_alpha > 0.0 else 0.0
        true_sd2 = np.nanstd(np.divide(D2_noisy, D2)) if noise_alpha > 0.0 else 0.0
        
        summary1 = CustomAnalysis(scenario1, D1_noisy, epsilon=epsilon)
        summary1.compute_statistics()
        summary1.write_summary_line(stats_file, sim_ID, "bias", 0.0, true_sd1, noise_alpha)
        
        summary2 = CustomAnalysis(scenario2, D2_noisy, epsilon=epsilon)
        summary2.compute_statistics()
        summary2.write_summary_line(stats_file, sim_ID+1, "bias", 0.0, true_sd2, noise_alpha)
        
        print(sim_ID/ (2*repeats), "alpha:", noise_alpha)
    
    sim_ID += 2
    
    print(sim_ID / (2*repeats))