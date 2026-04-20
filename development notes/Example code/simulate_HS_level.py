"""
This script simulates Symbiont trees inside of Host trees.

Parameters:
- host_n_species: a list of integers; the number of species in host trees.
                  The length of the list (N_h) determines the number of host trees to simulate.

- symbiont_DTL_rates: a list of list of floats; ech sublist contains duplication/loss/transfer rates for symbiont inside of a host.
                   Each of the N_r sublist will be used for all of the host trees.

- n_simulations: an integer (N_s); the number of independent simulations for a given set of parameters (host tree size and DTL rates).

- transfer_distance_bias: List; each of the N_db elements specifies whether closer related species have a higher probability to be the recipient species in a HGT event. (see https://david-schaller.de/docs/asymmetree/treeevolve/GeneTree.html)

- replace_prob: List; each of the N_rp elements indicates probability for replacing HGT events. (see https://david-schaller.de/docs/asymmetree/treeevolve/GeneTree.html)

In total, this script will generate N_h * N_r * N_s * N_db * N_rp evolutionary scenarios.
Note, there will be N_h unique host trees, each containing N_r * N_s * N_db * N_rp guest trees.
"""

from itertools import product
from pathlib import Path

from pandas import DataFrame, concat, Series
import random as py_random
import numpy.random as np_random
from tqdm import tqdm
tqdm.pandas()

from revolutionhtl.nhx_tools import get_nhx, read_nhx
from asymmetree.treeevolve import dated_gene_tree, prune_losses, species_tree_n


#########################
# Simulation parameters #
#########################
host_n_species= [5, 10, 30]

symbiont_DTL_rates= [[0.133, 0.266, 0.266],
                  [0.3, 0.6, 0.6],
                  [0.6, 1.2, 1.2]]

replace_prob= [0, 0.5, 1]

transfer_distance_bias= [False, 'exponential']

n_simulations= 5

output_dir= Path('../data/three_level_scenarios/')

#############
# Constants #
#############
seed= 13042026
R= lambda x: int(100*x)

total_simulations= len(host_n_species) * len(symbiont_DTL_rates) * len(replace_prob) * n_simulations * len(transfer_distance_bias)

#############
# Functions #
#############
def simulate_HG_scenario(speciens_tree, D, T, L, R_probability, T_bias):

    gene_tree = dated_gene_tree(speciens_tree,
                                dupl_rate= D,
                                loss_rate= L,
                                hgt_rate= T,
                                prohibit_extinction="per_family",
                                replace_prob= R_probability,
                                transfer_distance_bias=T_bias,
                                )

    prunned_gene_tree= prune_losses(gene_tree)

    return get_tralda_nhx(gene_tree), get_tralda_nhx(prunned_gene_tree)

def get_tralda_nhx(T):
    T, root= T.to_nx()
    return get_nhx(T, root, name_attr='label')

def get_str_vector(*X):
    return '_'.join(map(str, X))


###################
# Run simulations #
###################
py_random.seed(seed)
np_random.seed(seed)

df_scenarios= []
D_h_trees= dict()
D_repetitions= dict()

with tqdm(total=total_simulations, desc='Running simulations') as pbar:
    for n_Host in host_n_species:

        D_repetitions[n_Host]= D_repetitions.get(n_Host, 0) + 1
        hostTree_id= f'N{n_Host}.{D_repetitions[n_Host]}'

        D_h_trees[hostTree_id]= species_tree_n(n_Host, innovation= True)

        for (D, T, L), R_probability, T_bias, i in product(symbiont_DTL_rates,
                                                           replace_prob,
                                                           transfer_distance_bias,
                                                           range(n_simulations)):

            Tg_full, Tg_prunned= simulate_HG_scenario(D_h_trees[hostTree_id], D, T, L, R_probability, T_bias)
            parameter_vector= get_str_vector(hostTree_id, R(D), R(T), R(L), R(R_probability), T_bias, i)
            df_scenarios+= [[parameter_vector, Tg_full, Tg_prunned]]
            pbar.update(1)

# parameter vector fields: 'n_hosts', 'Dr', 'Tr', 'Lr', 'replacement_p', 'trasfer_distance_bias', 'idx'

# Get dataframe
#--------------
dfCols= ['parameters_vector', 'T_symbiont_unprunned', 'T_symbiont_prunned']
df_scenarios= DataFrame(df_scenarios, columns= dfCols).set_index('parameters_vector')

df_scenarios

D_h_trees= Series({id:get_tralda_nhx(T) for id,T in D_h_trees.items()})
D_h_trees

#############################
# Formatting recnciliations #
#############################
def improve_host_labels(T, preffix):
    T.old2new_label= dict()
    for X in T:
        oldlabel= T.nodes[X]['label']
        newlabel= f'{preffix}{oldlabel}'
        T.nodes[X]['label']= newlabel
        T.nodes[X]['reconc']= newlabel

        T.old2new_label[oldlabel]= newlabel

    return None

def improve_guest_labels(T_g, T_h, g_preffix):
    T_g.old2new_label= dict()
    for X_g in T_g:
        g_oldlabel= T_g.nodes[X_g]['label']
        g_map= T_g.nodes[X_g]['reconc']

        if g_map is not None:

            if g_map.startswith('('):
                g_map= g_map.split(', ')[1][:-1]
            h_label= T_h.old2new_label[g_map]
            g_newlabel= f'{g_preffix}{g_oldlabel}_{h_label}'

            T_g.nodes[X_g]['label']= g_newlabel
            T_g.nodes[X_g]['reconc']= h_label
            T_g.old2new_label[g_oldlabel]= g_newlabel

    return None

def add_reconciliation_label(T):
    for X in T:
        event= T.nodes[X]['event']
        g_label= T.nodes[X]['label']
        s_label= T.nodes[X]['reconc']
        dist= T.nodes[X]['dist']

        reclabel= f'{event}|{g_label}->{s_label}'
        reclabel_dated= f'{reclabel}:{dist}'

        T.nodes[X]['reclabel']= reclabel
        T.nodes[X]['reclabel_dated']= reclabel_dated

ignoreattrs= ['label', 'dist', 'transferred', 'reconc',
              'tstamp', 'event', 'sibling_nr', 'reclabel',
              'reclabel_dated']

def improve_row_labels(row, D_aux):

    hostTree_id= row.name.split('_')[0]
    T_host= D_aux[hostTree_id]

    improve_guest_labels(row.T_symbiont_unprunned, T_host, 'G')
    improve_guest_labels(row.T_symbiont_prunned, T_host, 'G')

    nhx_g_u= get_nhx(row.T_symbiont_unprunned, row.T_symbiont_unprunned.root, name_attr='label')
    nhx_g_p= get_nhx(row.T_symbiont_prunned, row.T_symbiont_unprunned.root, name_attr='label')

    add_reconciliation_label(row.T_symbiont_unprunned)
    add_reconciliation_label(row.T_symbiont_prunned)

    newick_g_u= get_nhx(row.T_symbiont_unprunned,
                        root= row.T_symbiont_unprunned.root,
                        name_attr= 'reclabel',
                        ignore_attrs= ignoreattrs)

    newick_g_p= get_nhx(row.T_symbiont_prunned,
                        root= row.T_symbiont_prunned.root,
                        name_attr= 'reclabel',
                        ignore_attrs= ignoreattrs)

    newick_g_u_dated= get_nhx(row.T_symbiont_unprunned,
                        root= row.T_symbiont_unprunned.root,
                        name_attr= 'reclabel_dated',
                        ignore_attrs= ignoreattrs)

    newick_g_p_dated= get_nhx(row.T_symbiont_prunned,
                        root= row.T_symbiont_prunned.root,
                        name_attr= 'reclabel_dated',
                        ignore_attrs= ignoreattrs)

    return Series(dict(T_symbiont_unprunned= nhx_g_u,
                       T_symbiont_prunned= nhx_g_p,
                       T_symbiont_unprunned_simple= newick_g_u,
                       T_symbiont_prunned_simple= newick_g_p,
                       T_symbiont_unprunned_simple_dated= newick_g_u_dated,
                       T_symbiont_prunned_simple_dated= newick_g_p_dated,
                       ))

# Load as nx trees
D_aux= D_h_trees.map(read_nhx)
df_aux= df_scenarios.map(read_nhx)

# Change gene and species names and simplify reconciliation format
D_aux.map(lambda T: improve_host_labels(T, 'H'))
df_aux= df_aux.apply(lambda row: improve_row_labels(row, D_aux), axis= 1)
df_aux

#####################
# Write simulations #
#####################

# original nhx format
ofile= output_dir / 'holobiont_scenarios.csv'
df_aux[['T_symbiont_unprunned', 'T_symbiont_prunned']].to_csv(ofile)
print(f'Written {ofile}')

# simplified dated format
ofile= output_dir / 'holobiont_scenarios_simple_dated.csv'
df_aux[['T_symbiont_unprunned_simple', 'T_symbiont_prunned_simple']].to_csv(ofile)
print(f'Written {ofile}')

# simplified, non dated format
ofile= output_dir / 'holobiont_scenarios_simple.csv'
df_aux[['T_symbiont_unprunned_simple_dated', 'T_symbiont_prunned_simple_dated']].to_csv(ofile)
print(f'Written {ofile}')
