# Generalizing AsymmeTree

Right now, AsymmeTree only allows to simulate two level evolution: species-gene reconciliation. In this project we will extend the framework to allow three level evolution: host-symbiont-gene reconciliation.

We will stick to the reconciliation simulation approach, but we add the host-symbiont reconciliation level. Furthermore, to model gene loss and gains, we will rely on game theory.

As starting point, we will simulate simple case: an epibiont evolving throughout genome reduction. We will consider only gene-to-gene interaction between homologs. and two types of interactions: redundancy and synergy.

# Current structure of AsymmeTree

This report is inside of the repository of Asymmetree, in a development branch. The main directory is `../`, there we can find the [main readme](../README.md) which includes links for [wiki](https://github.com/david-schaller/AsymmeTree/wiki/Manual) and [documentation](https://david-schaller.github.io/docs/asymmetree/). The following explanation is taken from the wiki.

The module `asymmetree.genome.GenomeSimulation` provides functions that combine the simulation of phylogenetic trees and sequences. In particular, the class `GenomeSimulator` combines multiple steps in order to conveniently simulate whole genomes/proteomes. The (optional) output directory contains serialized trees, fasta files, and the true alignments. The gene trees and the sequences are simulated in subsequent steps using the classes' functions (i) `simulate_gene_trees(n, **kwargs)` and (ii) `simulate_sequences(subst_model, **kwargs)`.

The first step (i) takes the same keyword parameters as input as the function `gene_trees()` where `n` is the number of gene families to be simulated. Thus, rates for the three event types (`dupl_rate`, `loss_rate`, `hgt_rate`), autocorrelation (`autocorr_variance`), the distribution of base rates (`base_rate_distr`) etc. can be specified.

The second step (ii) simulates the sequences along the pruned part (without loss branches) of the simulated gene trees.

After step (i), the `list`s of full and pruned gene trees are accessible via the attributes `true_gene_trees` and `pruned_gene_trees`, respectively. Moreover, the full gene trees are serialized into the directory 'true_gene_trees' if an output directory was specified. After step (ii), the `list`s of sequence `dict`s are accessible via the attribute `sequence_dicts`.

Example usage:

```python
from asymmetree.treeevolve import species_tree_n_age
from asymmetree.genome import GenomeSimulator
from asymmetree.seqevolve import SubstModel, IndelModel

# simulate the common species tree
S = species_tree_n_age(10, 1.0, model='yule')

# specify models for sequence evolution
subst_model = SubstModel('a', 'JTT')
indel_model = IndelModel(0.01, 0.01, length_distr=('zipf', 1.821))

# initialy GenomeSimulator instance
gs = GenomeSimulator(S, outdir='simulation_directory')

# simulate 50 gene trees along the species tree S (and write them to file)
gs.simulate_gene_trees(50, dupl_rate=1.0, loss_rate=0.5,
                       base_rate=('gamma', 1.0, 1.0),
                       prohibit_extinction='per_species')

# simulate sequences along the gene trees
gs.simulate_sequences(subst_model,
                      indel_model=indel_model,
                      het_model=None,
                      length_distr=('constant', 200))

# results have been written to directories 'true_gene_trees',
# 'fasta_files' and 'alignments' in 'path/to/genome_directory'
```

# Graphs and phylogenies preliminaries

In this section, we introduce all the graph concepts that we will need to describe both current AsymmeTree workflow and the extension we will purpose.

A *graph* $G$ is an ordered pair $G = (V, E)$, where $V$ is the set of *nodes* and $E \subseteq \{uv \mid u, v \in V, u\neq v\}$ is the set
of *edges*, representing connections between nodes.  If $uv\in E$ implies $vu\in E$ for all $u,v\in V$, we call $G$ *undirected* and, otherwise, *directed*.

The *degree* of a node $v \in V$ in an undirected graph is the number of edges incident to $v$. In a digraph, the *in-degree* and *out-degree* of $v$ refer to the number of incoming edges $uv$ and outgoing edges $vu$, denoted $\deg^-_G(v)$ and $\deg^+_G(v)$, respectively. A *subgraph* $H$ of $G$ is a graph such that $V(H)
\subseteq V(G)$ and $E(H) \subseteq E(G)$. The *induced subgraph* of
$G$ on a node set $V' \subseteq V(G)$, denoted $G[V']$, is the subgraph $H
= (V', E')$, where $E' = \{uv \in E(G) \mid u, v \in V'\}$ contains all those edges from $E(G)$ that connect pairs of nodes in $V'$.

A *directed acyclic graph (DAG)* is a directed graph without directed cycles.  A *rooted tree* $T$ is a DAG that does not contain nodes $v$ with $\deg^-_G(v)>1$ and that has a unique *root* $\rho_T\in V(T)$, i.e., a
vertex $\rho_T$ with $\deg^-_G(v)=0$ and from which all nodes $v$ are reachable via a directed path from $\rho_T$ to $v$.  By definition, if $|V|>1$, then for each $v \in V \setminus \{\rho_T\}$ there is a unique incoming edge $uv$ and we call $u$ the *parent* of $v$ and $v$ a *child* of $u$.  We collect in $\ch_T(u)$ all childen of $u$ in $T$. The nodes of a rooted tree can be classified as *leaves*, $L(T)$, which are terminal nodes (nodes with no children), and *internal nodes*, $V^0(T)$, which have at least one child.  If there is a directed path from $u$ to $v$, then $u$ is an *ancestor* of $v$ and $v$ a *descendant* of $u$, denoted by $v\preceq_T u$.  We write $u \parallel
v$ if neither $v\preceq_T u$ nor $u\preceq_T v$ holds and say that $u$ and $v$ are *incomparable* in $T$.  Furthermore, the $\preceq_T$ relationship has been extended to consider edges within $T$ as follows. Let $e=uv$ be an edge and $x$ be a node in $T$. Then, we put $x \prec_T e$ if and only if $x\preceq_T v$. Moreover, we put $e \prec_T x$ if and only if $u\preceq_T
x$. For edges $e=(u,v)$ and $f=(a,b)$ in $T$ we put $e\preceq_T f$ if and only if $v \preceq_T b$.

For any node $v \in V(T)$, the *subtree rooted at $v$*, denoted $T(v)$, is the subgraph of $T$ induced by $v$ and all its descendants.  For any node $v \in V(T)$ in a rooted tree $T$, the *cluster* of $v$, denoted $C(v)$, is the set of leaves in the subtree rooted at $v$, i.e., $C(v) = L(T(v))$.  The *last common ancestor* (LCA) of a set $X
\subseteq V(T)$, denoted $\lca_T(X)$, is the $\prec_T$-minimal node in $T$ that is an ancestor of all nodes in $X$. If $X =\{x,y\}$ we write $\lca_T(x,y)$ instead of
$\lca_T(\{x,y\})$.

A *phylogenetic tree* is a rooted tree where each internal node has at least two children. Since gene duplications, losses or other events may predate the root of a phylogenetic tree (in particular, in the context of tree reconciliations), we define a *planted tree* that is formed by adding a new node $0_T$ to a phylogenetic tree $T$ and an edge $0_T \rho_T$. In a planted tree, the node $0_T$ has degree one, with its only neighbor being $\rho_T$, and this property remains unchanged during any subsequent modifications, such as resolving 'polytomies' as explained next. A *polytomy* in a phylogenetic tree is a node with more than two children. Resolving a polytomy involves converting it into a binary tree, a tree where each internal node has exactly two children. This may introduce additional internal nodes to maintain the tree's structure. Given a gene tree $T$ and a subset of leaves $L'\subseteq L(T)$, the *restriction* $T|_{L'}$ of $T$ to the set $L'$ is obtained from the minimal subtree of $T$ connecting all the leaves in $L'$ by suppressing all the vertices with a single child, with the possible exception of the root $\rho_T$ or planted root $0_T$​.

# Two-level evolutionary scenarios

In this section we describe the reconciliation formalism used in asymmetree. In what follows we consider particular types of phylogenetic trees, namely species trees and gene trees.  A *species tree* $T_S$ represents the evolutionary history of a set of species, where the leaves $L(T_S) \subseteq V(T_S)$ correspond to extant species. A *gene tree* $T_G$ represents the evolutionary history of a set of genes, where the leaves $L(T_G) \subseteq V(T_G)$ correspond to the sampled genes.

As genes exist within species, we denote by $\sigma: L(T_G) \to
L(T_S)$ the map that assigns each leaf in the gene tree to the species in which it resides.

Moreover, we can assign evolutionary events or mechanisms to the nodes of a gene tree $T_G$ that act on genes through evolution. Specifically, the map $t\colon V(T_G) \to \{\bullet, \square, \times, \star, \odot \}$ classifies nodes in the gene tree based on evolutionary events: $t(x) = \bullet$ for (co-)speciation, $t(x) = \square$ for duplication, $t(x) = \times$ for gene loss, $t(x) = \star$ for gene transfer, and $t(x) = \odot$ for gene existing genes. We can extend the domain of $t$ to also include vertices of $T_S$, where all inner nodes corresponds to speciations and leaves are ether existing extinct species.

An *relaxed evolutionary  scenario* $(T_S,T_G,t,\sigma,\mu)$ consists of a species tree $T_S$, a labeled gene tree $(T_G,t,\sigma)$ and a *reconciliation map* $\mu: V(T_G)\rightarrow V(T_S)\cup E(T_S)$ that satisfies:
$$
\begin{align*}
\text{(U1)} & \text{ \textit{Gene Constraint}. If } t(x)=\odot \text{, then }  \mu(x)=\sigma(x) \in L(T_S) \\
\text{(U2)} & \text{ \textit{Speciation Constraint}. If } t(x)=\bullet \text{, then } \mu(x) = \lca_S(\sigma(C_G(x))) \in V^0(T_S) \text{, and} \\&  \text{ } \mu(y_0) \parallel_S \mu(y_1) \text{ for any two distinct children } y_0, y_1 \in \ch_G(x) \\
\text{(U3)} & \text{ \textit{Duplication Constraint}. If } t(x)=\square \text{, then } \new{\mu(x)=e \text{ for some edge } e\in E(T_S)} \\
\text{(U4)} & \text{ \textit{Ancestor Constraint}. If } x\prec_G y \text{, then } \mu(x) \preceq_S \mu(y)
\end{align*}
$$
Up to this point, a relaxed evolutionary scenario allow us model the relationships between nodes and branches of gene and species trees. We are only missing time information, for which we introduce the time maps $\tau_S:V(S)\rightarrow\mathbb{R}$ and  $\tau_G:V(S)\rightarrow\mathbb{R}$, this map indicates how much time ago an evolutionary event occur. The time maps have the following restrictions, given a node $x\in V(T_G)$ and $y=\mu(x)$:

- $\tau_G(x)=0$​ if and only if $t(x)=\odot$​,
- $\tau_G(x)>0$ if and only if $t(x)\not=\odot$,
- $\tau_G(x) = \tau_S(y)$ if $y\in V(T_S)$,
- $\tau_S(b)\leq\tau_G(x)\leq\tau_S(a)$ if $y=ab\in E(T_S)$.

When we have the maps  $\tau_G$ and $\tau_S$, we say that the trees $T_G$ and $T_S$​ are dated.

This allow us to define an *evolutionary  scenario* as the tuple $(T_S,T_G,t,\sigma,\mu, \tau_S, \tau_G)$​, fully describing the trajectory of genes across species and time. Note that this is a two level scenario.

By default, the planted roots of the trees simulated by asymmetree are dated at time $1$.

# AsymmeTree framework

AssymmeTree is a powerful tool for simularion of two-level evolutionary scenarios. The following text taken from the [paper of asymmetree](https://www.mdpi.com/2674-113X/1/3/13) expalins their implementation. For us, the most intersting part are steps (1) and (2).

Like many other tools for phylogenetic simulation [[25](https://www.mdpi.com/2674-113X/1/3/13#B25-software-01-00013),[26](https://www.mdpi.com/2674-113X/1/3/13#B26-software-01-00013),[40](https://www.mdpi.com/2674-113X/1/3/13#B40-software-01-00013),[48](https://www.mdpi.com/2674-113X/1/3/13#B48-software-01-00013)], `AsymmeTree` generates complex gene family histories in a multilevel process that consists of the following five steps. **(1)** A dated species tree S is simulated whose leaves correspond to extant species (thus having  time stamp 0) and possibly (depending on simulation models and  parameters) extinction events. **(2)** Along this species tree, one or  multiple gene trees T are simulated using a variant of a constant-rate birth-death process [[27](https://www.mdpi.com/2674-113X/1/3/13#B27-software-01-00013)] that considers speciations, gene duplications, and (additive) HGTs as  “birth events” and gene losses as “death events”. In addition, replacing HGT and gene conversion lead to a bifurcation in one gene lineage  while, at the same time, another lineage gets lost. **(3)** In a next step,  evolution rate heterogeneity between the branches is introduced in a  manner that accounts for both species effects as well as asymmetric  evolution of paralogs after gene duplication. **(4)** The pruned gene tree  is then obtained by removing all branches that lead to loss events only  and suppressing all vertices that only have a single child left. In a final step **(5)**, nucleotide or amino acid  sequences can be evolved along the gene tree using a continuous-time  Markov chain, which is the standard model for this purpose [[31](https://www.mdpi.com/2674-113X/1/3/13#B31-software-01-00013),[32](https://www.mdpi.com/2674-113X/1/3/13#B32-software-01-00013)]. `AsymmeTree` supports a variety of models for substitution, indels, and among-site  heterogeneity.

# Three-level scenarios

Similarly to an evolutionary scenario, we define an *holobiont scenario* as the tuple $(T_H, T_S, t', \sigma',\mu',\tau_H,\tau_S)$​ describing the evolution of symbiont species across host species and time. In this case

- $\sigma':L(T_S)\cup L(T_H)\rightarrow L(T_H)$​,
- $\mu':V(T_S)\cup V(T_H)\rightarrow V(T_H)\cup E(T_H)$,
- $t'\colon V(T_S)\cup V(T_H) \to \{\bullet, \square, \times, \star, \odot \}$.

with the restriction that for $z\in V(T_H)$ we have $\mu'(z)=z$​ and $t'(z)\in\{\bullet,\times,\odot\}$.

If we merge an evolutionary scenario and an holobiont scenario, we get a *hologenome scenario* $(T_H, T_S,T_G,t,\sigma,\mu, t', \sigma',\mu', \tau_H, \tau_S, \tau_G)$​, note that this is a three-level scenario; genes evolving inside species, inside other species.

To keep a simple notation, let expand the domain of $\mu'$ to also consider vertices of $T_G$ as follows: given $x\in V(T_G)$ then $\mu'(x)=\mu'(\mu(x))$. Note that this is well defined only when $\mu(x)\in V(T_S)$.

# Simulation of hologenomes with interactions

To model gene coexistence and interactions effectively, we would need to combine reconciliation with Evolutionary Game Theory (EGT). A good part of the reconciliation approach is already implemented in AsymmeTree, and we will build on the top of that.

In a first step, we simulate an holobiont scenario $(T_H, T_S, t', \sigma',\mu',\tau_H,\tau_S)$ using steps (1) and (2) of assymmetree, then we create an auxiliar tree $T_A$  by taking a copy of the trees $T_H, T_S$ and merging the planted roots $0_H, 0_S$ into a new root $\rho_A$ and adding a new planted root $0_A$ with the corresponding planted edge $0_A\rho_A$.

$\tau_A$ is constructed by coping $\tau_S$ and $\tau_H$ and adding $\tau_A(0_A)=1+\epsilon$. Note that $\tau_S(0_S) = \tau_H(0_H) = 1 = \tau_A(\rho_S)$. In a similar way, we construct a map $\mu_A$, with $\mu_A(\rho_A)= \rho_A$ and $\mu_A(0_A)= 0_A$.

Now we can use step (2) of AsymmeTree to generate again...



1. Phylogenetic Reconciliation Approach

In this model, "coexistence" is often represented by the mapping of  multiple gene trees into a single symbiont or species tree to infer  shared history. 

- **TALE Framework**: TALE uses a Monte Carlo algorithm to sample evolutionary scenarios  across three nested levels (Host, Symbiont, and Gene). It accounts for **horizontal gene transfer (HGT)** and **host switches**, which are the biological precursors to coexistence.
- **DTL Model**: The core mathematical engine is the **Duplication-Transfer-Loss (DTL)** model. Coexistence is implicitly modeled when two gene lineages map to  the same species branch without a speciation event, indicating they  occupied that "genomic niche" simultaneously.
- **Simulation vs. Inference**: TALE is primarily an *inference* tool. To *simulate* based on coexistence, you might look into [ALE (Amalgamated Likelihood Estimation)](https://pubmed.ncbi.nlm.nih.gov/24033262/) or Zombi, which can simulate gene trees within species trees with specific transfer and loss rates. 

- 
- Game Theory for Gene Interactions

Game theory provides the mathematical structure to explain why certain genes "coexist" (cooperate) or why some are lost (conflict).

- **Replicator Dynamics**: Used to track the frequency of genes or microbes within the holobiont  population. If two genes provide a synergistic benefit (e.g.,  cross-feeding in a metabolic pathway), their interaction is a **Cooperative Game**.
- **Nash Equilibrium**: In a holobiont, a stable state of gene coexistence can be viewed as a  Nash Equilibrium where neither the host nor the microbe can improve  their fitness by unilaterally losing or acquiring a gene.
- **Conflict and Cooperation**: Interactions between host and symbiont genes can be modeled as a **Prisoner's Dilemma**. For instance, a symbiont may "cheat" by losing a gene that benefits the host but is costly to maintain. Game theory models predict the  conditions (like vertical transmission) under which cooperation (gene  retention) is the stable strategy. 

Comparison of Models

Component 

|                        | Reconciliation (e.g., TALE)          | Game Theory (EGT)                   |
| ---------------------- | ------------------------------------ | ----------------------------------- |
| **Primary Goal**       | Reconstruct history (DTL events)     | Predict stable strategies (ESS)     |
| **Mathematical Basis** | Maximum Parsimony / Likelihood       | Payoff Matrices / Replicator Eqs    |
| **Coexistence View**   | Physical presence in the same branch | Functional stability and synergy    |
| **Key Limitation**     | Doesn't model "why" a gene is kept   | Hard to map to specific phylogenies |

----

**Deeper description**

This is a sophisticated setup. You are essentially modeling 

**reconstructive evolution** where the payoff matrix dictates the "transition probabilities" for your reconciliation model.

In an epibiont (a guest living on the surface of a host), the proximity is high, but the guest usually maintains a degree of metabolic autonomy.  To simulate this using game theory, we can define the payoffs based on 

Metabolic Cost (

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

**)** and 

Fitness Benefit (

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

**)**.

1. The Payoff Matrix (Gene Coexistence Game)

We can model the interaction between a host gene (

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

) and a guest homolog (

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

) as a symmetric or asymmetric game.

|      | Guest Keeps (![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==) |
| ---- | ------------------------------------------------------------ |
|      |                                                              |

| )    | Guest Loses (![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==) |
| ---- | ------------------------------------------------------------ |
|      |                                                              |

| )                                                            |
| ------------------------------------------------------------ |
| **Host Keeps (![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)** |

| **)** | **Scenario A (Redundancy):** ![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==) |
| ----- | ------------------------------------------------------------ |
|       |                                                              |

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)
**Scenario B (Synergy):**

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

|      | **Host Solo:** ![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==) |
| ---- | ------------------------------------------------------------ |
|      |                                                              |

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

|                                                              |
| ------------------------------------------------------------ |
| **Host Loses (![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)** |

| **)** | **Guest Solo:** ![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==) |
| ----- | ------------------------------------------------------------ |
|       |                                                              |

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

|      | **Null:** ![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==) |
| ---- | ------------------------------------------------------------ |
|      |                                                              |

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

- **Scenario A (Redundancy):** If the benefit 

  ![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

 is fixed, keeping both genes is inefficient because of the double cost (

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

). This drives the system toward a **Nash Equilibrium** at 

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

 or 

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

, leading to genome reduction.

**Scenario B (Synergy):** A synergy term (

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

) is added. If 

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

, the equilibrium shifts to 

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

- , and both genes are conserved.
- Implementation in the Simulation

To integrate this into a reconciliation framework (like a TALE-like simulation), you can use the payoffs to calculate the **Transfer/Loss probabilities** at each node of your tree:

(i) Gene lost in Host / (ii) Gene lost in Guest

Instead of a fixed rate 

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

, your loss probability 

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==) becomes a function of the payoff:

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

In Scenario A, the payoff for "Keep" is low relative to the cost, making 

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==) high. Genome reduction in the epibiont is usually faster because the guest has a smaller effective population size (

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

), making it more susceptible to **Muller's Ratchet** (stochastic gene loss).

(iii) Gene Transfer (HGT)

Transfer can be modeled as a move to reach a higher payoff state. If the Host loses a vital gene (

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

), the pressure for a Guest-to-Host transfer increases to restore the 

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

 benefit to the holobiont unit.

3. Suggested Workflow for Your Simulation

Since TALE is an inference tool, you might need a custom script to *generate* the data.

1. **Initialize:** Create a Host tree and a Guest (epibiont) tree with an initial set of homolog genes.

2. **Iterate:** At each internal branch, calculate the "Holobiont Fitness" based on your matrix.

3. **Decide:**

   - If **Scenario A**: Randomly trigger a Loss event in either tree (weighted by 

     ![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

).

If **Scenario B**: Apply a "Conservation Bonus" that lowers the Loss probability to near zero.

If **Transfer**: Use a **Poisson process** where the rate 

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

1. -  is scaled by the fitness deficit of the receiving partner.

**A quick tip for the "Naive" Matrix:**
Epibionts often lose genes related to **biosynthesis** (amino acids) but keep genes for **attachment and defense**. You might want to assign different 

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

 and 

![img](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

 values based on the *category* of the gene to see if your simulation produces realistic genome reduction patterns.

# Implementation

## Standards to follow
