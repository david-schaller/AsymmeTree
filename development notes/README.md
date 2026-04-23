# Generalizing AsymmeTree

Right now, AsymmeTree only allows to simulate two level evolution: species-gene reconciliation. In this project we will extend the framework to allow three level evolution: host-symbiont-gene reconciliation.

We will stick to the reconciliation simulation approach, but we add the host-symbiont reconciliation level. Furthermore, to model gene loss and gains, we will rely on game theory.

As starting point, we will simulate simple case: an epibiont evolving throughout genome reduction. We will consider only gene-to-gene interaction between homologs. and two types of interactions: redundancy and synergy.

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

[...]

Gene trees are simulated along a species tree using a constant-rate birth-death process [[27](https://www.mdpi.com/2674-113X/1/3/13#B27-software-01-00013)] where rates for the four types of events duplication, loss, HGT, and  gene conversion are user-defined. To this end, the user must specify  rates for the four event types which serve as parameters for exponential distributions from which waiting times until the next events are drawn. By setting the respective rate to zero, an event type is disabled  completely. The simulation starts with a single gene in the root of the species tree and proceeds stepwise in time towards the leaves by drawing waiting  times until the next event and lineages in which these events take place. At each point in the simulation, the total rate is given by the  sum of the four event types over the currently existing lineages in the  gene tree under construction. Hence, the simulation of all gene lineages progresses synchronously. This avoids difficulties arising by the fact  that some processes, such as replacing HGT as introduced below,  introduce dependencies between the gene lineages. In particular, with  this method, it is not necessary that simulated branches have to be  invalidated as a consequence of an event in a lineage that is processed  later in the simulation as it, e.g., occurs in [[26](https://www.mdpi.com/2674-113X/1/3/13#B26-software-01-00013)].

In addition to the randomly-occurring events, speciations and species  extinctions (which are determined by the species tree and therefore  fixed) are included as branching and loss events, respectively, and  affect all currently existing lineages in the respective species  lineage. In case of a speciation, each offspring species receives one  copy of each original gene lineage. On the other hand, the extinction of a species leads to loss of all of its gene lineages. If a waiting time  for a duplication, loss, or HGT event is drawn such that the next event  in the species tree occurs earlier, then this waiting time is discarded, the time is updated to the next species tree event, and the latter is  executed.

## Current structure of AsymmeTree code

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

Now we will generate a gene tree $T_G$ evolving along $T_A$ considering two interactions: (I) metabolic redundancy when both guest and host carry the same gene, and (II) co-symbiosis, which is required to have horizontal gene transfer. To perform this simulation, we will generate a modification of step (2) of asymmetree.

So the simulation first Initializes $T_G$ as a single node $0_G$ with $\mu(0_G)=0_A$, and add it to the list of 'growing branches', then the the branches grow in time and number by sampling 'waiting times' and evolutionary events from a distribution defined by user-provided rates and the number of existing growing branches.

The improvement we purpose is to update the rates of the events considering intersections: in a given iteration of the simulation we have an incomplete $T_G$ where the non-loss leaves $L^0=\{x\in L(T_G) \text{ such that } t(x)\not=\times\}$ are the collection of growing branches.

1. We will increase loss probability for genes residing in the same host; given two genes $x_0,x_1\in L^0$, we increase loss probability whenever

   - **(A0)** $\mu(x_0)=\mu(x_1)$      # Both genes are in the same species
   - **(A1)** $\mu'(x_0)=\mu'(x_1)$    # Genes are in different species, same host
   - **(A2)** $\mu(x_0)= \mu'(x_1)$     # One gene in in host, the other is in a symbiont of such a host.

   We should weight differentially those interactions.

   For the case (A2) we can introduce an asymmetry by setting a category for gene trees: (S-keeps) where $x_0$ is lost with higher probability, or (H-keeps)  where $x_1$ is lost with higher probability

2. We will decrease transfer probability between species residing in different hots

   Given a gene $x_0\in L^0$ together with the corresponding map $y_0=\mu(x_0)$, and a coexisting branch $y_1\in E(T_A)$, which is a possible recipient, we decrease transfer probability whenever:

   - **(A3)** $\mu'(y_0)\not=\mu'(y_1)$

# Implementation

The code should follow all AsymmeTree standards, in particular, have a look to `AsymmeTree/development notes/contributing.md`.

## Initialize `hologenome` module

- [x] Create a module `asymmetree.hologenome` inside of the directory `development notes/Example code/`
- [x] It should have the same structure as `asymmetree.genome`
- [x] It should be initialized with the code in `development notes/Example code/simulate_HS_level.py`. Until now, it simulates only the 'holobiont'.
- [x] Create the auxiliary tree $T_A$​ (One per each host-symbiont pair of trees).
- [x] Create an example script `development notes/Example code/example_simulations`. The result should be:
  - [x] One pandas series with Host trees trees (in general nhx format)
  - [x] One pandas dataframe with with symbiont trees  (in general nhx format)
- [x] Create a new copy of the whole package inside `development notes/Example code/` and add the new code.

### AI implementation notes

- The prototype now lives in `development notes/Example code/hologenome` with the same package
  pattern as `asymmetree.genome`: a thin `__init__.py` and a single simulation module.
- `HologenomeSimulator` currently covers the holobiont layer only and keeps the prototype logic
  close to `simulate_HS_level.py`, but removes the dependency on `revolutionhtl` by exporting
  trees with a local NHX serializer.
- The auxiliary tree uses collision-safe `H*`, `S*`, and `A*` labels and stores the host-side map
  in `reconc`, which gives us a clean hand-off for the later gene-level extension.
- `T_A` is built from the full symbiont tree rather than the pruned one so that later
  gene-level work still has access to extinct symbiont branches.
- The example script returns a host-tree `Series` and a symbiont-scenario `DataFrame`; the
  dataframe also includes the serialized auxiliary tree for each simulated host-symbiont pair.



## Add `holoevolve` module

- [x] Create a copy of `asymmetree.treeevolve` inside of the directory `development notes/Example code/asymmetree/`. Name this new module `holoevolve`
- [x] Add code to the `hologenome` module: simulate a single gene tree inside of the auxiliary tree $T_A$​, use the functions from `holoevolve` instead of `treeevolve`.
- [x] Make sure this new code is affectively called at the example script `development notes/Example code/example_simulations`. The result should be again a pandas dataframe.

### AI implementation notes

- The holobiont stage still uses `treeevolve`: host trees come from `treeevolve.species_tree_n()`
  and `HologenomeSimulator.simulate_symbiont_trees()` keeps using
  `treeevolve.dated_gene_tree()` / `treeevolve.prune_losses()`.
- `HologenomeSimulator.simulate_gene_trees()` simulates one gene tree per stored auxiliary tree and
  keeps both the unpruned and pruned versions alongside the host-symbiont results, using
  `holoevolve.dated_gene_tree()` / `holoevolve.prune_losses()`.
- `example_simulations.py` now runs the auxiliary-tree gene simulation for every scenario and adds
  `T_gene_unpruned` and `T_gene_pruned` to the returned scenario dataframe.
- The example currently reuses each scenario's DTL, replacement, and transfer-bias settings for
  the auxiliary-tree gene simulation so that every row remains self-contained.



## Add homolog-homolog interactions

We have to edit the `_get_branch_and_type` function. Right now it assumes that all the branches of the gene tree have exactly the same set of rates, which are the user-provided. The plan is to set a specific set of rates per branch. I'll deal with that later.

### AI implementation notes
