This file explains the mathematical framework to simulate two-level scenarios currently implmented in asymmetree, and the extension we are implementing to include three-level scenarios with interactions.

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

In a first step, we simulate an holobiont scenario $(T_H, T_S, t', \sigma',\mu',\tau_H,\tau_S)$ using steps (1) and (2) of assymmetree, then we create an auxiliar tree $T_A$  as follows:

1. Take a copy of the subtree of $T_H$ rooted at $\rho_H$ and the tree $T_S$ (rooted at $0_S$).
2. Add the roots $0_A$ and $\rho_A$ with the corresponding edges $0_A\rho_A$, $\rho_A\rho_H$, and $\rho_A0_S$.
3. Construct the time map $\tau_A$ by coping $\tau_S$ and $\tau_H$ and setting
   - $\tau_A(\rho_A)=\tau_A(0_S)+\epsilon_0$, and
   - $\tau_A(0_A)=\tau_A(0_S)+\epsilon_0+\epsilon_1$ .
4. In a similar way, construct a map $\mu_A$, with $\mu_A(\rho_A)= \rho_A$, $\mu_A(0_A)= 0_A$, and $\mu_A(0_S)= \rho_A\rho_H$.
5. For the evolutionary events:
   - $0_S$ becomes a transfer event with a transfer edge as the only descendant.
   - $\rho_A$ is an speciation
   - $0_A$  is the planted root.

Now we will generate a gene tree $T_G$ evolving along $T_A$ considering two interactions: (I) metabolic redundancy when both guest and host carry the same gene, and (II) co-symbiosis, which is required to have horizontal gene transfer. To perform this simulation, we will generate a modification of step (2) of asymmetree.

So the simulation first Initializes $T_G$ as a single node $0_G$ with $\mu(0_G)=0_A$, and add it to the list of 'growing branches', then the the branches grow in time and number by sampling 'waiting times' and evolutionary events from a distribution defined by user-provided rates and the number of existing growing branches.

The improvement we purpose is to update the rates of the events considering intersections: in a given iteration of the simulation we have an incomplete $T_G$ where the non-loss leaves $L^0=\{x\in L(T_G) \text{ such that } t(x)\not=\times\}$ are the collection of growing branches.

First, given a growing gene branch $g$, let's define the interactors, divided in several classes:

| Class    | Definition                                                   | Description       |
| -------- | ------------------------------------------------------------ | ----------------- |
| **(A0)** | $\{g_0\in L^0 \text{ such that } \mu(g)=\mu(g_0) \and \mu(g)\preceq \rho_H \and g\not=g_0 \}$ | intra-host        |
| **(A1)** | $\{g_0\in L^0 \text{ such that } \mu(g)=\mu(g_0) \and \mu(g)\preceq \rho_S \and g\not=g_0 \}$ | intra-symbiont    |
| **(A2)** | $\{g_0\in L^0 \text{ such that } \mu(g)=\mu(g_0) \and \mu(g)\not\preceq \rho_H \and \mu(g)\not\preceq \rho_S \and g\not=g_0 \}$ | intra-independent |
| **(B0)** | $\{g_0\in L^0 \text{ such that } \mu(g)\not=\mu(g_0) \and \mu'(g_0) = \mu(g) \}$ | host-to-symbiont  |
| **(B1)** | $\{g_0\in L^0 \text{ such that } \mu(g)\not=\mu(g_0) \and \mu(g_0) = \mu'(g) \}$ | symbiont-to-host  |
| **(B2)** | $\{g_0\in L^0 \text{ such that } \mu(g)\not=\mu(g_0) \and \mu'(g_0) = \mu'(g) \}$ | inter-symbiont    |



Before providing a quick algorithm to obtain the interactors, we need some of definitions. Let $T_I$ a base tree and $T_J$ a tree that is being simulated inside $T_I$, where the non-loss leaves $L^0=\{x\in L(T_J) \text{ such that } t(x)\not=\times\}$ are the collection of growing branches. Asymmetree, as it is, constructs the map $\mu:V(T_J)\rightarrow V(T_I)\cup E(T_I)$ during the simulation of $T_J$.

Now, we will also construct the inverse map $\mu^{(-1)} : V(T_I) \times \mathbb{R} \rightarrow 2^{V(T_J)}$. Note that this map include real numbers in the domain; this allow us to map not only nodes in $T_I$, but also a specific point in the parent edges of the nodes: given a node $a\in V(T_I)$ with parent $b$, a scalar $t\in \mathbb{R} $, and node $x\in V(T_J)$, then $x \in \mu^{(-1)}(a,t)$ if one of the two following conditions holds:

1. $\mu(x)=a \in V(T_I) \text{ and } t=0$, or
2. $\mu(x)=ba \in E(T_I) \text{ and } \tau_I(a) < t < \tau_I(b) $.

<u>when punning loss events in $T_J$...</u>



 Let $\mu^{(-1)}$ and $\mu'^{(-1)}$ be the inverse maps of $\mu$ and $\mu'$, now, we can obtain the interactors as follows:

1. Given a branch $g\in L^0$, initialize the empty sets $A_0,A_1,A_2,B_0,B_1,B_2$.
2. Set $s\leftarrow \mu(g)$ and $h\leftarrow \mu'(s)$.
3. If $s\preceq \rho_H$
   1. $A_0\leftarrow\mu^{(-1)}(s)$
   2. $B\leftarrow \mu'^{(-1)}(s) $
   3. $B_0\leftarrow \bigcup_{x\in B} \mu^{(-1)}(x)$
4. elif $s\preceq 0_S$
   1. $A_1\leftarrow\mu^{(-1)}(s)$
   2. $B\leftarrow \mu'^{(-1)}(h) \setminus \{s\} $
   3. $B_1\leftarrow \mu^{(-1)}(h)$
   4. $B_2\leftarrow \bigcup_{x\in B} \mu^{(-1)}(x)$

5. else
   1. $A_2\leftarrow\mu^{(-1)}(s)$

6. Return $A_0,A_1,A_2,B_0,B_1,B_2$,

---



We will increase loss probability for genes residing in the same host; given two genes $x_0,x_1\in L^0$, we increase loss probability whenever

- **(A0)** $\mu(x_0)=\mu(x_1)$      # Both genes are in the same species
- **(A1)** $\mu'(x_0)=\mu'(x_1)$    # Genes are in different species, same host
- **(A2)** $\mu(x_0)= \mu'(x_1)$     # One gene in in host, the other is in a symbiont of such a host.

We should weight differentially those interactions.

For the case (A2) we can introduce an asymmetry by setting a category for gene trees: (S-keeps) where $x_0$ is lost with higher probability, or (H-keeps)  where $x_1$ is lost with higher probability

We will decrease transfer probability between species residing in different hots

Given a gene $x_0\in L^0$ together with the corresponding map $y_0=\mu(x_0)$, and a coexisting branch $y_1\in E(T_A)$, which is a possible recipient, we decrease transfer probability whenever:

- **(A3)** $\mu'(y_0)\not=\mu'(y_1)$
