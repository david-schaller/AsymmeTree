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

For technical reasons, we assume that a host species is evolving inside of itself, which is ensured by the following restriction: for $z\in V(T_H)$ we have $\mu'(z)=z$ and $t'(z)\in\{\bullet,\times,\odot\}$.

If we merge an evolutionary scenario and an holobiont scenario, we get a *hologenome scenario* $(T_H, T_S,T_G,t,\sigma,\mu, t', \sigma',\mu', \tau_H, \tau_S, \tau_G)$​​, note that this is a three-level scenario; genes evolving inside species, inside other species.

we allow a slight abuse of notation: for and edge $uv\in E(T_S)\cup E(T_H)$, we set $\mu'(uv)= \mu'(v)$​.

Note that for a gene $g$, it resides in a symbiont if and only if $\mu(g)\not=\mu'(\mu(g))$, otherwise, $g$ resides directly into a host species.

# Simulation of hologenomes with interactions

## The big picture

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

So the simulation first Initializes $T_G$ as a single node $0_G$ with $\mu(0_G)=0_A$, then a 'growing branch' is created as descendant of of this node, and this growing branch is associated with the edge(s) descending from $0_A$. Furthermore, we initialize the remaining time $t=\tau_G(0_G)=\tau_A(0_A)$.  (Stress that the simulation is top-down, so time will decrease).

At each iteration of the simulation,The simulation proceeds by randomly obtaining a new time $t'< t$ at which an evolutionary event $\epsilon$ happens, affecting a set $\Gamma$ of growing branches; If $\epsilon$ is an speciation corresponding to a node $v\in V(T_A)$ with parent $u$, then $\Gamma$ contains every gene growing branches $g$ associated to the species branch $uv\in E(T_A)$. Otherwise, $\Gamma$ contains only one gene growing branch $g$ associated to a node $v\in V(T_A)$ with parent $u$​.

The evolutionary event $\epsilon$ is chosen based on user provided rates.

The evolutionary event $\epsilon$ implies the creation of a new gene node $x$ for each growing branch $g$ of $\Gamma$, thus $g$ becomes into an edge, and is not a growing branch anymore. Additionally, if $\epsilon$ is not a loss event, then new growing branches will be created below $x$.

Finally, update the remaining time: $t\leftarrow t'$ and iterate again until $t=0$​​​​.

### Note on branches position

Since **we are working on unpruned trees**, if $s$ is a symbiont branch, it is fully contained in a host branch. This is in general not true when we work with pruned symbiont trees:

Given a species bifurcation formed by the three edges $AB,BC,BD$. If another branch is evolving inside $AB$ and going through this speciation, it will include a bifurcation formed by the edges $ab,bc,bd$. If we prune $bc$, then the other two edges $ab$ and $bd$ will be merged into a new edge $ad$. By construction, the original edges $ab$ and $bd$ where fully contained in $AB$ and $BD$ correspondingly, nevertheless, the new edge $ad$ passes through more than one species edge.

## Map of growing branches and inverse

When simulating gene tree $T_G$ inside auxiliary tree $T_A$, in a given iteration of the simulation we have an incomplete $T_G$ with a set of growing branches $\Pi$. To track the relation between growing branches and base tree, we define the temporary **growing map** $\kappa : \Pi \rightarrow E(T_A)$. Each growing branch $g\in \Pi$ is associated with one species branch $s=uv\in E(T_A)$, so we write $\mu(g)=v$. Note that in the code $\kappa$ is implicitly stored as an attribute of the growing branches.

Recall, after simulation of symbiont tree $T_S$ inside of host tree $T_H$, we end up with the reconciliation map $\mu':V(T_S)\cup V(T_H)\rightarrow V(T_H)\cup E(T_H)$. Afterwards, when creating the auxiliary tree $T_A$, the same reconciliation map is inherited as $\mu' : V(T_A)\rightarrow V(T_A)\cup E(T_A)$​.

Note that for gene growing branch $g$, it resides in a symbiont if and only if $\kappa(g)\not=\mu'(\kappa(g))$, otherwise, $g$ resides directly into a host species.

During the simulation of $T_S$ we have to construct and return an **inverse map** $\gamma' : E(T_H) \rightarrow 2^{E(T_S)}$ of $\mu'$, where $\gamma'(ab)=\{ uv\in E(T_S) \mid b \preceq_H \mu'(v) \prec_H \mu'(u) \preceq_H a \}$ for $ab\in E(T_H)$​​.

Furthermore, we can **time-restrict the inverse map** to return only branches *touching* a time $t$ as follows:

$\gamma'_t(ab) = \{ uv \in \gamma(b) \mid \tau_S(v) \leq t \le \tau_S(u) \}$ for $ab\in E(T_H)$.

 In a similar way, we can define the inverse map $\gamma : E(T_A) \rightarrow 2^\Pi$ of $\kappa$ as $\gamma(ab)=\{ g\in \Pi \mid \kappa(g)=ab \}$ for $ab\in E(T_A)$.

We will use the maps $\kappa,\gamma,\text{ and } \gamma'$ during the simulation of a gene tree inside of the auxiliary tree to easily track gene-to-gene interactions.

## Track gene-gene interactions

> [!IMPORTANT]
>
> The six categories for interactions bellow are not being used in the implementation anymore. These interaction classes were useful for exploring the structure of overlaps, but they are no longer the right primitive for loss-rate updates.

Given a growing gene branch $g$, let $s=uv = \kappa(g)$ be the species where the gene is evolving, and $h=\mu'(s)$ be the host of the symbiont.

Recall that we assume that the hist branch is contained inside itself, so if $g$ is a gene evolving in the host, then $s=h\preceq \rho_H$, otherwise $s\preceq \rho_S \parallel h \preceq \rho_H$.

Let's define the interactors of $g$, divided in several classes:

| Class    | Definition                                                   | Description       |
| -------- | ------------------------------------------------------------ | ----------------- |
| **(A0)** | $\{g_0\in \Pi\setminus\{g\} \mid \kappa(g_0)=s \preceq \rho_H  \}$ | intra-host        |
| **(A1)** | $\{g_0\in \Pi\setminus\{g\} \mid \kappa(g_0)=s \preceq \rho_S  \}$ | intra-symbiont    |
| **(A2)** | $\{g_0\in \Pi\setminus\{g\} \mid \kappa(g_0)=s \and s\not\preceq \rho_H \and s\not\preceq \rho_S  \}$ | intra-independent |
| **(B0)** | $\{g_0\in \Pi \setminus\{g\} \mid \kappa(g_0)\not=  \mu'(\kappa(g_0)) = s \preceq \rho_H \}$ | host-to-symbiont  |
| **(B1)** | $\{g_0\in \Pi \setminus\{g\} \mid \kappa(g_0) = h \not= s \}$ | symbiont-to-host  |
| **(B2)** | $\{g_0\in \Pi \setminus\{g\} \mid s \not= \kappa(g_0)\not= \mu'(\kappa(g_0)) = h\not=s \}$ | inter-symbiont    |

The interactors can be quickly computed as follows:

1. Given a growing branch $g\in \Pi$, initialize the empty sets $A_0,A_1,A_2,B_0,B_1,B_2$.
2. $s\leftarrow \kappa(g)$ is the symbiont branch where $g$ is growing. Note $s=uv\in E(T_A)$​.
3. $h\leftarrow \mu'(s)$ is the species branch where $s$ is evolving.
4. If $s = h$
   1. $A_0\leftarrow\gamma(s)\setminus\{g\}$
   2. $B\leftarrow \gamma'(s) $
   3. $B_0\leftarrow \bigcup_{x\in B} \gamma(x)$
5. elif $s\preceq 0_S$
   1. $A_1\leftarrow\gamma(s)$
   2. $B\leftarrow \gamma'(h) \setminus \{s\} $
   3. $B_1\leftarrow \gamma(h)$
   4. $B_2\leftarrow \bigcup_{x\in B} \gamma(x)$

6. else
   1. $A_2\leftarrow\gamma(s)$
7. Return $A_0,A_1,A_2,B_0,B_1,B_2$,

## Update rates based on gene-gene interactions

The simulation starts with a base loss rate $r_l$. Whenever an holobiont system starts existing, the loss rate is updated for each of the involved species. Below we provide an expression to update the rates, this should be updated whenever gene content changes in a holobiont system, it my be due to duplication/loss/transfer of both genes and symbionts.

Let's say we have a host branch $h$ with $|\gamma'(h)|+1= N$, i.e. we have a total of $N$ species in the system, with a total of $M=|\gamma(h)|+\sum_{s\in\gamma'(h)}|\gamma(s)|$ genes. The loss rate for a gene existing in a symbiont $s$ will be:
$$
r_l^*= r_l + \alpha \frac{|\gamma(s)|}{M/N} = r_l + \alpha \frac{|\gamma(s)| N}{M}
$$
On the other hand, the loss rate for genes in the host species will be:
$$
r_l^*= r_l + \beta \frac{|\gamma(h)| N}{M}
$$
Where, $\alpha$ and $\beta$ are user-provided normalization factors.

The new loss rates will be higher in species with more genes than the average. Furthermore, the loss will be more likely to happen in a symbiont whenever $\alpha > \beta$ and conversely, loss will be more likely to happen in the host whenever $\beta>\alpha$.

This update should be interpreted at the level of **event rates** or **hazards**, not as a direct post-hoc probability correction. Increasing $r_l^*$ changes both the relative probability that the next event in a branch is a loss and the total event intensity used to sample the next event time. Therefore, in the implementation the effective loss rates should be recomputed from the current state of the holobiont system whenever gene content changes.

### Practical safeguards

The rate update should be activated only for true holobiont systems, i.e. when $N>1$. Otherwise a single host lineage with no active symbiont would still receive an interaction-driven loss penalty. Moreover, if duplications accumulate too many genes in one species, the term $\frac{|\gamma(\cdot)|N}{M}$ may become too aggressive. In that case it is reasonable to keep in reserve a capped crowding factor such as
$$
\operatorname{crowding}(x)= \min\left(\frac{|\gamma(x)|N}{M}, c_{\max}\right),
$$
with $c_{\max}$ around $2$ or $3$.

Independent auxiliary branches that are not part of any host-symbiont system should stay at the base loss rate $r_l$.

### Practical defaults for $\alpha$ and $\beta$

It is convenient to define $\alpha$ and $\beta$ as fractions of the base loss rate $r_l$, rather than as fixed absolute constants. This keeps the effect comparable across simulations run with different base loss rates.

A conservative starting point for a host-keeps regime is:
$$
\alpha = 0.25\,r_l, \qquad \beta = 0.10\,r_l.
$$

A neutral baseline is:
$$
\alpha = \beta = 0.15\,r_l.
$$

A stronger exploratory regime is:
$$
\alpha = 0.50\,r_l, \qquad \beta = 0.20\,r_l.
$$

The first option is the safest default for initial experiments, because around average occupancy $\frac{|\gamma(\cdot)|N}{M}\approx 1$, it increases symbiont losses by about $25\%$ of $r_l$ and host losses by about $10\%$ of $r_l$ without overwhelming the existing duplication, loss, transfer, and conversion dynamics.

## Update rates based in symbiont-host interactions

We will decrease transfer probability between species residing in different hots.

Given a gene growing branch $x_0\in \Pi$  at time $t$, together with the corresponding species $y_0=\kappa(x_0)$ and host $z_0=\mu'(y_0)$. The probability of transfer from $y_0$ to another branch $y_1\in E(T_A)$ decreases whenever 

- **(C0)** $y_1 \not\in \gamma'_t(z_0)$

for time-consistency purposes, the probability of transfer between $y_0$ and $y_1=ab$ is different from zero if only if $\tau_A(b) \leq t \le \tau_A(a)$.
