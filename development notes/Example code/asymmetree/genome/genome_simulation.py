"""Module for the simulation of genomes."""

from __future__ import annotations

from pathlib import Path

from tralda.datastructures import Tree
from tralda.datastructures import TreeNode

from asymmetree.io.sequences import write_alignment
from asymmetree.io.sequences import write_fasta
from asymmetree.treeevolve import gene_trees
from asymmetree.treeevolve import prune_losses
from asymmetree.seqevolve.alignment import AlignmentBuilder
from asymmetree.seqevolve.evolver import Evolver
from asymmetree.seqevolve.substitution import SubstModel
from asymmetree.utils.sampling import Sampler


class GenomeSimulator:
    """Class for the simulation of genomes."""

    def __init__(self, species_tree: Tree, outdir: str | Path = None) -> None:
        """Constructor for the GenomeSimulator class.

        Args:
            species_tree: The species tree along which gene trees are simulated.
            outdir: The path to a directory into which the results of the simulation (serialized
                trees, one fasta file per genome, one 'true' alignment per gene family). The default
                is None, in which case nothing is saved to file.
        """
        if not isinstance(species_tree, Tree):
            raise TypeError("species tree must be of type 'Tree'")

        self.S = species_tree
        self.outdir = Path(outdir) if outdir is not None else None

        if self.outdir:
            self._check_outdir()
            self.S.serialize(self.outdir / "species_tree.json", mode="json")

        self.true_gene_trees = []
        self.pruned_gene_trees = []

    def simulate_gene_trees(self, n: int, **kwargs):
        """Simulate gene trees.

        Simulates the 'true gene trees' which still contain the loss events as well as the 'pruned
        gene trees' in which the loss branches are removed.

        For the parameters of the simulation, see function 'gene_trees' in 'treeevolve' subpackage
        or documentation.

        Args:
            n: Number of gene trees to be simulated.
        """
        self.number_of_families = n

        self.true_gene_trees = gene_trees(self.S, n=n, **kwargs)

        self.pruned_gene_trees = [prune_losses(tree) for tree in self.true_gene_trees]

        # sequences should be emptied here if methods were called before
        if hasattr(self, "sequence_dicts"):
            self.sequence_dicts.clear()

        if self.outdir:
            for i in range(n):
                filename = self.outdir / "true_gene_trees" / f"gene_tree{i}.json"
                self.true_gene_trees[i].serialize(filename, mode="json")

        return self.true_gene_trees, self.pruned_gene_trees

    def simulate_sequences(
        self,
        subst_model: SubstModel,
        root_genome: list | None = None,
        length_distr: tuple | Sampler = ("constant", 200),
        min_length: int = 10,
        max_length: int | None = None,
        write_fastas: bool = True,
        write_alignments: bool = True,
        **kwargs,
    ) -> None:
        """Simulate sequences along the (pruned) gene trees.

        For the additional parameters of the simulation, see class 'Evolver' in the 'seqevolve'
        subpackage or documentation.

        Args:
            subst_model: Substitution model.
            root_genome: List of sequences for the roots of the gene trees, must contain the
                same number of str sequences as trees that were simulated. The default is None,
                in which case new sequences are generated at random. The sequences must be
                compatible with the specified substitution model (model_type='n'/'a').
            length_distr: Distribution of the length of the root sequences if root genome is not
                supplied. See documentation for available options.
            min_length: Minimal length at which the distribution of lengths is truncated. Must be
                less than the mean of this distribution.
            max_length: Maximal length at which the distribution of lengths is truncated. Must be
                less than the mean of this distribution. The default is None, in which case the
                distribution is not truncated.
            write_fastas: If True and an output directory was specified, write the sequences (one
                file per species) into the directory 'fasta_files' in the output directory.
            write_alignments: If True and an output directory was specified, write the true
                alignments (one file per gene tree) into the directory 'alignments' in the output
                directory.
        """
        self.subst_model = subst_model

        if hasattr(self, "sequence_dicts"):
            self.sequence_dicts.clear()
        else:
            self.sequence_dicts = []

        if root_genome:
            if len(root_genome) != self.number_of_families:
                raise ValueError(
                    "number of sequences in root genome does not match number of gene families"
                )
        else:
            if isinstance(length_distr, Sampler):
                self.sampler = length_distr
            else:
                self.sampler = Sampler(
                    length_distr, minimum=min_length, maximum=max_length, discrete=True
                )

        evolver = Evolver(subst_model, **kwargs)

        for i in range(self.number_of_families):
            PGT = self.pruned_gene_trees[i]

            if root_genome:
                evolver.evolve_along_tree(PGT, start_seq=root_genome[i].upper())
            else:
                evolver.evolve_along_tree(PGT, start_length=self.sampler())

            self.sequence_dicts.append(evolver.sequences)

            # write one alignment file per gene family
            if self.outdir and write_alignments:
                self._write_alignment(i)

        # write one fasta file per species
        if self.outdir and write_fastas:
            self._write_fastas(include_inner=False)

    def _check_outdir(self) -> None:
        """Check the output directory and create necessary subdirectories.

        Raises:
            FileExistsError: If the specified output directory exists but is not a directory.
        """
        if not self.outdir.exists():
            self.outdir.mkdir(parents=True, exist_ok=True)
        elif self.outdir.exists() and not self.outdir.is_dir():
            raise FileExistsError(f"'{self.outdir}' is not a directory")

        for directory in ("true_gene_trees", "fasta_files", "alignments"):
            path = self.outdir / directory
            if not path.exists():
                path.mkdir(parents=True, exist_ok=True)

    def _write_alignment(self, family_id: int) -> None:
        """Write the true alignment for a given gene family to file.

        Args:
            family_id: The id of the gene family for which the alignment should be written.

        Raises:
            RuntimeError: If no output directory was specified for the alignments.
        """
        if not self.outdir:
            raise RuntimeError("no output directory specified for alignments")

        alg_builder = AlignmentBuilder(
            self.pruned_gene_trees[family_id],
            self.sequence_dicts[family_id],
            self.subst_model.alphabet,
            include_inner=False,
        )

        alignment = []
        for node, sequence in alg_builder.build().items():
            label = self._compose_label(node, family_id)
            alignment.append((label, sequence))

        basename = "alignment{}.phylip".format(family_id)
        filename = self.outdir / "alignments" / basename
        write_alignment(filename, alignment, alignment_format="phylip")

    def _write_fastas(self, include_inner: bool = False) -> None:
        """Write the sequences for all species to fasta files.

        Args:
            include_inner: If True, also write sequences for inner nodes of the species tree. The
                default is False, in which case only sequences for the extant species (leaves of
                the species tree) are written.

        Raises:
            RuntimeError: If no output directory was specified for the fasta files.
        """
        if not self.outdir:
            raise RuntimeError("no output directory specified for fasta files")

        # list of all species labels
        if include_inner:
            species = [v.label for v in self.S.preorder()]
        else:
            species = [v.label for v in self.S.preorder() if not v.children]

        # labeled sequences sorted by species
        sorted_seqs = {s: [] for s in species}

        for i in range(self.number_of_families):
            for node, evoseq in self.sequence_dicts[i].items():
                # skip inner nodes
                if not include_inner and node.children:
                    continue

                label = self._compose_label(node, i)
                sequence = self.subst_model.to_sequence(evoseq)
                sorted_seqs[node.reconc].append((label, sequence))

        for spec, sequences in sorted_seqs.items():
            basename = f"{spec}.f{self.subst_model.model_type}a"
            filename = self.outdir / "fasta_files" / basename
            write_fasta(filename, sequences)

    def _compose_label(self, node: TreeNode, family_id: int) -> str:
        """Compose a label for a given node and gene family id.

        Args:
            node: The node for which the label should be composed.
            family_id: The id of the gene family to which the node belongs.

        Returns:
            The composed label for the given node and gene family id.
        """
        return f"fam{family_id}gene{node.label}spec{node.reconc}"
