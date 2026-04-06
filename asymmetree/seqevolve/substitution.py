"""Substitution models for nucleotide and amino acid evolution.

The module contains the class `SubstModel` that can be used to represent a substitution model
for nucleotide or amino acid evolution. The class provides methods to compute the transition
probability matrix for a given time and to convert sequences between their str representation and a
list of indices.

For amino acid substitution models, some empirical models are available that can be loaded via the
`EMPIRICAL_MODELS` object. The values in the matrices and the equilibrium frequencies are
taken from [2].

References:
    1. Z. Yang. Computational molecular evolution. Oxford series in ecology and evolution.
       Oxford University Press, 2006. ISBN 978-0-19-856699-1 978-0-19-856702-8.
    2. http://giphy.pasteur.fr/empirical-models-of-amino-acid-substitution/
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy import linalg

from asymmetree.io.substitution_model import parse_paml
from asymmetree.seqevolve.evolving_sequence import EvoSeq


_DATA_DIR: Path = Path(__file__).parent / "data"

# amino acid models available in the package data directory
_EMPIRICAL_MODEL_NAMES: list[str] = ["WAG", "JTT", "BLOSUM62", "LG", "DAYHOFF"]


class _EmpiricalModelLoader:
    """Helper class to load and store empirical amino acid models from PAML files.

    The models are loaded lazily, i.e., only when they are accessed for the first time.
    """

    def __init__(self):
        """Constructor for the _EmpiricalModelLoader class."""
        self.models: dict[str, tuple[np.ndarray, np.ndarray]] = {}

    def __contains__(self, item: str) -> bool:
        """Check if a model is available.

        Args:
            item: The name of the model (e.g. ``"WAG"``).

        Returns:
            True if the model is available, False otherwise.
        """
        return item in _EMPIRICAL_MODEL_NAMES

    def __getitem__(self, key: str) -> tuple[np.ndarray, np.ndarray]:
        """Get a model by name.

        Args:
            key: The name of the model (e.g. ``"WAG"``).

        Returns:
            The 20x20 exchangeability matrix.
            The vector of equilibrium frequencies.

        Raises:
            KeyError: If the model is not available.
        """
        if key not in self:
            raise KeyError(f"Model '{key}' is not available.")

        if key not in self.models:
            self.models[key] = self._load_model(key)

        return self.models[key]

    def _load_model(self, name: str) -> tuple[np.ndarray, np.ndarray]:
        """Load an empirical amino acid substitution model from a PAML file.

        Args:
            name: The name of the model (e.g. ``"WAG"``). A corresponding ``<name>.paml`` file must
                exist in the package data directory.

        Returns:
            The 20x20 exchangeability matrix.
            The vector of equilibrium frequencies.
        """
        matrix, freqs = parse_paml(_DATA_DIR / f"{name}.paml", model_type="a")

        return np.array(matrix), np.array(freqs)


EMPIRICAL_MODELS = _EmpiricalModelLoader()


class SubstModel:
    """Substitution model for nucleotide or amino acid evolution."""

    nuc_models = {"JC69", "K80", "GTR", "CUSTOM"}
    aa_models = {"JC69", "CUSTOM"}

    nucleotides = "ACGT"
    amino_acids = "ARNDCQEGHILKMFPSTWYV"

    def __init__(self, model_type: str, model_name: str, **kwargs) -> None:
        """Constructor for the SubstModel class.

        Some models require additional parameters, e.g. 'K80' requires the parameter 'kappa' and
        'GTR' requires 'abcdef' and 'f', see documentation for details. For a custom model, the
        keyword parameter 'filename' has to be specified with a path to a custom model in paml
        format. All of these parameters are passed as keyword arguments.

        Args:
            model_type: Available options are 'n' for nucleotide and 'a' for amino acid sequences.
            model_name: Available are e.g. 'JC69' (nuc./a.a.), 'K80' (nuc.), 'GTR' (nuc.), 'WAG'
                (a.a.), 'JTT' (a.a.), 'BLOSUM62' (a.a.), 'LG' (a.a.), and 'DAYHOFF' (a.a.). For the
                option 'CUSTOM', the keyword parameter 'filename' has to be specified with a path to
                a custom model in paml format.

        Raises:
            ValueError: If an invalid model was specified.
        """
        self.model_type = model_type.lower()
        self.model_name = model_name.upper()

        self._params = kwargs

        if (
            self.model_type in ("n", "nuc", "nucleotide")
            and self.model_name in SubstModel.nuc_models
        ):
            self.model_type = "n"
            self.alphabet = SubstModel.nucleotides

        elif self.model_type in ("a", "aa", "amino", "aminoacid", "protein") and (
            self.model_name in SubstModel.aa_models or self.model_name in EMPIRICAL_MODELS
        ):
            self.model_type = "a"
            self.alphabet = SubstModel.amino_acids

        else:
            raise ValueError(f"model '{model_type}', '{model_name}' is not available")

        self.alphabet_dict = {item: index for index, item in enumerate(self.alphabet)}

        self._load_exchangeability_and_freqs()
        self._build_rate_matrix()

    def _load_exchangeability_and_freqs(self):
        """Load the exchangeability matrix S and the stationary frequencies π.

        The exchangeability matrix S and the stationary frequencies π are loaded according to the
        specified model. For custom models, they are loaded from a provided PAML file.
        """
        # a custom model (via a paml file) was specified
        if self.model_name == "CUSTOM":
            if self._params is None or "filename" not in self._params:
                raise ValueError("custom model requires the parameter 'filename'")

            S, freqs = parse_paml(self._params["filename"], model_type=self.model_type)
            self.S, self.freqs = np.asarray(S), np.asarray(freqs)

        # non-empirical nucleotide models
        elif self.model_type == "n":
            if self.model_name == "JC69":
                self.S, self.freqs = _JC69_nuc()

            elif self.model_name == "K80":
                if self._params is None or "kappa" not in self._params:
                    raise ValueError("model 'K80' requires the parameter 'kappa'")

                self.S, self.freqs = _K80_nuc(self._params["kappa"])

            elif self.model_name == "GTR":
                if self._params is None or "abcdef" not in self._params or "f" not in self._params:
                    raise ValueError("model 'GTR' requires the parameters 'abcdef' and 'f'")

                self.S, self.freqs = _GTR_nuc(self._params["abcdef"], self._params["f"])

        # non-empirical and empirical amino acid models
        elif self.model_type == "a":
            if self.model_name == "JC69":
                self.S, self.freqs = _JC69_aa()
            elif self.model_name in EMPIRICAL_MODELS:
                self.S, self.freqs = EMPIRICAL_MODELS[self.model_name]

        # make sure stationary frequencies sum to 1
        self.freqs /= np.sum(self.freqs)

        # compute cumulative frequencies for evolver
        self.freqs_cumulative = np.cumsum(self.freqs)

    def _build_rate_matrix(self):
        """Compute the rate matrix Q.

        The rate matrix Q is computed from the exchangeability matrix S and the stationary
        frequencies π according to the formula Q[i, j] = S[i, j] * π[j] for i != j and
        Q[i, i] = -sum(Q[i, :]). The matrix is then normalized so that the expected number of
        substitutions per unit time is 1.
        """
        # rate matrix from exchangeability matrix and stationary frequencies
        Q = self.S @ np.diag(self.freqs)

        # compute diagonals, rows sum to 0
        for i in range(Q.shape[0]):
            Q[i, i] = -(np.sum(Q[i, :]) - Q[i, i])

        # normalize matrix
        k = -np.dot(self.freqs, np.diagonal(Q))
        Q /= k

        self.Q = Q

    def eigensystem(self):
        """Compute the eigensystem of the rate matrix Q."""
        if not hasattr(self, "eigenvals"):
            self.eigenvals, self.U, self.U_inv = diagonalize(self.Q, self.freqs)

        return self.eigenvals, self.U, self.U_inv

    def transition_prob_matrix(self, t: float) -> np.ndarray:
        """Calculate the transition probability matrix P(t).

        The transition probability matrix P(t) is computed from the eigensystem of the rate matrix Q
        according to the formula P(t) = U x e^(Λ * t) x U^(-1) where U is the matrix of eigenvectors
        and Λ is the diagonal matrix of eigenvalues.

        Args:
            t: The time / evolutionary distance for which to compute the transition probability
                matrix.

        Returns:
            The transition probability matrix P(t) as a `numpy.ndarray`.
        """
        # ensure that eigensystem has been computed
        self.eigensystem()

        # first multiplication element-wise, since corresponding matrix only has non-zero entries on
        # the diagonal

        return (self.U * np.exp(self.eigenvals * t)) @ self.U_inv

    def to_indices(self, sequence: str) -> list[int]:
        """Convert a sequence into a list of indices.

        The letters in the sequence are converted to their corresponding indices in the alphabet of
        the model.

        Args:
            sequence: The sequence to be converted.

        Returns:
            The list of indices corresponding to the letters in the sequence.

        Raises:
            ValueError: If the sequence contains a letter that is not in the alphabet of the model.
        """
        try:
            result = [self.alphabet_dict[letter] for letter in sequence]
        except KeyError:
            raise ValueError("invalid sequence for the specified model")

        return result

    def to_sequence(self, evoseq: EvoSeq) -> str:
        """Construct the str representation of a sequence.

        Args:
            evoseq: The sequence to be converted as an `EvoSeq` object.

        Returns:
            The nucleotide or amino acid sequence depending on the model.
        """
        return "".join(self.alphabet[x._value] for x in evoseq)


# --------------------------------------------------------------------------------------------------
#                                        Module functions
# --------------------------------------------------------------------------------------------------


def diagonalize(Q: np.ndarray, freqs: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Diagonalize a matrix.

    Diagonalize the rate matrix Q using the stationary frequencies freqs. If Q is symmetric, it is
    diagonalized directly. If Q is not symmetric but the model is time-reversible, it is first
    transformed to a symmetric matrix using the stationary frequencies and then diagonalized.

    Args:
        Q: Rate matrix to be diagonalized.
        freqs: The stationary frequencies π corresponding to the rate matrix Q.

    Returns:
        A 1-dimensional array containing the eigenvalues.
        A 2-dimensional array containing the eigenvectors as columns. The i-th column corresponds to
            the i-th eigenvalue.
        A 2-dimensional array containing the inverse of the eigenvector matrix.
    """
    # matrix is already symmetric
    if np.allclose(Q, Q.T):
        eigenvals, U = linalg.eigh(Q)
        U_inv = linalg.inv(U)

    # matrix is not symmetric but model is time-reversible
    else:
        Phi = np.diag(np.sqrt(freqs))
        Phi_inv = linalg.inv(Phi)

        B = Phi @ Q @ Phi_inv

        eigenvals, R = linalg.eigh(B)

        U = Phi_inv @ R
        U_inv = linalg.inv(R) @ Phi

    return eigenvals, U, U_inv


# --------------------------------------------------------------------------------------------------
#                                      NUCLEOTIDE MODELS
# --------------------------------------------------------------------------------------------------


def _JC69_nuc() -> tuple[np.ndarray, np.ndarray]:
    """Jukes-Cantor model for nucleotide evolution.

    Returns:
        The 4x4 exchangeability matrix.
        The vector of equilibrium frequencies.
    """
    S = np.array([[0, 1, 1, 1], [1, 0, 1, 1], [1, 1, 0, 1], [1, 1, 1, 0]])

    freqs = np.array([0.25, 0.25, 0.25, 0.25])

    return S, freqs


def _K80_nuc(kappa: float) -> tuple[np.ndarray, np.ndarray]:
    """Kimura 2-parameter model for nucleotide evolution.

    The Kimura 2-parameter model distinguishes between transitions (A <--> G and C <--> T) and
    transversions (all other substitutions) and has a parameter kappa that specifies the
    transition / transversion rate ratio.

    Args:
        kappa: The transition/transversion rate ratio.

    Returns:
        The 4x4 exchangeability matrix.
        The vector of equilibrium frequencies.
    """
    S = np.array([[0, 1, kappa, 1], [1, 0, 1, kappa], [kappa, 1, 0, 1], [1, kappa, 1, 0]])

    freqs = np.array([0.25, 0.25, 0.25, 0.25])

    return S, freqs


def _GTR_nuc(
    abcdef: tuple[float, float, float, float, float, float],
) -> tuple[np.ndarray, np.ndarray]:
    """Generalized time-reversible (GTR) model.

    Parameterization as in e.g. used in PAML and ALF:
    a:    C <--> T
    b:    A <--> T
    c:    G <--> T
    d:    A <--> C
    e:    C <--> G
    f:    A <--> G

    Args:
        abcdef: The six exchangeability parameters a, b, c, d, e, f as defined above.

    Returns:
        The 4x4 exchangeability matrix.
        The vector of equilibrium frequencies.
    """
    a, b, c, d, e, f = abcdef

    S = np.array([[0, d, f, b], [d, 0, e, a], [f, e, 0, c], [b, a, c, 0]])

    freqs = np.asarray(abcdef)
    freqs /= np.sum(freqs)

    return S, freqs


# --------------------------------------------------------------------------------------------------
#                                      AMINO ACID MODELS
# --------------------------------------------------------------------------------------------------


def _JC69_aa() -> tuple[np.ndarray, np.ndarray]:
    """Jukes-Cantor model for amino acid evolution.

    Returns:
        The 20x20 exchangeability matrix.
        The vector of equilibrium frequencies.
    """
    S = np.ones((20, 20))
    S -= np.identity(20)

    freqs = np.full((20,), 1 / 20)

    return S, freqs
