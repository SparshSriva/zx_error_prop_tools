from typing import List
import numpy as np
import pyzx as zx
from pyzx import VertexType
from numpy import sum, log2, abs, real, diag
from numpy.random import choice

def sample_bitstrings(prob_vector: np.ndarray, n_samples: int) -> List[str]:
    """
    Samples bitstrings based on the given probability vector.

    Args:
        prob_vector (np.ndarray): A vector where each entry represents the
                                  probability of a bitstring.
        n_samples (int): The number of samples to generate.

    Returns:
        List[str]: A list of sampled bitstrings.
    """

    prob_vector = prob_vector / sum(prob_vector)

    sampled_indices = choice(len(prob_vector), size=n_samples, p=prob_vector)

    bit_length = int(log2(len(prob_vector)))
    sampled_bitstrings = [format(idx, f'0{bit_length}b') for idx in sampled_indices]
    
    return sampled_bitstrings

def graphical_partial_trace(graph: zx.Graph, qubits: List[int]) -> zx.Graph:
    """
    Computes the graphical partial trace of a graph over specified qubits.

    Args:
        graph (zx.Graph): The input graph, must be a state.
        qubits (List[int]): List of qubit indices to trace out.

    Returns:
        zx.Graph: New graph after performing the partial trace.
    """
    gc = graph.copy().adjoint() + graph.copy()
    outs = gc.outputs()
    ins = gc.inputs()
    for q in qubits:
        gc.set_type(outs[q], VertexType.X)
        gc.set_type(ins[q], VertexType.X)
        gc.add_edge((outs[q], ins[q]))
    return gc

def sampler(graph: zx.Graph, qubits: List[int], n_samples: int) -> List[str]:
    """
    Samples bitstrings from the given graph.

    Args:
        graph (zx.Graph): The input graph, must be a state.
        qubits (List[int]): List of qubit indices to trace out.
                            Sample over remaining qubits.
        n_samples (int): The number of samples to generate.

    Returns:
        List[str]: A list of sampled bitstrings.
    """
    traced = graphical_partial_trace(graph, qubits)
    sample = abs(real(diag(traced.to_matrix(preserve_scalar=True))))
    return sample_bitstrings(diag(sample), n_samples)