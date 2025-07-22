from numpy import sum, log2
from numpy.random import choice
from pyzx import Graph, VertexType
from numpy import abs, real, diag
from typing import Dict, Any, Optional, Tuple, List

def sample_bitstrings(prob_vector, n_samples):
    """
    Samples bitstrings based on the given probability vector.

    Parameters:
        prob_vector (numpy.ndarray): A vector where each entry represents the probability of a bitstring.
        n_samples (int): The number of samples to generate.

    Returns:
        list: A list of sampled bitstrings.
    """

    prob_vector = prob_vector / sum(prob_vector)

    sampled_indices = choice(len(prob_vector), size=n_samples, p=prob_vector)

    bit_length = int(log2(len(prob_vector)))
    sampled_bitstrings = [format(idx, f'0{bit_length}b') for idx in sampled_indices]
    
    return sampled_bitstrings

def graphical_partial_trace(graph: Graph, qubits: list) -> Graph:
    """
    Computes the graphical partial trace of a graph over specified qubits.

    Parameters:
        graph (Graph): The input graph, must be a state.
        qubits (list): List of qubit indices to trace out.

    Returns:
        Graph: New graph after performing the partial trace.
    """
    gc = graph.copy().adjoint() + graph.copy()
    outs = gc.outputs()
    ins = gc.inputs()
    for q in qubits:
        gc.set_type(outs[q], VertexType.X)
        gc.set_type(ins[q], VertexType.X)
        gc.add_edge((outs[q], ins[q]))
    return gc

def sampler(graph: Graph, qubits: list, n_samples: int) -> list:
    """
    Samples bitstrings from the given graph.

    Parameters:
        graph (Graph): The input graph, must be a state.
        qubits (list): List of qubit indices to trace out. 
                       Sample over remaining qubits.
        n_samples (int): The number of samples to generate.

    Returns:
        list: A list of sampled bitstrings. 
    """
    traced = graphical_partial_trace(graph, qubits)
    sample = abs(real(diag(traced.to_matrix(preserve_scalar=True))))
    return sample_bitstrings(diag(sample), n_samples)