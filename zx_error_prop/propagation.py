from typing import List, Tuple, Dict, Any

import pyzx as zx
import networkx as nx
from pyzx.pauliweb import PauliWeb, compute_pauli_webs


def propagate_pauli_error(
    graph: zx.Graph,
    initial_errors: List[Tuple[Tuple[int, int], str]],
    add_initial_errors_to_web: bool = True
) -> Tuple[PauliWeb, List[str]]:
    """
    Propagates a list of initial errors through a circuit using gflow.

    This computes the full propagation and correction paths based on the
    Pauli webs of the circuit's ZX-diagram.

    Args:
        graph: The ZX-diagram of the circuit.
        initial_errors: A list of tuples, where each tuple defines an error
                        as ((vertex1, vertex2), error_type).
        add_initial_errors_to_web: If True, the output web includes the
                                   initial errors. If False, it includes only
                                   the resulting corrections.

    Returns:
        A tuple containing:
        - A PauliWeb object with the full error/correction paths.
        - A log of the corrections that were applied.
    """
    try:
        _, zwebs, xwebs = compute_pauli_webs(graph, backwards=True)
    except (ValueError, KeyError) as e:
        print(f"Error: Could not compute gflow for the graph: {e}")
        return None, []

    # Create a web to hold only the initial injected errors.
    initial_error_web = PauliWeb(graph)
    for half_edge, error_type in initial_errors:
        initial_error_web.add_half_edge(half_edge, error_type)

    # Create a separate web to accumulate the corrections.
    correction_web = PauliWeb(graph)
    log = []

    # Determine corrections by checking commutation of webs against the initial errors.
    output_edges = [(o, next(iter(graph.neighbors(o)))) for o in graph.outputs()]

    for o, n in output_edges:
        # Check Z-web. If it anti-commutes, an X correction is needed.
        if o in zwebs and not zwebs[o].commutes_with(initial_error_web):
            log.append(f"Z-web for output {o}: anti-commutes, added X")
            correction_web.add_half_edge((o, n), 'X')
        elif o in zwebs:
            log.append(f"Z-web for output {o}: commutes")

        # Check X-web. If it anti-commutes, a Z correction is needed.
        if o in xwebs and not xwebs[o].commutes_with(initial_error_web):
            log.append(f"X-web for output {o}: anti-commutes, added Z")
            correction_web.add_half_edge((o, n), 'Z')
        elif o in xwebs:
            log.append(f"X-web for output {o}: commutes")

    if add_initial_errors_to_web:
        # Combine initial errors and corrections for the full picture.
        final_web = initial_error_web * correction_web
    else:
        # Return only the corrections.
        final_web = correction_web

    return final_web, log


def get_output_errors(graph: zx.Graph, final_web: PauliWeb) -> Dict[int, str]:
    """
    Extracts the final Pauli errors at the circuit's outputs from a PauliWeb.

    Args:
        graph: The ZX-diagram of the circuit.
        final_web: The PauliWeb containing the full error and correction paths.

    Returns:
        A dictionary mapping each output vertex index to its resulting error
        ('X', 'Y', or 'Z'). Returns an empty dictionary if the web is None.
    """
    if not final_web:
        return {}

    output_errors = {}
    error_dict = final_web.es  # .es stores errors on each half-edge

    for v_out in graph.outputs():
        v_neighbor = next(iter(graph.neighbors(v_out)))
        half_edge_key = (v_out, v_neighbor)
        error = error_dict.get(half_edge_key, 'I')

        if error != 'I':
            output_errors[v_out] = error

    return output_errors


def get_qubit_errors(
    graph: zx.Graph,
    qubits_to_check: List[Dict[str, Any]],
    final_web: PauliWeb
) -> Dict[str, str]:
    """
    Extracts the final Pauli errors on a specified list of output qubits.

    Args:
        graph: The ZX-diagram of the circuit.
        qubits_to_check: A list of qubit info dictionaries, where each dict
                         contains the qubit's 'coord' and 'label'.
        final_web: The PauliWeb containing the full error and correction paths.

    Returns:
        A dictionary mapping each qubit's label to its resulting error
        ('X', 'Y', or 'Z').
    """
    if not final_web:
        return {}

    coord_to_vertex = {g.qubit_coordinates(v): v for v in graph.vertices()}
    output_errors = {}
    error_dict = final_web.es

    for qubit_info in qubits_to_check:
        qubit_coord = qubit_info['coord']
        qubit_label = qubit_info['label']

        v_idx = coord_to_vertex.get(qubit_coord)
        if v_idx is None or not graph.is_output(v_idx):
            continue

        v_neighbor = next(iter(graph.neighbors(v_idx)))
        half_edge_key = (v_idx, v_neighbor)
        error = error_dict.get(half_edge_key, 'I')

        if error != 'I':
            output_errors[qubit_label] = error

    return output_errors


def get_syndrome(
    graph: zx.Graph,
    ancillas: List[Dict[str, Any]],
    final_web: PauliWeb
) -> List[Dict[str, Any]]:
    """
    Extracts the error syndrome by identifying which ancilla measurements flipped.

    A flipped measurement is called a "defect".

    Args:
        graph: The ZX-diagram of the circuit.
        ancillas: A list of the ancilla qubit data (labels and coordinates).
        final_web: The PauliWeb containing the full error and correction paths.

    Returns:
        A list of ancilla qubit dictionaries that correspond to the defects.
    """
    ancilla_errors = get_qubit_errors(graph, ancillas, final_web)
    syndrome = []

    for ancilla_info in ancillas:
        label = ancilla_info['label']
        error = ancilla_errors.get(label, 'I')
        ancilla_type = label[0]

        is_defect = False
        # A Z-basis measurement is flipped by an X or Y error.
        if ancilla_type == 'Z' and error in ('X', 'Y'):
            is_defect = True
        # An X-basis measurement is flipped by a Z or Y error.
        elif ancilla_type == 'X' and error in ('Z', 'Y'):
            is_defect = True

        if is_defect:
            syndrome.append(ancilla_info)

    return syndrome


def create_matching_graph(syndrome: List[Dict[str, Any]]) -> nx.Graph:
    """
    Creates a graph for use with a Minimum Weight Perfect Matching (MWPM) decoder.

    The nodes in the graph are the defective ancillas from the syndrome. The
    edge weights are the Manhattan distance between the ancillas on the code's grid.

    Args:
        syndrome: A list of defective ancilla qubits.

    Returns:
        A networkx graph.
    """
    matching_graph = nx.Graph()
    for i in range(len(syndrome)):
        for j in range(i + 1, len(syndrome)):
            ancilla1 = syndrome[i]
            ancilla2 = syndrome[j]

            coord1 = ancilla1['coord']
            coord2 = ancilla2['coord']

            distance = abs(coord1[0] - coord2[0]) + abs(coord1[1] - coord2[1])
            matching_graph.add_edge(i, j, weight=distance)

    return matching_graph


def count_logical_errors(
    output_errors: Dict[int, str],
    z_logicals: List[List[int]],
    x_logicals: List[List[int]]
) -> Tuple[int, int]:
    """
    Counts the number of logical Z and X errors.

    A logical error occurs if a Pauli error string across the output qubits
    anti-commutes with a logical operator. This is checked by seeing if the
    error string crosses a logical line an odd number of times.

    Args:
        output_errors: A dictionary mapping output qubit indices to their
                       Pauli error ('X', 'Y', or 'Z').
        z_logicals: A list of Z-type logical operators, where each operator is
                    a list of qubit indices it acts on.
        x_logicals: A list of X-type logical operators.

    Returns:
        A tuple containing (number_of_z_errors, number_of_x_errors).
    """
    num_z_errors = 0
    num_x_errors = 0

    # A Z-logical operator is triggered by an X or Y error.
    for z_line in z_logicals:
        z_crossings = 0
        for qubit_idx in z_line:
            if qubit_idx in output_errors and output_errors[qubit_idx] in ('X', 'Y'):
                z_crossings += 1
        if z_crossings % 2 != 0:
            num_z_errors += 1

    # An X-logical operator is triggered by a Z or Y error.
    for x_line in x_logicals:
        x_crossings = 0
        for qubit_idx in x_line:
            if qubit_idx in output_errors and output_errors[qubit_idx] in ('Z', 'Y'):
                x_crossings += 1
        if x_crossings % 2 != 0:
            num_x_errors += 1

    return num_z_errors, num_x_errors
