import pyzx as zx
from pyzx.pauliweb import PauliWeb, compute_pauli_webs
from typing import List, Tuple, Dict
import matplotlib.pyplot as plt
import pprint
import re
import matplotlib.pyplot as plt
import pprint
import re
import networkx as nx

def generate_rotated_surface_code(d):
    """
    Generates a 2D qubit map for a rotated surface code of a given distance d,
    drawn on a square grid.

    Args:
        d (int): The distance of the surface code. Must be an odd integer >= 3.

    Returns:
        list[list]: A 2D list representing the qubit map. 0 indicates an
                    empty location. 'd' qubits are data qubits, 'X' and 'Z'
                    qubits are measurement qubits for the stabilizers.
    """
    if not isinstance(d, int) or d < 3 or d % 2 == 0:
        raise ValueError("Distance d must be an odd integer greater than or equal to 3.")

    # The grid size for this layout is (2d+1) x (2d+1)
    size = 2 * d + 1
    sc_map = [[0] * size for _ in range(size)]

    # --- Place Data Qubits ---
    # Data qubits are located at (row, col) where both row and col are odd.
    data_idx = 0
    for r in range(1, 2 * d, 2):
        for c in range(1, 2 * d, 2):
            sc_map[r][c] = f'd{data_idx}'
            data_idx += 1

    # --- Place Measurement Qubits (Stabilizers) ---
    # This layout places measurement qubits at (row, col) where both are even.
    x_idx = 0
    z_idx = 0

    # Place Z stabilizers
    for i in range(1, d):
        for j in range(d + 1):
            if (i + j) % 2 == 0:
                r, c = 2 * i, 2 * j
                if r < size and c < size:
                    sc_map[r][c] = f'Z{z_idx}'
                    z_idx += 1

    # Place X stabilizers
    for i in range(d + 1):
        for j in range(1, d):
            if (i + j) % 2 != 0:
                r, c = 2 * i, 2 * j
                if r < size and c < size:
                    sc_map[r][c] = f'X{x_idx}'
                    x_idx += 1

    return sc_map


def generate_pyzx_and_plot_rotated_code(d):
    """
    Generates a PyZX qubit map and plots the layout for a rotated surface code.

    Args:
        d (int): The distance of the surface code.

    Returns:
        list[tuple]: A list of (x, y) coordinates for PyZX, with data qubits
                     listed first, followed by positionally sorted ancilla qubits.
    """
    # 1. Generate the basic grid map
    sc_map = generate_rotated_surface_code(d)
    print(f"--- Generated Grid Map for d={d} ---")
    pprint.pprint(sc_map)

    # Dictionaries to store qubits before sorting them by their index
    data_qubits = {}
    x_ancillas = {}
    z_ancillas = {}
    
    # 2. Iterate through the map to extract qubit info
    size = len(sc_map)
    for r, row_list in enumerate(sc_map):
        for c, label in enumerate(row_list):
            if label != 0:
                match = re.match(r'([dXZ])(\d+)', label)
                if match:
                    q_type, q_index = match.groups()
                    q_index = int(q_index)
                    
                    # FIX: Convert matrix row 'r' to Cartesian y-coordinate
                    # The origin (0,0) is now at the bottom-left.
                    y_coord = (size - 1) - r
                    coordinate = (c, y_coord)
                    
                    if q_type == 'd':
                        data_qubits[q_index] = {'coord': coordinate, 'label': label}
                    elif q_type == 'X':
                        x_ancillas[q_index] = {'coord': coordinate, 'label': label}
                    elif q_type == 'Z':
                        z_ancillas[q_index] = {'coord': coordinate, 'label': label}

    # 3. Sort data qubits by their index
    sorted_data = [data_qubits[i] for i in sorted(data_qubits.keys())]
    
    # 4. Create a single list of all ancillas and sort by position
    all_ancillas = list(x_ancillas.values()) + list(z_ancillas.values())
    
    # FIX: Sort by y-coordinate descending (top-to-bottom), then x ascending (left-to-right)
    all_ancillas.sort(key=lambda item: (-item['coord'][1], item['coord'][0]))

    # 5. Create the final PyZX map
    pyzx_qubit_map = [item['coord'] for item in sorted_data] + \
                     [item['coord'] for item in all_ancillas]

    # 6. Prepare for plotting
    all_qubits_for_plot = sorted_data + all_ancillas
    
    # 7. Create the plot
    plt.figure(figsize=(10, 10))
    
    # Plot data qubits (blue circles)
    data_x = [c['coord'][0] for c in sorted_data]
    data_y = [c['coord'][1] for c in sorted_data]
    plt.scatter(data_x, data_y, s=200, facecolors='lightblue', edgecolors='black', label='Data Qubits')

    # Plot X ancillas (green squares)
    anc_x_x = [c['coord'][0] for c in x_ancillas.values()]
    anc_x_y = [c['coord'][1] for c in x_ancillas.values()]
    plt.scatter(anc_x_x, anc_x_y, s=220, marker='s', facecolors='lightcoral', edgecolors='black', label='X Ancillas')

    # Plot Z ancillas (yellow squares)
    anc_z_x = [c['coord'][0] for c in z_ancillas.values()]
    anc_z_y = [c['coord'][1] for c in z_ancillas.values()]
    plt.scatter(anc_z_x, anc_z_y, s=220, marker='s', facecolors='lightgreen', edgecolors='black', label='Z Ancillas')

    # Add text labels to each point
    for item in all_qubits_for_plot:
        x, y = item['coord']
        label = item['label']
        plt.text(x, y, label, fontsize=9, ha='center', va='center', weight='bold')
    
    plt.title(f'Surface Code Qubit Map (d={d})', fontsize=16)
    plt.xlabel('X-coordinate', fontsize=12)
    plt.ylabel('Y-coordinate', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    
    max_coord = 2 * d
    plt.xticks(range(max_coord + 1))
    plt.yticks(range(max_coord + 1))
    
    # FIX: The y-axis is no longer inverted.
    plt.gca().set_aspect('equal', adjustable='box')
    plt.legend()
    plt.show()

    return pyzx_qubit_map, all_ancillas

def PropagatePauliError(g: zx.Graph, initial_errors: List[Tuple[Tuple[int, int], str]]):
    """
    Takes a list of initial errors, calculates the full propagation and correction
    paths, and returns the final web and a log of corrections.

    Args:
        g: The ZX-diagram of the circuit.
        initial_errors: A list of tuples, where each tuple defines an error as
                        ((vertex1, vertex2), error_type).
    
    Returns:
        A tuple containing (final_web, log), where final_web is a PauliWeb object
        with the full error and correction paths, and log is a list of strings
        describing the corrections made.
    """
    # 1. Create the initial PauliWeb from the list of errors.
    err = PauliWeb(g)
    for half_edge, error_type in initial_errors:
        err.add_half_edge(half_edge, error_type)
        
    print("--- Visualizing Error Propagation and Correction ---")
    print("1. Initial errors introduced in the circuit:")
    zx.draw(g, labels=True, pauli_web=err)

    # 2. Compute the backward-propagated Pauli webs from the outputs.
    try:
        order, zwebs, xwebs = compute_pauli_webs(g, backwards=True)
    except (ValueError, KeyError) as e:
        print(f"Error: Could not compute gflow for the graph: {e}")
        return None, []

    # 3. Create a copy of the web to add corrections to.
    final_web = err.copy()
    
    # 4. Determine and apply corrections at the outputs.
    output_edges = [(o, next(iter(g.neighbors(o)))) for o in g.outputs()]
    
    log = []
    for o, n in output_edges:
        # Check Z-web commutation. If it anti-commutes, add an X correction.
        if o in zwebs and not zwebs[o].commutes_with(final_web):
            log.append(f"Z-web for output {o}: anti-commutes, added X")
            final_web.add_half_edge((o, n), 'X')
        else:
            if o in zwebs: log.append(f"Z-web for output {o}: commutes")
            
        # Check X-web commutation. If it anti-commutes, add a Z correction.
        # This is NOT an `elif` so that Y-errors (which anti-commute with both)
        # get both an X and a Z correction.
        if o in xwebs and not xwebs[o].commutes_with(final_web):
            log.append(f"X-web for output {o}: anti-commutes, added Z")
            final_web.add_half_edge((o, n), 'Z')
        else:
            if o in xwebs: log.append(f"X-web for output {o}: commutes")
            
    return final_web, log

def get_output_errors(g: zx.Graph, final_web: PauliWeb) -> Dict[int, str]:
    """
    Extracts the final Pauli errors at the circuit's outputs from a PauliWeb.

    Args:
        g: The ZX-diagram of the circuit.
        final_web: The PauliWeb containing the full error and correction paths.

    Returns:
        A dictionary mapping each output vertex to its resulting error ('X', 'Y', or 'Z').
        Returns an empty dictionary if the web is None.
    """
    if not final_web:
        return {}
        
    output_errors = {}
    # The .es dictionary stores the errors on each half-edge
    error_dict = final_web.es
    
    for v_out in g.outputs():
        # Find the neighbor of the output vertex that is inside the circuit
        v_neighbor = next(iter(g.neighbors(v_out)))
        
        # The key for the dictionary is the half-edge from the output *into* the circuit
        half_edge_key = (v_out, v_neighbor)
        
        # Get the error, defaulting to 'I' (Identity) if there is no error on that edge
        error = error_dict.get(half_edge_key, 'I')
        
        # We only care about non-Identity errors
        if error != 'I':
            output_errors[v_out] = error
            
    return output_errors

def generate_rotated_surface_code(d):
    """Generates a 2D qubit map for a rotated surface code."""
    if not isinstance(d, int) or d < 3 or d % 2 == 0:
        raise ValueError("Distance d must be an odd integer >= 3.")
    size = 2 * d + 1
    sc_map = [[0] * size for _ in range(size)]
    data_idx, x_idx, z_idx = 0, 0, 0
    for r in range(1, 2 * d, 2):
        for c in range(1, 2 * d, 2):
            sc_map[r][c] = f'd{data_idx}'; data_idx += 1
    for i in range(1, d):
        for j in range(d + 1):
            if (i + j) % 2 == 0:
                r, c = 2 * i, 2 * j
                if r < size and c < size: sc_map[r][c] = f'Z{z_idx}'; z_idx += 1
    for i in range(d + 1):
        for j in range(1, d):
            if (i + j) % 2 != 0:
                r, c = 2 * i, 2 * j
                if r < size and c < size: sc_map[r][c] = f'X{x_idx}'; x_idx += 1
    return sc_map

def parse_qubit_map(sc_map):
    """Parses the grid map to extract and sort qubit information."""
    data_qubits, x_ancillas, z_ancillas = {}, {}, {}
    size = len(sc_map)
    for r_matrix, row_list in enumerate(sc_map):
        for c, label in enumerate(row_list):
            if label != 0:
                match = re.match(r'([dXZ])(\d+)', label)
                if match:
                    q_type, q_index = match.groups()
                    y_coord = (size - 1) - r_matrix
                    info = {'label': label, 'coord': (c, y_coord), 'matrix_coord': (r_matrix, c)}
                    if q_type == 'd': data_qubits[int(q_index)] = info
                    elif q_type == 'X': x_ancillas[int(q_index)] = info
                    elif q_type == 'Z': z_ancillas[int(q_index)] = info
    
    sorted_data = [data_qubits[i] for i in sorted(data_qubits.keys())]
    
    all_ancillas = list(x_ancillas.values()) + list(z_ancillas.values())
    all_ancillas.sort(key=lambda item: (-item['coord'][1], item['coord'][0]))
    
    return sorted_data, all_ancillas

def find_neighboring_data_qubits(sc_map, r_matrix, c_matrix):
    """
    Finds data qubit neighbors for a given ancilla coordinate in a rotated layout.
    FIXED: This now correctly checks diagonal positions.
    """
    neighbors = []
    size = len(sc_map)
    # In this layout, neighbors are always diagonal to the ancilla
    for dr, dc in [(-1, -1), (-1, 1), (1, -1), (1, 1)]: # Check all four diagonal positions
        nr, nc = r_matrix + dr, c_matrix + dc
        # Check bounds and if there is a qubit at the location
        if 0 <= nr < size and 0 <= nc < size and sc_map[nr][nc] != 0:
            # Check if it's a data qubit
            if sc_map[nr][nc].startswith('d'):
                neighbors.append(sc_map[nr][nc])
    # Return sorted list of data qubit indices (e.g., [0, 1, 3, 4])
    return sorted([int(re.search(r'\d+', n).group()) for n in neighbors])

# --- Part 2: QASM Generation ---

def generate_surface_code_qasm(d, custom_cnot_orderings=None):
    """
    Generates a QASM string for a surface code circuit with configurable CNOT ordering.

    Args:
        d (int): The distance of the code.
        custom_cnot_orderings (dict): A dictionary to specify non-default CNOT orders.
            Example: {'Z0': [3, 2, 1, 0]} reverses the order for Z0.

    Returns:
        str: The complete QASM 2.0 circuit string.
    """
    # 1. Get the layout and qubit information
    sc_map = generate_rotated_surface_code(d)
    data_qubits, ancilla_qubits = parse_qubit_map(sc_map)

    num_data = len(data_qubits)
    num_ancilla = len(ancilla_qubits)

    # 2. Determine connectivity for each ancilla
    connectivity = {}
    for ancilla in ancilla_qubits:
        r_matrix, c_matrix = ancilla['matrix_coord']
        neighbor_indices = find_neighboring_data_qubits(sc_map, r_matrix, c_matrix)
        connectivity[ancilla['label']] = neighbor_indices
        
    # 3. Generate the QASM string
    qasm_header = f"// Surface Code d={d}\n"
    qasm_header += "OPENQASM 2.0;\n"
    qasm_header += 'include "qelib1.inc";\n'
    qasm_header += f"qreg q[{num_data}];\n"
    qasm_header += f"qreg a[{num_ancilla}];\n"
    
    qasm_body = ""
    for i, ancilla in enumerate(ancilla_qubits):
        ancilla_label = ancilla['label']
        ancilla_type = ancilla_label[0]
        
        qasm_body += f"\n// Stabilizer for {ancilla_label} (a[{i}])\n"
        
        connected_qubits = connectivity[ancilla_label]
        
        # Determine the CNOT order
        cnot_order = list(range(len(connected_qubits))) # Default order [0, 1, 2...]
        if custom_cnot_orderings and ancilla_label in custom_cnot_orderings:
            cnot_order = custom_cnot_orderings[ancilla_label]
            qasm_body += f"// Using custom CNOT order: {cnot_order}\n"

        # Apply the ordering to the list of connected qubits
        ordered_qubits = [connected_qubits[k] for k in cnot_order]
        
        # Generate the gates for this stabilizer
        if ancilla_type == 'Z':
            for data_q_idx in ordered_qubits:
                qasm_body += f"cx q[{data_q_idx}], a[{i}];\n"
        elif ancilla_type == 'X':
            qasm_body += f"h a[{i}];\n"
            for data_q_idx in ordered_qubits:
                qasm_body += f"cx a[{i}], q[{data_q_idx}];\n"
            qasm_body += f"h a[{i}];\n"   
    return qasm_header + qasm_body

def get_qubit_errors(g: zx.Graph, qubits_to_check: list, final_web: PauliWeb) -> dict:
    """
    Extracts the final Pauli errors on a specified list of qubits that are outputs.

    Args:
        g: The ZX-diagram of the circuit.
        qubits_to_check: A list of qubit info dicts (containing coords and labels).
        final_web: The PauliWeb containing the full error and correction paths.

    Returns:
        A dictionary mapping each qubit's label to its resulting error ('X', 'Y', or 'Z').
    """
    if not final_web:
        return {}

    # Create a reverse mapping from qubit coordinates to vertex indices
    coord_to_vertex = {coord: v for v, coord in g.qubit_map.items()}
    output_errors = {}
    error_dict = final_web.es

    for qubit_info in qubits_to_check:
        qubit_coord = qubit_info['coord']
        qubit_label = qubit_info['label']

        v_idx = coord_to_vertex.get(qubit_coord)
        # Ensure the qubit is a graph output before checking for an error
        if v_idx is None or v_idx not in g.outputs():
            continue

        v_neighbor = next(iter(g.neighbors(v_idx)))
        half_edge_key = (v_idx, v_neighbor)
        error = error_dict.get(half_edge_key, 'I')

        if error != 'I':
            output_errors[qubit_label] = error

    return output_errors

def get_syndrome(g: zx.Graph, ancillas: list, final_web: PauliWeb) -> list:
    """
    Extracts the error syndrome by identifying which ancilla measurements flipped.

    Args:
        g: The ZX-diagram of the circuit.
        ancillas: A list of the ancilla qubit data (labels and coordinates).
        final_web: The PauliWeb containing the full error and correction paths.

    Returns:
        A list of ancilla qubits that form the syndrome (i.e., "defects").
    """
    # Get all non-identity errors on the ancilla qubits
    ancilla_errors = get_qubit_errors(g, ancillas, final_web)

    syndrome = []
    for ancilla_info in ancillas:
        ancilla_label = ancilla_info['label']
        error = ancilla_errors.get(ancilla_label, 'I')

        # Determine if the ancilla measurement is flipped.
        # A Z-stabilizer measures in the Z-basis, so an X or Y error flips it.
        # An X-stabilizer measures in the X-basis, so a Z or Y error flips it.
        ancilla_type = ancilla_label[0]

        is_defect = False
        if ancilla_type == 'Z' and error in ('X', 'Y'):
            is_defect = True
        elif ancilla_type == 'X' and error in ('Z', 'Y'):
            is_defect = True

        if is_defect:
            syndrome.append(ancilla_info)

    return syndrome

def create_matching_graph(syndrome: list) -> nx.Graph:
    """
    Creates a graph for the MWPM decoder.

    Args:
        syndrome: A list of defective ancilla qubits.

    Returns:
        A networkx graph where nodes are defective ancillas and edge weights are
        the Manhattan distance between them.
    """
    matching_graph = nx.Graph()
    for i in range(len(syndrome)):
        for j in range(i + 1, len(syndrome)):
            ancilla1 = syndrome[i]
            ancilla2 = syndrome[j]
            
            coord1 = ancilla1['coord']
            coord2 = ancilla2['coord']
            
            # Calculate Manhattan distance
            distance = abs(coord1[0] - coord2[0]) + abs(coord1[1] - coord2[1])
            
            matching_graph.add_edge(i, j, weight=distance)
            
    return matching_graph