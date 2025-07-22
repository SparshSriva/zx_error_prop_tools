import pyzx as zx
from pyzx.pauliweb import PauliWeb, compute_pauli_webs
from typing import List, Tuple, Dict
import matplotlib.pyplot as plt
import pprint
import re

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
    sc_map = generate_rotated_surface_code(d)
    print(f"--- Generated Grid Map for d={d} ---")
    pprint.pprint(sc_map)

    data_qubits = {}
    x_ancillas = {}
    z_ancillas = {}

    size = len(sc_map)
    for r, row_list in enumerate(sc_map):
        for c, label in enumerate(row_list):
            if label != 0:
                match = re.match(r'([dXZ])(\d+)', label)
                if match:
                    q_type, q_index = match.groups()
                    q_index = int(q_index)

                    y_coord = (size - 1) - r
                    coordinate = (c, y_coord)
                    
                    if q_type == 'd':
                        data_qubits[q_index] = {'coord': coordinate, 'label': label}
                    elif q_type == 'X':
                        x_ancillas[q_index] = {'coord': coordinate, 'label': label}
                    elif q_type == 'Z':
                        z_ancillas[q_index] = {'coord': coordinate, 'label': label}

    sorted_data = [data_qubits[i] for i in sorted(data_qubits.keys())]

    all_ancillas = list(x_ancillas.values()) + list(z_ancillas.values())

    all_ancillas.sort(key=lambda item: (-item['coord'][1], item['coord'][0]))

    pyzx_qubit_map = [item['coord'] for item in sorted_data] + \
                     [item['coord'] for item in all_ancillas]

    all_qubits_for_plot = sorted_data + all_ancillas
    
    plt.figure(figsize=(10, 10))
    
    data_x = [c['coord'][0] for c in sorted_data]
    data_y = [c['coord'][1] for c in sorted_data]
    plt.scatter(data_x, data_y, s=200, facecolors='lightblue', edgecolors='black', label='Data Qubits')

    anc_x_x = [c['coord'][0] for c in x_ancillas.values()]
    anc_x_y = [c['coord'][1] for c in x_ancillas.values()]
    plt.scatter(anc_x_x, anc_x_y, s=220, marker='s', facecolors='lightcoral', edgecolors='black', label='X Ancillas')

    anc_z_x = [c['coord'][0] for c in z_ancillas.values()]
    anc_z_y = [c['coord'][1] for c in z_ancillas.values()]
    plt.scatter(anc_z_x, anc_z_y, s=220, marker='s', facecolors='lightgreen', edgecolors='black', label='Z Ancillas')

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

    plt.gca().set_aspect('equal', adjustable='box')
    plt.legend()
    plt.show()

    return pyzx_qubit_map, all_ancillas

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
    sc_map = generate_rotated_surface_code(d)
    data_qubits, ancilla_qubits = parse_qubit_map(sc_map)

    num_data = len(data_qubits)
    num_ancilla = len(ancilla_qubits)

    connectivity = {}
    for ancilla in ancilla_qubits:
        r_matrix, c_matrix = ancilla['matrix_coord']
        neighbor_indices = find_neighboring_data_qubits(sc_map, r_matrix, c_matrix)
        connectivity[ancilla['label']] = neighbor_indices

    qasm_header = f"// Surface Code d={d}\n"
    qasm_header += "OPENQASM 2.0;\n"
    qasm_header += 'include "qelib1.inc";\n'
    qasm_header += f"qreg q[{num_data}];\n"
    qasm_header += f"qreg qq[{num_ancilla}];\n"
    
    qasm_body = ""
    for i, ancilla in enumerate(ancilla_qubits):
        ancilla_label = ancilla['label']
        ancilla_type = ancilla_label[0]
        
        qasm_body += f"\n// Stabilizer for {ancilla_label} (qq[{i}])\n"
        
        connected_qubits = connectivity[ancilla_label]

        cnot_order = list(range(len(connected_qubits)))
        if custom_cnot_orderings and ancilla_label in custom_cnot_orderings:
            cnot_order = custom_cnot_orderings[ancilla_label]
            qasm_body += f"// Using custom CNOT order: {cnot_order}\n"

        ordered_qubits = [connected_qubits[k] for k in cnot_order]

        if ancilla_type == 'Z':
            for data_q_idx in ordered_qubits:
                qasm_body += f"cx qq[{i}], q[{data_q_idx}];\n"
        elif ancilla_type == 'X':
            qasm_body += f"h qq[{i}];\n"
            for data_q_idx in ordered_qubits:
                qasm_body += f"cx q[{data_q_idx}], qq[{i}];\n"
            qasm_body += f"h qq[{i}];\n"   
    return qasm_header + qasm_body

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