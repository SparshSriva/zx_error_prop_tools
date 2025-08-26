import re
import pprint
from typing import List, Dict, Any, Tuple

import pyzx as zx
import matplotlib.pyplot as plt


def generate_rotated_surface_code(distance: int) -> List[List[Any]]:
    """
    Generates a 2D qubit map for a rotated surface code of a given distance.

    This layout is drawn on a square grid where data qubits are at odd-numbered
    coordinates and measurement qubits are at even-numbered coordinates.

    Args:
        distance (int): The distance of the surface code. Must be an odd
                        integer >= 3.

    Returns:
        A 2D list representing the qubit map. '0' indicates an empty
        location, 'd' qubits are data qubits, and 'X'/'Z' qubits are
        measurement qubits for the stabilizers.
    """
    if not isinstance(distance, int) or distance < 3 or distance % 2 == 0:
        raise ValueError("Distance d must be an odd integer >= 3.")

    size = 2 * distance + 1
    sc_map = [[0] * size for _ in range(size)]

    # Place data qubits at (r, c) where r and c are odd
    data_idx = 0
    for r in range(1, 2 * distance, 2):
        for c in range(1, 2 * distance, 2):
            sc_map[r][c] = f'd{data_idx}'
            data_idx += 1

    # Place measurement qubits (stabilizers)
    x_idx, z_idx = 0, 0
    # Place Z stabilizers
    for i in range(1, distance):
        for j in range(distance + 1):
            if (i + j) % 2 == 0:
                r, c = 2 * i, 2 * j
                if r < size and c < size:
                    sc_map[r][c] = f'Z{z_idx}'
                    z_idx += 1
    # Place X stabilizers
    for i in range(distance + 1):
        for j in range(1, distance):
            if (i + j) % 2 != 0:
                r, c = 2 * i, 2 * j
                if r < size and c < size:
                    sc_map[r][c] = f'X{x_idx}'
                    x_idx += 1
    return sc_map


def parse_qubit_map(sc_map: List[List[Any]]) -> Tuple[List[Dict], List[Dict]]:
    """
    Parses the grid map to extract and sort qubit information.

    Args:
        sc_map: The 2D list representing the surface code map.

    Returns:
        A tuple containing two lists:
        - A list of data qubit information dictionaries.
        - A list of ancilla qubit information dictionaries, sorted by their
          position (top-to-bottom, then left-to-right).
    """
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
                    if q_type == 'd':
                        data_qubits[int(q_index)] = info
                    elif q_type == 'X':
                        x_ancillas[int(q_index)] = info
                    elif q_type == 'Z':
                        z_ancillas[int(q_index)] = info

    sorted_data = [data_qubits[i] for i in sorted(data_qubits.keys())]

    all_ancillas = list(x_ancillas.values()) + list(z_ancillas.values())
    all_ancillas.sort(key=lambda item: (-item['coord'][1], item['coord'][0]))

    return sorted_data, all_ancillas


def find_neighboring_data_qubits(sc_map: List[List[Any]], r_matrix: int, c_matrix: int) -> List[int]:
    """
    Finds data qubit neighbors for a given ancilla coordinate in a rotated layout.

    In the rotated layout, data qubits are diagonal to the ancilla qubits.

    Args:
        sc_map: The 2D list representing the surface code map.
        r_matrix: The row index of the ancilla in the map.
        c_matrix: The column index of the ancilla in the map.

    Returns:
        A sorted list of integer indices for the neighboring data qubits.
    """
    neighbors = []
    size = len(sc_map)
    # Check all four diagonal positions for data qubits
    for dr, dc in [(-1, -1), (-1, 1), (1, -1), (1, 1)]:
        nr, nc = r_matrix + dr, c_matrix + dc
        if 0 <= nr < size and 0 <= nc < size and isinstance(sc_map[nr][nc], str):
            if sc_map[nr][nc].startswith('d'):
                neighbors.append(sc_map[nr][nc])
    # Return sorted list of data qubit indices (e.g., [0, 1, 3, 4])
    return sorted([int(re.search(r'\d+', n).group()) for n in neighbors])


def generate_surface_code_qasm(
    distance: int,
    custom_cnot_orderings: Dict[str, List[int]] = None
) -> str:
    """
    Generates a QASM string for a surface code circuit.

    Allows for configurable CNOT ordering for each stabilizer, which is key
    to studying the impact of gate ordering on hook errors.

    Args:
        distance: The distance of the code.
        custom_cnot_orderings: A dictionary to specify non-default CNOT orders.
            The keys are ancilla labels (e.g., 'Z0') and values are the
            desired order of application to neighbor qubits (e.g., [3, 2, 1, 0]).

    Returns:
        The complete QASM 2.0 circuit string.
    """
    sc_map = generate_rotated_surface_code(distance)
    data_qubits, ancilla_qubits = parse_qubit_map(sc_map)

    num_data = len(data_qubits)
    num_ancilla = len(ancilla_qubits)

    connectivity = {}
    for ancilla in ancilla_qubits:
        r_matrix, c_matrix = ancilla['matrix_coord']
        neighbor_indices = find_neighboring_data_qubits(sc_map, r_matrix, c_matrix)
        connectivity[ancilla['label']] = neighbor_indices

    qasm_header = f"// Surface Code d={distance}\n"
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

        # Use the default order [0, 1, 2...] unless a custom one is provided
        cnot_order = list(range(len(connected_qubits)))
        if custom_cnot_orderings and ancilla_label in custom_cnot_orderings:
            cnot_order = custom_cnot_orderings[ancilla_label]
            qasm_body += f"// Using custom CNOT order: {cnot_order}\n"

        ordered_qubits = [connected_qubits[k] for k in cnot_order]

        if ancilla_type == 'Z':
            # For Z-stabilizers, CNOT from data to ancilla
            for data_q_idx in ordered_qubits:
                qasm_body += f"cx q[{data_q_idx}], a[{i}];\n"
        elif ancilla_type == 'X':
            # For X-stabilizers, H on ancilla, then CNOT from ancilla to data
            qasm_body += f"h a[{i}];\n"
            for data_q_idx in ordered_qubits:
                qasm_body += f"cx a[{i}], q[{data_q_idx}];\n"
            qasm_body += f"h a[{i}];\n"

    return qasm_header + qasm_body


def plot_rotated_code(distance: int) -> None:
    """
    Generates and plots the layout for a rotated surface code.

    Args:
        distance: The distance of the surface code.
    """
    sc_map = generate_rotated_surface_code(distance)
    print(f"--- Generated Grid Map for d={distance} ---")
    pprint.pprint(sc_map)

    data_qubits, all_ancillas = parse_qubit_map(sc_map)
    x_ancillas = [anc for anc in all_ancillas if anc['label'].startswith('X')]
    z_ancillas = [anc for anc in all_ancillas if anc['label'].startswith('Z')]

    all_qubits_for_plot = data_qubits + all_ancillas

    plt.figure(figsize=(10, 10))

    # Plot data qubits
    data_x = [c['coord'][0] for c in data_qubits]
    data_y = [c['coord'][1] for c in data_qubits]
    plt.scatter(data_x, data_y, s=200, facecolors='lightblue', edgecolors='black', label='Data Qubits')

    # Plot X ancillas
    anc_x_x = [c['coord'][0] for c in x_ancillas]
    anc_x_y = [c['coord'][1] for c in x_ancillas]
    plt.scatter(anc_x_x, anc_x_y, s=220, marker='s', facecolors='lightcoral', edgecolors='black', label='X Ancillas')

    # Plot Z ancillas
    anc_z_x = [c['coord'][0] for c in z_ancillas]
    anc_z_y = [c['coord'][1] for c in z_ancillas]
    plt.scatter(anc_z_x, anc_z_y, s=220, marker='s', facecolors='lightgreen', edgecolors='black', label='Z Ancillas')

    for item in all_qubits_for_plot:
        x, y = item['coord']
        label = item['label']
        plt.text(x, y, label, fontsize=9, ha='center', va='center', weight='bold')

    plt.title(f'Surface Code Qubit Map (d={distance})', fontsize=16)
    plt.xlabel('X-coordinate', fontsize=12)
    plt.ylabel('Y-coordinate', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)

    max_coord = 2 * distance
    plt.xticks(range(max_coord + 1))
    plt.yticks(range(max_coord + 1))

    plt.gca().set_aspect('equal', adjustable='box')
    plt.legend()
    plt.show()
