import stim
import itertools
import statistics
import time
from typing import List, Dict, Tuple

# === Constants for d=3 Surface Code ===
# These are hardcoded for the specific d=3 topology used in the notebooks.
# A more general implementation would generate these from the surface_code module.

CONNECTIVITY = {
    'X0': [0, 1],
    'Z0': [0, 1, 3, 4],
    'X1': [1, 2, 4, 5],
    'Z1': [2, 5],
    'X2': [3, 4, 6, 7],
    'Z3': [4, 5, 7, 8],
    'Z2': [3, 6],
    'X3': [7, 8],
}
STAB_TYPES = {k: k[0] for k in CONNECTIVITY.keys()}
DATA_COUNT = 9
# This order is important for the analysis
ANC_ORDER = ['Z0', 'X1', 'X2', 'Z3']
ANC_INDEX = {a: DATA_COUNT + i for i, a in enumerate(ANC_ORDER)}


def build_memory_circuit(
    ordering: List[int],
    rounds: int = 10,
    p2: float = 0.001
) -> stim.Circuit:
    """
    Builds a multi-round memory experiment circuit in Stim.

    Args:
        ordering: The CNOT permutation to use for the weight-4 stabilizers.
        rounds: The number of syndrome extraction rounds.
        p2: The 2-qubit depolarizing error rate to apply after each CNOT.

    Returns:
        A stim.Circuit object for the memory experiment.
    """
    perm_map = {stab: ordering for stab in ANC_ORDER}
    c = stim.Circuit()
    data_qubits = list(range(DATA_COUNT))

    c.append('R', data_qubits)
    c.append('TICK')

    for r in range(rounds):
        # Reset ancillas
        for a in ANC_ORDER:
            q = ANC_INDEX[a]
            c.append('R', [q])
            if STAB_TYPES[a] == 'X':
                c.append('H', [q])
        c.append('TICK')

        # Apply stabilizer checks
        for a in ANC_ORDER:
            base = CONNECTIVITY[a]
            permuted_dqs = [base[i] for i in perm_map[a]]
            anc_q = ANC_INDEX[a]
            if STAB_TYPES[a] == 'Z':
                for dq in permuted_dqs:
                    c.append('CNOT', [dq, anc_q])
                    if p2 > 0:
                        c.append('DEPOLARIZE2', [dq, anc_q], p2)
                    c.append('TICK')
            else:  # X stabilizer
                for dq in permuted_dqs:
                    c.append('CNOT', [anc_q, dq])
                    if p2 > 0:
                        c.append('DEPOLARIZE2', [anc_q, dq], p2)
                    c.append('TICK')

        # Measure ancillas
        for a in ANC_ORDER:
            q = ANC_INDEX[a]
            if STAB_TYPES[a] == 'X':
                c.append('H', [q])
            c.append('M', [q])
        c.append('TICK')

    # Final data measurement for logical error check
    c.append('M', data_qubits)
    return c


def get_logical_failure_rate(circuit: stim.Circuit, shots: int = 10000) -> float:
    """
    Samples a circuit and computes the logical Z failure rate.

    Assumes the logical Z operator is the parity of the top row of data qubits.

    Args:
        circuit: The stim circuit to simulate.
        shots: The number of times to sample the circuit.

    Returns:
        The fraction of shots that resulted in a logical Z flip.
    """
    sampler = circuit.compile_sampler()
    measurements = sampler.sample(shots=shots)

    data_block = measurements[:, -DATA_COUNT:]
    top_row = data_block[:, [0, 1, 2]]

    # Parity of the top row indicates a logical Z flip from the |0>_L state
    parity = top_row.sum(axis=1) & 1
    return parity.mean()


def run_stim_memory_experiment(
    rounds: int = 6,
    p2: float = 0.00002,
    shots: int = 5000,
    repeats: int = 10
) -> Dict[Tuple[int, ...], Dict[str, float]]:
    """
    Runs a full memory experiment for all CNOT permutations and aggregates stats.

    Args:
        rounds: Number of syndrome rounds per circuit.
        p2: 2-qubit depolarizing noise level.
        shots: Number of shots per simulation repeat.
        repeats: Number of times to repeat the simulation for statistics.

    Returns:
        A dictionary mapping each permutation to its statistics (mean, std, stderr).
    """
    permutations = list(itertools.permutations([0, 1, 2, 3]))
    accumulated_results = {p: [] for p in permutations}

    print(f"Running memory experiment for {len(permutations)} permutations.")
    print(f"Params: {rounds} rounds, p2={p2}, {shots} shots, {repeats} repeats.")
    start_time = time.time()

    for i in range(repeats):
        for perm in permutations:
            circuit = build_memory_circuit(list(perm), rounds=rounds, p2=p2)
            fail_rate = get_logical_failure_rate(circuit, shots=shots)
            accumulated_results[perm].append(fail_rate)

        if (i + 1) % max(1, repeats // 5) == 0:
            print(f"  Completed {i+1}/{repeats} repeats in {time.time()-start_time:.1f}s")

    stats = {}
    for perm, values in accumulated_results.items():
        mean = statistics.fmean(values)
        sd = statistics.pstdev(values) if len(values) > 1 else 0.0
        stderr = (sd / (len(values)**0.5)) if len(values) > 1 else 0.0
        stats[perm] = {'mean': mean, 'std_dev': sd, 'std_err': stderr}

    return stats


# ===== Classical Pauli Propagation Simulator =====
# This section contains functions to run a noise-free, classical simulation
# of Pauli error propagation, mimicking the behavior of the PyZX gflow analysis.
# This is useful for validating results and for analyses where the exact
# error propagation path is more important than performance.

import re
from .surface_code import generate_rotated_surface_code, parse_qubit_map, find_neighboring_data_qubits, generate_surface_code_qasm

def _p_encode(p: str) -> Tuple[int, int]:
    """Encodes a Pauli operator string to a binary vector."""
    return {'I': (0,0), 'X': (1,0), 'Z': (0,1), 'Y': (1,1)}[p]

def _p_decode(xz: Tuple[int, int]) -> str:
    """Decodes a binary vector to a Pauli operator string."""
    x, z = xz
    if (x,z) == (0,0): return 'I'
    if (x,z) == (1,0): return 'X'
    if (x,z) == (0,1): return 'Z'
    return 'Y'

def _qasm_to_schedule(qasm_str: str, expected_data: int) -> Tuple:
    """Parses a QASM string into a structured schedule for simulation."""
    lines = [ln.strip() for ln in qasm_str.splitlines()]
    m_q = next((re.match(r"qreg q\\[(\\d+)\\];", ln) for ln in lines if ln.startswith('qreg q[')), None)
    m_a = next((re.match(r"qreg a\\[(\\d+)\\];", ln) for ln in lines if ln.startswith('qreg a[')), None)
    num_data = int(m_q.group(1)) if m_q else 0
    num_anc = int(m_a.group(1)) if m_a else 0
    if num_data != expected_data:
        print(f"WARNING: QASM header reports num_data={num_data} (expected {expected_data})")

    anc_label_to_aidx = {}
    per_anc = {}
    order = []
    cur = None
    for ln in lines:
        if ln.startswith('// Stabilizer for '):
            m = re.match(r"// Stabilizer for\\s+([XZ]\\d+)\\s+\\(a\\[(\\d+)\\]\\)", ln)
            if not m: continue
            cur = m.group(1)
            aidx = int(m.group(2))
            anc_label_to_aidx[cur] = aidx
            per_anc[cur] = {'H_pre': False, 'cnots': [], 'H_post': False}
            order.append(cur)
        elif ln.startswith('h a['):
            if cur is None: continue
            if not per_anc[cur]['H_pre']:
                per_anc[cur]['H_pre'] = True
            else:
                per_anc[cur]['H_post'] = True
        elif ln.startswith('cx '):
            if cur is None: continue
            m1 = re.match(r"cx\\s+q\\[(\\d+)\\],\\s*a\\[(\\d+)\\];", ln)
            m2 = re.match(r"cx\\s+a\\[(\\d+)\\],\\s*q\\[(\\d+)\\];", ln)
            if m1:
                dq = int(m1.group(1))
                per_anc[cur]['cnots'].append(('D->A', dq))
            elif m2:
                dq = int(m2.group(2))
                per_anc[cur]['cnots'].append(('A->D', dq))
    return num_data, num_anc, anc_label_to_aidx, per_anc, order

def _simulate_hook_paulis_qasm(
    distance: int,
    ordering: List[int],
    target_stab: str,
    hook_step: int,
    error_type: str,
    weight4_stabs: List[str]
) -> List[str]:
    """Simulates a single hook error propagation for a given ordering."""
    perm4 = list(ordering)
    custom_orders = {lab: perm4 for lab in weight4_stabs}
    qasm = generate_surface_code_qasm(distance, custom_cnot_orderings=custom_orders)
    num_data, num_anc, anc_label_to_aidx, per_anc, anc_order = _qasm_to_schedule(qasm, distance*distance)
    paulis = {q: (0,0) for q in range(num_data + num_anc)}

    def H_on(q):
        x,z = paulis[q]; paulis[q] = (z,x)
    def CNOT_(c,t):
        xc,zc = paulis[c]; xt,zt = paulis[t]
        xt ^= xc; zc ^= zt
        paulis[c] = (xc,zc); paulis[t] = (xt,zt)
    def add_pauli(q, p):
        ax,az = paulis[q]; bx,bz = _p_encode(p)
        paulis[q] = (ax^bx, az^bz)

    for lab in anc_order:
        steps = per_anc[lab]
        aidx = anc_label_to_aidx[lab]
        aq = num_data + aidx
        if steps['H_pre']: H_on(aq)
        for i, (kind, dq) in enumerate(steps['cnots']):
            if kind=='D->A': CNOT_(dq, aq)
            else: CNOT_(aq, dq)
            if lab == target_stab and i == hook_step:
                add_pauli(aq, error_type)
        if steps['H_post']: H_on(aq)

    return [_p_decode(paulis[q]) for q in range(num_data)]

def run_classical_hook_simulation(distance: int, hook_step: int) -> Dict:
    """
    Runs a classical simulation to find robust CNOT orderings for a given code distance.

    This analysis identifies which CNOT permutations are robust to hook errors by
    simulating the exact Pauli propagation for each case and checking the final
    weight of error chains on the logical operator lines.

    Args:
        distance: The distance of the surface code.
        hook_step: The CNOT step (0-indexed) after which to inject the hook error.

    Returns:
        A dictionary containing the analysis results, including the robust permutations.
    """
    # 1. Build code connectivity
    sc_map = generate_rotated_surface_code(distance)
    _, ancillas = parse_qubit_map(sc_map)
    connectivity = {anc['label']: find_neighboring_data_qubits(sc_map, *anc['matrix_coord']) for anc in ancillas}
    ancilla_types = {anc['label']: anc['label'][0] for anc in ancillas}
    weight4_stabs = sorted([lab for lab, neigh in connectivity.items() if len(neigh) == 4], key=lambda s: (s[0], int(s[1:])))

    # 2. Define logical lines and error types
    D2 = distance * distance
    Z_LINES = [list(range(r*distance, r*distance + D2)) for r in range(distance)]
    X_LINES = [list(range(c, D2, distance)) for c in range(distance)]
    ERROR_TYPE_BY_STAB_TYPE = {'X': 'X', 'Z': 'Z'}

    # 3. Loop through all permutations and simulate
    permutations = list(itertools.permutations([0, 1, 2, 3]))
    results = {}
    print(f"Evaluating {len(permutations)} permutations for d={distance} at hook_step={hook_step}")
    for p in permutations:
        max_z_weight, max_x_weight = 0, 0
        for stab in weight4_stabs:
            err_type = ERROR_TYPE_BY_STAB_TYPE[ancilla_types[stab]]
            data_paulis = _simulate_hook_paulis_qasm(distance, p, stab, hook_step, err_type, weight4_stabs)

            z_weights = [sum(1 for q in row if data_paulis[q] in ('Z','Y')) for row in Z_LINES]
            x_weights = [sum(1 for q in col if data_paulis[q] in ('X','Y')) for col in X_LINES]

            max_z_weight = max(max_z_weight, max(z_weights) if z_weights else 0)
            max_x_weight = max(max_x_weight, max(x_weights) if x_weights else 0)
        results[p] = (max_z_weight, max_x_weight)

    # 4. Identify robust permutations (max line weight < 2 is a good heuristic)
    robust_perms = sorted([p for p, (mz, mx) in results.items() if mz <= 1 and mx <= 1])

    return {
        "distance": distance,
        "hook_step": hook_step,
        "weight4_stabilizers": weight4_stabs,
        "all_results": results,
        "robust_permutations": robust_perms
    }
