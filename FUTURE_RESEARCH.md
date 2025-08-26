# Future Research Directions

This document outlines potential future research directions that build upon the tools and analyses in this repository. The core methodology of analyzing CNOT ordering effects on error propagation can be extended in several exciting ways.

## 1. Deeper Analysis of the Surface Code

The current analysis can be extended even for the surface code itself.

### Mixed CNOT Orderings

The current notebooks apply the same CNOT permutation to all weight-4 stabilizers. A key question from the original presentation was: *what is the effect of different orderings between plaquettes?*

A future experiment could explore this by assigning different permutations to different stabilizers. The `generate_surface_code_qasm` function in the library already supports this, as it accepts a dictionary mapping stabilizer labels to their CNOT orders. An analysis could search the space of mixed orderings to see if a combination of "good" and "bad" orderings can yield even better results.

### Different Error Models

The analysis currently focuses on single Pauli errors injected at specific times (hook errors). This could be expanded to:
-   **Two-qubit errors:** Injecting correlated errors on two qubits after a CNOT gate.
-   **More complex noise channels:** Using the noisy Stim simulator to model more realistic noise, such as coherent errors or leakage.

## 2. Extension to Other Quantum Error Correcting Codes

The methodology developed here is not limited to the surface code. In principle, it can be applied to any QECC with a syndrome extraction circuit that can be translated to the ZX-calculus.

Some promising candidates for future analysis include:

-   **Color Codes:** These are another family of topological codes with interesting properties, such as the ability to perform some logical gates transversally. Analyzing their CNOT ordering robustness would be a valuable comparison to the surface code.
-   **Concatenated Codes:** Codes like the Steane code concatenated with another code are important for fault-tolerance. The ZX-calculus is well-suited to representing these hierarchical structures.
-   **Graph Codes:** As this entire analysis is graph-based, extending it to the broader family of graph codes is a natural next step. The literature suggests ZX-calculus is a powerful tool for constructing and understanding these codes.

A proof-of-concept for analyzing a different code could involve writing a new function (e.g., `generate_color_code_qasm`) and then applying the same simulation notebooks to it.

## 3. Deeper Integration with Stim

Stim is an extremely fast simulator, and the current analysis only scratches the surface of its capabilities.

### Gflow in Stim

The presentation noted that Craig Gidney has his own definition of gflow. A major research direction would be to investigate if this gflow definition can be used to perform the same kind of deterministic error propagation directly within Stim. This would combine the analytical power of the gflow approach with the performance of an industry-grade tool, potentially allowing for the analysis of much larger and more complex codes.

### Advanced Stim Features

Stim's `detector` and `observable` features were used in a basic way in the memory experiment. A deeper integration could involve more complex detector setups to analyze different error types or to design more efficient decoding schemes based on the results of the adversarial error analysis.

## 4. Analytical Results with PyZX

While Stim is fast for sampling, PyZX is a research tool that allows for graph-based reasoning and rewriting. A long-term goal could be to move beyond simulation and towards *analytical* results. For example, it might be possible to use the graph rewriting rules of the ZX-calculus to *prove* that certain CNOT orderings are robust, rather than just observing it through simulation. This would be a powerful result and a significant contribution to the field.
