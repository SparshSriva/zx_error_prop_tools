# ZX-Calculus Tools for Analyzing Hook Errors in Surface Codes

This project provides a set of tools and analyses for studying the effect of CNOT gate ordering on the propagation of "hook errors" in surface codes. It uses a combination of PyZX for ZX-calculus-based error propagation and Stim for noisy simulations and validation.

The core idea is that the order of CNOTs in a stabilizer measurement circuit can affect how physical errors propagate into logical errors. This toolkit was developed to identify CNOT orderings that are more robust against specific, problematic hook errors.

## Methodology

Two primary methods of analysis are used in this project:

1.  **PyZX `gflow` Analysis:** This is a deterministic, graph-based method. An error is injected into a ZX-diagram of a circuit, and its propagation is tracked using the `gflow` algorithm. This allows for a precise understanding of the error chain, but can be computationally intensive for large codes. The main analysis is in `hook errors.ipynb`.

2.  **Stim Simulation:** Stim is used for two purposes:
    *   **Noisy Memory Experiment:** A full memory experiment is simulated with depolarizing noise to see how different CNOT orderings perform in a more realistic, noisy environment. This is found in `stim_d3_sc.ipynb`.
    *   **Classical Pauli Propagation:** A custom classical simulator was developed to mimic the exact Pauli propagation of the PyZX analysis. This is used for validation and for analyzing larger codes like the d=5 case in `stim_d5_sc.ipynb`.

## Repository Structure

The repository has been refactored into a small Python library and a set of example notebooks.

-   `zx_error_prop/`: This directory contains the core Python library.
    -   `surface_code.py`: Functions for generating surface code layouts and QASM circuits.
    -   `propagation.py`: Functions for running the PyZX `gflow` analysis and counting logical errors.
    -   `stim_interface.py`: Functions for running the Stim-based memory experiments and classical simulations.
    -   `utils.py`: General helper functions.

-   `*.ipynb`: Jupyter notebooks that serve as the entry points for the analyses.
    -   `hook errors.ipynb`: The main PyZX-based analysis for the d=3 code.
    -   `stim_d3_sc.ipynb`: The noisy memory experiment for the d=3 code.
    -   `stim_d5_sc.ipynb`: The classical hook error simulation for the d=5 code.

## Getting Started

### 1. Setup Environment

It is recommended to use a virtual environment.

```bash
python -m venv venv
source venv/bin/activate  # On Windows, use `venv\\Scripts\\activate`
```

Install the required packages using the `requirements.txt` file:

```bash
pip install -r requirements.txt
```

### 2. Run the Analyses

To run the analyses, open the Jupyter notebooks and execute the cells.

-   **For the main PyZX analysis:**
    Open and run `hook errors.ipynb`.

-   **For the noisy Stim simulation (d=3):**
    Open and run `stim_d3_sc.ipynb`.

-   **For the classical simulation (d=5):**
    Open and run `stim_d5_sc.ipynb`.

You can modify the parameters at the top of each notebook to explore different code distances, noise levels, or hook error configurations.
