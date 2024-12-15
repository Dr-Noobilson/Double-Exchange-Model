# Double-Exchange-Model

Monte Carlo Simulation for Simplified Double Exchange Model

## Overview

This repository contains the implementation of a Monte Carlo simulation for the simplified Double Exchange (DE) model, focusing on the interaction between localized spins and itinerant electrons and verifying the ferromagnetic order in the lowest energy configuration. The project computes key physical quantities such as magnetization, spin-spin correlations, and phase diagrams across different lattice sizes.

## Files in the Repository

### 1. **mcmc_de.ipynb**
   - Implements the Metropolis algorithm for simulating the double exchange model.
   - Handles temperature sweeps, thermalization, state updates and plotting of results.

### 2. **utils.py**
   - Contains utility functions for key computations:
     - **`generate_lattice`**: Generates the spin lattice with specified dimensions.
     - **`calculate_mag`**: Computes localized spin magnetization of the system.
     - **`metropolis`**: Runs the Metropolis algorithm for energy state updates.
     - **`calculate_hamiltonian`**: Computes the Hamiltonian of the lattice.
     - **`compute_eigenvalues`**: Determines eigenvalues and eigenstates for energy computation.
     - **`calculate_dq0`**: Calculates spin-spin correlation for given lattice points in reciprocal space.
     - **`quantum_mag`**: Estimates hoping electron magnetization.

### 3. **plotter.py**
   - Provides functions to visualize the results:
     - Visualization of the lattice spins
     - Heatmaps for spin-spin correlation in reciprocal lattice space.

## Results

Simulations have been conducted for **4×4** and **6×6** lattices. The key findings include:
- **Localized Spin Magnetization vs. Temperature**: Displays the behavior of localized spins as a function of temperature at different filling fractions of itinerant electrons.
- **Hopping Electron Magnetization vs. Temperature**: Investigates the magnetization of itinerant electrons vs temperature and their filling fractions.
- **Spin-Spin Correlation vs. Temperature**: Captures correlations between spins at q=0 (large order) vs Temperature at different filling fractions.
- **Heatmaps**: Visualizes spin-spin correlations for various reciprocal lattice vectors \(q\).
- **Phase Diagram**: Heatmap plot for temperaure, filling fraction and magnetization that identifies the different phases of the system.

## Requirements

- Python 3.x
- Required libraries:
  - `numpy`
  - `matplotlib`
  - `scipy`
  - `tqdm`

Install the dependencies using:
```bash
pip install -r requirements.txt
```

## How to Use

```bash
# Clone the repository:
git clone https://github.com/<your-username>/Double-Exchange-Model.git
cd Double-Exchange-Model
```

# Run the Monte Carlo simulation using the Jupyter Notebook:
```
jupyter notebook mcmc_de.ipynb
```
