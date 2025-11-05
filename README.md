# Evolutionary Dynamics in Subdivided Populations

This repository contains code and notebooks to simulate **the evolutionary dynamics of subdividied populations** using a Gillespie algorithm and to compute the **fixation probability)**, **the mean number of local fixations**, **the mean number of dispersal events**, and **the fixation time**.

This repository accompanies the following paper:
> **L. Marrec**, *Impact of genotype-dependent dispersal on mutation fixation in subdivided populations*, *Journal of Evolutionary Biology*, in press (2025).  
> BioRxiv preprint: [doi.org/10.1101/2023.11.29.569213](https://doi.org/10.1101/2023.11.29.569213)

If you find this code useful, **please cite the paper above**.

---

## Repository Structure

### 1. `C/`
Contains **individual-based simulations** (birth, death, dispersal) in C.

Subfolders:
- **Island/**
- **SteppingStone/**

Each subfolder contains:
- `main.c` – Implements the full simulation based on the algorithm described in the paper.

**Usage:**
```bash
# Navigate to the subfolder
cd C/Island  # or C/SteppingStone

# Compile the program
gcc main.c -o main

# Run the simulation
./main
```
### 2. `Matlab/`
This folder contains **deme-level simulations** using a Gillespie algorithm.  
These scripts simulate the dynamics at the level of demes rather than individuals, which allows **faster simulations** and **more stochastic replicates**. The results are consistent with the C simulations.

Subfolders:
- **Island/**
- **SteppingStone/**

Each subfolder contains:

| File | Description |
|------|-------------|
| `RunFixationAnalysis.m` | Main script to run the simulation and generate plots. |
| `GillespieDemeFixation.m` | Gillespie algorithm simulating fixation probability (mutants can either go extinct or fix). |
| `GillespieDemeFixationConditioned.m` | Gillespie algorithm conditioned on mutant fixation, used to calculate fixation time, number of dispersal events, and number of local fixations. |

**Usage:**

1. Open Matlab.
2. Navigate to the folder corresponding to the model:
```matlab
cd 'Matlab/Island'       % or 'Matlab/SteppingStone'

