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


