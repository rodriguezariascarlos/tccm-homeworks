# Project: Hartree-Fock and MP2 Energy Calculation

## Overview
This project implements a program in C to compute the Hartree-Fock (HF) energy and the second-order Møller-Plesset perturbation theory (MP2) correlation energy using molecular orbital integrals stored in a TREXIO file.

## Features
- Reads one-electron and two-electron integrals from a TREXIO file.
- Computes nuclear repulsion energy.
- Calculates the Hartree-Fock energy based on occupied molecular orbitals.
- Performs MP2 energy calculations by summing over occupied and virtual molecular orbitals.
- Efficiently handles memory allocation and deallocation for 4D arrays.

## File Structure
/project
│
├── src/            
│   └── main.c      # Main program implementing Hartree-Fock and MP2 calculations
│
├── tests/          
│   └── *.h5        # Test molecular files in TREXIO format
│
├── AUTHORS         # List of contributors
├── LICENSE         # License information
├── README.md       # Project documentation
└── INSTALL.md      # Installation instructions

## Usage
1. Ensure that the `tests` folder contains valid TREXIO files with molecular orbital integrals.
2. Compile the program following the steps in `INSTALL.md`.
3. Run the program by providing the path to a TREXIO file as input:
   ```bash
   ./main <file.h5>
