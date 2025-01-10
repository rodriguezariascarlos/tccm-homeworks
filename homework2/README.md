# Project: Sparse and Dense Matrix Multiplication in Fortran

## Overview
This project implements different methods for multiplying sparse and dense matrices using Fortran. The goal is to compare the computational efficiency of each method and analyze their performance based on the size and filling degree of the matrices.

## Features
- Reads sparse matrices from files in the `tests/` folder.
- Converts sparse matrices to dense format.
- Implements matrix multiplication using the following methods:
  - BLAS (DGEMM)
  - Manual algorithms in Fortran
  - Direct multiplication in sparse format
- Measures execution time and computational efficiency for performance analysis.

## File Structure
/project
│
├── src/            
│   ├── sparse.f90              # Handles operations on sparse matrices
│   ├── dense_operations.f90    # Functions for dense matrix operations
│   └── main.f90                # Main program to perform matrix multiplications
│
├── tests/          
│   └── matrix_*                # Test matrices for performance evaluation
│
├── AUTHORS                     # List of contributors
├── LICENSE                     # License information
├── README.md                   # Project documentation
└── INSTALL.md                  # Installation instructions

## Usage
1. Edit the file paths for the input matrices directly in the `main.f90` source code.
2. Compile the program following the steps in `INSTALL.md`.
3. Run the program and analyze the generated results.
