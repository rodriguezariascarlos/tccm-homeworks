# Installation Instructions

## Requirements
- Fortran compiler (gfortran or equivalent)
- BLAS library installed (optional but recommended for DGEMM)

## Compilation
1. Navigate to the main project directory.
2. Run the following command: gfortran -o matrix_mult src/sparse.f90 src/dense_operations.f90 src/main.f90 -lblas

## Execution
1. Run the program: ./matrix_mult
2. The program will read the matrices from tests/... and display the results of the multiplication and efficiency analysis.
