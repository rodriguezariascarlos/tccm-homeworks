program main
  use sparse                ! Use the 'sparse' module, which likely contains functions related to sparse matrices
  use dense_operations      ! Use the 'dense_operations' module, which likely contains functions for dense matrix operations
  implicit none             ! Disables implicit declaration of variables, requiring all variables to be explicitly declared
    
  ! Declare variables
  type(sparse_matrix) :: A, B, C_sparse
  real(c_double), allocatable :: A_dense(:,:), B_dense(:,:), C_dense(:,:)
  real :: start, finish ! Variables to measure execution time
  integer :: actual_mults, theoretical_mults  ! Counters for multiplications (real vs. theoretical)
  integer :: i, repetitions ! Loop counter and number of repetitions
  character(len=100) :: fileA, fileB ! File paths for input matrices
 
  repetitions = 1000

  fileA = 'tests/matrix'
  fileB = 'tests/matrix'

  ! Load matrices from files
  call read_matrix(fileA, A)
  call read_matrix(fileB, B)

  ! Convert sparse matrices to dense
  call sparse_to_dense(A, A_dense)
  call sparse_to_dense(B, B_dense)

  ! Allocate memory for the result matrix C (dense version)
  allocate(C_dense(size(A_dense, 1), size(B_dense, 2)))

  ! Dense matrix multiplication using BLAS
  call cpu_time(start)
  do i = 1, repetitions
  call multiply_dense_blas(A_dense, B_dense, C_dense)
  end do
  call cpu_time(finish)
  print *, "Multiplication with Dense Matrices (BLAS): Average Time =", (finish - start)/repetitions

  ! Dense matrix multiplication with manual code (without BLAS)
  call cpu_time(start)
  do i = 1, repetitions
  call multiply_dense_manual(A_dense, B_dense, C_dense)
  end do
  call cpu_time(finish)
  print *, "Multiplication with Dense Matrices (hand-written): Average Time =", (finish - start)/repetitions

  ! Sparse matrix multiplication
  call cpu_time(start)
  do i = 1, repetitions
  call multiply_matrices(A, B, C_sparse, actual_mults)
  end do
  call cpu_time(finish)
  theoretical_mults = A%rows**3
  print *, "Multiplication with Sparse Matrices: Average Time =", (finish - start)/repetitions

  ! Print the results and performance metrics
  print *, "Multiplications performed:", actual_mults
  print *, "Theorical maximum of multiplications:", theoretical_mults
  print *, "Filling rate of the resulting matrix: ", C_sparse%non_zero / real(C_sparse%rows * C_sparse%cols) * 100.0, "%"
  print *, "Efficiency of sparse multiplication: ", actual_mults / real(theoretical_mults) * 100.0, "%"

  print *, "Result:"
  call print_matrix(C_sparse)

end program main
