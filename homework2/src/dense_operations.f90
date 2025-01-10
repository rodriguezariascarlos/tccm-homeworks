module dense_operations
  use sparse, only: sparse_matrix  ! Use the 'sparse' module to access the sparse_matrix type and its functions
  use iso_c_binding, only: c_double  ! Use the 'iso_c_binding' module to define 'c_double' for compatibility with C
  implicit none  ! Disables implicit declaration of variables, requiring all variables to be explicitly declared

contains

 ! Subroutine to convert a sparse matrix to a dense matrix
 subroutine sparse_to_dense(sparse, dense)
    implicit none
    type(sparse_matrix), intent(in) :: sparse   ! Input: sparse matrix
    real(c_double), allocatable, intent(out) :: dense(:,:)  ! Output: dense matrix
    integer :: i  ! Loop variable for non-zero elements in the sparse matrix

    ! Allocate memory for the dense matrix (matching the size of the sparse matrix)
    allocate(dense(sparse%rows, sparse%cols))
    
    ! Initialize the dense matrix to zero
    dense = 0.0_c_double

    ! Loop over the non-zero elements in the sparse matrix and populate the dense matrix
    do i = 1, sparse%non_zero
      ! Fill both (row, col) and (col, row) positions in the dense matrix (for symmetric matrices)  
      dense(sparse%row(i), sparse%col(i)) = sparse%val(i)
      if (sparse%row(i) /= sparse%col(i)) then
        dense(sparse%col(i), sparse%row(i)) = sparse%val(i)
      end if
    end do
  end subroutine sparse_to_dense

  ! Subroutine to multiply dense matrices using BLAS (high-performance optimized library)
  subroutine multiply_dense_blas(A, B, C)
    implicit none
    real(c_double), intent(in) :: A(:,:), B(:,:)  ! Input dense matrices A and B
    real(c_double), intent(out) :: C(:,:)  ! Output dense matrix C (result of A * B)
    real(c_double) :: alpha, beta  ! Scalars for matrix multiplication (used by dgemm)
    integer :: m, n, k, lda, ldb, ldc  ! Matrix dimensions and leading dimensions for the BLAS call

    ! Set matrix dimensions based on the input matrices A and B
    m = size(A, 1)  ! Number of rows of A
    k = size(A, 2)  ! Number of columns of A (also the number of rows of B)
    n = size(B, 2)  ! Number of columns of B

    ! Set leading dimensions (for memory storage optimization in BLAS)
    lda = m
    ldb = k
    ldc = m

    ! Set scalar values for the multiplication (alpha = 1 for multiplication, beta = 0 for no addition)
    alpha = 1.0_c_double
    beta = 0.0_c_double
 
    ! Call BLAS's dgemm function to perform matrix multiplication
    call dgemm('N', 'N', m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
  end subroutine multiply_dense_blas

  ! Subroutine to multiply dense matrices manually (without using BLAS)
  subroutine multiply_dense_manual(A, B, C)
    implicit none
    real(c_double), intent(in) :: A(:,:), B(:,:)  ! Input dense matrices A and B
    real(c_double), intent(out) :: C(:,:)  ! Output dense matrix C (result of A * B)
    integer :: i, j, k  ! Loop variables for matrix multiplication

    ! Initialize matrix C to zero
    C = 0.0_c_double

    ! Loop over rows of A, columns of B, and the common dimension (k) for the multiplication
    do i = 1, size(A, 1)  ! Iterate over rows of A
      do j = 1, size(B, 2)  ! Iterate over columns of B
        do k = 1, size(A, 2)  ! Iterate over the common dimension (columns of A and rows of B)
          C(i, j) = C(i, j) + A(i, k) * B(k, j)  ! Perform the multiplication and accumulate the result in C
        end do
      end do
    end do
  end subroutine multiply_dense_manual
end module dense_operations
