module sparse
  implicit none

  ! Define the sparse_matrix type which will represent a sparse matrix
  type :: sparse_matrix
     integer :: rows, cols, non_zero   ! Number of rows, columns, and non-zero elements
     integer, allocatable :: row(:), col(:)   ! Arrays to store row and column indices of non-zero elements
     real, allocatable :: val(:)        ! Array to store the values of the non-zero elements
  end type sparse_matrix

contains

  ! Subroutine to read a sparse matrix from a file
  subroutine read_matrix(filename, matrix)
    implicit none
    character(len=*), intent(in) :: filename  ! Input filename to read the matrix from
    type(sparse_matrix), intent(out) :: matrix  ! Output sparse matrix
    integer :: i, row, col, max_dim, io_status  ! Loop index, max dimension, I/O status
    real :: val  ! Value of a non-zero element
    integer, allocatable :: temp_row(:), temp_col(:)  ! Temporary arrays for rows and columns
    real, allocatable :: temp_val(:)  ! Temporary array for values of non-zero elements
    integer :: n_entries  ! Number of non-zero entries in the matrix

    ! Open the file for reading
    print *, "Opening file:", filename
    open(unit=10, file=filename, status='old', action='read', iostat=io_status)
    if (io_status /= 0) then  ! Check if file opening was successful
      print *, "Error, unable to open file:", filename
      stop
    end if

    ! Initialize counters and variables
    n_entries = 0
    max_dim = 0

    ! Read through the file and determine the matrix size and the number of non-zero elements
    do
      read(10, *, iostat=io_status) row, col, val
      if (io_status /= 0) exit  ! Exit if end of file or error
      n_entries = n_entries + 1
      if (col > max_dim) max_dim = col  ! Update the maximum column index
    end do
    rewind(10)  ! Rewind the file pointer to the beginning

    ! Allocate temporary arrays to store the matrix data
    allocate(temp_row(n_entries))
    allocate(temp_col(n_entries))
    allocate(temp_val(n_entries))

    ! Read all non-zero elements into temporary arrays
    do i = 1, n_entries
      read(10, *) temp_row(i), temp_col(i), temp_val(i)
    end do
    close(10)  ! Close the file

    ! Set up the output sparse matrix structure
    matrix%rows = max_dim
    matrix%cols = max_dim
    matrix%non_zero = n_entries
    allocate(matrix%row(n_entries))
    allocate(matrix%col(n_entries))
    allocate(matrix%val(n_entries))

    ! Copy the temporary arrays into the matrix structure
    matrix%row = temp_row
    matrix%col = temp_col
    matrix%val = temp_val

    ! Deallocate temporary arrays
    deallocate(temp_row, temp_col, temp_val)

    ! Print matrix information
    print *, "File read successfully"
    print *, "Matrix dimensions: ", matrix%rows, "x", matrix%cols
    print *, "Number of non-zero elements: ", matrix%non_zero
  end subroutine read_matrix

  ! Subroutine to print a sparse matrix
  subroutine print_matrix(matrix)
    implicit none
    type(sparse_matrix), intent(in) :: matrix  ! Input sparse matrix
    integer :: i  ! Loop index for non-zero elements

    ! Print matrix information
    print *, 'Rows:', matrix%rows, 'Cols:', matrix%cols, 'Non-zeros:', matrix%non_zero
    ! Print the row, column, and value of each non-zero element
    do i = 1, matrix%non_zero
      print '(I5, I5, F10.4)', matrix%row(i), matrix%col(i), matrix%val(i)
    end do
  end subroutine print_matrix

  ! Subroutine to resize a sparse matrix to a new size
  subroutine resize_matrix(matrix, new_size)
    implicit none
    type(sparse_matrix), intent(inout) :: matrix  ! Input/output sparse matrix
    integer, intent(in) :: new_size  ! New size for the matrix
    integer, allocatable :: temp_row(:), temp_col(:)  ! Temporary arrays for rows and columns
    real, allocatable :: temp_val(:)  ! Temporary array for values

    ! If the new size is greater than the current size, allocate new arrays
    if (new_size > size(matrix%row)) then
      allocate(temp_row(new_size))
      allocate(temp_col(new_size))
      allocate(temp_val(new_size))

      ! Copy current data into temporary arrays
      temp_row(:size(matrix%row)) = matrix%row(:)
      temp_col(:size(matrix%col)) = matrix%col(:)
      temp_val(:size(matrix%val)) = matrix%val(:)

      ! Deallocate the old arrays
      deallocate(matrix%row, matrix%col, matrix%val)

      ! Assign the new arrays to the matrix
      matrix%row = temp_row
      matrix%col = temp_col
      matrix%val = temp_val
    end if
  end subroutine resize_matrix

  ! Subroutine to multiply two sparse matrices A and B, and store the result in C
  subroutine multiply_matrices(A, B, C, actual_mults)
    implicit none
    type(sparse_matrix), intent(in) :: A, B  ! Input sparse matrices A and B
    type(sparse_matrix), intent(out) :: C  ! Output sparse matrix C (result of A * B)
    integer, intent(out) :: actual_mults  ! Output: number of multiplications performed
    integer :: i, j, nz  ! Loop variables and non-zero element counter
    real :: sum  ! Sum of products for each non-zero element

    ! Initialize multiplication counter and non-zero element counter
    actual_mults = 0
    nz = 0

    ! Set matrix dimensions for the result matrix C
    C%rows = A%rows
    C%cols = B%cols

    ! Allocate memory for result matrix C (in the worst case, every non-zero in A multiplies with every non-zero in B)
    allocate(C%row(A%non_zero * B%non_zero))
    allocate(C%col(A%non_zero * B%non_zero))
    allocate(C%val(A%non_zero * B%non_zero))

    ! Perform matrix multiplication
    do i = 1, A%non_zero
      do j = 1, B%non_zero
        if (A%col(i) == B%row(j)) then  ! If the column of A matches the row of B, multiply the elements
          sum = A%val(i) * B%val(j)
          actual_mults = actual_mults + 1  ! Increment multiplication counter
          
          ! Check if the non-zero element already exists in the result matrix C
          if (nz == 0 .or. (C%row(nz) /= A%row(i) .or. C%col(nz) /= B%col(j))) then
            nz = nz + 1
            C%row(nz) = A%row(i)
            C%col(nz) = B%col(j)
            C%val(nz) = sum  ! Store the result
          else
            C%val(nz) = C%val(nz) + sum  ! Accumulate the result if the position already exists
          end if
        end if
      end do
    end do

    ! If there are more allocated slots than needed, resize the matrix to fit the actual number of non-zero elements
    if (nz < size(C%row)) then
      call resize_matrix(C, nz)
    end if

    ! Set the final number of non-zero elements in the result matrix
    C%non_zero = nz
  end subroutine multiply_matrices

end module sparse

