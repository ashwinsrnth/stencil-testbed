subroutine halo_swap(cartcomm, f, halo_left, halo_right, halo_bottom, halo_top, &
                halo_temp_lr, halo_temp_bt, nx, ny)
        ! Swaps the left-right and bottom-top halos
        ! for all processes
        !
        ! Parameters:
        !   cartcomm:   2-D Cartesian communicator
        !   f:  2-D array from which the halos take values
        !   halo_{left, right}: Space for left and right halos [sized (1,N)]
        !   halo_{bottom, top}: Space for bottom and top halos [sized (N,1)]
        !   halo_temp_{lr, bt}: Scratch space sized (1,N) and (N,1) respectively
        !   nx, ny: Size of array f
    use mpi
    implicit none

    integer, intent(in) :: cartcomm, nx, ny
    double precision, intent(inout) :: f(nx,ny), halo_left(1,ny), halo_right(1,ny), halo_bottom(nx,1), halo_top(nx,1)
    double precision, intent(inout) :: halo_temp_lr(1,ny), halo_temp_bt(nx,1)
    integer :: rank, left_rank, right_rank, bottom_rank, top_rank, ierr, i, istatus(MPI_STATUS_SIZE)
    call MPI_Comm_rank(cartcomm, rank, ierr)
    call MPI_Cart_shift(cartcomm, 1, 1, left_rank, right_rank, ierr)
    call MPI_Cart_shift(cartcomm, 0, 1, bottom_rank, top_rank, ierr)

    ! Send left halo and receive on right
    do i = 1,ny
        halo_left(1,i) = f(1,i)
        halo_temp_lr(1,i) = f(nx,i)
    end do

    call MPI_Sendrecv(halo_left, ny, MPI_DOUBLE_PRECISION, left_rank, &
        10, &
        halo_right, ny, MPI_DOUBLE_PRECISION, right_rank, &
        10, &
        cartcomm, istatus, ierr)

    ! Send right halo and receive on left
    call MPI_Sendrecv(halo_temp_lr, ny, MPI_DOUBLE_PRECISION, right_rank, &
        20, &
        halo_left, ny, MPI_DOUBLE_PRECISION, left_rank, &
        20, &
        cartcomm, istatus, ierr)

    ! Send bottom halo and receive on top
    do i = 1,nx
        halo_bottom(i,1) = f(i,1)
        halo_temp_bt(1,i) = f(i,ny)
    end do
    
    call MPI_Sendrecv(halo_bottom, nx, MPI_DOUBLE_PRECISION, bottom_rank, &
        30, &
        halo_top, nx, MPI_DOUBLE_PRECISION, top_rank, &
        30, &
        cartcomm, istatus, ierr)

    ! Send top halo and receive on bottom
    call MPI_Sendrecv(halo_temp_bt, nx, MPI_DOUBLE_PRECISION, top_rank, &
        40, &
        halo_bottom, nx, MPI_DOUBLE_PRECISION, bottom_rank, &
        40, &
        cartcomm, istatus, ierr)
end subroutine

subroutine print_array2d(f, nx, ny)
    implicit none
    double precision, intent(in) :: f(nx, ny)
    integer, intent(in) :: nx, ny
    integer i, j
    do j=1,ny
        do i=1,nx
            write (*,"(f6.2)",advance="no") f(i,j)
        end do
        write (*,*)
    end do
end subroutine

program stencil
    use mpi
    implicit none

    integer :: ierr, rank, cartcomm
    integer :: dims(2)
    logical :: periods(2), reorder
    integer, parameter :: N = 5
    double precision :: f(N, N)
    double precision :: halo_left(1, N), halo_right(1, N), halo_bottom(N, 1), halo_top(N, 1)
    double precision :: halo_temp_lr(1, N), halo_temp_bt(N, 1)
    integer :: i, j

    call MPI_Init(ierr)
    call MPI_comm_rank(MPI_Comm_world, rank, ierr)

    dims(1) = 3
    dims(2) = 3
    periods(1) = .false.
    periods(2) = .false.
    reorder = .false.
    call MPI_Cart_create(MPI_Comm_world, 2, dims, periods, reorder, cartcomm, ierr)

    ! Initialize array
    do j=1, N
        do i=1,N
            f(i,j) = rank 
        end do
    end do

    ! Perform halo swaps
    call halo_swap(cartcomm, f, halo_left, halo_right, halo_bottom, halo_top, &
             halo_temp_lr, halo_temp_bt, N, N)
    
    if (rank.eq.4) then
        call print_array2d(halo_left, 1, N)
        call print_array2d(halo_right, 1, N)
        call print_array2d(halo_bottom, N, 1)
        call print_array2d(halo_top, N, 1)
    end if

    call MPI_Finalize(ierr)
end program
