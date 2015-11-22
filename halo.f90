module halo

contains

subroutine halo_swap(cartcomm, nx, ny, f, halo_left, halo_right, halo_bottom, halo_top, &
                halo_temp_lr, halo_temp_bt)
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

    integer, intent(in) :: cartcomm
    integer(kind=8), intent(in) :: nx, ny
    real(kind=8), intent(in) :: f(nx,ny)
    real(kind=8), intent(inout) :: halo_left(ny,1), halo_right(ny,1), halo_bottom(nx,1), halo_top(nx,1)
    real(kind=8), intent(inout) :: halo_temp_lr(ny,1), halo_temp_bt(nx,1)
    integer :: rank, left_rank, right_rank, bottom_rank, top_rank
    integer :: ierr, istatus(MPI_STATUS_SIZE)
    integer :: i,j

    call MPI_Comm_rank(cartcomm, rank, ierr)
    call MPI_Cart_shift(cartcomm, 1, 1, left_rank, right_rank, ierr)
    call MPI_Cart_shift(cartcomm, 0, 1, bottom_rank, top_rank, ierr)

    ! copy from function to halos
    do i = 1,ny
        halo_left(i,1) = f(1,i)
        halo_temp_lr(i,1) = f(nx,i)
    end do

    do i = 1,nx
        halo_bottom(i,1) = f(i,1)
        halo_temp_bt(i,1) = f(i,ny)
    end do
    
    call MPI_Barrier(cartcomm, ierr)

    ! swap halos
    call MPI_Sendrecv(halo_left, ny, MPI_DOUBLE_PRECISION, left_rank, &
        10, &
        halo_right, ny, MPI_DOUBLE_PRECISION, right_rank, &
        10, &
        cartcomm, istatus, ierr)

    call MPI_Sendrecv(halo_temp_lr, ny, MPI_DOUBLE_PRECISION, right_rank, &
        20, &
        halo_left, ny, MPI_DOUBLE_PRECISION, left_rank, &
        20, &
        cartcomm, istatus, ierr)
    
    call MPI_Sendrecv(halo_bottom, nx, MPI_DOUBLE_PRECISION, bottom_rank, &
        30, &
        halo_top, nx, MPI_DOUBLE_PRECISION, top_rank, &
        30, &
        cartcomm, istatus, ierr)
    
    call MPI_Sendrecv(halo_temp_bt, nx, MPI_DOUBLE_PRECISION, top_rank, &
        40, &
        halo_bottom, nx, MPI_DOUBLE_PRECISION, bottom_rank, &
        40, &
        cartcomm, istatus, ierr)
end subroutine

end module
