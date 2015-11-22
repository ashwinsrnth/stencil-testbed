module halo

contains

subroutine halo_swap(cartcomm, nx, ny, f, stencil_width, &
                halo_left, halo_right, halo_bottom, halo_top, &
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
    integer, intent(in) :: stencil_width
    real(kind=8), intent(inout) :: halo_left(ny,1), halo_right(ny,1), halo_bottom(nx,1), halo_top(nx,1)
    real(kind=8), intent(inout) :: halo_temp_lr(ny,1), halo_temp_bt(nx,1)
    integer :: rank, left_rank, right_rank, bottom_rank, top_rank
    integer :: ierr, istatus(MPI_STATUS_SIZE)
    integer :: i,j

    call MPI_Comm_rank(cartcomm, rank, ierr)
    call MPI_Cart_shift(cartcomm, 1, 1, left_rank, right_rank, ierr)
    call MPI_Cart_shift(cartcomm, 0, 1, bottom_rank, top_rank, ierr)

    ! copy from function to halos
    do i=1,ny
        do j=1,stencil_width
            halo_left(i,j) = f(j,i)
            halo_temp_lr(i,j) = f(nx-(j-1),i)
        end do
    end do

    do i=1,nx
        do j=1,stencil_width
            halo_bottom(i,j) = f(i,j)
            halo_temp_bt(i,j) = f(i,ny-(j-1))
        end do
    end do
    
    call MPI_Barrier(cartcomm, ierr)

    ! swap halos
    call MPI_Sendrecv(halo_left, stencil_width*ny, MPI_DOUBLE_PRECISION, left_rank, &
        10, &
        halo_right, stencil_width*ny, MPI_DOUBLE_PRECISION, right_rank, &
        10, &
        cartcomm, istatus, ierr)

    call MPI_Sendrecv(halo_temp_lr, stencil_width*ny, MPI_DOUBLE_PRECISION, right_rank, &
        20, &
        halo_left, stencil_width*ny, MPI_DOUBLE_PRECISION, left_rank, &
        20, &
        cartcomm, istatus, ierr)
    
    call MPI_Sendrecv(halo_bottom, stencil_width*nx, MPI_DOUBLE_PRECISION, bottom_rank, &
        30, &
        halo_top, stencil_width*nx, MPI_DOUBLE_PRECISION, top_rank, &
        30, &
        cartcomm, istatus, ierr)
    
    call MPI_Sendrecv(halo_temp_bt, stencil_width*nx, MPI_DOUBLE_PRECISION, top_rank, &
        40, &
        halo_bottom, stencil_width*nx, MPI_DOUBLE_PRECISION, bottom_rank, &
        40, &
        cartcomm, istatus, ierr)
end subroutine

end module
