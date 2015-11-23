module central2
contains

subroutine calcdfdx(cartcomm, nx, ny, f, dx, dfdx, &
        halo_left, halo_right, halo_bottom, halo_top, &
        halo_temp_lr, halo_temp_bt)
    use mpi
    use halo
    implicit none

        ! Compute derivative of function f (2-D array)
        ! using 2th order central finite difference
        !
        ! cartcomm: Cartesian communicator
        ! nx, ny: Size of array
        ! f: Array of function values
        ! dx: Spacing between grid points
        ! dfdx: Space for output
        ! halo*: Space for halos
    
    integer, intent(in) :: cartcomm
    integer(kind=8), intent(in) :: nx, ny
    real(kind=8), intent(in) :: f(nx,ny), dx
    real(kind=8), intent(inout) :: dfdx(nx,ny)
    real(kind=8), intent(inout) :: halo_left(ny,1), &
        halo_right(ny,1), halo_bottom(nx,1), halo_top(nx,1)
    real(kind=8), intent(inout) :: halo_temp_lr(ny,1), &
        halo_temp_bt(nx,1)
    integer :: left_rank, right_rank
    integer :: ierr, istatus(MPI_STATUS_SIZE)
    integer :: i,j
    

    call MPI_Cart_shift(cartcomm, 1, 1, left_rank, right_rank, ierr)

    !$omp parallel shared(dfdx,f,dx) private(i,j)
    !$omp do
    do j=1,nx
        do i=2,ny-1
            dfdx(i,j) = (f(i+1,j) - f(i-1,j))/(2*dx)
        end do
    end do
    !$omp end do
    !$omp end parallel

    call MPI_Barrier(cartcomm, ierr)

    ! Perform halo swaps
    call halo_swap(cartcomm, nx, ny, f, 1, &
        halo_left, halo_right, halo_bottom, halo_top, &
        halo_temp_lr, halo_temp_bt)
   
    ! Handle the boundaries

    ! left
    if (left_rank.eq.MPI_PROC_NULL) then
        do i=1,ny
            dfdx(1,i) = (f(2,i) - f(1,i))/dx
        end do
    else
        do i=1,ny
            dfdx(1,i) = (f(2,i) - halo_left(i,1))/(2*dx)
        end do
    end if

    ! right
    if (right_rank.eq.MPI_PROC_NULL) then
        do i=1,ny
            dfdx(nx,i) = (f(nx,i) - f(nx-1,i))/dx
        end do
    else
        do i=1,ny
            dfdx(ny,i) = (halo_right(i,1) - f(ny-1,i))/(2*dx)
        end do
    end if
end subroutine

subroutine calcdfdy(cartcomm, nx, ny, f, dy, dfdy, &
        halo_left, halo_right, halo_bottom, halo_top, &
        halo_temp_lr, halo_temp_bt)
    use mpi
    use halo
    implicit none
    
    integer, intent(in) :: cartcomm
    integer(kind=8), intent(in) :: nx, ny
    real(kind=8), intent(in) :: f(nx,ny), dy
    real(kind=8), intent(inout) :: dfdy(nx,ny)
    real(kind=8), intent(inout) :: halo_left(ny,1), &
        halo_right(ny,1), halo_bottom(nx,1), halo_top(nx,1)
    real(kind=8), intent(inout) :: halo_temp_lr(ny,1), &
        halo_temp_bt(nx,1)
    integer :: left_rank, right_rank
    integer :: ierr, istatus(MPI_STATUS_SIZE)
    integer :: i,j
    integer :: bottom_rank, top_rank

    call MPI_Cart_shift(cartcomm, 0, 1, bottom_rank, top_rank, ierr)

    !$omp parallel shared(f, dfdy, dy) private(i,j)
    !$omp do
    do j=2,nx-1
        do i=1,ny
            dfdy(i,j) = (f(i,j+1) - f(i,j-1))/(2*dy)
        end do
    end do
    !$omp end do
    !$omp end parallel

    call MPI_Barrier(cartcomm, ierr)

    ! Perform halo swaps
    call halo_swap(cartcomm, nx, ny, f, 1, &
        halo_left, halo_right, halo_bottom, halo_top, &
        halo_temp_lr, halo_temp_bt)
    
    ! bottom
    if (bottom_rank.eq.MPI_PROC_NULL) then
        do i=1,nx
            dfdy(i,1) = (f(i,2) - f(i,1))/dy
        end do
    else
        do i=1,nx
            dfdy(i,1) = (f(i,2) - halo_bottom(i,1))/(2*dy)
        end do
    end if
    
    ! top
    if (top_rank.eq.MPI_PROC_NULL) then
        do i=1,nx
            dfdy(i,ny) = (f(i,ny) - f(i,ny-1))/dy
        end do
    else
        do i=1,nx
            dfdy(i,ny) = (halo_top(i,1) - f(i,ny-1))/(2*dy)
        end do
    end if
end subroutine

end module 
