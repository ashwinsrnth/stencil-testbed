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
    integer :: rank, left_rank, right_rank, bottom_rank, top_rank
    integer :: ierr, istatus(MPI_STATUS_SIZE)
    integer :: i,j

    call MPI_Comm_rank(cartcomm, rank, ierr)
    call MPI_Cart_shift(cartcomm, 1, 1, left_rank, right_rank, ierr)
    call MPI_Cart_shift(cartcomm, 0, 1, bottom_rank, top_rank, ierr)

    ! copy from function to halos
    do i = 1,ny
        halo_left(1,i) = f(1,i)
        halo_temp_lr(1,i) = f(nx,i)
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
    integer :: left_rank, right_rank, bottom_rank, top_rank
    integer :: dims(2), coords(2)
    logical :: periods(2), reorder
    integer, parameter :: N = 256
    double precision :: dx
    double precision :: x(N, N), y(N, N), f(N, N), dfdx(N, N), dfdy(N, N)
    double precision :: halo_left(1, N), halo_right(1, N), halo_bottom(N, 1), halo_top(N, 1)
    double precision :: halo_temp_lr(1, N), halo_temp_bt(N, 1)
    double precision :: error = 0
    integer :: i, j

    call MPI_Init(ierr)

    dims(1) = 3
    dims(2) = 3
    periods(1) = .false.
    periods(2) = .false.
    reorder = .false.
    call MPI_Cart_create(MPI_Comm_world, 2, dims, periods, reorder, cartcomm, ierr)
    call MPI_Comm_rank(cartcomm, rank, ierr)
    call MPI_Cart_shift(cartcomm, 1, 1, left_rank, right_rank, ierr)
    call MPI_Cart_shift(cartcomm, 0, 1, bottom_rank, top_rank, ierr)
    call MPI_Cart_get(cartcomm, 2, dims, periods, coords, ierr)

    ! Initialize array

    dx = 1./(dims(2)*N-1)
    do j=1,N
        do i=1,N
            x(i,j) = (coords(2)*N + i - 1)*dx
            y(i,j) = (coords(1)*N + j - 1)*dx
            f(i,j) = sin(x(i,j))
        end do
    end do

    call halo_swap(cartcomm, x, halo_left, halo_right, halo_bottom, halo_top, &
             halo_temp_lr, halo_temp_bt, N, N)

    ! Calculate derivatives in the interior
    do j=1,N
        do i=2,N-1
            dfdx(i,j) = (f(i+1,j) - f(i-1,j))/(2*dx)
        end do
    end do

    call MPI_Barrier(cartcomm, ierr)

    ! Perform halo swaps
    call halo_swap(cartcomm, f, halo_left, halo_right, halo_bottom, halo_top, &
             halo_temp_lr, halo_temp_bt, N, N)
   
    ! Handle the boundaries

    ! left
    if (left_rank.eq.MPI_PROC_NULL) then
        do i=1,N
            dfdx(1,i) = (f(2,i) - f(1,i))/dx
        end do
    else
        do i=1,N
            dfdx(1,i) = (f(2,i) - halo_left(1,i))/(2*dx)
        end do
    end if

    ! right
    if (right_rank.eq.MPI_PROC_NULL) then
        do i=1,N
            dfdx(N,i) = (f(N,i) - f(N-1,i))/dx
        end do
    else
        do i=1,N
            dfdx(N,i) = (halo_right(1,i) - f(N-1,i))/(2*dx)
        end do
    end if
    
    call MPI_Barrier(cartcomm, ierr)
    ! bottom
    if (bottom_rank.eq.MPI_PROC_NULL) then
        do i=1,N
            dfdy(i,1) = (f(i,2) - f(i,1))/dx
        end do
    else
        do i=1,N
            dfdy(i,1) = (f(i,2) - halo_bottom(i,1))/(2*dx)
        end do
    end if
    
    ! top
    if (top_rank.eq.MPI_PROC_NULL) then
        do i=1,N
            dfdy(i,N) = (f(i,N) - f(i,N-1))/dx
        end do
    else
        do i=1,N
            dfdy(i,N) = (halo_top(i,1) - f(i,N-1))/(2*dx)
        end do
    end if

    ! compute error
    if (rank.eq.4) then
    do j=1,N
        do i=1,N
            error = error + dabs(dfdx(i,j) - cos(x(i,j)))
        end do
    end do
    end if

    if (rank.eq.4) then
        write (*,*), error/(N*N)
    end if

    call MPI_Finalize(ierr)
end program
