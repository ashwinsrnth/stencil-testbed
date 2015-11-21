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
    real(kind=8), intent(inout) :: f(nx,ny), halo_left(ny,1), halo_right(ny,1), halo_bottom(nx,1), halo_top(nx,1)
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

subroutine print_array2d(f, nx, ny)
    implicit none
    real(kind=8), intent(in) :: f(nx, ny)
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

    integer :: ierr, rank, nprocs, cartcomm
    integer :: left_rank, right_rank, bottom_rank, top_rank
    integer :: dims(2), coords(2)
    logical :: periods(2), reorder
    integer :: N, N_global
    real(kind=8) :: dx
    real(kind=8), allocatable :: x(:,:), y(:,:), f(:,:), dfdx(:,:), dfdy(:,:)
    real(kind=8), allocatable :: halo_left(:,:), halo_right(:,:), halo_bottom(:,:), halo_top(:,:)
    real(kind=8), allocatable :: halo_temp_lr(:,:), halo_temp_bt(:,:)
    real(kind=8) :: error, global_error
    integer :: i, j, step
    real(kind=8) :: t1, t2

    call MPI_Init(ierr)
    call MPI_Comm_size(MPI_comm_world, nprocs, ierr)
    dims(1) = int(sqrt(real(nprocs)))
    dims(2) = int(sqrt(real(nprocs)))
    N_global = 8192
    N = N_global/dims(2)
    periods(1) = .false.
    periods(2) = .false.
    reorder = .false.
    call MPI_Cart_create(MPI_Comm_world, 2, dims, periods, reorder, cartcomm, ierr)
    call MPI_Comm_rank(cartcomm, rank, ierr)
    call MPI_Cart_shift(cartcomm, 1, 1, left_rank, right_rank, ierr)
    call MPI_Cart_shift(cartcomm, 0, 1, bottom_rank, top_rank, ierr)
    call MPI_Cart_get(cartcomm, 2, dims, periods, coords, ierr)

    allocate(f(N,N))
    allocate(x(N,N))
    allocate(y(N,N))
    allocate(dfdx(N,N))
    allocate(dfdy(N,N))
    allocate(halo_left(1,N))
    allocate(halo_right(1,N))
    allocate(halo_bottom(N,1))
    allocate(halo_top(N,1))
    allocate(halo_temp_lr(1,N))
    allocate(halo_temp_bt(N,1))

    ! Initialize array
    dx = 1./(N_global-1)
    do j=1,N
        do i=1,N
            x(i,j) = (coords(2)*N + i - 1)*dx
            y(i,j) = (coords(1)*N + j - 1)*dx
            f(i,j) = sin(x(i,j)) + sin(y(i,j))
        end do
    end do
    
    call MPI_Barrier(cartcomm, ierr)
    t1 = MPI_Wtime()
    
    do step=1,50
        ! Calculate derivatives in the interior
        
        !$omp parallel shared(dfdx,f,dx) private(i,j)
        !$omp do
        do j=1,N
            do i=2,N-1
                dfdx(i,j) = (f(i+1,j) - f(i-1,j))/(2*dx)
            end do
        end do
        !$omp end do
        !$omp end parallel

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
                dfdx(1,i) = (f(2,i) - halo_left(i,1))/(2*dx)
            end do
        end if

        ! right
        if (right_rank.eq.MPI_PROC_NULL) then
            do i=1,N
                dfdx(N,i) = (f(N,i) - f(N-1,i))/dx
            end do
        else
            do i=1,N
                dfdx(N,i) = (halo_right(i,1) - f(N-1,i))/(2*dx)
            end do
        end if
    end do

    call MPI_Barrier(cartcomm, ierr)
    t2 = MPI_Wtime()

    if (rank.eq.0) then
        write (*,*) "Time for dfdx: ", t2-t1
    end if

    call MPI_Barrier(cartcomm, ierr)
    t1 = MPI_Wtime()
    
    do step=1,50

        !$omp parallel shared(f, dfdy, dx) private(i,j)
        !$omp do
        do j=2,N-1
            do i=1,N
                dfdy(i,j) = (f(i,j+1) - f(i,j-1))/(2*dx)
            end do
        end do
        !$omp end do
        !$omp end parallel

        call MPI_Barrier(cartcomm, ierr)

        ! Perform halo swaps
        call halo_swap(cartcomm, f, halo_left, halo_right, halo_bottom, halo_top, &
                 halo_temp_lr, halo_temp_bt, N, N)
        
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
    end do

    call MPI_Barrier(cartcomm, ierr)
    t2 = MPI_Wtime()

    if (rank.eq.0) then
        write (*,*) "Time for dfdy: ", t2-t1
    end if

    ! compute error
    error = 0.0
    do j=1,N
        do i=1,N
            error = error + dabs(dfdx(i,j) - cos(x(i,j)))
        end do
    end do

    call MPI_Reduce(error, global_error, 1, MPI_DOUBLE_PRECISION, &
        MPI_SUM, 0, cartcomm)

    if (rank.eq.0) then
        write (*,*), "err in dfdx: ", global_error/(N_global*N_global)
    end if

    error = 0.0
    do j=1,N
        do i=1,N
            error = error + dabs(dfdy(i,j) - (cos(y(i,j))))
        end do
    end do

    call MPI_Reduce(error, global_error, 1, MPI_DOUBLE_PRECISION, &
        MPI_SUM, 0, cartcomm)

    if (rank.eq.0) then
        write (*,*), "err in dfdy: ", global_error/(N_global*N_global)
    end if

    deallocate(halo_left)
    deallocate(halo_right)
    deallocate(halo_bottom)
    deallocate(halo_top)
    deallocate(halo_temp_lr)
    deallocate(halo_temp_bt)
    deallocate(f)
    deallocate(x)
    deallocate(y)
    deallocate(dfdx)
    deallocate(dfdy)

    call MPI_Finalize(ierr)
end program
