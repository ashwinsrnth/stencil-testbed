
subroutine calcdfdx(cartcomm, nx, ny, f, dx, dfdx, &
        halo_left, halo_right, halo_bottom, halo_top, &
        halo_temp_lr, halo_temp_bt)
    use mpi
    use halo
    implicit none
    
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
    call halo_swap(cartcomm, nx, ny, f, &
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
    call halo_swap(cartcomm, nx, ny, f, &
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




program stencil
    use omp_lib
    use mpi
    implicit none

    integer :: ierr, rank, nprocs, nthreads, threadID, cartcomm
    integer :: left_rank, right_rank, bottom_rank, top_rank
    integer :: dims(2), coords(2)
    logical :: periods(2), reorder
    integer(kind=8) :: N, N_global
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
    N_global = 1024
    N = N_global/dims(2)
    periods(1) = .false.
    periods(2) = .false.
    reorder = .false.
    call MPI_Cart_create(MPI_Comm_world, 2, dims, periods, reorder, cartcomm, ierr)
    call MPI_Comm_rank(cartcomm, rank, ierr)
    call MPI_Cart_shift(cartcomm, 1, 1, left_rank, right_rank, ierr)
    call MPI_Cart_shift(cartcomm, 0, 1, bottom_rank, top_rank, ierr)
    call MPI_Cart_get(cartcomm, 2, dims, periods, coords, ierr)
    
    if (rank.eq.0) then
        write (*,*), "--------------------------------------------------------"
        write (*,*), "Number of MPI processes: ", nprocs
        write (*,*), "Array size per dimension per process: ", N
    end if
    
    if (rank.eq.0) then    
        !$omp parallel private(threadID)
        threadID = omp_get_thread_num()
        nthreads = omp_get_num_threads()
        if (threadID.eq.0) then
            write (*,*), "Number of OpenMP threads per process: ", nthreads
        end if
        !$omp end parallel
    end if

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
    
    do step=1,5
        call calcdfdx(cartcomm, N, N, f, dx, dfdx, &
        halo_left, halo_right, halo_bottom, halo_top, &
        halo_temp_lr, halo_temp_bt)
    end do

    call MPI_Barrier(cartcomm, ierr)
    t2 = MPI_Wtime()

    if (rank.eq.0) then
        write (*,*) "Time for dfdx: ", t2-t1
    end if

    call MPI_Barrier(cartcomm, ierr)
    t1 = MPI_Wtime()
    
    do step=1,5
        call calcdfdy(cartcomm, N, N, f, dx, dfdy, &
        halo_left, halo_right, halo_bottom, halo_top, &
        halo_temp_lr, halo_temp_bt)
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
        write (*,*), "--------------------------------------------------------"
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
