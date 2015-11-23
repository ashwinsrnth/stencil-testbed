
program stencil
    use omp_lib
    use mpi
    use central8

    implicit none

    integer :: ierr, rank, nprocs, nthreads, threadID, cartcomm
    integer :: left_rank, right_rank, bottom_rank, top_rank
    integer :: dims(2), coords(2)
    logical :: periods(2), reorder
    integer(kind=8) :: N, N_global
    real(kind=8) :: dx
    real(kind=8), allocatable :: x(:,:), y(:,:), f(:,:), dfdx(:,:)
    real(kind=8), allocatable :: halo_left(:,:), halo_right(:,:), halo_bottom(:,:), halo_top(:,:)
    real(kind=8), allocatable :: halo_temp_lr(:,:), halo_temp_bt(:,:)
    real(kind=8) :: error, global_error
    integer :: i, j, step
    real(kind=8) :: t1, t2
    character(32) :: N_global_string
    
    call MPI_Init(ierr)
    call MPI_Comm_size(MPI_comm_world, nprocs, ierr)

    call getarg(1, N_global_string)
    read(N_global_string, *), N_global
    
    dims(1) = int(sqrt(real(nprocs)))
    dims(2) = int(sqrt(real(nprocs)))
    periods(1) = .false.
    periods(2) = .false.
    reorder = .false.
    N = N_global/dims(2)
    call MPI_Cart_create(MPI_Comm_world, 2, dims, periods, reorder, cartcomm, ierr)
    call MPI_Comm_rank(cartcomm, rank, ierr)
    call MPI_Cart_shift(cartcomm, 1, 1, left_rank, right_rank, ierr)
    call MPI_Cart_shift(cartcomm, 0, 1, bottom_rank, top_rank, ierr)
    call MPI_Cart_get(cartcomm, 2, dims, periods, coords, ierr)
    
    allocate(f(N,N))
    allocate(x(N,N))
    allocate(y(N,N))
    allocate(dfdx(N,N))
    allocate(halo_left(4,N))
    allocate(halo_right(4,N))
    allocate(halo_bottom(4,N))
    allocate(halo_top(4,N))
    allocate(halo_temp_lr(4,N))
    allocate(halo_temp_bt(4,N))

    ! Initialize function values
    dx = 1.00d+00/(N_global-1)
    do j=1,N
        do i=1,N
            x(i,j) = (coords(2)*N+i-1.00d+00)*dx
            y(i,j) = (coords(1)*N+j-1.00d+00)*dx
            f(i,j) = dsin(x(i,j))
        end do
    end do
    
    ! Evaluate derivative (and time the evaluation)
    call MPI_Barrier(cartcomm, ierr)
    t1 = MPI_Wtime()
    
    do step=1,50
        call calcdfdx(cartcomm, N, N, f, dx, dfdx, &
        halo_left, halo_right, halo_bottom, halo_top, &
        halo_temp_lr, halo_temp_bt)
    end do

    call MPI_Barrier(cartcomm, ierr)
    t2 = MPI_Wtime()

    ! Compute error
    error = 0.0d+00
    do j=1,N
        do i=1,N
            error = error + dabs(dfdx(i,j) - dcos(x(i,j)))
        end do
    end do

    call MPI_Reduce(error, global_error, 1, MPI_DOUBLE_PRECISION, &
        MPI_SUM, 0, cartcomm)
    
     
    ! Write report
    if (rank.eq.0) then
        print '(A33, I11)', "Problem size per direction: ", N_global
        !print '(A33, E11.4)', "Global error: ", (global_error/(N_global*N_global))
        print '(A33, I11)', "Number of MPI processes: ", nprocs

        !$omp parallel private(threadID)
        threadID = omp_get_thread_num()
        nthreads = omp_get_num_threads()
        if (threadID.eq.0) then
            print '(A33, I11)', "OMP threads per process: ", nthreads
        end if
        !$omp end parallel
        print '(A33, F9.4, A2)', "Time for dfdx: ", (t2-t1)/50.0d+00, "s"
        write (*,*), global_error/(N_global*N_global)
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

    call MPI_Finalize(ierr)
end program
