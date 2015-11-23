module central8
contains

subroutine calcdfdx(cartcomm, nx, ny, f, dx, dfdx, &
        halo_left, halo_right, halo_bottom, halo_top, &
        halo_temp_lr, halo_temp_bt)
    use mpi
    use halo
    implicit none

        ! Compute derivative of function f (2-D array)
        ! using 8-th ordered central finite difference
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
    real(kind=8), intent(inout) :: halo_left(4,ny), &
        halo_right(4,ny), halo_bottom(4,nx), halo_top(4,nx)
    real(kind=8), intent(inout) :: halo_temp_lr(4,ny), &
        halo_temp_bt(4,nx)
    integer :: left_rank, right_rank
    integer :: ierr, istatus(MPI_STATUS_SIZE)
    integer :: i,j

    call MPI_Cart_shift(cartcomm, 1, 1, left_rank, right_rank, ierr)
    
    ! dfdx at inner points (this is the tight loop):

    !$omp parallel shared(dfdx,f,dx) private(i,j)
    !$omp do
    do j=1,ny
        do i=5,nx-4
            dfdx(i,j) = ((4.00d+00/5.00d+00)*(f(i+1,j)-f(i-1,j)) + &
                        (-1.00d+00/5.00d+00)*(f(i+2,j)-f(i-2,j)) + &
                        (4.00d+00/105.00d+00)*(f(i+3,j)-f(i-3,j)) + &
                        (-1.00d+00/280.00d+00)*(f(i+4,j)-f(i-4,j)))/dx
        end do
    end do
    !$omp end do
    !$omp end parallel

    call MPI_Barrier(cartcomm, ierr)

    ! Perform halo swaps
    call halo_swap(cartcomm, nx, ny, f, 4, &
        halo_left, halo_right, halo_bottom, halo_top, &
        halo_temp_lr, halo_temp_bt)

    ! Handle the boundaries

    ! left
    if (left_rank.eq.MPI_PROC_NULL) then
        do j=1,ny
            dfdx(1,j) = (-375150.00d+00*f(1,j) + &
                1565172.00d+00*f(2,j) - &
                4039140.00d+00*f(3,j) + &
                7180075.00d+00*f(4,j) - &
                8503950.00d+00*f(5,j) + &
                6676740.00d+00*f(6,j) - &
                3378844.00d+00*f(7,j) + &
                1024890.00d+00*f(8,j) - &
                157500.00d+00*f(9,j) + &
                8400.00d+00*f(10,j) - &
                840.00d+00*f(11,j) + &
                147.00d+00*f(12,j))/(105000.00d+00*dx)

            dfdx(2,j) = (56400.00d+00*f(1,j) - &
                725676.00d+00*f(2,j) + &
                2331945.00d+00*f(3,j) - &
                4322850.00d+00*f(4,j) + &
                5296900.00d+00*f(5,j) - &
                4232970.00d+00*f(6,j) + &
                2147502.00d+00*f(7,j) - &
                641320.00d+00*f(8,j) + &
                94500.00d+00*f(9,j) - &
                5250.00d+00*f(10,j) + &
                945.00d+00*f(11,j) - &
                126.00d+00*f(12,j))/(105000.00d+00*dx)

            dfdx(3,j) = (8186*f(1,j) - &
                66262.00d+00*f(2,j) + &
                181370.00d+00*f(3,j) - &
                336385.00d+00*f(4,j) + &
                407120.00d+00*f(5,j) - &
                289352.00d+00*f(6,j) + &
                114674.00d+00*f(7,j) - &
                18070.00d+00*f(8,j) - &
                1890.00d+00*f(9,j) + &
                630.00d+00*f(10,j) - &
                84.00d+00*f(11,j) + &
                63.00d+00*f(12,j))/(21000.00d+00*dx)

            dfdx(4,j) = (16480*f(1,j) - &
                121338*f(2,j) + &
                348810.00d+00*f(3,j) - &
                944475.00d+00*f(4,j) + &
                1234800.00d+00*f(5,j) - &
                792120.00d+00*f(6,j) + &
                323876.00d+00*f(7,j) - &
                77310.00d+00*f(8,j) + &
                16800.00d+00*f(9,j) - &
                7350.00d+00*f(10,j) + &
                1890.00d+00*f(11,j) - &
                63.00d+00*f(12,j))/(210000.00d+00*dx)
        end do
    else
        do j=1,ny
            i = 1
            dfdx(i,j) = ((4.00d+00/5.00d+00)*(f(i+1,j)-halo_left(1,j)) + &
                (-1.00d+00/5.00d+00)*(f(i+2,j)-halo_left(2,j)) + &
                (4.00d+00/105.00d+00)*(f(i+3,j)-halo_left(3,j)) + &
                (-1.00d+00/280.00d+00)*(f(i+4,j)-halo_left(4,j)))/dx

            i = 2
            dfdx(i,j) = ((4.00d+00/5.00d+00)*(f(i+1,j)-f(i-1,j)) + &
                (-1.00d+00/5.00d+00)*(f(i+2,j)-halo_left(1,j)) + &
                (4.00d+00/105.00d+00)*(f(i+3,j)-halo_left(2,j)) + &
                (-1.00d+00/280.00d+00)*(f(i+4,j)-halo_left(3,j)))/dx

            i = 3
            dfdx(i,j) = ((4.00d+00/5.00d+00)*(f(i+1,j)-f(i-1,j)) + &
                (-1.00d+00/5.00d+00)*(f(i+2,j)-f(i-2,j)) + &
                (4.00d+00/105.00d+00)*(f(i+3,j)-halo_left(1,j)) + &
                (-1.00d+00/280.00d+00)*(f(i+4,j)-halo_left(2,j)))/dx
            
            i = 4
            dfdx(i,j) = ((4.00d+00/5.00d+00)*(f(i+1,j)-f(i-1,j)) + &
                (-1.00d+00/5.00d+00)*(f(i+2,j)-f(i-2,j)) + &
                (4.00d+00/105.00d+00)*(f(i+3,j)-f(i-3,j)) + &
                (-1.00d+00/280.00d+00)*(f(i+4,j)-halo_left(1,j)))/dx
        end do
    end if

    if (right_rank.eq.MPI_PROC_NULL) then
        do j=1,ny
            dfdx(nx,j) = (-375150.00d+00*f(nx,j) + &
                1565172.00d+00*f(nx-1,j) - &
                4039140.00d+00*f(nx-2,j) + &
                7180075.00d+00*f(nx-3,j) - &
                8503950.00d+00*f(nx-4,j) + &
                6676740.00d+00*f(nx-5,j) - &
                3378844.00d+00*f(nx-6,j) + &
                1024890.00d+00*f(nx-7,j) - &
                157500.00d+00*f(nx-8,j) + &
                8400.00d+00*f(nx-9,j) - &
                840.00d+00*f(nx-10,j) + &
                147.00d+00*f(nx-11,j))/(-105000.00d+00*dx)

            dfdx(nx-1,j) = (56400.00d+00*f(nx,j) - &
                725676*f(nx-1,j) + &
                2331945.00d+00*f(nx-2,j) - &
                4322850.00d+00*f(nx-3,j) + &
                5296900.00d+00*f(nx-4,j) - &
                4232970.00d+00*f(nx-5,j) + &
                2147502.00d+00*f(nx-6,j) - &
                641320.00d+00*f(nx-7,j) + &
                94500.00d+00*f(nx-8,j) - &
                5250.00d+00*f(nx-9,j) + &
                945.00d+00*f(nx-10,j) - &
                126.00d+00*f(nx-11,j))/(-105000.00d+00*dx)

            dfdx(nx-2,j) = (8186.00d+00*f(nx,j) - &
                66262.00d+00*f(nx-1,j) + &
                181370.00d+00*f(nx-2,j) - &
                336385.00d+00*f(nx-3,j) + &
                407120.00d+00*f(nx-4,j) - &
                289352.00d+00*f(nx-5,j) + &
                114674.00d+00*f(nx-6,j) - &
                18070.00d+00*f(nx-7,j) - &
                1890.00d+00*f(nx-8,j) + &
                630.00d+00*f(nx-9,j) - &
                84.00d+00*f(nx-10,j) + &
                63.00d+00*f(nx-11,j))/(-21000.00d+00*dx)

            dfdx(nx-3,j) = (16480.00d+00*f(nx,j) - &
                121338.00d+00*f(nx-1,j) + &
                348810.00d+00*f(nx-2,j) - &
                944475.00d+00*f(nx-3,j) + &
                1234800.00d+00*f(nx-4,j) - &
                792120.00d+00*f(nx-5,j) + &
                323876.00d+00*f(nx-6,j) - &
                77310.00d+00*f(nx-7,j) + &
                16800.00d+00*f(nx-8,j) - &
                7350.00d+00*f(nx-9,j) + &
                1890.00d+00*f(nx-10,j) - &
                63.00d+00*f(nx-11,j))/(-210000.00d+00*dx)
        end do
    else
        do j=1,ny
            i = nx
            dfdx(i,j) = ((4.00d+00/5.00d+00)*(halo_right(1,j)-f(i-1,j)) + &
                (-1.00d+00/5.00d+00)*(halo_right(2,j)-f(i-2,j)) + &
                (4.00d+00/105.00d+00)*(halo_right(3,j)-f(i-3,j)) + &
                (-1.00d+00/280.00d+00)*(halo_right(4,j)-f(i-4,j)))/dx


            i = nx-1
            dfdx(i,j) = ((4.00d+00/5.00d+00)*(f(i+1,j)-f(i-1,j)) + &
                (-1.00d+00/5.00d+00)*(halo_right(1,j)-f(i-2,j)) + &
                (4.00d+00/105.00d+00)*(halo_right(2,j)-f(i-3,j)) + &
                (-1.00d+00/280.00d+00)*(halo_right(3,j)-f(i-4,j)))/dx

            i = nx-2
            dfdx(i,j) = ((4.00d+00/5.00d+00)*(f(i+1,j)-f(i-1,j)) + &
                (-1.00d+00/5.00d+00)*(f(i+2,j)-f(i-2,j)) + &
                (4.00d+00/105.00d+00)*(halo_right(1,j)-f(i-3,j)) + &
                (-1.00d+00/280.00d+00)*(halo_right(2,j)-f(i-4,j)))/dx

            i = nx-3
            dfdx(i,j) = ((4.00d+00/5.00d+00)*(f(i+1,j)-f(i-1,j)) + &
                (-1.00d+00/5.00d+00)*(f(i+2,j)-f(i-2,j)) + &
                (4.00d+00/105.00d+00)*(f(i+3,j)-f(i-3,j)) + &
                (-1.00d+00/280.00d+00)*(halo_right(1,j)-f(i-4,j)))/dx
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
    do j=2,ny-1
        do i=1,nx
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
