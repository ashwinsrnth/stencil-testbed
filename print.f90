module print

contains

subroutine print_array2d(f, nx, ny)
    implicit none
    real(kind=8), intent(in) :: f(nx, ny)
    integer(kind=8) :: nx, ny
    integer i, j
    do j=1,ny
        do i=1,nx
            write (*,"(f6.2)",advance="no") f(i,j)
        end do
        write (*,*)
    end do
end subroutine

end module 
