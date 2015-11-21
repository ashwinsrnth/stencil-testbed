default:
	mpif90 stencil.f90

openmp:
	mpif90 stencil.f90 -fopenmp
