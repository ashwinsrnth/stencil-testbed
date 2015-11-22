default:
	mpif90 -c halo.f90 -O2 -fopenmp 
	mpif90 -c stencil.f90 -O2 -fopenmp
	mpif90 -o a.out stencil.o halo.o -fopenmp
	rm *.o *.mod 
