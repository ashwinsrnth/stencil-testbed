
export OMP_NUM_THREADS=1
export MCAFLAGS="--mca btl_openib_warn_nonexistent_if 0 --mca btl_openib_warn_no_device_params_found 0"
echo
echo
make
for size in 1024 2048 4096 8192 16384 32768
do
    echo mpiexec $MCAFLAGS -n 256 ./a.out $size
    mpiexec $MCAFLAGS -n 256 ./a.out $size
    echo '______________________________'
    echo
    echo
done
module purge
