export OMP_NUM_THREADS=16
export MCAFLAGS="--mca btl_openib_warn_nonexistent_if 0 --mca btl_openib_warn_no_device_params_found 0"
make
echo
echo
for size in 1024 2048 4096 8192 16384 32768 65536
do
    echo
    mpiexec $MCAFLAGS -n 16 ./a.out $size
    echo '________________________________________________'
    echo
done
module purge
