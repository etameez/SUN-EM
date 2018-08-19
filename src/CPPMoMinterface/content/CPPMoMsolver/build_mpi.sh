#!/bin/bash

cd build_mpi/
cmake -DMPI=ON -DATLAS=ON .. -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
make
time mpiexec -np 1 -genv OMP_NUM_THREADS=1 ./mom_mpi
time mpiexec -np 2 -genv OMP_NUM_THREADS=1 ./mom_mpi
time mpiexec -np 3 -genv OMP_NUM_THREADS=1 ./mom_mpi
time mpiexec -np 4 -genv OMP_NUM_THREADS=1 ./mom_mpi
time mpiexec -np 5 -genv OMP_NUM_THREADS=1 ./mom_mpi
time mpiexec -np 6 -genv OMP_NUM_THREADS=1 ./mom_mpi
time mpiexec -np 7 -genv OMP_NUM_THREADS=1 ./mom_mpi
time mpiexec -np 8 -genv OMP_NUM_THREADS=1 ./mom_mpi

