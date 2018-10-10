#!/bin/bash

cd build_mpi/
cmake -DMPI=ON -DATLAS=ON .. -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
make
time mpiexec -np 4 -genv OMP_NUM_THREADS=1 ./mom_mpi ../../../../../examples/example-11/sphere.mom

