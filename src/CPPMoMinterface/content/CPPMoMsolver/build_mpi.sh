#!/bin/bash

cd build_mpi/
cmake -DMPI=ON .. -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
make
mpiexec -np 4 -oversubscribe --bind-to none ./mom_mpi

