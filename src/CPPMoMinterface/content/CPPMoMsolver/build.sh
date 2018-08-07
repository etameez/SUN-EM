#!/bin/bash

cd .build/
cmake .. -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
make
mpiexec -np 4 ./mom
