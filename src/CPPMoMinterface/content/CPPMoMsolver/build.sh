#!/bin/bash

cd build/
cmake .. -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
make
time ./mom ../../../../../examples/example-10/pec_plate.mom
