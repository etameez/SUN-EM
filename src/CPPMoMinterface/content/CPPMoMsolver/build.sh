#!/bin/bash

cd build/
cmake .. -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
make
time ./mom pec_plate.mom
