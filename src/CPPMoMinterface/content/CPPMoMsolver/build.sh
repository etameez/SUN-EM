#!/bin/bash

cd ~/git/SUN-EM/src/CPPMoMinterface/content/CPPMoMsolver/.build/
cmake .. -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
make
./mom
