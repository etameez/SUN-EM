# Method of Moments C++

## Overview

This project will solve for the Method of Moments using C++.

The will be both a serial and parallel version. The parallel version will use both MPI and OpenMP.

Written By: Tameez Ebrahim  
Contact:    etameez@gmail.com  
Date:       Q3 2018  

## Dependencies

* C++11
* CMake Version 3.1 or higher
* An MPI implementaion (MPICH, IntelMPI, OpenMPI) 
* OpenMP

This was tested on Arch Linux 4.17.11 and macOS 10.12.6 using MPICH3 for the MPI implementation. Windows can only run the serial version.

## Setup

### Quick Setup
1. cd content/CPPMoMsolver/
2. mkdir build                      // For serial 
3. mkdir build_mpi                  // For parallel
4. ./build.sh or ./build_mpi.sh     // Serial and Parallel respectively  

### Self Compilation
1. cmake -DMPI=ON <path_to_build_directory> // For parallel  
   cmake <path_to_build_directory>          // For serial  
2. make                                     // Linux and macOS  
   nmake                                    // Windows if compiler is MSVC   
   mingw32-make                             // Windows if compiler is MinGW  
3. ./mom or ./mom_mpi                       // Linux and macOS  
   mom or mom_mpi                           // Windows  

If you are using Windows use either -G "MinGW Makefiles" or -G "NMake Makefiles" as an extra CMake compile command
e.g. cmake <path_to_build_directory> -G "MinGW Makefiles". Only MinGW was tested so YMMV with nmake.

