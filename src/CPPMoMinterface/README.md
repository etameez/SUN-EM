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
* Scalapack

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

## In Depth Setup For Parallel With Scalapack

There are 3 options  
1. Using the refernce libraries
2. Using ATLAS
3. Using ATLAS threaded libraries explicitly

### Using the reference libraries
The libraries are located in lib/scalapack
Have  to have libgfortran in /usr/lib/

Either
1. mkdir build_mpi
2. ./build_mpi.sh

Or
1. mkdir folder
2. cd folder
3. cmake -DMPI=ON ..
4. make
5. mpiexec -np num_procs ./mom_mpi

### Using ATLAS
1. Disable hyperthreading in bios
2. Install CpuPower (pacman -S cpupower)
3. cpupower frequency-info
4. If you see intel_pstate then,
 	Edit /etc/default/grub and append GRUB_CMDLINE_LINUX_DEFAULT = "intel_pstate=disable"
 	Rebuild grub config
5. cpupower frequency-set -g performance
6. Install ATLAS (yaourt -S atlas-lapack)
7. Make sure the libaries are in or symlinked in /usr/bin
8. Download Scalapack from the scalapack website
8. Unzip it
9. Copy the resulting folder to /usr/bin/ and rename the folder to scalapack
10, In folder scalapack
11. mkdir build
12. cd build
13. cmake ..
14. make
15. make test 
16. If the tests pass, scalapack is installed properly
17. Go back to the project directory
18. mkdir folder
19. cd folder
20. cmake -DMPI=ON -DATLAS=ON ..
21. make
22. mpiexec -np num_procs ./mom_mpi 

###Using ATLAS threaded libraries explicitly
1. Do the same as above, but stop at step 20
2. cmake -DMPI=ON -DPARALLEL=ON ..
3. make
4. mpiexec -np num_procs ./mom_mpi 

















