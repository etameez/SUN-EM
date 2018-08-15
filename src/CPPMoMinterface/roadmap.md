# Roadmap

## MoM

- [X] Output File For C++
- [X] Read .mom file in C++
- [X] Store the relevant information in C++ classes
- [X] MoM Solution using FEKO data
- [ ] Calculate missing data
- [X] MoM solution using calculated data
- [X] Compare diagonal of Zmn to FEKO Zmn matrix
- [ ] Add better singularity handling to fix Zmn diagonal discrepencies
- [X] Calculate Vrhs internally (Vm)
- [X] Calculate the surface current In
- [ ] Interface C++ with Matlab
- [ ] Test on 3D sphere

## HPC
- [X] Integrate MPI
- [X] Use MPI to hasten Zmn calculations
- [X] Use OpenMP to hasten Ipq calculations
- [ ] Use OpenMP to hasten Amn and Phimn calculations
- [X] Use OpenMP to hasten Eigen matrix calculation of In
- [X] Use SCALAPACK to hasten In calculation with MPI