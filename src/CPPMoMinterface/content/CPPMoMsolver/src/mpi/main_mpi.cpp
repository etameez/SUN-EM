/**
 *   \file main.cpp
 *   \brief C++ MoM Solver Entrypoint
 *
 *  Detailed description
 *  The entrypoint for a MoM solver that will be ran in parallel using MPI and OpenMP
 *  This will be used to read the .mom file that is created by MATLAB.
 *
 *  Author:    Tameez Ebrahim
 *  Created:   09 August 2018 
 *
 */

#include "../mom_file_reader.h"
#include "mom_solver_mpi.h"
#include "../timer.h"
#include <mpi.h>
#include <omp.h>

int main()
{
    int provided;
    MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &provided);

    int size;
    int rank;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    std::string path = "../../../../../examples/example-10/pec_plate_fine_mesh.mom";
    MoMFileReader reader(path); 
    MoMSolverMPI solver(reader.getNodes(), reader.getTriangles(), reader.getEdges(), reader.getVrhs(), reader.getConstMap());

    solver.calculateZmnByFaceMPI();
    if(rank == 0)
    {
        solver.calculateVrhsInternally();
    }
    solver.calculateJMatrixSCALAPACK();

    MPI_Finalize();
    return 0;
}
