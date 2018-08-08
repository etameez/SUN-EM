/**
 *   \file main.cpp
 *   \brief C++ MoM Solver Entrypoint
 *
 *  Detailed description
 *  The entrypoint for a MoM solver. This will be used to read the .mom
 *  file that is created by MATLAB. This is standalone and serves as an
 *  alternative to the mex file.
 *
 *  Author:    Tameez Ebrahim
 *  Created:   25 July 2018 
 *
 */

#include "../mom_file_reader.h"
#include "mom_solver_mpi.h"
#include "../timer.h"
#include <mpi.h>

int main()
{
    MPI_Init(NULL, NULL);

    int size;
    int rank;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    std::string path = "../../../../../examples/example-10/pec_plate.mom";
    MoMFileReader reader(path); 
    MoMSolverMPI solver(reader.getNodes(), reader.getTriangles(), reader.getEdges(), reader.getVrhs(), reader.getConstMap());

    Timer t;
    if(rank == 0)
    {
        t.startTimer();
    }

    solver.calculateZmnByFaceMPI();
    if(rank == 0)
    {
        t.endTimer();
        std::cout << "The MPI ZMN TIME: " << std::endl;
        t.printTime();
        std::cout << std::endl << std::endl;

        solver.calculateVrhsInternally();
        solver.calculateJMatrix();
    }

    MPI_Finalize();
  
    return 0;
}
