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
#include "../mom_file_writer.h"
#include <mpi.h>
#include <omp.h>

int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        std::cout << "ERROR: Invalid number of arguments" << std::endl;
    }
    else
    {
        int provided;
        MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &provided);

        int size;
        int rank;

        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        std::string path = argv[1];
        MoMFileReader reader(path);
    
        MoMSolverMPI solver(reader.getNodes(), reader.getTriangles(), reader.getEdges(), reader.getVrhs(), reader.getConstMap());
    
        solver.calculateZmnByFaceMPI();
        solver.calculateVrhsInternally();

        solver.calculateJMatrixSCALAPACK();

        if(rank == 0)
        {
            MoMFileWriter file_writer;
            std::string file_name = path.substr(0 , path.size() - 3) + "sol";
            file_writer.writeIlhsToFile(file_name, solver.getIlhs());
        }

        MPI_Finalize();
    }
    return 0;
}
