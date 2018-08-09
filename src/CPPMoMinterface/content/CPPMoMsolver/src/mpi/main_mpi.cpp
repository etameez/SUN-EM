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

        solver.calculateVrhsInternally();
        solver.calculateJMatrix();
    }

    // OpenMP Test
    if(rank == 0)
    {
        int n;
        int id;

        omp_set_num_threads(4);
        #pragma omp parallel private(n, id)
        {
            id = omp_get_thread_num();
            std::cout << "Thread no: " << id << std::endl;

            if(id == 0)
            {
                n = omp_get_num_threads();
                std::cout << "Num threads: " << n << std::endl;
            }
        } 
    }

    MPI_Finalize();
    return 0;
}
