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
#include <chrono>
 
#define TIMING
 
#ifdef TIMING
#define INIT_TIMER auto start = std::chrono::high_resolution_clock::now();
#define START_TIMER  start = std::chrono::high_resolution_clock::now();
#define STOP_TIMER(name)  std::cout << "RUNTIME of " << name << ": " << \
    std::chrono::duration_cast<std::chrono::milliseconds>( \
            std::chrono::high_resolution_clock::now()-start \
    ).count() << " ms " << std::endl; 
#else
#define INIT_TIMER
#define START_TIMER
#define STOP_TIMER(name)
#endif

#include "mom_file_reader.h"
#include "mom_solver.h"
#include "mom_file_writer.h"

int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        std::cout << "ERROR: Invalid number of arguments" << std::endl;
    }
    else
    {
        INIT_TIMER
        // START_TIMER
        std::string path = argv[1];
        MoMFileReader reader(path);
        // STOP_TIMER("FILE READ") 
        MoMSolver solver(reader.getNodes(), reader.getTriangles(), reader.getEdges(), reader.getVrhs(), reader.getConstMap());

        START_TIMER
        solver.calculateZmnByEdge();
        // solver.calculateZmnByFace();
        // solver.printZmnToFile("pec_plate.zmn");
        STOP_TIMER("ZMN")

        // START_TIMER
        solver.calculateVrhsInternally();
        // STOP_TIMER("VRHS")

        // START_TIMER
        solver.calculateJMatrixLAPACK();
        // STOP_TIMER("J")

        // START_TIMER
        MoMFileWriter file_writer;
        std::string file_name = path.substr(0 , path.size() - 3) + "sol";
        file_writer.writeIlhsToFile(file_name, solver.getIlhs());
        // STOP_TIMER("FWRITER")

        // TEST
    }
    return 0;
}
