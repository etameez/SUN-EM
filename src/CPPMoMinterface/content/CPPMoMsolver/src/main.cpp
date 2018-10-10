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
        std::string path = argv[1];
        MoMFileReader reader(path); 
        MoMSolver solver(reader.getNodes(), reader.getTriangles(), reader.getEdges(), reader.getVrhs(), reader.getConstMap());

        //solver.calculateZmnByEdge();
        solver.calculateZmnByFace();
        // solver.printZmnToFile("pec_plate.zmn");
        solver.calculateVrhsInternally();
        solver.calculateJMatrixLAPACK();

        MoMFileWriter file_writer;
        std::string file_name = path.substr(0 , path.size() - 3) + "sol";
        file_writer.writeIlhsToFile(file_name, solver.getIlhs());

        // TEST
    }
    return 0;
}
