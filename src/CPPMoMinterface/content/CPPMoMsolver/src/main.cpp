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

int main()
{
    // std::string path = "../../../../../examples/example-10/pec_plate_super_fine_mesh.mom";
    std::string path = "../../../../../examples/example-10/pec_plate.mom";
    MoMFileReader reader(path); 
    MoMSolver solver(reader.getNodes(), reader.getTriangles(), reader.getEdges(), reader.getVrhs(), reader.getConstMap());

    solver.calculateZmnByFace();
    solver.calculateVrhsInternally();
    solver.calculateJMatrixLAPACK();

    MoMFileWriter file_writer;
    std::string file_name = path.substr(0 , path.size() - 3) + "sol";
    file_writer.writeIlhsToFile(file_name, solver.getIlhs());
    return 0;
}
