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
#include "timer.h"

int main()
{
    std::string path = "../../../../../examples/example-10/pec_plate_fine_mesh.mom";
    MoMFileReader reader(path); 
    MoMSolver solver(reader.getNodes(), reader.getTriangles(), reader.getEdges(), reader.getVrhs(), reader.getConstMap());

    solver.calculateZmnByFace();
    solver.calculateVrhsInternally();
    solver.calculateJMatrix();
    return 0;
}
