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
    std::string path = "../../../../../examples/example-10/pec_plate.mom";
    MoMFileReader reader(path); 
    MoMSolver solver(reader.getNodes(), reader.getTriangles(), reader.getEdges(), reader.getVrhs(), reader.getConstMap());
    std::cout << "After Solver" << std::endl;

    solver.calculateZmnByFace();
    std::cout << "After ZMN" << std::endl;
    solver.calculateVrhsInternally();
    std::cout << "After VRHS" << std::endl;
    // solver.calculateJMatrix();
    solver.calculateJMatrixLAPACK();
    std::cout << "After J" << std::endl;
    return 0;
}
