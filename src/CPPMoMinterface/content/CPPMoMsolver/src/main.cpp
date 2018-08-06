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
#include <Eigen/Dense>
using Eigen::MatrixXd;

int main()
{
  // std::string path = "/home/tameez/Dropbox/pec_plate.mom";
  // std::string path = "C:\\Users\\Tameez\\Dropbox\\pec_plate.mom";
  std::string path = "/Users/tameez/Dropbox/pec_plate.mom";
  MoMFileReader reader(path);
  MoMSolver solver(reader.getNodes(), reader.getTriangles(), reader.getEdges(), reader.getVrhs(), reader.getConstMap());
  solver.calculateZmnByFace();
  solver.calculateVrhsInternally();
  solver.calculateJMatrix();
  

  // for(int i = 0; i < 1000; i++)
  // { 
  //   solver.calculateZmnByFace();
  // }

  solver.timeProfiler(1);


  return 0;
}
