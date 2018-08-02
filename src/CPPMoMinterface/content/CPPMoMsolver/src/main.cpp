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
#include <mpi.h>

int main(int argc, char *argv[])
{
  // std::string path = "/home/tameez/Dropbox/pec_plate.mom";
  // // std::string path = "C:\\Users\\Tameez\\Dropbox\\pec_plate.mom";
  // MoMFileReader reader(path);
  // MoMSolver solver(reader.getNodes(), reader.getTriangles(), reader.getEdges(), reader.getVrhs(), reader.getConstMap());
  //solver.calculateZmnByFace();

  // Lets test MPI
  int id;
  int ierr;
  int p;
  double wtime;

  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &p);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &id);

  if(id == 0)
  {
    std::cout << "Master. Num processes is " << p << std::endl;
  }

  std::cout << "Process: " << id << std::endl;

  MPI_Finalize();
  
    
  return 0;
}
