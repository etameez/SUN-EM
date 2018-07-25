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

int main()
{
  std::string path = "/Users/tameez/Dropbox/dipoles.mom";
  MoMFileReader reader(path);
  return 0;
}
