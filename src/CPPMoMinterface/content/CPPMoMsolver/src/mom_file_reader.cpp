#include "mom_file_reader.h"

/**
 *   \file mom_file_reader.cpp
 *   \brief Read the .mom file
 *
 *  Detailed description
 *  Read the .mom file and assign the information to all the relevant
 *  data structures.
 *
 *  Author:  Tameez Ebrahim
 *  Created: 25 July 2018
 *
 */

MoMFileReader::MoMFileReader(std::string file_path)
{
  using namespace std;

  ifstream file(file_path);
  string str;

  if(file.is_open())
  {
    while(getline(file, str))
    {
      cout << str << endl;
    }
  }
  else
  {
    cout << "ERROR: The file cannot be opened." << endl;
  }
}

