#ifndef MOM_FILE_READER
#define MOM_FILE_READER

#include <string>
#include <fstream>
#include <iostream>
#include <map>
#include <algorithm>
#include <cstring>
#include <vector>
#include "node.h"

class MoMFileReader
{
 	public:
    	MoMFileReader(std::string file_path);

    	std::map<std::string, std::string> getConstMap();
  
 	protected:
    	std::map<std::string, std::string> const_map;
    	std::vector<Node> node_vector;

    	std::vector<std::string> constLineReader(std::string line);
    	std::vector<std::string> numberLineReader(std::string line, int num_values);
};

#endif
