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

  // Lets open the file
  std::ifstream file(file_path);

  // Lets create a string to use for the received lines
  std::string str;

  // Lets create an in to hold the number of fields
  int num_fields = 0;

  // Lets create a string vector to store key-value pairs from MoMFileReader::constLineReader
  // This will be used for Const entries and entry numbers
  std::vector<std::string> line_vector; 

  if(file.is_open())
  {

    // Lets read first three lines
    // Basically header stuff
    getline(file, str);
    getline(file, str);
    getline(file, str);
    
    // Lets read Const
    getline(file, str);
    if(str == "CONST START")
    {
      // Lets get the number of fields
      getline(file, str);
      line_vector = this->constLineReader(str);
      num_fields = std::stoi(line_vector[1]);

      // Lets read the empty line
      getline(file, str);
      
      // Lets read all the Const entries now
      for(int i = 0; i < num_fields; i++)
      {
        getline(file, str);
        line_vector = this->constLineReader(str);
        
        // Now lets put the information from Const into the map (Dictionary) created in the class
        this->const_map[line_vector[0]] = line_vector[1];
      }
    }
    else
    {
      std::cout << "ERROR: THERE IS SOMETHING WRONG WITH THE FILE" << std::endl;
      std::cout << "ERROR: CANNOT FIND CONST ENTRIES" << std::endl;
    }

    // Lets check that all Const entries are read
    getline(file, str);

    if(!(str == "CONST END"))
    {
      std::cout << "ERROR: THERE IS SOMETHING WRONG WITH THE FILE" << std::endl;
      std::cout << "ERROR: THE NUMBER OF CONST ENTRIES ARE WRONG" << std::endl;
    }

    // Lets start getting the FEKO_data
    // Lets first read the empty line
    getline(file, str);

    // Now lets check for FEKO_data
    getline(file, str);

    if(str == "FEKO_DATA START")
    {
      // Lets first read the nodes
      getline(file, str);

      if(str == "NODES START")
      {
        // Lets get the number of fields i.e the number of nodes
        getline(file, str);
        line_vector = this->constLineReader(str);
        num_fields = std::stoi(line_vector[1]);

        // Lets read the node header
        getline(file, str);

        // Now lets read the nodes
        for(int i = 0; i < num_fields; i++)
        {
          getline(file, str);
          line_vector = this->numberLineReader(str, 3);

          // Lets add the values to the Node class and push it to node_vector
          // Also have to change the string to float
          this->node_vector.push_back(Node(std::stof(line_vector[0]),
                                           std::stof(line_vector[1]),
                                           std::stof(line_vector[2])));
        }
      }
      else
      {
        std::cout << "ERROR: THERE IS SOMETHING WRONG WITH THE FILE" << std::endl;
        std::cout << "ERROR: CANNOT FIND NODES" << std::endl;
      }

      // Lets check if all the nodes have been read
      getline(file, str);

      if(!(str == "NODES END"))
      {
        std::cout << "ERROR: THERE IS SOMETHING WRONG WITH THE FILE" << std::endl;
        std::cout << "ERROR: THE NUMBER OF NODE ENTRIES ARE WRONG" << std::endl;
      }

      // Lets start reading the triangles
      // First read the empty line
      getline(file, str);

      // Lets check for the triangles
      getline(file, str);

      if(str == "TRIANGLES START")
      {
        // First lets get the number of triangles
        getline(file, str);
        line_vector = this->constLineReader(str);
        num_fields = std::stoi(line_vector[1]);

        // Now lets read the empty line and the header
        getline(file, str);
        getline(file, str);

        // Now lets read the triangle data
        // The data is of the form VERTEX-VERTEX-VERTEX-CENTREX-CENTREY-CENTREZ-AREA
        // The vertices refer to the index of the nodes
        for(int i = 0; i < num_fields; i++)
        {
          getline(file, str);

          line_vector = this->numberLineReader(str, 7);

          // Lets read to the triangle class
          // Lets first make a node to store the centre
          // Dont forget to convert from string
          Node centre_node(std::stof(line_vector[3]),
                           std::stof(line_vector[4]),
                           std::stof(line_vector[5]));

          // Now lets read to a triangle
          Triangle triangle(std::stoi(line_vector[0]),
                            std::stoi(line_vector[1]),
                            std::stoi(line_vector[2]),
                            centre_node,
                            std::stof(line_vector[6]));

          // Finally, lets push to vector
          this->triangles.push_back(triangle);

        }
      }
      else
      {
        std::cout << "ERROR: THERE IS SOMETHING WRONG WITH THE FILE" << std::endl;
        std::cout << "ERROR: CANNOT FIND TRIANGLE DATA" << std::endl;
      }

      // Lets check that all the triangle data was read
      getline(file, str);
      if(!(str == "TRIANGLES END"))
      {
        std::cout << "ERROR: THERE IS SOMETHING WRONG WITH THE FILE" << std::endl;
        std::cout << "ERROR: THE NUMBER OF TRAINGLE ENTRIES ARE WRONG" << std::endl;
      }

      // Lets read the edge data
      // Start by reading the empty line
      getline(file, str);

      // Now lets check for the edges
      getline(file, str);

      if(str == "EDGES START")
      {
        // Lets get the number of edges
        getline(file, str);
        line_vector = this->constLineReader(str);
        num_fields = std::stoi(line_vector[1]);

        // Now lets read the empty line and the header
        getline(file, str);
        getline(file, str);

        // Lets read the edge data
        // The data is of the form VERTEX-VERTEX-CENTRE_X-CENTRE_Y-CENTRE_Z-LENGTH-TRIANGLE_MINUS-
        // TRIANGLE_PLUS-MINUS_FREE_VERTEX-PLUS_FREE_VERTEX-RHO_C_MINUS_X-RHO_C_MINUS_Y-
        // RHO_C_MINUS_Z-RHO_C_PLUS_X-RHO_C_PLUS_Y-RHO_C_PLUS_Z
        for(int i = 0; i < num_fields; i++)
        {
          getline(file, str);

          //TODO ADD STRUCT
          line_vector = this->numberLineReader(str, 16);
        }
      }
      else
      {
        std::cout << "ERROR: THERE IS SOMETHING WRONG WITH THE FILE" << std::endl;
        std::cout << "ERROR: CANNOT FIND EDGE DATA" << std::endl;
      }

      // Lets check that all the edges were read
      getline(file, str);

      if(!(str == "EDGES END"))
      {
        std::cout << "ERROR: THERE IS SOMETHING WRONG WITH THE FILE" << std::endl;
        std::cout << "ERROR: THE NUMBER OF TRAINGLE ENTRIES ARE WRONG" << std::endl;
      }
    }
    else
    {
      std::cout << "ERROR: THERE IS SOMETHING WRONG WITH THE FILE" << std::endl;
      std::cout << "ERROR: CANNOT FIND FEKO_DATA" << std::endl;
    }

    // Lets check that all the FEKO_data was read
    getline(file, str);

    if(!(str == "FEKO_DATA END"))
    {
      std::cout << "ERROR: THERE IS SOMETHING WRONG WITH THE FILE" << std::endl;
      std::cout << "ERROR: FEKO_DATA HAS NOT ENDED" << std::endl;
    }

    // Lets finally read the Vrhs data
    // Lets check if the data is available to read
    // Lets first start reading the empty line
    getline(file, str);
    getline(file, str);

    if(str == "VRHS START")
    {
      // Lets first get the number of values
      getline(file, str);
      line_vector = this->constLineReader(str);
      num_fields = std::stoi(line_vector[1]);

      // Now lets read the data
      // First lets read the empty line and the header
      getline(file, str);
      getline(file, str);

      // The data is of the form VALUE
      // No need for any processing
      for(int i = 0; i < num_fields; i++)
      {
        getline(file, str);

        //TODO ADD STRUCT
        std::cout << str << std::endl;
      }
    }
    else
    {
      std::cout << "ERROR: THERE IS SOMETHING WRONG WITH THE FILE" << std::endl;
      std::cout << "ERROR: VRHS CANNOT BE FOUND" << std::endl;
    }

    // Lets finally check that all Vrhs data was read
    getline(file, str);

    if(!(str == "VRHS END"))
    {
      std::cout << "ERROR: THERE IS SOMETHING WRONG WITH THE FILE" << std::endl;
      std::cout << "ERROR: THE NUMBER OF VRHS ENTRIES ARE WRONG" << std::endl;
    }

  }
  else
  {
    std::cout << "ERROR: THE FILE CANNOT BE OPENED" << std::endl;
  }
}

std::map<std::string, std::string> MoMFileReader::getConstMap()
{
  // Return the map of all the Const data
  return this->const_map;
}

std::vector<Triangle> MoMFileReader::getTriangles()
{
  return this->triangles;
}

std::vector<std::string> MoMFileReader::constLineReader(std::string line)
{
  // line is of the form
  // KEY..........--VALUE
  // The .'s are spaces and the -'s are tabs
  // The file is tab delimited but the spaces are introduced -
  // by matlab due to formatting reasons   


  // Lets get the index till the first space
  // Then cut the string till the first space to get the starting word
  int first_whitespace_char_index = std::strcspn(line.c_str(), " ");
  std::string key = line.substr(0, first_whitespace_char_index);

  // Lets reverse the line
  // Then lets get the index of the first tab
  // It looks like eulav--....yek 
  std::reverse(line.begin(), line.end());
  int last_whitespace_char_index = line.length() - std::strcspn(line.c_str(), "\t");

  // Lets then reverse back to normal to cut
  // Lets cut from the back till the index
  std::reverse(line.begin(), line.end());
  std::string value = line.substr(last_whitespace_char_index, line.length() - 1);

  // Lets return the values as a string vector
  std::vector<std::string> key_value_vector;
  key_value_vector.push_back(key);
  key_value_vector.push_back(value);
  return key_value_vector;
}

std::vector<std::string> MoMFileReader::numberLineReader(std::string line, int num_values)
{
  // line is of the form 
  // VALUE-VALUE-VALUE etc.
  // The -'s are tabs
  // num_values provides the number of values


  // Lets create a vector to return the values
  std::vector<std::string> value_vector;

  // Lets get the values with a tab find
  // Then lets add them to value_vector
  // Let us first create an index
  int index = 0;

  // The loop needs to be one less than the values because we don't need to erase the last one
  for(int i = 0; i < num_values - 1; i++)
  {
    // Lets find the first tab
    index = line.find('\t');

    // Now lets add the string till the first tab
    value_vector.push_back(line.substr(0, index));

    // Then lets erase the the added value as well as the tab
    line.erase(0, index + 1);
  }

  // Lets add the last value
  value_vector.push_back(line);

  // TODO DELETE THIS
  for(int i =0; i < value_vector.size(); i++)
  {
    std::cout << value_vector[i] << std::endl;
  }
  std::cout << std::endl;

  // Lets return the values
  return value_vector;
}