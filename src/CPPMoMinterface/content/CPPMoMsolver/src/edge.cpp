#include "edge.h"

/**
 *   \file edge.cpp
 *   \brief Store edge data
 *
 *  Detailed description
 *  A data structure(class) to store a single edge. The data is read from FEKO through a .mom file
 *  generate by MATLAB.
 *
 *  Author:  Tameez Ebrahim
 *  Created: 27 July 2018
 *
 */

// This is the constructor for Edge when reading from the .mom file
// There are no methods called, just assigning of arguments to internal variables
Edge::Edge(int vertex_1,
           int vertex_2,
           Node centre,
           float length,
           int minus_triangle_index,
           int plus_triangle_index,
           int minus_free_vertex,
           int plus_free_vertex,
           Node rho_c_minus,
           Node rho_c_plus)
{
    // The two vertices are indices to a vector<Node>
    // They define the endpoints of the edge
    // o------o
    // This is illustrated above, with the two o's the vertices
    // and the -'s the edge
    this->vertex_1 = vertex_1;
    this->vertex_2 = vertex_2;

    // The centre of the edge is of yet uneeded and may be removed
    this->centre = centre;

    this->length = length;

    // The minus and plus triangle indices are each indexes to a vector<Triangle>
    // They define the positivity of the triangles TODO FINISH
    this->minus_triangle_index = minus_triangle_index;
    this->plus_triangle_index = plus_triangle_index;

    // The minus and plus free vertices refer to the vertex which is not part of the
    // edge in the triangle. This is illustrated below
    // O\
    // | >  <-- This is a triangle
    // O/
    // In the above drawing where the O's are the edge vertices and the | is the edge,
    // the > is the free vertex
    // This is necessary to calculate equation 32 in RWG80 and is the term ri
    this->minus_free_vertex = minus_free_vertex;
    this->plus_free_vertex = plus_free_vertex;

    // The two rho's are again associated with the minus/plus triangles
    // They form part of equation 17 in RWG80 where depending on the method,
    // one or both will be used
    // Remember that these are vectors and therefore contain the x, y and z points
    // Thus, they are placed into class Node
    this->rho_c_minus = rho_c_minus;
    this->rho_c_plus = rho_c_plus;  
}

// This is the constructor for Edge when reading directly from a mesh
// Just the two vertices are given and assigned internally
// All the other internal variables will be calculated and assigned
// externally in another class
// See above for the other internal variables
Edge::Edge(int vertex_1, int vertex_2)
{
    this->vertex_1 = vertex_1;
    this->vertex_2 = vertex_2;
}

int Edge::getVertex1()
{
    return this->vertex_1;
}

int Edge::getVertex2()
{
    return this->vertex_2;
}

Node Edge::getCentre()
{
    return this->centre;
}

float Edge::getLength()
{
    return this->length;
}

int Edge::getMinusTriangleIndex()
{
    return this->minus_triangle_index;
}

int Edge::getPlusTriangleIndex()
{
    return this->plus_triangle_index;
}

int Edge::getMinusFreeVertex()
{
    return this->minus_free_vertex;
}

int Edge::getPlusFreeVertex()
{
    return this->plus_free_vertex;
}

Node Edge::getRhoCMinus()
{
    return this->rho_c_minus;
}

Node Edge::getRhoCPlus()
{
    return this->rho_c_plus;
}



