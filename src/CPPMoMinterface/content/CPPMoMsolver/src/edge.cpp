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
        this->vertex_1 = vertex_1;
        this->vertex_2 = vertex_2;
        this->centre = centre;
        this->length = length;
        this->minus_triangle_index = minus_triangle_index;
        this->plus_triangle_index = plus_triangle_index;
        this->minus_free_vertex = minus_free_vertex;
        this->plus_free_vertex = plus_free_vertex;
        this->rho_c_minus = rho_c_minus;
        this->rho_c_plus = rho_c_plus;  
}

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



