#include "node.h"

/**
 *   \file node.cpp
 *   \brief Store node data
 *
 *  Detailed description
 *  This class will be used to store the nodes of the mesh. A node consists of x, y and z
 *  co-ordinates,
 *
 *  Author:  Tameez Ebrahim
 *  Created: 27 July 2018
 *
 */

Node::Node()
{
    
}

Node::Node(double x_coord, double y_coord, double z_coord)
{
    this->x = x_coord;
    this->y = y_coord;
    this->z = z_coord;  
    this->isComplex = false;
}

Node::Node(std::complex<double> x_coord, std::complex<double> y_coord, std::complex<double> z_coord)
{
    this->x_complex = x_coord;
    this->y_complex = y_coord;
    this->z_complex = z_coord;  
    this->isComplex = true;
}

double Node::getXCoord()
{
    return this->x;
}

double Node::getYCoord()
{
    return this->y;
}

double Node::getZCoord()
{
    return this->z;
}

std::complex<double> Node::getXComplexCoord()
{
    return this->x_complex;
}

std::complex<double> Node::getYComplexCoord()
{
    return this->y_complex;
}

std::complex<double> Node::getZComplexCoord()
{
    return this->z_complex;
}

double Node::getDistanceTo(Node node)
{
    // Lets get the distance
    // Distance = sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)
    return std::sqrt(std::pow((node.getXCoord() - this->x), 2) + 
                     std::pow((node.getYCoord() - this->y), 2) + 
                     std::pow((node.getZCoord() - this->z), 2));        
}

Node Node::getDifferenceBetween(Node node)
{
    // Lets get the difference between two nodes
    // Simply (x, y, z) = (x1 - x2, y1 - y2, z1 - z2)
    Node difference_node(this->x - node.getXCoord(),
                         this->y - node.getYCoord(),
                         this->z - node.getZCoord());

    return difference_node;
}

double Node::getNorm()
{
    // Lets get the norm of a node
    // norm = sqrt(x^2 + y^2 + z^2)
    return std::sqrt(std::pow(this->x, 2) + std::pow(this->y, 2) + std::pow(this->z, 2));   
}

Node Node::getScalarMultiply(double scalar)
{
    Node return_node;

    if(this->isComplex)
    {
        return_node = Node(scalar * this->x_complex, scalar * this->y_complex, scalar * this->z_complex);
    }
    else
    {
        return_node = Node(scalar * this->x, scalar * this->y, scalar * this->z);
    }
    return return_node;
}


Node Node::getScalarMultiply(std::complex<double> scalar)
{
    Node return_node;

    if(this->isComplex)
    {
        return_node = Node(scalar * this->x_complex, scalar * this->y_complex, scalar * this->z_complex);
    }
    else
    {
        return_node = Node(scalar * this->x, scalar * this->y, scalar * this->z);
    }
    return return_node;
}

Node Node::getAddNode(Node node)
{
    Node return_node;

    if(this->isComplex)
    {
        if(node.getIsComplex())
        {
            return_node = Node(this->x_complex + node.getXComplexCoord(), 
                               this->y_complex + node.getYComplexCoord(), 
                               this->z_complex + node.getZComplexCoord());
        }
        else
        {
            return_node = Node(this->x_complex + node.getXCoord(),
                               this->y_complex + node.getYCoord(),
                               this->z_complex + node.getZCoord());
        }
    }
    else
    {
        if(node.getIsComplex())
        {
            return_node = Node(this->x + node.getXComplexCoord(), 
                               this->y + node.getYComplexCoord(), 
                               this->z + node.getZComplexCoord());
        }
        else
        {
            return_node = Node(this->x + node.getXCoord(),
                               this->y + node.getYCoord(),
                               this->z + node.getZCoord());
        }
    }
    return return_node;
}

Node Node::getSubtractComplexNode(Node node)
{
    Node return_node;

    if(this->isComplex)
    {
        return_node = Node(this->x_complex - node.getXComplexCoord(), 
                           this->y_complex - node.getYComplexCoord(), 
                           this->z_complex - node.getZComplexCoord());
    }
    else
    {
        return_node = Node(this->x - node.getXComplexCoord(), 
                           this->y - node.getYComplexCoord(), 
                           this->z - node.getZComplexCoord());
    }
    return return_node;
}

std::complex<double> Node::getDot(Node node)
{
    if(this->isComplex)
    {
        return (this->x_complex * node.getXCoord()) +
               (this->y_complex * node.getYCoord()) + 
               (this->z_complex * node.getZCoord());
    }
    else
    {
        return (this->x * node.getXCoord()) +
               (this->y * node.getYCoord()) + 
               (this->z * node.getZCoord());
    }            
}

double Node::getDotNoComplex(Node node)
{
    return (this->x * node.getXCoord()) +
           (this->y * node.getYCoord()) + 
           (this->z * node.getZCoord());
}

bool Node::getIsComplex()
{
    return this->isComplex;
}