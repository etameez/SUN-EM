#include "node.h"

/**
 *   \file node.cpp
 *   \brief Store node data
 *
 *  Detailed description
 *  This class will be used to store the nodes of the mesh. A node consists of x, y and z
 *  co-ordinates. This class will also be used to store vector data created by intermediate
 *  calculations.
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
    // The constructor for a non-complex node
    this->x = x_coord;
    this->y = y_coord;
    this->z = z_coord;  
    this->isComplex = false;
}

Node::Node(std::complex<double> x_coord, std::complex<double> y_coord, std::complex<double> z_coord)
{
    // The constructor for a complex node
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
    // Lets get the scalar multiplication of a node by a real number
    // Node = a * (x, y, z)
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
    // Lets get the scalar multiplication of a node by a complex number
    // Node = (a + bi) * (x, y, z) 
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
    // Lets add two nodes together
    // The if statements are needed to take into account the complexity of both nodes
    // It is important to remember that each node could be either complex or real
    // sum = (x1 + x2, y1 + y2, z1 + z2)
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
    // Lets subtract this node from another
    // The if statements are needed to take into account the complexity of both nodes
    // It is important to remember that the current node is either real or complex, 
    // but the input node is always complex
    // sum = (x1 + x2, y1 + y2, z1 + z2)

    // TODO: Rename 

    Node return_node;

    if(this->isComplex)
    {
        if(node.getIsComplex())
        {
            return_node = Node(this->x_complex - node.getXComplexCoord(), 
                               this->y_complex - node.getYComplexCoord(), 
                               this->z_complex - node.getZComplexCoord());
        }
        else
        {
            return_node = Node(this->x_complex - node.getXCoord(),
                               this->y_complex - node.getYCoord(),
                               this->z_complex - node.getZCoord());
        }
    }
    else
    {
        if(node.getIsComplex())
        {
            return_node = Node(this->x - node.getXComplexCoord(), 
                               this->y - node.getYComplexCoord(), 
                               this->z - node.getZComplexCoord());
        }
        else
        {
            return_node = Node(this->x - node.getXCoord(),
                               this->y - node.getYCoord(),
                               this->z - node.getZCoord());
        }
    }
    return return_node;
}

std::complex<double> Node::getDot(Node node)
{
    // Lets compute the dot product of this node by another
    // It is important to remember that this node can be either complex or real,
    // but the input node is always real
    // dot_product = (x1 * x2) + (y1 * y2) + (z1 * z2)

    // TODO: change to be more general. Input node could be complex or real

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
    // No need for this I think
    // TODO: check how this is used and change to getDot(Node node)
    // Used in Vrhs calculations
    // Needs to return a double
    // But, vrhs needs to be complex aswell for calculation reasons
    // TODO: change vrhs to complex
    return (this->x * node.getXCoord()) +
           (this->y * node.getYCoord()) + 
           (this->z * node.getZCoord());
}

bool Node::getIsComplex()
{
    // Lets check if the node is complex
    return this->isComplex;
}

Node Node::getCrossProduct(Node node)
{
    // TODO: Make general
    
    // Lets get the cross product of two vectors
    // http://tutorial.math.lamar.edu/Classes/CalcII/CrossProduct.aspx

    return Node(this->y * node.getZCoord() - this->z * node.getYCoord(),
                this->z * node.getXCoord() - this->x * node.getZCoord(),
                this->x * node.getYCoord() - this->y * node.getXCoord());           
}

Node Node::getScalarDivide(double scalar)
{
   // Lets divide a node by a real number
   
   Node return_node;

   if(this->isComplex)
   {
        return_node = Node(this->x_complex / std::complex<double>(scalar, 0),
                           this->y_complex / std::complex<double>(scalar, 0),
                           this->z_complex / std::complex<double>(scalar, 0));
   } 
   else
   {
        return_node = Node(this->x / scalar,
                           this->y / scalar,
                           this->z / scalar);
   }
   return return_node;
}

void Node::printNode()
{
    if(this->isComplex)
    {
        std::cout << "[" << this->x_complex << " "
                         << this->y_complex << " "
                         << this->z_complex << "]"
                         << std::endl;
    }
    else
    {
        std::cout << "[" << this->x << " "
                         << this->y << " "
                         << this->z << "]"
                         << std::endl;
    }
}

Node Node::getHat()
{
    return this->getScalarDivide(this->getNorm());
}


