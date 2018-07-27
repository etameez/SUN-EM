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

Node::Node(float x_coord, float y_coord, float z_coord)
{
	this->x = x_coord;
	this->y = y_coord;
	this->z = z_coord;	
}

float Node::getXCoord()
{
	return this->x;
}

float Node::getYCoord()
{
	return this->y;
}

float Node::getZCoord()
{
	return this->z;
}