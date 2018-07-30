#include "triangle.h"

/**
 *   \file triangle.cpp
 *   \brief Data structure to store triangles
 *
 *  Detailed description
 *  A data structure(class) to store a mesh triangle and it's relevant data. It will be able to 
 *  store the data read from FEKO through the .mom file(created using MATLAB). It can also calculate
 *  the relevant values(centre and area) on its own;
 *
 *  Author:  Tameez Ebrahim
 *  Created: 27 July 2018
 *
 */

Triangle::Triangle(int vertex_1, int vertex_2, int vertex_3, Node centre, float area)
{
	this->vertex_1 = vertex_1;	
	this->vertex_2 = vertex_2;	
	this->vertex_3 = vertex_3;
	this->centre = centre;
	this->area = area;	
}

Triangle::Triangle(int vertex_1, int vertex_2, int vertex_3,
						Node node_vertex_1, Node node_vertex_2, Node node_vertex_3)
{
	this->vertex_1 = vertex_1;	
	this->vertex_2 = vertex_2;	
	this->vertex_3 = vertex_3;

	// Lets calculate the centre and area since it's not given
	this->centre = this->calculateCentre(node_vertex_1, node_vertex_2, node_vertex_2);
	this->area = this->calculateArea(node_vertex_1, node_vertex_2, node_vertex_2);	
}

int Triangle::getVertex1()
{
	return this->vertex_1;
}

int Triangle::getVertex2()
{
	return this->vertex_2;
}

int Triangle::getVertex3()
{
	return this->vertex_3;
}

std::vector<int> Triangle::getVertices()
{
	// Lets create a vector of the vertices
	std::vector<int> vertices_vector;
	vertices_vector.push_back(this->vertex_1);
	vertices_vector.push_back(this->vertex_2);
	vertices_vector.push_back(this->vertex_3);

	return vertices_vector;
}

std::vector<int> Triangle::getEdges()
{
	return this->edge_indices;
}

float Triangle::getArea()
{
	return this->area;
}

Node Triangle::getCentre()
{
	return this->centre;
}

void Triangle::setEdgeIndex(int index)
{
	this->edge_indices.push_back(index);
}

float Triangle::calculateArea(Node node_vertex_1, Node node_vertex_2, Node node_vertex_3)
{
	// Lets calculate the area
	// Easiest way is to use Heron's formula
	// S = (A + B + C) / 2 where A, B and C are the side lengths
	// Area = sqrt(S(S-A)(S-B)(S-C))

	// First lets get A, B and C
	float a = node_vertex_1.getDistanceTo(node_vertex_2);
	float b = node_vertex_1.getDistanceTo(node_vertex_3);
	float c = node_vertex_2.getDistanceTo(node_vertex_3);
	
	// Then lets get S
	float s = (a + b + c) / 2;

	// Finally, lets get the area
	return std::sqrt(s * (s - a) * (s - b) * (s - c)); 	
}

Node Triangle::calculateCentre(Node node_vertex_1, Node node_vertex_2, Node node_vertex_3)
{
	// Lets calculate the centre
	// The centre is the average of the three nodes

	float x_coord = (node_vertex_1.getXCoord() + 
					 node_vertex_2.getXCoord() + 
					 node_vertex_3.getXCoord()) / 3;

	float y_coord = (node_vertex_1.getYCoord() + 
					 node_vertex_2.getYCoord() + 
					 node_vertex_3.getYCoord()) / 3;

	float z_coord = (node_vertex_1.getZCoord() + 
					 node_vertex_2.getZCoord() + 
					 node_vertex_3.getZCoord()) / 3;

	Node centre_node(x_coord, y_coord, z_coord);
	return centre_node;
}


























