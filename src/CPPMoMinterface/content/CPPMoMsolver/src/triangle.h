#ifndef TRIANGLE
#define TRIANGLE

#include <vector>
#include "node.h"

class Triangle
{
	public:
		Triangle(int vertex_1, int vertex_2, int vertex_3, Node centre, float area);
		Triangle(int vertex_1, int vertex_2, int vertex_3,
				 	Node node_vertex_1, Node node_vertex_2, Node node_vertex_3);

		int getVertex1();
		int getVertex2();
		int getVertex3();
		std::vector<int> getVertices();

		float getArea();
		Node getCentre();

	protected:
		int vertex_1;
		int vertex_2;
		int vertex_3;
		float area;
		Node centre;


		float calculateArea(Node node_vertex_1, Node node_vertex_2, Node node_vertex_3);
		Node calculateCentre(Node node_vertex_1, Node node_vertex_2, Node node_vertex_3);
};

#endif