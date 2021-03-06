#ifndef TRIANGLE
#define TRIANGLE

#include <vector>
#include "node.h"

class Triangle
{
    public:
        Triangle(int vertex_1, int vertex_2, int vertex_3, Node centre, double area);
        Triangle(int vertex_1, int vertex_2, int vertex_3,
                    Node node_vertex_1, Node node_vertex_2, Node node_vertex_3);

        int getVertex1();
        int getVertex2();
        int getVertex3();
        std::vector<int> getVertices();
        std::vector<int> getEdges();

        double getArea();
        Node getCentre();

        void setEdgeIndex(int index);

    protected:
        int vertex_1;
        int vertex_2;
        int vertex_3;
        double area;
        Node centre;
        std::vector<int> edge_indices;

        double calculateArea(Node node_vertex_1, Node node_vertex_2, Node node_vertex_3);
        Node calculateCentre(Node node_vertex_1, Node node_vertex_2, Node node_vertex_3);
};

#endif