#ifndef EDGE
#define EDGE

#include "node.h"

class Edge
{
    public:
        Edge(int vertex_1,
             int vertex_2,
             Node centre,
             float length,
             int minus_triangle_index,
             int plus_triangle_index,
             int minus_free_vertex,
             int plus_free_vertex,
             Node rho_c_minus,
             Node rho_c_plus);
        Edge(int vertex_1, int vertex_2);

        int getVertex1();
        int getVertex2();
        Node getCentre();
        float getLength();
        int getMinusTriangleIndex();
        int getPlusTriangleIndex();
        int getMinusFreeVertex();
        int getPlusFreeVertex();
        Node getRhoCMinus();
        Node getRhoCPlus();

    protected:
        int vertex_1;
        int vertex_2;
        Node centre;
        float length;
        int minus_triangle_index;
        int plus_triangle_index;
        int minus_free_vertex;
        int plus_free_vertex;
        Node rho_c_minus;
        Node rho_c_plus;
    
};

#endif