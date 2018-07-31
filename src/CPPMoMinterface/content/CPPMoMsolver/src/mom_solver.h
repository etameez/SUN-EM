#ifndef MOM_SOLVER
#define MOM_SOLVER

#include <vector>
#include <map>
#include <string>
#include <complex>
#include <array>
#include <cmath>
#include <iostream>
#include "node.h"
#include "edge.h"
#include "triangle.h"
#include "quadrature.h"

class MoMSolver
{
    public:
        MoMSolver(std::vector<Node> nodes,
                  std::vector<Triangle> triangles,
                  std::vector<Edge> edges,
                  std::vector<float> vrhs,
                  std::map<std::string, std::string> const_map);

        void calculateZmnByFace();

    protected:
        std::vector<Node> nodes;
        std::vector<Triangle> triangles;
        std::vector<Edge> edges;
        std::vector<float> vrhs;
        std::map<std::string, std::string> const_map;
        std::vector<std::array<float, 4>> quadrature_weights_values;
        std::complex<float> j;
        float k; 
        float frequency;
        float omega;
        float lambda;
        float c;

        std::vector<Node> calculateAAndPhi(int p, int q);
        std::vector<std::complex<float>> calculateIpq(int p, int q);


    
};
#endif