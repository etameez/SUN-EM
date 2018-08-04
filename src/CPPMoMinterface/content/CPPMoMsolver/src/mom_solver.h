#ifndef MOM_SOLVER
#define MOM_SOLVER

#include <vector>
#include <map>
#include <string>
#include <complex>
#include <array>
#include <cmath>
#include <iostream>
#include <chrono>
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
                  std::vector<double> vrhs,
                  std::map<std::string, std::string> const_map);

        void calculateZmnByFace();
        void calculateJMatrix();

    protected:
        std::vector<Node> nodes;
        std::vector<Triangle> triangles;
        std::vector<Edge> edges;
        std::vector<double> vrhs;
        std::map<std::string, std::string> const_map;
        std::vector<std::array<double, 4>> quadrature_weights_values;
        std::vector<std::vector<std::complex<double>>> z_mn; 
        std::complex<double> j;
        double k; 
        double frequency;
        double omega;
        double lambda;
        double c;

        std::vector<Node> calculateAAndPhi(int p, int q);
        std::vector<std::complex<double>> calculateIpq(int p, int q);


    
};
#endif