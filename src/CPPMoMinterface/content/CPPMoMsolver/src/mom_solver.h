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
#include <Eigen/Dense>
#include <fstream>
#include "node.h"
#include "edge.h"
#include "triangle.h"
#include "quadrature.h"
#include "timer.h"

extern "C"
{
    void zgetrf_(int*, int*, std::complex<double>*, int*, int*, int*);
    void zgetrs_(const char*, int*, int*, std::complex<double>*, int*, int*, std::complex<double>*, int*, int*);
}

class MoMSolver
{
    public:
        MoMSolver(std::vector<Node> nodes,
                  std::vector<Triangle> triangles,
                  std::vector<Edge> edges,
                  std::vector<double> vrhs,
                  std::map<std::string, std::string> const_map);

        void calculateZmnByFace();
        void calculateVrhsInternally();
        void calculateJMatrixLAPACK();
        std::vector<std::complex<double>> getIlhs();

    protected:
        std::vector<Node> nodes;
        std::vector<Triangle> triangles;
        std::vector<Edge> edges;
        std::vector<double> vrhs;
        std::map<std::string, std::string> const_map;
        std::vector<std::array<double, 4>> quadrature_weights_values;
        std::complex<double> *z_mn;
        std::vector<std::complex<double>> vrhs_internal;
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