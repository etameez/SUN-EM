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
        void printZmnToFile(std::string path);
        std::vector<double> getMatrixXVector(std::vector<std::array<double, 3>> matrix,
                                             std::vector<double> vec_tor);
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

        int num_points_j;
        int num_points_i;
        std::vector<std::array<double, 2>> quadrature_weights_values_j;
        std::vector<std::array<double, 2>> quadrature_weights_values_i;

        std::vector<Node> calculateAAndPhi(int p, int q);
        std::vector<std::complex<double>> calculateIpq(int p, int q);

        std::vector<std::complex<double>> getIpqSING(int p);


};
#endif