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
#include "node.h"
#include "edge.h"
#include "triangle.h"
#include "quadrature.h"
#include "timer.h"

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
        void calculateJMatrix();

        // Time profiling
        void timeProfiler(int num_iter);


    protected:
        std::vector<Node> nodes;
        std::vector<Triangle> triangles;
        std::vector<Edge> edges;
        std::vector<double> vrhs;
        std::vector<double> vrhs_internal;
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

        // Time profiling
        Timer z_mn_timer;
        Timer i_timer;
        Timer a_phi_timer;
        Timer j_timer;
        double i_time;
        double a_phi_time;
        double z_mn_time;
        double j_time;
};
#endif