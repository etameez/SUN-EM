#ifndef MOM_SOLVER_MPI
#define MOM_SOLVER_MPI

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
#include <mpi.h>
#include "../node.h"
#include "../edge.h"
#include "../triangle.h"
#include "../quadrature.h"
#include "../timer.h"

class MoMSolverMPI
{
    public:
        MoMSolverMPI(std::vector<Node> nodes,
                  std::vector<Triangle> triangles,
                  std::vector<Edge> edges,
                  std::vector<double> vrhs,
                  std::map<std::string, std::string> const_map);

        void calculateVrhsInternally();
        void calculateJMatrix();
        void calculateZmnByFaceMPI();
        int numValuesMPI(int num_procs, int rank, int data_length);
        std::vector<double> workMPI(std::vector<int> p_values); // RENAME

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

};
#endif