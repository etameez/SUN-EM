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
#include <fstream>
#include <mpi.h>
#include <omp.h>
#include "../node.h"
#include "../edge.h"
#include "../triangle.h"
#include "../quadrature.h"
#include "../timer.h"

extern "C"
{
    void Cblacs_pinfo(int*, int*); 
    void Cblacs_get(int , int, int*);
    void Cblacs_gridinit(int*, const char*, int, int);
    void Cblacs_gridinfo(int, int*, int*, int*,int*);
    void Cblacs_barrier(int, const char*);

    void Czgesd2d(int, int, int, std::complex<double>*, int, int, int);
    void Czgerv2d(int, int, int, std::complex<double>*, int, int, int);

    int numroc_(int*, int*, int*, int*, int*);
    void descinit_(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*);

    void pzgetrf_(int*, int*, std::complex<double>*, int*, int*, int*, int*, int*);
    void pzgetrs_(const char*, int*, int*, std::complex<double>*, int*, int*,int* ,int*, std::complex<double>*, int*, int*, int*, int*);

}

class MoMSolverMPI
{
    public:
        MoMSolverMPI(std::vector<Node> nodes,
                  std::vector<Triangle> triangles,
                  std::vector<Edge> edges,
                  std::vector<double> vrhs,
                  std::map<std::string, std::string> const_map);

        void calculateVrhsInternally();
        void calculateZmnByFaceMPI();
        void calculateJMatrixSCALAPACK();

        std::vector<std::complex<double>> getIlhs();
        void printPartialZMN();

    protected:
        std::vector<Node> nodes;
        std::vector<Triangle> triangles;
        std::vector<Edge> edges;
        std::vector<double> vrhs;
        std::vector<std::complex<double>> vrhs_internal;
        std::map<std::string, std::string> const_map;
        std::vector<std::array<double, 4>> quadrature_weights_values;
        std::complex<double> *zmn;
        std::complex<double> j;
        std::vector<std::complex<double>> ilhs;
        double k; 
        double frequency;
        double omega;
        double lambda;
        double c;

        std::vector<Node> calculateAAndPhi(int p, int q);
        std::vector<std::complex<double>> calculateIpq(int p, int q);
        int numValuesMPI(int num_procs, int rank, int data_length);
        std::vector<double> workMPIMP(std::vector<int> p_values); // RENAME

        std::vector<int> getProcessGrid(int num_procs);

        int num_threads;

        std::vector<std::complex<double>> getIpqSING(int p);
        std::vector<double> getMatrixXVector(std::vector<std::array<double, 3>> matrix,
                                             std::vector<double> vec_tor);
        int num_points_j;
        int num_points_i;
        std::vector<std::array<double, 2>> quadrature_weights_values_j;
        std::vector<std::array<double, 2>> quadrature_weights_values_i;

};
#endif