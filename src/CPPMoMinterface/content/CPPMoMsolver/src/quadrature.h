#ifndef QUADRATURE
#define QUADRATURE

#include <vector>
#include <array>
#include <iostream>
#include <cmath>

std::vector<std::array<double, 4>> getQuadratureWeightsAndValues(int num_quadrature_points);
std::vector<std::array<double, 2>> getGaussLegendreQuadratureWeightsAndValues(int num_quadrature_points);

#endif