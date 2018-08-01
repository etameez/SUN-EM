#ifndef QUADRATURE
#define QUADRATURE

#include <vector>
#include <array>
#include <iostream>

std::vector<std::array<double, 4>> getQuadratureWeightsAndValues(int num_quadrature_points);

#endif