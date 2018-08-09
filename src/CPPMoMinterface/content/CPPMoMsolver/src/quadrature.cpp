#include "quadrature.h"

// TODO add in 3 and 12 quadrature points

std::vector<std::array<double, 4>> getQuadratureWeightsAndValues(int num_quadrature_points)
{
    // Lets return the quadrature weights and values
    // It is of the form [weight lambda_1 lambda_2 lambda_3]
    // with the number of rows equalling the number of quadrature points

    // Lets create the return vector
    std::vector<std::array<double, 4>> quadrature_weights_values;    
    
    switch(num_quadrature_points)
    {
        case 6:
        {
            std::array<double, 4> point_1 = {0.223381589678011, 0.108103018168070, 0.445948490915965, 0.445948490915965};  
            std::array<double, 4> point_2 = {0.223381589678011, 0.445948490915965, 0.108103018168070, 0.445948490915965};  
            std::array<double, 4> point_3 = {0.223381589678011, 0.445948490915965, 0.445948490915965, 0.108103018168070};  
            std::array<double, 4> point_4 = {0.109951743655322, 0.816847572980459, 0.091576213509771, 0.091576213509771};  
            std::array<double, 4> point_5 = {0.109951743655322, 0.091576213509771, 0.816847572980459, 0.091576213509771};  
            std::array<double, 4> point_6 = {0.109951743655322, 0.091576213509771, 0.091576213509771, 0.816847572980459};

            quadrature_weights_values.push_back(point_1);   
            quadrature_weights_values.push_back(point_2);   
            quadrature_weights_values.push_back(point_3);   
            quadrature_weights_values.push_back(point_4);   
            quadrature_weights_values.push_back(point_5);   
            quadrature_weights_values.push_back(point_6);   
            break;
        }

        default:
            std::cout << "ERROR: THERE IS NO DATA FOR THIS NUMBER OF QUADRATURE POINTS" << std::endl;
    }   

    return quadrature_weights_values;
}   