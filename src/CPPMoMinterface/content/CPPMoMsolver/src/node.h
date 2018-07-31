#ifndef NODE
#define NODE

#include <cmath>
#include <complex>

class Node
{
    public:
        Node();
        Node(float x_coord, float y_coord, float z_coord);
        Node(std::complex<float> x_coord, std::complex<float> y_coord, std::complex<float> z_coord);

        float getXCoord();
        float getYCoord();
        float getZCoord();

        std::complex<float> getXComplexCoord();
        std::complex<float> getYComplexCoord();
        std::complex<float> getZComplexCoord();

        float getDistanceTo(Node node);
        Node getDifferenceBetween(Node node);
        float getNorm();
        Node getScalarMultiply(float scalar);
        Node getScalarMultiply(std::complex<float> scalar);
        Node getAddComplexNode(Node node);
        std::complex<float> getDot(Node node);

    protected:
        float x;
        float y;
        float z;

        bool isComplex;
        std::complex<float> x_complex;
        std::complex<float> y_complex;
        std::complex<float> z_complex;
};

#endif