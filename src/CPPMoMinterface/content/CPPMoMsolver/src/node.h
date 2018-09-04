#ifndef NODE
#define NODE

#include <cmath>
#include <complex>

class Node
{
    public:
        Node();
        Node(double x_coord, double y_coord, double z_coord);
        Node(std::complex<double> x_coord, std::complex<double> y_coord, std::complex<double> z_coord);

        double getXCoord();
        double getYCoord();
        double getZCoord();

        std::complex<double> getXComplexCoord();
        std::complex<double> getYComplexCoord();
        std::complex<double> getZComplexCoord();

        double getDistanceTo(Node node);
        Node getDifferenceBetween(Node node);
        double getNorm();
        Node getScalarMultiply(double scalar);
        Node getScalarMultiply(std::complex<double> scalar);
        Node getAddNode(Node node);
        Node getSubtractComplexNode(Node node);
        std::complex<double> getDot(Node node);
        double getDotNoComplex(Node node);
        bool getIsComplex();
        Node getCrossProduct(Node node);
        Node getScalarDivide(double scalar);

    protected:
        double x;
        double y;
        double z;

        bool isComplex;
        std::complex<double> x_complex;
        std::complex<double> y_complex;
        std::complex<double> z_complex;
};

#endif