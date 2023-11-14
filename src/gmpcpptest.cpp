#include <iostream>
#include <cstdlib>
#include "lrsinterface.hpp"

FILE *lrs_ifp;
FILE *lrs_ofp;

int main(int argc, char const **argv)
{
    
    Eigen::MatrixXd Ad(4,3);
    Ad <<   1, 1, 0,
            1, 0, 1,
            1, -1, 0,
            1, 0, -1;

    std::cout << Ad << std::endl;

    Polytope poly;
    std::cout << "Got to " << __LINE__ << " in " << __FILE__ << std::endl;
    poly.setVertices(fromMatrixXd(Ad));

    std::cout << "Got to " << __LINE__ << " in " << __FILE__ << std::endl;

    poly.getA();

    std::cout << "Got to " << __LINE__ << " in " << __FILE__ << std::endl;

    Eigen::MatrixXd Vmat = fromMatrixXq(poly.getA());

    std::cout << "Got to " << __LINE__ << " in " << __FILE__ << std::endl;

    std::cout << std::endl << Vmat << std::endl;

    return 0;
}
