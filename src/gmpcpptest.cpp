#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <numeric>
#include "lrsinterface.hpp"

FILE *lrs_ifp;
FILE *lrs_ofp;

int main(int argc, char const **argv)
{
    
    Eigen::MatrixXd Vd(4,2);
    Vd <<   1, 0,
            0, 1,
            -1, 0,
            0, -1;

    Polytope poly;
    Vrep vrep;
    Hrep hrep;
    auto V = fromMatrixXd(Vd);
    std::vector<bool> isVertex(V.rows(), true);

    poly.setVrep(V, isVertex);
    poly.getHrep(hrep);

    MatrixXq Hh(hrep.rows.size(), hrep.nDim + 1);
    for (std::size_t iRow = 0; iRow < hrep.rows.size(); iRow++) {
        Hh.row(iRow) = hrep.rows[iRow];
    }

    Eigen::MatrixXd H = fromMatrixXq(Hh);

    std::cout << std::endl << H << std::endl;

    Eigen::MatrixXd Ad(4,2);
    Ad <<   1, 0,
            0, 1,
            -1, 0,
            0, -1;
    Eigen::VectorXd bd(4);
    bd << 1, 1, 1, 1;

    Polytope poly2;
    auto A = fromMatrixXd(Ad);
    auto b = fromVectorXd(bd);

    poly2.setHrep(A, b);
    poly2.getVrep(vrep);

    std::size_t nVert = std::accumulate(vrep.isVertex.begin(), vrep.isVertex.end(), 0, [](auto a, auto b) { return a + b; });
    MatrixXq Vv(nVert, vrep.nDim + 1);
    MatrixXq Rv(vrep.rows.size() - nVert, vrep.nDim + 1);
    std::size_t iVert{0},iRay{0};
    for (std::size_t iRow = 0; iRow < vrep.rows.size(); iRow++) {
        if (vrep.isVertex[iRow])
            Vv.row(iVert++) = vrep.rows[iRow]; 
        else
            Rv.row(iRay++) = vrep.rows[iRow]; 
    }

    std::cout << std::endl << fromMatrixXq(Vv) << std::endl;

    // std::cout << "Got to " << __LINE__ << " in " << __FILE__ << std::endl;

    // Eigen::MatrixXd Vmat = fromMatrixXq(poly.getA());

    // std::cout << "Got to " << __LINE__ << " in " << __FILE__ << std::endl;

    // std::cout << std::endl << Vmat << std::endl;

    return 0;
}
