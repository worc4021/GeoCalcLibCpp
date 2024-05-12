#pragma once
#include <list>
#include <vector>
#include <algorithm>
#if defined(GMPMATRIX)
#include "gmpMatrix.hpp"
typedef mpz_class zType;
typedef mpq_class qType;
#else
#include "mpqMatrix.hpp"
typedef mpz zType;
typedef mpq qType;
#endif

extern "C" {
#include "stdio.h"
#include "lrsrestart.h"
#include "lrslib.h"
}

struct Hrep {
    std::size_t nDim{}; //< Dimension of vertices as well as x in b+A*x>=0
    std::vector<RowXq> rows{}; //< rows of the halfspace representation b_i, -A_i
    std::vector<std::size_t> linearities{}; //< indices of linear constraints


    /**
     * @brief Append the latest return value from lrs to the H-representation buffer.
     *
     * @param output The lrs_mp_vector to be added as a new row.
     * @param bLinearity A boolean flag indicating whether the row is a linearity.
     */
    void push_back(const lrs_mp_vector output, bool bLinearity) {
        RowXq retVal(nDim + 1);
        if (bLinearity) 
            linearities.push_back(rows.size());
        
        if (zero(output[0])) {
            // ray 
            RowXz ray(nDim + 1);
            for (auto i = 0; i < ray.cols(); i++)
                ray(i) = zType(output[i]);

            zType norm = -ray.squaredNorm();

            for (auto i = 0; i < retVal.cols(); i++)
                retVal(i) = qType(zType(output[i]),norm);
        } else {
            // handle hyperplanes with negative sign
            retVal(0) = qType((mpz_cmp_si(output[0],0L)+1) ? 1: -1);
            zType normaliser = -abs(zType(output[0]));
            for (auto i = 1; i <= nDim; i++) {
                retVal(i) = qType(zType(output[i]), normaliser);
            }
        }
        rows.push_back(retVal);
    }
};

struct Vrep {
    std::size_t nDim{};
    std::vector<RowXq> rows{};
    std::vector<bool> isVertex{};
    std::vector<std::size_t> linearities{};

    void push_back(const lrs_mp_vector output, bool bLinearity) {
        RowXq retVal(nDim + 1);
        if (bLinearity) 
            linearities.push_back(rows.size());

        isVertex.push_back(!zero(output[0]));
        if (zero(output[0])) {
            // ray 
            RowXz ray(nDim + 1);
            for (auto i = 0; i < ray.cols(); i++)
                ray(i) = zType(output[i]);

            zType norm(ray.squaredNorm());
            for (auto i = 0; i < retVal.cols(); i++)
                retVal(i) = qType(zType(output[i]),norm);
        } else {
            // handle hyperplanes with negative sign
            retVal(0) = qType((mpz_cmp_si(output[0],0L)+1) ? 1: -1);
            zType normaliser = abs(zType(output[0]));
            for (auto i = 1; i <= nDim; i++) {
                retVal(i) = qType(zType(output[i]), normaliser);
            }
        }
        rows.push_back(retVal);
    }
};


struct Polytope {
    
    Hrep hRep{};
    Vrep vRep{};

    void vertexEnumerate(void);
    void facetEnumerate(void);
    double volume;

    void setHrep(const MatrixXq& A, const MatrixXq& b, const std::vector<std::size_t>& linearities = std::vector<std::size_t>()) {
        vRep.nDim = A.cols();
        hRep.nDim = A.cols();
        MatrixXq Rows(A.rows(), A.cols() + 1);
        Rows << b, A;
        for (auto i = 0; i < A.rows(); i++) {
            hRep.rows.push_back(Rows.row(i));
        }
        hRep.linearities.resize(linearities.size());
        std::transform(linearities.begin(), linearities.end(), hRep.linearities.begin(), [](auto i) { return i; });
    }

    void setVrep(const MatrixXq& V, const std::vector<bool>& isVertex, const std::vector<std::size_t>& linearities = std::vector<std::size_t>()) {
        hRep.nDim = V.cols();
        vRep.nDim = V.cols();
        RowXq row(V.cols() + 1);
        const qType one_t(1L);
        const qType zero_t(0L);
        for (auto i = 0; i < V.rows(); i++) {
            row << (isVertex[i] ? one_t : zero_t), V.row(i);
            vRep.rows.push_back(row);
            vRep.isVertex.push_back(isVertex[i]);
        }
        vRep.linearities.resize(linearities.size());
        std::transform(linearities.begin(), linearities.end(), vRep.linearities.begin(), [](auto i) { return i; });
    }

    void getHrep(Hrep& _hRep) {
        if (!hRep.rows.size())
            facetEnumerate();
        _hRep.nDim = hRep.nDim;
        _hRep.rows = hRep.rows;
        _hRep.linearities.resize(hRep.linearities.size());
        std::transform(hRep.linearities.begin(), hRep.linearities.end(), _hRep.linearities.begin(), [](auto i) { return i; });
    }

    void getVrep(Vrep& _vRep) {
        if (!vRep.rows.size())
            vertexEnumerate();
        _vRep.nDim = vRep.nDim;
        _vRep.rows = vRep.rows;
        _vRep.isVertex = vRep.isVertex;
        _vRep.linearities.resize(vRep.linearities.size());
        std::transform(vRep.linearities.begin(), vRep.linearities.end(), _vRep.linearities.begin(), [](auto i) { return i; });
    }

        // void setHrepIneq(const MatrixXq& A, const MatrixXq& b) {
        //     HrepIneq = MatrixXq(A.rows(), A.cols()+1);
        //     HrepIneq << b,A;
        //     if (!inHrep)
        //         HrepEq = MatrixXq(0,HrepIneq.cols());
        //     inHrep = true;
            
        // }

        // void setHrepEq(const MatrixXq& A, const MatrixXq& b) {
        //     HrepEq = MatrixXq(A.rows(), A.cols()+1);
        //     HrepEq << b,A;
        //     if (!inHrep)
        //         HrepIneq = MatrixXq(0,HrepEq.cols());
        //     inHrep = true;
        // }

        // void setVertices(const MatrixXq& V) {
        //     VrepV = MatrixXq(V.rows(), V.cols()+1);
        //     VrepV << MatrixXq::Constant(V.rows(),1,1.), V;
        //     if (!inVrep)
        //         VrepR = MatrixXq(0,VrepV.cols());
        //     inVrep = true;
        // }

        // void setRays(const MatrixXq& R) {
        //     VrepR = MatrixXq(R.rows(), R.cols()+1);
        //     VrepR << MatrixXq::Zero(R.rows(),1),R;
        //     if (!inVrep)
        //         VrepV = MatrixXq(0,VrepR.cols());
        //     inVrep = true;
        // }

        // MatrixXq getVertices(void) {
        //     if (!inVrep)
        //         vertexEnumerate();
        //     return VrepV.block(0,1,VrepV.rows(), VrepV.cols()-1);
        // };

        // MatrixXq getRays(void) {
        //     if (!inVrep)
        //         vertexEnumerate();
        //     return VrepR.block(0,1,VrepR.rows(), VrepR.cols()-1);
        // }

        // double getVolume(void) {
        //     if (!inVrep)
        //         vertexEnumerate();
        //     return volume;
        // }

        // MatrixXq getA(void) {
        //     if (!inHrep)
        //         facetEnumerate();
        //     return HrepIneq.block(0,1,HrepIneq.rows(), HrepIneq.cols()-1);
        // }

        // MatrixXq getB(void) {
        //     if (!inHrep)
        //         facetEnumerate();
        //     return HrepIneq.col(0);
        // }

        // MatrixXq getAeq(void) {
        //     if (!inHrep)
        //         facetEnumerate();
        //     return HrepEq.block(0, 1, HrepEq.rows(), HrepEq.cols() - 1);
        // }

        // MatrixXq getBeq(void) {
        //     if (!inHrep)
        //         facetEnumerate();
        //     return HrepEq.col(0);
        // }

        
};