#pragma once
#include <list>
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

class Polytope {
    private:
        MatrixXq HrepIneq;
        MatrixXq HrepEq;
        MatrixXq VrepV;
        MatrixXq VrepR;
        bool inHrep, inVrep;

        void vertexEnumerate(void);
        void facetEnumerate(void);
        double volume;

    public:
        Polytope() : HrepIneq(), HrepEq(), VrepV(), VrepR(), inHrep(false), inVrep(false), volume(0.) {}     

        void setHrepIneq(const MatrixXq& A, const MatrixXq& b) {
            HrepIneq = MatrixXq(A.rows(), A.cols()+1);
            HrepIneq << b,A;
            if (!inHrep)
                HrepEq = MatrixXq(0,HrepIneq.cols());
            inHrep = true;
            
        }

        void setHrepEq(const MatrixXq& A, const MatrixXq& b) {
            HrepEq = MatrixXq(A.rows(), A.cols()+1);
            HrepEq << b,A;
            if (!inHrep)
                HrepIneq = MatrixXq(0,HrepEq.cols());
            inHrep = true;
        }

        void setVertices(const MatrixXq& V) {
            VrepV = MatrixXq(V.rows(), V.cols()+1);
            VrepV << MatrixXq::Constant(V.rows(),1,1.), V;
            if (!inVrep)
                VrepR = MatrixXq(0,VrepV.cols());
            inVrep = true;
        }

        void setRays(const MatrixXq& R) {
            VrepR = MatrixXq(R.rows(), R.cols()+1);
            VrepR << MatrixXq::Zero(R.rows(),1),R;
            if (!inVrep)
                VrepV = MatrixXq(0,VrepR.cols());
            inVrep = true;
        }

        MatrixXq getVertices(void) {
            if (!inVrep)
                vertexEnumerate();
            return VrepV.block(0,1,VrepV.rows(), VrepV.cols()-1);
        };

        MatrixXq getRays(void) {
            if (!inVrep)
                vertexEnumerate();
            return VrepR.block(0,1,VrepR.rows(), VrepR.cols()-1);
        }

        double getVolume(void) {
            if (!inVrep)
                vertexEnumerate();
            return volume;
        }

        MatrixXq getA(void) {
            if (!inHrep)
                facetEnumerate();
            return HrepIneq.block(0,1,HrepIneq.rows(), HrepIneq.cols()-1);
        }

        MatrixXq getB(void) {
            if (!inHrep)
                facetEnumerate();
            return HrepIneq.col(0);
        }

        
};