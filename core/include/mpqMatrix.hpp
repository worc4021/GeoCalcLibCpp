#pragma once
#include <Eigen/Core>
#include "mpq.hpp"


namespace Eigen {
    template<> struct NumTraits<mpz> :
        NumTraits<signed long int> 
        {
            typedef mpz Real;
            typedef mpz Nested;
            typedef mpq NonInteger;
            typedef mpz Literal;
            enum {
                IsComplex = 0,
                IsInteger = 1,
                IsSigned = 1,
                RequireInitialization = 1,
                ReadCost = 1,
                AddCost = 3,
                MulCost = 3
            };
        };

    template<> struct NumTraits<mpq> :
        NumTraits<double>
        {
            typedef mpq Real;
            typedef mpq Nested;
            typedef mpq NonInteger;
            typedef mpq Literal;
            enum {
                IsComplex = 0,
                IsInteger = 0,
                IsSigned = 1,
                RequireInitialization = 1,
                ReadCost = 1,
                AddCost = 3,
                MulCost = 3
            };
        };
};




typedef Eigen::Matrix<mpq, Eigen::Dynamic, Eigen::Dynamic> MatrixXq;
typedef Eigen::Matrix<mpz, Eigen::Dynamic, Eigen::Dynamic> MatrixXz;
typedef Eigen::Matrix<mpq, 1, Eigen::Dynamic> RowXq;
typedef Eigen::Matrix<mpz, 1, Eigen::Dynamic> RowXz;

void canonicalise(MatrixXq& A);
// RowXz getDenominator(const MatrixXq& A, Eigen::Index row);
void getDenominator(const MatrixXq& A, Eigen::Index row, RowXz& currentRow);
// RowXz getNumerator(const MatrixXq& A, Eigen::Index row);
void getNumerator(const MatrixXq& A, Eigen::Index row, RowXz& currentRow);
void getRow(const MatrixXq& A, Eigen::Index row, RowXz& numerator, RowXz& denominator);
void setRow(MatrixXq& A, Eigen::Index row, const RowXq& newRow);
MatrixXq fromMatrixXd(const Eigen::MatrixXd& A);
Eigen::MatrixXd fromMatrixXq(const MatrixXq& A);