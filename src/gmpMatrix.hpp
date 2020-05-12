#pragma once
#include <Eigen/Dense>
#include <gmpxx.h>


typedef Eigen::Matrix<mpq_class, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXq;
typedef Eigen::Matrix<mpz_class, 1, Eigen::Dynamic, Eigen::RowMajor> RowXz;
typedef Eigen::Matrix<mpq_class, 1, Eigen::Dynamic, Eigen::RowMajor> RowXq;


void canonicalise(MatrixXq& A);
RowXz getDenominator(const MatrixXq& A, Eigen::Index row);
void getDenominator(const MatrixXq& A, Eigen::Index row, RowXz& currentRow);
RowXz getNumerator(const MatrixXq& A, Eigen::Index row);
void getNumerator(const MatrixXq& A, Eigen::Index row, RowXz& currentRow);
void getRow(const MatrixXq& A, Eigen::Index row, RowXz& numerator, RowXz& denominator);
void setRow(MatrixXq& A, Eigen::Index row, const RowXq& newRow);
MatrixXq fromMatrixXd(const Eigen::MatrixXd& A);
Eigen::MatrixXd fromMatrixXq(const MatrixXq& A);