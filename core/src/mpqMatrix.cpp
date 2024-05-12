#include "mpqMatrix.hpp"


void canonicalise(MatrixXq& A) {
    for (auto i = 0; i<A.rows(); i++){
        for (auto j = 0; j<A.cols(); j++)
            A(i,j).canonicalise();
    }
}

void getDenominator(const MatrixXq& A, Eigen::Index row, RowXz& currentRow) {
    assert(currentRow.cols() == A.cols());
    for (auto i = 0; i<A.cols(); i++)
        currentRow[i] = A(row, i).get_den_mpz_t();
}


void getNumerator(const MatrixXq& A, Eigen::Index row, RowXz& currentRow) {
    assert(currentRow.cols() == A.cols());
    for (auto i = 0; i<A.cols(); i++)
        currentRow[i] = A(row, i).get_num_mpz_t();
}

void getRow(const MatrixXq& A, Eigen::Index row, RowXz& numerator, RowXz& denominator) {
    assert(numerator.cols() == A.cols());
    assert(denominator.cols() == A.cols());
    for (auto i = 0; i < A.cols(); i++) {
        numerator[i] = A(row,i).get_num_mpz_t();
        denominator[i] = A(row,i).get_den_mpz_t();
    }
}

void setRow(MatrixXq& A, Eigen::Index row, const RowXq& newRow) {
    assert(A.cols() == newRow.cols());
    for (auto i = 0; i < A.cols(); i++)
        A(row, i) = newRow(i);
}


MatrixXq fromMatrixXd(const Eigen::MatrixXd& A) {
    MatrixXq retVal(A.rows(), A.cols());
    for (auto i = 0; i < A.rows(); i++)
        for (auto j = 0; j < A.cols(); j++)
            retVal(i,j) = mpq(A(i,j));

    return retVal;
}

Eigen::MatrixXd fromMatrixXq(const MatrixXq& A) {
    Eigen::MatrixXd retVal(A.rows(), A.cols());
    for (auto i = 0; i < A.rows(); i++)
        for (auto j = 0; j < A.cols(); j++)
            retVal(i,j) = A(i,j);

    return retVal;
}

VectorXq fromVectorXd(const Eigen::VectorXd& A) {
    VectorXq retVal(A.rows());
    for (auto i = 0; i < A.rows(); i++)
        retVal(i) = mpq(A(i));

    return retVal;
}