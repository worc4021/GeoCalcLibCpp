#include "gmpMatrix.hpp"

void canonicalise(MatrixXq& A) {
    for (auto i = 0; i<A.rows(); i++){
        for (auto j = 0; j<A.cols(); j++)
            A(i,j).canonicalize();
    }
}

RowXz getDenominator(const MatrixXq& A, Eigen::Index row) {
    RowXz retVal(A.cols());
    for (auto i = 0; i<A.cols(); i++)
        retVal[i] = A(row, i).get_den();
    return retVal;
}

void getDenominator(const MatrixXq& A, Eigen::Index row, RowXz& currentRow) {
    assert(currentRow.cols() == A.cols());
    for (auto i = 0; i<A.cols(); i++)
        currentRow[i] = A(row, i).get_den();
}

RowXz getNumerator(const MatrixXq& A, Eigen::Index row) {
    RowXz retVal(A.cols());
    for (auto i = 0; i<A.cols(); i++)
        retVal[i] = A(row, i).get_num();
    return retVal;
}

void getNumerator(const MatrixXq& A, Eigen::Index row, RowXz& currentRow) {
    assert(currentRow.cols() == A.cols());
    for (auto i = 0; i<A.cols(); i++)
        currentRow[i] = A(row, i).get_num();
}

void getRow(const MatrixXq& A, Eigen::Index row, RowXz& numerator, RowXz& denominator) {
    assert(numerator.cols() == A.cols());
    assert(denominator.cols() == A.cols());
    for (auto i = 0; i < A.cols(); i++) {
        numerator[i] = A(row,i).get_num();
        denominator[i] = A(row,i).get_den();
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
            retVal(i,j) = mpq_class(A(i,j));

    return retVal;
}

Eigen::MatrixXd fromMatrixXq(const MatrixXq& A) {
    Eigen::MatrixXd retVal(A.rows(), A.cols());
    for (auto i = 0; i < A.rows(); i++)
        for (auto j = 0; j < A.cols(); j++)
            retVal(i,j) = A(i,j).get_d();

    return retVal;
}