#include "mpq.hpp"
#include <string>
#include <iostream>

std::ostream& operator<<(std::ostream& strm, const mpq& x) {
    // mpq_canonicalize(x._data);
    mpz num, den;
    mpq_get_num(num.get_mpz_t(), x._data[0]);
    mpq_get_den(den.get_mpz_t(), x._data[0]);
    auto length = num.strlength() + den.strlength() + 3;
    char *buffer = new char[length];
    mpq_get_str(buffer, 10, x._data[0]);
    std::string sbuffer(buffer, length);
    delete[] buffer;
    return strm << sbuffer;
    
}

mpq::mpq(void) : _data(new mpq_t[1]) {
    mpq_init(_data[0]);
}

mpq::mpq(const mpq& x) : _data(new mpq_t[1]) {
    mpq_init(_data[0]);
    mpq_set(_data[0], x._data[0]);
}

mpq::mpq(double x) : _data(new mpq_t[1]) {
    mpq_init(_data[0]);
    mpq_set_d(_data[0], x);
}

mpq::mpq(unsigned long int num, unsigned long int den) : _data(new mpq_t[1]) {
    mpq_init(_data[0]);
    mpq_set_ui(_data[0], num, den);
    mpq_canonicalize(_data[0]);
}

mpq::mpq(signed long int num, signed long int den) : _data(new mpq_t[1]) {
    mpq_init(_data[0]);
    mpq_set_si(_data[0], num, den);
    mpq_canonicalize(_data[0]);
}

mpq::mpq(const mpz& num, const mpz& den) : _data(new mpq_t[1]) {
    mpq_init(_data[0]);
    mpq_set_num(_data[0], num.get_mpz_t());
    mpq_set_den(_data[0], den.get_mpz_t());
    mpq_canonicalize(_data[0]);
}

mpq::~mpq() {
    mpq_clear(_data[0]);
    delete[] _data;
}

mpq& mpq::operator=(const mpq& x) {
    mpq_set(this->_data[0], x._data[0]);
    return *this;
}

mpq::operator double() const {
    return mpq_get_d(_data[0]);
}

bool mpq::operator<(const mpq& x) const {
    return (0>mpq_cmp(_data[0], x._data[0]));
}

bool mpq::operator==(const mpq& x) const {
    return (0 != mpq_equal(_data[0], x._data[0]));
}

mpq& mpq::operator+=(const mpq& x) {
    mpq_t tempVal;
    mpq_init(tempVal);
    mpq_add(tempVal, this->_data[0], x._data[0]);
    mpq_swap(tempVal, this->_data[0]);
    mpq_clear(tempVal);
    return *this;
}

mpq& mpq::operator-=(const mpq& x) {
    mpq_t tempVal;
    mpq_init(tempVal);
    mpq_sub(tempVal, this->_data[0], x._data[0]);
    mpq_swap(tempVal, this->_data[0]);
    mpq_clear(tempVal);
    return *this;
}

mpq& mpq::operator*=(const mpq& x) {
    mpq_t tempVal;
    mpq_init(tempVal);
    mpq_mul(tempVal, this->_data[0], x._data[0]);
    mpq_swap(tempVal, this->_data[0]);
    mpq_clear(tempVal);
    return *this;
}

mpq& mpq::operator/=(const mpq& x) {
    mpq_t tempVal;
    mpq_init(tempVal);
    mpq_div(tempVal, this->_data[0], x._data[0]);
    mpq_swap(tempVal, this->_data[0]);
    mpq_clear(tempVal);
    return *this;
}

void mpq::canonicalise(void) {
    mpq_canonicalize(_data[0]);
}

const mpz_srcptr mpq::get_num_mpz_t(void) const {
    return mpq_numref(_data[0]);
}

mpz_srcptr mpq::get_num_mpz_t(void) {
    return const_cast<mpz_srcptr>(const_cast<const mpq*>(this)->get_num_mpz_t());
}

const mpz_srcptr mpq::get_den_mpz_t(void) const {
    return mpq_denref(_data[0]);
}

mpz_srcptr mpq::get_den_mpz_t(void) {
    return const_cast<mpz_srcptr>(const_cast<const mpq*>(this)->get_den_mpz_t());
}

void mpq::set_num(const mpz& num) {
    mpq_set_num(_data[0], num.get_mpz_t());
}

void mpq::set_den(const mpz& den) {
    mpq_set_den(_data[0], den.get_mpz_t());
}

double mpq::getValue(void) const {
    return mpq_get_d(_data[0]);
}