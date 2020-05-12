#include "mpz.hpp"
#include <string>

mpz::~mpz(void) {
    mpz_clear(_data[0]);
    delete[] _data;
}

mpz::mpz() : _data(new mpz_t[1]) {
    mpz_init(_data[0]);
}

mpz::mpz(const mpz_t& c) : _data(new mpz_t[1]) {
    mpz_init_set(_data[0], c);
}

mpz::mpz(const mpz& c) : _data(new mpz_t[1]) {
    mpz_init_set(_data[0], c._data[0]);
}

mpz::mpz(int x) : _data(new mpz_t[1]) {
    std::string buffer = std::to_string(x);
    mpz_init_set_str(_data[0], buffer.data(), 10);
}

mpz::mpz(unsigned long int x) : _data(new mpz_t[1]) {
    mpz_init_set_ui(_data[0], x);
}

mpz::mpz(signed long int x) : _data(new mpz_t[1]) {
    mpz_init_set_si(_data[0], x);
}

mpz::mpz(double x) : _data(new mpz_t[1]) {
    mpz_init_set_d(_data[0], x);
}

mpz& mpz::operator=(const mpz& x) {
    mpz_set(this->_data[0], x._data[0]);
    return *this;
}

mpz& mpz::operator=(const mpz_t x) {
    mpz_set(this->_data[0], x);
    return *this;
}

bool mpz::operator<(const mpz& x) const {
    return (0>mpz_cmp(_data[0],x._data[0]));
}

bool mpz::operator==(const mpz& x) const {
    return (0==mpz_cmp(_data[0], x._data[0]));
}

mpz& mpz::operator+=(const mpz& x) {
    mpz_t tmp;
    mpz_init(tmp);
    mpz_add(tmp, this->_data[0], x._data[0]);
    mpz_swap(tmp, this->_data[0]);
    mpz_clear(tmp);
    return *this;
}

mpz& mpz::operator-=(const mpz& x) {
    mpz_t tmp;
    mpz_init(tmp);
    mpz_sub(tmp, this->_data[0], x._data[0]);
    mpz_swap(tmp, this->_data[0]);
    mpz_clear(tmp);
    return *this;
}

mpz& mpz::operator*=(const mpz& x){
    mpz_t tmp;
    mpz_init(tmp);
    mpz_mul(tmp, this->_data[0], x._data[0]);
    mpz_swap(tmp, this->_data[0]);
    mpz_clear(tmp);
    return *this;
}

mpz& mpz::operator/=(const mpz& x) {
    mpz_t tmp;
    mpz_init(tmp);
    mpz_tdiv_q(tmp, this->_data[0], x._data[0]);
    mpz_swap(tmp, this->_data[0]);
    mpz_clear(tmp);
    return *this;
}
    
bool mpz::operator<(signed long int x) const {
    return (0>mpz_cmp_si(_data[0], x));
}

bool mpz::operator>(signed long int x) const {
    return (0<mpz_cmp_si(_data[0], x));
}

bool mpz::operator==(signed long int x) const {
    return (0==mpz_cmp_si(_data[0], x));
}

const mpz_ptr mpz::get_mpz_t(void) const {
    return _data[0];
}

mpz_ptr mpz::get_mpz_t(void) {
    return const_cast<mpz_ptr>(const_cast<const mpz*>(this)->get_mpz_t());
}

size_t mpz::strlength(void) const {
    return mpz_sizeinbase(_data[0], 10);
}

std::ostream& operator<<(std::ostream &strm, const mpz &a) {
    auto length = a.strlength();
    
    char *buffer = new char[length];

    mpz_get_str(buffer, 10, a._data[0]);

    std::string sbuffer(buffer, length);
    delete[] buffer;

    return strm << sbuffer;
}