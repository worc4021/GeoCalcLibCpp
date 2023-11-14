#include "mpz.hpp"

class mpq
    :   boost::less_than_comparable<mpq>,
        boost::equality_comparable<mpq>,
        boost::addable<mpq>,
        boost::subtractable<mpq>,
        boost::multipliable<mpq>,
        boost::dividable<mpq>
{
private:
    mpq_t *_data;
    friend std::ostream& operator<<(std::ostream&, const mpq&);
public:
    mpq(void);
    mpq(const mpq&);
    mpq(double);
    mpq(unsigned long int, unsigned long int);
    mpq(signed long int, signed long int);
    mpq(const mpz&, const mpz&);
    ~mpq();

    mpq& operator=(const mpq&);

    operator double() const;

    bool operator<(const mpq&) const;
    bool operator==(const mpq&) const;
    mpq& operator+=(const mpq&);
    mpq& operator-=(const mpq&);
    mpq& operator*=(const mpq&);
    mpq& operator/=(const mpq&);

    void canonicalise(void);

    
    mpz_srcptr get_num_mpz_t(void);
    const mpz_srcptr get_num_mpz_t(void) const;
    mpz_srcptr get_den_mpz_t(void);
    const mpz_srcptr get_den_mpz_t(void) const;

    void set_num(const mpz&);
    void set_den(const mpz&);

    double getValue(void) const;
};

