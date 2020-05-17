#include <gmp.h>
#include <boost/operators.hpp>
#include <ostream>

class mpz 
    :   boost::less_than_comparable<mpz>,
        boost::equality_comparable<mpz>,
        boost::addable<mpz>,
        boost::subtractable<mpz>,
        boost::multipliable<mpz>,
        boost::dividable<mpz>,
        boost::less_than_comparable<mpz, signed long int>,
        boost::equality_comparable<mpz, signed long int>
     
{
private:
    mpz_t *_data;

    friend std::ostream& operator<<(std::ostream&, const mpz&);
    friend mpz abs(const mpz&);
public:
    ~mpz();

    mpz();
    mpz(const mpz_t&);
    mpz(const mpz&);
    mpz(int);
    mpz(unsigned long int);
    mpz(signed long int);
    mpz(double);

    mpz& operator=(const mpz&);
    mpz& operator=(const mpz_t);
    
    bool operator<(const mpz& x) const;
    bool operator==(const mpz& x) const;
    mpz& operator+=(const mpz& x);
    mpz& operator-=(const mpz& x);
    mpz& operator*=(const mpz& x);
    mpz& operator/=(const mpz& x);

    bool operator<(signed long int) const;
    bool operator>(signed long int) const;
    bool operator==(signed long int) const;

    mpz_ptr get_mpz_t(void);
    const mpz_ptr get_mpz_t(void) const;

    size_t strlength(void) const;

    mpz& operator-(void);
};
