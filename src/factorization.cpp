#include "factorization.h"
#include <iostream>
#include "boost/math/ccmath/abs.hpp"

enum
{
    DEF_C_VALUE = 3,
    DEF_FREQ_VALUE = 10,
    DEF_MAX_ITER_VALUE = 10000000
};

namespace lpn
{
Factor BasicFactorization(long_int a)
{
    Factor factor;
    if (auto counter = FullDiv(a, 2); counter > 0)
        factor[2] += counter;

    long_int div = 3;
    while (div * div <= a)
    {
        if (auto counter = FullDiv(a, div); counter > 0)
            factor[div] += counter;
        div += 2;
    }
    if (a > 1)
        factor[a]++;
    return factor;
}

RhoFactorization::RhoFactorization()
{
    c = DEF_C_VALUE;
    max = DEF_MAX_ITER_VALUE;
    frequency = DEF_FREQ_VALUE;
}

RhoFactorization::RhoFactorization(const long_int & c, const size_t & max, const size_t & frequency) : c(c), max(max), frequency(frequency)
{
}

size_t RhoFactorization::check_gcd(long_int & n, const long_int & g)
{
    return g > 1 ? FullDiv(n, g) : 0;
}

Factor RhoFactorization::factorize(long_int n, long_int x1)
{
    const std::array<long_int, 6> primes = {2, 3, 5, 7, 11, 13};
    Factor factor;
    if (PseudoPrimeTest(n, primes))
    {
        std::cerr << "Warning: Value " << n << " is Pseudoprime\n";
        factor[n] += 1;
        return factor;
    }
    long_int x2 = next(x1, n);
    size_t range = 1;
    size_t terms = 0;
    while (terms <= max && !PseudoPrimeTest(n, primes))
    {
        long_int product = 1;
        for (size_t j = 0; j < range; ++j)
        {
            x2 = next(x2, n);
            product = (product * abs(x1 - x2)) % n;
            if (product == 0)
                product = 1;
            terms++;
            if (terms % frequency == 0 || j + 1 == range)
            {
                auto g = lpn::gcd(n, product);
                if (g > 1)
                    factor[g] += check_gcd(n, g);
                product = 1;
            }
        }
        x1 = x2;
        range <<= 1;
        for (size_t j = 0; j < range; ++j)
            x2 = next(x2, n);
    }
    if (n > 1)
        factor[n]++;
    if (factor.size() == 1)
        std::cerr << "No denominators were found. Please try again with different parameters.\n";
    return factor;
}

long_int RhoFactorization::next(const long_int & x2, const long_int & n)
{
    return (x2 * x2 + c) % n;
}
};
