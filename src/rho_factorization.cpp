#include "rho_factorization.hpp"
#include "basic.hpp"

#include <boost/random.hpp>
#include <iostream>

enum
{
    DEF_C_VALUE = 3,
    DEF_FREQ_VALUE = 10,
    DEF_MAX_ITER_VALUE = 10000000,
    MAX_ATTEMP = 100,
};

namespace lpn
{

RhoFactorization::RhoFactorization()
{
    c_ = DEF_C_VALUE;
    max_iter_amount_ = DEF_MAX_ITER_VALUE;
    frequency_ = DEF_FREQ_VALUE;
}

RhoFactorization::RhoFactorization(const long_int & c, const size_t & max, const size_t & frequency)
    : c_(c), max_iter_amount_(max), frequency_(frequency)
{
}

Factor RhoFactorization::Factorize(long_int n, long_int x1)
{
    std::cout << "Start to factorize number " << n << '\n';
    Factor factor;

    for (size_t i = 0; i < MAX_ATTEMP && n > 1; ++i)
    {
        auto divisor = FindNewDivisor(n, x1);
        factor[divisor] += FullDiv(n, divisor);
        NotifyUser(divisor);
    }
    PrintResult(factor);
    return factor;
}

long_int RhoFactorization::FindNewDivisor(long_int & n, long_int & x1)
{
    const std::array<long_int, 7> primes = {2, 3, 5, 7, 11, 13, 61631};
    if (PseudoPrimeTest(n, primes))
    {
        std::cerr << "Warning: Value " << n << " is Pseudoprime\n";
        return n;
    }
    long_int x2 = Next(x1, n);
    size_t range = 1;
    size_t terms = 0;
    while (terms <= max_iter_amount_ && !PseudoPrimeTest(n, primes))
    {
        long_int product = 1;
        for (size_t j = 0; j < range; ++j)
        {
            x2 = Next(x2, n);
            product = (product * abs(x1 - x2)) % n;
            if (product == 0)
            {
                product = 1;
            }
            terms++;
            if (terms % frequency_ == 0 || j + 1 == range)
            {
                auto g = lpn::gcd(n, product);
                if (g > 1)
                {
                    return g;
                }
            }
        }
        x1 = x2;
        range <<= 1;
        for (size_t j = 0; j < range; ++j)
        {
            x2 = Next(x2, n);
        }
    }
    std::cerr << "No denominators were found. Please try again with different parameters.\n";
    return n;
}

long_int RhoFactorization::Next(const long_int & x2, const long_int & n) { return (x2 * x2 + c_) % n; }

void RhoFactorization::PrintResult(const Factor & factor)
{
    std::cout << "RESULT:\n";
    for (const auto & [div, deg] : factor)
    {
        std::cout << "Divisor " << div << " "
                  << "Degree " << deg << '\n';
    }
}

void RhoFactorization::NotifyUser(const long_int & div) { std::cout << "New divisor has been found " << div << '\n'; }

};  // namespace lpn