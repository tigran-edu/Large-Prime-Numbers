#include "factorization.h"

#include <iostream>

#include "basic.h"

enum
{
    DEF_C_VALUE = 3,
    DEF_FREQ_VALUE = 10,
    DEF_MAX_ITER_VALUE = 10000000,
    MAX_ATTEMP = 100,
};

namespace lpn
{
Factor BasicFactorization(long_int a)
{
    Factor factor;
    if (auto counter = FullDiv(a, 2); counter > 0)
    {
        factor[2] += counter;
    }

    long_int div = 3;
    while (div * div <= a)
    {
        if (auto counter = FullDiv(a, div); counter > 0)
        {
            factor[div] += counter;
        }
        div += 2;
    }
    if (a > 1)
    {
        factor[a]++;
    }
    return factor;
}

RhoFactorization::RhoFactorization()
{
    c = DEF_C_VALUE;
    max_iter_amount = DEF_MAX_ITER_VALUE;
    frequency = DEF_FREQ_VALUE;
}

RhoFactorization::RhoFactorization(const long_int & c, const size_t & max, const size_t & frequency)
    : c(c), max_iter_amount(max), frequency(frequency)
{
}

Factor RhoFactorization::factorize(long_int n, long_int x1)
{
    std::cout << "Start to factorize number " << n << '\n';
    Factor factor;

    for (size_t i = 0; i < MAX_ATTEMP && n > 1; ++i)
    {
        auto divisor = find_new_divisor(n, x1);
        factor[divisor] += FullDiv(n, divisor);
        notify_user(divisor);
    }
    print_result(factor);
    return factor;
}

long_int RhoFactorization::find_new_divisor(long_int & n, long_int & x1)
{
    const std::array<long_int, 7> primes = {2, 3, 5, 7, 11, 13, 61631};
    if (PseudoPrimeTest(n, primes))
    {
        std::cerr << "Warning: Value " << n << " is Pseudoprime\n";
        return n;
    }
    long_int x2 = next(x1, n);
    size_t range = 1;
    size_t terms = 0;
    while (terms <= max_iter_amount && !PseudoPrimeTest(n, primes))
    {
        long_int product = 1;
        for (size_t j = 0; j < range; ++j)
        {
            x2 = next(x2, n);
            product = (product * abs(x1 - x2)) % n;
            if (product == 0)
            {
                product = 1;
            }
            terms++;
            if (terms % frequency == 0 || j + 1 == range)
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
            x2 = next(x2, n);
        }
    }
    std::cerr << "No denominators were found. Please try again with different parameters.\n";
    return n;
}

long_int RhoFactorization::next(const long_int & x2, const long_int & n) { return (x2 * x2 + c) % n; }

void RhoFactorization::print_result(const Factor & factor)
{
    std::cout << "RESULT:\n";
    for (const auto & [div, deg] : factor)
    {
        std::cout << "Divisor " << div << " "
                  << "Degree " << deg << '\n';
    }
}

void RhoFactorization::notify_user(const long_int & div) { std::cout << "New divisor had been found " << div << '\n'; }

};  // namespace lpn
