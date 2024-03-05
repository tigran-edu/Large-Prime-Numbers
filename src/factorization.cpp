#include "factorization.h"
#include "basic.h"

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

long_int QuadraticCongruences::SolvingQuadraticCongruences(const long_int & n, const long_int & p)
{
    if (p % 4 == 3)
    {
        return FastExponentiationWithMod(n, (p + 1) / 4, p);
    }
    if (p % 8 == 5 && FastExponentiationWithMod(n, (p - 1) / 4, p) == 1)
    {
        return FastExponentiationWithMod(n, (p + 3) / 8, p);
    }
    if (p % 8 == 5 && FastExponentiationWithMod(n, (p - 1) / 4, p) == p - 1)
    {
        return (FastExponentiationWithMod(4 * n, (p + 3) / 8, p) * (p + 1) / 2) % p;
    }
    long_int h = FindStartValue(n, p);
    return (Solve(n, h, p) * (p + 1) / 2) % p;
}

long_int QuadraticCongruences::FindStartValue(const long_int & n, const long_int & p)
{
    boost::random::mt19937 gen{42};
    while (true)
    {
        long_int h = gen();
        if (LegendreSymbol::Compute(h * h - 4 * n, p) == -1)
        {
            return h;
        }
    }
}

long_int QuadraticCongruences::Solve(const long_int & n, const long_int & h, const long_int & p)
{
    long_int m = n;
    long_int v = h;
    long_int w = (h * h - 2 * n) % p;
    auto bitset = Binary((p + 1) / 2);

    for (size_t i = 1; i < bitset.size(); ++i)
    {
        long_int x = (v * w - h * m) % p;
        v = (v * v - 2 * m) % p;
        w = (w * w - 2 * n * m) % p;
        m = (m * m) % p;
        if (bitset[i])
        {
            w = x;
        }
        else
        {
            v = x;
            m = (n * m) % p;
        }
    }
    return (p + v) % p;
}

std::vector<bool> QuadraticCongruences::Binary(long_int val)
{
    std::vector<bool> bitset;
    while (val > 0)
    {
        bitset.push_back(val % 2 == 0);
        val /= 2;
    }
    std::reverse(bitset.begin(), bitset.end());
    return std::move(bitset);
}

};  // namespace lpn
