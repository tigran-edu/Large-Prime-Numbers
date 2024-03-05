#pragma once

#include <iostream>
#include <map>

#include "boost/multiprecision/cpp_int.hpp"

namespace lpn
{
using long_int = boost::multiprecision::cpp_int;  // NOLINT

long_int FastExponentiation(long_int a, long_int b);
long_int FastExponentiationWithMod(long_int a, long_int b, const long_int & m);

size_t FullDiv(long_int & a, const long_int & b);

template <typename T>
T gcd(T a, T b)  // NOLINT
{
    while (b != 0)
    {
        a = std::exchange<T, T>(b, a % b);
    }
    return a;
}

template <size_t N>
bool PseudoPrimeTest(const long_int & p, const std::array<long_int, N> & primes)
{
    for (auto & prime : primes)
    {
        if (p % prime == 0 || FastExponentiationWithMod(prime, p - 1, p) != 1)
        {
            return false;
        }
    }
    return true;
}

template <size_t N>
bool StrongPseudoPrimeTest(const long_int & p, const std::array<long_int, N> & primes)
{
    if (p % 2 == 0)
    {
        return false;
    }
    for (auto & prime : primes)
    {
        long_int t = p - 1;
        size_t a = 0;
        while (t % 2 == 0)
        {
            t /= 2;
            a += 1;
        }
        long_int test = FastExponentiationWithMod(prime, t, p);
        if (test == 1 || test == p - 1)
        {
            return true;
        }
        for (size_t i = 1; i < a - 1; ++i)
        {
            test = (test * test) % p;
            if (test == p - 1)
            {
                return true;
            }
        }
    }
    return false;
}

class LegendreSymbol
{
   public:
    static int Compute(long_int n, long_int p);

   private:
    static void PullTwos(long_int & n, int & legendre, const long_int & p);
};

};  // namespace lpn
