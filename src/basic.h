#pragma once

#include <map>
#include "boost/multiprecision/cpp_int.hpp"

namespace lpn
{
using long_int = boost::multiprecision::cpp_int;

long_int FastExponentiation(long_int a, long_int b);
long_int FastExponentiationWithMod(long_int a, long_int b, const long_int & m);

size_t FullDiv(long_int & a, const long_int & b);

template <typename T>
T gcd(T a, T b)
{
    while (b != 0)
        a = std::exchange(b, a % b);
    return a;
}

template <size_t N>
bool PseudoPrimeTest(const long_int & p, const std::array<long_int, N> & primes)
{
    for (auto & prime : primes)
        if (lpn::gcd(p, prime) == 1 && FastExponentiationWithMod(prime, p - 1, p) != 1)
            return false;
    return true;
}
};
