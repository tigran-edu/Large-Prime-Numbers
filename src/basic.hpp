#pragma once

#include "aliases.hpp"
#include <optional>

namespace lpn
{

void MergeFactorSets(FactorSet & first, const FactorSet & second);

long_int Eval(const FactorSet & factor);

long_int FastExponentiation(long_int a, long_int b);
long_int FastExponentiationWithMod(long_int a, long_int b, const long_int & m);

size_t ExtractPower(long_int & a, const long_int & b);
size_t ExtractPowerFast(long_int & a, const long_int & b);

std::vector<bool> ToBinaryFormat(long_int val);

template <typename T>
T gcd(T a, T b)  // NOLINT
{
    while (b != 0)
    {
        a = std::exchange<T, T>(b, a % b);
    }
    return a;
}

template <typename Container>
bool IsPseudoPrime(const long_int & p, const Container & primes)
{
    for (const auto & prime : primes)
    {
        if (p % prime == 0 || FastExponentiationWithMod(prime, p - 1, p) != 1)
        {
            return false;
        }
    }
    return true;
}

template <typename Container>
bool IsStrongPseudoPrime(const long_int & p, const Container & primes)
{
    if (p % 2 == 0)
    {
        return false;
    }
    for (auto & prime : primes)
    {
        long_int t = p - 1;
        size_t a = ExtractPowerFast(t, 2);
        long_int test = FastExponentiationWithMod(prime, t, p);
        if (test == 1 || test == p - 1)
        {
            return true;
        }
        for (size_t i = 1; i < a; ++i)
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

template <typename Container>
static std::optional<FactorSet> TryToDecompose(const Container & primes, long_int n)
{
    FactorSet factor;
    for (size_t i = 0; i < primes.size() && n > 1; ++i)
    {
        if (n % primes[i] == 0)
        {
            factor[primes[i]] += ExtractPowerFast(n, primes[i]);
        }
    }
    return n == 1 ? std::optional<FactorSet>(factor) : std::nullopt;
}

FactorSet FactorizeBasic(long_int a);

bool IsPrimeBasic(const long_int & a);

long_int PowBasic(const long_int & a, long_int b);

};  // namespace lpn
