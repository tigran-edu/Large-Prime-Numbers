#include "congruence.hpp"

#include <climits>
#include <random>

namespace lpn
{

int ComputeLegendreSymbol(long_int n, long_int p)
{
    int legendre = 1;
    n = n % p;

    if (n == 0)
    {
        return 0;
    }

    if (n < 0)
    {
        n *= -1;
        if (p % 4 == 3)
        {
            legendre *= -1;
        }
    }

    while (n > 1)
    {
        size_t deg = ExtractPower(n, 2);
        if (deg % 2 == 1 && (p * p - 1) % 16 == 8)
        {
            legendre *= -1;
        }
        if ((n - 1) * (p - 1) % 8 == 4)
        {
            legendre *= -1;
        }
        p = std::exchange<long_int, long_int>(n, p % n);
    }
    return legendre;
}

long_int QuadraticCongruences::Solve(const long_int & n, const long_int & p)
{
    assert(p % 2 != 0);

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
    return (SolveCongruence(n, h, p) * (p + 1) / 2) % p;
}

long_int QuadraticCongruences::FindStartValue(const long_int & n, const long_int & p)
{
    static std::mt19937 gen{42};
    static std::uniform_int_distribution<size_t> distrib(1, UINT64_MAX);
    long_int h = 1;
    while (ComputeLegendreSymbol(h * h - 4 * n, p) != -1)
    {
        h = distrib(gen);
    }
    return h;
}

long_int QuadraticCongruences::SolveCongruence(const long_int & n, const long_int & h, const long_int & p)
{
    long_int m = n;
    long_int v = h;
    long_int w = (h * h - 2 * n) % p;
    auto bitset = ToBinaryFormat((p + 1) / 2);

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

std::vector<bool> QuadraticCongruences::ToBinaryFormat(long_int val)
{
    std::vector<bool> bitset;
    while (val > 0)
    {
        bitset.push_back(val % 2 == 0);
        val >>= 1;
    }
    std::reverse(bitset.begin(), bitset.end());
    return bitset;
}

};  // namespace lpn
