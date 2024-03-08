#include "qs.hpp"
#include <boost/random.hpp>
#include <stdexcept>

namespace lpn
{

long_int QuadraticCongruences::SolvingQuadraticCongruences(const long_int & n, const long_int & p)
{
    if (p % 2 == 0)
    {
        throw std::runtime_error(
            "SolvingQuadraticCongruences available only for odd prime numbers. "
            "Entered value: " +
            p.str());
    }

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

std::vector<long_int> Sieve::Sieving(long_int n)
{
    Config cf(n);
    std::vector<long_int> potential_divs;
    std::vector<cpp_bin_float_100> logs(2 * cf.m, 0);
    long_int r = boost::multiprecision::sqrt(n) - cf.m;

    for (size_t pos = 0; pos < cf.primes.size(); ++pos)
    {
        size_t p = cf.primes[pos];
        // n = x^2 mod p; f(r) = r^2 - n; target f(r) == 0 mod p
        size_t i = (size_t)((cf.congruences[pos] + p - (r % p)) % p);
        Fill(i, p, logs);
        size_t j = (size_t)((2 * p - cf.congruences[pos] - (r % p)) % p);
        if (j != i)
        {
            Fill(i, p, logs);
        }
    }

    for (size_t i = 0; i < logs.size(); ++i)
    {
        if (logs[i] > cf.closeness && IsDecomposed(cf, i, n))
        {
            potential_divs.push_back(Eval(r, n, i));
        }
    }
    return potential_divs;
}

void Sieve::Fill(size_t i, size_t p, std::vector<cpp_bin_float_100> & logs)
{
    cpp_bin_float_100 p_log = log10(cpp_bin_float_100(p));
    while (i < logs.size())
    {
        logs[i] += p_log;
        i += p;
    }
}

bool Sieve::IsDecomposed(Config & cf, size_t i, const long_int & n)
{
    long_int r = boost::multiprecision::sqrt(n) - cf.m;
    long_int value = Eval(r, n, i);
    for (auto & p : cf.primes)
    {
        FullDiv(value, p);
    }
    return value == 1;
}

long_int Sieve::Eval(const long_int & r, const long_int & n, size_t i)
{
    return boost::multiprecision::abs((r + i) * (r + i) - n);
}

};  // namespace lpn
