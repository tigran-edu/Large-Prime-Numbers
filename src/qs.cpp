#include "qs.hpp"
#include <boost/random.hpp>
#include "basic.hpp"
#include <stdexcept>

namespace lpn
{

FactorSet QuadraticSieve::Factorize(const long_int & n)
{
    FactorSet factor;
    SieveResult result = Sieve::Sieving(n);
    GaussianBasic gs = GaussianBasic(result.factors, result.primes);
    auto solutions = gs.Solve();

    for (auto & sol : solutions)
    {
        if (sol.mask.none())
        {
            long_int check = CheckResults(result, sol, n);
            if (check != 1 && check != n)
            {
                factor[check] = 1;
                factor[n / check] = 1;
                return factor;
            }
        }
    }
    return factor;
}

long_int QuadraticSieve::CheckResults(SieveResult & result, const GaussianBasic::Bitset & solution, const long_int & n)
{
    long_int x = 1;
    FactorSet factor;
    for (size_t i = 0; i < result.factors.size(); ++i)
    {
        if (solution.participants[i])
        {
            x = (x * (result.r + result.positions[i])) % n;
            MergeFactorSets(factor, result.factors[i]);
        }
    }
    long_int y = 1;
    for (auto & [key, value] : factor)
    {
        assert(value % 2 == 0);
        while (value > 0)
        {
            y = (y * key) % n;
            value -= 2;
        }
    }
    return lpn::gcd<long_int>(boost::multiprecision::abs(x - y), n);
}

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
        if (ComputeLegendreSymbol(h * h - 4 * n, p) == -1)
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

SieveResult Sieve::Sieving(const long_int & n)
{
    SieveResult result;
    Config cf(n);
    std::vector<float> logs(2 * cf.m, 0.0);
    long_int r = boost::multiprecision::sqrt(n) - cf.m;
    result.primes = cf.primes;
    result.r = r;

    for (size_t pos = 0; pos < cf.primes.size(); ++pos)
    {
        size_t p = cf.primes[pos];
        // n = x^2 mod p; f(r) = r^2 - n; target f(r) == 0 mod p
        size_t i = (size_t)((cf.congruences[pos] + p - (r % p)) % p);
        size_t j = (size_t)((2 * p - cf.congruences[pos] - (r % p)) % p);
        CacheSaveFill(i, j, p, logs);
    }

    for (size_t i = 0; i < logs.size(); ++i)
    {
        if (cf.closeness <= logs[i] && logs[i] <= cf.target && result.positions.size() < 1.15 * cf.factor_size)
        {
            auto factor = IsDecomposed(cf, i, n);
            if (factor != std::nullopt)
            {
                result.positions.push_back(i);
                result.factors.push_back(std::move(factor.value()));
            }
        }
    }
    return result;
}

void Sieve::Fill(size_t i, size_t p, std::vector<float> & logs)
{
    float p_log = log10(float(p));
    while (i < logs.size())
    {
        logs[i] += p_log;
        i += p;
    }
}

void Sieve::CacheSaveFill(size_t i, size_t j, size_t p, std::vector<float> & logs)
{
    if (j < i)
    {
        std::swap(i, j);
    }
    else if (i == j)
    {
        Fill(i, p, logs);
        return;
    }

    float p_log = log10(float(p));
    while (j < logs.size())
    {
        logs[i] += p_log;
        logs[j] += p_log;
        i += p;
        j += p;
    }
    if (i < logs.size())
    {
        logs[i] += p_log;
    }
}

std::optional<FactorSet> Sieve::IsDecomposed(Config & cf, size_t i, const long_int & n)
{
    long_int r = boost::multiprecision::sqrt(n) - cf.m;
    long_int value = F(r, n, i);
    FactorSet factor;
    for (auto p : cf.primes)
    {
        if (value % p == 0)
        {
            factor[p] += ExtractPowerFast(value, p);
        }
    }
    return value == 1 ? std::optional<FactorSet>(factor) : std::nullopt;
}

long_int Sieve::F(const long_int & r, const long_int & n, size_t i)
{
    return boost::multiprecision::abs((r + i) * (r + i) - n);
}

};  // namespace lpn
