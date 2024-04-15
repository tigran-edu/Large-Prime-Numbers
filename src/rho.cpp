#include "rho.hpp"

namespace lpn
{

RhoFactorization::RhoFactorization()
    : c_(BasicConfig::kDefCValue),
      max_iter_(BasicConfig::kDefMaxIterValue),
      frequency_(BasicConfig::kDefFreqValue),
      max_attemp_(BasicConfig::kDefMaxAttempValue)
{
}

RhoFactorization::RhoFactorization(const long_int & c, size_t max_iter, size_t frequency, size_t max_attemp)
    : c_(c), max_iter_(max_iter), frequency_(frequency), max_attemp_(max_attemp)
{
}

FactorSet RhoFactorization::Factorize(const long_int & n, size_t starting_point)
{
    FactorSet factor;
    long_int x1 = starting_point;

    for (size_t i = 0; i < max_attemp_ && n > 1; ++i)
    {
        long_int divisor = FindDivisor(n, x1);
        if (divisor != 1 && divisor != n)
        {
            factor[divisor] += 1;
            factor[n / divisor] += 1;
            return factor;
        }
    }
    factor[n] = 1;
    return factor;
}

size_t RhoFactorization::ComputeNextRange(size_t range) const { return range << 1; }

void RhoFactorization::UpdateX1(const long_int & x2, long_int & x1) { x1 = x2; }

void RhoFactorization::UpdateX2(const long_int & n, size_t range, long_int & x2)
{
    for (size_t j = 0; j < range; ++j)
    {
        x2 = Next(x2, n);
    }
}

std::optional<long_int> RhoFactorization::TryToFindDivisor(const long_int & n, size_t range, size_t terms,
                                                           long_int & x1, long_int & x2)
{
    long_int product = 1;
    long_int divisor = 1;
    for (size_t j = 0; j < range; ++j)
    {
        x2 = Next(x2, n);
        product = (product * abs(x1 - x2)) % n;
        if (product == 0)
        {
            product = 1;
        }
        terms++;
        if ((terms % frequency_ == 0 || j + 1 == range) && ((divisor = lpn::gcd(n, product)) > 1))
        {
            return divisor;
        }
    }
    return std::nullopt;
}

long_int RhoFactorization::FindDivisor(const long_int & n, long_int & x1)
{
    long_int x2 = Next(x1, n);
    size_t range = 1;
    size_t terms = 0;
    while (terms <= max_iter_ && !IsPseudoPrime(n, BasicConfig::kPrimes))
    {
        auto divisor = TryToFindDivisor(n, range, terms, x1, x2);
        if (divisor.has_value())
        {
            return divisor.value();
        }
        terms += range;
        range = ComputeNextRange(range);
        UpdateX1(x2, x1);
        UpdateX2(n, range, x2);
    }
    return n;
}

long_int RhoFactorization::Next(const long_int & x2, const long_int & n) const { return (x2 * x2 + c_) % n; }

};  // namespace lpn
