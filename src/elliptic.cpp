#include "elliptic.hpp"
#include "aliases.hpp"
#include "basic.hpp"

namespace lpn
{

EllipticCurveFactorization::EllipticCurveFactorization(long_int a, long_int b, long_int x, long_int z)
    : a_(a), b_(b), x_(x), z_(z)
{
}

long_int EllipticCurveFactorization::ComputeB(const long_int & x, const long_int & y, const long_int & a,
                                              const long_int & n)
{
    return (n + (y * y - x * (x * x - a)) % n) % n;
}

FactorSet EllipticCurveFactorization::Factorize(const long_int & n, const long_int & a, const long_int & x,
                                                const long_int & y, size_t max_iter)
{
    return EllipticCurveFactorization(a, ComputeB(x, y, a, n), x, 1).FindFactor(n, max_iter);
}

long_int EllipticCurveFactorization::ComputeX(const long_int & r, const long_int & s, const long_int & n) const
{
    long_int s2 = s * s;
    long_int term = (r * r - a_ * s2) % n;
    return (n + (term * term - 8 * b_ * r * s2 * s) % n) % n;
}

long_int EllipticCurveFactorization::ComputeZ(const long_int & r, const long_int & s, const long_int & n) const
{
    return (n + (4 * s * (r * r * r + s * s * (a_ * r + b_ * s))) % n) % n;
}

long_int EllipticCurveFactorization::ComputeU(const long_int & r, const long_int & s, const long_int & u,
                                              const long_int & v, const long_int & n) const
{
    long_int term1 = (r * u - a_ * s * v) % n;
    long_int term2 = (b_ * s * v * (r * v + s * u)) % n;
    return (n + (z_ * (term1 * term1 - 4 * term2)) % n) % n;
}

long_int EllipticCurveFactorization::ComputeW(const long_int & r, const long_int & s, const long_int & u,
                                              const long_int & v, const long_int & n) const
{
    long_int term = (u * s - r * v) % n;
    return (x_ * term * term) % n;
}

FactorSet EllipticCurveFactorization::FindFactor(const long_int & n, size_t max_iter)
{
    FactorSet factor;
    long_int g = lpn::gcd<long_int>(4 * a_ * a_ * a_ + 27 * b_ * b_, n);
    if (g != 1)
    {
        factor[g] = 1;
        factor[n / g] = 1;
        return factor;
    }
    for (size_t k = 2; k <= max_iter;)
    {
        for (size_t i = 0; i < 10; ++i, ++k)
        {
            ComputeNextValues(k, n);
        }
        g = lpn::gcd<long_int>(z_, n);
        if (g != 1)
        {
            factor[g] = 1;
            factor[n / g] = 1;
            return factor;
        }
    }

    factor[n] = 1;
    return factor;
}

void EllipticCurveFactorization::ComputeNextValues(size_t k, const long_int & n)
{
    auto bitset = ToBinaryFormat(k);
    long_int x1 = x_;
    long_int z1 = z_;
    long_int x2 = ComputeX(x_, z_, n);
    long_int z2 = ComputeZ(x_, z_, n);
    for (size_t i = 1; i < bitset.size(); ++i)
    {
        auto u = ComputeU(x1, z1, x2, z2, n);
        auto w = ComputeW(x1, z1, x2, z2, n);
        if (!bitset[i])
        {
            auto tmp = ComputeX(x1, z1, n);
            z1 = ComputeZ(x1, z1, n);
            x1 = std::move(tmp);
            x2 = std::move(u);
            z2 = std::move(w);
        }
        else
        {
            auto tmp = ComputeX(x2, z2, n);
            z2 = ComputeZ(x2, z2, n);
            x2 = std::move(tmp);
            x1 = std::move(u);
            z1 = std::move(w);
        }
    }
    x_ = std::move(x1);
    z_ = std::move(z1);
}

};  // namespace lpn
