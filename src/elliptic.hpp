#pragma once

#include "aliases.hpp"

namespace lpn
{

class EllipticCurveFactorization
{
   public:
    static FactorSet Factorize(const long_int & n);
    static FactorSet Factorize(const long_int & n, const long_int & a, const long_int & x, const long_int & y,
                               size_t max_iter);

   private:
    EllipticCurveFactorization(long_int a, long_int b, long_int x, long_int z);
    static long_int ComputeB(const long_int & x, const long_int & y, const long_int & a, const long_int & n);

   private:
    long_int ComputeX(const long_int & r, const long_int & s, const long_int & n) const;
    long_int ComputeZ(const long_int & r, const long_int & s, const long_int & n) const;
    long_int ComputeU(const long_int & r, const long_int & s, const long_int & u, const long_int & v,
                      const long_int & n) const;
    long_int ComputeW(const long_int & r, const long_int & s, const long_int & u, const long_int & v,
                      const long_int & n) const;
    void ComputeNextValues(size_t k, const long_int & n);
    FactorSet FindFactor(const long_int & n, size_t max_iter);

   private:
    const long_int a_;
    const long_int b_;
    long_int x_;
    long_int z_;
};

};  // namespace lpn
