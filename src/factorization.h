#pragma once

#include <map>
#include <vector>

#include "basic.h"

namespace lpn
{
using Factor = std::unordered_map<long_int, size_t>;

class RhoFactorization
{
   public:
    RhoFactorization();
    RhoFactorization(const long_int & c, const size_t & max, const size_t & frequency);
    Factor Factorize(long_int n, long_int x1);

   private:
    long_int FindNewDivisor(long_int & n, long_int & x1);
    long_int Next(const long_int & x2, const long_int & n);
    void PrintResult(const Factor & factor);
    void NotifyUser(const long_int & div);
    long_int c_;
    size_t max_iter_amount_;
    size_t frequency_;
};

class QuadraticSieve
{
   public:
    static Factor Factorize(long_int n);
};

class QuadraticCongruences
{
   public:
    static long_int SolvingQuadraticCongruences(const long_int & n, const long_int & p);

   private:
    static long_int Solve(const long_int & n, const long_int & h, const long_int & p);
    static long_int FindStartValue(const long_int & n, const long_int & p);
    static std::vector<bool> Binary(long_int val);
};

Factor BasicFactorization(long_int a);
};  // namespace lpn
