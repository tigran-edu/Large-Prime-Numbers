#pragma once

#include <map>
#include <vector>

#include "basic.hpp"

namespace lpn
{

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

};  // namespace lpn
