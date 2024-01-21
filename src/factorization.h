#pragma once

#include <map>
#include "basic.h"

namespace lpn
{
using Factor = std::unordered_map<long_int, size_t>;

struct RhoFactorization
{
    RhoFactorization();
    RhoFactorization(const long_int & c, const size_t & max, const size_t & frequency);
    Factor factorize(long_int n, long_int x1);

private:
    void compute_diff(size_t & range, long_int & x1, long_int & x2);
    long_int next(const long_int & x2, const long_int & n);
    long_int compute_diff();
    size_t check_gcd(long_int & n, const long_int & g);
    long_int c;
    size_t max;
    size_t frequency;
};

Factor BasicFactorization(long_int a);
};
