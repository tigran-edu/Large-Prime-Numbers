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
    long_int find_new_divisor(long_int & n, long_int & x1);
    long_int next(const long_int & x2, const long_int & n);
    void print_result(const Factor & factor);
    void notify_user(const long_int & div);
    long_int c;
    size_t max_iter_amount;
    size_t frequency;
};

Factor BasicFactorization(long_int a);
};
