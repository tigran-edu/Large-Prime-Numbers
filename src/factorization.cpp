#include "factorization.h"
#include <iostream>
#include "boost/math/ccmath/abs.hpp"


namespace lpn
{

std::unordered_map<cpp_int, size_t> RhoFactorization(cpp_int n, cpp_int c, cpp_int max)
{
    std::unordered_map<cpp_int, size_t> factor;
    if (PseudoprimeTest(n) == 1)
    {
        std::cerr << "Warning: Value " << n << " is Pseudoprime\n";
        factor[n] += 1;
        return factor;
    }
    cpp_int x1 = 2;
    cpp_int x2 = 4 + c;
    cpp_int range = 1;
    cpp_int product = 1;
    cpp_int terms = 0;

    while (terms <= max && PseudoprimeTest(n) != 1)
    {
        for (cpp_int j = 0; j < range; ++j)
        {
            x2 = (x2 * x2 + c) % n;
            product = (product * abs(x1 - x2)) % n;
            terms++;
            if (terms % 10 == 0 || j + 1 == range)
            {
                auto g = gcd(n, product);
                if (g > 1)
                {
                    while (n % g == 0)
                    {
                        factor[g]++;
                        n /= g;
                    }
                }
                product = 1;
            }
        }
        x1 = x2;
        range = 2 * range;
        for (cpp_int j = 0; j < range; ++j)
            x2 = (x2 * x2 + c) % n;
    }
    if (n > 1)
        factor[n]++;
    return factor;
}


};