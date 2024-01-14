#include "basic.h"


namespace lpn
{

std::vector<cpp_int> primes = {2, 3, 5, 7, 11, 13};


std::unordered_map<cpp_int, size_t> BasicFactorization(cpp_int a)
{
    std::unordered_map<cpp_int, size_t> factor;
    cpp_int div = 2;
    while (div * div <= a)
    {
        if (a % div == 0)
        {
            a /= div;
            factor[div]++;
        }
        else
        {
            div++;
        }
    }
    if (a > 1)
        factor[a]++;
    return factor;
}


cpp_int gcd(cpp_int a, cpp_int b)
{
    while (b != 0)
    {
        auto tmp = b;
        b = a % b;
        a = tmp;
    }
    return a;
}


cpp_int FastExponentiationWithMod(cpp_int a, cpp_int b, cpp_int m)
{
    cpp_int n = 1;
    while (b != 0)
    {
        if (b % 2 != 0)
            n = (n * a) % m;
        b = b / 2;
        a = (a * a) % m;
    }
    return n;
}

cpp_int FastExponentiation(cpp_int a, cpp_int b)
{
    cpp_int n = 1;
    while (b != 0)
    {
        if (b % 2 != 0)
            n = n * a;
        b = b / 2;
        a = (a * a);
    }
    return n;
}


bool PseudoprimeTest(cpp_int p)
{
    for (auto & prime : primes)
        if (gcd(p, prime) == 1 && FastExponentiationWithMod(prime, p - 1, p) != 1)
            return false;
    return true;
}
};