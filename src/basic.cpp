#include "basic.hpp"
#include <cassert>

namespace lpn
{

void MergeFactorSets(FactorSet & first, const FactorSet & second)
{
    for (const auto & [key, value] : second)
    {
        first[key] += value;
    }
}

long_int Eval(const FactorSet & factor)
{
    long_int result = 1;
    for (const auto & [div, amount] : factor)
    {
        result *= FastExponentiation(div, amount);
    }
    return result;
}

long_int FastExponentiationWithMod(long_int a, long_int b, const long_int & m)
{
    assert(b > 0);

    long_int n = 1;
    while (b != 0)
    {
        if (b % 2 != 0)
        {
            n = (n * a) % m;
        }
        b = b / 2;
        a = (a * a) % m;
    }
    return n;
}

long_int FastExponentiation(long_int a, long_int b)
{
    assert(b > 0);

    long_int n = 1;
    while (b != 0)
    {
        if (b % 2 != 0)
        {
            n = n * a;
        }
        b = b / 2;
        a = (a * a);
    }
    return n;
}

size_t ExtractPower(long_int & a, const long_int & b)
{
    assert(b > 1);

    size_t power = 0;
    while (a % b == 0)
    {
        a /= b;
        power++;
    }
    return power;
}

size_t ExtractPowerFast(long_int & a, const long_int & b)
{
    assert(b > 1);

    if (a % b != 0)
    {
        return 0;
    }

    size_t power = ExtractPowerFast(a, b * b) * 2;
    if (a % b == 0)
    {
        power += 1;
        a /= b;
    }
    return power;
}

FactorSet FactorizeBasic(long_int a)
{
    FactorSet factor;
    if (auto counter = ExtractPowerFast(a, 2); counter > 0)
    {
        factor[2] += counter;
    }

    long_int sqrt_a = boost::multiprecision::sqrt(a);
    for (long_int div = 3; div <= sqrt_a; div += 2)
    {
        if (auto counter = ExtractPowerFast(a, div); counter > 0)
        {
            factor[div] += counter;
        }
    }
    if (a > 1)
    {
        factor[a]++;
    }
    return factor;
}

bool IsPrimeBasic(const long_int & a)
{
    if (a <= 3)
    {
        return true;
    }
    if (a % 2 == 0)
    {
        return false;
    }
    long_int sqrt_a = boost::multiprecision::sqrt(a);
    for (long_int div = 3; div <= sqrt_a; div += 2)
    {
        if (a % div == 0)
        {
            return false;
        }
    }

    return true;
}

long_int PowBasic(const long_int & a, long_int b)
{
    long_int answer = 1;
    while (b > 0)
    {
        answer *= a;
        --b;
    }
    return answer;
}

};  // namespace lpn
