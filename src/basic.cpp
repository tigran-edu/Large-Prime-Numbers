#include "basic.hpp"

namespace lpn
{

void MergeFactors(Factor & first, const Factor & second)
{
    for (auto & [key, value] : second)
    {
        first[key] += value;
    }
}

long_int Eval(const Factor & factor)
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

size_t FullDiv(long_int & a, const long_int & b)
{
    if (b < 2)
    {
        return 0;
    }
    size_t counter = 0;
    while (a % b == 0)
    {
        a /= b;
        counter++;
    }
    return counter;
}

size_t FullDivFast(long_int & a, const long_int & b, size_t i)
{
    if (b < 2 || a % b != 0)
    {
        return 0;
    }

    size_t counter = FullDivFast(a, b * b, 2 * i);
    if (a % b == 0)
    {
        counter += i;
        a /= b;
    }
    return counter;
}

int LegendreSymbol::Compute(long_int n, long_int p)
{
    int legendre = 1;
    n = n % p;

    if (n == 0)
    {
        return 0;
    }

    if (n < 0)
    {
        n *= -1;
        if (p % 4 == 3)
        {
            legendre *= -1;
        }
    }

    PullTwos(n, legendre, p);

    while (n > 1)
    {
        if ((n - 1) * (p - 1) % 8 == 4)
        {
            legendre *= -1;
        }
        p = std::exchange<long_int, long_int>(n, p % n);
        PullTwos(n, legendre, p);
    }
    return legendre;
}

void LegendreSymbol::PullTwos(long_int & n, int & legendre, const long_int & p)
{
    int count = 0;
    while (n % 2 == 0)
    {
        n /= 2;
        count = 1 - count;
    }
    if ((p * p - 1) * count % 16 == 8)
    {
        legendre *= -1;
    }
}

Factor BasicFactorization(long_int a)
{
    Factor factor;
    if (auto counter = FullDiv(a, 2); counter > 0)
    {
        factor[2] += counter;
    }

    long_int div = 3;
    while (div * div <= a)
    {
        if (auto counter = FullDiv(a, div); counter > 0)
        {
            factor[div] += counter;
        }
        div += 2;
    }
    if (a > 1)
    {
        factor[a]++;
    }
    return factor;
}

bool BasicIsPrime(const long_int & a)
{
    if (a <= 3)
    {
        return true;
    }
    if (a % 2 == 0)
    {
        return false;
    }
    long_int div = 3;
    while (div * div <= a)
    {
        if (a % div == 0)
        {
            return false;
        }
        div++;
    }

    return true;
}

};  // namespace lpn
