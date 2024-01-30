#include "basic.h"

namespace lpn
{
long_int FastExponentiationWithMod(long_int a, long_int b, const long_int & m)
{
    long_int n = 1;
    while (b != 0)
    {
        if (b % 2 != 0)
            n = (n * a) % m;
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
            n = n * a;
        b = b / 2;
        a = (a * a);
    }
    return n;
}

size_t FullDiv(long_int & a, const long_int & b)
{
    if (b < 2)
        return 0;
    size_t counter = 0;
    while (a % b == 0)
    {
        a /= b;
        counter++;
    }
    return counter;
}
};
