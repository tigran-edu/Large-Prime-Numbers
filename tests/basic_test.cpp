#include <boost/math/special_functions/pow.hpp>
#include <gtest/gtest.h>
#include "basic.h"

namespace
{

using long_int = boost::multiprecision::cpp_int;

long_int basicPow(const long_int & a, long_int b)
{
    long_int answer = 1;
    while (b > 0)
    {
        answer *= a;
        --b;
    }
    return answer;
}

long_int basicPowWithMod(const long_int & a, const long_int & b, const long_int & mod)
{
    return lpn::FastExponentiation(a, b) % mod;
}

TEST(Basic, fastExp)
{
    for (long_int b = 2; b < 100; ++b)
        for (long_int a = 2; a < 100; ++a)
            ASSERT_EQ(lpn::FastExponentiation(a, b), basicPow(a, b));
}

TEST(Basic, fastExpWithMod)
{
    long_int mod = 41887;
    long_int a = 389;
    long_int b = 563;

    for (size_t i = 1; i <= 20; ++i)
    {
        ASSERT_EQ(lpn::FastExponentiationWithMod(a, b, mod), basicPowWithMod(a, b, mod));
        mod *= mod + 1;
        a += i;
        b += mod % i;
    }
}

TEST(Basic, fullDiv)
{
    long_int a = lpn::FastExponentiation(2, 10);
    ASSERT_EQ(lpn::FullDiv(a, 2), 10);
    ASSERT_EQ(a, 1);
}

TEST(Basic, gcd)
{
    ASSERT_EQ(lpn::gcd(1, 10), 1);
    ASSERT_EQ(lpn::gcd(2, 10), 2);
    ASSERT_EQ(lpn::gcd(119, 3456712), 119);
    ASSERT_EQ(lpn::gcd(34567890, 8827), 679);
}

TEST(Basic, PseudoPrimeTest)
{
    const std::array<long_int, 6> primes = {2, 3, 5, 7, 11, 13};
    // Large Prime Number
    long_int value("987011932680753357835689163067");
    ASSERT_EQ(lpn::PseudoPrimeTest(value, primes), true);

    value = 561; // 561 = 3 x 11 x 17 and 561 is a Carmichael number.
    ASSERT_EQ(lpn::PseudoPrimeTest(value, primes), true);

    value = 341;
    ASSERT_EQ(lpn::PseudoPrimeTest(value, primes), false);
}

};
