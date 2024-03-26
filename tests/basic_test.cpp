#include "basic.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/pow.hpp>
#include <gtest/gtest.h>

namespace
{

using namespace lpn;  // NOLINT

long_int BasicPow(const long_int & a, long_int b)
{
    long_int answer = 1;
    while (b > 0)
    {
        answer *= a;
        --b;
    }
    return answer;
}

long_int BasicPowWithMod(const long_int & a, const long_int & b, const long_int & mod)
{
    return FastExponentiation(a, b) % mod;
}

TEST(Basic, fastExp)
{
    for (long_int b = 2; b < 100; ++b)
    {
        for (long_int a = 2; a < 100; ++a)
        {
            ASSERT_EQ(FastExponentiation(a, b), BasicPow(a, b));
        }
    }
}

TEST(Basic, fastExpWithMod)
{
    long_int mod = 41887;
    long_int a = 389;
    long_int b = 563;

    for (size_t i = 1; i <= 20; ++i)
    {
        ASSERT_EQ(FastExponentiationWithMod(a, b, mod), BasicPowWithMod(a, b, mod));
        mod *= mod + 1;
        a += i;
        b += mod % i;
    }
}

TEST(Basic, fullDiv)
{
    for (size_t i = 2; i < 5; ++i)
    {
        for (size_t j = 10; j < 32; ++j)
        {
            long_int a = FastExponentiation(i, j);
            ASSERT_EQ(FullDiv(a, i), j);
            ASSERT_EQ(a, 1);
        }
    }
}

TEST(Basic, fullDivFast)
{
    for (size_t i = 2; i < 5; ++i)
    {
        for (size_t j = 10; j < 32; ++j)
        {
            long_int a = FastExponentiation(i, j);
            ASSERT_EQ(FullDivFast(a, i), j);
            ASSERT_EQ(a, 1);
        }
    }
}

TEST(Basic, IsPrime)
{
    ASSERT_TRUE(BasicIsPrime(2));
    ASSERT_TRUE(BasicIsPrime(3));
    ASSERT_TRUE(BasicIsPrime(5));
    ASSERT_TRUE(BasicIsPrime(7));
    ASSERT_FALSE(BasicIsPrime(128));
    ASSERT_FALSE(BasicIsPrime(119));
    ASSERT_FALSE(BasicIsPrime(282112));
}

TEST(Basic, gcd)
{
    ASSERT_EQ(gcd(1, 10), 1);
    ASSERT_EQ(gcd(2, 10), 2);
    ASSERT_EQ(gcd(119, 3456712), 119);
    ASSERT_EQ(gcd(34567890, 8827), 679);
}

TEST(Basic, PseudoPrimeTest)
{
    const std::array<long_int, 5> primes = {2, 3, 5, 7, 11};
    // Large Prime Number
    long_int value("987011932680753357835689163067");
    ASSERT_EQ(PseudoPrimeTest(value, primes), true);

    value = 561;  // 561 = 3 x 11 x 17 and 561 is a Carmichael number.
    ASSERT_EQ(PseudoPrimeTest(value, primes), false);

    value = 341;
    ASSERT_EQ(PseudoPrimeTest(value, primes), false);
}

TEST(Basic, StrongPseudoPrimeTest)
{
    const std::array<long_int, 4> primes = {2, 3, 5, 7};
    // Large Prime Number
    long_int value("55850466220760803896422741537242561024552133147337966990409402967679250928454649");
    ASSERT_EQ(StrongPseudoPrimeTest(value, primes), true);

    value += 1;  // value % 2 == 0
    ASSERT_EQ(StrongPseudoPrimeTest(value, primes), false);

    value = 3215031751;  // 3215031751 = 151 x 751 x 28351.
    ASSERT_EQ(StrongPseudoPrimeTest(value, primes), true);

    value = 1;
    for (size_t i = 1; i < 100; ++i)
    {
        value *= i;
        value++;
    }
    ASSERT_EQ(StrongPseudoPrimeTest(value, primes), false);
}

TEST(Basic, LegendreSymbol)
{
    ASSERT_EQ(LegendreSymbol::Compute(10, 5), 0);
    ASSERT_EQ(LegendreSymbol::Compute(1003, 1151), -1);
    ASSERT_EQ(LegendreSymbol::Compute(2, 7), 1);
    ASSERT_EQ(LegendreSymbol::Compute(32, 7), 1);
    ASSERT_EQ(LegendreSymbol::Compute(5, 7), -1);
    ASSERT_EQ(LegendreSymbol::Compute(144, 7), 1);
    ASSERT_EQ(LegendreSymbol::Compute(145, 7), -1);
    ASSERT_EQ(LegendreSymbol::Compute(453, 13), -1);
}

TEST(Factorization, BasicPrime)
{
    long_int value("9111657031");
    auto factor = BasicFactorization(value);
    ASSERT_TRUE(factor.size() == 1 && factor[value] == 1);
}

TEST(Factorization, BasicComplex)
{
    long_int value("1307674368000");
    auto factor = BasicFactorization(value);
    ASSERT_TRUE(factor.size() == 6);
    long_int result = 1;
    for (const auto & [div, amount] : factor)
    {
        result *= FastExponentiation(div, amount);
    }
    ASSERT_TRUE(result == value);
}

};  // namespace
