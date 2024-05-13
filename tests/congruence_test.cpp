#include "congruence.hpp"
#include "boost/random.hpp"

#include <gtest/gtest.h>
#include <fstream>

namespace
{

using namespace lpn;  // NOLINT

bool Check(const long_int & x, const long_int & p, const size_t i) { return (x * x) % p == (i * i) % p; }

TEST(Basic, LegendreSymbol)
{
    ASSERT_EQ(ComputeLegendreSymbol(10, 5), 0);
    ASSERT_EQ(ComputeLegendreSymbol(1003, 1151), -1);
    ASSERT_EQ(ComputeLegendreSymbol(2, 7), 1);
    ASSERT_EQ(ComputeLegendreSymbol(32, 7), 1);
    ASSERT_EQ(ComputeLegendreSymbol(5, 7), -1);
    ASSERT_EQ(ComputeLegendreSymbol(144, 7), 1);
    ASSERT_EQ(ComputeLegendreSymbol(145, 7), -1);
    ASSERT_EQ(ComputeLegendreSymbol(453, 13), -1);
}

TEST(Basic, LegendreSymbolRandom)
{
    boost::mt19937 rng(42);
    {
        long_int p = 51439;
        for (size_t i = 1; i < 20000; ++i)
        {
            long_int value = p * rng() + i * i;
            ASSERT_TRUE(ComputeLegendreSymbol(value, p) == 1);
        }
    }

    {
        long_int p = 56287631;
        for (size_t i = 1; i < 20000; ++i)
        {
            long_int value = p * rng() + p - 1;
            ASSERT_TRUE(ComputeLegendreSymbol(value, p) == -1);
        }
    }

    {
        long_int p("850089435441534261939703688213");
        for (size_t i = 1; i < 20000; ++i)
        {
            long_int value = p * rng() + i * i;
            ASSERT_TRUE(ComputeLegendreSymbol(value, p) == 1);
        }
    }
}

TEST(Basic, LegendreSymbolRandomLarge)
{
    long_int value("1136576893960171706831330319665099465299");
    long_int p("706546294168342682401061556019");
    ASSERT_TRUE(ComputeLegendreSymbol(value, p) == 1);
}

TEST(Factorization, BasicPrime)
{
    long_int value("9111657031");
    auto factor = FactorizeBasic(value);
    ASSERT_TRUE(factor.size() == 1 && factor[value] == 1);
}

TEST(Sieve, QuadraticCongruencesBasic)
{
    long_int p = 41;
    long_int value = 64;
    ASSERT_TRUE(ComputeLegendreSymbol(value, p) == 1);
    ASSERT_TRUE(Check(QuadraticCongruences::Solve(value, p), p, 8));
}

TEST(Sieve, QuadraticCongruencesRandom)
{
    boost::mt19937 rng(42);
    {
        long_int p("706546294168342682401061556019");  // p % 4 == 3
        for (size_t i = 1; i < 20000; ++i)
        {
            long_int value = p * rng() + i * i;
            ASSERT_TRUE(ComputeLegendreSymbol(value, p) == 1);
            ASSERT_TRUE(Check(QuadraticCongruences::Solve(value, p), p, i));
        }
    }

    {
        long_int p("850089435441534261939703688213");  // p % 8 == 5
        for (size_t i = 1; i < 20000; ++i)
        {
            long_int value = p * rng() + i * i;
            ASSERT_TRUE(ComputeLegendreSymbol(value, p) == 1);
            ASSERT_TRUE(Check(QuadraticCongruences::Solve(value, p), p, i));
        }
    }

    {
        std::vector<long_int> primes = {128484928526173, 468655121175827, 343338799356211,
                                        383943794884159, 638267933343401, 455249100334669};
        for (long_int & p : primes)
        {
            for (size_t i = 1; i < 2000; ++i)
            {
                long_int value = p * rng() + i * i;
                ASSERT_TRUE(ComputeLegendreSymbol(value, p) == 1);
                ASSERT_TRUE(Check(QuadraticCongruences::Solve(value, p), p, i));
            }
        }
    }
}

TEST(Sieve, QuadraticCongruencesHeavy)
{
    std::ifstream file("primes.txt", std::ios_base::in);
    long_int p;
    boost::mt19937 rng(42);
    while (file >> p)
    {
        for (size_t i = 1'000'000; i < 1'002'000; ++i)
        {
            long_int value = p * rng() + i * i;
            ASSERT_TRUE(Check(QuadraticCongruences::Solve(value, p), p, i));
        }
    }
}

}  // namespace
