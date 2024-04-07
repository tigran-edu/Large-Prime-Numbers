#include "qs.hpp"
#include "boost/random.hpp"

#include <gtest/gtest.h>
#include <fstream>

namespace
{

using namespace lpn;  // NOLINT

bool Check(const long_int & x, const long_int & p, const size_t i) { return (x * x) % p == (i * i) % p; }

TEST(Sieve, QuadraticCongruencesBasic)
{
    long_int p = 41;
    long_int value = 64;
    ASSERT_TRUE(ComputeLegendreSymbol(value, p) == 1);
    ASSERT_TRUE(Check(QuadraticCongruences::SolvingQuadraticCongruences(value, p), p, 8));
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
            ASSERT_TRUE(Check(QuadraticCongruences::SolvingQuadraticCongruences(value, p), p, i));
        }
    }

    {
        long_int p("850089435441534261939703688213");  // p % 8 == 5
        for (size_t i = 1; i < 20000; ++i)
        {
            long_int value = p * rng() + i * i;
            ASSERT_TRUE(ComputeLegendreSymbol(value, p) == 1);
            ASSERT_TRUE(Check(QuadraticCongruences::SolvingQuadraticCongruences(value, p), p, i));
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
                ASSERT_TRUE(Check(QuadraticCongruences::SolvingQuadraticCongruences(value, p), p, i));
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
            ASSERT_TRUE(Check(QuadraticCongruences::SolvingQuadraticCongruences(value, p), p, i));
        }
    }
}

TEST(Sieve, QuadraticSieve)
{
    long_int n("59469489332848408438249254427481121839977");  // 338555568168236555657 * 175656509371887105761
    FactorSet factor = QuadraticSieve::Factorize(n);
    ASSERT_EQ(factor.size(), 2);
    ASSERT_EQ(Eval(factor), n);
}

TEST(Sieve, QuadraticSieveWithConfig)
{
    long_int n("4482406424966880742829846540605971439398287609");  // 86738535685150523290199 * 51677220390685710220591
    auto config = Sieve::CreateConfig(n, 500'000'000, 5000, 1);
    FactorSet factor = QuadraticSieve::Factorize(n, config);
    ASSERT_EQ(factor.size(), 2);
    ASSERT_EQ(Eval(factor), n);
}

};  // namespace
