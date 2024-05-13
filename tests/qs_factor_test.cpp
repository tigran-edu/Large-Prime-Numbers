#include "qs.hpp"

#include <gtest/gtest.h>

namespace
{

using namespace lpn;  // NOLINT

TEST(Sieve, QuadraticSieve)
{
    long_int n("59469489332848408438249254427481121839977");  // 338555568168236555657 * 175656509371887105761
    FactorSet factor = QuadraticSieveFactorization::Factorize(n);
    ASSERT_EQ(factor.size(), 2);
    ASSERT_EQ(Eval(factor), n);
}

TEST(Sieve, QuadraticSieveWithConfig)
{
    long_int n("4482406424966880742829846540605971439398287609");  // 86738535685150523290199 * 51677220390685710220591
    auto config = Sieve::CreateConfig(n, 260'000'000, 8000, 1);
    FactorSet factor = QuadraticSieveFactorization::Factorize(n, config);
    ASSERT_EQ(factor.size(), 2);
    ASSERT_EQ(Eval(factor), n);
}

};  // namespace
