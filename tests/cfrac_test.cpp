#include "cfrac.hpp"

#include <gtest/gtest.h>

namespace
{

using namespace lpn;  // NOLINT

TEST(CFRAC, ContinuedFractions)
{
    long_int n("59469489332848408438249254427481121839977");  // 338555568168236555657 * 175656509371887105761
    FactorSet factor = ContinuedFractionsFactorization::Factorize(n);
    ASSERT_EQ(factor.size(), 2);
    ASSERT_EQ(Eval(factor), n);
}

TEST(CFRAC, ContinuedFractionsWithConfig)
{
    long_int n("275496260473012513310855787717477881849");  // 23074736938723048343 * 11939302328976341743
    FactorSet factor = ContinuedFractionsFactorization::Factorize(n, 500);
    ASSERT_EQ(factor.size(), 2);
    ASSERT_EQ(Eval(factor), n);
}

};  // namespace
