#include "rho.hpp"

#include <gtest/gtest.h>

namespace
{
using namespace lpn;  // NOLINT

long_int ComputeNumber(size_t i) { return i * (i + 1) * (i + 2); }

TEST(RhoFactorization, RhoBasic)
{
    long_int value("5465458763478567834678638");
    auto factor = RhoFactorization::Factorize(value, 2);

    ASSERT_TRUE(factor.size() == 2);
    ASSERT_TRUE(Eval(factor) == value);
}

TEST(RhoFactorization, RhoLargeComplexValue)
{
    long_int value("71636382152868291931");

    for (size_t i = 2; i <= 20; ++i)
    {
        value *= i;
    }
    // value now is 39 digit complex number
    auto factor1 = RhoFactorization::Factorize(value, 2);
    auto factor2 = RhoFactorization::Factorize(value, 5);

    // starting point makes difference
    ASSERT_TRUE(factor1 != factor2);
    ASSERT_TRUE(Eval(factor1) == Eval(factor2));
}

TEST(RhoFactorization, Heavy)
{
    for (size_t i = 1000; i < 4000; ++i)
    {
        long_int value = ComputeNumber(i);
        auto factor = RhoFactorization::Factorize(value, 2);

        ASSERT_TRUE(factor.size() == 2);
        ASSERT_TRUE(Eval(factor) == value);
    }
}
};  // namespace
