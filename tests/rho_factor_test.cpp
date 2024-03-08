#include "rho_factorization.hpp"

#include <gtest/gtest.h>

namespace
{
using namespace lpn;  // NOLINT

long_int Eval(const Factor & factor)
{
    long_int result = 1;
    for (const auto & [div, amount] : factor)
    {
        result *= FastExponentiation(div, amount);
    }
    return result;
}

TEST(Factorization, RhoBasic)
{
    long_int value("5465458763478567834678638");
    auto factor = RhoFactorization().Factorize(value, 2);

    ASSERT_TRUE(factor.size() == 4);
    ASSERT_TRUE(Eval(factor) == value);
}

TEST(Factorization, RhoLargeComplexValue)
{
    long_int value("71636382152868291931");

    for (size_t i = 2; i <= 20; ++i)
    {
        value *= i;
    }
    // value now is 39 digit complex number
    auto factor1 = RhoFactorization().Factorize(value, 2);
    auto factor2 = RhoFactorization().Factorize(value, 5);

    // starting point makes difference
    ASSERT_TRUE(factor1 != factor2);
    ASSERT_TRUE(Eval(factor1) == Eval(factor2));
}
};  // namespace
