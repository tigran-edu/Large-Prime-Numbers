#include "factorization.h"

#include <gtest/gtest.h>

#include "basic.h"

namespace
{
using long_int = boost::multiprecision::cpp_int;

long_int eval(const lpn::Factor & factor)
{
    long_int result = 1;
    for (const auto & [div, amount] : factor)
    {
        result *= lpn::FastExponentiation(div, amount);
    }
    return result;
}

TEST(Factorization, BasicPrime)
{
    long_int value("9111657031");
    auto factor = lpn::BasicFactorization(value);
    ASSERT_TRUE(factor.size() == 1 && factor[value] == 1);
}

TEST(Factorization, BasicComplex)
{
    long_int value("1307674368000");
    auto factor = lpn::BasicFactorization(value);
    ASSERT_TRUE(factor.size() == 6);
    long_int result = 1;
    for (const auto & [div, amount] : factor)
    {
        result *= lpn::FastExponentiation(div, amount);
    }
    ASSERT_TRUE(result == value);
}

TEST(Factorization, RhoBasic)
{
    long_int value("5465458763478567834678638");
    auto factor = lpn::RhoFactorization().factorize(value, 2);

    ASSERT_TRUE(factor.size() == 4);
    ASSERT_TRUE(eval(factor) == value);
}

TEST(Factorization, RhoLargeComplexValue)
{
    long_int value("71636382152868291931");

    for (size_t i = 2; i <= 20; ++i)
    {
        value *= i;
    }
    // value now is 39 digit complex number
    auto factor1 = lpn::RhoFactorization().factorize(value, 2);
    auto factor2 = lpn::RhoFactorization().factorize(value, 5);

    // starting point makes difference
    ASSERT_TRUE(factor1 != factor2);
    ASSERT_TRUE(eval(factor1) == eval(factor2));
}
};  // namespace
