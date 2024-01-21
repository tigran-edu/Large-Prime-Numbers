#include <gtest/gtest.h>
#include "basic.h"
#include "factorization.h"

TEST(Factorization, Basic)
{
    using boost::multiprecision::cpp_int;
    cpp_int value1("4615021121");
    auto factor = lpn::BasicFactorization(value1);
    ASSERT_TRUE(factor.size() == 1 && factor[value1] == 1);

    cpp_int value2("4615021122");
    factor = lpn::BasicFactorization(value2);
    ASSERT_TRUE(factor.size() == 4);
    cpp_int result = 1;
    for (auto [div, amount] : factor)
        result *= div;
    ASSERT_TRUE(result == value2);
}

TEST(Factorization, Rho)
{
    using boost::multiprecision::cpp_int;

    cpp_int value("5465458763478567834678638");
    auto factor = lpn::RhoFactorization().factorize(value, 2);

    ASSERT_TRUE(factor.size() == 4);

    cpp_int result = 1;
    for (auto [div, amount] : factor)
        result *= div;

    ASSERT_TRUE(result == value);

    value <<= 20;
    value += 1;
    std::cout << value << '\n';
    factor = lpn::RhoFactorization().factorize(value, 2);

    ASSERT_TRUE(factor.size() == 2);

    result = 1;
    for (auto [div, amount] : factor)
        result *= div;

    ASSERT_TRUE(result == value);

}
