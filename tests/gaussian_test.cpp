#include "gaussian.hpp"
#include "basic.hpp"
#include "gtest/gtest.h"
#include "boost/random.hpp"

namespace
{

using namespace lpn;  // NOLINT

TEST(GaussianBasic, mask)
{
    std::vector<size_t> primes = {2, 3, 5, 7, 11, 13};
    std::vector<Factor> factors;
    size_t test_size = primes.size() - 1;
    for (size_t i = 0; i < test_size; ++i)
    {
        factors.push_back(std::move(BasicFactorization(primes[i] * primes[i + 1])));
    }
    auto gaus = GaussianBasic(factors, primes);
    auto matrix = std::move(gaus.Solve());
    for (size_t i = 0; i < matrix.size(); ++i)
    {
        boost::dynamic_bitset<> answer(6);
        answer[i] = true;
        answer[5] = true;
        ASSERT_TRUE(answer == matrix[i].mask);
    }
}

TEST(GaussianBasic, participants)
{
    std::vector<size_t> primes = {2, 3, 5, 7, 11, 13};
    std::vector<Factor> factors;
    size_t test_size = primes.size() - 1;
    for (size_t i = 0; i < test_size; ++i)
    {
        factors.push_back(std::move(BasicFactorization(primes[i] * primes[i + 1])));
    }
    auto gaus = GaussianBasic(factors, primes);
    auto matrix = std::move(gaus.Solve());
    for (size_t i = 0; i < matrix.size(); ++i)
    {
        boost::dynamic_bitset<> answer(5);
        for (size_t j = i; j < 5; ++j)
        {
            answer[j] = true;
        }
        ASSERT_TRUE(answer == matrix[i].participants);
    }
}

TEST(GaussianBasic, random)
{
    boost::mt19937 rng(42);
    std::vector<size_t> primes = {2, 3, 5, 7, 11, 13};
    std::vector<Factor> factors;
    size_t test_size = 10;
    for (size_t i = 0; i < test_size; ++i)
    {
        std::vector<size_t> tmp;
        std::sample(primes.begin(), primes.end(), std::back_inserter(tmp), 4, rng);
        factors.push_back(std::move(BasicFactorization(tmp[0] * tmp[0] * tmp[1] * tmp[1] * tmp[2] * tmp[3])));
    }
    auto gaus = GaussianBasic(factors, primes);
    auto matrix = std::move(gaus.Solve());
    for (auto & line : matrix)
    {
        if (line.mask.none())
        {
            return;
        }
    }

    ASSERT_FALSE(true);
}

};  // namespace
