#include "qs.hpp"
#include "boost/random.hpp"

#include <gtest/gtest.h>
#include <vector>
#include <fstream>

namespace
{

using namespace lpn;  // NOLINT

bool Check(const long_int & x, const long_int & p, const size_t i) { return (x * x) % p == i * i; }

TEST(Sieve, QuadraticCongruencesRandom)
{
    boost::mt19937 rng(42);
    {
        long_int p("706546294168342682401061556019");  // p % 4 == 3
        for (size_t i = 1; i < 20000; ++i)
        {
            long_int value = p * rng() + i * i;
            ASSERT_TRUE(Check(QuadraticCongruences::SolvingQuadraticCongruences(value, p), p, i));
        }
    }

    {
        long_int p("850089435441534261939703688213");  // p % 8 == 5
        for (size_t i = 1; i < 20000; ++i)
        {
            long_int value = p * rng() + i * i;
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

// TEST(Sieve, QuadraticSieve)
// {
//     const std::array<long_int, 4> primes = {391708601274539, 136564985324737, 261009627419033, 630274293935959};
//     for (const long_int & p1 : primes)
//     {
//         for (const long_int & p2 : primes)
//         {
//             long_int value = p1 * p2;
//             Factor factor = QuadraticSieve::Factorize(value);
//             for (auto & [key, value] : factor)
//             {
//                 std::cout << key << " " << value << std::endl;
//             }
//             std::cout << factor[p1] << " " << factor[p2] << std::endl;
//             ASSERT_TRUE(factor[p1] == 1 && factor[p2] == 1 && factor.size() == 2);
//         }
//     }
// }

};  // namespace