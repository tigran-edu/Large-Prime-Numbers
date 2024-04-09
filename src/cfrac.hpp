#pragma once

#include "aliases.hpp"
#include "basic.hpp"
#include "gaussian.hpp"

namespace lpn
{

class ContinuedFractions
{
   private:
    struct State
    {
        explicit State(const long_int & n);

        void Update();

        long_int n;
        long_int sqrt_n;
        long_int a;
        long_int b0;
        long_int b1;
        long_int c0;
        long_int c1;
        long_int p0;
        long_int p1;
    };

   public:
    struct Solution
    {
        FactorSets factors;
        std::vector<long_int> xs;
        std::vector<size_t> primes;
    };

   public:
    static Solution Solve(const long_int & n, size_t factor_size);
};

class ContinuedFractionsFactorization
{
   private:
    using Matrix = GaussianBasic::Matrix;
    using Line = GaussianBasic::Line;
    using Solution = ContinuedFractions::Solution;

   public:
    static FactorSet Factorize(const long_int & n);
    static FactorSet Factorize(const long_int & n, size_t factor_size);

   private:
    static long_int ComputeX(const Solution & solution, const std::vector<size_t> & positions, const long_int & n);
    static long_int ComputeY(const Solution & solution, const std::vector<size_t> & positions, const long_int & n);
    static bool IsPerfectSquare(const Line & line);
    static FactorSet FindFactor(const Solution & solution, const Matrix & matrix, const long_int & n);

   private:
    static constexpr size_t kBasicFactorSize = 4000;
};

};  // namespace lpn
