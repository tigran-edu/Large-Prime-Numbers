#pragma once

#include "aliases.hpp"
#include "base.hpp"
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
    static Solution Solve(const long_int & n, size_t factor_size);

   private:
    static void AddFactor(const std::optional<FactorSet> & factor, const State & state, Solution & solution);
};

class ContinuedFractionsFactorization : private FactorizationBase
{
   public:
    static FactorSet Factorize(const long_int & n);
    static FactorSet Factorize(const long_int & n, size_t factor_size);

   private:
    static constexpr size_t kBasicFactorSize = 4000;
};

};  // namespace lpn
