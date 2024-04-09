#pragma once

#include <optional>

#include "gaussian.hpp"
#include "base.hpp"

namespace lpn
{

class Sieve
{
   private:
    struct BasicConfig
    {
        static constexpr size_t kDefSegmentSize = 50'000'000;
        static constexpr size_t kDefFactorSize = 2000;
        static constexpr float kExpansionRate = 1.5;
    };

   public:
    struct Config
    {
        friend Sieve;

       private:
        Config(size_t segment_size, size_t factor_size);

        void ComputeCloseness(float expansion);

        void Reserve();

        void ComputeTarget(const long_int & n);

        void ComputeSmallPrimes(const long_int & n);

        size_t segment_size_;
        size_t factor_size_;
        float target_;
        float closeness_;
        std::vector<size_t> primes_;
        std::vector<long_int> congruences_;
    };

    static Solution Solve(const long_int & n, const Config & config);
    static Config CreateConfig(const long_int & n);
    static Config CreateConfig(const long_int & n, size_t segment_size, size_t factor_size, float expansion_rate);

   private:
    static void AddLogAtSpecificPositions(size_t i, size_t j, size_t p, std::vector<float> & logs);
    static long_int ComputeTargetFunction(const long_int & r, const long_int & n, size_t i);
};

class QuadraticSieveFactorization : private FactorizationBase
{
   public:
    static FactorSet Factorize(const long_int & n);
    static FactorSet Factorize(const long_int & n, const Sieve::Config & config);
};

};  // namespace lpn
